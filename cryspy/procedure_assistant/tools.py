"""Tools of CrysPy for AI assistants and other programs.

Every tool takes the CrysPy object and keyword arguments with plain values
(str, float, bool) and returns text, so it can be called by a language
model through JSON. CrysPy Editor and any other front end use the same
tools:

    >>> import cryspy
    >>> obj = cryspy.load_file("main.rcif")
    >>> cryspy.call_assistant_tool(obj, "list_parameters", {"pattern": "cell"})
    {'text': '...', 'flag_object_changed': False}

assistant_tools() describes the tools in the JSON schema format used by
chat models ("tools" of OpenAI and llama.cpp).

A tool that lacks information only the user knows (unit cell, wavelength,
...) does not fail: the result has the key 'missing', a list of the
arguments to ask the user (as "elicitation" of the Model Context Protocol).
The program asks them, checks each answer with convert_assistant_answer and
calls the tool again with the answers added to the arguments:

    >>> result = cryspy.call_assistant_tool(obj, name, arguments)
    >>> for d_missing in result["missing"]:
    ...     arguments[d_missing["name"]] = cryspy.convert_assistant_answer(
    ...         d_missing, input(d_missing["question"]))
    >>> result = cryspy.call_assistant_tool(obj, name, arguments)

A required argument is asked from the user when its description has
"question"; a tool asks for other values by raising ArgumentsNeeded.
"""
import re
from typing import Callable, Dict, List, Optional

import numpy

from cryspy.B_parent_classes.cl_3_data import DataN
from cryspy.B_parent_classes.cl_4_global import GlobalN
from cryspy.C_item_loop_classes.cl_1_atom_site import AtomSite
from cryspy.E_data_classes.cl_1_crystal import Crystal

from cryspy.procedure_assistant.help import get_help
from cryspy.procedure_assistant.procedures import (
    describe_procedure, describe_procedures, find_cif_objects, take_procedure_parameters,
    take_procedures)
from cryspy.procedure_assistant.history import (
    form_metrics, format_metrics, keep_checkpoint, keep_snapshot, restore_snapshot, take_history)
from cryspy.procedure_assistant.recipes import recipes_as_text, take_recipe
from cryspy.procedure_assistant.parameters import (
    find_parameter_records, suggest_parameter_names, take_parameter_record,
    take_parameter_records)

MAX_PARAMETERS_IN_LIST = 200
CORRELATION_LIMIT = 0.7
MISSING_KEYS = ("type", "items", "minItems", "maxItems", "enum", "aliases",
                "minimum", "maximum", "hint", "proposal", "s_proposal")
NULL_WORDS = ("none", "null", "identity", "-")
RE_NUMBER = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?")


class ArgumentsNeeded(Exception):
    """Raised by a tool for the values to ask the user.

    l_missing: dicts with "name" (argument of the tool), "question" (what
    is asked, e.g. "the unit cell parameters") and optional JSON schema
    keys of the value ("type", "items", "minItems", "maxItems", "enum"),
    "hint" (format of the answer), "proposal" (a value known already, the
    user confirms it), "s_proposal" (its text)."""

    def __init__(self, l_missing: List[dict], message: str = ""):
        super().__init__(message)
        self.l_missing = l_missing
        self.message = message


class AssistantTool:
    """A function with the description of its arguments for a chat model."""

    def __init__(self, name: str, description: str, function: Callable,
                 parameters: Optional[Dict[str, dict]] = None,
                 required: tuple = (), flag_object_changed: bool = False,
                 flag_object_needed: bool = True,
                 flag_extra_arguments: bool = False):
        self.name = name
        self.description = description
        self.function = function
        self.parameters = parameters if parameters is not None else {}
        self.required = required
        self.flag_object_changed = flag_object_changed
        self.flag_object_needed = flag_object_needed
        # arguments not in parameters are given to the function as they are
        self.flag_extra_arguments = flag_extra_arguments

    def to_dict(self) -> dict:
        """Description in the 'tools' format of chat completion APIs."""
        return {
            "type": "function",
            "function": {
                "name": self.name,
                "description": self.description,
                "parameters": {
                    "type": "object",
                    "properties": self.parameters,
                    "required": list(self.required),
                },
            },
        }

    def to_short_text(self) -> str:
        """One line for the prompt of a small model."""
        l_argument = []
        for name, d_param in self.parameters.items():
            s_type = d_param.get("type", "string")
            if name in self.required:
                l_argument.append(f"{name:}: {s_type:}")
            else:
                l_argument.append(f"{name:}?: {s_type:}")
        return f"- {self.name:}({', '.join(l_argument):}): {self.description:}"


def format_value(value) -> str:
    if isinstance(value, (float, numpy.floating)):
        return f"{float(value):.6g}"
    return str(value)


def format_dictionary_name(name) -> str:
    """Name of a parameter from the dictionary pipeline as text.

    ('crystal_pbso4', 'atom_fract_xyz', (0, 1)) -> 'crystal_pbso4.atom_fract_xyz[0, 1]'
    """
    if isinstance(name, tuple) and len(name) == 3 and isinstance(name[0], str):
        index = ", ".join(str(int(i)) for i in numpy.atleast_1d(name[2]))
        return f"{name[0]:}.{name[1]:}[{index:}]"
    return str(name)


def form_refined_parameter_table(cryspy_object) -> str:
    l_record = [record for record in take_parameter_records(cryspy_object)
                if record.flag_refinement]
    ls_out = []
    for record in l_record:
        if record.sigma != 0.:
            ls_out.append(f"{record.path:} = {format_value(record.value):} "
                          f"+- {format_value(record.sigma):}")
        else:
            ls_out.append(f"{record.path:} = {format_value(record.value):}")
    return "\n".join(ls_out)


def format_rhochi_result(result: dict, cryspy_object) -> str:
    """Chi-square, refined parameters and strong correlations as text."""
    ls_out = []
    if "chi_sq" in result:
        s_chi = f"chi_sq = {float(result['chi_sq']):.5g}"
        if "n_point" in result and float(result["n_point"]) > 0:
            s_chi += (f", points = {int(result['n_point']):}, chi_sq/n = "
                      f"{float(result['chi_sq']) / float(result['n_point']):.4g}")
        ls_out.append(s_chi)
    table = form_refined_parameter_table(cryspy_object)
    if table != "":
        ls_out.append("Refined parameters:\n" + table)
    names = result.get("parameter_name", [])
    matrix = result.get("correlation_matrix", None)
    if matrix is not None and len(names) == len(matrix):
        l_pair = []
        for ind_1 in range(len(names)):
            for ind_2 in range(ind_1 + 1, len(names)):
                value = float(matrix[ind_1][ind_2])
                if abs(value) >= CORRELATION_LIMIT:
                    l_pair.append((abs(value), value, ind_1, ind_2))
        l_pair.sort(reverse=True)
        if len(l_pair) > 0:
            ls_out.append(f"Correlations |r| >= {CORRELATION_LIMIT:}:")
            for _, value, ind_1, ind_2 in l_pair[:20]:
                ls_out.append(
                    f"{format_dictionary_name(names[ind_1]):} - "
                    f"{format_dictionary_name(names[ind_2]):}: {value:.2f}")
        else:
            ls_out.append(f"No correlations with |r| >= {CORRELATION_LIMIT:}.")
    return "\n".join(ls_out)


def tool_get_help(cryspy_object, query: str) -> str:
    return get_help(query)


def tool_summarize_object(cryspy_object) -> str:
    ls_out = []
    if isinstance(cryspy_object, GlobalN):
        l_item = cryspy_object.items
    else:
        l_item = [cryspy_object]
    for item in l_item:
        if isinstance(item, DataN):
            l_block = [block.get_name() for block in item.items]
            ls_out.append(f"data block '{item.data_name:}' "
                          f"({type(item).__name__:}): {', '.join(l_block):}")
        else:
            ls_out.append(f"item {item.get_name():}")
    if len(ls_out) == 0:
        ls_out.append("The object is empty.")
    l_record = take_parameter_records(cryspy_object)
    n_refined = sum(1 for record in l_record if record.flag_refinement)
    ls_out.append(f"Refinable parameters: {len(l_record):}, "
                  f"marked for refinement: {n_refined:}.")
    return "\n".join(ls_out)


def form_no_match_error(cryspy_object, pattern: str) -> str:
    """Error text with the names the model can use in the next call."""
    l_name = suggest_parameter_names(cryspy_object, pattern)
    s_name = ", ".join(l_name) if len(l_name) > 0 else "none"
    return (f"Error: no parameters match '{pattern:}'. Parameters of this "
            f"object look like: {s_name:}.")


def tool_list_parameters(cryspy_object, pattern: str = "",
                         refined_only: bool = False) -> str:
    l_record = find_parameter_records(
        cryspy_object, pattern, flag_refined_only=refined_only)
    if len(l_record) == 0:
        if not refined_only or len(find_parameter_records(
                cryspy_object, pattern)) == 0:
            return form_no_match_error(cryspy_object, pattern)
        if pattern.strip() == "":
            return "No parameters are marked for refinement."
        return "No parameters matching the pattern are marked for refinement."
    ls_out = ["path | value | sigma | refined | constrained"]
    for record in l_record[:MAX_PARAMETERS_IN_LIST]:
        ls_out.append(
            f"{record.path:} | {format_value(record.value):} | "
            f"{format_value(record.sigma):} | {record.flag_refinement:} | "
            f"{record.flag_constraint:}")
    if len(l_record) > MAX_PARAMETERS_IN_LIST:
        ls_out.append(f"... {len(l_record) - MAX_PARAMETERS_IN_LIST:} more, "
                      "give a narrower pattern.")
    return "\n".join(ls_out)


def tool_refine_parameters(cryspy_object, pattern: str) -> str:
    if pattern.strip() == "":
        return "Error: give a pattern of the parameters to refine."
    l_record = find_parameter_records(cryspy_object, pattern)
    if len(l_record) == 0:
        return form_no_match_error(cryspy_object, pattern)
    l_free, l_constrained = [], []
    for record in l_record:
        if record.set_refinement(True):
            l_free.append(record.path)
        else:
            l_constrained.append(record.path)
    ls_out = []
    if len(l_free) > 0:
        ls_out.append("Marked for refinement:\n" + "\n".join(l_free))
    if len(l_constrained) > 0:
        ls_out.append("Not refinable because of constraints:\n"
                      + "\n".join(l_constrained))
    if len(l_free) == 0:
        ls_out.insert(0, "Error: none of the matching parameters can be "
                         "refined.")
    return "\n".join(ls_out)


def tool_fix_parameters(cryspy_object, pattern: str = "") -> str:
    if pattern.strip() in ("", "*", "all"):
        l_record = find_parameter_records(
            cryspy_object, "", flag_refined_only=True)
        n_variable = len(cryspy_object.get_variable_names())
        # also the refined parameters without a path (e.g. of nested loops)
        cryspy_object.fix_variables()
        if n_variable == 0:
            return "No parameters were marked for refinement: all are fixed."
        ls_out = [f"All {n_variable:} refined parameters are fixed now."]
        ls_out.extend(record.path for record in l_record)
        return "\n".join(ls_out)
    l_match = find_parameter_records(cryspy_object, pattern)
    if len(l_match) == 0:
        return form_no_match_error(cryspy_object, pattern)
    l_record = [record for record in l_match if record.flag_refinement]
    for record in l_record:
        record.set_refinement(False)
    if len(l_record) == 0:
        return (f"The {len(l_match):} parameters matching '{pattern:}' "
                "were not refined: nothing to fix.")
    return "Fixed:\n" + "\n".join(record.path for record in l_record)


def tool_set_parameter(cryspy_object, path: str, value: float) -> str:
    record = take_parameter_record(cryspy_object, path)
    if record is None:
        l_record = find_parameter_records(cryspy_object, path)
        if len(l_record) == 1:
            record = l_record[0]
        elif len(l_record) > 1:
            return ("Error: the path is ambiguous, give one of:\n"
                    + "\n".join(record.path for record in l_record[:20]))
        else:
            return form_no_match_error(cryspy_object, path)
    if record.flag_constraint:
        return (f"Error: {record.path:} is constrained by symmetry and "
                "cannot be set.")
    old_value = record.value
    record.set_value(value)
    return (f"{record.path:}: {format_value(old_value):} -> "
            f"{format_value(record.value):}")


def tool_calculate(cryspy_object) -> dict:
    from cryspy.procedure_rhochi.rhochi import rhochi_no_refinement
    history = take_history(cryspy_object)
    result = rhochi_no_refinement(cryspy_object)
    d_metrics = form_metrics(result, cryspy_object, chi_sq_previous=history.chi_sq)
    keep_checkpoint(cryspy_object, d_metrics.get("chi_sq"))
    return {"text": form_calculation_text("Calculation finished.", result, cryspy_object, d_metrics),
            "metrics": d_metrics}


def form_calculation_text(title: str, result: dict, cryspy_object, d_metrics: dict) -> str:
    # the change of chi_sq and the warnings first: a long table of parameters may be cut
    ls_out = [title, format_metrics(d_metrics), format_rhochi_result(result, cryspy_object)]
    return "\n".join(line for line in ls_out if line != "")


def tool_refine(cryspy_object) -> dict:
    from cryspy.procedure_rhochi.rhochi import rhochi_rietveld_refinement
    if len(cryspy_object.get_variable_names()) == 0:
        return {"text": "Error: no parameters are marked for refinement. Use refine_parameters first."}
    history = take_history(cryspy_object)
    d_last = {}  # the last iteration of the minimizer: chi_sq and the number of points
    keep_snapshot(cryspy_object)
    try:
        result = rhochi_rietveld_refinement(cryspy_object, callback=d_last.update)
    except Exception:
        history.l_snapshot.pop()  # nothing is changed: nothing to undo
        raise
    d_metrics = form_metrics(result, cryspy_object, n_point=d_last.get("n_point"),
                             chi_sq_previous=history.chi_sq)
    if "n_point" in d_metrics:
        result = dict(result, n_point=d_metrics["n_point"])
    keep_checkpoint(cryspy_object, d_metrics.get("chi_sq"))
    return {"text": form_calculation_text("Refinement finished.", result, cryspy_object, d_metrics),
            "metrics": d_metrics}


def tool_undo_refinement(cryspy_object) -> dict:
    from cryspy.procedure_rhochi.rhochi import rhochi_no_refinement
    snapshot = restore_snapshot(cryspy_object)
    if snapshot is None:
        return {"text": "Error: there is no refinement to undo."}
    # the calculated data of the object are those of the undone refinement
    result = rhochi_no_refinement(cryspy_object)
    d_metrics = form_metrics(result, cryspy_object)
    keep_checkpoint(cryspy_object, d_metrics.get("chi_sq"))
    n_refined = sum(1 for _, flag in snapshot["flags"] if flag)
    text = (f"The parameters and their refinement flags are returned to the state before the "
            f"last refinement ({n_refined:} marked for refinement).\n"
            + format_rhochi_result(result, cryspy_object))
    if d_metrics["warnings"]:
        text += "\n" + "\n".join(f"Warning: {warning:}" for warning in d_metrics["warnings"])
    return {"text": text, "metrics": d_metrics}


def tool_list_recipes(cryspy_object) -> str:
    return ("Recipes (typical algorithms of analysis; get_recipe gives the steps):\n"
            + recipes_as_text())


def tool_get_recipe(cryspy_object, name: str) -> str:
    recipe = take_recipe(name)
    if recipe is None:
        return f"Error: recipe '{name:}' is not known. Known recipes:\n{recipes_as_text():}"
    return recipe.to_text()


def take_data_block(cryspy_object, name: str):
    name_sh = name.strip().lower()
    for item in cryspy_object.items:
        if isinstance(item, DataN):
            if name_sh in (str(item.data_name).lower(), item.get_name()):
                return item
    return None


def tool_show_model(cryspy_object, block: str = "",
                    max_length: int = 8000) -> str:
    if block.strip() == "":
        text = cryspy_object.to_cif()
    else:
        data = take_data_block(cryspy_object, block)
        if data is None:
            return f"Error: data block '{block:}' is not found."
        text = data.to_cif()
    if len(text) > max_length:
        text = text[:max_length] + "\n[... the rest is cut ...]"
    return text


def tool_add_atom(cryspy_object, crystal: str, label: str, type_symbol: str,
                  fract_x: float, fract_y: float, fract_z: float,
                  occupancy: float = 1., adp_type: str = "",
                  b_iso_or_equiv: float = 0.) -> str:
    data = take_data_block(cryspy_object, crystal)
    if not isinstance(data, Crystal):
        return f"Error: crystal '{crystal:}' is not found."
    atom_site = data.atom_site
    if label in [str(item.label) for item in atom_site.items]:
        return f"Error: atom '{label:}' already exists in '{crystal:}'."
    if adp_type == "":
        # the same kind of isotropic parameter as the other atoms
        l_adp_type = [item.adp_type for item in atom_site.items
                      if item.is_attribute("adp_type")]
        adp_type = "Uiso" if "Uiso" in l_adp_type else "Biso"
    if adp_type not in ("Biso", "Uiso"):
        return ("Error: adp_type should be 'Biso' or 'Uiso'; "
                "anisotropic parameters are set in atom_site_aniso.")
    atom = AtomSite(label=label, type_symbol=type_symbol, fract_x=fract_x,
                    fract_y=fract_y, fract_z=fract_z, occupancy=occupancy,
                    adp_type=adp_type)
    if adp_type == "Biso":
        atom.b_iso_or_equiv = b_iso_or_equiv
    else:
        atom.u_iso_or_equiv = b_iso_or_equiv
    # keep the same columns as the other atoms, so the loop stays printable
    for name in ("b_iso_or_equiv", "u_iso_or_equiv"):
        if not atom.is_attribute(name) and all(
                item.is_attribute(name) for item in atom_site.items):
            setattr(atom, name, 0.)
    atom_site.items.append(atom)
    data.apply_constraints()
    return (f"Atom {label:} ({type_symbol:}) is added to {crystal:} at "
            f"({fract_x:}, {fract_y:}, {fract_z:}).")


def take_hkl(h, k, l, hkl: str):
    """numpy array (3, 1) of the indices given as h, k, l or as '5 0 3'."""
    if hkl.strip() != "":
        l_index = re.findall(r"-?\d+(?:\.\d*)?", hkl)
        if len(l_index) != 3:
            raise ValueError(f"hkl '{hkl:}' should have three indices.")
        h, k, l = (float(index) for index in l_index)
    if h is None or k is None or l is None:
        raise ValueError("give the indices h, k, l (or hkl as '5 0 3').")
    return numpy.array([[h], [k], [l]], dtype=float)


def tool_calculate_structure_factor(
        cryspy_object, h: float = None, k: float = None, l: float = None,
        hkl: str = "", crystal: str = "") -> str:
    try:
        index_hkl = take_hkl(h, k, l, hkl)
    except ValueError as error:
        return f"Error: {error:}"
    if crystal.strip() == "":
        l_crystal = [item for item in cryspy_object.items
                     if isinstance(item, Crystal)]
    else:
        l_crystal = [take_data_block(cryspy_object, crystal)]
        if not isinstance(l_crystal[0], Crystal):
            return f"Error: crystal '{crystal:}' is not found."
    if len(l_crystal) == 0:
        return "Error: the object has no crystal."
    s_hkl = " ".join(f"{index:g}" for index in index_hkl[:, 0])
    ls_out = [f"Nuclear structure factor F_N of the reflection ({s_hkl:}), "
              "in 10^-12 cm (per unit cell, with the Debye-Waller factor):"]
    for item in l_crystal:
        f_nucl = complex(numpy.asarray(item.calc_f_nucl(index_hkl)).ravel()[0])
        sthovl = float(numpy.asarray(item.cell.calc_sthovl(
            *index_hkl[:, 0])).ravel()[0])
        ls_out.append(
            f"{item.data_name:}: F_N = {f_nucl.real:.5f} {f_nucl.imag:+.5f}i, "
            f"|F_N| = {abs(f_nucl):.5f}, |F_N|^2 = {abs(f_nucl)**2:.5f}, "
            f"sin(theta)/lambda = {sthovl:.5f} 1/Å")
    return "\n".join(ls_out)


def tool_add_crystal_from_cod(cryspy_object, formula: str = "",
                              cod_id: str = "", name: str = "") -> str:
    from cryspy.H_functions_global.cod import add_crystal_from_cod
    if formula.strip() == "" and str(cod_id).strip() == "":
        return "Error: give the chemical formula or the COD number."
    try:
        crystal = add_crystal_from_cod(
            cryspy_object, formula=formula, cod_id=str(cod_id),
            data_name=name)
    except ValueError as error:
        return f"Error: {error:}"
    cell = crystal.cell
    return (f"Crystal {crystal.data_name:} is added from COD: space group "
            f"{crystal.space_group.name_hm_alt:}, a, b, c = {cell.length_a:g}, "
            f"{cell.length_b:g}, {cell.length_c:g} Å, "
            f"{len(crystal.atom_site.items):} atoms.")


def tool_list_procedures(cryspy_object) -> str:
    ls_out = []
    for d_procedure in describe_procedures():
        ls_param = []
        for d_param in d_procedure["parameters"]:
            if d_param["kind"] == "global":
                continue
            s_type = d_param.get("type") or "|".join(d_param.get("classes", []))
            s_default = f"={d_param['default']!r:}" if "default" in d_param else ""
            ls_param.append(f"{d_param['name']:}: {s_type:}{s_default:}")
        ls_out.append(f"- {d_procedure['name']:}({', '.join(ls_param):}): {d_procedure['summary']:}")
    return "\n".join(ls_out)


def form_procedure_result(name: str, result, cryspy_object) -> str:
    from cryspy.B_parent_classes.cl_1_item import ItemN
    from cryspy.B_parent_classes.cl_2_loop import LoopN
    if isinstance(result, dict) and "chi_sq" in result:
        return f"{name:} finished.\n" + format_rhochi_result(
            result, cryspy_object)
    if result is None:
        return f"{name:} finished."
    if isinstance(result, (ItemN, LoopN, DataN, GlobalN)):
        text = result.to_cif()
    elif isinstance(result, (list, tuple)) and len(result) > 0 and all(
            isinstance(item, (ItemN, LoopN, DataN)) for item in result):
        text = "\n\n".join(item.to_cif() for item in result)
    else:
        text = str(result)
    if len(text) > 4000:
        text = text[:4000] + "\n[... the rest is cut ...]"
    return f"{name:} finished.\n{text:}"


def tool_run_procedure(cryspy_object, name: str, **arguments) -> str:
    """Run a procedure: the object, and the items of the object found by
    their classes, are given at once; values without a default and items
    found more than once are asked from the user (ArgumentsNeeded) with
    the questions made from the docstring of the procedure."""
    d_function = take_procedures()
    if name not in d_function:
        return (f"Error: procedure '{name:}' is not known. Available:\n"
                + tool_list_procedures(cryspy_object))
    function = d_function[name]
    d_param_type = {param: types for param, _, types, _ in take_procedure_parameters(function)}
    kwargs, l_missing = {}, []
    for d_param in describe_procedure(function)["parameters"]:
        param, kind = d_param["name"], d_param["kind"]
        d_ask = {"name": param, "question": d_param.get("question", param)}
        if kind == "global":
            kwargs[param] = cryspy_object
        elif kind == "basic":
            if param in arguments:
                kwargs[param] = convert_argument(arguments[param], d_param)
            elif "default" not in d_param:
                l_missing.append(dict(d_ask, type=d_param["type"], hint=param))
        elif kind == "cif":
            d_found = find_cif_objects(cryspy_object, d_param_type[param])
            s_type = " or ".join(d_param["classes"])
            if param in arguments and str(arguments[param]) in d_found:
                kwargs[param] = d_found[str(arguments[param])]
            elif len(d_found) == 1:
                kwargs[param] = list(d_found.values())[0]
            elif len(d_found) == 0:
                return f"Error: {name:} needs {s_type:} ({param:}), the object has none."
            else:
                l_missing.append(dict(d_ask, type="string", enum=list(d_found), hint=", ".join(d_found)))
    if l_missing:
        raise ArgumentsNeeded(l_missing)
    return form_procedure_result(name, function(**kwargs), cryspy_object)


S_PATH = ("Parameter path or pattern: 'crystal_pbso4.cell.length_a', "
          "'atom_site[Fe].fract_x', 'b_iso', '*.pd_background[*].intensity'. "
          "Without '*' any path containing the text matches; the labels "
          "may be left out ('atom_site.fract_x'). Several patterns are "
          "separated by ','. Words work too: 'atomic positions', "
          "'scale factor', 'polarization', 'unit cell', 'adp'.")

L_ASSISTANT_TOOLS: List[AssistantTool] = [
    AssistantTool(
        "get_help",
        "Search CrysPy help: how to describe something in rcif "
        "(thermal vibrations, multipoles, extinction, twins, powder, ...) "
        "or the description of a category.",
        tool_get_help,
        {"query": {"type": "string",
                   "description": "Question, topic or category name."}},
        required=("query", ), flag_object_needed=False),
    AssistantTool(
        "summarize_object",
        "Short overview of the loaded object: data blocks, their categories "
        "and the number of refined parameters.",
        tool_summarize_object),
    AssistantTool(
        "list_parameters",
        "List refinable parameters with values, sigmas and flags.",
        tool_list_parameters,
        {"pattern": {"type": "string", "description": S_PATH},
         "refined_only": {"type": "boolean",
                          "description": "Only the parameters marked for refinement."}}),
    AssistantTool(
        "refine_parameters",
        "Mark the parameters matching the pattern for refinement.",
        tool_refine_parameters,
        {"pattern": {"type": "string", "description": S_PATH}},
        required=("pattern", ), flag_object_changed=True),
    AssistantTool(
        "fix_parameters",
        "Switch off refinement of the parameters matching the pattern "
        "(all parameters if the pattern is empty).",
        tool_fix_parameters,
        {"pattern": {"type": "string", "description": S_PATH}},
        flag_object_changed=True),
    AssistantTool(
        "set_parameter",
        "Set the value of one parameter.",
        tool_set_parameter,
        {"path": {"type": "string", "description": S_PATH},
         "value": {"type": "number", "description": "New value."}},
        required=("path", "value"), flag_object_changed=True),
    AssistantTool(
        "calculate",
        "Calculate the model without refinement: chi-square and calculated data.",
        tool_calculate, flag_object_changed=True),
    AssistantTool(
        "refine",
        "Refine the parameters marked for refinement (RhoChi); returns "
        "chi-square and its change since the previous calculation, refined "
        "values, strong correlations and warnings (unphysical values).",
        tool_refine, flag_object_changed=True),
    AssistantTool(
        "undo_refinement",
        "Return the parameters and their refinement flags to the state "
        "before the last refinement (e.g. when chi-square got worse or a "
        "value is unphysical).",
        tool_undo_refinement, flag_object_changed=True),
    AssistantTool(
        "list_recipes",
        "List the recipes: typical algorithms of analysis (e.g. the "
        "step-by-step structure refinement).",
        tool_list_recipes, flag_object_needed=False),
    AssistantTool(
        "get_recipe",
        "Steps of a recipe: for each step the actions, what counts as "
        "success and what to do otherwise. Follow them one action at a time.",
        tool_get_recipe,
        {"name": {"type": "string", "description": "Name of the recipe from list_recipes."}},
        required=("name", ), flag_object_needed=False),
    AssistantTool(
        "show_model",
        "Show the object, or one data block, in rcif format.",
        tool_show_model,
        {"block": {"type": "string",
                   "description": "Name of a data block; empty for all."}}),
    AssistantTool(
        "add_atom",
        "Add an atom to the atom_site loop of a crystal.",
        tool_add_atom,
        {"crystal": {"type": "string", "description": "Name of the crystal data block."},
         "label": {"type": "string"},
         "type_symbol": {"type": "string", "description": "Element or ion, e.g. Fe3+."},
         "fract_x": {"type": "number"},
         "fract_y": {"type": "number"},
         "fract_z": {"type": "number"},
         "occupancy": {"type": "number"},
         "adp_type": {"type": "string", "enum": ["Biso", "Uiso"],
                      "description": "By default as the other atoms."},
         "b_iso_or_equiv": {"type": "number",
                            "description": "Isotropic B (or U for Uiso), Å²."}},
        required=("crystal", "label", "type_symbol", "fract_x", "fract_y",
                  "fract_z"),
        flag_object_changed=True),
    AssistantTool(
        "calculate_structure_factor",
        "Nuclear structure factor F_N (complex, |F_N|, |F_N|^2) of one "
        "reflection hkl for the crystals of the object.",
        tool_calculate_structure_factor,
        {"h": {"type": "number"}, "k": {"type": "number"},
         "l": {"type": "number"},
         "hkl": {"type": "string",
                 "description": "The indices as text, e.g. '5 0 3'; "
                                "instead of h, k, l."},
         "crystal": {"type": "string",
                     "description": "Name of the crystal; empty for all."}}),
    AssistantTool(
        "add_crystal_from_cod",
        "Download a crystal structure from the Crystallography Open "
        "Database (COD, needs internet) and add it to the object as a "
        "crystal data block.",
        tool_add_crystal_from_cod,
        {"formula": {"type": "string",
                     "description": "Chemical formula, e.g. 'Al', 'TiO2'; "
                                    "the first COD entry is taken."},
         "cod_id": {"type": "string",
                    "description": "COD number, e.g. '9008460'; instead of "
                                   "the formula."},
         "name": {"type": "string",
                  "description": "Name of the new data block; by default "
                                 "the formula."}},
        flag_object_changed=True),
    AssistantTool(
        "list_procedures",
        "List CrysPy procedures that can be run on the object "
        "(MEM, simulations, reports, ...).",
        tool_list_procedures),
    AssistantTool(
        "run_procedure",
        "Run a CrysPy procedure from list_procedures by its name; its "
        "parameters are further arguments (e.g. \"iteration_max\": 500). "
        "Missing values are asked from the user; results that are rcif "
        "objects are given in rcif format.",
        tool_run_procedure,
        {"name": {"type": "string"}}, required=("name", ),
        flag_object_changed=True, flag_extra_arguments=True),
]


def take_assistant_tool(name: str) -> Optional[AssistantTool]:
    for tool in L_ASSISTANT_TOOLS:
        if tool.name == name:
            return tool
    return None


def assistant_tools() -> List[dict]:
    """Descriptions of the tools for a chat model (JSON schema format)."""
    return [tool.to_dict() for tool in L_ASSISTANT_TOOLS]


def assistant_tools_as_text() -> str:
    """Compact description of the tools for the prompt of a small model."""
    return "\n".join(tool.to_short_text() for tool in L_ASSISTANT_TOOLS)


def convert_argument(value, d_param: dict):
    """Value of the argument of the type given by its JSON schema; raises
    ValueError or TypeError when it does not fit.

    Keys of d_param: "type" (a list with "null" accepts None, also given as
    'none', 'null', 'identity' or '-'), "items", "minItems", "maxItems"
    (arrays), "minimum", "maximum" (numbers), "enum" and "aliases" ({item
    of enum: [words]}: the item whose word is in the answer)."""
    s_type = d_param.get("type", "string")
    if isinstance(s_type, list):
        if "null" in s_type and (value is None or (
                isinstance(value, str) and value.strip().lower() in NULL_WORDS)):
            return None
        l_type = [item for item in s_type if item != "null"]
        d_param = dict(d_param, type=l_type[0] if l_type else "string")
        return convert_argument(value, d_param)
    if s_type == "array":
        if isinstance(value, str):
            value = [word for word in re.split(r"[\s,;]+", value.strip(" []()"))
                     if word != ""]
        l_value = [convert_argument(item, d_param.get("items", {}))
                   for item in value]
        n_min = d_param.get("minItems", 0)
        n_max = d_param.get("maxItems", len(l_value))
        if not n_min <= len(l_value) <= n_max:
            s_number = (f"{n_min:}" if n_min == n_max
                        else f"from {n_min:} to {n_max:}")
            raise ValueError(f"{s_number:} values are expected, "
                             f"{len(l_value):} are given.")
        return l_value
    if s_type in ("number", "integer"):
        s_kind = "a number" if s_type == "number" else "an integer"
        if isinstance(value, str):  # '2.36 A': the only number with its unit
            l_number = RE_NUMBER.findall(value)
            if len(l_number) == 1:
                value = l_number[0]
        try:
            out = float(value) if s_type == "number" else int(value)
        except ValueError:
            raise ValueError(f"'{value:}' is not {s_kind:}.") from None
        if "minimum" in d_param and out < d_param["minimum"]:
            raise ValueError(f"{out:g} is less than {d_param['minimum']:g}.")
        if "maximum" in d_param and out > d_param["maximum"]:
            raise ValueError(f"{out:g} is more than {d_param['maximum']:g}.")
    elif s_type == "boolean":
        if isinstance(value, str):
            return value.strip().lower() in ("true", "yes", "1")
        return bool(value)
    else:
        out = str(value)
    if "enum" in d_param:
        d_enum = {str(item).lower(): item for item in d_param["enum"]}
        if str(out).lower() not in d_enum:
            s_lower = str(out).lower()
            for item, l_word in d_param.get("aliases", {}).items():
                if any(re.search(rf"(?<!\w){re.escape(str(word).lower()):}(?!\w)", s_lower)
                       for word in l_word):
                    return item
            raise ValueError(f"'{out:}' is not one of: "
                             f"{', '.join(str(item) for item in d_param['enum']):}.")
        out = d_enum[str(out).lower()]
    return out


def convert_assistant_answer(d_missing: dict, text: str):
    """Value of the answer of the user to a question of 'missing' (see
    call_assistant_tool). Raises ValueError with a message for the user
    when the answer does not fit."""
    try:
        return convert_argument(text.strip(), d_missing)
    except (TypeError, ValueError) as error:
        message = str(error)
        if d_missing.get("hint"):
            message = f"{message:} Expected: {d_missing['hint']:}."
        raise ValueError(message) from error


def form_missing(name: str, d_param: dict) -> dict:
    """The dict of 'missing' for the argument of a tool."""
    d_missing = {"name": name,
                 "question": d_param.get("question", f"the {name:}")}
    d_missing.update({key: d_param[key] for key in MISSING_KEYS
                      if key in d_param})
    return d_missing


def form_missing_text(name: str, l_missing: List[dict],
                      message: str = "") -> str:
    s_question = "; ".join(d_missing["question"] for d_missing in l_missing)
    s_text = f"{name:} needs to know: {s_question:}."
    return f"{message:}\n{s_text:}" if message else s_text


def call_assistant_tool(cryspy_object, name: str,
                        arguments: Optional[dict] = None) -> dict:
    """Run a tool by name with arguments from a chat model.

    Never raises: errors are returned as text, so that the model can read
    them and correct the call.

    Returns dict with keys 'text', 'flag_object_changed', 'flag_error'
    (True if the tool did not do what was asked: the text starts with
    "Error"; a chat program must not report such a call as done) and
    'missing' (the arguments to ask the user, see the module docstring;
    an empty list when nothing is missing), 'metrics' (numbers of the
    result: "chi_sq", "n_point", "chi_sq_per_n", "chi_sq_previous",
    "chi_sq_change", "n_refined", "max_correlation", "warnings"; see
    history.form_metrics; empty for tools without numbers) and, when the
    tool raised an exception, 'error_type' (its class name).
    """
    result = run_assistant_tool(cryspy_object, name, arguments)
    result["flag_error"] = str(result["text"]).startswith("Error")
    result.setdefault("missing", [])
    result.setdefault("metrics", {})
    return result


def run_assistant_tool(cryspy_object, name: str,
                       arguments: Optional[dict] = None) -> dict:
    tool = take_assistant_tool(name)
    if tool is None:
        return {"text": f"Error: unknown tool '{name:}'. Known tools: "
                        + ", ".join(tool.name for tool in L_ASSISTANT_TOOLS),
                "flag_object_changed": False}
    if tool.flag_object_needed and cryspy_object is None:
        return {"text": "Error: no CrysPy object is loaded.",
                "flag_object_changed": False}
    if not isinstance(arguments, dict):
        arguments = {}
    kwargs = {}
    try:
        for key, value in arguments.items():
            if key in tool.parameters:
                kwargs[key] = convert_argument(value, tool.parameters[key])
            elif tool.flag_extra_arguments:
                kwargs[key] = value
    except (TypeError, ValueError) as error:
        return {"text": f"Error in arguments of {name:}: {error:}",
                "flag_object_changed": False}
    l_missing = [key for key in tool.required if key not in kwargs]
    l_ask = [key for key in l_missing if "question" in tool.parameters[key]]
    if len(l_missing) > len(l_ask):
        return {"text": f"Error: {name:} needs arguments: {', '.join(l_missing):}",
                "flag_object_changed": False}
    if len(l_ask) > 0:
        l_missing = [form_missing(key, tool.parameters[key]) for key in l_ask]
        return {"text": form_missing_text(name, l_missing),
                "flag_object_changed": False, "missing": l_missing}
    try:
        text = tool.function(cryspy_object, **kwargs)
    except ArgumentsNeeded as error:
        return {"text": form_missing_text(name, error.l_missing, error.message),
                "flag_object_changed": False, "missing": error.l_missing}
    except Exception as error:
        return {"text": f"Error in {name:}: {type(error).__name__:}: {error:}",
                "flag_object_changed": tool.flag_object_changed,
                "error_type": type(error).__name__}
    if isinstance(text, dict):  # a tool with numbers: {"text", "metrics"}
        return {"text": text["text"], "flag_object_changed": tool.flag_object_changed,
                "metrics": text.get("metrics", {})}
    return {"text": text, "flag_object_changed": tool.flag_object_changed}
