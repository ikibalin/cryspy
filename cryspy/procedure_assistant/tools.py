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
"question"; a tool asks for other values by raising ArgumentsNeeded. When
the user calls a tool directly (e.g. '@set_parameter' typed in a chat)
call_assistant_tool(..., ask_required=True) asks every required argument
that is not given; a chat model gets an error instead and corrects its call.
"""
import re
from typing import Callable, Dict, List, Optional

import numpy

from cryspy.B_parent_classes.cl_3_data import DataN
from cryspy.B_parent_classes.cl_4_global import GlobalN
from cryspy.C_item_loop_classes.cl_1_atom_site import AtomSite
from cryspy.C_item_loop_classes.cl_1_inversed_hessian import inversed_hessian_to_correlation
from cryspy.E_data_classes.cl_1_crystal import Crystal

from cryspy.procedure_assistant.attributes import (
    describe_item_class, explain_allowed_value, find_attribute_records, find_item_classes,
    take_attribute_record, take_block, take_object_item_classes)
from cryspy.procedure_assistant.help import get_help
from cryspy.procedure_assistant.procedures import (
    describe_procedure, describe_procedures, find_cif_objects, take_procedure_parameters,
    take_procedures)
from cryspy.procedure_assistant.history import (
    form_metrics, format_metrics, keep_checkpoint, keep_snapshot, restore_snapshot, take_history)
from cryspy.procedure_assistant.recipes import recipes_as_text, suggest_recipes, take_recipe
from cryspy.procedure_assistant import experiments
from cryspy.procedure_assistant.constraints import (
    add_constraint, take_constraints, take_marked_parameters)
from cryspy.procedure_assistant.parameters import (
    explain_value_words, find_parameter_records, suggest_parameter_names, take_parameter_record,
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
    names = result.get("parameter_path", result.get("parameter_name", []))
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
    s_explanation = explain_value_words(pattern)
    if s_explanation == "":
        # '_' stays in a word: 'a' of 'length_a' is not the value 'A'
        l_word = [word for word in re.split(r"[^A-Za-z0-9_+\-]+", pattern) if word != ""]
        s_explanation = " ".join(dict.fromkeys(
            explain_allowed_value(word) for word in l_word if explain_allowed_value(word)))
    s_explanation = f" {s_explanation:}" if s_explanation else ""
    return (f"Error: no parameters match '{pattern:}'.{s_explanation:} Parameters of this "
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


def parse_mark_paths(parameters: str) -> Dict[str, str]:
    """'x=path_1, y=path_2' (or 'x: path_1; y: path_2') -> {'x': path_1, ...}."""
    d_mark_path = {}
    for part in re.split(r"\s*[,;\n]\s*", parameters.strip()):
        if part == "":
            continue
        l_mark_path = re.split(r"\s*[=:]\s*", part, maxsplit=1)
        if len(l_mark_path) != 2:
            raise ValueError(f"'{part:}' is not 'mark=path'")
        d_mark_path[l_mark_path[0]] = l_mark_path[1]
    return d_mark_path


def tool_list_constraints(cryspy_object) -> str:
    l_expression = take_constraints(cryspy_object)
    l_marked = [d_par for d_par in take_marked_parameters(cryspy_object) if d_par["mark"] != ""]
    ls_out = ["Constraints:"] + [f"  {expression:}" for expression in l_expression]
    if len(l_expression) == 0:
        ls_out = ["No constraints."]
    if len(l_marked) > 0:
        ls_out.append("Marked parameters:")
        ls_out.extend(f"  {d_par['mark']:} = {d_par['path']:} ({format_value(d_par['value']):})"
                      for d_par in l_marked)
    return "\n".join(ls_out)


def tool_add_constraint(cryspy_object, expression: str, parameters: str = "",
                        sigma: Optional[float] = None) -> str:
    try:
        d_mark_path = parse_mark_paths(parameters)
        item = add_constraint(cryspy_object, expression, d_mark_path, sigma=sigma)
    except ValueError as error:
        return f"Error: {error:}"
    return f"Constraint {item.id:} added: {item.expression:}\n" + tool_list_constraints(cryspy_object)


def tool_calculate(cryspy_object, flag_model_change: bool = False) -> dict:
    """flag_model_change: the model was changed on purpose before the
    calculation (set_chi_type): a worse chi_sq/n is expected, the next
    refinement is compared with this calculation."""
    from cryspy.procedure_rhochi.rhochi import rhochi_no_refinement
    history = take_history(cryspy_object)
    result = rhochi_no_refinement(cryspy_object)
    d_metrics = form_metrics(result, cryspy_object, d_previous=history.d_metrics,
                             flag_model_change=flag_model_change)
    keep_checkpoint(cryspy_object, d_metrics)
    text = form_calculation_text("Calculation finished.", result, cryspy_object, d_metrics)
    if flag_model_change:
        text += ("\nThe model is changed, not refined: the change of chi_sq/n comes from the change; "
                 "do not undo it. Refine the parameters of the new model next.")
    return {"text": text, "metrics": d_metrics}


def form_calculation_text(title: str, result: dict, cryspy_object, d_metrics: dict) -> str:
    # the change of chi_sq/n, R_F2 and the warnings first: a long table of parameters may be cut
    ls_out = [title, format_metrics(d_metrics), format_rhochi_result(result, cryspy_object)]
    return "\n".join(line for line in ls_out if line != "")


def take_numerical_correlation(cryspy_object, result: dict) -> dict:
    """The result with the correlations of the inversed Hessian estimated by
    numerical derivatives at the minimum (the sigmas of the refined values
    are from it too); the matrix of BFGS is accumulated along the path of
    the minimizer and may show |r| = 1 for independent parameters."""
    if not cryspy_object.is_attribute("inversed_hessian"):
        return result
    matrix = getattr(cryspy_object.inversed_hessian, "matrix", None)
    if matrix is None or len(matrix) != len(result.get("parameter_name", [])):
        return result
    matrix = numpy.array(matrix, dtype=float)
    # the matrix and the result may list the parameters in other orders
    l_label = [str(label) for label in getattr(cryspy_object.inversed_hessian, "label", [])]
    l_path = list(result.get("parameter_path", []))
    if len(l_path) == len(l_label) and sorted(l_path) == sorted(l_label):
        l_index = [l_label.index(path) for path in l_path]
        matrix = matrix[l_index, :][:, l_index]
    correlation_matrix, _ = inversed_hessian_to_correlation(matrix)
    return dict(result, correlation_matrix=correlation_matrix)


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
    result = take_numerical_correlation(cryspy_object, result)
    d_metrics = form_metrics(result, cryspy_object, n_point=d_last.get("n_point"),
                             d_previous=history.d_metrics)
    if "n_point" in d_metrics:
        result = dict(result, n_point=d_metrics["n_point"])
    keep_checkpoint(cryspy_object, d_metrics)
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
    keep_checkpoint(cryspy_object, d_metrics)
    n_refined = sum(1 for _, flag in snapshot["flags"] if flag)
    text = (f"The parameters and their refinement flags are returned to the state before the "
            f"last refinement ({n_refined:} marked for refinement).\n"
            + format_rhochi_result(result, cryspy_object))
    s_metrics = format_metrics(d_metrics, flag_undo_hint=False)
    if s_metrics != "":
        text = f"{text:}\n{s_metrics:}"
    return {"text": text, "metrics": d_metrics}


def tool_list_recipes(cryspy_object) -> str:
    return ("Recipes (typical algorithms of analysis; get_recipe gives the steps):\n"
            + recipes_as_text())


def tool_get_recipe(cryspy_object, name: str) -> str:
    recipe = take_recipe(name)
    if recipe is None:
        return f"Error: recipe '{name:}' is not known. Known recipes:\n{recipes_as_text():}"
    return recipe.to_text()


def tool_suggest_recipes(cryspy_object) -> str:
    d_kind = experiments.take_object_kinds(cryspy_object)
    if len(d_kind) == 0:
        return "Error: the object has no experiments."
    ls_out, d_name = [], {}
    for name, l_kind in d_kind.items():
        l_recipe = suggest_recipes(l_kind)
        s_recipe = ", ".join(recipe.name for recipe in l_recipe) if l_recipe else "no recipe"
        ls_out.append(f"{name:} ({', '.join(l_kind):}): {s_recipe:}")
        if l_recipe:
            d_name.setdefault(l_recipe[0].name, []).append(name)
    if d_name:
        ls_out.append("Run: " + "; ".join(f"{recipe:} for {', '.join(l_name):}"
                                          for recipe, l_name in d_name.items())
                      + " (get_recipe gives the steps).")
    return "\n".join(ls_out)


def tool_check_flip_ratio_model(cryspy_object, experiment: str = "",
                                field: Optional[float] = None,
                                polarization: Optional[float] = None,
                                magnetic_atoms: Optional[List[str]] = None,
                                lande: Optional[float] = None) -> dict:
    l_experiment = experiments.take_flip_ratio_experiments(cryspy_object, experiment)
    if len(l_experiment) == 0:
        return {"text": "Error: no experiment with flip ratios (diffrn_refln_fr)"
                        + (f" named '{experiment:}'." if experiment.strip() else ".")}
    ls_change = []
    for diffrn in l_experiment:
        if field is not None and experiments.take_field(diffrn) == 0.:
            experiments.set_field(diffrn, field)
            ls_change.append(f"{diffrn.data_name:}: magnetic field = {field:} T")
        if polarization is not None and experiments.take_polarization(diffrn) == 0.:
            experiments.set_polarization(diffrn, polarization)
            ls_change.append(f"{diffrn.data_name:}: polarization = {polarization:}")
        for crystal in experiments.take_experiment_crystals(cryspy_object, diffrn):
            if magnetic_atoms and len(experiments.take_susceptible_atoms(crystal)) == 0:
                l_label = [str(item.label) for item in crystal.atom_site.items]
                for label in magnetic_atoms:
                    if label in l_label:
                        ls_change.append(experiments.add_susceptibility(
                            crystal, label, 2. if lande is None else lande))
    l_missing = []
    for diffrn in l_experiment:
        if experiments.take_field(diffrn) == 0.:
            l_missing.append({
                "name": "field", "question": f"the magnetic field of {diffrn.data_name:} in tesla",
                "type": "number", "hint": "a number in T, e.g. 1.0"})
        if experiments.take_polarization(diffrn) == 0.:
            l_missing.append({
                "name": "polarization", "question": f"the polarization of the beam of {diffrn.data_name:}",
                "type": "number", "minimum": -1., "maximum": 1., "hint": "a number from -1 to 1, e.g. 0.9"})
        for crystal in experiments.take_experiment_crystals(cryspy_object, diffrn):
            if len(experiments.take_susceptible_atoms(crystal)) == 0:
                l_label = [str(item.label) for item in crystal.atom_site.items]
                l_missing.append({
                    "name": "magnetic_atoms",
                    "question": f"which atoms of {crystal.data_name:} are magnetic",
                    "type": "array", "items": {"type": "string", "enum": l_label}, "minItems": 1,
                    "hint": "labels of atoms: " + ", ".join(l_label)})
                l_missing.append({
                    "name": "lande", "question": "the Lande factor of the magnetic atoms",
                    "type": "number", "proposal": 2., "s_proposal": "2.0 (spin only)",
                    "hint": "a number, e.g. 2.0 for spin only"})
    l_missing = list({d_missing["name"]: d_missing for d_missing in l_missing}.values())
    if len(l_missing) > 0:
        l_problem = []
        for diffrn in l_experiment:
            l_problem.extend(experiments.find_flip_ratio_problems(cryspy_object, diffrn))
        raise ArgumentsNeeded(l_missing, "The model gives no magnetic scattering: "
                                         + "; ".join(l_problem) + ".")
    # new tensors and zero ones (FR = 1, no derivative step) get the best start value
    l_item, l_label = [], []
    for diffrn in l_experiment:
        for crystal in experiments.take_experiment_crystals(cryspy_object, diffrn):
            for item in crystal.atom_site_susceptibility.items:
                if (str(item.label) in (magnetic_atoms or []) and item.chi_type == "Ciso") or all(
                        float(getattr(item, f"chi_{ij:}", 0.)) == 0. for ij in ("11", "22", "33")):
                    if item.chi_type != "Ciso":
                        item.chi_type = "Ciso"
                    if item not in l_item:
                        l_item.append(item)
                        l_label.append(f"{crystal.data_name:}: {item.label:}")
    if len(l_item) > 0:
        chi_start = experiments.scan_chi_start(cryspy_object, l_item)
        ls_change.append(f"start value chi_11 = {chi_start:g} mu_B/T (Ciso), the best of "
                         f"{', '.join(f'{chi:g}' for chi in experiments.L_CHI_START):} for "
                         + ", ".join(l_label))
    result = tool_calculate(cryspy_object)
    ls_out = ["Changed:\n" + "\n".join(ls_change)] if ls_change else []
    flag_magnetic = True
    for diffrn in l_experiment:
        deviation = experiments.calc_flip_ratio_deviation(diffrn)
        if deviation is None or deviation < experiments.FLIP_RATIO_TOLERANCE:
            flag_magnetic = False
            ls_out.append(f"{diffrn.data_name:}: all calculated flip ratios are equal to 1 although "
                          "the field, the polarization and the magnetic atoms are given: check the "
                          "magnetic atoms and the reflections.")
        else:
            ls_out.append(f"{diffrn.data_name:}: the calculated flip ratios differ from 1 "
                          f"(max |FR_calc - 1| = {deviation:.3g}): the model gives magnetic scattering.")
    if not flag_magnetic:
        ls_out.insert(0, "Error: the model gives no magnetic scattering.")
    ls_out.append(result["text"])
    return {"text": "\n".join(ls_out), "metrics": result["metrics"]}


def tool_set_chi_type(cryspy_object, chi_type: str, atoms: str = "") -> dict:
    l_changed = experiments.set_chi_type(cryspy_object, chi_type, atoms)
    if len(l_changed) == 0:
        if not any(crystal.is_attribute("atom_site_susceptibility")
                   for crystal in experiments.take_crystals(cryspy_object)):
            return {"text": "Error: no atom has a susceptibility tensor (atom_site_susceptibility)."}
        return {"text": f"The susceptibility tensors are {chi_type:} already."}
    # the model is changed: the next refinement is compared with this calculation
    result = tool_calculate(cryspy_object, flag_model_change=True)
    return {"text": f"Susceptibility type {chi_type:} is set for:\n" + "\n".join(l_changed)
                    + "\n" + result["text"],
            "metrics": result["metrics"]}


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


def tool_list_attributes(cryspy_object, pattern: str = "") -> str:
    l_record = find_attribute_records(cryspy_object, pattern)
    if len(l_record) == 0:
        s_explanation = explain_allowed_value(pattern)
        if s_explanation:
            return f"Error: no attributes match '{pattern:}'. {s_explanation:}"
        return (f"Error: no attributes match '{pattern:}'. Use summarize_object to see "
                "the data blocks and describe_category for the attributes of a category.")
    ls_out = ["path | value | allowed values"]
    for record in l_record[:MAX_PARAMETERS_IN_LIST]:
        allowed = record.allowed_values
        s_allowed = ""
        if allowed is not None and len(allowed) <= 30:
            s_allowed = ", ".join(str(_) for _ in allowed)
        elif allowed is not None:
            s_allowed = f"{len(allowed):} values"
        ls_out.append(f"{record.path:} | {format_value(record.value):} | {s_allowed:}")
    if len(l_record) > MAX_PARAMETERS_IN_LIST:
        ls_out.append(f"... {len(l_record) - MAX_PARAMETERS_IN_LIST:} more, "
                      "give a narrower pattern.")
    return "\n".join(ls_out)


def tool_set_attribute(cryspy_object, path: str, value: str) -> str:
    record = take_attribute_record(cryspy_object, path)
    if record is None:
        l_record = find_attribute_records(cryspy_object, path, flag_defined_only=False)
        if len(l_record) == 1:
            record = l_record[0]
        elif len(l_record) > 1:
            return ("Error: the path is ambiguous, give one of:\n"
                    + "\n".join(record.path for record in l_record[:20]))
        else:
            return (f"Error: no attribute '{path:}'. Use list_attributes to find the path.")
    if record.flag_constraint:
        return f"Error: {record.path:} is constrained by symmetry and cannot be set."
    old_value = record.value
    try:
        record.set_value(value)
    except ValueError as error:
        return f"Error: {error:}"
    return f"{record.path:}: {format_value(old_value):} -> {format_value(record.value):}"


def tool_show_item(cryspy_object, path: str, max_length: int = 8000) -> str:
    block = take_block(cryspy_object, path)
    s_value = ""
    if block is None:  # the path of one attribute: its value and its item
        record = take_attribute_record(cryspy_object, path)
        if record is None:
            l_record = find_attribute_records(cryspy_object, path)
            record = l_record[0] if len(l_record) == 1 else None
        if record is None:
            l_name = [str(data.get_name()) for data in getattr(cryspy_object, "items", [])]
            return (f"Error: '{path:}' is not found. Give a data block, a category, an item of a "
                    "loop or an attribute: 'crystal_fe', 'crystal_fe.cell', "
                    "'crystal_fe.atom_site[Fe]', 'diffrn_exp1.setup.field'. The object has: "
                    f"{', '.join(l_name):}. list_attributes finds the paths of attributes.")
        block = record.item
        s_value = f"{record.path:} = {format_value(record.value):}\n\n"
    text = s_value + block.to_cif()
    if len(text) > max_length:
        text = text[:max_length] + "\n[... the rest is cut ...]"
    return text


def tool_describe_category(cryspy_object, name: str = "") -> str:
    if name.strip() == "":
        l_class = take_object_item_classes(cryspy_object) if cryspy_object is not None else []
        if len(l_class) == 0:
            return "Error: give the name of a category, e.g. 'cell' or 'atom_site'."
    else:
        l_class = find_item_classes(name)
        if len(l_class) == 0:
            return f"Error: there is no category '{name:}'."
        if len(l_class) > 10:
            return ("Several categories match, give one of: "
                    + ", ".join(item_class.PREFIX for item_class in l_class))
    return "\n\n".join(describe_item_class(item_class) for item_class in l_class)


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


def take_named_crystals(cryspy_object, crystal: str = "") -> list:
    """Crystals of the object, or the one named; ValueError if not found."""
    if crystal.strip() == "":
        l_crystal = [item for item in cryspy_object.items if isinstance(item, Crystal)]
        if len(l_crystal) == 0:
            raise ValueError("the object has no crystal.")
        return l_crystal
    item = take_data_block(cryspy_object, crystal)
    if not isinstance(item, Crystal):
        raise ValueError(f"crystal '{crystal:}' is not found.")
    return [item]


def tool_list_atoms_in_unit_cell(cryspy_object, crystal: str = "",
                                 max_atoms: int = 500) -> str:
    try:
        l_crystal = take_named_crystals(cryspy_object, crystal)
    except ValueError as error:
        return f"Error: {error:}"
    ls_out = []
    for item in l_crystal:
        atoms = item.calc_atoms_in_unit_cell()
        moment = atoms["moment_crystalaxis"]
        n_atom = len(atoms["label"])
        ls_out.append(f"Crystal {item.data_name:}: {n_atom:} atoms in the unit cell "
                      "(generated by the symmetry elements from atom_site).")
        s_head = "| atom | from | x | y | z |"
        s_line = "|---|---|---|---|---|"
        if moment is not None:
            s_head += " M_a, M_b, M_c (μB, along a/|a|, b/|b|, c/|c|) |"
            s_line += "---|"
        ls_out += ["", s_head, s_line]
        for i_out in range(min(n_atom, max_atoms)):
            xyz = atoms["fract_xyz"][:, i_out]
            line = (f"| {atoms['label'][i_out]:} | "
                    f"{item.atom_site.items[atoms['atom_index'][i_out]].label:} | "
                    f"{xyz[0]:.5f} | {xyz[1]:.5f} | {xyz[2]:.5f} |")
            if moment is not None:
                m = moment[:, i_out]
                line += (" |" if numpy.any(numpy.isnan(m)) else
                         " " + ", ".join(f"{v:.4f}" for v in m) + " |")
            ls_out.append(line)
        if n_atom > max_atoms:
            ls_out.append(f"[... {n_atom - max_atoms:} more atoms are cut ...]")
        ls_out.append("")
    return "\n".join(ls_out)


def tool_convert_to_p1(cryspy_object, crystal: str = "", name: str = "") -> str:
    try:
        l_crystal = take_named_crystals(cryspy_object, crystal)
    except ValueError as error:
        return f"Error: {error:}"
    if len(l_crystal) > 1:
        return ("Error: the object has several crystals ("
                + ", ".join(item.data_name for item in l_crystal)
                + "); give the crystal to convert.")
    item = l_crystal[0]
    data_name = name.strip() if name.strip() != "" else f"{item.data_name}_p1"
    if take_data_block(cryspy_object, data_name) is not None:
        return f"Error: the data block '{data_name:}' exists already; give another name."
    crystal_p1 = item.to_p1(data_name)
    cryspy_object.add_items([crystal_p1])
    ls_out = [f"Crystal {data_name:} in P1 is added: {len(crystal_p1.atom_site.items):} atoms "
              f"(from {len(item.atom_site.items):} atoms of {item.data_name:})."]
    for loop_name in ("atom_site_aniso", "atom_site_moment", "atom_site_susceptibility",
                      "atom_site_scat"):
        if crystal_p1.is_attribute(loop_name):
            ls_out.append(f"{loop_name:}: {len(getattr(crystal_p1, loop_name).items):} atoms.")
    ls_out.append("Multipoles and local axes are not converted. Experiments still use "
                  f"{item.data_name:}; the new crystal is not in their phases.")
    return "\n".join(ls_out)


def tool_report_ellipsoids(cryspy_object, kind: str = "ADP", crystal: str = "",
                           unit_cell: bool = False, field: float = 1.) -> str:
    try:
        l_crystal = take_named_crystals(cryspy_object, crystal)
        ls_out = []
        for item in l_crystal:
            ls_out.append(f"Crystal {item.data_name:}.")
            ls_out.append(item.report_ellipsoids(kind, flag_unit_cell=unit_cell, field=field))
            ls_out.append("")
    except ValueError as error:
        return f"Error: {error:}"
    return "\n".join(ls_out)


def tool_add_crystal_from_cod(cryspy_object, formula: str = "",
                              cod_id: str = "", name: str = "") -> str:
    from cryspy.H_functions_global.cod import add_crystal_from_cod
    if formula.strip() == "" and str(cod_id).strip() == "":
        raise ArgumentsNeeded([{
            "name": "formula", "type": "string",
            "question": "the chemical formula of the crystal or its COD number",
            "hint": "e.g. TiO2 or 9008460"}])
    if formula.strip().isdigit() and str(cod_id).strip() == "":  # the COD number as the answer
        formula, cod_id = "", formula.strip()
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
          "'scale factor', 'polarization', 'unit cell', 'adp'. "
          "Write the path only, without the value or '()' of rcif.")

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
        "Mark the parameters matching the pattern for refinement (as '()' "
        "after a value in rcif); to try one parameter give its path from "
        "list_parameters, e.g. 'crystal_fe.cell.length_a'.",
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
        "list_attributes",
        "List all attributes of the object (also text ones such as "
        "atom_site_susceptibility.chi_type or extinction.model) with their values "
        "and allowed values.",
        tool_list_attributes,
        {"pattern": {"type": "string", "description": (
            "Path or pattern: 'chi_type', '*.extinction.model', "
            "'atom_site[Fe]'; empty for all.")}}),
    AssistantTool(
        "set_attribute",
        "Set the value of one attribute of any category (text, integer or number), "
        "e.g. 'exp.extinction.model' to 'gauss'; text values must be allowed values.",
        tool_set_attribute,
        {"path": {"type": "string", "description": "Path from list_attributes."},
         "value": {"type": "string", "description": "New value."}},
        required=("path", "value"), flag_object_changed=True),
    AssistantTool(
        "list_constraints",
        "List the constraints between refined parameters (x+y=1, hard or soft) and "
        "the marks of the parameters used in them.",
        tool_list_constraints),
    AssistantTool(
        "add_constraint",
        "Add a constraint between refined parameters, e.g. 'x+y=1' "
        "to keep the sum of two occupancies; the parameters are named by "
        "marks, given in 'parameters' as 'x=path_1, y=path_2' (they are "
        "refined then). Without sigma the equation is linear and holds "
        "exactly; with sigma it is soft (restraint, any formula of the marks).",
        tool_add_constraint,
        {"expression": {"type": "string",
                        "description": "Linear equation of the marks: 'x+y=1', '2*a-b=0'."},
         "parameters": {"type": "string",
                        "description": ("Marks of the parameters: 'x=crystal_fe.atom_site[Fe1].occupancy, "
                                        "y=crystal_fe.atom_site[Co1].occupancy'; empty if "
                                        "the parameters are marked already.")},
         "sigma": {"type": "number",
                   "description": "Allowed deviation of a soft constraint; omit for a hard one."}},
        required=("expression", ), flag_object_changed=True),
    AssistantTool(
        "show_item",
        "Show in rcif format one data block, category (item or loop) or one item of a loop; "
        "for the path of an attribute its value and its item.",
        tool_show_item,
        {"path": {"type": "string", "description": (
            "'crystal_fe', 'crystal_fe.cell', 'crystal_fe.atom_site[Fe]' or an attribute "
            "'diffrn_exp1.setup.field'; the data block may be given by its data name ('exp1').")}},
        required=("path", )),
    AssistantTool(
        "describe_category",
        "Attributes of a category (any CrysPy item or loop class) with types, "
        "refinable ones and allowed values; without a name the categories of the object.",
        tool_describe_category,
        {"name": {"type": "string", "description": (
            "Category, e.g. 'cell', 'atom_site_susceptibility', 'extinction'.")}},
        flag_object_needed=False),
    AssistantTool(
        "calculate",
        "Calculate the model without refinement: chi-square and calculated data.",
        tool_calculate, flag_object_changed=True),
    AssistantTool(
        "refine",
        "Refine the parameters marked for refinement (RhoChi); returns "
        "chi-square, the change of chi_sq/n (and of R_F2 for single crystals) "
        "since the previous calculation, refined values, strong correlations "
        "and warnings (unphysical values).",
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
        "suggest_recipes",
        "Kinds of the experiments of the object (single crystal, powder, flip ratios, "
        "polarized or not, nuclear or magnetic crystals) and the recipes for them.",
        tool_suggest_recipes),
    AssistantTool(
        "check_flip_ratio_model",
        "Calculate the model of the flip ratio experiments and check that it gives "
        "magnetic scattering (not all calculated flip ratios equal to 1); asks the user for "
        "a zero magnetic field or polarization and for the magnetic atoms when no atom has a "
        "susceptibility tensor.",
        tool_check_flip_ratio_model,
        {"experiment": {"type": "string",
                        "description": "Name of the experiment; empty for all with flip ratios."},
         "field": {"type": "number", "description": "Magnetic field, T."},
         "polarization": {"type": "number", "description": "Polarization of the beam."},
         "magnetic_atoms": {"type": "array", "items": {"type": "string"},
                            "description": "Labels of the magnetic atoms."},
         "lande": {"type": "number", "description": "Lande factor of the magnetic atoms."}},
        flag_object_changed=True),
    AssistantTool(
        "set_chi_type",
        "Set the type of the susceptibility tensors: 'Ciso' (isotropic, only chi_11 is "
        "refined) or 'Cani' (anisotropic, the components allowed by the symmetry of the site).",
        tool_set_chi_type,
        {"chi_type": {"type": "string", "enum": ["Ciso", "Cani"]},
         "atoms": {"type": "string",
                   "description": "Labels of the atoms; empty for all atoms with a susceptibility."}},
        required=("chi_type", ), flag_object_changed=True),
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
        "list_atoms_in_unit_cell",
        "All atoms of the unit cell made by the symmetry elements (the "
        "structure in P1): fractional coordinates and ordered moments of "
        "every atom; nothing is changed.",
        tool_list_atoms_in_unit_cell,
        {"crystal": {"type": "string",
                     "description": "Name of the crystal; empty for all."}}),
    AssistantTool(
        "convert_to_p1",
        "Add a copy of a crystal in space group P1: every atom of the unit "
        "cell becomes a separate atom ('Fe1_1', 'Fe1_2', ...) with rotated "
        "anisotropic displacement parameters, moments and susceptibility "
        "tensors.",
        tool_convert_to_p1,
        {"crystal": {"type": "string",
                     "description": "Name of the crystal; may be empty when there is one."},
         "name": {"type": "string",
                  "description": "Data name of the new crystal; by default the old one with '_p1'."}},
        flag_object_changed=True),
    AssistantTool(
        "report_ellipsoids",
        "Principal axes of the atomic ellipsoids: displacement ellipsoids "
        "(kind 'ADP': U_i, rms displacements) or magnetization ellipsoids "
        "of the susceptibility tensors (kind 'susceptibility': |χ_i| and the "
        "moments induced by a field), with their directions in Cartesian and "
        "crystal axes; for the atoms of atom_site or of the whole unit cell.",
        tool_report_ellipsoids,
        {"kind": {"type": "string", "enum": ["ADP", "susceptibility"]},
         "crystal": {"type": "string",
                     "description": "Name of the crystal; empty for all."},
         "unit_cell": {"type": "boolean",
                       "description": "For every atom of the unit cell (orientations "
                                      "of the equivalent atoms)."},
         "field": {"type": "number",
                   "description": "Magnetic field, T, for the induced moments; 1 by default."}}),
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
                 "question": d_param.get("question", f"the {name.replace('_', ' '):}")}
    d_missing.update({key: d_param[key] for key in MISSING_KEYS
                      if key in d_param})
    if "hint" not in d_missing and d_param.get("description"):
        d_missing["hint"] = d_param["description"]
    return d_missing


def form_missing_text(name: str, l_missing: List[dict],
                      message: str = "") -> str:
    s_question = "; ".join(d_missing["question"] for d_missing in l_missing)
    s_text = f"{name:} needs to know: {s_question:}."
    return f"{message:}\n{s_text:}" if message else s_text


def call_assistant_tool(cryspy_object, name: str,
                        arguments: Optional[dict] = None, *,
                        ask_required: bool = False) -> dict:
    """Run a tool by name with arguments from a chat model.

    ask_required: the user called the tool (not a model), so every
    required argument that is not given is asked ('missing'), not only
    those with "question".

    Never raises: errors are returned as text, so that the model can read
    them and correct the call.

    Returns dict with keys 'text', 'flag_object_changed', 'flag_error'
    (True if the tool did not do what was asked: the text starts with
    "Error"; a chat program must not report such a call as done) and
    'missing' (the arguments to ask the user, see the module docstring;
    an empty list when nothing is missing), 'metrics' (numbers of the
    result: "chi_sq", "n_point", "chi_sq_per_n", "chi_sq_per_n_previous",
    "chi_sq_per_n_change", "r_f2" and "r_f2_previous" ({data block: R_F2}
    of single crystals), "n_refined", "max_correlation", "warnings"; see
    history.form_metrics; empty for tools without numbers) and, when the
    tool raised an exception, 'error_type' (its class name).
    """
    result = run_assistant_tool(cryspy_object, name, arguments, ask_required=ask_required)
    result["flag_error"] = str(result["text"]).startswith("Error")
    result.setdefault("missing", [])
    result.setdefault("metrics", {})
    return result


def run_assistant_tool(cryspy_object, name: str,
                       arguments: Optional[dict] = None, ask_required: bool = False) -> dict:
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
    l_ask = [key for key in l_missing if ask_required or "question" in tool.parameters[key]]
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
