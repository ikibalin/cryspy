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
"""
import inspect
from typing import Callable, Dict, List, Optional

import numpy

from cryspy.B_parent_classes.cl_3_data import DataN
from cryspy.B_parent_classes.cl_4_global import GlobalN
from cryspy.C_item_loop_classes.cl_1_atom_site import AtomSite
from cryspy.E_data_classes.cl_1_crystal import Crystal

from cryspy.procedure_assistant.help import get_help
from cryspy.procedure_assistant.parameters import (
    find_parameter_records, take_parameter_record, take_parameter_records)

MAX_PARAMETERS_IN_LIST = 200
CORRELATION_LIMIT = 0.7


class AssistantTool:
    """A function with the description of its arguments for a chat model."""

    def __init__(self, name: str, description: str, function: Callable,
                 parameters: Optional[Dict[str, dict]] = None,
                 required: tuple = (), flag_object_changed: bool = False,
                 flag_object_needed: bool = True):
        self.name = name
        self.description = description
        self.function = function
        self.parameters = parameters if parameters is not None else {}
        self.required = required
        self.flag_object_changed = flag_object_changed
        self.flag_object_needed = flag_object_needed

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


def tool_list_parameters(cryspy_object, pattern: str = "",
                         refined_only: bool = False) -> str:
    l_record = find_parameter_records(
        cryspy_object, pattern, flag_refined_only=refined_only)
    if len(l_record) == 0:
        return f"No parameters match '{pattern:}'."
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
        return "Give a pattern of the parameters to refine."
    l_record = find_parameter_records(cryspy_object, pattern)
    if len(l_record) == 0:
        return f"No parameters match '{pattern:}'."
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
    return "\n".join(ls_out)


def tool_fix_parameters(cryspy_object, pattern: str = "") -> str:
    l_record = find_parameter_records(
        cryspy_object, pattern, flag_refined_only=True)
    for record in l_record:
        record.set_refinement(False)
    if len(l_record) == 0:
        return "No refined parameters match the pattern."
    return "Fixed:\n" + "\n".join(record.path for record in l_record)


def tool_set_parameter(cryspy_object, path: str, value: float) -> str:
    record = take_parameter_record(cryspy_object, path)
    if record is None:
        l_record = find_parameter_records(cryspy_object, path)
        if len(l_record) == 1:
            record = l_record[0]
        elif len(l_record) > 1:
            return ("The path is ambiguous, give one of:\n"
                    + "\n".join(record.path for record in l_record[:20]))
        else:
            return f"Parameter '{path:}' is not found."
    if record.flag_constraint:
        return f"{record.path:} is constrained by symmetry and cannot be set."
    old_value = record.value
    record.set_value(value)
    return (f"{record.path:}: {format_value(old_value):} -> "
            f"{format_value(record.value):}")


def tool_calculate(cryspy_object) -> str:
    from cryspy.procedure_rhochi.rhochi import rhochi_no_refinement
    result = rhochi_no_refinement(cryspy_object)
    return "Calculation finished.\n" + format_rhochi_result(result, cryspy_object)


def tool_refine(cryspy_object) -> str:
    from cryspy.procedure_rhochi.rhochi import rhochi_rietveld_refinement
    if len(cryspy_object.get_variable_names()) == 0:
        return ("No parameters are marked for refinement. "
                "Use refine_parameters first.")
    result = rhochi_rietveld_refinement(cryspy_object)
    return "Refinement finished.\n" + format_rhochi_result(result, cryspy_object)


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
            return f"Data block '{block:}' is not found."
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
        return f"Crystal '{crystal:}' is not found."
    atom_site = data.atom_site
    if label in [str(item.label) for item in atom_site.items]:
        return f"Atom '{label:}' already exists in '{crystal:}'."
    if adp_type == "":
        # the same kind of isotropic parameter as the other atoms
        l_adp_type = [item.adp_type for item in atom_site.items
                      if item.is_attribute("adp_type")]
        adp_type = "Uiso" if "Uiso" in l_adp_type else "Biso"
    if adp_type not in ("Biso", "Uiso"):
        return ("adp_type should be 'Biso' or 'Uiso'; "
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


def take_procedures() -> Dict[str, Callable]:
    """CrysPy procedures that need only the object (as in CrysPy Editor menus)."""
    import cryspy
    d_procedure = {}
    for function in cryspy.L_FUNCTION + cryspy.L_FUNCTION_ADD:
        try:
            signature = inspect.signature(function)
        except (TypeError, ValueError):
            continue
        l_param = list(signature.parameters.values())
        if len(l_param) == 0:
            continue
        if l_param[0].annotation is not GlobalN:
            continue
        if any(param.default is inspect.Parameter.empty
               for param in l_param[1:]):
            continue
        d_procedure[function.__name__] = function
    return d_procedure


def tool_list_procedures(cryspy_object) -> str:
    ls_out = []
    for name, function in sorted(take_procedures().items()):
        doc = (function.__doc__ or "").strip().split("\n")[0]
        ls_out.append(f"- {name:}: {doc:}")
    return "\n".join(ls_out)


def tool_run_procedure(cryspy_object, name: str) -> str:
    d_procedure = take_procedures()
    if name not in d_procedure:
        return (f"Procedure '{name:}' is not known. Available:\n"
                + tool_list_procedures(cryspy_object))
    result = d_procedure[name](cryspy_object)
    if isinstance(result, dict) and "chi_sq" in result:
        return f"{name:} finished.\n" + format_rhochi_result(
            result, cryspy_object)
    if result is None:
        return f"{name:} finished."
    text = str(result)
    if len(text) > 4000:
        text = text[:4000] + "\n[... the rest is cut ...]"
    return f"{name:} finished.\n{text:}"


S_PATH = ("Parameter path or pattern: 'crystal_pbso4.cell.length_a', "
          "'atom_site[Fe].fract_x', 'b_iso', '*.pd_background[*].intensity'. "
          "Without '*' any path containing the text matches.")

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
        "chi-square, refined values and strong correlations.",
        tool_refine, flag_object_changed=True),
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
        "list_procedures",
        "List CrysPy procedures that can be run on the object "
        "(MEM, simulations, reports, ...).",
        tool_list_procedures),
    AssistantTool(
        "run_procedure",
        "Run a CrysPy procedure from list_procedures by its name.",
        tool_run_procedure,
        {"name": {"type": "string"}}, required=("name", ),
        flag_object_changed=True),
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
    s_type = d_param.get("type", "string")
    if s_type == "number":
        return float(value)
    if s_type == "integer":
        return int(value)
    if s_type == "boolean":
        if isinstance(value, str):
            return value.strip().lower() in ("true", "yes", "1")
        return bool(value)
    return str(value)


def call_assistant_tool(cryspy_object, name: str,
                        arguments: Optional[dict] = None) -> dict:
    """Run a tool by name with arguments from a chat model.

    Never raises: errors are returned as text, so that the model can read
    them and correct the call.

    Returns dict with keys 'text' and 'flag_object_changed'.
    """
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
    except (TypeError, ValueError) as error:
        return {"text": f"Error in arguments of {name:}: {error:}",
                "flag_object_changed": False}
    l_missing = [key for key in tool.required if key not in kwargs]
    if len(l_missing) > 0:
        return {"text": f"Error: {name:} needs arguments: {', '.join(l_missing):}",
                "flag_object_changed": False}
    try:
        text = tool.function(cryspy_object, **kwargs)
    except Exception as error:
        return {"text": f"Error in {name:}: {type(error).__name__:}: {error:}",
                "flag_object_changed": tool.flag_object_changed}
    return {"text": text, "flag_object_changed": tool.flag_object_changed}
