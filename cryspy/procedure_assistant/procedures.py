"""Procedures of CrysPy described from their signatures and docstrings.

The procedures are the functions of cryspy.L_FUNCTION and L_FUNCTION_ADD
whose parameters are the CrysPy object (GlobalN), basic values (float,
int, bool, str) or items, loops and data blocks of the object (as Cell).
Their description is taken from the type hints and from the docstring in
the numpydoc format, so that a function needs no other description to be
shown in a menu, in a form of its parameters or by an AI assistant:

    def mempy_cycle_with_parameters(obj: GlobalN, n_cycle: int = 10):
        '''Reconstruct the density and refine susceptibilities in cycles.

        Parameters
        ----------
        obj : GlobalN
            The object with the crystal, the MEM parameters and the
            polarized neutron experiments.
        n_cycle : int, optional
            Number of cycles of density reconstruction and refinement.
        '''

describe_procedure(function) gives the title, the group of the menu, the
summary (first line), the description and, for every parameter, its kind,
type, default, description and the question to ask it from the user.
"""
import inspect
import re
from typing import Callable, Dict, List, Optional

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN
from cryspy.B_parent_classes.cl_3_data import DataN
from cryspy.B_parent_classes.cl_4_global import GlobalN

D_BASIC_TYPE = {float: "number", int: "integer", bool: "boolean", str: "string"}
# words of the titles of the menus
D_TITLE_WORD = {"Mempy": "MEMPy", "Rhochi": "RhoChi", "Calc": "Calculate"}
PARAMETER_HEADERS = ("parameters", "arguments", "args")
RETURN_HEADERS = ("returns", "return", "output", "outputs")


def take_annotation_types(annotation) -> tuple:
    """Classes of an annotation: (A, ) for A, (A, B) for Union[A, B]."""
    l_type = getattr(annotation, "__args__", None)
    if l_type is None or getattr(annotation, "__origin__", None) is None:
        return (annotation, )
    return tuple(l_type)


def take_procedure_parameters(function: Callable) -> Optional[List[tuple]]:
    """(name, kind, types, default) of the parameters of a CrysPy function,
    None if it cannot be run as a procedure. kind: "global" (the object),
    "basic" (number, integer, boolean, string), "cif" (an item, loop or
    data block of the object), "skip" (keyword with a default: callback,
    d_info)."""
    try:
        signature = inspect.signature(function)
    except (TypeError, ValueError):
        return None
    l_out = []
    for param in signature.parameters.values():
        if param.kind in (param.VAR_POSITIONAL, param.VAR_KEYWORD):
            return None
        types = take_annotation_types(param.annotation)
        if GlobalN in types:
            kind = "global"
        elif len(types) == 1 and types[0] in D_BASIC_TYPE:
            kind = "basic"
        elif all(isinstance(item, type) and issubclass(item, (ItemN, LoopN, DataN))
                 and item not in (ItemN, LoopN, DataN) for item in types):
            kind = "cif"
        elif param.default is not inspect.Parameter.empty:
            kind = "skip"
        else:
            return None
        l_out.append((param.name, kind, types, param.default))
    return l_out


def take_procedures() -> Dict[str, Callable]:
    """CrysPy procedures (of the menus of CrysPy Editor): their parameters
    are the object, numbers, strings, booleans, or items, loops and data
    blocks of the object (as Cell)."""
    import cryspy
    d_procedure = {}
    for function in cryspy.L_FUNCTION + cryspy.L_FUNCTION_ADD:
        name = getattr(function, "__name__", "")
        if name.startswith(("check_", "plot_", "get_")):
            continue
        l_param = take_procedure_parameters(function)
        if l_param is None or len(l_param) == 0:
            continue
        l_kind = [kind for _, kind, _, _ in l_param]
        if "global" not in l_kind and "cif" not in l_kind:
            continue
        d_procedure[name] = function
    return {name: d_procedure[name] for name in sorted(d_procedure)}


def parse_docstring(doc: str) -> dict:
    """Summary, description, parameters and returns of a docstring in the
    numpydoc format (sections "Parameters" and "Returns" underlined with
    "-"); a docstring without sections is the summary and description.

    Returns {"summary": str, "description": str, "returns": str,
    "parameters": {name: {"type": str, "optional": bool,
    "description": str}}}."""
    l_line = inspect.cleandoc(doc or "").splitlines()
    d_section = {"": []}
    section = ""
    i_line = 0
    while i_line < len(l_line):
        line = l_line[i_line]
        if (i_line + 1 < len(l_line) and line.strip() != ""
                and re.fullmatch(r"\s*-{3,}\s*", l_line[i_line + 1])):
            section = line.strip().rstrip(":").lower()
            d_section[section] = []
            i_line += 2
            continue
        d_section[section].append(line)
        i_line += 1

    l_text = "\n".join(d_section[""]).strip().split("\n\n", 1)
    summary = " ".join(l_text[0].split())
    description = l_text[1].strip() if len(l_text) > 1 else ""

    d_parameter = {}
    l_param_line = next((d_section[key] for key in PARAMETER_HEADERS if key in d_section), [])
    name = None
    for line in l_param_line:
        if line.strip() == "":
            continue
        if not line.startswith((" ", "\t")):
            # 'name : type' or 'Na, Nb, Nc : int' (several with one description)
            match = re.match(r"\s*-?\s*(\w+(?:\s*,\s*\w+)*)\s*(?::\s*(.*))?$", line)
            if match is None:
                continue
            l_name = [word.strip() for word in match.group(1).split(",")]
            s_type = (match.group(2) or "").strip()
            flag_optional = bool(re.search(r"\boptional\b", s_type))
            s_type = re.sub(r",?\s*optional\b", "", s_type).strip(" ,")
            d_description = {"type": s_type, "optional": flag_optional, "description": ""}
            for name in l_name:
                d_parameter[name] = d_description
        elif name is not None:
            text = d_parameter[name]["description"]
            d_parameter[name]["description"] = f"{text:} {line.strip():}".strip()
    returns = "\n".join(next((d_section[key] for key in RETURN_HEADERS if key in d_section), [])).strip()
    return {"summary": summary, "description": description,
            "parameters": d_parameter, "returns": returns}


def form_title(name: str) -> str:
    """Title of a procedure in menus: 'rhochi_no_refinement' -> 'RhoChi No
    Refinement'."""
    return " ".join(D_TITLE_WORD.get(word, word) for word in name.replace("_", " ").title().split())


def form_question(description: str, default_text: str) -> str:
    """What is asked from the user: the first sentence (up to '.', ';' or
    ':') of the description of the parameter ('the starting value of
    lambda')."""
    sentence = re.split(r"(?<=[^.]\.)\s|[;:]\s", description.strip(), maxsplit=1)[0].rstrip(".").strip()
    if sentence == "":
        return default_text
    if not (len(sentence) > 1 and sentence[1].isupper()):  # keep 'MEM', 'UB'
        sentence = sentence[0].lower() + sentence[1:]
    if sentence.split()[0] not in ("the", "a", "an"):
        sentence = f"the {sentence:}"
    return sentence


def describe_procedure(function: Callable) -> Optional[dict]:
    """Description of a procedure from its signature and its docstring;
    None if it cannot be run as a procedure.

    Returns {"name", "title", "group" (first word of the title),
    "summary", "description", "returns", "parameters": [{"name", "kind"
    ("global", "basic", "cif"), "type" (JSON schema type for "basic"),
    "classes" (class names for "cif"), "default" (when there is one),
    "description", "question"}, ...]}."""
    l_param = take_procedure_parameters(function)
    if l_param is None:
        return None
    name = function.__name__
    d_doc = parse_docstring(function.__doc__)
    title = form_title(name)
    l_out = []
    for param, kind, types, default in l_param:
        if kind == "skip":
            continue
        d_param_doc = d_doc["parameters"].get(param, {})
        description = d_param_doc.get("description", "")
        d_param = {"name": param, "kind": kind, "description": description}
        if kind == "basic":
            d_param["type"] = D_BASIC_TYPE[types[0]]
            d_param["question"] = form_question(description, f"the {param.replace('_', ' '):}")
        elif kind == "cif":
            d_param["classes"] = [item.__name__ for item in types]
            d_param["question"] = form_question(
                description, f"the {' or '.join(d_param['classes']):} for {param:}")
        if default is not inspect.Parameter.empty:
            d_param["default"] = default
        l_out.append(d_param)
    return {"name": name, "title": title, "group": title.split(" ")[0],
            "summary": d_doc["summary"], "description": d_doc["description"],
            "returns": d_doc["returns"], "parameters": l_out}


def describe_procedures() -> List[dict]:
    """Descriptions of all procedures (describe_procedure), e.g. to build
    the menus of a program."""
    return [describe_procedure(function) for function in take_procedures().values()]


def find_cif_objects(cryspy_object, types: tuple) -> Dict[str, object]:
    """{label: object} of the items, loops and data blocks of the object
    of the given classes; labels as 'crystal_name.cell'."""
    d_found = {}
    if cryspy_object is None:
        return d_found
    for item in cryspy_object.items:
        if isinstance(item, DataN):
            if isinstance(item, types):
                d_found[item.data_name] = item
            for sub_item in item.items:
                if isinstance(sub_item, types):
                    d_found[f"{item.data_name:}.{sub_item.get_name():}"] = sub_item
        elif isinstance(item, types):
            d_found[item.get_name()] = item
    return d_found
