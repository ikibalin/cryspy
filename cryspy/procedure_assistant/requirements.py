"""CIF categories that a procedure or a method needs and adds.

A function that adds a description to the object (multipole form factor,
texture by spherical harmonics, magnetic broadening, ...) names in its
numpydoc docstring the categories it needs in the object and the ones it
adds, with the column given by one of its parameters:

    Needs
    -----
    atom_site_scat : label = atom_label
        The magnetic ion of the atom.

    Adds
    ----
    atom_rho_multipole : atom_label = atom_label
        The multipole coefficients of the atom.
    atom_local_axes : atom_label = atom_label
        The local axes, when the atom has none.

From these sections alone, with no list of features anywhere,
:func:`check_requirements` looks at the object and finds what is there
already and what is missing, and the possible values of the parameters:
the values of the first needed category (atom_site.label), the ones that
lack a needed category (dummy atoms without atom_site_scat) and the ones
that have the added category already. An AI assistant (and any program)
then adds only what is missing, proposes the possible values and asks the
user only when it is not clear which one is meant:

    >>> cryspy.call_assistant_tool(obj, "check_requirements",
    ...                            {"name": "add_multipole_form_factor"})
"""
from typing import Callable, Dict, List, Optional, Tuple

from cryspy.A_functions_base.function_1_strings import split_list_text as split_values
from cryspy.B_parent_classes.cl_2_loop import LoopN
from cryspy.B_parent_classes.cl_3_data import DataN
from cryspy.B_parent_classes.cl_4_global import GlobalN
from cryspy.procedure_assistant.procedures import parse_docstring

# a parameter whose description says "one or several" takes several values
# separated by commas ('Mn, F')
SEVERAL_WORDS = ("one or several", "several labels", "several values", "separated by commas")
MAX_VALUES_IN_TEXT = 12
# more possible values than this are given in the hint only, not as buttons
MAX_OPTIONS = 5


def take_requirements(function: Callable) -> Dict[str, List[dict]]:
    """{"needs": [...], "adds": [...]} of the docstring of the function
    (procedures.parse_requirement_lines); empty lists when it has none.

    The categories needed by the rows a function adds come from their
    classes (dependencies.py): an added row keyed by a parameter needs the
    referenced row (atom_site : label = atom_label) and the rows of its
    REQUIRED_CATEGORIES that the function does not add itself
    (atom_site_scat for atom_rho_multipole), so the section Needs lists
    only what the classes cannot say."""
    d_doc = parse_docstring(getattr(function, "__doc__", "") or "")
    l_need = list(d_doc["needs"])
    l_add_category = [d_item["category"] for d_item in d_doc["adds"]]
    for d_add in d_doc["adds"]:
        if not d_add["parameter"]:
            continue
        for d_need in take_class_needs(d_add):
            if d_need["category"] in l_add_category or any(
                    d_item["category"] == d_need["category"] for d_item in l_need):
                continue
            l_need.append(d_need)
    return {"needs": l_need, "adds": d_doc["adds"]}


def take_class_needs(d_add: dict) -> List[dict]:
    """Categories needed by an added keyed category, from its class."""
    from cryspy.procedure_assistant.dependencies import split_reference, take_companion_key
    from cryspy.procedure_assistant.items import take_category_class
    try:
        item_class = take_category_class(d_add["category"])
    except ValueError:
        return []
    reference = getattr(item_class, "D_REFERENCES", {}).get(d_add["column"])
    if reference is None:
        return []
    category, column = split_reference(reference)
    l_out = [{"category": category, "column": column, "parameter": d_add["parameter"],
              "description": f"referenced by {d_add['category']}.{d_add['column']}"}]
    for required in getattr(item_class, "REQUIRED_CATEGORIES", ()):
        required_class = take_category_class(required)
        l_out.append({"category": required,
                      "column": take_companion_key(required_class, reference),
                      "parameter": d_add["parameter"],
                      "description": f"needed by the rows of {d_add['category']}"})
    return l_out


def has_requirements(function: Callable) -> bool:
    d_requirement = take_requirements(function)
    return bool(d_requirement["needs"] or d_requirement["adds"])


def takes_several_values(description: str) -> bool:
    return any(word in description.lower() for word in SEVERAL_WORDS)


def take_data_blocks(cryspy_object, block=None) -> List[DataN]:
    """Data blocks where the categories are looked for: the data block of
    the method (the block itself or the one that holds the loop or item),
    else every data block of the object (procedures)."""
    if isinstance(block, DataN):
        return [block]
    l_data = []
    if isinstance(cryspy_object, GlobalN):
        l_data = [item for item in cryspy_object.items if isinstance(item, DataN)]
    elif isinstance(cryspy_object, DataN):
        l_data = [cryspy_object]
    if block is not None:
        l_holder = [data for data in l_data if any(item is block for item in data.items)]
        if l_holder:
            return l_holder
    return l_data


def take_column_values(cif_object, column: str) -> List[str]:
    """Values of a column of a loop (or of an item) as text."""
    l_item = cif_object.items if isinstance(cif_object, LoopN) else [cif_object]
    l_value = []
    for item in l_item:
        value = getattr(item, column, None)
        if value is not None and str(value) not in l_value:
            l_value.append(str(value))
    return l_value


def find_category(l_data: List[DataN], d_requirement: dict) -> Dict[str, Optional[List[str]]]:
    """{data block name: values of the column (or [] without a column)} of
    the data blocks that have the category."""
    d_found = {}
    for data in l_data:
        if not data.is_attribute(d_requirement["category"]):
            continue
        cif_object = getattr(data, d_requirement["category"])
        d_found[data.get_name()] = (take_column_values(cif_object, d_requirement["column"])
                                    if d_requirement["column"] else [])
    return d_found


def contains_value(l_value: List[str], value: str) -> bool:
    return value.lower() in [item.lower() for item in l_value]


def check_requirements(cryspy_object, function: Callable, block=None,
                       arguments: Optional[dict] = None) -> dict:
    """What the object has of the categories the function needs and adds.

    Parameters
    ----------
    cryspy_object : GlobalN or DataN
        The object.
    function : callable
        The procedure or the method (its docstring has the sections Needs
        and Adds).
    block : DataN, LoopN or ItemN, optional
        The block of a method; its data block is searched. None for a
        procedure: every data block is searched.
    arguments : dict, optional
        The values given already: each value of a keyed parameter is
        checked.

    Returns
    -------
    dict
        {"needs": [requirement + {"blocks": {data name: values}}], "adds":
        [...], "choices": {parameter: {"values", "eligible", "lacking"
        {value: [category]}, "done" {value: [category]}, "several"}},
        "errors": [text]}. "values" are the values of the first needed
        category with this parameter; "eligible" the ones that have every
        needed category and not yet the first added one.
    """
    l_data = take_data_blocks(cryspy_object, block)
    d_requirement = take_requirements(function)
    d_doc_param = parse_docstring(getattr(function, "__doc__", "") or "")["parameters"]
    d_out = {"needs": [], "adds": [], "choices": {}, "errors": []}
    for key in ("needs", "adds"):
        for d_item in d_requirement[key]:
            d_out[key].append(dict(d_item, blocks=find_category(l_data, d_item)))
    l_param = []
    for d_item in d_out["needs"] + d_out["adds"]:
        if d_item["parameter"] and d_item["parameter"] not in l_param:
            l_param.append(d_item["parameter"])
    for param in l_param:
        description = d_doc_param.get(param, {}).get("description", "")
        d_out["choices"][param] = take_choices(d_out, param, takes_several_values(description))
    for param, value in (arguments or {}).items():
        if param not in d_out["choices"] or not d_out["choices"][param]["has_values"]:
            continue
        d_out["errors"].extend(check_given_values(d_out, param, value))
    return d_out


def take_choices(d_check: dict, param: str, several: bool) -> dict:
    """Possible values of a keyed parameter (check_requirements)."""
    l_need = [d_item for d_item in d_check["needs"] if d_item["parameter"] == param]
    l_add = [d_item for d_item in d_check["adds"] if d_item["parameter"] == param]
    l_value = []
    for values in (l_need[0]["blocks"].values() if l_need else []):
        l_value.extend(value for value in values if not contains_value(l_value, value))
    d_lacking, d_done = {}, {}
    for value in l_value:
        for d_item in l_need[1:]:
            if not any(contains_value(values, value) for values in d_item["blocks"].values()):
                d_lacking.setdefault(value, []).append(d_item["category"])
        for d_item in l_add:
            if any(contains_value(values, value) for values in d_item["blocks"].values()):
                d_done.setdefault(value, []).append(d_item["category"])
    main_add = l_add[0]["category"] if l_add else ""
    l_eligible = [value for value in l_value
                  if value not in d_lacking and main_add not in d_done.get(value, [])]
    return {"values": l_value, "eligible": l_eligible, "lacking": d_lacking,
            "done": d_done, "several": several, "has_values": len(l_need) > 0}


def check_given_values(d_check: dict, param: str, value) -> List[str]:
    """Errors for the values of a keyed parameter: unknown values and values
    that lack a needed category (to be added first)."""
    d_choice = d_check["choices"][param]
    l_value = split_values(value)
    if len(l_value) == 0:
        return []  # empty: the function takes every value (as phase_label='')
    l_error = []
    l_need = [d_item for d_item in d_check["needs"] if d_item["parameter"] == param]
    l_unknown = [item for item in l_value if not contains_value(d_choice["values"], item)]
    if l_unknown:
        d_first = l_need[0]
        l_error.append(
            f"no {', '.join(l_unknown)} in {d_first['category']}.{d_first['column']} "
            f"({param}); the values are {format_values(d_choice['values'])}.")
    if len(l_value) > 1 and not d_choice["several"]:
        l_error.append(f"{param} takes one value, not {', '.join(l_value)}; call it for each.")
    d_lower = {key.lower(): categories for key, categories in d_choice["lacking"].items()}
    for item in l_value:
        for category in d_lower.get(item.lower(), []):
            d_need = next(d_need for d_need in l_need if d_need["category"] == category)
            reason = f" ({d_need['description']})" if d_need["description"] else ""
            l_error.append(
                f"{item} has no row of {category}{reason}; add it first "
                f"(add_item {{\"category\": \"{category}\"}}).")
    return l_error


def format_values(l_value: List[str]) -> str:
    if len(l_value) == 0:
        return "none"
    s_value = ", ".join(l_value[:MAX_VALUES_IN_TEXT])
    return s_value + (", ..." if len(l_value) > MAX_VALUES_IN_TEXT else "")


def add_requirement_choices(d_check: dict, l_missing: List[dict]) -> Tuple[List[dict], str]:
    """Questions of the missing keyed parameters with the possible values
    as options and, when it is clear, a proposal: the only eligible value,
    or all of them for a parameter that takes several. Returns the
    questions and an error text when no value is possible ('' otherwise)."""
    l_out = []
    for d_missing in l_missing:
        d_choice = d_check["choices"].get(d_missing["name"])
        if d_choice is None or not d_choice["has_values"]:
            l_out.append(d_missing)
            continue
        l_eligible = d_choice["eligible"]
        if len(l_eligible) == 0:
            return l_missing, form_no_value_text(d_check, d_missing["name"])
        d_missing = dict(d_missing, hint=form_choice_text(d_check, d_missing["name"]))
        if len(l_eligible) == 1:
            d_missing.update(proposal=l_eligible[0], s_proposal=l_eligible[0])
        elif len(l_eligible) > MAX_OPTIONS:
            pass
        elif d_choice["several"]:
            d_missing["options"] = [", ".join(l_eligible)] + l_eligible
        else:
            d_missing["options"] = list(l_eligible)
        l_out.append(d_missing)
    return l_out, ""


def form_choice_text(d_check: dict, param: str) -> str:
    """'Mn, F can get it; D1, D2 lack atom_site_scat; Cu has atom_rho_multipole already'."""
    d_choice = d_check["choices"][param]
    ls_out = [f"{format_values(d_choice['eligible'])} "
              + ("(one or several)" if d_choice["several"] else "(one of them)")]
    d_group = {}
    for value, categories in d_choice["lacking"].items():
        d_group.setdefault(f"lack {', '.join(categories)}", []).append(value)
    for value, categories in d_choice["done"].items():
        d_group.setdefault(f"have {', '.join(categories)} already", []).append(value)
    for text, l_value in d_group.items():
        ls_out.append(f"{format_values(l_value)} {text}")
    return "; ".join(ls_out)


def form_no_value_text(d_check: dict, param: str) -> str:
    d_choice = d_check["choices"][param]
    text = f"Error: no value of {param} is possible"
    if d_choice["values"]:
        return f"{text}: " + form_choice_text(d_check, param).split("; ", 1)[-1] + "."
    l_need = [d_item for d_item in d_check["needs"] if d_item["parameter"] == param]
    return (f"{text}: the object has no {l_need[0]['category']}; add it first "
            f"(add_item {{\"category\": \"{l_need[0]['category']}\"}}).")


def format_requirement(d_item: dict) -> str:
    key = f" ({d_item['column']} = {d_item['parameter']})" if d_item["column"] else ""
    return f"{d_item['category']}{key}"


def form_requirements_text(d_check: dict, title: str = "", flag_short: bool = False) -> str:
    """Markdown of check_requirements: what is in the object, what the
    function adds, and the possible values of its parameters."""
    ls_out = [f"Requirements of {title}:"] if title else []
    for key, word in (("needs", "needs"), ("adds", "adds")):
        for d_item in d_check[key]:
            d_block = d_item["blocks"]
            if len(d_block) == 0:
                state = "not in the object" + (" (will be added)" if key == "adds" else
                                               " (add it first)")
            elif d_item["column"]:
                l_value = [value for values in d_block.values() for value in values]
                state = f"rows {format_values(l_value)}"
            else:
                state = f"present in {', '.join(d_block)}"
            description = "" if flag_short or not d_item["description"] \
                else f" - {d_item['description']}"
            ls_out.append(f"- {word} {format_requirement(d_item)}: {state}{description}")
    for param, d_choice in d_check["choices"].items():
        if d_choice["has_values"]:
            ls_out.append(f"- {param}: {form_choice_text(d_check, param)}")
    return "\n".join(ls_out)
