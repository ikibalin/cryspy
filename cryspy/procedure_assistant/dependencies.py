"""Add rows of any category together with the rows they depend on.

No list of features is kept anywhere: the class of a category says how its
rows depend on the other categories of the data block (ItemN.D_REFERENCES,
REQUIRED_CATEGORIES, D_PROPOSALS, propose_values). For a request as
"add multipoles for Mn and F" (category atom_rho_multipole) the assistant

1. finds the key of a row (atom_label, a value of atom_site.label) and
   the possible keys: the atoms without a row yet, first those that have
   the first required category (the magnetic atoms, with atom_site_scat);
2. keeps the rows that exist already and plans the missing ones, the
   rows of the required categories first (atom_local_axes of the atom);
3. fills the values from the definitions: the given ones ('Mn: 3d, F:
   2p' for several keys), the ones of the rows with the same key
   (atom_type_symbol from atom_site_scat), the proposals of the class
   (the orbital of the ion, the nearest atoms for the local axes) and the
   defaults;
4. asks, in one round, the values that are not sure, with the possible
   values and the reasons, and adds nothing before all are known.

The procedures and methods whose docstring section ``Adds`` names a
category (requirements.py) are given with it: they compute rows that a
user would have to type (dummy atoms of the local axes, the allowed
texture coefficients).
"""
from typing import Dict, List, Optional, Tuple

from cryspy.A_functions_base.function_1_strings import parse_values_by_label, split_list_text
from cryspy.B_parent_classes.attributes import take_allowed_values
from cryspy.B_parent_classes.cl_2_loop import LoopN
from cryspy.B_parent_classes.cl_3_data import DataN
from cryspy.exceptions import ArgumentsNeeded
from cryspy.procedure_assistant import items

MAX_OPTIONS = 6


def take_key_attribute(item_class: type) -> str:
    """Attribute that names the row (the first of D_REFERENCES), '' when the
    rows of the category do not refer to another one."""
    return next(iter(getattr(item_class, "D_REFERENCES", {}) or {}), "")


def split_reference(reference: str) -> Tuple[str, str]:
    category, _, column = reference.partition(".")
    return category, column


def take_rows(data: DataN, category: str) -> list:
    if not data.is_attribute(category):
        return []
    obj = getattr(data, category)
    return list(obj.items) if isinstance(obj, LoopN) else [obj]


def take_column(data: DataN, reference: str) -> List[str]:
    category, column = split_reference(reference)
    l_value = []
    for row in take_rows(data, category):
        value = getattr(row, column, None) if row.is_attribute(column) else None
        if value is not None and str(value) not in l_value:
            l_value.append(str(value))
    return l_value


def take_companion_key(companion_class: type, reference: str) -> str:
    """Attribute of the companion category with the same reference as the
    key ('label' of atom_site_scat for atom_site.label)."""
    category, column = split_reference(reference)
    if companion_class.PREFIX == category:
        return column
    for name, value in getattr(companion_class, "D_REFERENCES", {}).items():
        if value == reference:
            return name
    return take_key_attribute(companion_class)


def find_row(data: DataN, item_class: type, key_name: str, key: str):
    for row in take_rows(data, item_class.PREFIX):
        if row.is_attribute(key_name) and str(getattr(row, key_name)).lower() == key.lower():
            return row
    return None


def take_given_values(d_value: dict, prefix: str, item_class: type, l_key: List[str]) -> Dict[str, dict]:
    """{attribute: {key: value}} of the values given for a category: by the
    name of the attribute (main category) or 'category.attribute'
    (required category); 'Mn: 3d, F: 2p' gives each key its value, a plain
    value is for every key ("plain": True)."""
    d_out = {}
    for name, value in d_value.items():
        category, dot, attribute = str(name).rpartition(".")
        if (dot and category != prefix) or (not dot and prefix != ""):
            continue
        attribute_name = items.take_attribute_name(item_class, attribute)
        if attribute_name is None:
            if dot or name in ("category", "block", "values"):
                continue
            raise ValueError(f"{item_class.PREFIX:} has no attribute {attribute:}; its attributes: "
                             + ", ".join(item_class.ATTR_CIF) + ".")
        text = str(value)
        # a value without keys is "plain" only for several keys (one key: it is meant for it)
        plain = ":" not in text and "=" not in text and len(l_key) > 1
        d_out[attribute_name] = {"plain": plain, "values": parse_values_by_label(text, l_key)}
    return d_out


def find_methods_adding(cryspy_object, category: str) -> List[str]:
    """'add_multipole_form_factor of crystal_mnf2' for the methods and the
    procedures whose docstring section Adds names the category."""
    from cryspy.procedure_assistant.methods import iterate_method_blocks, take_methods
    from cryspy.procedure_assistant.procedures import take_procedures
    from cryspy.procedure_assistant.requirements import take_requirements
    l_out = []
    for block_path, block in iterate_method_blocks(cryspy_object):
        for name, function in take_methods(type(block)).items():
            if any(d_item["category"] == category for d_item in take_requirements(function)["adds"]):
                text = f"{name} of {block_path}"
                if not any(item.startswith(f"{name} of ") for item in l_out):
                    l_out.append(text)
    for name, function in take_procedures().items():
        if any(d_item["category"] == category for d_item in take_requirements(function)["adds"]):
            l_out.append(f"procedure {name}")
    return l_out


def take_proposed_from_rows(data: DataN, item_class: type, key: str, d_row: dict,
                            d_planned: Optional[dict] = None) -> dict:
    """Values of D_PROPOSALS for the attributes without a value: from the
    rows of other categories with the same key, or fixed values; returns
    the ones taken from rows (data the user may doubt)."""
    reference = item_class.D_REFERENCES[take_key_attribute(item_class)]
    d_data_value = {}
    for name, source in getattr(item_class, "D_PROPOSALS", {}).items():
        if name in d_row:
            continue
        if not (isinstance(source, str) and "." in source):
            d_row[name] = source
            continue
        category, column = split_reference(source)
        companion_class = items.take_category_class(category)
        row = find_row(data, companion_class, take_companion_key(companion_class, reference), key)
        d_planned_row = (d_planned or {}).get((category, key), {})
        if row is not None and row.is_attribute(column):
            d_row[name] = d_data_value[name] = str(getattr(row, column))
        elif column in d_planned_row:  # the row added in the same call
            d_row[name] = d_data_value[name] = str(d_planned_row[column])
    return d_data_value


def take_options(data: DataN, item_class: type, name: str) -> List[str]:
    """Possible values of an attribute: the referenced column or the
    allowed values of the class."""
    if name in item_class.D_REFERENCES:
        return take_column(data, item_class.D_REFERENCES[name])
    allowed = take_allowed_values(item_class, name)
    if allowed is not None and len(allowed) <= 30:
        return [str(value) for value in allowed]
    return []


def plan_row(data: DataN, item_class: type, key: str, d_given: Dict[str, dict],
             d_planned: Optional[dict] = None) -> Tuple[dict, dict]:
    """Values of a new row and the values to ask {attribute: {"proposal",
    "options", "reason"}}.

    Explicit values ('Mn: 3d') come first, then the values of rows with
    the same key and the proposals of the class; a proposal that is not
    sure is asked. A plain value given for several keys ('2p', the answer
    to a question about F) fills only what would be asked, or what has no
    value at all, so it does not replace the sure 3d of Mn."""
    key_name = take_key_attribute(item_class)
    d_row = {key_name: key}
    d_plain = {}
    for name, d_given_item in d_given.items():
        if key in d_given_item["values"] and name != key_name:
            (d_plain if d_given_item["plain"] else d_row)[name] = d_given_item["values"][key]
    l_explicit = list(d_row)
    d_data_value = take_proposed_from_rows(data, item_class, key, d_row, d_planned)
    d_open = {}
    for name, d_proposal in item_class.propose_values(data, d_row).items():
        if name in l_explicit:
            continue
        if d_proposal.get("sure", False) and d_proposal.get("proposal", "") != "":
            if name not in d_row or name in d_data_value:
                d_row[name] = d_proposal["proposal"]
        else:
            d_row.pop(name, None)
            d_open[name] = d_proposal
    for name in item_class.ATTR_MANDATORY_NAMES:
        if name not in d_row and name not in d_open:
            d_open[name] = {"proposal": "", "options": take_options(data, item_class, name),
                            "reason": ""}
    for name, value in d_plain.items():
        if name in d_open or name not in d_row:
            d_row[name] = value
            d_open.pop(name, None)
    return d_row, d_open


def plan_rows(cryspy_object, data: DataN, item_class: type, l_key: List[str],
              d_value: dict) -> dict:
    """Rows to add (required categories first), rows kept and questions.

    Returns
    -------
    dict
        {"rows": [(category, key, values)], "kept": {category: [key]},
        "open": {(category, attribute): {key: proposal dict}}}.
    """
    key_name = take_key_attribute(item_class)
    reference = item_class.D_REFERENCES[key_name]
    d_plan = {"rows": [], "kept": {}, "open": {}}

    def plan(klass: type, prefix: str, l_key_class: List[str]):
        d_given = take_given_values(d_value, prefix, klass, l_key_class)
        d_planned = {(category, key): d_row for category, key, d_row in d_plan["rows"]}
        for key in l_key_class:
            d_row, d_open = plan_row(data, klass, key, d_given, d_planned)
            # a wrong value (a word for a number) is an error before any question
            items.form_item(klass, f"{data.get_name()}.{klass.PREFIX}", d_row)
            d_plan["rows"].append((klass.PREFIX, key, d_row))
            for name, d_proposal in d_open.items():
                d_plan["open"].setdefault((klass.PREFIX, name), {})[key] = d_proposal

    l_new = []
    for key in l_key:
        if find_row(data, item_class, key_name, key) is not None:
            d_plan["kept"].setdefault(item_class.PREFIX, []).append(key)
        else:
            l_new.append(key)
    for category in getattr(item_class, "REQUIRED_CATEGORIES", ()):
        companion_class = items.take_category_class(category)
        if not items.accepts(data, companion_class):
            continue
        companion_key = take_companion_key(companion_class, reference)
        l_missing = [key for key in l_new if find_row(data, companion_class, companion_key, key) is None]
        d_plan["kept"].setdefault(category, []).extend(key for key in l_new if key not in l_missing)
        if l_missing:
            plan(companion_class, category, l_missing)
    plan(item_class, "", l_new)
    return d_plan


def form_questions(cryspy_object, item_class: type, d_open: dict) -> List[dict]:
    """One question per attribute for all its keys; the answer is 'key:
    value' for each of them (a plain value is for all)."""
    l_question = []
    for (category, name), d_key in d_open.items():
        arg_name = name if category == item_class.PREFIX else f"{category}.{name}"
        l_key = list(d_key)
        l_reason, l_proposal = [], []
        for key, d_proposal in d_key.items():
            s_option = ", ".join(str(option) for option in d_proposal.get("options", [])[:MAX_OPTIONS])
            text = f"{key}: {s_option}" if s_option else f"{key}"
            if d_proposal.get("reason"):
                text += f" ({d_proposal['reason']})"
            l_reason.append(text)
            if d_proposal.get("proposal", "") != "":
                l_proposal.append(f"{key}: {d_proposal['proposal']}")
        hint = "; ".join(l_reason)
        if len(l_key) > 1:
            hint += "; answer as 'key: value' for each, or one value for all"
        if category != item_class.PREFIX:
            l_method = find_methods_adding(cryspy_object, category)
            if l_method:
                hint += f"; {category} can also be made by {', '.join(l_method)}"
        d_question = {"name": arg_name, "type": "string", "hint": hint,
                      "question": f"{name} of {category} of {', '.join(l_key)}"}
        if len(l_proposal) == len(l_key):
            d_question["proposal"] = ", ".join(l_proposal)
            d_question["s_proposal"] = d_question["proposal"]
        options = d_key[l_key[0]].get("options", [])
        if len(l_key) == 1 and 1 < len(options) <= MAX_OPTIONS:
            d_question["options"] = [f"{l_key[0]}: {option}" for option in options]
        l_question.append(d_question)
    return l_question


def take_possible_keys(data: DataN, item_class: type) -> Tuple[List[str], List[str], Dict[str, str]]:
    """(preferred, other, reason of the other) keys of new rows: values of
    the referenced column without a row; preferred: those with a row of
    the first required category (the magnetic atoms for multipoles), all
    of them when the block has no row of it yet (texture_sh is added with the first
    coefficients)."""
    key_name = take_key_attribute(item_class)
    reference = item_class.D_REFERENCES[key_name]
    l_value = [value for value in take_column(data, reference)
               if find_row(data, item_class, key_name, value) is None]
    l_required = [category for category in getattr(item_class, "REQUIRED_CATEGORIES", ())
                  if items.accepts(data, items.take_category_class(category))]
    if not l_required:
        return l_value, [], {}
    first_class = items.take_category_class(l_required[0])
    first_key = take_companion_key(first_class, reference)
    l_preferred = [value for value in l_value
                   if find_row(data, first_class, first_key, value) is not None]
    if len(take_rows(data, l_required[0])) == 0:  # added with the new rows
        return l_value, [], {}
    l_other = [value for value in l_value if value not in l_preferred]
    return l_preferred, l_other, {value: f"no {l_required[0]}" for value in l_other}


def describe_category_rows(cryspy_object, data: DataN, item_class: type) -> str:
    """What the data block has of the category and of the categories it
    depends on, the possible keys of new rows and the procedures and
    methods that compute rows."""
    prefix = item_class.PREFIX
    key_name = take_key_attribute(item_class)
    ls_out = [f"Rows of {prefix} in {data.get_name()}:"]
    if key_name:
        reference = item_class.D_REFERENCES[key_name]
        l_present = [str(getattr(row, key_name)) for row in take_rows(data, prefix)
                     if row.is_attribute(key_name)]
        ls_out.append(f"- one row per {reference} ({key_name}); present: "
                      + (", ".join(l_present) or "none"))
        for category in getattr(item_class, "REQUIRED_CATEGORIES", ()):
            companion_class = items.take_category_class(category)
            companion_key = take_companion_key(companion_class, reference)
            l_row = [str(getattr(row, companion_key)) for row in take_rows(data, category)
                     if row.is_attribute(companion_key)]
            ls_out.append(f"- needs a row of {category} for the same {key_name}; present: "
                          + (", ".join(l_row) or "none") + " (the missing ones are added with it)")
        l_preferred, l_other, d_reason = take_possible_keys(data, item_class)
        ls_out.append(f"- new rows possible for: {', '.join(l_preferred) or 'none'}"
                      + (f"; {', '.join(l_other)}: {d_reason[l_other[0]]}" if l_other else ""))
    l_method = find_methods_adding(cryspy_object, prefix)
    if l_method:
        ls_out.append(f"- also made by {', '.join(l_method)} (computes the rows)")
    return "\n".join(ls_out)


def add_rows(cryspy_object, category: str, d_value: dict, block: str = "") -> Optional[str]:
    """Add rows of the category for one or several keys ('Mn, F') with the
    missing rows of the required categories; None for a category whose
    rows refer to no other one (items.add_item adds it)."""
    item_class = items.take_category_class(category)
    key_name = take_key_attribute(item_class)
    if key_name == "":
        return None
    data = items.choose_data_block(cryspy_object, item_class, block)
    reference = item_class.D_REFERENCES[key_name]
    l_all = take_column(data, reference)
    if len(l_all) == 0:
        ref_category = split_reference(reference)[0]
        raise ValueError(f"{data.get_name()} has no {ref_category}: add it first.")
    given = d_value.get(key_name, d_value.get(items.take_attribute_name(item_class, key_name), ""))
    l_key = split_list_text(given)
    if len(l_key) == 0:
        l_preferred, l_other, d_reason = take_possible_keys(data, item_class)
        if len(l_preferred) == 0:
            raise ValueError(
                f"no new row of {item_class.PREFIX} is possible in {data.get_name()}:\n"
                + describe_category_rows(cryspy_object, data, item_class))
        d_question = {
            "name": key_name, "type": "string",
            "question": f"{key_name} of the new rows of {item_class.PREFIX} (one or several)",
            "hint": ", ".join(l_preferred)
            + (f"; {', '.join(l_other)}: {d_reason[l_other[0]]}" if l_other else "")}
        if len(l_preferred) == 1:
            d_question.update(proposal=l_preferred[0], s_proposal=l_preferred[0])
        elif len(l_preferred) < MAX_OPTIONS:
            d_question["options"] = [", ".join(l_preferred)] + l_preferred
        raise ArgumentsNeeded([d_question])
    l_unknown = [key for key in l_key if key.lower() not in [value.lower() for value in l_all]]
    if l_unknown:
        raise ValueError(f"no {', '.join(l_unknown)} in {reference}; the values are {', '.join(l_all)}.")
    d_lower = {value.lower(): value for value in l_all}
    l_key = [d_lower[key.lower()] for key in l_key]
    d_plan = plan_rows(cryspy_object, data, item_class, l_key, d_value)
    l_kept = d_plan["kept"].get(item_class.PREFIX, [])
    if len(d_plan["rows"]) == 0:
        return f"{', '.join(l_kept)}: {item_class.PREFIX} is there already, nothing is added."
    if d_plan["open"]:
        ls_plan = [f"{category} of {key}" for category, key, _ in d_plan["rows"]]
        raise ArgumentsNeeded(form_questions(cryspy_object, item_class, d_plan["open"]),
                              "To add: " + ", ".join(ls_plan) + ".")
    ls_out = []
    for row_category, key, d_row in d_plan["rows"]:
        path, _ = items.add_item(cryspy_object, row_category, d_row, data.get_name())
        ls_out.append(f"{path} is added: " + ", ".join(f"{name} {value}" for name, value in d_row.items()))
    for row_category, l_key_kept in d_plan["kept"].items():
        if l_key_kept:
            ls_out.append(f"{row_category} of {', '.join(l_key_kept)}: there already, kept.")
    return "\n".join(ls_out)
