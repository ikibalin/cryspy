"""Description of the attributes of ItemN and LoopN classes.

Everything is taken from the class definition (ATTR_NAMES, ATTR_TYPES,
ATTR_CIF, ATTR_REF, D_CONSTRAINTS, PREFIX), so it works for every class
derived from ItemN or LoopN, also for the classes added later:

    >>> import cryspy
    >>> cryspy.take_value_choices(cryspy.AtomSiteSusceptibilityL)
    {'_atom_site_susceptibility_chi_type': ('Ciso', 'Cani')}
    >>> cryspy.next_value_choice(('Ciso', 'Cani'), "ciso")
    'Cani'

Front ends (CrysPy Editor, AI assistants) use these functions instead of
lists of their own.
"""
from typing import List, Optional

from cryspy.B_parent_classes.cl_1_item import ItemN

# more allowed values than this are not a choice to switch by a click
# (e.g. the space group names); sets have no order and are no choice either
MAX_VALUE_CHOICES = 12


def take_item_class(obj) -> Optional[type]:
    """Class derived from ItemN of an ItemN or LoopN object or class;
    None when it is not known (e.g. LoopN without ITEM_CLASS)."""
    if isinstance(obj, type):
        if issubclass(obj, ItemN):
            item_class = obj
        else:
            item_class = getattr(obj, "ITEM_CLASS", None)
    elif isinstance(obj, ItemN):
        item_class = type(obj)
    else:  # no getattr: LoopN.__getattr__ recurses without ITEM_CLASS
        item_class = getattr(obj, "__dict__", {}).get(
            "ITEM_CLASS", getattr(type(obj), "ITEM_CLASS", None))
    if not isinstance(item_class, type) or item_class is ItemN:
        return None
    if not isinstance(getattr(item_class, "ATTR_NAMES", None), (tuple, list)):
        return None
    return item_class


def normalize_cif_tag(tag: str) -> str:
    """'_Atom_site_susceptibility.chi_type' -> '_atom_site_susceptibility_chi_type'."""
    return tag.strip().lower().replace(".", "_")


def take_cif_tag(item_class: type, name: str) -> str:
    """Normalized CIF tag of the attribute name of the class."""
    l_name = list(item_class.ATTR_NAMES)
    l_cif = list(getattr(item_class, "ATTR_CIF", l_name))
    name_cif = l_cif[l_name.index(name)] if name in l_name else name
    return normalize_cif_tag(f"_{item_class.PREFIX:}_{name_cif or name:}")


def take_allowed_values(item_class: type, name: str) -> Optional[tuple]:
    """All allowed values of the attribute (D_CONSTRAINTS) or None."""
    d_constraints = getattr(item_class, "D_CONSTRAINTS", {})
    if not isinstance(d_constraints, dict) or name not in d_constraints:
        return None
    return tuple(d_constraints[name])


def take_attribute_choices(item_class: type, name: str) -> Optional[tuple]:
    """Allowed values of a text attribute when they are a short ordered
    list (e.g. chi_type: Ciso, Cani), else None."""
    d_constraints = getattr(item_class, "D_CONSTRAINTS", {})
    if not isinstance(d_constraints, dict):
        return None
    allowed = d_constraints.get(name)
    if not isinstance(allowed, (tuple, list)):
        return None
    if not 1 < len(allowed) <= MAX_VALUE_CHOICES:
        return None
    if not all(isinstance(value, str) for value in allowed):
        return None
    return tuple(allowed)


def take_value_choices(obj) -> dict:
    """Normalized CIF tag -> allowed values (a short ordered list) for the
    text attributes of an ItemN or LoopN object or class."""
    item_class = take_item_class(obj)
    if item_class is None:
        return {}
    d_choice = {}
    for name in item_class.ATTR_NAMES:
        choices = take_attribute_choices(item_class, name)
        if choices is not None:
            d_choice[take_cif_tag(item_class, name)] = choices
    return d_choice


def take_refinable_cif_tags(obj) -> set:
    """Normalized CIF tags of the refinable attributes (ATTR_REF) of an
    ItemN or LoopN object or class."""
    item_class = take_item_class(obj)
    if item_class is None:
        return set()
    attr_ref = getattr(item_class, "ATTR_REF", ())
    if not isinstance(attr_ref, (tuple, list)):  # a property on a class
        return set()
    return {take_cif_tag(item_class, name) for name in item_class.ATTR_NAMES
            if name in attr_ref}


def next_value_choice(choices, value: str) -> str:
    """Value after value in choices, the first after the last; the
    comparison ignores case. The first value when value is not there."""
    l_choice = list(choices)
    value_sh = str(value).strip().strip("'\"").lower()
    for index, choice in enumerate(l_choice):
        if choice.lower() == value_sh:
            return l_choice[(index + 1) % len(l_choice)]
    return l_choice[0]


def describe_attributes(obj) -> List[dict]:
    """Attributes of an ItemN or LoopN object or class.

    Each dict has 'name', 'cif_tag', 'type' (name of the type),
    'mandatory', 'refinable' and 'allowed_values' (a tuple or None; a long
    list or a set is given as it is defined)."""
    item_class = take_item_class(obj)
    if item_class is None:
        return []
    l_type = list(getattr(item_class, "ATTR_TYPES", ()))
    attr_mandatory = getattr(item_class, "ATTR_MANDATORY_NAMES", ())
    attr_ref = getattr(item_class, "ATTR_REF", ())
    if not isinstance(attr_ref, (tuple, list)):
        attr_ref = ()
    l_description = []
    for index, name in enumerate(item_class.ATTR_NAMES):
        val_type = l_type[index] if index < len(l_type) else str
        allowed = take_allowed_values(item_class, name)
        if allowed is not None and not isinstance(
                item_class.D_CONSTRAINTS[name], (tuple, list)):
            allowed = tuple(sorted(allowed, key=str))
        l_description.append({
            "name": name,
            "cif_tag": take_cif_tag(item_class, name),
            "type": getattr(val_type, "__name__", str(val_type)),
            "mandatory": name in attr_mandatory,
            "refinable": name in attr_ref,
            "allowed_values": allowed})
    return l_description


def take_item_classes() -> List[type]:
    """All classes derived from ItemN that are imported (also those of
    other packages), each once, sorted by PREFIX."""
    l_class, l_new = [], list(ItemN.__subclasses__())
    while len(l_new) > 0:
        item_class = l_new.pop(0)
        if item_class in l_class:
            continue
        l_class.append(item_class)
        l_new.extend(item_class.__subclasses__())
    l_class = [item_class for item_class in l_class
               if isinstance(getattr(item_class, "PREFIX", None), str)]
    return sorted(l_class, key=lambda item_class: item_class.PREFIX)
