"""Add and remove items of any category in the data blocks of an object.

A row of a loop (atom_rho_multipole, atom_local_axes, atom_site_scat, ...)
or a single item (extinction, setup, ...) is added to a data block by the
name of its category; the loop is created when the data block has none.
Only the class definitions are used (ATTR_NAMES, ATTR_MANDATORY_NAMES,
ATTR_CIF, D_DEFAULT, CLASSES of the data blocks, ITEM_CLASS of the loops),
so new categories need no change here:

    >>> add_item(obj, "atom_rho_multipole", {"atom_label": "Fe",
    ...          "atom_type_symbol": "Fe3+", "orbital": "3d"})
    >>> remove_item(obj, "crystal_fe.atom_rho_multipole[Fe]")
"""
from typing import List, Optional, Tuple

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN
from cryspy.B_parent_classes.cl_3_data import DataN
from cryspy.B_parent_classes.cl_4_global import GlobalN
from cryspy.procedure_assistant.attributes import (
    AttributeRecord, find_item_classes, take_block)
from cryspy.procedure_assistant.parameters import take_item_label

class BlockNeeded(Exception):
    """The data block should be chosen; args[0]: names of the data blocks."""


class ValuesNeeded(Exception):
    """Mandatory attributes have no value; args: their names."""


# attributes which hold the label of an atom of atom_site
ATOM_LABEL_ATTRIBUTES = ("atom_label", )


def take_data_blocks(cryspy_object) -> List[DataN]:
    if isinstance(cryspy_object, GlobalN):
        return [item for item in cryspy_object.items if isinstance(item, DataN)]
    return [cryspy_object] if isinstance(cryspy_object, DataN) else []


def take_loop_class(data: DataN, item_class: type) -> Optional[type]:
    """Class of the loop of the data block whose rows are item_class."""
    for block_class in data.CLASSES:
        if issubclass(block_class, LoopN) and getattr(block_class, "ITEM_CLASS", None) is item_class:
            return block_class
    return None


def accepts(data: DataN, item_class: type) -> bool:
    """True when the data block may have the item (or a loop of them)."""
    return item_class in data.CLASSES or take_loop_class(data, item_class) is not None


def take_category_class(category: str) -> type:
    """ItemN class of the category ('atom_rho_multipole', 'AtomRhoMultipole',
    '_atom_rho_multipole'). Raises ValueError with the explanation."""
    l_class = find_item_classes(category.strip().lstrip("_"))
    if len(l_class) == 0:
        raise ValueError(f"there is no category '{category:}'. describe_category lists the categories "
                         "of the object.")
    if len(l_class) > 1:
        raise ValueError("several categories match, give one of: "
                         + ", ".join(item_class.PREFIX for item_class in l_class[:20]))
    return l_class[0]


def take_attribute_name(item_class: type, key: str) -> Optional[str]:
    """Attribute of the class given by its name, its CIF name or its CIF
    tag ('coeff_P1m1', 'coeff_P1-1', '_atom_rho_multipole_coeff_P1-1')."""
    key_sh = key.strip().lstrip("_")
    prefix = item_class.PREFIX
    if key_sh.lower().startswith(prefix.lower() + "_"):
        key_sh = key_sh[len(prefix) + 1:]
    for name, name_cif in zip(item_class.ATTR_NAMES, item_class.ATTR_CIF):
        if key_sh in (name, name_cif):
            return name
    for name, name_cif in zip(item_class.ATTR_NAMES, item_class.ATTR_CIF):
        if key_sh.lower() in (name.lower(), name_cif.lower()):
            return name
    return None


def form_item(item_class: type, path: str, values: dict) -> Tuple[ItemN, List[str]]:
    """Item with the values (converted to the types of the attributes) and
    the names of the mandatory attributes without a value."""
    item = item_class()
    l_unknown = []
    for key, value in values.items():
        name = take_attribute_name(item_class, key)
        if name is None:
            l_unknown.append(key)
            continue
        AttributeRecord(f"{path:}.{name:}", item, name).set_value(value)
    if l_unknown:
        raise ValueError(
            f"{item_class.PREFIX:} has no attribute {', '.join(l_unknown):}; its attributes: "
            + ", ".join(item_class.ATTR_CIF) + ".")
    l_missing = [name for name in item_class.ATTR_MANDATORY_NAMES if not item.is_attribute(name)]
    return item, l_missing


def take_atom_labels(data: DataN) -> Optional[List[str]]:
    atom_site = next((obj for obj in data.items if obj.get_name() == "atom_site"), None)
    if not isinstance(atom_site, LoopN):
        return None
    return [str(atom.label) for atom in atom_site.items]


def check_atom_labels(data: DataN, item: ItemN):
    l_label = take_atom_labels(data)
    if l_label is None:
        return
    for name in ATOM_LABEL_ATTRIBUTES:
        if item.is_attribute(name) and str(getattr(item, name)) not in l_label:
            raise ValueError(f"atom '{getattr(item, name):}' is not in atom_site of {data.get_name():}; "
                             f"the atoms are {', '.join(l_label):}.")


def fill_columns(loop: LoopN, item: ItemN):
    """Rows of a loop are printed with the columns of the first row: an
    attribute defined in some rows is given to the others with its default
    value (0 for numbers without a default, as a missing coefficient), so
    no column is lost in the rcif file."""
    l_item = loop.items + [item]
    for name, value_type in zip(item.ATTR_NAMES, item.ATTR_TYPES):
        l_defined = [row for row in l_item if row.is_attribute(name)]
        if len(l_defined) in (0, len(l_item)):
            continue
        if name in item.D_DEFAULT:
            value = item.D_DEFAULT[name]
        elif value_type in (float, int):
            value = value_type(0)
        else:
            continue  # printed as '.'
        for row in l_item:
            if not row.is_attribute(name):
                setattr(row, name, value)


def add_item(cryspy_object, category: str, values: dict, block: str = "") -> Tuple[str, object]:
    """Add an item of the category to a data block; for a loop category the
    item is a new row (the loop is created when there is none).

    Parameters
    ----------
    cryspy_object : GlobalN or DataN
        The object.
    category : str
        Name of the category, e.g. 'atom_rho_multipole'.
    values : dict
        Values of the attributes by their names or CIF names.
    block : str
        Name of the data block; empty when only one data block may have
        the category.

    Returns
    -------
    (path, block) : the path and the loop (or item) where the item is.

    Raises ValueError with the explanation, BlockNeeded when the data
    block should be chosen and ValuesNeeded when mandatory attributes have
    no value.
    """
    item_class = take_category_class(category)
    l_data = [data for data in take_data_blocks(cryspy_object) if accepts(data, item_class)]
    if block.strip() != "":
        data = take_block(cryspy_object, block.strip())
        if not isinstance(data, DataN):
            raise ValueError(f"there is no data block '{block:}'.")
        if data not in l_data:
            raise ValueError(f"the data block {data.get_name():} ({type(data).__name__:}) cannot have "
                             f"{item_class.PREFIX:}.")
    elif len(l_data) == 1:
        data = l_data[0]
    elif len(l_data) == 0:
        raise ValueError(f"no data block of the object can have {item_class.PREFIX:}.")
    else:
        raise BlockNeeded([data.get_name() for data in l_data])
    prefix = item_class.PREFIX
    item, l_missing = form_item(item_class, f"{data.get_name():}.{prefix:}", values)
    if l_missing:
        raise ValuesNeeded(*l_missing)
    check_atom_labels(data, item)
    loop_class = take_loop_class(data, item_class)
    container = next((obj for obj in data.items if isinstance(obj, loop_class or item_class)), None)
    path = f"{data.get_name():}.{prefix:}"
    if loop_class is None:
        if container is not None:
            raise ValueError(f"{path:} exists already; change its values with set_attribute.")
        data.items.append(item)
        container = item
    else:
        label = take_item_label(item, len(getattr(container, "items", [])))
        if container is None:
            container = loop_class()
            data.items.append(container)
        elif label in [take_item_label(row, index) for index, row in enumerate(container.items)]:
            raise ValueError(f"{path:}[{label:}] exists already; change its values with set_attribute.")
        fill_columns(container, item)
        container.items.append(item)
        path = f"{path:}[{label:}]"
    if hasattr(data, "apply_constraints"):
        data.apply_constraints()
    return path, container


def take_parent(cryspy_object, block) -> Tuple[Optional[DataN], Optional[LoopN]]:
    """Data block and loop (None for a single item) which hold the block."""
    for data in take_data_blocks(cryspy_object):
        for obj in data.items:
            if obj is block:
                return data, None
            if isinstance(obj, LoopN) and any(row is block for row in obj.items):
                return data, obj
    return None, None


def remove_item(cryspy_object, path: str) -> str:
    """Remove a row of a loop ('crystal_fe.atom_rho_multipole[Fe]'), a
    whole loop or a single item ('crystal_fe.atom_rho_multipole',
    'exp.extinction') from its data block. A loop left without rows is
    removed too; mandatory categories of a data block are kept. Returns
    the removed path; raises ValueError with the explanation."""
    block = take_block(cryspy_object, path)
    if block is None or isinstance(block, DataN):
        raise ValueError(f"'{path:}' is not a row, loop or item of a data block. Give e.g. "
                         "'crystal_fe.atom_rho_multipole[Fe]' or 'crystal_fe.atom_rho_multipole'.")
    data, loop = take_parent(cryspy_object, block)
    if data is None:
        raise ValueError(f"'{path:}' is not in a data block.")
    if loop is not None:
        index = next(index for index, row in enumerate(loop.items) if row is block)
        label = take_item_label(block, index)
        loop.items.pop(index)
        if len(loop.items) > 0:
            return f"{data.get_name():}.{loop.get_name():}[{label:}]"
        block = loop
    if type(block) in data.CLASSES_MANDATORY:
        raise ValueError(f"{block.get_name():} is mandatory in {data.get_name():} and is kept.")
    data.items.remove(block)
    return f"{data.get_name():}.{block.get_name():}"
