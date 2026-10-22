"""All attributes of a CrysPy object addressed by readable paths.

The paths are those of parameters.py (refinable parameters), but every
attribute of every ItemN is reached, also text ones such as
atom_site_susceptibility.chi_type or extinction.model:

    crystal.atom_site_susceptibility[Fe].chi_type
    exp.extinction.model

Only the class definitions (ATTR_NAMES, ATTR_TYPES, D_CONSTRAINTS, ...)
are used, so the objects of new classes derived from DataN, ItemN and
LoopN are reached without any change here.
"""
import re
from typing import Iterator, List, Optional, Tuple

from cryspy.B_parent_classes.attributes import (
    describe_attributes, take_allowed_values, take_item_class, take_item_classes)
from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN
from cryspy.B_parent_classes.cl_3_data import DataN
from cryspy.B_parent_classes.cl_4_global import GlobalN
from cryspy.procedure_assistant.parameters import (
    form_path_regex, remove_labels, split_pattern, take_item_label)


class AttributeRecord:
    """One attribute of one item: its path and where it is stored."""

    def __init__(self, path: str, item: ItemN, attribute: str):
        self.path = path
        self.item = item
        self.attribute = attribute

    @property
    def value(self):
        if not self.item.is_attribute(self.attribute):
            return None
        return getattr(self.item, self.attribute)

    @property
    def allowed_values(self) -> Optional[tuple]:
        return take_allowed_values(type(self.item), self.attribute)

    @property
    def value_type(self) -> type:
        ind = self.item.ATTR_NAMES.index(self.attribute)
        return self.item.ATTR_TYPES[ind]

    @property
    def flag_constraint(self) -> bool:
        if self.attribute not in self.item.ATTR_REF:
            return False
        return bool(getattr(self.item, f"{self.attribute:}_constraint"))

    def convert_value(self, value):
        """Value of the type of the attribute; an allowed value is found
        ignoring case. Raises ValueError with the explanation."""
        allowed = self.allowed_values
        if allowed is not None:
            value_sh = str(value).strip().strip("'\"").lower()
            for allowed_value in allowed:
                if str(allowed_value).lower() == value_sh:
                    return allowed_value
            s_allowed = ", ".join(str(_) for _ in list(allowed)[:30])
            if len(allowed) > 30:
                s_allowed += ", ..."
            raise ValueError(
                f"'{value:}' is not allowed for {self.path:}; allowed values: {s_allowed:}.")
        value_type = self.value_type
        try:
            if value_type is bool and isinstance(value, str):
                return value.strip().lower() in ("true", "yes", "1")
            return value_type(value)
        except (TypeError, ValueError) as error:
            raise ValueError(
                f"'{value:}' is not {value_type.__name__:} as {self.path:} needs.") from error

    def set_value(self, value):
        setattr(self.item, self.attribute, self.convert_value(value))


def iterate_items(cryspy_object) -> Iterator[Tuple[str, ItemN]]:
    """(path, item) of every ItemN of the object, also of the loops:
    'crystal.cell', 'crystal.atom_site[Fe]', ..."""
    if isinstance(cryspy_object, GlobalN):
        for obj in cryspy_object.items:
            yield from iterate_items(obj)
    elif isinstance(cryspy_object, DataN):
        data_name = cryspy_object.get_name()
        for block in cryspy_object.items:
            for path, item in iterate_items(block):
                yield f"{data_name:}.{path:}".strip("."), item
    elif isinstance(cryspy_object, LoopN):
        name = cryspy_object.get_name()
        for index, item in enumerate(cryspy_object.items):
            yield f"{name:}[{take_item_label(item, index):}]", item
    elif isinstance(cryspy_object, ItemN):
        yield cryspy_object.get_name(), cryspy_object


def iterate_blocks(cryspy_object) -> Iterator[Tuple[str, object]]:
    """(path, object) of the data blocks, loops and items of the object."""
    for path, item in iterate_items(cryspy_object):
        yield path, item
    if isinstance(cryspy_object, GlobalN):
        l_obj = cryspy_object.items
    else:
        l_obj = [cryspy_object]
    for obj in l_obj:
        if isinstance(obj, DataN):
            yield obj.get_name(), obj
            for block in obj.items:
                if isinstance(block, LoopN):
                    yield f"{obj.get_name():}.{block.get_name():}", block
        elif isinstance(obj, LoopN):
            yield obj.get_name(), obj


def take_attribute_records(cryspy_object, flag_defined_only: bool = True
                           ) -> List[AttributeRecord]:
    """All attributes of the object; with flag_defined_only only those
    which have a value."""
    l_record = []
    for path, item in iterate_items(cryspy_object):
        for attribute in item.ATTR_NAMES:
            if flag_defined_only and not item.is_attribute(attribute):
                continue
            l_record.append(AttributeRecord(f"{path:}.{attribute:}", item, attribute))
    return l_record


def find_attribute_records(cryspy_object, pattern: str = "",
                           flag_defined_only: bool = True) -> List[AttributeRecord]:
    """Attributes whose path matches the pattern as in
    find_parameter_records ('chi_type', '*.extinction.model', ...)."""
    l_regex = [form_path_regex(part) for part in split_pattern(pattern)]
    l_record = []
    for record in take_attribute_records(cryspy_object, flag_defined_only):
        path = record.path.lower()
        path_no_label = remove_labels(path)
        if any(regex.search(path) is not None or regex.search(path_no_label) is not None
               for regex in l_regex):
            l_record.append(record)
    return l_record


def take_attribute_record(cryspy_object, path: str,
                          flag_defined_only: bool = False) -> Optional[AttributeRecord]:
    path_sh = path.strip().lower()
    for record in take_attribute_records(cryspy_object, flag_defined_only):
        if record.path.lower() == path_sh:
            return record
    return None


def take_block(cryspy_object, path: str):
    """Data block, loop or item with this path (labels may be left out
    when there is only one such item), or None."""
    path_sh = path.strip().lower()
    l_data = cryspy_object.items if isinstance(cryspy_object, GlobalN) else [cryspy_object]
    # the data block 'crystal_fe' may be called by its data name 'fe'
    d_data_name = {data.get_name().lower(): str(data.data_name).lower()
                   for data in l_data if isinstance(data, DataN)}
    l_match = []
    for block_path, block in iterate_blocks(cryspy_object):
        block_path = block_path.lower()
        l_head_tail = block_path.split(".", 1)
        l_path = [block_path]
        if l_head_tail[0] in d_data_name:
            l_path.append(".".join([d_data_name[l_head_tail[0]]] + l_head_tail[1:]))
        if path_sh in l_path:
            return block
        if path_sh in [remove_labels(_) for _ in l_path] + [block_path.split(".", 1)[-1]]:
            l_match.append(block)
    if len(l_match) == 1:
        return l_match[0]
    return None


def find_item_classes(name: str) -> List[type]:
    """ItemN classes (also of loops) whose PREFIX or class name is the
    name, or contains it when no name is equal."""
    name_sh = re.sub(r"[\s\-\.]+", "_", name.strip().lower()).strip("_")
    l_class = take_item_classes()
    l_equal = [item_class for item_class in l_class
               if name_sh in (item_class.PREFIX.lower(), item_class.__name__.lower())]
    if len(l_equal) > 0:
        return l_equal
    return [item_class for item_class in l_class
            if name_sh in item_class.PREFIX.lower() or name_sh in item_class.__name__.lower()]


def explain_allowed_value(word: str) -> str:
    """'X is a value of category.attribute ...' when the word is an allowed
    value of some attribute (e.g. 'gauss' of extinction.model), else ''."""
    word_sh = word.strip().lower()
    if word_sh == "":
        return ""
    l_place = []
    for item_class in take_item_classes():
        for d_attribute in describe_attributes(item_class):
            allowed = d_attribute["allowed_values"]
            if allowed is None or len(allowed) > 30:
                continue
            if any(str(value).lower() == word_sh for value in allowed):
                l_place.append(
                    f"{item_class.PREFIX:}.{d_attribute['name']:} "
                    f"(allowed: {', '.join(str(_) for _ in allowed):})")
    l_place = list(dict.fromkeys(l_place))
    if len(l_place) == 0:
        return ""
    return (f"'{word:}' is a value of {'; '.join(l_place[:5]):}, not a parameter: "
            "set it with set_attribute.")


def describe_item_class(item_class: type) -> str:
    """Text description of a category: first line of the docstring and the
    attributes with types and allowed values."""
    s_doc = (item_class.__doc__ or "").strip()
    s_doc = s_doc.split("\n\n")[0].replace("\n", " ").strip()
    ls_out = [f"{item_class.PREFIX:} ({item_class.__name__:}): {s_doc:}".rstrip(": ")]
    for d_attribute in describe_attributes(item_class):
        l_flag = [d_attribute["type"]]
        if d_attribute["mandatory"]:
            l_flag.append("mandatory")
        if d_attribute["refinable"]:
            l_flag.append("refinable")
        s_line = f"  {d_attribute['cif_tag']:} ({', '.join(l_flag):})"
        allowed = d_attribute["allowed_values"]
        if allowed is not None:
            s_allowed = ", ".join(str(_) for _ in allowed[:30])
            if len(allowed) > 30:
                s_allowed += f", ... ({len(allowed):} values)"
            s_line += f": {s_allowed:}"
        ls_out.append(s_line)
    return "\n".join(ls_out)


def take_object_item_classes(cryspy_object) -> List[type]:
    l_class = []
    for _, block in iterate_blocks(cryspy_object):
        item_class = take_item_class(block)
        if item_class is not None and item_class not in l_class:
            l_class.append(item_class)
    return l_class
