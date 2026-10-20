"""Refinable parameters of a CrysPy object addressed by readable paths.

A path names one refinable parameter as text, so that a person or a
language model can refer to it without knowing the internal tuple names:

    ho2ti2o7.cell.length_a
    ho2ti2o7.atom_site[Ho].fract_x
    exp_1.pd_background[3].intensity

The first part is the data block name, then the item (loop items get the
label, atom_label or id in brackets, or their index when there is no
label), then the attribute name as in ATTR_REF.
"""
import re
from typing import List, Optional

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN
from cryspy.B_parent_classes.cl_3_data import DataN
from cryspy.B_parent_classes.cl_4_global import GlobalN

L_LABEL_ATTRIBUTES = ("label", "atom_label", "id")


class ParameterRecord:
    """One refinable parameter: its path and where it is stored."""

    def __init__(self, path: str, item: ItemN, attribute: str):
        self.path = path
        self.item = item
        self.attribute = attribute

    @property
    def value(self):
        return getattr(self.item, self.attribute)

    @property
    def sigma(self) -> float:
        return float(getattr(self.item, f"{self.attribute:}_sigma"))

    @property
    def flag_refinement(self) -> bool:
        return bool(getattr(self.item, f"{self.attribute:}_refinement"))

    @property
    def flag_constraint(self) -> bool:
        """True if the parameter is fixed by symmetry or other constraint."""
        return bool(getattr(self.item, f"{self.attribute:}_constraint"))

    def set_refinement(self, flag: bool) -> bool:
        """Switch refinement on or off; constrained parameters stay fixed.

        Returns True if the flag has the requested value afterwards.
        """
        if flag and self.flag_constraint:
            return False
        setattr(self.item, f"{self.attribute:}_refinement", bool(flag))
        return True

    def set_value(self, value: float):
        setattr(self.item, self.attribute, float(value))


def take_item_label(item: ItemN, index: int) -> str:
    for attribute in L_LABEL_ATTRIBUTES:
        if attribute in item.ATTR_NAMES and item.is_attribute(attribute):
            label = getattr(item, attribute)
            if label is not None:
                return str(label)
    return str(index)


def form_item_records(item: ItemN, prefix: str) -> List[ParameterRecord]:
    l_record = []
    for attribute in item.ATTR_REF:
        if not item.is_attribute(attribute):
            continue
        if getattr(item, attribute) is None:
            continue
        l_record.append(
            ParameterRecord(f"{prefix:}.{attribute:}", item, attribute))
    return l_record


def form_block_records(block, prefix: str) -> List[ParameterRecord]:
    """Records of an item or loop placed in a data block or global."""
    l_record = []
    if isinstance(block, LoopN):
        name = f"{prefix:}.{block.get_name():}".strip(".")
        for index, item in enumerate(block.items):
            label = take_item_label(item, index)
            l_record.extend(form_item_records(item, f"{name:}[{label:}]"))
    elif isinstance(block, ItemN):
        name = f"{prefix:}.{block.get_name():}".strip(".")
        l_record.extend(form_item_records(block, name))
    return l_record


def take_parameter_records(cryspy_object) -> List[ParameterRecord]:
    """All refinable parameters which have a value in the object."""
    l_record = []
    if isinstance(cryspy_object, (GlobalN, DataN)):
        if isinstance(cryspy_object, DataN):
            l_data = [cryspy_object]
        else:
            l_data = cryspy_object.items
        for data in l_data:
            if isinstance(data, DataN):
                data_name = data.get_name()
                for block in data.items:
                    l_record.extend(form_block_records(block, data_name))
            else:
                l_record.extend(form_block_records(data, ""))
    else:
        l_record.extend(form_block_records(cryspy_object, ""))
    return l_record


def form_path_regex(pattern: str):
    """Compile a path pattern: '*' is any text, otherwise a substring.

    Square brackets are literal, so 'atom_site[Fe]' works as written.
    """
    pattern_sh = pattern.strip().lower()
    if "*" in pattern_sh or "?" in pattern_sh:
        s_regex = re.escape(pattern_sh).replace(r"\*", ".*").replace(r"\?", ".")
        return re.compile(f"^{s_regex:}$")
    return re.compile(re.escape(pattern_sh))


def find_parameter_records(
        cryspy_object, pattern: str = "",
        flag_refined_only: bool = False) -> List[ParameterRecord]:
    """Parameters whose path matches the pattern (case-insensitive).

    Without '*' the pattern is searched as a part of the path, so
    'fract_x' finds the x coordinate of every atom. With '*' the whole
    path must match: '*.atom_site[*].b_iso_or_equiv'.
    """
    regex = form_path_regex(pattern)
    l_record = []
    for record in take_parameter_records(cryspy_object):
        if flag_refined_only and not record.flag_refinement:
            continue
        if regex.search(record.path.lower()) is not None:
            l_record.append(record)
    return l_record


def take_parameter_record(cryspy_object, path: str) -> Optional[ParameterRecord]:
    """Parameter with exactly this path, or None."""
    path_sh = path.strip().lower()
    for record in take_parameter_records(cryspy_object):
        if record.path.lower() == path_sh:
            return record
    return None
