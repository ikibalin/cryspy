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
import difflib
import re
from typing import List, Optional

import numpy

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN
from cryspy.B_parent_classes.cl_3_data import DataN
from cryspy.B_parent_classes.cl_4_global import GlobalN, find_changed_dict_name

L_LABEL_ATTRIBUTES = ("label", "atom_label", "id")

# words of a request (as a language model writes them) -> parts of the paths
D_PATTERN_ALIASES = {
    "position": ("fract_x", "fract_y", "fract_z"),
    "atomic_position": ("fract_x", "fract_y", "fract_z"),
    "atom_position": ("fract_x", "fract_y", "fract_z"),
    "coordinate": ("fract_x", "fract_y", "fract_z"),
    "atomic_coordinate": ("fract_x", "fract_y", "fract_z"),
    "fract": ("fract_x", "fract_y", "fract_z"),
    "fract_xyz": ("fract_x", "fract_y", "fract_z"),
    "xyz": ("fract_x", "fract_y", "fract_z"),
    "scale_factor": ("scale", ),
    "scale_parameter": ("scale", ),
    "polarisation": ("polarization", ),
    "beam_polarization": ("polarization", ),
    "flipper_efficiency": ("efficiency", ),
    "unit_cell": ("cell.length", "cell.angle"),
    "cell_parameter": ("cell.length", "cell.angle"),
    "lattice_parameter": ("cell.length", "cell.angle"),
    "adp": ("iso_or_equiv", "atom_site_aniso"),
    "thermal_parameter": ("iso_or_equiv", "atom_site_aniso"),
    "atomic_displacement": ("iso_or_equiv", "atom_site_aniso"),
    "b_iso": ("b_iso_or_equiv", ),
    "u_iso": ("u_iso_or_equiv", ),
    "occupancies": ("occupancy", ),
    "moment": ("atom_site_moment", ),
    "magnetic_moment": ("atom_site_moment", ),
    "susceptibility": ("atom_site_susceptibility", ),
}

# values of non-refinable attributes that a model writes as parameters:
# word of the pattern -> explanation for the error
D_VALUE_WORDS = {
    "ciso": "Ciso and Cani are values of atom_site_susceptibility.chi_type of each atom, not "
            "parameters: set them with set_chi_type. With Ciso only the scalar "
            "atom_site_susceptibility.chi_11 is refined, with Cani the components of the tensor "
            "(atom_site_susceptibility.chi_11, chi_22, chi_33, chi_12, chi_13, chi_23).",
}
D_VALUE_WORDS["cani"] = D_VALUE_WORDS["ciso"]


def explain_value_words(pattern: str) -> str:
    """Explanation when the pattern names a value of an attribute (e.g.
    'atom_site_susceptibility.Ciso'), else an empty string."""
    l_word = re.split(r"[^a-z0-9]+", pattern.lower())
    return " ".join(dict.fromkeys(D_VALUE_WORDS[word] for word in l_word if word in D_VALUE_WORDS))


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



def is_number(value) -> bool:
    return isinstance(value, (int, float, numpy.integer, numpy.floating)) and \
        not isinstance(value, bool)


def take_dict_paths(cryspy_object, dict_global: dict) -> dict:
    """Paths of the refined parameters by their ways in the dictionary.

    Parameters
    ----------
    cryspy_object : GlobalN
        Object with the crystals and the experiments.
    dict_global : dict
        Its dictionary (``cryspy_object.get_dictionary()``).

    Returns
    -------
    dict
        Way (block name, key, index) -> path
        ('crystal_fe.atom_site[Fe1].occupancy'), for the refinement log.
    """
    d_path = {}
    if not isinstance(cryspy_object, GlobalN):
        return d_path
    for data in cryspy_object.items:
        if not isinstance(data, DataN) or data.get_name() not in dict_global:
            continue
        data_name = data.get_name()
        for block in data.items:
            for record in form_block_records(block, data_name):
                # loops unknown to CrysPy keep their values as text
                if not record.flag_refinement or not is_number(record.value):
                    continue
                way = find_changed_dict_name(
                    data, data_name, dict_global[data_name], record.set_value,
                    record.value)
                if way is not None:
                    d_path[way] = record.path
    return d_path

def form_path_regex(pattern: str):
    """Compile a path pattern: '*' is any text, otherwise a substring.

    Square brackets are literal, so 'atom_site[Fe]' works as written.
    """
    pattern_sh = pattern.strip().lower()
    if "*" in pattern_sh or "?" in pattern_sh:
        s_regex = re.escape(pattern_sh).replace(r"\*", ".*").replace(r"\?", ".")
        return re.compile(f"^{s_regex:}$")
    return re.compile(re.escape(pattern_sh))


def remove_value_marks(path: str) -> str:
    """Path without the rcif marks of refinement written after a name:
    'cell.length_a()' or 'cell.length_a(3)' -> 'cell.length_a' (a model
    copies '()' from the rcif text into the path)."""
    # only after a name: the label 'O(1)' in 'atom_site[O(1)]' stays
    return re.sub(r"(?<=\w)\([^()]*\)(?=$|[.,;|\s])", "", path)


def remove_labels(path: str) -> str:
    """Path without the labels of loop items: 'b.atom_site[Fe].fract_x' ->
    'b.atom_site.fract_x' (a model often writes the path this way)."""
    return re.sub(r"\[[^\]]*\]", "", path)


def normalize_pattern(pattern: str) -> str:
    pattern_sh = remove_value_marks(pattern).strip().lower()
    pattern_sh = re.sub(r"[\s\-]+", "_", pattern_sh).strip("_")
    for prefix in ("all_", "the_"):
        if pattern_sh.startswith(prefix):
            pattern_sh = pattern_sh[len(prefix):]
    return pattern_sh


def take_alias_patterns(pattern: str) -> List[str]:
    """Patterns of the paths for a word of a request ('atomic positions',
    'scale factor'); the pattern itself when it is not an alias."""
    pattern_sh = normalize_pattern(pattern)
    for key in (pattern_sh, pattern_sh.rstrip("s")):
        if key in D_PATTERN_ALIASES:
            return list(D_PATTERN_ALIASES[key])
    return [pattern_sh]


def split_pattern(pattern: str) -> List[str]:
    """Several patterns given at once: 'fract_x, fract_y' or 'scale | b_iso'."""
    l_part = re.split(r"\s*(?:[,;|]|\s+or\s+|\s+and\s+)\s*", pattern.strip())
    l_pattern = []
    for part in l_part:
        if part.strip() == "":
            continue
        l_pattern.extend(take_alias_patterns(part))
    return l_pattern if len(l_pattern) > 0 else [""]


def find_parameter_records(
        cryspy_object, pattern: str = "",
        flag_refined_only: bool = False) -> List[ParameterRecord]:
    """Parameters whose path matches the pattern (case-insensitive).

    Without '*' the pattern is searched as a part of the path, so
    'fract_x' finds the x coordinate of every atom. With '*' the whole
    path must match: '*.atom_site[*].b_iso_or_equiv'. The path is also
    tried without the labels of loop items ('atom_site.fract_x').
    Several patterns are separated by ',' or '|'; words such as
    'atomic positions' or 'scale factor' are understood (D_PATTERN_ALIASES).
    """
    l_regex = [form_path_regex(part) for part in split_pattern(pattern)]
    l_record = []
    for record in take_parameter_records(cryspy_object):
        if flag_refined_only and not record.flag_refinement:
            continue
        path = record.path.lower()
        path_no_label = remove_labels(path)
        if any(regex.search(path) is not None
               or regex.search(path_no_label) is not None
               for regex in l_regex):
            l_record.append(record)
    return l_record


def take_parameter_names(cryspy_object) -> List[str]:
    """Paths without data block names and labels, each once:
    'cell.length_a', 'atom_site.fract_x', ..."""
    l_name = []
    for record in take_parameter_records(cryspy_object):
        name = remove_labels(record.path)
        name = name.split(".", 1)[1] if name.count(".") > 1 else name
        if name not in l_name:
            l_name.append(name)
    return l_name


def suggest_parameter_names(cryspy_object, pattern: str,
                            n_max: int = 15) -> List[str]:
    """Names of the parameters of the object close to the pattern; all
    names (at most n_max) when nothing is close."""
    l_name = take_parameter_names(cryspy_object)
    l_suggestion = []
    for part in split_pattern(pattern):
        l_word = [word for word in re.split(r"[^a-z0-9]+", part) if len(word) > 2]
        for name in l_name:
            if any(word in name.lower() for word in l_word):
                l_suggestion.append(name)
        l_suggestion.extend(difflib.get_close_matches(part, l_name, n=5, cutoff=0.5))
    l_suggestion = list(dict.fromkeys(l_suggestion))
    if len(l_suggestion) == 0:
        l_suggestion = l_name
    return l_suggestion[:n_max]


def take_data_name_aliases(cryspy_object) -> dict:
    """Name of each data block in the paths -> its data name: the block
    'diffrn_exp1' of the class Diffrn is often called 'exp1'."""
    l_data = cryspy_object.items if isinstance(cryspy_object, GlobalN) else [cryspy_object]
    return {data.get_name().lower(): str(data.data_name).lower()
            for data in l_data if isinstance(data, DataN)}


def form_path_aliases(path: str, d_data_name: dict) -> List[str]:
    """The path (lower case) and the same path with the data name of its
    data block: 'diffrn_exp1.setup.field' -> also 'exp1.setup.field'."""
    path = path.lower()
    l_head_tail = path.split(".", 1)
    l_path = [path]
    if l_head_tail[0] in d_data_name:
        l_path.append(".".join([d_data_name[l_head_tail[0]]] + l_head_tail[1:]))
    return l_path


def take_parameter_record(cryspy_object, path: str) -> Optional[ParameterRecord]:
    """Parameter with exactly this path (the data block may be given by its
    data name), or None."""
    path_sh = remove_value_marks(path).strip().lower()
    d_data_name = take_data_name_aliases(cryspy_object)
    for record in take_parameter_records(cryspy_object):
        if path_sh in form_path_aliases(record.path, d_data_name):
            return record
    return None
