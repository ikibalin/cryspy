"""Parent class DataN."""
import logging
import os
import numpy
import os.path
from warnings import warn
from typing import Union, NoReturn
from pycifstar import Global, to_global

from cryspy.A_functions_base.function_1_markdown import md_to_html
from cryspy.A_functions_base.linear_constraints import parse_soft_constraint
from cryspy.A_functions_base.function_1_objects import \
    get_functions_of_objet, get_table_html_for_variables

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN
from cryspy.exceptions import CrysPyWarning, RcifFileNotFoundError, RcifFormatError
from cryspy.B_parent_classes.cl_3_data import DataN

logger = logging.getLogger(__name__)



def find_changed_dict_entry(block, block_name: str, d_before: dict, set_value, value: float,
                            flag_only_refined: bool = True):
    """Way (block name, key, index) of a parameter of the data block in its
    dictionary d_before together with the linear map from the value of the
    object to the value of the dictionary, None if it is not found.

    set_value(new_value) changes the parameter a little, and the entry of
    the dictionary of the block that changes with it is the way; the value
    is restored afterwards. The dictionary keeps some parameters in other
    units (b_iso = 8 pi^2 u_iso, cell angles in radians):
    ``value_dict = slope * value + offset``.

    Returns
    -------
    (way, slope, offset) or None
    """
    # both directions: a value at its limit (occupancy 1) cannot grow
    for sign in (1., -1.):
        delta = sign * 0.01 * (abs(value) + 1.)
        try:
            set_value(value + delta)
            d_after = block.get_dictionary()
        finally:
            set_value(value)
        # the refined entries first; a parameter which is not refined has
        # the flag False
        l_pass = (True, ) if flag_only_refined else (True, False)
        for flag_refined in l_pass:
            for key, flags in d_before.items():
                if not key.startswith("flags_") or key[len("flags_"):] not in d_before:
                    continue
                key_value = key[len("flags_"):]
                try:
                    array_before = numpy.asarray(d_before[key_value], dtype=float)
                    array_after = numpy.asarray(d_after[key_value], dtype=float)
                    flag_changed = array_before != array_after
                    if flag_refined:
                        flag_changed = numpy.logical_and(
                            numpy.asarray(flags, dtype=bool), flag_changed)
                except (KeyError, TypeError, ValueError):
                    continue
                l_index = numpy.argwhere(flag_changed)
                if len(l_index) > 0:
                    index = tuple(int(hh) for hh in l_index[0])
                    value_before = float(array_before[index])
                    # rounding errors of the finite difference are removed:
                    # the maps are linear (1, 8 pi^2, pi/180)
                    slope = float(f"{(float(array_after[index]) - value_before) / delta:.8g}")
                    offset = value_before - slope * float(value)
                    if abs(offset) <= 1e-7 * (abs(value_before) + 1.):
                        offset = 0.
                    return (block_name, key_value, index), slope, offset
    return None


def find_changed_dict_name(block, block_name: str, d_before: dict, set_value, value: float):
    """Way (block name, key, index) of a refined parameter of the data
    block in its dictionary d_before, None if it is not found.

    set_value(new_value) changes the parameter a little, and the refined
    entry of the dictionary of the block that changes with it is the way;
    the value is restored afterwards.
    """
    entry = find_changed_dict_entry(block, block_name, d_before, set_value, value)
    if entry is None:
        return None
    return entry[0]


def take_marked_variable_names(obj) -> list:
    """Names (as get_variable_names gives them) and marks of all marked
    parameters of the object, refined or not.

    Returns
    -------
    list
        ``[(variable_name, mark, flag_refinement), ...]``.
    """
    l_out = []
    if isinstance(obj, ItemN):
        for name, name_flag in zip(obj.ATTR_REF, obj.ATTR_REF_FLAG):
            mark = getattr(obj, f"{name:}_mark", "")
            if isinstance(mark, str) and mark != "":
                l_out.append((((obj.PREFIX, None), (name, None)), mark,
                              bool(getattr(obj, name_flag, False))))
    elif isinstance(obj, LoopN):
        loop_name = obj.loop_name
        if isinstance(loop_name, str):
            loop_name = loop_name.lower()
        for ind, item in enumerate(obj.items):
            for name, mark, flag in take_marked_variable_names(item):
                l_out.append((((item.PREFIX, loop_name), (name[1][0], ind)), mark, flag))
    elif isinstance(obj, (DataN, GlobalN)):
        if isinstance(obj, DataN):
            name = obj.data_name
        else:
            name = obj.global_name
        if isinstance(name, str):
            name = name.lower()
        for item in obj.items:
            for var_name, mark, flag in take_marked_variable_names(item):
                l_out.append((((obj.PREFIX, name), ) + var_name, mark, flag))
    return l_out


def take_active_marks(l_mark_flag: list, l_equation_labels: list) -> set:
    """Labels of marks which are refined.

    A mark is active when a parameter with it is refined, or when it enters
    a hard constraint together with an active mark (x+y=1 with only x
    refined makes y refined too).

    Parameters
    ----------
    l_mark_flag : list
        ``[(mark, flag_refinement), ...]`` of the marked parameters.
    l_equation_labels : list
        Labels of each hard constraint.
    """
    active = set(mark.rstrip("-+") for mark, flag in l_mark_flag if flag)
    flag_changed = True
    while flag_changed:
        flag_changed = False
        for l_label in l_equation_labels:
            if active.intersection(l_label) and not active.issuperset(l_label):
                active.update(l_label)
                flag_changed = True
    return active


class MarkedParameter:
    """Marked parameter of an object: its name, mark, refinement and the
    entry (way, slope, offset) in the dictionary, None if it does not
    change the dictionary."""

    def __init__(self, variable_name: tuple, mark: str, flag_refinement: bool,
                 entry, value):
        self.variable_name = variable_name
        self.mark = mark
        self.flag_refinement = flag_refinement
        self.entry = entry
        self.value = value

    @property
    def label(self) -> str:
        return self.mark.rstrip("-+")

    @property
    def sign(self) -> float:
        return -1. if self.mark.endswith("-") else 1.


def form_dict_equation(l_coeff_marked: list, const: float):
    """Equation ``sum coeff_i x_i = const`` on the values of the object
    written for the values of the dictionary
    ``p_i = slope_i x_i + offset_i``.

    Returns
    -------
    ([(coeff, way), ...], const)
    """
    l_term = []
    const_dict = float(const)
    for coeff, marked in l_coeff_marked:
        way, slope, offset = marked.entry
        l_term.append((coeff / slope, way))
        const_dict += coeff * offset / slope
    return l_term, const_dict


def form_variable_path(variable_name: tuple) -> str:
    """Readable name of a parameter: crystal_pbcn.atom_site[1].u_iso_or_equiv."""
    l_part = []
    for hh in variable_name[1:-1]:
        name = "_".join(str(h) for h in hh if h not in (None, ""))
        if name != "":
            l_part.append(name)
    attribute, index = variable_name[-1][0], variable_name[-1][1]
    if index is None:
        l_part.append(attribute)
    else:
        l_part[-1] = f"{l_part[-1]:}[{index:}]"
        l_part.append(attribute)
    return ".".join(l_part)


class GlobalN(object):
    """Data container of data bloks, loops and items."""

    global_name = None  # name after 'global_' in the file, if any

    def __repr__(self):
        """
        Magic method print() is redefined.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        ls_out = [f"# Object '{self.get_name():}':"]
        for item in self.items:
            if isinstance(item, LoopN):
                ls_out.append(f"{4*' ':}.{item.get_name():} (loop)")
            else:
                ls_out.append(f"{4*' ':}.{item.get_name():}")
            if isinstance(item, DataN):
                for i_d in item.items:
                    if isinstance(i_d, ItemN):
                        ls_out.append(f"{8*' ':}.{i_d.get_name():}")
                    else:
                        ls_out.append(
                            f"{8*' ':}.{i_d.get_name():} (loop)")

        method = self.methods_html()
        if method != "":
            ls_out.append(f"\n# Methods:\n{method:}\n")                
        return "\n".join(ls_out)

    def _repr_html_(self):
        """Representation in HTML format."""
        ls_html = [f"<h1>Object '{self.get_name():}'</h1>"]
        ls_html.append(self.attributes_to_html())

        ls_html.append(get_table_html_for_variables(self))

        report = self.report_html()
        if report != "":
            ls_html.append(f"<h1>Description </h1> {report:}")

        ls_html.append("<h1>Classes and methods</h1>")
        try:
            names = sorted([obj.__name__ for obj in self.CLASSES_MANDATORY])
            if len(names) != 0:
                ls_html.append("<b>Mandatory classes: </b>")
                ls_html.append(f"{', '.join(names):}.<br>")
        except AttributeError:
            pass
        try:
            names = sorted([obj.__name__ for obj in self.CLASSES_OPTIONAL])
            if len(names) != 0:
                ls_html.append("<b>Optional classes: </b>")
                ls_html.append(f"{', '.join(names):}.<br>")
        except AttributeError:
            pass

        method = self.methods_html()
        if method != "":
            ls_html.append(f"<b>Methods: </b> {method:}")
        return " ".join(ls_html)

    def attributes_to_html(self) -> str:
        """Representation of defined parameters in HTML format.
        """
        ls_html = ["<table>"]
        ls_html.append("<tr><th>Attribute</th><th> Note </th></tr>")
        items_sorted = sorted(self.items, key=lambda item: item.get_name())
        for item in items_sorted:
            item_type = item.__doc__.strip().split("\n")[0]
            ls_html.append(f"<tr><td>.<b>{item.get_name():}</b></td>\
<td><b>{item_type:}</b></td></tr>")
            if isinstance(item, DataN):
                items_s_d = sorted(item.items, key=lambda x: x.get_name())
                for item_d in items_s_d:
                    item_type_d = item_d.__doc__.strip().split("\n")[0]
                    ls_html.append(f"<tr><td>.{item.get_name():}.<b>{item_d.get_name():}</b></td>\
<td>{item_type_d:}</td></tr>")
        ls_html.append("</table>")
        return " ".join(ls_html)

    def methods_html(self):
        ls_html = [f".{func_name}" for func_name in
                   get_functions_of_objet(self)]
        return ", ".join(ls_html)+"."

    def __str__(self):
        """
        Magic method str() is redefined.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        return self.to_cif()

    def take_item_by_name(self, name: str):
        """Data block (or item, loop) of the container by its name
        (``rcif.crystal_fe``); None when there is none."""
        name = name.lower()
        for item in self.__dict__.get("items", ()):
            if name == item.get_name():
                return item
        return None

    def __getattr__(self, name):
        """The data blocks are reached by their names (``rcif.crystal_fe``);
        an unknown name raises AttributeError (``rcif.is_attribute(name)``
        or ``take_item_by_name`` check for a block that may be absent). Only
        names that are not attributes of the object come here; the blocks
        are read from __dict__, so an object without them (made by copy or
        pickle) never comes back here."""
        if name.startswith("__") and name.endswith("__"):
            raise AttributeError(name)
        item = self.take_item_by_name(name)
        if item is None:
            raise AttributeError(
                f"There is no data block or item '{name:}'")
        return item

    def __setattr__(self, name, value) -> NoReturn:
        """
        Rules to set attribute.

        Parameters
        ----------
        name : TYPE
            DESCRIPTION.
        value : TYPE
            DESCRIPTION.

        Returns
        -------
        NoReturn
            DESCRIPTION.

        """
        flag_items, flag_direct = False, True
        if name == "global_name":
            flag_direct = False
            val_new = str(value).strip()
        elif name == "items":
            flag_items = True
            self.add_items(value)
        else:
            cls_value = type(value)
            if cls_value in self.CLASSES:
                l_name = [item.get_name() for item in self.items]
                name_new = value.get_name()
                if name_new in l_name:
                    self.items.pop(l_name.index(name))
                self.items.append(value)
                flag_items, flag_direct = True, False
                if name_new != name:
                    warn(f"Access to variable by '{name_new:}'.", CrysPyWarning)

        if flag_items:
            pass
        elif flag_direct:
            self.__dict__[name] = value
        else:
            self.__dict__[name] = val_new

    def is_attribute(self, name):
        """Temporary construction.

        Better to use:

        try:
            obj = self.attribute_name
        except AttributeError as e:
            obj = ...
        """
        for item in self.items:
            if name.lower() == item.get_name():
                return True
        return False

    def add_items(self, items: list):
        """Add items."""
        items = [hh for hh in items if hh is not None]
        l_name = [item.get_name() for item in items]
        s_name = set(l_name)
        if len(s_name) != len(l_name):
            warn("Double items were given.", CrysPyWarning)
            items_unique = [items[l_name.index(name)] for name in s_name]
        else:
            items_unique = items
        l_ind_del = []
        for ind_item, item in enumerate(self.items):
            if item.get_name() in s_name:
                l_ind_del.append(ind_item)
        l_ind_del.reverse()
        for ind in l_ind_del:
            self.items.pop(ind)
        for item in items_unique:
            if isinstance(item, self.CLASSES):
                self.items.append(item)
            elif type(self) is GlobalN:
                if issubclass(type(item), (ItemN, LoopN, DataN)):
                    self.CLASSES = self.CLASSES + (type(item), )
                    self.CLASSES_OPTIONAL = self.CLASSES_OPTIONAL + (type(item), )
                    self.items.append(item)

    @classmethod
    def make_container(cls, cls_mandatory, cls_optional, prefix):
        """Create GlobalN object as a container for items."""
        if cls is not GlobalN:
            warn("The method 'make_container' is used only for GlobalN class.")
            return
        obj = cls()
        obj.__dict__["CLASSES_MANDATORY"] = cls_mandatory
        obj.__dict__["CLASSES_OPTIONAL"] = cls_optional
        obj.__dict__["CLASSES"] = cls_mandatory+cls_optional
        obj.__dict__["PREFIX"] = prefix
        obj.__dict__["D_DEFAULT"] = {}
        obj.__dict__["items"] = []
        obj.__dict__["data_name"] = ""
        return obj

    def __getitem__(self, name: Union[int, str]):
        """
        Get item by index or predefined index.

        Parameters
        ----------
        name : TYPE
            DESCRIPTION.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        if isinstance(name, int):
            return self.items[name]
        elif isinstance(name, str):
            for item in self.items:
                if name.lower() == item.get_name():
                    return item
        return None

    def get_name(self) -> str:
        """Name of object."""
        name = self.PREFIX
        global_name = self.global_name
        if global_name is not None:
            name = f"{name:}_{global_name:}"
        return name.lower()

    def get_variable_names(self) -> list:
        """
        Get names of variable as a list.

        (((#prefix, #NAME), (#prefix, #NAME), (#attribute, #index))

        Returns
        -------
        list
            List of names of variable.

        """
        prefix = self.PREFIX
        global_name = self.global_name
        if isinstance(global_name, str):
            global_name = global_name.lower()
        l_var = []
        for item in self.items:
            l_var.extend(item.get_variable_names())
        l_var_out = [((prefix, global_name), ) + var for var in l_var]
        return l_var_out

    def is_variables(self) -> bool:
        """Define is there variables or not."""
        flag = False
        for item in self.items:
            if item.is_variables():
                flag = True
                break
        return flag

    def get_variable_by_name(self, name: tuple) -> Union[float, int, str]:
        """
        Get variable given by name.

        Parameters
        ----------
        name : tuple
            (((#prefix, #data_name), (#prefix, #loop_name),
               (#attribute, #index_item))

        Returns
        -------
        Union[float, int, str]
            DESCRIPTION.

        """
        prefix = self.PREFIX
        global_name = self.global_name
        if isinstance(global_name, str):
            global_name = global_name.lower()

        prefix_g, prefix_n = name[0], name[1]
        if isinstance(prefix_g[1], str):
            prefix_g_2 = prefix_g[1].lower()
        else:
            prefix_g_2 = prefix_g[1]

        if (prefix_g[0], prefix_g_2) != (prefix, global_name):
            return None

        name_sh = tuple(name[1:])
        for item in self.items:
            if isinstance(item, ItemN):
                prefix = item.PREFIX
            elif isinstance(item, LoopN):
                prefix = item.ITEM_CLASS.PREFIX
            elif isinstance(item, DataN):
                prefix = item.PREFIX
            else:
                raise AttributeError(
                    f"Unknown type object '{type(item).__name__:}'")
            if prefix == prefix_n[0]:
                res = item.get_variable_by_name(name_sh)
                if res is not None:
                    return res
        return None

    def set_variable_by_name(self, name: tuple, value) -> NoReturn:
        """
        Set value to variable given by name.

        Parameters
        ----------
        name : tuple
            DESCRIPTION.
        value : TYPE
            DESCRIPTION.

        Returns
        -------
        NoReturn
            DESCRIPTION.

        """
        prefix = self.PREFIX
        global_name = self.global_name

        prefix_g, prefix_n = name[0], name[1]

        if prefix_g != (prefix, global_name):
            return

        name_sh = tuple(name[1:])
        for item in self.items:
            if isinstance(item, ItemN):
                prefix = item.PREFIX
            elif isinstance(item, LoopN):
                prefix = item.ITEM_CLASS.PREFIX
            elif isinstance(item, DataN):
                prefix = item.PREFIX
            else:
                raise AttributeError(
                    f"Unknown type object '{type(item).__name__:}'")
            if prefix == prefix_n[0]:
                item.set_variable_by_name(name_sh, value)

    def is_defined(self) -> bool:
        """
        If all mandatory attributes is defined.

        Returns
        -------
        bool
            DESCRIPTION.

        """
        flag = True
        for item in self.items:
            if not(item.is_defined()):
                flag = False
                if isinstance(item, (ItemN, DataN)):
                    warn(f"Item {item.PREFIX:} is not fully described.",
                         CrysPyWarning)
                    break
                elif isinstance(item, LoopN):
                    warn(f"Item {item.ITEM_CLASS.PREFIX:} is not fully described.",
                         CrysPyWarning)
                    break
        if flag:
            cls_items = [type(item) for item in self.items]
            for cls_mand in self.CLASSES_MANDATORY:
                if not(cls_mand in cls_items):
                    flag = False
                    warn(f"The object of {cls_mand.__name__:} is not defined.",
                         CrysPyWarning)
                    break
        return flag

    def form_object(self):
        """Form object."""
        pass

    def to_cif(self, separator="_") -> str:
        """
        Print information about object in string in STAR format.

        Arguments
        ---------
            prefix: prefix in front of label of attribute
            separator: separator between prefix and attribute ("_" or ".")
            flag: for undefined attribute "." will be printed
            flag_minimal if it's True the minimal set of object will be printed

        Returns
        -------
            A string in STAR/CIF format
        """
        ls_out = []
        if self.global_name is None:
            ls_out.append("global_\n")
        else:
            ls_out.append(f"global_{self.global_name:}\n")
        l_item = self.items
        l_s_itemn = [item.to_cif(separator=separator)+"\n"
                     for item in l_item if isinstance(item, ItemN)]
        l_s_loopn = [item.to_cif(separator=separator)+"\n"
                     for item in l_item if isinstance(item, LoopN)]
        l_s_datan = [item.to_cif(separator=separator)+"\n"
                     for item in l_item if isinstance(item, DataN)]
        if l_s_loopn != []:
            n_max_loop = max([len(_) for _ in l_s_loopn])
            if n_max_loop < 1000:
                n_max_loop = 1000
        else:
            n_max_loop = 10000
        l_n_max_item = [len(_) for _ in l_s_itemn]

        ls_out.extend([_1 for _1, _2 in zip(l_s_itemn, l_n_max_item)
                       if _2 <= n_max_loop])
        ls_out.extend([_ for _ in l_s_loopn])
        ls_out.extend([_1 for _1, _2 in zip(l_s_itemn, l_n_max_item)
                       if _2 > n_max_loop])
        ls_out.extend(sorted(l_s_datan, key=len))
        return "\n".join(ls_out)

    @classmethod
    def from_cif(cls, string: str):
        """Generate object from string of CIF format."""
        cif_global = Global()
        flag = cif_global.take_from_string(string)

        cif_items = cif_global.items
        cif_loops = cif_global.loops
        cif_datas = cif_global.datas

        items = []
        flag = True
        n_mandatory = len(cls.CLASSES_MANDATORY)
        for i_cls, cls_ in enumerate(cls.CLASSES):
            flag = i_cls >= n_mandatory
            if issubclass(cls_, ItemN):
                prefix_cls = cls_.PREFIX
                if cif_items.is_prefix(prefix_cls):
                    cif_items_prefix = cif_items[prefix_cls]
                    cif_string = str(cif_items_prefix)
                    obj_prefix = cls_.from_cif(cif_string)
                    if obj_prefix is not None:
                        items.append(obj_prefix)
                        flag = True
            elif issubclass(cls_, LoopN):
                prefix_cls = cls_.ITEM_CLASS.PREFIX
                for cif_loop in cif_loops:
                    if cif_loop.is_prefix("_"+prefix_cls):
                        cif_string = str(cif_loop)
                        obj_prefix = cls_.from_cif(cif_string)
                        if obj_prefix is not None:
                            items.append(obj_prefix)
                            flag = True
            elif issubclass(cls_, DataN):
                prefix_cls = cls_.PREFIX
                for cif_data in cif_datas:
                    cif_string = str(cif_data)
                    obj_prefix = cls_.from_cif(cif_string)
                    if obj_prefix is not None:
                        items.append(obj_prefix)
                        flag = True
            if not(flag):
                warn(f"Mandatory class: '{cls_.__name__:}' is not given.",
                     CrysPyWarning)
                break

        if not(flag):
            return None

        global_name = cif_global.name
        obj = cls(global_name=global_name, items=items)
        obj.form_object()
        return obj

    @classmethod
    def from_cif_file(cls, f_name: str):
        """Read from cif file."""
        if not(os.path.isfile(f_name)):
            raise RcifFileNotFoundError(f"File {f_name:} is not found.")
            return None
        str_from_cif = str(to_global(f_name))
        obj = cls.from_cif(str_from_cif)
        obj.file_input = f_name
        return obj

    def report(self):
        return ""

    def report_html(self):
        return md_to_html(self.report())

    def plots(self):
        l_res = []
        for item in self.items:
            for plot in item.plots():
                if plot is not None:
                    l_res.append(plot)
        return l_res

    def fix_variables(self):
        """Fix variables."""
        for item in self.items:
            item.fix_variables()

    def set_variable(self, name: str, index=None):
        """Set refinement for variable given by name.
        
        Index parameters is used only for objects given as a matrix.
        """
        name_sh = name.strip(".").lower()
        l_name = name_sh.split(".")
        name_1 = l_name[0]
        for item in self.items:
            if name_1 == item.get_name(): 
                if len(l_name) == 0:
                    pass
                else:
                    item.set_variable(".".join(l_name[1:]), index=index)

    def get_dictionary(self):
        self.refine_constrained_parameters()
        dict_out = {}
        for item in self.items:
            if "get_dictionary" in dir(item):
                name = item.get_name()
                dict_out[name] = item.get_dictionary()

        variable_names = self.get_variable_names()
        dict_names = transfer_cif_names_to_dict_names(variable_names)

        # the table of transfer_cif_names_to_dict_names does not know every
        # parameter (occupancy, cell lengths, ...) and the dictionary keeps
        # some of them in other units (b_iso = 8 pi^2 u_iso, angles in
        # radians): marked parameters, refined or not, are found in the
        # dictionary of their data block together with the units
        l_marked = []
        for variable_name, mark, flag_refinement in take_marked_variable_names(self):
            entry = self.find_dict_entry(variable_name, dict_out, flag_only_refined=flag_refinement)
            l_marked.append(MarkedParameter(
                variable_name, mark, flag_refinement, entry,
                self.get_variable_by_name(variable_name)))
        d_marked = {marked.variable_name: marked for marked in l_marked}
        for ind, variable_name in enumerate(variable_names):
            marked = d_marked.get(variable_name)
            if marked is not None and marked.entry is not None:
                dict_names[ind] = marked.entry[0]

        # equal and opposite values: marks a, a-
        l_refined = [marked for marked in l_marked if marked.flag_refinement]
        linear_constraints, linear_equations = [], []
        for ind_1, marked_1 in enumerate(l_refined):
            for marked_2 in l_refined[ind_1+1:]:
                if marked_1.label != marked_2.label:
                    continue
                coeff_1 = -1. if marked_1.mark.endswith("-") else 1.
                coeff_2 = 1. if marked_2.mark.endswith("-") else -1.
                l_term, const = form_dict_equation(
                    [(coeff_1, marked_1), (coeff_2, marked_2)], 0.)
                if abs(const) < 1e-12:
                    linear_constraints.append(l_term)
                else:
                    linear_equations.append((l_term, const))
        if len(linear_constraints) > 0:
            dict_out["linear_constraints"] = linear_constraints
        linear_equations_loop, soft_constraints = self.get_constraints(l_marked)
        linear_equations.extend(linear_equations_loop)
        if len(linear_equations) > 0:
            dict_out["linear_equations"] = linear_equations
        if len(soft_constraints) > 0:
            dict_out["soft_constraints"] = soft_constraints
        # restraints of distances and angles (restr_distance, ... of crystals)
        restraints = [restraint for item in self.items
                      if hasattr(type(item), "get_restraints")
                      for restraint in item.get_restraints()]
        if len(restraints) > 0:
            dict_out["restraints"] = restraints
        try:
            expression = self.punishment.function
            d_mark_dict_name = {}
            for rcif_name, dict_name in zip(variable_names, dict_names):
                rcif_name_mark = rcif_name[:-1]+ ((rcif_name[-1][0]+"_mark", rcif_name[-1][1]),)
                mark = self.get_variable_by_name(rcif_name_mark)
                if mark != "":
                    d_mark_dict_name[mark] = dict_name
            dict_out["punishment_function"] = (expression, d_mark_dict_name)
        except Exception:
            pass
        return dict_out

    def find_dict_entry(self, variable_name: tuple, dict_global: dict,
                        flag_only_refined: bool = True):
        """Way (block name, key, index) of a parameter in the dictionary of
        the object with the linear map of its units, None if it is not
        found (see find_changed_dict_entry)."""
        block_name = "_".join(str(hh) for hh in variable_name[1] if hh not in (None, ""))
        if block_name not in dict_global:
            return None
        l_block = [item for item in self.items if item.get_name() == block_name]
        if len(l_block) != 1:
            return None
        value = self.get_variable_by_name(variable_name)
        try:
            value = float(value)
        except (TypeError, ValueError):
            return None
        return find_changed_dict_entry(
            l_block[0], block_name, dict_global[block_name],
            lambda new_value: self.set_variable_by_name(variable_name, new_value), value,
            flag_only_refined=flag_only_refined)

    def find_dict_name(self, variable_name: tuple, dict_global: dict):
        """Way (block name, key, index) of a refined parameter in the
        dictionary of the object, None if it is not found.

        The value is changed a little, and the refined entry of the
        dictionary of its data block that changes with it is the way.
        """
        entry = self.find_dict_entry(variable_name, dict_global)
        if entry is None:
            return None
        return entry[0]

    def take_hard_constraint_labels(self) -> list:
        """Labels of the marks of each hard constraint of the constraint
        loops (wrong expressions are skipped here and reported by
        get_constraints)."""
        l_labels = []
        for item in self.items:
            if not callable(type(item).__dict__.get("get_linear_equations")):
                continue
            for constraint in item.items:
                try:
                    if not constraint.is_soft():
                        l_labels.append(list(constraint.get_linear_equation()[0]))
                except ValueError:
                    continue
        return l_labels

    def refine_constrained_parameters(self) -> list:
        """Switch on the refinement of the marked parameters tied to refined
        ones by hard constraints (x+y=1 with only x refined) or by equal
        marks (a, a-).

        Returns
        -------
        list
            Names of the parameters whose refinement was switched on.
        """
        l_marked = take_marked_variable_names(self)
        if len(l_marked) == 0:
            return []
        active = take_active_marks([(mark, flag) for name, mark, flag in l_marked],
                                   self.take_hard_constraint_labels())
        l_name = []
        for variable_name, mark, flag_refinement in l_marked:
            if flag_refinement or mark.rstrip("-+") not in active:
                continue
            name_flag = variable_name[:-1] + (
                (f"{variable_name[-1][0]:}_refinement", variable_name[-1][1]), )
            self.set_variable_by_name(name_flag, True)
            l_name.append(variable_name)
            logger.info(f"{form_variable_path(variable_name):} [{mark:}] is refined "
                        "because of a constraint.")
        return l_name

    def get_constraints(self, l_marked: list):
        """Constraints of the constraint loops.

        Parameters
        ----------
        l_marked : list of MarkedParameter
            All marked parameters, refined or not.

        Returns
        -------
        linear_equations : list
            Hard constraints (x+y=1) as ``[([(coeff, dict_name), ...], const), ...]``
            in the units of the dictionary.
        soft_constraints : list
            Soft constraints as
            ``[(expression, sigma, {mark: (sign, dict_name, slope, offset)}), ...]``:
            the value of the mark is ``sign*(value_dict - offset)/slope``.

        A mark ``a`` names the first refined value marked ``a``; when only
        ``a-`` is written, the value is the opposite of ``a``. Constraints
        without refined values are skipped.
        """
        d_mark = {}
        for marked in l_marked:
            label = marked.label
            if label not in d_mark:
                d_mark[label] = marked
                continue
            old = d_mark[label]
            if (marked.flag_refinement, not marked.mark.endswith("-")) > \
                    (old.flag_refinement, not old.mark.endswith("-")):
                d_mark[label] = marked

        def check_marks(constraint, l_label):
            l_unknown = [label for label in l_label if label not in d_mark]
            if len(l_unknown) > 0:
                raise RcifFormatError(
                    f"_constraint_expression '{constraint.expression:}': "
                    f"no value is marked {', '.join(l_unknown):} "
                    "(write the mark after the parentheses, e.g. 0.5()x)")
            for label in l_label:
                if d_mark[label].entry is None:
                    raise RcifFormatError(
                        f"_constraint_expression '{constraint.expression:}': "
                        f"the value marked {label:} "
                        f"({form_variable_path(d_mark[label].variable_name):}) "
                        "does not change the calculation. Is it used by the model "
                        "(for example _atom_site_adp_type Uiso for u_iso_or_equiv)?")

        linear_equations, soft_constraints = [], []
        for item in self.items:
            if not callable(type(item).__dict__.get("get_linear_equations")):
                continue
            for constraint in item.items:
                try:
                    if constraint.is_soft():
                        if not constraint.sigma > 0.:
                            raise ValueError("sigma must be positive")
                        l_label = parse_soft_constraint(constraint.expression)[1]
                    else:
                        d_coeff, const = constraint.get_linear_equation()
                        l_label = list(d_coeff)
                except ValueError as error:
                    raise RcifFormatError(
                        f"_constraint_expression '{constraint.expression:}': {error:}") from None
                l_unknown = [label for label in l_label if label not in d_mark]
                if len(l_unknown) == 0 and not any(
                        d_mark[label].flag_refinement for label in l_label):
                    # nothing is refined: the constraint does not act
                    logger.info(f"Constraint '{constraint.expression:}' is not used: "
                                "none of its values is refined.")
                    continue
                check_marks(constraint, l_label)
                if constraint.is_soft():
                    soft_constraints.append((
                        constraint.expression, float(constraint.sigma),
                        {label: (d_mark[label].sign, ) + tuple(d_mark[label].entry)
                         for label in l_label}))
                else:
                    linear_equations.append(form_dict_equation(
                        [(coeff * d_mark[label].sign, d_mark[label])
                         for label, coeff in d_coeff.items()], const))
        return linear_equations, soft_constraints

    def take_parameters_from_dictionary(self, dict_global, l_parameter_name: list=None, l_sigma: list=None):
        keys = dict_global.keys()

        if (l_parameter_name is None) or (l_sigma is None):
            l_parameter_name = []
            l_sigma = []

        for item in self.items:
            item_name = item.get_name().lower()
            if (item_name in keys) and ("take_parameters_from_dictionary" in dir(item)):
                l_parameter_name_data, l_sigma_data = [], []
                for parameter_name, sigma in zip(l_parameter_name, l_sigma):
                    if parameter_name[0].lower() == item_name.lower():
                        l_parameter_name_data.append(parameter_name[1:])
                        l_sigma_data.append(sigma)
                dict_data = dict_global[item_name]
                item.take_parameters_from_dictionary(dict_data, l_parameter_name=l_parameter_name_data, l_sigma=l_sigma_data)


def transfer_cif_names_to_dict_names(variable_names):
    dict_names = []
    for cif_name in variable_names:
        data_block = f"{cif_name[1][0]:}_{cif_name[1][1]:}"
        name = cif_name[-1][0]
        index = cif_name[-1][1]
        if index is None:
            index = (0,)
        if name == "polarization":
            name = "beam_polarization"
        elif name == "scale":
            name = "phase_scale"
        elif name == "efficiency":
            name = "flipper_efficiency"
        elif name == "radius":
            name = "extinction_radius"
        elif name == "mosaicity":
            name = "extinction_mosaicity"
        elif name == "chi_11":
            name = "atom_para_susceptibility"
            index = (0, index)
        elif name == "chi_22":
            name = "atom_para_susceptibility"
            index = (1, index)
        elif name == "chi_33":
            name = "atom_para_susceptibility"
            index = (2, index)
        elif name == "chi_12":
            name = "atom_para_susceptibility"
            index = (3, index)
        elif name == "chi_13":
            name = "atom_para_susceptibility"
            index = (4, index)
        elif name == "chi_23":
            name = "atom_para_susceptibility"
            index = (5, index)
        elif name == "u":
            name = "resolution_parameters"
            index = (0, )
        elif name == "v":
            name = "resolution_parameters"
            index = (1, )
        elif name == "w":
            name = "resolution_parameters"
            index = (2, )
        elif name == "x":
            name = "resolution_parameters"
            index = (3, )
        elif name == "y":
            name = "resolution_parameters"
            index = (4, )
        elif name == "fract_x":
            name = "atom_fract_xyz"
            index = (0, index)
        elif name == "fract_y":
            name = "atom_fract_xyz"
            index = (1, index)
        elif name == "fract_z":
            name = "atom_fract_xyz"
            index = (2, index)
        if isinstance(index, int):
            index = (index, )
        dict_names.append((data_block.lower(), name.lower(), index))
    return dict_names