"""Parent class DataN."""
import os
import os.path
from warnings import warn
from typing import Union, NoReturn
from pycifstar import Global, to_global

from cryspy.A_functions_base.function_1_markdown import md_to_html
from cryspy.A_functions_base.function_1_objects import \
    get_functions_of_objet, get_table_html_for_variables

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN
from cryspy.B_parent_classes.cl_3_data import DataN


class GlobalN(object):
    """Data container of data bloks, loops and items."""

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

    def __getattr__(self, name):
        """
        Magic method __getattr__ is slightly changed for special attributes.

        Parameters
        ----------
        name : TYPE
            DESCRIPTION.

        Raises
        ------
        AttributeError
            DESCRIPTION.

        Returns
        -------
        res : TYPE
            DESCRIPTION.

        """
        for item in self.items:
            if name.lower() == item.get_name():
                return item
        return None

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
                    warn(f"Access to variable by '{name_new:}'.", UserWarning)

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
            warn("Double items were given.", UserWarning)
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
                item_cls = item.ITEM_CLASS
                if item_cls is ItemN:
                    prefix = item[0].PREFIX
                else:
                    prefix = item_cls.PREFIX
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
                item_cls = item.ITEM_CLASS
                if item_cls is ItemN:
                    prefix = item[0].PREFIX
                else:
                    prefix = item_cls.PREFIX
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
                         UserWarning)
                    break
                elif isinstance(item, LoopN):
                    warn(f"Item {item.ITEM_CLASS.PREFIX:} is not fully described.",
                         UserWarning)
                    break
        if flag:
            cls_items = [type(item) for item in self.items]
            for cls_mand in self.CLASSES_MANDATORY:
                if not(cls_mand in cls_items):
                    flag = False
                    warn(f"The object of {cls_mand.__name__:} is not defined.",
                         UserWarning)
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
                     UserWarning)
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
            raise UserWarning(f"File {f_name:} is not found.")
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
        dict_out = {}
        for item in self.items:
            if "get_dictionary" in dir(item):
                name = item.get_name()
                dict_out[name] = item.get_dictionary()

        variable_names = self.get_variable_names()
        dict_names = transfer_cif_names_to_dict_names(variable_names)

        mark_names = [hh[:-1] + ((f"{hh[-1][0]:}_mark", hh[-1][1]), ) for hh in variable_names]
        marks = [self.get_variable_by_name(hh) for hh in mark_names]
        constraint_labels = [hh.rstrip("-+") for hh in marks]
        linear_constraints = []
        for constraint_label in set(constraint_labels).difference([""]):
            ind_1 = constraint_labels.index(constraint_label)
            for i_hh, hh in enumerate(constraint_labels[ind_1+1:]):
                if hh == constraint_label:
                    ind_2 = i_hh+ind_1+1
                    name_1 = dict_names[ind_1]
                    if marks[ind_1].endswith("-"):
                        coeff_1 = -1.
                    else:
                        coeff_1 = 1.
                    
                    name_2 = dict_names[ind_2]
                    if marks[ind_2].endswith("-"):
                        coeff_2 = 1.
                    else:
                        coeff_2 = -1.
                    linear_constraints.append([(coeff_1, name_1), (coeff_2, name_2)])
        if len(linear_constraints) > 0:
            dict_out["linear_constraints"] = linear_constraints
        try:
            expression = self.punishment.function
            d_marked = {}
            for rcif_name, dict_name in zip(variable_names, dict_names):
                rcif_name_mark = rcif_name[:-1]+ ((rcif_name[-1][0]+"_mark", rcif_name[-1][1]),)
                mark = self.get_variable_by_name(rcif_name_mark)
                if mark != "":
                    d_marked[mark] = dict_name
            dict_out["punishment_function"] = (expression, d_marked)
        except:
            pass
        return dict_out



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
        if isinstance(index, int):
            index = (index, )
        dict_names.append((data_block.lower(), name.lower(), index))
    return dict_names