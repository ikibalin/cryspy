"""
Describing of parent class LoopN.

It should be used ass parent to describe new loop classes.
"""
import os
import os.path
import copy
from typing import NoReturn, Union
import numpy
from pycifstar import Data, to_data

from cryspy.A_functions_base.function_1_markdown import md_to_html
from cryspy.A_functions_base.function_1_strings import find_prefix, \
    string_to_value_error_mark
from cryspy.A_functions_base.function_1_objects import \
    get_functions_of_objet, get_table_html_for_variables

from cryspy.B_parent_classes.cl_1_item import ItemN


class LoopN(object):
    """Loop data.

    It is internal class of cryspy library.
    You should use it only to create your own classes.
    """

    def __repr__(self):
        """
        Magic method print() is redefined.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        ls_out = [f"# Object '{self.get_name():}':\n"]
        ls_out.append(str(self))

        method = self.methods_html()
        if method != "":
            ls_out.append(f"\n# Methods:\n{method:}\n")                
        return "\n".join(ls_out)

    def _repr_html_(self):
        """Representation in HTML format."""
        ls_html = [f"<h3>Object '{self.get_name():}'</h3>"]
        ls_html.append(self.attributes_to_html())

        ls_html.append(get_table_html_for_variables(self))
        
        report = self.report_html()
        if report != "":
            ls_html.append(f"<h3>Description </h3> {report:}")

        ls_html.append(f"<h3>Methods</h3>")
        method = self.methods_html()
        if method != "":
            ls_html.append(f"<b>Methods: </b> {method:}")

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
        item_class = self.ITEM_CLASS
        if item_class is ItemN:
            item_class = self.items[0]
        if name.startswith("numpy_"):
            name_sh = name[6:]
            if name_sh in (item_class.ATTR_NAMES + item_class.ATTR_SIGMA +
                           item_class.ATTR_CONSTR_FLAG +
                           item_class.ATTR_REF_FLAG + item_class.ATTR_CONSTR_MARK):
                l_val = [getattr(item, name_sh) for item in self.items]
                if name_sh in item_class.ATTR_NAMES:
                    type_array = item_class.ATTR_TYPES[
                        item_class.ATTR_NAMES.index(name_sh)]
                elif name_sh in (item_class.ATTR_CONSTR_FLAG +
                                 item_class.ATTR_REF_FLAG):
                    type_array = bool
                elif name_sh in item_class.ATTR_CONSTR_MARK:
                    type_array = str
                else:
                    type_array = float
                res = numpy.array(l_val, dtype=type_array)
                self.__dict__[name] = res
                return res
            elif name_sh in (item_class.ATTR_INT_NAMES +
                             item_class.ATTR_INT_PROTECTED_NAMES):
                l_val = [getattr(item, name_sh) for item in self.items]
                type_array = type(l_val[0])
                res = numpy.array(l_val, dtype=type_array)
                self.__dict__[name] = res
                return res
        elif name in (item_class.ATTR_NAMES + item_class.ATTR_SIGMA +
                      item_class.ATTR_CONSTR_FLAG + item_class.ATTR_REF_FLAG +
                      item_class.ATTR_CONSTR_MARK):
            res = [getattr(item, name) for item in self.items]
            return res
        elif name in (item_class.ATTR_INT_NAMES +
                      item_class.ATTR_INT_PROTECTED_NAMES):
            res = [getattr(item, name) for item in self.items]
            return res
        raise AttributeError(
            f"'{type(self).__name__:}' object has no attribute '{name:}'")

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
        flag_direct = True
        if name == "items":
            flag = all([isinstance(val, self.ITEM_CLASS) for val in value])
            if not(flag):
                flag_direct = False
                val_new = [val for val in value if isinstance(val,
                                                              self.ITEM_CLASS)]
        elif name == "loop_name":
            flag_direct = False
            val_new = str(value).strip()

        if flag_direct:
            self.__dict__[name] = value
        else:
            self.__dict__[name] = val_new

    def is_attribute(self, name: str):
        """Give True if all attributes are defined."""
        flag = all([item.is_attribute(name) for item in self.items])
        return flag

    def numpy_to_items(self) -> NoReturn:
        """
        Transform all internal numpy arrays to elements of items.

        Returns
        -------
        NoReturn
            DESCRIPTION.

        """
        loop_dict = self.__dict__
        item_class = self.ITEM_CLASS
        keys = loop_dict.keys()
        del_name = []
        for name in keys:
            if name.startswith("numpy_"):
                name_sh = name[6:]
                if name_sh in (item_class.ATTR_NAMES + item_class.ATTR_SIGMA +
                               item_class.ATTR_CONSTR_FLAG +
                               item_class.ATTR_REF_FLAG + 
                               item_class.ATTR_CONSTR_MARK):
                    numpy_val = self.__dict__[name]
                    if len(self.items) == 0:
                        items = [item_class() for val in numpy_val]
                        self.items = items
                    for item, val in zip(self.items, numpy_val):
                        setattr(item, name_sh, val)
                    del_name.append(name)
        for name in del_name:
            del self.__dict__[name]

    def form_object(self):
        pass

    @classmethod
    def define_by_items(cls, l_item: list, loop_name: str = ""):
        """Define LoopN by defined items.
        
        (it's for the future)
        """
        if len(l_item) == 0:
            raise ValueError("Number of items is zero")
        if not(isinstance(l_item[0], ItemN)):
            raise ValueError("type of items is not ItemN")
        item_class = type(l_item[0])
        flag = all([type(item) is item_class for item in l_item])
        if not(flag):
            raise ValueError("There are different types of given items")
        obj = cls()
        if cls is LoopN:
            obj.__dict__["ITEM_CLASS"] = item_class
            obj.__dict__["ATTR_INDEX"] = None
            obj.__dict__["items"] = l_item
            obj.__dict__["loop_name"] = loop_name
        elif obj.ITEM_CLASS is item_class:
            obj.items = l_item
        else:
            raise ValueError("Type of items does not correspond to loop class")            
        return obj

            
    @classmethod
    def from_cif(cls, string: str):
        """
        Create object of Cell class from string.

        Parameters
        ----------
        cls : TYPE
            DESCRIPTION.
        string : str
            DESCRIPTION.

        Returns
        -------
        obj : TYPE
            DESCRIPTION.

        """
        cif_data = Data()
        cif_data.take_from_string(string)
        obj = None
        if cls is LoopN:
            loop = cif_data.loops[0]
            l_name = loop.names
            prefix = find_prefix(*tuple(l_name))
            l_attr = [name[len(prefix):].strip("_").strip(".") for name in l_name]
            l_cif_attr = [_.lower() for _ in l_attr]
            if "" in l_attr:
                l_attr[l_attr.index("")] = "label"
            prefix = prefix.strip("_").strip("_")
            obj = cls()
            obj.__dict__["ITEM_CLASS"] = ItemN
            obj.__dict__["ATTR_INDEX"] = None
            item_class = ItemN
            obj.__dict__["items"] = []
            obj.__dict__["loop_name"] = loop.name
        else:
            item_class = cls.ITEM_CLASS
            l_attr = item_class.ATTR_NAMES
            l_cif_attr = item_class.ATTR_CIF

            l_cif_attr = [_.lower() for _ in l_cif_attr]
            prefix = item_class.PREFIX

        for cif_loop in cif_data.loops:
            # cond = ("_" + prefix.lower()) == cif_loop.prefix.lower()
            # print("cif_loop.prefix: ", cif_loop.prefix)
            # print(("_" + prefix.lower()) == cif_loop.prefix.lower(),
            #       cif_loop.prefix.lower().startswith("_" + prefix.lower()))
            # if ("_" + prefix.lower()) == cif_loop.prefix.lower():
            if (cif_loop.prefix.lower().startswith("_" + prefix.lower())):
                prefix = f"_{prefix:}"
                # prefix = cif_loop.prefix
                l_name = cif_loop.names
                loop_name = cif_loop.name
                l_name_short = [_[(len(prefix) + 1):] for _ in l_name]

                if obj is None:
                    obj = cls(loop_name=loop_name)
                _i = 0
                for _name, _name_short in zip(l_name, l_name_short):
                    if _name_short.lower() in l_cif_attr:
                        _name_short_obj = l_attr[l_cif_attr.index(_name_short)]
                    else:
                        _name_short_obj = _name_short.lower()

                    if _name_short_obj.strip("_").strip(".") in l_attr:
                        if _i == 0:
                            items = []
                            for _val in cif_loop[_name]:
                                if item_class is ItemN:
                                    item = item_class()
                                    item._base_attributes_itemn(
                                        prefix.strip("_").strip("."),
                                        tuple(l_cif_attr), tuple(l_attr))
                                else:
                                    item = item_class()

                                if _name_short_obj in item.ATTR_REF:
                                    value, error, mark = string_to_value_error_mark(_val)
                                    
                                    setattr(item, _name_short_obj, value)
                                    if error is not None:
                                        setattr(item, f"{_name_short_obj:}_sigma", error)
                                        setattr(item, f"{_name_short_obj:}_refinement", True)
                                        setattr(item, f"{_name_short_obj:}_mark", mark)
                                    else:
                                        setattr(
                                            item,
                                            f"{_name_short_obj:}_refinement",
                                            False)
                                else:
                                    setattr(item, _name_short_obj, _val)
                                items.append(item)
                        else:
                            for _val, item in zip(cif_loop[_name], items):
                                if _name_short_obj in item.ATTR_REF:
                                    value, error, mark = string_to_value_error_mark(_val)
                                    setattr(item, _name_short_obj, value)
                                    if error is not None:
                                        setattr(item, f"{_name_short_obj:}_sigma", error)
                                        setattr(item,f"{_name_short_obj:}_refinement", True)
                                        setattr(item,f"{_name_short_obj:}_mark", mark)
                                    else:
                                        setattr(item, f"{_name_short_obj:}_refinement", False)
                                else:
                                    setattr(item, _name_short_obj, _val)
                        _i += 1
                    else:
                        pass

                obj.items = items
                break

        return obj

    def to_cif(self, separator: str = "_") -> str:
        """
        Give data in CIF/STAR format.

        Parameters
        ----------
        separator : str, optional
            DESCRIPTION. The default is "_".

        Returns
        -------
        str
            DESCRIPTION.

        """
        ls_out = []
        item_class = self.ITEM_CLASS
        if self.loop_name is None:
            ls_out.append("loop_")
        else:
            ls_out.append(f"loop_{self.loop_name:}")

        if len(self.items) == 0:
            if item_class == ItemN:
                return ""
            prefix = item_class.PREFIX
            for name_cif in item_class.ATTR_CIF:
                ls_out.append(f"_{prefix:}{separator:}{name_cif:}")
            ls_out.append(" ".join(len(item_class.ATTR_CIF)*["."]))
            return "\n".join(ls_out)

        item_0 = self.items[0]
        prefix = item_0.PREFIX
        ls_out_2 = []
        for name, name_cif in zip(item_0.ATTR_NAMES, item_0.ATTR_CIF):
            flag_value = item_0.is_attribute(name)
            if flag_value:
                if name_cif == "":
                    ls_out.append(f"_{prefix:}")
                else:
                    ls_out.append(f"_{prefix:}{separator:}{name_cif:}")
                list_value = [getattr(item, f"{name}_as_string")
                              for item in self.items]
                if len(ls_out_2) == 0:
                    n_max = max(map(len, [str(h) for h in list_value])) + 2
                    for s_val in list_value:
                        if ((len(s_val.split(" ")) > 1) |
                                s_val.startswith("_")):
                            ls_out_2.append([("\""+s_val+"\"").rjust(n_max)])
                        else:
                            ls_out_2.append([s_val.rjust(n_max)])
                else:
                    n_max = max(map(len, [str(h) for h in list_value])) + 2
                    for list_out, s_val in zip(ls_out_2, list_value):
                        if ((len(s_val.split(" ")) > 1) | 
                                s_val.startswith("_")):
                            list_out.append(("\""+s_val+"\"").rjust(n_max))
                        else:
                            list_out.append(s_val.rjust(n_max))

        for list_out in ls_out_2:
            s_line = " ".join(list_out)
            ls_out.append(s_line)
        return "\n".join(ls_out)

    def attributes_to_html(self) -> str:
        """Representation of defined parameters in HTML format.
        """
        ls_html = ["<table>"]
        item_class = self.ITEM_CLASS
        if item_class is ItemN:
            if len(self.items) != 0:
                item_n = self.items[0]
                attribute_names = item_n.ATTR_NAMES
            else:
                return ""
        else:
            attribute_names = item_class.ATTR_NAMES
        flags = [False for name in attribute_names]
        for item in self.items:
            for ind_name, name in enumerate(attribute_names):
                flag_value = item.is_attribute(name)
                if flag_value:
                    flags[ind_name] = True

        ls_html.append(
            "<tr><th> </th>"+"".join([f"<th>.{name:}</th>" for name, flag in zip(
                attribute_names, flags)  if flag])+"</tr>")

        for ind_item, item in enumerate(self.items):
            ls_line = []
            for name, flag in zip(attribute_names, flags):
                if flag:
                    flag_value = item.is_attribute(name)
                    if flag_value:
                        value = getattr(item, f"{name}_as_string")
                        ls_line.append(f"<td>{value:}</td>")
                    else:
                        ls_line.append("<td>.</td>")
            ls_html.append(f"<tr><th>[{ind_item}]</th>"+"".join(ls_line)+"</tr>")
            if ind_item > 10:
                ls_line = ["<td>...</td>" for flag in flags if flag ]
                ls_html.append(f"<tr><th>[..., {len(self.items):}]</th>"+"".join(ls_line)+
                               "</tr>")
                break

        ls_html.append("</table>")
        return " ".join(ls_html)

    def save_to_csv(self, file_name: str = "output.csv", delimiter: str = ";")\
            -> NoReturn:
        """
        Save loop into csv format.

        Parameters
        ----------
        file_name : str, optional
            Output file name. The default is "output.csv".
        delimiter : str, optional
            Delimeter. The default is ";".

        Returns
        -------
        NoReturn

        """
        item_class = self.ITEM_CLASS
        l_attr_print = item_class.ATTR_NAMES
        l_cif_attr = item_class.ATTR_NAMES_CIF
        """
        l_flag = [all([_item.is_defined_attribute(_attr) for _item in self.item]) for _attr in l_attr_print]
        l_attr_print = [_ for _, _flag in zip(l_attr_print, l_flag) if _flag]
        l_cif_attr = [_ for _, _flag in zip(l_cif_attr, l_flag) if _flag]
        ls_out = [(f"{delimiter:}").join([f"{_attr:}" for _attr in l_cif_attr])]
        for _item in self.item:
            ls_out.append(_item.print_attribute(l_attr_print, delimiter=delimiter))
        s_cont = "\n".join(ls_out)
        with open(file_name, "w") as fid:
            fid.write(s_cont)
        """

    def __getitem__(self, name):
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
        attr_index = self.ATTR_INDEX
        if attr_index is None:
            if isinstance(name, int):
                return self.items[name]
            else:
                return None
        else:
            for item in self.items:
                if name == getattr(item, attr_index):
                    return item
            if isinstance(name, int):
                return self.items[name]

    def get_name(self) -> str:
        """Get name."""
        if self.ITEM_CLASS is ItemN:
            if len(self.items) > 0:
                name = self.items[0].PREFIX
            else:
                return ""
        else:
            name = self.ITEM_CLASS.PREFIX
        loop_name = self.loop_name
        if ((loop_name is None) | (loop_name == "")):
            name = f"{name:}"
        else:
            name = f"{name:}_{loop_name:}"
        return name.lower()

    def get_variable_names(self) -> list:
        """
        Get names of variable as a list.

        (((#prefix, #NAME), (#attribute, #index))

        Returns
        -------
        list
            List of names of variable.

        """
        item_class = self.ITEM_CLASS
        if item_class is ItemN:
            if len(self.items) != 0:
                prefix = self.items[0].PREFIX
            else:
                return []
        else:
            prefix = item_class.PREFIX
        loop_name = self.loop_name
        if isinstance(loop_name, str):
            loop_name = loop_name.lower()
        l_var = []
        for ind, item in enumerate(self.items):
            l_var.extend([((prefix, loop_name), (name[1][0], ind))
                          for name in item.get_variable_names()])
        return l_var

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
            (((#prefix, #loop_name), (#attribute, #index_item))

        Returns
        -------
        Union[float, int, str]
            DESCRIPTION.

        """
        item_class = self.ITEM_CLASS
        if item_class is ItemN:
            if len(self.items) != 0:
                prefix = self.items[0].PREFIX
            else:
                return
        else:
            prefix = item_class.PREFIX
        prefix_t = name[0]
        if isinstance(prefix_t[1], str):
            prefix_t_2 = prefix_t[1].lower()
        else:
            prefix_t_2 = prefix_t[1]

        loop_name = self.loop_name
        if isinstance(loop_name, str):
            loop_name = loop_name.lower()

        if (prefix_t[0], prefix_t_2) != (prefix, loop_name):
            return None

        if len(name) == 1:
            return self
        attr_t = name[1]
        attr_name = attr_t[0]
        if len(attr_t) == 1:
            return getattr(self.items, attr_name)
        item_index = attr_t[1]
        if attr_name is None:
            return self.items[item_index]
        return getattr(self.items[item_index], attr_name)

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
        item_class = self.ITEM_CLASS
        if item_class is ItemN:
            if len(self.items) != 0:
                prefix = self.items[0].PREFIX
            else:
                return
        else:
            prefix = item_class.PREFIX
        
        prefix_t, attr_t = name
        if prefix_t != (prefix, self.loop_name):
            return
        attr_name, item_index = attr_t
        setattr(self.items[item_index], attr_name, value)

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
                break
        return flag

    @classmethod
    def from_cif_file(cls, f_name: str):
        """Read from cif file."""
        if not(os.path.isfile(f_name)):
            raise UserWarning(f"File {f_name:} is not found.")
            return None
        str_from_cif = str(to_data(f_name))
        obj = cls.from_cif(str_from_cif)
        obj.file_input = f_name
        return obj

    def copy_from(self, obj):
        """Copy attributes from obj to self."""
        
        if type(obj) is not type(self):
            return
        self.items.clear()
        self.items = [type(item).from_cif(str(item))  for item in obj.items]
        
        

    def report(self):
        return ""

    def report_html(self):
        return md_to_html(self.report())

    def plots(self):
        return []

    def fix_variables(self):
        """Fix variables."""
        for item in self.items:
            item.fix_variables()

    def refine_all_variables(self):
        for item in self.items:
            item.refine_all_variables()

    def set_variable(self, name: str, index=None):
        """Set refinement for variable given by name.
        
        Index parameters is used only for objects given as a matrix.
        """
        name_sh = name.strip(".").lower()
        item_class = self.ITEM_CLASS
        if item_class is ItemN:
            if len(self.items) != 0:
                att_ref = self.items[0].ATTR_REF
            else:
                return
        else:
            att_ref = item_class.ATTR_REF

        if name_sh in att_ref:
            if index is None:
                for item in self.items:
                    item.set_variable(name_sh, index=index)
            else:
                item = self[index]
                item.set_variable(name_sh, index=None)
    
    def get_dictionary(self):
        res = {}
        return res
    
    def take_parameters_from_dictionary(self, ddict_diffrn, l_parameter_name: list=None, l_sigma: list=None):
        """
        """
        pass
        return None

def get_prefix_of_loop(item: LoopN):
    if not(isinstance(item, LoopN)):
        return ""
    if item.ITEM_CLASS is ItemN:
        if len(item.items) > 0:
            prefix = item.items[0].PREFIX
        else:
            return ""
    else:
        prefix = item.ITEM_CLASS.PREFIX
    return prefix