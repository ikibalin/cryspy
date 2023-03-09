"""
Describing of parent class ItemN.

It should be used ass parent to describe new item classes.
"""
import os
import os.path
from typing import NoReturn, Union
from pycifstar import Data, to_data

from cryspy.A_functions_base.function_1_markdown import md_to_html
from cryspy.A_functions_base.function_1_strings import find_prefix, \
    string_to_value_error_mark, value_error_mark_to_string
from cryspy.A_functions_base.function_1_objects import get_functions_of_objet


def string_to_fixed_or_refined(s_value: str):
    flag_refined = False
    try:
        value = float(s_value.split("(")[0])
        flag = (("(" in s_value) and (")" in s_value))
        if flag:
            flag_refined = True
    except ValueError:
        value = s_value
    except AttributeError:
        value = s_value
    return value, flag_refined


class ItemN(object):
    """Items data.

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
        res = str(self)
        if res == "":
            res = self.to_cif(flag_all_attributes=True)
        ls_out.append(res)

        method = self.methods_html()
        if method != "":
            ls_out.append(f"\n# Methods:\n{method:}\n")                
        return "\n".join(ls_out)

    def _repr_html_(self):
        """Representation in HTML format."""
        ls_html = [f"<h3>Object '{self.get_name():}'</h3>"]
        ls_html.append(self.attributes_to_html())

        report = self.report_html()
        if report != "":
            ls_html.append(f"<h3>Description </h3> {report:}")

        ls_html.append(f"<h3>Attributes and methods</h3>")
        try:
            names = sorted(self.ATTR_MANDATORY_NAMES)
            if len(names) != 0:
                ls_html.append("<b>Mandatory attributes: </b>")
                ls_html.append(f".{', .'.join(names):}<br>")
        except AttributeError:
            pass
        try:
            names = sorted(self.ATTR_OPTIONAL_NAMES)
            if len(names) != 0:
                ls_html.append("<b>Optional attributes: </b>")
                ls_html.append(f".{', .'.join(names):}<br>")
        except AttributeError:
            pass
        try:
            names = sorted(self.ATTR_INT_NAMES)
            if len(names) != 0:
                ls_html.append("<b>Internal attributes: </b>")
                ls_html.append(f".{', .'.join(names):}<br>")
        except AttributeError:
            pass
        try:
            names = sorted(self.ATTR_INT_PROTECTED_NAMES)
            if len(names) != 0:
                ls_html.append("<b>Internal protected attributes: </b>")
                ls_html.append(f".{', .'.join(names):}<br>")
        except AttributeError:
            pass

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
            Name of the attribute.

        Raises
        ------
        AttributeError
            Standard attribute error when it's necessary.

        Returns
        -------
        None.

        """
        if name.endswith("_as_string"):
            name_sh = name[:-(len("_as_string"))]
            try:
                val = getattr(self, name_sh)
            except AttributeError:
                return "."
            if val is None:
                return "."
            flag_ref = False
            keys = self.__dict__.keys()
            mark = ""
            if ((f"{name_sh:}_refinement" in keys) &
                    (f"{name_sh:}_sigma" in keys)):
                flag_ref = self.__dict__[f"{name_sh:}_refinement"]
                if (f"{name_sh:}_mark" in keys):
                    mark = self.__dict__[f"{name_sh:}_mark"]
            
            if flag_ref:
                val_sig = self.__dict__[f"{name_sh:}_sigma"]
                res = value_error_mark_to_string(val, val_sig, mark)
                
            else:
                s_format = None
                if (("D_FORMATS" in self.__dict__.keys()) |
                        ("D_FORMATS" in type(self).__dict__.keys())):
                    if name_sh in self.D_FORMATS.keys():
                        s_format = self.D_FORMATS[name_sh]
                if s_format is None:
                    res = f"{val:}"
                else:
                    res = s_format.format(val)
            return res
        if name in (self.ATTR_NAMES + self.ATTR_SIGMA
                    + self.ATTR_INT_PROTECTED_NAMES):
            raise AttributeError(f"Attribute '{name:}' is not defined in \
'{type(self).__name__:}'")
            # return None
        elif name in self.ATTR_INT_NAMES:
            if self.is_defined():
                # print("is_defined: ", self.is_defined())
                self.form_object()
                if name in  self.__dict__.keys():
                    return self.__dict__[name]
                else:
                    raise AttributeError(f"Attribute '{name:}' is not defined in \
'{type(self).__name__:}'")
            else:
                raise AttributeError(f"Attribute '{name:}' is not defined in \
'{type(self).__name__:}'")
                # return None
        else:
            raise AttributeError(f"Attribute '{name:}' is not defined in \
'{type(self).__name__:}'")

    def __setattr__(self, name: str, value) -> NoReturn:
        """
        Rules to set attribute.

        Parameters
        ----------
        name : str
            Name of attribute.
        value : TYPE
            Value of attribute.

        Returns
        -------
        None.

        """
        flag, flag_write = False, False
        flag_none_string = False
        if isinstance(value, str):
            flag_none_string = ((value == ".") | (value == "?"))

        if (value is None):
            flag, flag_write = True, True
            val_new = None
        elif flag_none_string:
            flag, flag_write = True, True
            val_new = None
        elif name in self.ATTR_NAMES:
            ind = self.ATTR_NAMES.index(name)
            val_type = self.ATTR_TYPES[ind]
            if ((val_type is bool) and (isinstance(value, str))):
                val_new = value == "True"
            else:
                val_new = val_type(value)
            flag, flag_write = True, True

            if name in self.D_CONSTRAINTS.keys():
                if val_new not in self.D_CONSTRAINTS[name]:
                    flag_write = False
            if name in self.D_MIN.keys():
                if val_new < self.D_MIN[name]:
                    val_new = self.D_MIN[name]
            if name in self.D_MAX.keys():
                if val_new > self.D_MAX[name]:
                    val_new = self.D_MAX[name]
            self.delete_internal_parameters()
        elif name in self.ATTR_SIGMA:
            val_new = float(value)
            flag, flag_write = True, True
        elif name in (self.ATTR_CONSTR_FLAG+self.ATTR_REF_FLAG):
            val_new = bool(value)
            flag, flag_write = True, True
        elif name in self.ATTR_CONSTR_MARK:
            val_new = str(value).strip()
            flag, flag_write = True, True

        if flag and flag_write:
            self.__dict__[name] = val_new
        elif not(flag):
            self.__dict__[name] = value

    def delete_internal_parameters(self):
        """
        Delete all internal parameters except protected one.

        Returns
        -------
        None.

        """
        l_del_name = []
        for key in self.__dict__.keys():
            flag = True
            if ((key == "D_MIN") | (key == "D_MAX")):
                flag = False
            elif key in (self.ATTR_NAMES + self.ATTR_SIGMA +
                         self.ATTR_CONSTR_FLAG + self.ATTR_REF_FLAG + self.ATTR_CONSTR_MARK +
                         self.ATTR_INT_PROTECTED_NAMES):
                flag = False
            if flag:
                l_del_name.append(key)
        for key in l_del_name:
            del self.__dict__[key]

    def is_attribute(self, name: str):
        """Give True if attribute is defined.
        
        Attention if parameter defined as None it's considered as defined.
        """
        flag = True
        try:
            getattr(self, name)
            # val = getattr(self, name)
            # if val is None:  # temporary solution
            #     flag = False
        except AttributeError or KeyError:
            flag = False
        return flag

    def is_defined(self) -> bool:
        """
        If all mandatory attributes is defined.

        Returns
        -------
        bool
            DESCRIPTION.

        """
        flag = True
        for name in self.ATTR_MANDATORY_NAMES:
            flag = self.is_attribute(name)
            if not(flag):
                break
        if (len(self.ATTR_MANDATORY_NAMES) == 0):
            flag = False
            for name in self.ATTR_OPTIONAL_NAMES:
                # print(name, getattr(self, name))
                if name in self.__dict__.keys():
                    flag = self.is_attribute(name)
                    if flag:
                        break
        return flag

    def form_object(self) -> NoReturn:
        """
        Calculate internal parameters for the object.

        Returns
        -------
        NoReturn
        """
        if ((type(self) == ItemN) and (len(self.ATTR_REF)==0)):
            l_ATTR_REF, l_ATTR_SIGMA = [], []
            l_ATTR_CONSTR_FLAG, l_ATTR_REF_FLAG = [], []
            l_ATTR_CONSTR_MARK = []
            variable_names = []
            for name in self.ATTR_NAMES:
                s_value = getattr(self, name)
                value, flag_refined = string_to_fixed_or_refined(s_value)
                if isinstance(value, float):
                    value, error, mark = string_to_value_error_mark(s_value)
                    name_sigma = f"{name:}_sigma"
                    name_constraint = f"{name:}_constraint"
                    name_refinement = f"{name:}_refinement"
                    name_mark = f"{name:}_mark"

                    self.__dict__[name] = value
                    self.__dict__[name_sigma] = error
                    self.__dict__[name_constraint] = False
                    self.__dict__[name_refinement] = flag_refined
                    self.__dict__[name_mark] = mark

                    l_ATTR_REF.append(name)
                    l_ATTR_SIGMA.append(name_sigma)
                    l_ATTR_CONSTR_FLAG.append(name_constraint)
                    l_ATTR_REF_FLAG.append(name_refinement)
                    l_ATTR_CONSTR_MARK.append(name_mark)
            self.__dict__["ATTR_REF"] = tuple(l_ATTR_REF)
            self.__dict__["ATTR_SIGMA"] = tuple(l_ATTR_SIGMA)
            self.__dict__["ATTR_CONSTR_FLAG"] = tuple(l_ATTR_CONSTR_FLAG)
            self.__dict__["ATTR_REF_FLAG"] = tuple(l_ATTR_REF_FLAG)
            self.__dict__["ATTR_CONSTR_MARK"] = tuple(l_ATTR_CONSTR_MARK)

    def to_cif(self, separator: str = "_", flag_all_attributes:bool=False) -> str:
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
        prefix = self.PREFIX
        for name, name_cif in zip(self.ATTR_NAMES, self.ATTR_CIF):
            if name_cif == "":
                name_cif = name
            flag_value = self.is_attribute(name)
            if (flag_value | flag_all_attributes):
                if flag_value:
                    s_val = str(getattr(self, name))
                else:
                    s_val = "."
                l_s_val = s_val.split("\n")
                if len(l_s_val) > 1:
                    ls_out.append(f"_{prefix:}{separator:}{name_cif:}")
                    ls_out.append(";")
                    ls_out.extend(l_s_val)
                    ls_out.append(";")
                else:
                    if len(s_val.split(" ")) > 1:
                        ls_out.append(
                            f"_{prefix:}{separator:}{name_cif:} \"{s_val:}\"")
                    else:
                        s_val = getattr(self, f"{name:}_as_string")
                        ls_out.append(
                            f"_{prefix:}{separator:}{name_cif:} {s_val:}")
        return "\n".join(ls_out)

    def attributes_to_html(self) -> str:
        """Representation of defined parameters in HTML format.
        """
        ls_html = ["<table>"]
        ls_html.append("<tr><th>Attributes</th><th> Value </th></tr>")
        for name in self.ATTR_NAMES:
            flag_value = self.is_attribute(name)
            if flag_value:
                if name in self.ATTR_REF:
                    s_val = str(getattr(self, f"{name:}_as_string"))
                else:
                    s_val = str(getattr(self, name))
                if len(s_val.split("\n")) > 1:
                    ls_html.append(f"<tr><td>.{name:}</td><td> \
<i>multiline expression</i> </td></tr>")
                elif len(s_val)>40:                    
                    ls_html.append(f"<tr><td>.{name:}</td><td>{s_val[:40]:}...\
</td></tr>")
                else:
                    ls_html.append(f"<tr><td>.{name:}</td><td>{s_val:}\
</td></tr>")
        ls_html.append("</table>")
        return " ".join(ls_html)

    def get_name(self) -> str:
        """Get name."""
        return self.PREFIX.lower()

    def get_variable_names(self) -> list:
        """
        Get names of variable as a list.

        (((#prefix, #NAME), (#attribute, #index))

        Returns
        -------
        list
            List of names of variable.

        """
        if type(self) == ItemN:
            self.form_object()
        atr_ref = self.ATTR_REF
        atr_ref_flag = self.ATTR_REF_FLAG
        prefix = self.PREFIX
        variable_names = [((prefix, None), (name, None)) 
            for name, name_flag in
                zip(atr_ref, atr_ref_flag) if getattr(self, name_flag)]
        return variable_names

    def is_variables(self) -> bool:
        """Define is there variables or not."""
        atr_ref = self.ATTR_REF
        atr_ref_flag = self.ATTR_REF_FLAG
        flag = False
        for name, name_flag in zip(atr_ref, atr_ref_flag):
            if getattr(self, name_flag):
                flag = True
                break
        return flag

    def get_variable_by_name(self, name: tuple) -> Union[float, int, str]:
        """
        Get variable given by name.

        Parameters
        ----------
        name : tuple
            (((#prefix, ), (#attribute, ))

        Returns
        -------
        Union[float, int, str]
            Value.

        """
        
        prefix_t = name[0]
        if prefix_t[0] != self.PREFIX:
            return None
        if len(name) == 1:
            return self
        attr_t = name[1]
        attr_name = attr_t[0]
        return getattr(self, attr_name)

    def set_variable_by_name(self, name: tuple, value) -> NoReturn:
        """
        Set value to variable given by name.

        Parameters
        ----------
        name : tuple
            (((#prefix, ), (#attribute, ))
        value : TYPE
            DESCRIPTION.

        Returns
        -------
        NoReturn

        """
        prefix_t, attr_t = name
        if prefix_t[0] != self.PREFIX:
            return
        attr_name = attr_t[0]
        setattr(self, attr_name, value)

    def _base_attributes_itemn(self, prefix, name_opt_cif, name_opt=None):
        name_mand, name_mand_cif = (), ()
        attr_ref = ()
        if name_opt is None:
            name_opt = name_opt_cif
        self.__dict__["PREFIX"] = prefix
        self.__dict__["ATTR_MANDATORY_CIF"] = name_mand_cif
        self.__dict__["ATTR_OPTIONAL_CIF"] = name_opt_cif
        self.__dict__["ATTR_REF"] = attr_ref
        self.__dict__["ATTR_CIF"] = name_mand_cif+name_opt_cif
        self.__dict__["ATTR_MANDATORY_NAMES"] = name_mand
        self.__dict__["ATTR_OPTIONAL_NAMES"] = name_opt
        self.__dict__["ATTR_NAMES"] = name_mand+name_opt
        types_mand = tuple(len(name_mand)*[str])
        types_opt = tuple(len(name_opt)*[str])
        self.__dict__["ATTR_MANDATORY_TYPES"] = types_mand
        self.__dict__["ATTR_OPTIONAL_TYPES"] = types_opt
        self.__dict__["ATTR_TYPES"] = types_mand + types_opt
        self.__dict__["ATTR_REF"] = ()
        self.__dict__["ATTR_SIGMA"] = ()
        self.__dict__["ATTR_CONSTR_FLAG"] = ()
        self.__dict__["ATTR_REF_FLAG"] = ()
        self.__dict__["ATTR_CONSTR_MARK"] = ()
        self.__dict__["D_FORMATS"] = {}
        self.__dict__["D_CONSTRAINTS"] = {}
        self.__dict__["D_DEFAULT"] = {}
        self.__dict__["ATTR_INT_NAMES"] = ()
        self.__dict__["D_MIN"] = {}
        self.__dict__["D_MAX"] = {}
        self.__dict__["ATTR_INT_PROTECTED_NAMES"] = (
            "PREFIX", "ATTR_MANDATORY_CIF", "ATTR_OPTIONAL_CIF",
            "ATTR_MANDATORY_NAMES", "ATTR_OPTIONAL_NAMES",
            "ATTR_REF", "ATTR_NAMES", "ATTR_CIF",
            "ATTR_MANDATORY_TYPES", "ATTR_TYPES", "ATTR_REF",
            "ATTR_SIGMA", "ATTR_CONSTR_FLAG", "ATTR_REF_FLAG", "ATTR_CONSTR_MARK",
            "D_FORMATS", "D_CONSTRAINTS", "D_DEFAULT",
            "ATTR_INT_NAMES", "ATTR_INT_NAMES",
            "ATTR_INT_PROTECTED_NAMES", "D_MIN", "D_MAX")

    @classmethod
    def from_cif(cls, string: str):
        """
        Create object of Cell class from string.

        Parameters
        ----------
        cls : TYPE
            DESCRIPTION.
        string : str
            string with descritption of the file.

        Returns
        -------
        ItemN
            DESCRIPTION.

        """
        cif_data = Data()
        flag = cif_data.take_from_string(string)
        item = None
        if cls is ItemN:
            l_name = [item.name for item in cif_data.items]
            if len(l_name) == 1:
                name = l_name[0]
                if name.find(".") != -1:
                    prefix = name.split(".")[0]
                else:
                    prefix = "_"+name.split("_")[1]
            else:
                prefix = find_prefix(*tuple(l_name))
            l_attr = [name[len(prefix)+1:] for name in l_name]

            prefix = prefix.strip("_").strip(".")
            name_opt_cif = tuple(l_attr)
            item = cls()
            item._base_attributes_itemn(prefix, name_opt_cif)
            prefix = item.PREFIX
            name_mand_cif = item.ATTR_MANDATORY_CIF
            name_opt_cif = item.ATTR_OPTIONAL_CIF
            name_mand = item.ATTR_MANDATORY_NAMES
            name_opt = item.ATTR_OPTIONAL_NAMES
            attr_ref = item.ATTR_REF
        else:
            prefix = cls.PREFIX
            name_mand_cif = cls.ATTR_MANDATORY_CIF
            name_opt_cif = cls.ATTR_OPTIONAL_CIF
            name_mand = cls.ATTR_MANDATORY_NAMES
            name_opt = cls.ATTR_OPTIONAL_NAMES
            attr_ref = cls.ATTR_REF

        name = name_mand + name_opt
        name_cif = name_mand_cif + name_opt_cif

        separator = "_"
        flag = all([cif_data.is_value("_" + prefix + separator + _attr) for
                    _attr in name_mand_cif])
        if not (flag):
            separator = "."
            flag = all(
                [cif_data.is_value("_" + prefix + separator + _attr) for _attr
                 in name_mand_cif])
        if flag:
            separator = "_"
            flag_2 = any([cif_data.is_value("_" + prefix + separator + _attr)
                          for _attr in name_cif])
            if not (flag_2):
                separator = "."
                flag_2 = any([cif_data.is_value("_" + prefix + separator +
                                                _attr) for _attr in name_cif])
            if flag_2:
                if item is None:
                    item = cls()
                for _attr, _cif_attr in zip(name, name_cif):
                    _cif_attr_full = "_" + prefix + separator + _cif_attr
                    if cif_data.is_value(_cif_attr_full):
                        if _attr in attr_ref:
                            value, error, mark = \
                                string_to_value_error_mark(
                                    cif_data[_cif_attr_full].value)
                            setattr(item, _attr, value)
                            if error is not None:
                                setattr(item, f"{_attr:}_sigma", error)
                                setattr(item, f"{_attr:}_refinement", True)
                                setattr(item, f"{_attr:}_mark", mark)
                            else:
                                setattr(item, f"{_attr:}_refinement", False)
                        else:
                            setattr(item, _attr,
                                    cif_data[_cif_attr_full].value)
            else:
                item = None
        else:
            item = None
        if item is not None:
            if item.is_defined():
                item.form_object()
        return item

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
        for attr in (self.ATTR_NAMES + self.ATTR_SIGMA +
                     self.ATTR_CONSTR_FLAG + self.ATTR_REF_FLAG + self.ATTR_CONSTR_MARK):
            if obj.is_attribute(attr):
                setattr(self, attr, getattr(obj, attr))

    def report(self):
        return ""

    def report_html(self):
        return md_to_html(self.report())

    def plots(self):
        return []
    
    def fix_variables(self):
        """Fix variables."""
        for name in self.get_variable_names():
            name_ref = name[:-1] + ((f"{name[-1][0]:}_refinement",
                                     name[-1][1]), )
            self.set_variable_by_name(name_ref, False)

    def refine_all_variables(self):
        """Fix variables."""
        for name in self.ATTR_REF:
            name_constr = f"{name:}_constraint"
            constr = getattr(self, name_constr)
            if not(constr):
                name_ref = f"{name:}_refinement"
                setattr(self, name_ref, True)

    def set_variable(self, name: str, index=None):
        """Set refinement for variable given by name.
        
        Index parameters is used only for objects given as a matrix.
        """
        name_sh = name.strip(".").lower()
        if name_sh in self.ATTR_REF:
            flag_1 = self.is_attribute(name_sh)
            flag_2 = not(getattr(self, f"{name_sh:}_constraint"))
            if flag_1 & flag_2:
                setattr(self, f"{name_sh:}_refinement", True)
        