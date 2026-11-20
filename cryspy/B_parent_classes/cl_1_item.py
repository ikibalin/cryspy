"""
Describing of parent class ItemN.

It should be used ass parent to describe new item classes.
"""
import copyreg
import os
import re
import os.path
from typing import NoReturn, Union
from pycifstar import Data, to_data

from cryspy.A_functions_base.function_1_markdown import md_to_html
from cryspy.A_functions_base.function_1_strings import find_prefix, \
    string_to_value_error_mark, value_error_mark_to_string
from cryspy.A_functions_base.function_1_objects import get_functions_of_objet
from cryspy.exceptions import RcifFileNotFoundError, RcifFormatError
from cryspy.B_parent_classes.cif_attributes import install_cif_descriptors


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


INTEGER_LITERAL = re.compile(r"[+-]?(0|[1-9][0-9]*)")
CODE_LITERAL = re.compile(r"[+-]?0[0-9]+")  # '007': a code, not a number


def take_value_type(l_value) -> type:
    """Type of an attribute of a category unknown to CrysPy from its values
    in the file: int when every value is an integer (a year, a number of
    a list; '007' stays a string), float when they are numbers with or
    without sigma (a refinable parameter), str otherwise. '.' and '?' are
    not counted."""
    value_type = None
    for value in l_value:
        if not isinstance(value, str) or value in (".", "?"):
            continue
        if CODE_LITERAL.fullmatch(value) is not None:
            return str
        if INTEGER_LITERAL.fullmatch(value) is not None:
            type_value = int
        elif isinstance(string_to_fixed_or_refined(value)[0], float):
            type_value = float
        else:
            return str
        if value_type is not float:
            value_type = type_value
    return str if value_type is None else value_type


def take_float_format(l_value) -> Union[str, None]:
    """Format of the numbers of an attribute of a category unknown to
    CrysPy that writes them back as in the file: as many decimals as the
    longest of them ('0.00006' stays '0.00006', not '6e-05'); None when a
    value is in exponential form. Values with sigma are written by their
    sigma, '.' and '?' are not counted."""
    n_decimal = None
    for value in l_value:
        if not isinstance(value, str) or value in (".", "?") or "(" in value:
            continue
        if "e" in value.lower():
            return None
        n_value = len(value.split(".")[1]) if "." in value else 0
        n_decimal = n_value if n_decimal is None else max(n_decimal, n_value)
    return None if n_decimal is None else f"{{:.{n_decimal:}f}}"


def value_error_mark_from_cif(string: str, cif_name: str):
    """Value, sigma and constraint mark of a refinable rcif value.

    "." and "?" give None. A value that is not a number raises
    RcifFormatError naming the tag, instead of being silently dropped.
    """
    try:
        value, error, mark = string_to_value_error_mark(string)
    except ValueError as error_value:
        raise RcifFormatError(
            f"Value '{string:}' of '{cif_name:}' is not a number.") from error_value
    if value is None and string.strip() not in (".", "?"):
        raise RcifFormatError(f"Value '{string:}' of '{cif_name:}' is not a number.")
    return value, error, mark


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

        ls_html.append("<h3>Attributes and methods</h3>")
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
        # a special name (__setstate__ asked by copy.deepcopy and pickle on
        # an object made without __init__) is never a CIF attribute: without
        # this check the lookups below come here again until RecursionError,
        # a stack overflow in a thread (CrysPy Editor closed on macOS)
        if name.startswith("__") and name.endswith("__"):
            raise AttributeError(name)
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

    def __init_subclass__(cls, **kwargs):
        """Every attribute of the ATTR_* tuples of a new item class becomes
        a descriptor of the class (cif_attributes.py): the type, the
        allowed values and the limits are checked when it is set."""
        super().__init_subclass__(**kwargs)
        install_cif_descriptors(cls)

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
        pass

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
        if cls is ItemN:  # a category unknown to CrysPy: a class for it
            return take_item_class_of_cif_data(cif_data).from_cif(string)
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
                            value, error, mark = value_error_mark_from_cif(
                                cif_data[_cif_attr_full].value, _cif_attr_full)
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
            raise RcifFileNotFoundError(f"File {f_name:} is not found.")
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
    
    def get_dictionary(self):
        res = {}
        return res
    
    def take_parameters_from_dictionary(self, ddict_diffrn, l_parameter_name: list=None, l_sigma: list=None):
        """
        """
        pass
        return None

# classes of the CIF/STAR categories unknown to CrysPy, made at first use
GENERIC_ITEM_CLASSES = {}


class GenericCategoryType(type):
    """Type of the classes made for the categories unknown to CrysPy
    (make_item_class, make_loop_class). Such a class cannot be imported by
    its name, so pickle saves its description and makes it again
    (copyreg below): the objects are pickled (multiprocessing) and copied
    as the objects of the other classes."""


def reduce_generic_class(cls):
    """(function, arguments) making the class again."""
    return cls.__dict__["REBUILD"]


copyreg.pickle(GenericCategoryType, reduce_generic_class)


def is_generic(obj) -> bool:
    """True for an item or loop (object or class) of a CIF/STAR category
    that CrysPy does not know (made by make_item_class, make_loop_class)."""
    cls = obj if isinstance(obj, type) else type(obj)
    return bool(cls.__dict__.get("IS_GENERIC", False))


def make_item_class(prefix: str, attr_cif: tuple, attr_names: tuple = None,
                    attr_ref: tuple = (), attr_int: tuple = (),
                    formats: tuple = ()) -> type:
    """Class derived from ItemN for a CIF/STAR category that CrysPy does
    not know.

    Parameters
    ----------
    prefix : str
        Category, without the leading '_' ('publ_author').
    attr_cif : tuple
        Names of the attributes in the CIF tags ('name' of _publ_author_name).
    attr_names : tuple
        Names of the attributes in Python; attr_cif when None.
    attr_ref : tuple
        Attributes with numerical values (refinable, float).
    attr_int : tuple
        Attributes with integer values (int, not refinable); the others are
        strings.
    formats : tuple
        Pairs (attribute, format) for writing the numbers (D_FORMATS).

    Returns
    -------
    type
        The class, the same object for the same arguments; it is described
        by the ATTR_* tuples as every item class (its attributes are
        descriptors) and is pickled by this description.
    """
    attr_cif = tuple(attr_cif)
    attr_names = attr_cif if attr_names is None else tuple(attr_names)
    attr_ref = tuple(name for name in attr_names if name in attr_ref)
    attr_int = tuple(name for name in attr_names if name in attr_int and name not in attr_ref)
    formats = tuple((name, s_format) for name, s_format in formats
                    if name in attr_ref and s_format is not None)
    key = (prefix, attr_cif, attr_names, attr_ref, attr_int, formats)
    item_class = GENERIC_ITEM_CLASSES.get(key)
    if item_class is not None:
        return item_class
    attr_types = tuple(float if name in attr_ref else (int if name in attr_int else str)
                       for name in attr_names)
    namespace = {
        "__doc__": f"Item of the CIF category _{prefix:} (not described in CrysPy).",
        "__module__": __name__,
        "IS_GENERIC": True, "GENERIC_KEY": key,
        "REBUILD": (make_item_class, key),
        "PREFIX": prefix,
        "ATTR_MANDATORY_NAMES": (), "ATTR_MANDATORY_TYPES": (), "ATTR_MANDATORY_CIF": (),
        "ATTR_OPTIONAL_NAMES": attr_names, "ATTR_OPTIONAL_TYPES": attr_types,
        "ATTR_OPTIONAL_CIF": attr_cif,
        "ATTR_NAMES": attr_names, "ATTR_TYPES": attr_types, "ATTR_CIF": attr_cif,
        "ATTR_INT_NAMES": (), "ATTR_INT_PROTECTED_NAMES": (),
        "ATTR_REF": attr_ref,
        "ATTR_SIGMA": tuple(f"{name:}_sigma" for name in attr_ref),
        "ATTR_CONSTR_FLAG": tuple(f"{name:}_constraint" for name in attr_ref),
        "ATTR_REF_FLAG": tuple(f"{name:}_refinement" for name in attr_ref),
        "ATTR_CONSTR_MARK": tuple(f"{name:}_mark" for name in attr_ref),
        "D_FORMATS": dict(formats), "D_CONSTRAINTS": {}, "D_DEFAULT": {}, "D_MIN": {}, "D_MAX": {},
    }
    class_name = "ItemN_" + "".join(char if char.isalnum() else "_" for char in prefix)
    item_class = GenericCategoryType(class_name, (ItemN,), namespace)
    GENERIC_ITEM_CLASSES[key] = item_class
    return item_class


def take_item_class_of_cif_data(cif_data) -> type:
    """Class of the items of cif_data (pycifstar Data) for a category that
    CrysPy does not know: the category is the common prefix of the tags
    (the first part of a single tag), a value that is a number makes its
    attribute refinable."""
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
    d_type = {attr: take_value_type([cif_data[name].value]) for attr, name in zip(l_attr, l_name)}
    return make_item_class(
        prefix.strip("_").strip("."), tuple(l_attr),
        attr_ref=tuple(attr for attr in l_attr if d_type[attr] is float),
        attr_int=tuple(attr for attr in l_attr if d_type[attr] is int),
        formats=tuple((attr, take_float_format([cif_data[name].value]))
                      for attr, name in zip(l_attr, l_name)))
