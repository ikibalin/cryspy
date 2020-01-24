from pycifstar import Data
from pycifstar.item import print_string
from typing import List, Tuple
import warnings

def val_to_str(val):
    if isinstance(val, str):
        s_val = print_string(val)
    elif val is None:
        s_val = "."
    else:
        s_val = f"{val:}"
    return s_val

class ItemConstr(object):
    def __init__(self, mandatory_attribute = (), optional_attribute = (), internal_attribute = (), prefix=""):
        super(ItemConstr, self).__init__()
        self.__mandatory_attribute = mandatory_attribute
        self.__optional_attribute = optional_attribute
        self.__internal_attribute = internal_attribute
        self.__prefix = prefix
        self.clean_attribute
        self.__flag_renewed = True

    def __repr__(self) -> str:
        ls_out = []
        ls_out.append("ItemConstr: ")
        ls_out.append(f"{str(self):}")
        return "\n".join(ls_out)

    def __str__(self) -> str:
        ls_out = []
        for _attr in self.__mandatory_attribute:
            s_attr = f"__{_attr:}"
            _val = getattr(self, s_attr)
            if _val is not None:
                s_val = val_to_str(_val)
                ls_out.append(f"{_attr:}: {s_val:}")
        for _attr in self.__optional_attribute:
            s_attr = f"__{_attr:}"
            _val = getattr(self, s_attr)
            if _val is not None:
                s_val = val_to_str(_val)
                ls_out.append(f"{_attr:}: {s_val:}")
        for _attr in self.__internal_attribute:
            s_attr = f"__{_attr:}"
            _val = getattr(self, s_attr)
            if _val is not None:
                s_val = val_to_str(_val)
                ls_out.append(f"{_attr:}: {s_val:}")
        return "\n".join(ls_out)

    def to_cif(self, separator="_", flag=False) -> str: 
        """
Save information about object in string in STAR format

Args:
    prefix: prefix in front of label of attribute
    separator: separator between prefix and attribute ("_" or ".")
    flag: for undefined attribute "." will be printed

Returns:
    A string in STAR/CIF format
        """
        prefix = self.prefix
        ls_out = []

        mandatory_attribute = self.MANDATORY_ATTRIBUTE
        optional_attribute = self.OPTIONAL_ATTRIBUTE
        try:
            related_cif_mandatory_attribute = self.RELATED_CIF_MANDATORY_ATTRIBUTE
            related_cif_optional_attribute = self.RELATED_CIF_OPTIONAL_ATTRIBUTE
        except:
            related_cif_mandatory_attribute = mandatory_attribute
            related_cif_optional_attribute = optional_attribute
        for _attr, cif_attr in zip(mandatory_attribute, related_cif_mandatory_attribute):
            #s_attr = f"__{_attr:}"
            _val = getattr(self, _attr)
            if _val is not None:
                s_val = val_to_str(_val)
                ls_out.append(f"_{prefix:}{separator:}{cif_attr:} {s_val:}")
            elif flag:
                ls_out.append(f"_{prefix:}{separator:}{cif_attr:} .")
        for _attr, cif_attr in zip(optional_attribute, related_cif_optional_attribute):
            #s_attr = f"__{_attr:}"
            _val = getattr(self, _attr)
            if _val is not None:
                s_val = val_to_str(_val)
                ls_out.append(f"_{prefix:}{separator:}{cif_attr:} {s_val:}")
            elif flag:
                ls_out.append(f"_{prefix:}{separator:}{cif_attr:} .")
        return "\n".join(ls_out)
    def print_attribute(self, l_attr=()) -> str:
        """
        Save attributes in one string
        Args:
            l_attr is a list/tuple of attributes to print
        Returns:
            A string of values (if None it will be ".")
        """
        ls_out = []
        for _attr in l_attr:
            s_attr = f"__{_attr:}"
            _val = getattr(self, s_attr)
            if _val is not None:
                s_val = val_to_str(_val)
                ls_out.append(f"{s_val:}")
            else:
                ls_out.append(f".")
        return " ".join(ls_out)

    @property
    def mandatory_attribute(self) -> Tuple[str]:
        return self.__mandatory_attribute
    @property
    def optional_attribute(self) -> Tuple[str]:
        return self.__optional_attribute
    @property
    def internal_attribute(self) -> Tuple[str]:
        return self.__internal_attribute
    @property
    def prefix(self) -> str:
        return self.__prefix

    @property
    def flag_renewed(self)->bool:
        return self.__flag_renewed

    @flag_renewed.setter
    def flag_renewed(self, x:bool):
        self.__flag_renewed = bool(x)

    @property
    def clean_attribute(self):
        for _ in self.__mandatory_attribute:
            setattr(self, f"__{_:}", None)
        for _ in self.__optional_attribute:
            setattr(self, f"__{_:}", None)
        for _ in self.__internal_attribute:
            setattr(self, f"__{_:}", None)

    def is_attribute_mandatory(self, attr:str) -> bool:
        flag = False
        if attr in self.__mandatory_attribute:
            flag = True
        return flag

    def is_attribute_optional(self, attr:str) -> bool:
        flag = False
        if attr in self.__optional_attribute:
            flag = True
        return flag

    def is_attribute_internal(self, attr:str) -> bool:
        flag = False
        if attr in self.__internal_attribute:
            flag = True
        return flag

    def is_attribute(self, attr:str) -> bool:
        flag_1 = self.is_attribute_mandatory(attr)
        flag_2 = self.is_attribute_optional(attr)
        flag_3 = self.is_attribute_internal(attr)
        return (flag_1 | flag_2 | flag_3)

    @property
    def is_defined(self) -> bool:
        flag = True
        for _ in self.mandatory_attribute:
            s_attr = f"__{_:}"
            val = getattr(self, s_attr)
            if val is None:
                flag = False
        return flag

    @property
    def form_object(self)->bool:
        return True

    def is_defined_attribute(self, _attr:str) -> bool:
        s_attr = f"__{_attr:}"
        flag = hasattr(self, s_attr)
        if flag:
            val = getattr(self, s_attr)
            if val is None:
                flag = False
        return flag

    @property
    def is_variable(self) -> bool:
        """
Output: True if there is any refined parameter
        """
        return False

    def get_variables(self) -> List:
        """
Output: the list of the refined parameters
        """
        return []

    @classmethod
    def from_cif(cls, string: str):
        cif_data = Data()
        flag = cif_data.take_from_string(string)
        mandatory_attribute = cls.MANDATORY_ATTRIBUTE
        optional_attribute = cls.OPTIONAL_ATTRIBUTE
        try:
            related_cif_mandatory_attribute = cls.RELATED_CIF_MANDATORY_ATTRIBUTE
            related_cif_optional_attribute = cls.RELATED_CIF_OPTIONAL_ATTRIBUTE
        except:
            related_cif_mandatory_attribute = mandatory_attribute
            related_cif_optional_attribute = optional_attribute
        prefix = cls.PREFIX
        separator = "_"
        flag = all([cif_data.is_value("_"+prefix+separator+_attr) for _attr in related_cif_mandatory_attribute])
        if not (flag):
            separator = "."
            flag = all([cif_data.is_value("_"+prefix+separator+_attr) for _attr in related_cif_mandatory_attribute])
        if flag:
            l_attr = list(mandatory_attribute)+list(optional_attribute)
            l_cif_attr = list(related_cif_mandatory_attribute)+list(related_cif_optional_attribute)
            separator = "_"
            flag_2 = any([cif_data.is_value("_"+prefix+separator+_attr) for _attr in l_cif_attr])
            if not(flag_2):
                separator = "."
                flag_2 = any([cif_data.is_value("_"+prefix+separator+_attr) for _attr in l_cif_attr])
            if flag_2:
                _item = cls()
                for _attr, _cif_attr in zip(l_attr, l_cif_attr):
                    _cif_attr_full = "_"+prefix+separator+_cif_attr
                    if cif_data.is_value(_cif_attr_full):
                        setattr(_item, _attr, cif_data[_cif_attr_full].value)
            else:
                _item = None
        else:
            _item = None
        if _item is not None:
            if _item.is_defined:
                _item.form_object
        return _item        

    def _show_message(self, s_out: str):
        warnings.warn("***  Error ***\n"+s_out, UserWarning, stacklevel=2)


