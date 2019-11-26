
from typing import List, Tuple

class ItemConstr(object):
    def __init__(self, mandatory_attribute = (), optional_attribute = (), internal_attribute = ()):
        super(ItemConstr, self).__init__()
        self.__mandatory_attribute = mandatory_attribute
        self.__optional_attribute = optional_attribute
        self.__internal_attribute = internal_attribute
        self.clean_attribute

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
                ls_out.append(f"{_attr:}: {_val:}")
        for _attr in self.__optional_attribute:
            s_attr = f"__{_attr:}"
            _val = getattr(self, s_attr)
            if _val is not None:
                ls_out.append(f"{_attr:}: {_val:}")
        for _attr in self.__internal_attribute:
            s_attr = f"__{_attr:}"
            _val = getattr(self, s_attr)
            if _val is not None:
                ls_out.append(f"{_attr:}: {_val:}")
        return "\n".join(ls_out)

    def to_cif(self, prefix="", separator="_", flag=False) -> str: 
        """
        Save information about object in string in STAR format

        Args:
            prefix: prefix in front of label of attribute
            separator: separator between prefix and attribute ("_" or ".")
            flag: for undefined attribute "." will be printed

        Returns:
            A string in STAR/CIF format
        """
        ls_out = []
        for _attr in self.__mandatory_attribute:
            s_attr = f"__{_attr:}"
            _val = getattr(self, s_attr)
            if _val is not None:
                ls_out.append(f"{prefix:}{separator:}{_attr:} {_val:}")
            elif flag:
                ls_out.append(f"{prefix:}{separator:}{_attr:} .")
        for _attr in self.__optional_attribute:
            s_attr = f"__{_attr:}"
            _val = getattr(self, s_attr)
            if _val is not None:
                ls_out.append(f"{prefix:}{separator:}{_attr:} {_val:}")
            elif flag:
                ls_out.append(f"{prefix:}{separator:}{_attr:} .")
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
                ls_out.append(f"{_val:}")
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
        flag_1 = self.is_attribute_mandatory
        flag_2 = self.is_attribute_optional
        flag_3 = self.is_attribute_internal
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

    def is_defined_attribute(self, _attr:str) -> bool:
        s_attr = f"__{_attr:}"
        flag = hasattr(self, s_attr)
        if flag:
            val = getattr(self, s_attr)
            if val is None:
                flag = False
        return flag
