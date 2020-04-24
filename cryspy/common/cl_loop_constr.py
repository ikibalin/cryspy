import warnings
from cryspy.common.cl_item_constr import ItemConstr
from typing import List, Tuple
from pycifstar import Data


class LoopConstr(object):
    def __init__(self, category_key=(), item_class=ItemConstr, loop_name="", internal_attribute=()):
        super(LoopConstr, self).__init__()
        self.__item = []
        self.__category_key = category_key
        self.__item_class = item_class
        self.__internal_attribute_loop = internal_attribute
        self.clean_attribute

        self.__mandatory_attribute = item_class.MANDATORY_ATTRIBUTE
        self.__optional_attribute = item_class.OPTIONAL_ATTRIBUTE
        self.__internal_attribute = item_class.INTERNAL_ATTRIBUTE

        try:
            self.__related_cif_mandatory_attribute = item_class.RELATED_CIF_MANDATORY_ATTRIBUTE
            self.__related_cif_optional_attribute = item_class.RELATED_CIF_OPTIONAL_ATTRIBUTE
        except:
            self.__related_cif_mandatory_attribute = item_class.MANDATORY_ATTRIBUTE
            self.__related_cif_optional_attribute = item_class.OPTIONAL_ATTRIBUTE
        self.__loop_name = ""

        self.loop_name = loop_name

    @property
    def clean_attribute(self):
        for _ in self.__internal_attribute_loop:
            setattr(self, f"__{_:}", None)

    def __repr__(self) -> str:
        ls_out = [f"{type(self).__name__:}: ", f"{str(self):}"]
        return "\n".join(ls_out)

    def __str__(self) -> str:
        ls_out = []
        ls_out.append(f"prefix: {self.prefix:}")
        ls_out.append(f"category_key: {', '.join(self.__category_key):}")
        ls_out.append(f"mandatory_attribute: {', '.join(self.__mandatory_attribute):}")
        ls_out.append(f"optional_attribute: {', '.join(self.__optional_attribute):}")
        for _item in self.item:
            ls_out.append(" ".join(
                [str(getattr(_item, _attr)) if (getattr(_item, _attr) is not None) else "." for _attr in
                 self.__mandatory_attribute]) + "  " +
                          " ".join(
                              [str(getattr(_item, _attr)) if (getattr(_item, _attr) is not None) else "." for _attr in
                               self.__optional_attribute]))

        return "\n".join(ls_out)

    def to_cif(self, separator="_", flag=False, flag_minimal=True) -> str:
        """
Print information about object in string in STAR format

Args:
    prefix is a prefix in front of label of attribute
    separator is a separator between prefix and attribute ("_" or ".")
    flag if it's True the value "." will be printed for undefined attributes
    flag_minimal if it's True the minimal set of object will be printed

Returns:
    A string in STAR/CIF format
        """
        prefix = self.prefix
        if (not (prefix.startswith("_")) & (prefix != "")):
            prefix = f"_{prefix:}"
        ls_out = []
        l_attr_print = list(self.__mandatory_attribute) + list(self.__optional_attribute)
        l_cif_attr = list(self.__related_cif_mandatory_attribute) + list(self.__related_cif_optional_attribute)
        if not (flag):
            l_flag = [all([_item.is_defined_attribute(_attr) for _item in self.item]) for _attr in l_attr_print]
            l_attr_print = [_ for _, _flag in zip(l_attr_print, l_flag) if _flag]
            l_cif_attr = [_ for _, _flag in zip(l_cif_attr, l_flag) if _flag]
        ls_out.append(f"loop_{self.loop_name}")
        ls_out.append("\n".join([f"{prefix:}{separator:}{_attr:}" for _attr in l_cif_attr]))
        for _item in self.item:
            ls_out.append(_item.print_attribute(l_attr_print))
        return "\n".join(ls_out)

    @property
    def is_defined(self) -> bool:
        flag_1 = True
        flag_2 = all([_obj.is_defined for _obj in self.item])
        flag = flag_1 & flag_2
        return flag

    @property
    def form_object(self) -> bool:
        flag_1 = True
        flag_2 = all([_obj.form_object for _obj in self.item if _obj.is_defined])
        flag = flag_1 & flag_2
        return flag

    @classmethod
    def from_cif(cls, string: str) -> List:
        cif_data = Data()
        flag = cif_data.take_from_string(string)
        _item_class = cls.ITEM_CLASS
        l_attr = list(_item_class.MANDATORY_ATTRIBUTE) + list(_item_class.OPTIONAL_ATTRIBUTE)
        try:
            l_cif_attr = list(_item_class.RELATED_CIF_MANDATORY_ATTRIBUTE) + list(
                _item_class.RELATED_CIF_OPTIONAL_ATTRIBUTE)
        except:
            l_cif_attr = l_attr
        l_cif_attr = [_.lower() for _ in l_cif_attr]
        prefix = _item_class.PREFIX
        l_obj = []
        for cif_loop in cif_data.loops:
            if ("_" + prefix.lower()) == cif_loop.prefix.lower():
                flag = False
                prefix = cif_loop.prefix
                l_name = cif_loop.names
                loop_name = cif_loop.name
                l_name_short = [_[(len(prefix) + 1):] for _ in l_name]

                _obj = cls(loop_name=loop_name)
                _i = 0
                for _name, _name_short in zip(l_name, l_name_short):
                    if _name_short.lower() in l_cif_attr:
                        _name_short_obj = l_attr[l_cif_attr.index(_name_short)]
                    else:
                        _name_short_obj = _name_short.lower()
                    if _i == 0:
                        item = []
                        for _val in cif_loop[_name]:
                            _item = _item_class()
                            setattr(_item, _name_short_obj, _val)
                            item.append(_item)
                    else:
                        for _val, _item in zip(cif_loop[_name], item):
                            setattr(_item, _name_short_obj, _val)
                    _i += 1
                _obj.item = item
                l_obj.append(_obj)
        return l_obj

    @property
    def category_key(self) -> Tuple[str]:
        return self.__category_key

    @property
    def prefix(self) -> str:
        return self.__item_class.PREFIX

    @property
    def loop_name(self) -> str:
        return getattr(self, "__loop_name")

    @loop_name.setter
    def loop_name(self, x) -> str:
        if x is None:
            x_in = ""
        else:
            x_in = str(x)
        setattr(self, "__loop_name", x_in)

    @property
    def item_class(self):
        return self.__item_class

    @property
    def items(self):
        _attrs = self.category_key
        l_id = [tuple([getattr(_, _attr) for _attr in _attrs]) for _ in self.__item]
        return l_id

    def is_item(self, _id) -> bool:
        flag = False
        l_id = self.items
        if _id in l_id: flag = True
        return flag

    def __getitem__(self, *args):
        _id = tuple(args)
        _attr = self.category_key
        if self.is_item(_id):
            _ind = self.items.index(_id)
            return self.__item[_ind]
        else:
            s_out = " ".join([f"{_:}" for _ in _id])
            warnings.warn(f"Element '{s_out:}' is not found in list", UserWarning, stacklevel=2)
            return None

    @property
    def item(self):
        return self.__item

    @item.setter
    def item(self, l_x):
        flag = False
        if all([isinstance(x, self.item_class) for x in l_x]):
            flag = True
        if flag:
            self.__item = l_x
            self.form_object
        else:
            warnings.warn(f"One of the introduced object does not correspond to class {self.item_class:}", UserWarning,
                          stacklevel=2)

    def __getattr__(self, attr):
        if attr in self.__mandatory_attribute:
            res = [getattr(_item, attr) for _item in self.__item]
        elif attr in self.__optional_attribute:
            res = [getattr(_item, attr) for _item in self.__item]
        elif attr in self.__internal_attribute:
            res = [getattr(_item, attr) for _item in self.__item]
        elif attr in self.__internal_attribute_loop:
            res = getattr(self, f"__{attr:}")
        else:
            res = None
            warnings.warn(f"Attribute '{attr:}' is not defined", UserWarning, stacklevel=2)
        return res

    @property
    def is_variable(self) -> bool:
        """
Output: True if there is any refined parameter
        """
        return any([_item.is_variable for _item in self.__item])

    def get_variables(self) -> List:
        """
Output: the list of the refined parameters
        """
        l_res = []
        for _item in self.__item:
            l_res.extend(_item.get_variables())
        return l_res

    def _show_message(self, s_out: str):
        warnings.warn("***  Error ***\n" + s_out, UserWarning, stacklevel=2)
