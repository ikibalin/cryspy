from cryspy.f_common.cl_item_constr import ItemConstr
from typing import List, Tuple


class LoopConstr(object):
    def __init__(self, category_key = (), item_class=ItemConstr, prefix="", label=""):
        super(LoopConstr, self).__init__()
        self.__item = []
        self.__category_key = category_key   
        self.__item_class = item_class
        self.__prefix = prefix
        self.__mandatory_attribute = item_class.MANDATORY_ATTRIBUTE
        self.__optional_attribute = item_class.OPTIONAL_ATTRIBUTE
        self.__internal_attribute = item_class.INTERNAL_ATTRIBUTE
        self.__label = ""

        self.label = label

    def __repr__(self) -> str:
        ls_out = []
        ls_out.append("LoopConstr: ")
        ls_out.append(f"{str(self):}")
        return "\n".join(ls_out)

    def __str__(self) -> str:
        ls_out = []
        ls_out.append(f"prefix: {', '.join(self.__prefix):}")
        ls_out.append(f"category_key: {', '.join(self.__category_key):}")
        ls_out.append(f"mandatory_attribute: {', '.join(self.__mandatory_attribute):}")
        ls_out.append(f"optional_attribute: {', '.join(self.__optional_attribute):}")
        for _item in self.item:
            ls_out.append(" ".join([str(getattr(_item, _attr)) for _attr in self.__mandatory_attribute]))
            ls_out.append(" ".join([str(getattr(_item, _attr)) for _attr in self.__optional_attribute]))

        return "\n".join(ls_out)

    def to_cif(self, separator="_", flag=False) -> str: 
        """
        Print information about object in string in STAR format

        Args:
            prefix: prefix in front of label of attribute
            separator: separator between prefix and attribute ("_" or ".")
            flag: for undefined attribute "." will be printed

        Returns:
            A string in STAR/CIF format
        """
        prefix = self.prefix
        if (not(prefix.startswith("_")) & (prefix != "")):
            prefix = f"_{prefix:}"
        ls_out = []
        l_attr_print = list(self.__mandatory_attribute)+list(self.__optional_attribute)
        if not(flag):
            l_flag = [all([_item.is_defined_attribute(_attr) for _item in self.item]) for _attr in l_attr_print]
            l_attr_print = [_ for _, _flag in zip(l_attr_print, l_flag) if _flag]
        ls_out.append(f"loop_{self.label}")
        ls_out.append("\n".join([f"{prefix:}{separator:}{_attr:}" for _attr in l_attr_print]))
        for _item in self.item:
            ls_out.append(_item.print_attribute(l_attr_print))
        return "\n".join(ls_out)


    @property
    def category_key(self) -> Tuple[str]:
        return self.__category_key

    @property
    def prefix(self) -> str:
        return self.__prefix

    @property
    def label(self) -> str:
        return getattr(self, "__label")
    @label.setter
    def label(self, x) -> str:
        if x is None:
            x_in = ""
        else:
            x_in = str(x)
        setattr(self, "__label", x_in)

    @property
    def item_class(self):
        return self.__item_class
        

    @property
    def items(self):
        _attrs = self.category_key
        l_id = [tuple([getattr(_,_attr) for _attr in _attrs]) for _ in self.__item]
        return l_id

    def is_item(self,  _id) -> bool:
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
            print(f"Element '{s_out:}' is not found in list")
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
        else:
            print(f"One of the introduced object does not correspond to class {self.item_class:}")

    def __getattr__(self, attr):
        if attr in self.__mandatory_attribute:
            res = [getattr(_item, attr) for _item in self.__item]
        elif attr in self.__optional_attribute:
            res = [getattr(_item, attr) for _item in self.__item]
        elif attr in self.__internal_attribute:
            res = [getattr(_item, attr) for _item in self.__item]
        else:
            res = None
            print(f"Attribute '{attr:}' is not defined")
        return res

