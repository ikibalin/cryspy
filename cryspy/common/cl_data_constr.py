import warnings
from cryspy.common.cl_item_constr import ItemConstr
from cryspy.common.cl_loop_constr import LoopConstr
from typing import List, Tuple
from pycifstar import Data


class DataConstr(object):
    def __init__(self, data_name="",
                 mandatory_classes = (),
                 optional_classes = (),
                 internal_classes = ()):
        super(DataConstr, self).__init__()
        setattr(self, "__data_name", data_name)
        setattr(self, "__mandatory_classes", mandatory_classes)
        setattr(self, "__optional_classes", optional_classes)
        setattr(self, "__internal_classes", internal_classes)

        setattr(self, "__mandatory_objs", [])
        setattr(self, "__optional_objs", [])
        setattr(self, "__internal_objs", [])

    def __repr__(self) -> str:
        ls_out = []
        ls_out.append("DataConstr: ")
        ls_out.append(f"{str(self):}")
        return "\n".join(ls_out)

    def __str__(self) -> str:
        ls_out = []
        ls_out.append("Mandatory classes:")
        ls_out.extend([str(_) for _ in self.mandatory_classes])
        ls_out.append("\nOptional classes:")
        ls_out.extend([str(_) for _ in self.optional_classes])
        ls_out.append("\nMandatory objects:")
        ls_out.extend([str(_)+"\n" for _ in self.mandatory_objs])
        ls_out.append("\nOptional objects:")
        ls_out.extend([str(_)+"\n" for _ in self.optional_objs])
        ls_out.append("\nInternal objects:")
        ls_out.extend([str(_)+"\n" for _ in self.internal_objs])
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
        ls_out = []
        ls_out.append(f"data_{self.data_name:}\n")
        ls_out.extend([_.to_cif(separator=separator, flag=flag)+"\n" for _ in self.mandatory_objs])
        ls_out.extend([_.to_cif(separator=separator, flag=flag)+"\n" for _ in self.optional_objs])
        return "\n".join(ls_out)

    @classmethod
    def from_cif(cls, string: str):
        cif_data = Data()
        flag = cif_data.take_from_string(string)

        cif_items = cif_data.items
        cif_loops = cif_data.loops

        mandatory_objs = []
        optional_objs = []

        flag = True
        for _cls in cls.MANDATORY_CLASSES:
            flag = False
            if issubclass(_cls, ItemConstr):
                prefix_cls = _cls.PREFIX
                if cif_items.is_prefix(prefix_cls):
                    cif_items_prefix = cif_items[prefix_cls]
                    cif_string = str(cif_items_prefix)
                    _obj_prefix = _cls.from_cif(cif_string)
                    mandatory_objs.append(_obj_prefix)
                    flag = True
            elif issubclass(_cls, LoopConstr):
                prefix_cls = _cls.ITEM_CLASS.PREFIX
                for cif_loop in cif_loops:
                    if cif_loop.is_prefix("_"+prefix_cls):
                        cif_string = str(cif_loop)
                        _obj_prefix = _cls.from_cif(cif_string)
                        mandatory_objs.extend(_obj_prefix)
                        flag = True
            if not(flag):
                #warnings.warn(f"unknown class type : '{_cls:}'", UserWarning, stacklevel=2)
                break
        
        if not(flag):
            return None

        for _cls in cls.OPTIONAL_CLASSES:
            if issubclass(_cls, ItemConstr):
                prefix_cls = _cls.PREFIX
                if cif_items.is_prefix(prefix_cls):
                    cif_items_prefix = cif_items[prefix_cls]
                    cif_string = str(cif_items_prefix)
                    _obj_prefix = _cls.from_cif(cif_string)
                    optional_objs.append(_obj_prefix)
            elif issubclass(_cls, LoopConstr):
                prefix_cls = _cls.ITEM_CLASS.PREFIX
                for cif_loop in cif_loops:
                    if cif_loop.is_prefix("_"+prefix_cls):
                        cif_string = str(cif_loop)
                        _obj_prefix = _cls.from_cif(cif_string)
                        optional_objs.extend(_obj_prefix)
            else:
                warnings.warn(f"unknown class type : '{_cls:}'", UserWarning, stacklevel=2)
        data_name = cif_data.name
        _obj = cls(data_name=data_name)
        _obj.mandatory_objs = mandatory_objs
        _obj.optional_objs = optional_objs
        
        return _obj


    @property
    def data_name(self) -> str:
        return getattr(self, "__data_name")
    @data_name.setter
    def data_name(self, x) -> str:
        if x is None:
            x_in = ""
        else:
            x_in = str(x)
        setattr(self, "__data_name", x_in)

    @property
    def mandatory_classes(self):
        return getattr(self, "__mandatory_classes")
    @property
    def optional_classes(self):
        return getattr(self, "__optional_classes")

    @property
    def mandatory_objs(self):
        return getattr(self, "__mandatory_objs")
    @mandatory_objs.setter
    def mandatory_objs(self, l_x: List) -> str:
        l_x_in = []
        for x in l_x:
            if any([isinstance(x, _cls) for _cls in self.mandatory_classes]):
                l_x_in.append(x)
        return setattr(self, "__mandatory_objs", l_x_in)

    @property
    def optional_objs(self) -> str:
        return getattr(self, "__optional_objs")
    @optional_objs.setter
    def optional_objs(self, l_x: List) -> str:
        l_x_in = []
        for x in l_x:
            if any([isinstance(x, _cls) for _cls in self.optional_classes]):
                l_x_in.append(x)
        return setattr(self, "__optional_objs", l_x_in)

    @property
    def internal_objs(self) -> str:
        return getattr(self, "__internal_objs")

    def is_class(self,  _cls) -> bool:
        flag = any([isinstance(_obj, _cls) for _obj in 
                  (self.mandatory_objs + self.optional_objs + self.internal_objs)])
        return flag

    def __getitem__(self, _cls):
        l_res = []
        if self.is_class(_cls):
            for _obj in (self.mandatory_objs + self.optional_objs + self.internal_objs):
                if isinstance(_obj, _cls):
                    l_res.append(_obj)
        return tuple(l_res)
    @property
    def is_defined(self):
        flag_1 = all([self.is_class(_cls) for _cls in getattr(self, "mandatory_classes")])
        flag_2 = all([_obj.is_defined for _obj in getattr(self, "mandatory_objs")])
        flag = flag_1 & flag_2
        return flag
    @property
    def form_object(self):
        return True

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


    #def __getattr__(self, attr):
    #    if attr in self.__mandatory_attribute:
    #        res = [getattr(_item, attr) for _item in self.__item]
    #    elif attr in self.__optional_attribute:
    #        res = [getattr(_item, attr) for _item in self.__item]
    #    elif attr in self.__internal_attribute:
    #        res = [getattr(_item, attr) for _item in self.__item]
    #    else:
    #        res = None
    #        print(f"Attribute '{attr:}' is not defined")
    #    return res

