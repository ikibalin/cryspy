"""Parent class DataN."""
import os
import os.path
from warnings import warn
from typing import Union, NoReturn
from pycifstar import Global, to_global

from cryspy.A_functions_base.function_1_strings import ciftext_to_html

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN
from cryspy.B_parent_classes.cl_3_data import DataN


class GlobalN(object):
    """GlobalN class."""

    def __repr__(self):
        """
        Magic method print() is redefined.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        ls_out = [self.get_name()]
        for item in self.items:
            if isinstance(item, LoopN):
                ls_out.append(f"{4*' ':}{item.get_name():} (loop)")
            else:
                ls_out.append(f"{4*' ':}{item.get_name():}")
            if isinstance(item, DataN):
                for i_d in item.items:
                    if isinstance(i_d, ItemN):
                        ls_out.append(f"{8*' ':}{i_d.get_name():}")
                    else:
                        ls_out.append(
                            f"{8*' ':}{i_d.get_name():} (loop)")
        return "\n".join(ls_out)

    def _repr_html_(self):
        """Representation in HTML format."""
        return ciftext_to_html(self.__repr__())

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

    def add_items(self, items: list):
        """Add items."""
        l_name = [item.get_name() for item in items]
        s_name = set(l_name)
        if len(s_name) != l_name:
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

        prefix_g, prefix_n = name[0], name[1]

        if prefix_g != (prefix, global_name):
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
                warn(f"Item {item.PREFIX:} is not fully described.",
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
        ls_out.extend(l_s_datan)
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
