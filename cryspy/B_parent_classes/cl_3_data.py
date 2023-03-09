"""Parent class DataN."""
import os
import os.path
from warnings import warn
from typing import Union, NoReturn
from pycifstar import Data, to_data

from cryspy.A_functions_base.function_1_markdown import md_to_html
from cryspy.A_functions_base.function_1_objects import \
    get_functions_of_objet, get_table_html_for_variables

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN


class DataN(object):
    """Data container of loops and items."""

    def __repr__(self):
        """
        Magic method print() is redefined.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        ls_out = [f"# Object '{self.get_name():}'"]
        for item in self.items:
            if isinstance(item, ItemN):
                ls_out.append(f"{4*' ':}.{item.get_name():}")
            else:
                ls_out.append(f"{4*' ':}.{item.get_name():} (loop)")
                
        method = self.methods_html()
        if method != "":
            ls_out.append(f"\n# Methods:\n{method:}\n")                
        return "\n".join(ls_out)

    def _repr_html_(self):
        """Representation in HTML format."""
        ls_html = [f"<h2>Object '{self.get_name():}'</h2>"]
        ls_html.append(self.attributes_to_html())


        ls_html.append(get_table_html_for_variables(self))

        report = self.report_html()
        if report != "":
            ls_html.append(f"<h2>Description </h2> {report:}")

        ls_html.append(f"<h2>Classes and methods</h2>")
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

    def methods_html(self):
        ls_html = [f".{func_name}" for func_name in
                   get_functions_of_objet(self)]
        return ", ".join(ls_html)+"."
    
    def attributes_to_html(self) -> str:
        """Representation of defined parameters in HTML format.
        """
        ls_html = ["<table>"]
        ls_html.append("<tr><th>Attribute</th><th> Note </th></tr>")
        items_sorted = sorted(self.items, key=lambda item: item.get_name())
        for item in items_sorted:
            item_type = item.__doc__.strip().split("\n")[0]
            ls_html.append(f"<tr><td>.{item.get_name():}</td>\
<td>{item_type:}</td></tr>")
        ls_html.append("</table>")
        return " ".join(ls_html)

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
        for item in self.__getattribute__('items'): # self.items gives an error during the parallelization process
            if name.lower() == item.get_name():
                return item
        raise AttributeError(f"Attribute '{name:}' is not defined")

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
        if name == "data_name":
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
            # if isinstance(item, self.CLASSES):
            if isinstance(item, (ItemN, LoopN)):
                self.items.append(item)
            elif type(self) is DataN:
                if issubclass(type(item), (ItemN, LoopN)):
                    self.CLASSES = self.CLASSES + (type(item), )
                    self.CLASSES_OPTIONAL = self.CLASSES_OPTIONAL + (type(item), )
                    self.items.append(item)

    @classmethod
    def make_container(cls, cls_mandatory, cls_optional, prefix):
        """Create DataN object as a container for items."""
        if cls is not DataN:
            warn("The method 'make_container' is used only for DataN class.")
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

    @classmethod
    def get_mandatory_attributes(cls, separator: str = "_"):
        """Get a list of mandatory attributes from mandatory classes."""
        l_res = []
        for cls_obj in cls.CLASSES_MANDATORY:
            if issubclass(cls_obj, ItemN):
                cls_item = cls_obj
            else: #LoopN
                cls_item = cls_obj.ITEM_CLASS
            l_res.extend([f"{cls_item.PREFIX:}{separator:}{name_cif:}" 
                              for name_cif in cls_item.ATTR_MANDATORY_CIF])
        return l_res

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
        data_name = self.data_name
        if data_name is not None:
            name = f"{name:}_{data_name:}"
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
        data_name = self.data_name
        if isinstance(data_name, str):
            data_name = data_name.lower()
        l_var = []
        for item in self.items:
            l_var.extend(item.get_variable_names())
        l_var_out = [((prefix, data_name), ) + var for var in l_var]
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
        data_name = self.data_name
        if isinstance(data_name, str):
            data_name = data_name.lower()

        prefix_d, prefix_n = name[0], name[1]
        if isinstance(prefix_d[1], str):
            prefix_d_2 = prefix_d[1].lower()
        else:
            prefix_d_2 = prefix_d[1]

        if (prefix_d[0], prefix_d_2) != (prefix, data_name):
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
        data_name = self.data_name

        prefix_d, prefix_n = name[0], name[1]

        if prefix_d != (prefix, data_name):
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
                if isinstance(item, ItemN):
                    warn(f"{item.PREFIX:} is not fully described.",
                         UserWarning)
                    break
                elif isinstance(item, LoopN):
                    warn(f"{item.ITEM_CLASS.PREFIX:} is not fully described.",
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
        """Print information about object in string in STAR format.

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
        if self.data_name is None:
            ls_out.append("data_\n")
        else:
            ls_out.append(f"data_{self.data_name:}\n")
        l_item = self.items
        l_s_itemn = [item.to_cif(separator=separator)+"\n"
                     for item in l_item if isinstance(item, ItemN)]
        l_s_loopn = [item.to_cif(separator=separator)+"\n"
                     for item in l_item if isinstance(item, LoopN)]
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
        return "\n".join(ls_out)

    @classmethod
    def from_cif(cls, string: str):
        """Generate object from string of CIF format."""
        cif_data = Data()
        flag = cif_data.take_from_string(string)

        cif_items = cif_data.items
        cif_loops = cif_data.loops

        items = []
        flag = True
        n_mandatory = len(cls.CLASSES_MANDATORY)

        # FIXME: if ItemN is not defined in data class it will be loosed. It should be fixed.        
        flag = len(cls.CLASSES) >= n_mandatory
        for cls_ in cls.CLASSES:
            if issubclass(cls_, ItemN):
                prefix_cls = cls_.PREFIX
                if cif_items.is_prefix(prefix_cls):
                    cif_items_prefix = cif_items[prefix_cls]
                    cif_string = str(cif_items_prefix)
                    obj_prefix = cls_.from_cif(cif_string)
                    if obj_prefix is not None:
                        items.append(obj_prefix)
                        flag = True
        
        for cif_loop in cif_loops:
            flag_loop = True
            for i_cls, cls_ in enumerate(cls.CLASSES):
                if issubclass(cls_, LoopN):
                    prefix_cls = cls_.ITEM_CLASS.PREFIX
                    if cif_loop.is_prefix("_"+prefix_cls):
                        cif_string = str(cif_loop)
                        obj_prefix = cls_.from_cif(cif_string)
                        if obj_prefix is not None:
                            items.append(obj_prefix)
                            flag = True
                            flag_loop = False
            if flag_loop:
                loopn = LoopN.from_cif(str(cif_loop))
                if not(loopn is None):
                    items.append(loopn)
                    

        if not(flag):
            return None

        data_name = cif_data.name
        obj = cls(data_name=data_name, items=items)
        obj.form_object()
        return obj

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

    def copy(self, data_name: str = ""):
        """Deep copy of object with new data name."""
        s_cif = self.to_cif()
        obj_new = type(self).from_cif(s_cif)
        obj_new.data_name = data_name
        return obj_new

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
                if len(l_name) == 1:
                    attr_refs = []
                    if isinstance(item, ItemN):
                        attr_refs = item.ATTR_REF
                    elif isinstance(item, LoopN):
                        item_class = item.ITEM_CLASS
                        if item_class is ItemN:
                            if len(self.items) != 0:
                                attr_refs = item.items[0].ATTR_REF
                        else:
                            attr_refs = item_class.ATTR_REF
                    for attr_ref in attr_refs:
                        item.set_variable(attr_ref, index=index)
                        
                else:
                    item.set_variable(".".join(l_name[1:]), index=index)

    def get_dictionary(self):
        self.form_object()
        res = {}
        res["name"] = getattr(self, "data_name")
        res["type_name"] = self.get_name()
        for item in self.items:
            if isinstance(item, (LoopN, ItemN)):
                dict = item.get_dictionary()
                res.update(dict)
        return res
    
    def take_parameters_from_dictionary(self, ddict_diffrn, l_parameter_name: list=None, l_sigma: list=None):
        """
        """
        for item in self.items:
            if isinstance(item, (LoopN, ItemN)):
                item.take_parameters_from_dictionary(ddict_diffrn, l_parameter_name=l_parameter_name, l_sigma=l_sigma)
        return None