# calculate d, sthovl based on hkl and cell for objects that supprot it
from pyparsing.common import Union

import numpy
from cryspy.B_parent_classes.cl_2_loop import LoopN
from cryspy.B_parent_classes.cl_1_item import ItemN

from cryspy.A_functions_base.unit_cell import calc_inv_d_by_unit_cell_parameters, calc_sthovl_by_unit_cell_parameters


def calc_d_sthovl_for_hkl(cell:ItemN, obj_hkl:Union[LoopN, ItemN]):
    """The object cell should have method 'get_unit_cell_parameters' and the object obj_hkl should have attributes 'index_h', 'index_k', 'index_l'
    """
    try:
        index_hkl = numpy.array([obj_hkl.index_h,obj_hkl.index_k,obj_hkl.index_l], dtype=float)
    except AttributeError:
        # Handle the case where the object doesn't have these attributes
        return
    try:
        unit_cell_parameters=cell.get_unit_cell_parameters()
    except AttributeError:
        # Handle the case where the object doesn't have the required method
        return
    inv_d = calc_inv_d_by_unit_cell_parameters(index_hkl=index_hkl, unit_cell_parameters=unit_cell_parameters)[0]
    sthovl = calc_sthovl_by_unit_cell_parameters(index_hkl=index_hkl, unit_cell_parameters=unit_cell_parameters)[0]
    d =1 / inv_d
    flag_loop = isinstance(obj_hkl, LoopN)
    try:
        if flag_loop:
            obj_hkl.numpy_d_spacing = d
        else:
            obj_hkl.d_spacing = d
    except AttributeError:
        # Handle the case where the object doesn't have these attributes
        pass
    try:
        if flag_loop:
            obj_hkl.numpy_sintlambda = sthovl
        else:
            obj_hkl.sintlambda = sthovl
    except AttributeError:
        # Handle the case where the object doesn't have these attributes
        pass
    if flag_loop:
        obj_hkl.numpy_to_items()
    return
