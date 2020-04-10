import pytest
import math
import numpy
from cryspy.common.cl_fitable import Fitable

from cryspy.corecif.cl_cell import Cell

STR_FROM_CIF_1 = """
     _cell_length_a                     5.959(1)
     _cell_length_b                     14.956(1)
     _cell_length_c                     19.737(3)
     _cell_angle_alpha                  90
     _cell_angle_beta                   90
     _cell_angle_gamma                  90
    """

def test_init():
    try:
        _object = Cell()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif():
    try:
        _object = Cell()
        _str = _object.to_cif
        flag = True
    except:
        flag = False
    assert flag



def test_from_cif():
    _obj = Cell.from_cif(STR_FROM_CIF_1)
    assert math.isclose(float(_obj.length_a), 5.959, rel_tol =0.001, abs_tol=0.001)
    cell=Cell(9.)
    cell.apply_constraint("cP", "")
    assert math.isclose(float(cell.length_b), 9.000, rel_tol =0.001, abs_tol=0.001)
    
