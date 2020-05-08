import pytest
import math
import numpy
from cryspy.common.cl_fitable import Fitable

from cryspy.corecif.cl_refln import Refln, ReflnL

STR_FROM_CIF_1 = """
 loop_
  _refln_index_h
  _refln_index_k
  _refln_index_l
  _refln_d_spacing
  _refln_A_calc
  _refln_B_calc
  0 0 2 2.315 3.25  1.232
  2 2 0 4.213 5.00 -4.05
    """

def test_init():
    try:
        _object = ReflnL()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif():
    try:
        _object = ReflnL()
        _str = _object.to_cif
        flag = True
    except:
        flag = False
    assert flag

def test_init_el():
    try:
        _object = Refln()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif_el():
    try:
        _object = Refln()
        _str = _object.to_cif
        flag = True
    except:
        flag = False
    assert flag


def test_from_cif():
    rel_tol, abs_tol = 0.0001, 0.0001
    l_obj = ReflnL.from_cif(STR_FROM_CIF_1)
    assert len(l_obj) == 1
    _obj = l_obj[0]
    assert _obj.index_h == [0, 2]
    #assert math.isclose(_obj[0,0,2].a_calc, 3.25, rel_tol=rel_tol, abs_tol=abs_tol)
    val = _obj.to_cif(separator=".")

