import pytest
import math
import numpy
from cryspy.common.cl_fitable import Fitable

from cryspy.cif_like.cl_diffrn_refln import DiffrnRefln, DiffrnReflnL

STR_FROM_CIF_1 = """
 loop_
 _diffrn_refln_index_h
 _diffrn_refln_index_k
 _diffrn_refln_index_l
 _diffrn_refln_fr
 _diffrn_refln_fr_sigma
     0    0    8   0.64545   0.01329 
     2    0    6   1.75682   0.0454  
    """

def test_init():
    try:
        _object = DiffrnReflnL()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif():
    try:
        _object = DiffrnReflnL()
        _str = _object.to_cif
        flag = True
    except:
        flag = False
    assert flag

def test_init_el():
    try:
        _object = DiffrnRefln()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif_el():
    try:
        _object = DiffrnRefln()
        _str = _object.to_cif
        flag = True
    except:
        flag = False
    assert flag


def test_from_cif():
    rel_tol, abs_tol = 0.0001, 0.0001
    l_obj = DiffrnReflnL.from_cif(STR_FROM_CIF_1)
    assert len(l_obj) == 1
    _obj = l_obj[0]
    assert _obj.index_h == [0, 2]
    #assert math.isclose(_obj[0,0,2].a_calc, 3.25, rel_tol=rel_tol, abs_tol=abs_tol)
    val = _obj.to_cif(separator=".")

