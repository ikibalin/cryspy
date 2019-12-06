import pytest
import math
import numpy
from cryspy.common.cl_fitable import Fitable

from cryspy.magneticcif.cl_refln_susceptibility import ReflnSusceptibility, ReflnSusceptibilityL

STR_FROM_CIF_1 = """
loop_
_refln_susceptibility.index_h    
_refln_susceptibility.index_k    
_refln_susceptibility.index_l    
_refln_susceptibility.d_spacing       
_refln_susceptibility.sintlambda      
_refln_susceptibility.chi_11_calc     
_refln_susceptibility.chi_12_calc     
_refln_susceptibility.chi_13_calc     
_refln_susceptibility.chi_21_calc     
_refln_susceptibility.chi_22_calc     
_refln_susceptibility.chi_23_calc     
_refln_susceptibility.chi_31_calc     
_refln_susceptibility.chi_32_calc     
_refln_susceptibility.chi_33_calc     
_refln_susceptibility.moment_11_calc  
_refln_susceptibility.moment_12_calc  
_refln_susceptibility.moment_13_calc  
_refln_susceptibility.moment_21_calc  
_refln_susceptibility.moment_22_calc  
_refln_susceptibility.moment_23_calc  
_refln_susceptibility.moment_31_calc  
_refln_susceptibility.moment_32_calc  
_refln_susceptibility.moment_33_calc  
2 0 0 4.52 0.123 0+0j 0+0j 0+0j 0+0j 0+0j 0+0j 0+0j 0+0j 0+0j  0+0j 0+0j 0+0j 0+0j 0+0j 0+0j 0+0j 0+0j 0+0j
0 2 0 4.52 0.123 0+0j 0+0j 0+0j 0+0j 0+0j 0+0j 0+0j 0+0j 0+0j  0+0j 0+0j 0+0j 0+0j 0+0j 0+0j 0+0j 0+0j 0+0j
    """

def test_init():
    try:
        _object = ReflnSusceptibilityL()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif():
    try:
        _object = ReflnSusceptibilityL()
        _str = _object.to_cif
        flag = True
    except:
        flag = False
    assert flag

def test_init_el():
    try:
        _object = ReflnSusceptibility()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif_el():
    try:
        _object = ReflnSusceptibility()
        _str = _object.to_cif
        flag = True
    except:
        flag = False
    assert flag


def test_from_cif():
    rel_tol, abs_tol =0.001, 0.0001
    l_obj = ReflnSusceptibilityL.from_cif(STR_FROM_CIF_1)
    assert len(l_obj) == 1
    _obj = l_obj[0]
    assert _obj.index_h == [2, 0]
    assert _obj.chi_11_calc == [complex(0), complex(0)]
    val = _obj.to_cif(separator=".")

