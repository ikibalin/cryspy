import pytest
import numpy
from cryspy.common.cl_fitable import Fitable

from cryspy.pd2dcif_like.cl_pd2d_peak import Pd2dPeak, Pd2dPeakL

STR_FROM_CIF_1 = """
 loop_
 _pd2d_peak_index_h
 _pd2d_peak_index_k
 _pd2d_peak_index_l
 _pd2d_peak_index_mult
 _pd2d_peak_ttheta
 _pd2d_peak_f_nucl_sq
 _pd2d_peak_f_m_p_sin_sq
 _pd2d_peak_f_m_p_cos_sq
 _pd2d_peak_cross_sin
 _pd2d_peak_width_ttheta
  2  2  0  4 17.2 100.0 101.2   90.0  87.4 2.3
  2  1  3  2 19.1  78.6 101.0   92.0  82.7 5.7
"""

def test_init():
    try:
        _object = Pd2dPeakL()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif():
    try:
        _object = Pd2dPeakL()
        _str = _object.to_cif()
        flag = True
    except:
        flag = False
    assert flag

def test_init_el():
    try:
        _object = Pd2dPeak()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif_el():
    try:
        _object = Pd2dPeak()
        _str = _object.to_cif
        flag = True
    except:
        flag = False
    assert flag


def test_from_cif():
    l_obj = Pd2dPeakL.from_cif(STR_FROM_CIF_1)
    assert len(l_obj) == 1
    _obj = l_obj[0]
    assert _obj.index_h == [2, 2]
    assert _obj.index_k == [2, 1]
    assert _obj.index_l == [0, 3]
    val = _obj.to_cif(separator=".")

