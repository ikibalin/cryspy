import pytest
import numpy
from cryspy.common.cl_fitable import Fitable

from cryspy.pd1dcif_like.cl_pd_peak import PdPeak, PdPeakL

STR_FROM_CIF_1 = """
 loop_
 _pd_peak_index_h
 _pd_peak_index_k
 _pd_peak_index_l
 _pd_peak_index_mult
 _pd_peak_2theta
 _pd_peak_intensity_up
 _pd_peak_intensity_down
 _pd_peak_width_2theta
  2  2  0  4 17.2 100.0  90.0  2.3
"""

def test_init():
    try:
        _object = PdPeakL()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif():
    try:
        _object = PdPeakL()
        _str = _object.to_cif()
        flag = True
    except:
        flag = False
    assert flag

def test_init_el():
    try:
        _object = PdPeak()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif_el():
    try:
        _object = PdPeak()
        _str = _object.to_cif
        flag = True
    except:
        flag = False
    assert flag


def test_from_cif():
    l_obj = PdPeakL.from_cif(STR_FROM_CIF_1)
    assert len(l_obj) == 1
    _obj = l_obj[0]
    assert _obj.index_h == [2]
    assert _obj.index_k == [2]
    assert _obj.index_l == [0]
    val = _obj.to_cif(separator=".")

