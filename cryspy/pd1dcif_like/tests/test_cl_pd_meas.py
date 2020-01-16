import pytest
import numpy
from cryspy.common.cl_fitable import Fitable

from cryspy.pd1dcif_like.cl_pd_meas import PdMeas, PdMeasL

STR_FROM_CIF_1 = """
 loop_
 _pd_meas_2theta
 _pd_meas_intensity_up
 _pd_meas_intensity_up_sigma
 _pd_meas_intensity_down
 _pd_meas_intensity_down_sigma
  4.00   465.80000   128.97000   301.88000   129.30000
  4.20   323.78000   118.22000   206.06000   120.00000
"""

def test_init():
    try:
        _object = PdMeasL()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif():
    try:
        _object = PdMeasL()
        _str = _object.to_cif
        flag = True
    except:
        flag = False
    assert flag

def test_init_el():
    try:
        _object = PdMeas()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif_el():
    try:
        _object = PdMeas()
        _str = _object.to_cif
        flag = True
    except:
        flag = False
    assert flag


def test_from_cif():
    l_obj = PdMeasL.from_cif(STR_FROM_CIF_1)
    assert len(l_obj) == 1
    _obj = l_obj[0]
    assert _obj.ttheta == [4.00, 4.20]
    assert _obj[4.20].intensity_up == 323.78000
    val = _obj.to_cif(separator=".")

