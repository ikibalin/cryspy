import pytest
import numpy
from cryspy.common.cl_fitable import Fitable

from cryspy.pd2dcif_like.cl_pd2d_background import Pd2dBackground

STR_FROM_CIF_1 = """
 _pd2d_background_2theta_phi_intensity
;
      2    4.5     40.0     80.0
 -3.000 -350.0()   -350.0   -400.0
 41.000 -351.0()   -350.0   -400.0
;
"""

def test_init():
    try:
        _object = Pd2dBackground()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif():
    try:
        _object = Pd2dBackground()
        _str = _object.to_cif
        flag = True
    except:
        flag = False
    assert flag

def test_init_el():
    try:
        _object = Pd2dBackground()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif_el():
    try:
        _object = Pd2dBackground()
        _str = _object.to_cif
        flag = True
    except:
        flag = False
    assert flag


def test_from_cif():
    _obj = Pd2dBackground.from_cif(STR_FROM_CIF_1)
    print(_obj)
    assert _obj.phi[0] == -3.00
    assert _obj.phi[1] == 41.00
    assert _obj.ttheta[1] == 40.00
    assert _obj.intensity[0][1].refinement
    assert _obj.intensity[0][1].value == -351.00

