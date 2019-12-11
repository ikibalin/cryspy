import pytest
import numpy
from cryspy.common.cl_fitable import Fitable

from cryspy.pd1dcif_like.cl_pd_background import PdBackground, PdBackgroundL

STR_FROM_CIF_1 = """
 loop_
 _pd_background_ttheta
 _pd_background_intensity
  4.5  256.0
  40.0  158.0
  80.0  65.0
"""

def test_init():
    try:
        _object = PdBackgroundL()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif():
    try:
        _object = PdBackgroundL()
        _str = _object.to_cif
        flag = True
    except:
        flag = False
    assert flag

def test_init_el():
    try:
        _object = PdBackground()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif_el():
    try:
        _object = PdBackground()
        _str = _object.to_cif
        flag = True
    except:
        flag = False
    assert flag


def test_from_cif():
    l_obj = PdBackgroundL.from_cif(STR_FROM_CIF_1)
    assert len(l_obj) == 1
    _obj = l_obj[0]
    assert _obj.ttheta == [4.5, 40., 80.]
    assert float(_obj[40.].intensity) == 158.0
    val = _obj.to_cif(separator=".")

