import pytest
import math
import numpy
from cryspy.common.cl_fitable import Fitable

from cryspy.cif_like.cl_setup import Setup

STR_FROM_CIF_1 = """
 _setup_wavelength   0.84
 _setup_field        1.00
 _setup_offset_ttheta -0.385
    """

def test_init():
    try:
        _object = Setup()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif():
    try:
        _object = Setup()
        _str = _object.to_cif
        flag = True
    except:
        flag = False
    assert flag


def test_from_cif():
    rel_tol, abs_tol =0.001, 0.001
    _obj = Setup.from_cif(STR_FROM_CIF_1)

    assert math.isclose(float(_obj.wavelength), 0.84, rel_tol=rel_tol, abs_tol=abs_tol)
    assert math.isclose(float(_obj.field), 1.00, rel_tol=rel_tol, abs_tol=abs_tol)
    assert math.isclose(float(_obj.offset_ttheta), -0.385, rel_tol=rel_tol, abs_tol=abs_tol)
    assert _obj.is_defined
    
