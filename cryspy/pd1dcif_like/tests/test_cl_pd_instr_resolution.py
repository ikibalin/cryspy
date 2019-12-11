import pytest
import math
import numpy
from cryspy.common.cl_fitable import Fitable

from cryspy.pd1dcif_like.cl_pd_instr_resolution import PdInstrResolution

STR_FROM_CIF_1 = """
 _pd_instr_resolution_u 16.9776
 _pd_instr_resolution_v -2.8357
 _pd_instr_resolution_w  0.5763
 _pd_instr_resolution_x  0.1000
 _pd_instr_resolution_y  0.0000
    """

def test_init():
    try:
        _object = PdInstrResolution()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif():
    try:
        _object = PdInstrResolution()
        _str = _object.to_cif
        flag = True
    except:
        flag = False
    assert flag


def test_from_cif():
    rel_tol, abs_tol =0.001, 0.001
    _obj = PdInstrResolution.from_cif(STR_FROM_CIF_1)

    assert math.isclose(float(_obj.u), 16.9776, rel_tol=rel_tol, abs_tol=abs_tol)
    assert math.isclose(float(_obj.v), -2.8357, rel_tol=rel_tol, abs_tol=abs_tol)
    assert math.isclose(float(_obj.w), 0.5763, rel_tol=rel_tol, abs_tol=abs_tol)
    assert math.isclose(float(_obj.x), 0.1, rel_tol=rel_tol, abs_tol=abs_tol)
    assert math.isclose(float(_obj.y), 0.0, rel_tol=rel_tol, abs_tol=abs_tol)
    assert _obj.is_defined
    
