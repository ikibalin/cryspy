import pytest
import math
import numpy
from cryspy.common.cl_fitable import Fitable

from cryspy.cif_like.cl_range import Range

STR_FROM_CIF_1 = """
 _range_ttheta_min     4.000
 _range_ttheta_max    80.000
    """

def test_init():
    try:
        _object = Range()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif():
    try:
        _object = Range()
        _str = _object.to_cif
        flag = True
    except:
        flag = False
    assert flag


def test_from_cif():
    rel_tol, abs_tol =0.001, 0.001
    _obj = Range.from_cif(STR_FROM_CIF_1)

    assert math.isclose(float(_obj.ttheta_min),  4.0, rel_tol=rel_tol, abs_tol=abs_tol)
    assert math.isclose(float(_obj.ttheta_max), 80.0, rel_tol=rel_tol, abs_tol=abs_tol)
    assert _obj.is_defined
    
