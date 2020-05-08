import pytest
import math
import numpy
from cryspy.common.cl_fitable import Fitable

from cryspy.cif_like.cl_chi2 import Chi2

STR_FROM_CIF_1 = """
 _chi2_sum  True
 _chi2_diff False
 _chi2_up   False
 _chi2_down False
    """

def test_init():
    try:
        _object = Chi2()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif():
    try:
        _object = Chi2()
        _str = _object.to_cif
        flag = True
    except:
        flag = False
    assert flag


def test_from_cif():
    rel_tol, abs_tol =0.001, 0.001
    _obj = Chi2.from_cif(STR_FROM_CIF_1)

    assert _obj.sum
    assert not(_obj.up)
    assert not(_obj.down)
    assert not(_obj.diff)
    
