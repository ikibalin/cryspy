import pytest
import math
import numpy
from cryspy.common.cl_fitable import Fitable

from cryspy.cif_like.cl_extinction import Extinction

STR_FROM_CIF_1 = """
 _extinction_mosaicity 100.0
 _extinction_radius    50.0
 _extinction_model     gauss
    """

def test_init():
    try:
        _object = Extinction()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif():
    try:
        _object = Extinction()
        _str = _object.to_cif
        flag = True
    except:
        flag = False
    assert flag


def test_from_cif():
    rel_tol, abs_tol =0.001, 0.001
    _obj = Extinction.from_cif(STR_FROM_CIF_1)

    assert math.isclose(float(_obj.mosaicity), 100., rel_tol=rel_tol, abs_tol=abs_tol)
    assert math.isclose(float(_obj.radius), 50., rel_tol=rel_tol, abs_tol=abs_tol)
    assert _obj.model == "gauss"
    assert _obj.is_defined
    
