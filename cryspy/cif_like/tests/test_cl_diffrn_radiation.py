import pytest
import math
import numpy
from cryspy.common.cl_fitable import Fitable

from cryspy.cif_like.cl_diffrn_radiation import DiffrnRadiation

STR_FROM_CIF_1 = """
 _diffrn_radiation_efficiency    1.0
 _diffrn_radiation_polarization -0.87
    """

def test_init():
    try:
        _object = DiffrnRadiation()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif():
    try:
        _object = DiffrnRadiation()
        _str = _object.to_cif
        flag = True
    except:
        flag = False
    assert flag


def test_from_cif():
    rel_tol, abs_tol =0.001, 0.001
    _obj = DiffrnRadiation.from_cif(STR_FROM_CIF_1)

    assert math.isclose(float(_obj.polarization), -0.87, rel_tol=rel_tol, abs_tol=abs_tol)
    assert math.isclose(float(_obj.efficiency), 1.000, rel_tol=rel_tol, abs_tol=abs_tol)

    assert _obj.is_defined
    
