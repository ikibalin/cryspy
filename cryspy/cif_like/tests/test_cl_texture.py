import pytest
import math
import numpy
from cryspy.common.cl_fitable import Fitable

from cryspy.cif_like.cl_texture import Texture

STR_FROM_CIF_1 = """
 _texture_g_1 0.1239
 _texture_g_2 0.94211
 _texture_h_ax -0.66119
 _texture_k_ax -0.0541
 _texture_l_ax 3.0613
    """

def test_init():
    try:
        _object = Texture()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif():
    try:
        _object = Texture()
        _str = _object.to_cif
        flag = True
    except:
        flag = False
    assert flag


def test_from_cif():
    rel_tol, abs_tol =0.001, 0.001
    _obj = Texture.from_cif(STR_FROM_CIF_1)

    assert math.isclose(float(_obj.g_1), 0.1239, rel_tol=rel_tol, abs_tol=abs_tol)
    assert math.isclose(float(_obj.g_2), 0.94211, rel_tol=rel_tol, abs_tol=abs_tol)
    assert math.isclose(float(_obj.h_ax), -0.66119, rel_tol=rel_tol, abs_tol=abs_tol)
    assert _obj.is_defined
    
