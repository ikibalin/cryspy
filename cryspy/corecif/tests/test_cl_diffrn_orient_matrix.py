import pytest
import math
import numpy
from cryspy.common.cl_fitable import Fitable

from cryspy.corecif.cl_diffrn_orient_matrix import DiffrnOrientMatrix

STR_FROM_CIF_1 = """
_diffrn_orient_matrix_UB_11           -0.04170
_diffrn_orient_matrix_UB_12           -0.01429
_diffrn_orient_matrix_UB_13           -0.02226
_diffrn_orient_matrix_UB_21           -0.00380
_diffrn_orient_matrix_UB_22           -0.05578
_diffrn_orient_matrix_UB_23           -0.05048
_diffrn_orient_matrix_UB_31            0.00587
_diffrn_orient_matrix_UB_32           -0.13766
_diffrn_orient_matrix_UB_33            0.02277
    """

def test_init_el():
    try:
        _object = DiffrnOrientMatrix()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif_el():
    try:
        _object = DiffrnOrientMatrix()
        _str = _object.to_cif
        flag = True
    except:
        flag = False
    assert flag


def test_from_cif():
    rel_tol, abs_tol = 0.0001, 0.0001
    _obj = DiffrnOrientMatrix.from_cif(STR_FROM_CIF_1)
    assert math.isclose(_obj.ub_11, -0.04170, rel_tol=rel_tol, abs_tol=abs_tol)
    val = _obj.to_cif(separator=".")

