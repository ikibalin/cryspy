import pytest
import os
import math
import numpy
from cryspy.common.cl_fitable import Fitable

from cryspy.scripts.cl_rhochi import RhoChi
f_diffrn_1 = os.path.join(os.path.dirname(__file__), "input_diffrn_1.rcif")
f_pd_1 = os.path.join(os.path.dirname(__file__), "input_pd_1.rcif")
f_pd2d_1 = os.path.join(os.path.dirname(__file__), "input_pd2d_1.rcif")

with open(f_diffrn_1, "r") as fid:
    STR_FROM_CIF_1 = fid.read()

with open(f_pd_1, "r") as fid:
    STR_FROM_CIF_2 = fid.read()

with open(f_pd2d_1, "r") as fid:
    STR_FROM_CIF_3 = fid.read()

def test_init():
    try:
        _object = RhoChi()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif():
    try:
        _object = RhoChi()
        _str = _object.to_cif
        flag = True
    except:
        flag = False
    assert flag


def test_from_cif_mono():
    rel_tol, abs_tol =0.001, 0.001
    _obj = RhoChi.from_cif(STR_FROM_CIF_1)
    assert _obj.crystals[0].space_group.it_number == 227
    assert len(_obj.get_variables()) == 2
    assert _obj.apply_constraint()
    chi_sq, n_res = _obj.calc_chi_sq()
    assert n_res == 10
    assert math.isclose(chi_sq, 75.2,rel_tol=rel_tol, abs_tol=abs_tol)
    

def test_from_cif_powder_1d():
    rel_tol, abs_tol =0.001, 0.001
    _obj = RhoChi.from_cif(STR_FROM_CIF_2)
    assert _obj.crystals[0].space_group.it_number == 227
    assert _obj.apply_constraint()
    chi_sq, n_res = _obj.calc_chi_sq()
    assert n_res == 381
    assert math.isclose(chi_sq, 1141.4532,rel_tol=rel_tol, abs_tol=abs_tol)
    


def test_from_cif_powder_2d():
    rel_tol, abs_tol =0.001, 0.001
    _obj = RhoChi.from_cif(STR_FROM_CIF_3)
    assert _obj.crystals[0].space_group.it_number == 227
    assert _obj.apply_constraint()
    chi_sq, n_res = _obj.calc_chi_sq()
    assert n_res == 25553
    assert math.isclose(chi_sq, 55600.84, rel_tol=rel_tol, abs_tol=abs_tol)
    