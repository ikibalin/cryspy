import pytest
import numpy
from cryspy.common.cl_fitable import Fitable

from cryspy.pd1dcif_like.cl_pd_phase import PdPhase, PdPhaseL

STR_FROM_CIF_1 = """

 loop_
 _pd_phase_label
 _pd_phase_scale
 _pd_phase_igsize
  Fe3O4 0.02381 0.0
"""

def test_init():
    try:
        _object = PdPhaseL()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif():
    try:
        _object = PdPhaseL()
        _str = _object.to_cif
        flag = True
    except:
        flag = False
    assert flag

def test_init_el():
    try:
        _object = PdPhase()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif_el():
    try:
        _object = PdPhase()
        _str = _object.to_cif
        flag = True
    except:
        flag = False
    assert flag


def test_from_cif():
    l_obj = PdPhaseL.from_cif(STR_FROM_CIF_1)
    assert len(l_obj) == 1
    _obj = l_obj[0]
    assert _obj.label == ["Fe3O4"]
    assert isinstance(_obj["Fe3O4"].igsize, Fitable)
    val = _obj.to_cif(separator=".")

