import pytest
import math
import numpy
from cryspy.common.cl_fitable import Fitable

from cryspy.magneticcif.cl_atom_site_moment import AtomSiteMoment, AtomSiteMomentL

STR_FROM_CIF_1 = """
loop_                                     
_atom_site_moment_label
_atom_site_moment_crystalaxis_x
_atom_site_moment_crystalaxis_y
_atom_site_moment_crystalaxis_z
Fe3A 4.8  0.0  0.0
Fe3B 0.0 -4.5  0.0
    """

def test_init():
    try:
        _object = AtomSiteMomentL()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif():
    try:
        _object = AtomSiteMomentL()
        _str = _object.to_cif
        flag = True
    except:
        flag = False
    assert flag

def test_init_el():
    try:
        _object = AtomSiteMoment()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif_el():
    try:
        _object = AtomSiteMoment()
        _str = _object.to_cif
        flag = True
    except:
        flag = False
    assert flag


def test_from_cif():
    rel_tol, abs_tol =0.001, 0.0001
    l_obj = AtomSiteMomentL.from_cif(STR_FROM_CIF_1)
    assert len(l_obj) == 1
    _obj = l_obj[0]
    assert _obj.label == ["Fe3A", "Fe3B"]
    assert math.isclose(_obj["Fe3A"].crystalaxis_x, 4.8, rel_tol=rel_tol, abs_tol=abs_tol)
    val = _obj.to_cif(separator=".")

