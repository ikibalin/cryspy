import pytest
import numpy
from cryspy.common.cl_fitable import Fitable

from cryspy.corecif.cl_atom_site_aniso import AtomSiteAniso, AtomSiteAnisoL

STR_FROM_CIF_1 = """
loop_
_atom_site_aniso_label
_atom_site_aniso_U_11
_atom_site_aniso_U_22
_atom_site_aniso_U_33
_atom_site_aniso_U_12
_atom_site_aniso_U_13
_atom_site_aniso_U_23
 O1   .071(1) .076(1) .0342(9) .008(1)   .0051(9) -.0030(9) 
 C2   .060(2) .072(2) .047(1)  .002(2)   .013(1)  -.009(1)  
    """

def test_init():
    try:
        _object = AtomSiteAnisoL()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif():
    try:
        _object = AtomSiteAnisoL()
        _str = _object.to_cif
        flag = True
    except:
        flag = False
    assert flag

def test_init_el():
    try:
        _object = AtomSiteAniso()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif_el():
    try:
        _object = AtomSiteAniso()
        _str = _object.to_cif
        flag = True
    except:
        flag = False
    assert flag


def test_from_cif():
    l_obj = AtomSiteAnisoL.from_cif(STR_FROM_CIF_1)
    assert len(l_obj) == 1
    _obj = l_obj[0]
    assert _obj.label == ["O1", "C2"]
    assert isinstance(_obj["O1"].u_11, Fitable)
    val = _obj.to_cif(separator=".")

