import pytest
import numpy
from cryspy.f_common.cl_fitable import Fitable

from cryspy.f_crystal.corecif.cl_atom_site import AtomSite, AtomSiteL

STR_FROM_CIF_1 = """
    loop_                                     
    _atom_site_label          
    _atom_site_type_symbol   
    _atom_site_fract_x       
    _atom_site_fract_y       
    _atom_site_fract_z       
    _atom_site_adp_type       
    _atom_site_B_iso_or_equiv
    _atom_site_occupancy     
     Fe3A   Fe  0.12500 0.12500 0.12500  Uani   0.0   1.0
     Fe3B   Fe  0.50000 0.50000 0.50000  Uani   0.0   1.0
     O1     O   0.25521 0.25521 0.25521  Uiso   0.0   1.0
    """

def test_init():
    try:
        _object = AtomSiteL()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif():
    try:
        _object = AtomSiteL()
        _str = _object.to_cif
        flag = True
    except:
        flag = False
    assert flag

def test_init_el():
    try:
        _object = AtomSite()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif_el():
    try:
        _object = AtomSite()
        _str = _object.to_cif
        flag = True
    except:
        flag = False
    assert flag


def test_from_cif():
    l_obj = AtomSiteL.from_cif(STR_FROM_CIF_1)
    assert len(l_obj) == 1
    _obj = l_obj[0]
    assert _obj.label == ["Fe3A", "Fe3B", "O1"]
    assert _obj.type_symbol == ["Fe", "Fe", "O"]
    assert isinstance(_obj["Fe3A"].fract_x, Fitable)
    val = _obj.to_cif(separator=".")

