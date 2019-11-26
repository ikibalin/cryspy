import pytest
from cryspy.f_crystal.cl_atom_site import AtomSite

def test_atom_site_init():
    try:
        _object = AtomSite()
        flag = True
    except:
        flag = False
    assert flag

def test_atom_site_empty_to_cif():
    try:
        _object = AtomSite()
        _str = _object.to_cif
        flag = True
    except:
        flag = False
    assert flag
