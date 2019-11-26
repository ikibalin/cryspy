import pytest
from cryspy.f_crystal.cl_space_group_new import SpaceGroup

def test_init():
    try:
        _object = SpaceGroup()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif():
    try:
        _object = SpaceGroup()
        _str = _object.to_cif
        flag = True
    except:
        flag = False
    assert flag

    
