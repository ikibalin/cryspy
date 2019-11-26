import pytest
from cryspy.f_crystal.cl_space_group_2 import SpaceGroup2

def test_init():
    try:
        _object = SpaceGroup2()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif():
    try:
        _object = SpaceGroup2()
        _str = _object.to_cif
        flag = True
    except:
        flag = False
    assert flag

    
