import pytest
from cryspy.f_crystal.cl_space_group_wyckoff import SpaceGroupWyckoff

def test_init():
    try:
        _object = SpaceGroupWyckoff()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif():
    try:
        _object = SpaceGroupWyckoff()
        _str = _object.to_cif
        flag = True
    except:
        flag = False
    assert flag

    
