import pytest
from cryspy.f_crystal.cl_space_group_symop import SpaceGroupSymopEl, SpaceGroupSymop

STR_FROM_CIF_1 = """
loop_
_space_group_symop.id
_space_group_symop.operation_xyz
_space_group_symop.operation_description
  1    x,y,z              'identity mapping'
  2    -x,-y,-z           'inversion'
  3    -x,1/2+y,1/2-z
              '2-fold screw rotation with axis in (0,y,1/4)'
  4    x,1/2-y,1/2+z
            'c glide reflection through the plane (x,1/4,y)'
    """

STR_FROM_CIF_2 = """
loop_
_space_group_symop.id
_space_group_symop.generator_xyz
  1    x,1/2-y,1/2+z
    """

def test_init():
    try:
        _object = SpaceGroupSymop()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif():
    try:
        _object = SpaceGroupSymop()
        _str = _object.to_cif
        flag = True
    except:
        flag = False
    assert flag

def test_init_el():
    try:
        _object = SpaceGroupSymopEl()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif_el():
    try:
        _object = SpaceGroupSymopEl()
        _str = _object.to_cif
        flag = True
    except:
        flag = False
    assert flag

    
def test_init_2():
    i_1 = SpaceGroupSymopEl(id=1, operation_xyz="x,y,z")
    i_2 = SpaceGroupSymopEl(id="2", operation_xyz="-x,-y,-z")
    _object = SpaceGroupSymop(item=[i_1, i_2])
    assert all([_1==_2 for _1, _2 in  zip(_object.id, ["1", "2"])])
    assert all([_1==_2 for _1, _2 in  zip(_object.operation_xyz, ["x,y,z", "-x,-y,-z"])])
    assert _object["1"].operation_xyz == "x,y,z"
    assert _object["2"].operation_xyz == "-x,-y,-z"
    assert _object["1"].id == "1"
    assert _object.prefix == "space_group_symop"

"""
def test_from_cif_1():
    try: 
        _object = SpaceGroupSymop.from_cif(STR_FROM_CIF_1)
        flag = True
    except:
        flag = False
    assert flag
    assert _object.item("1").operation_xyz == "x,y,z"
    assert _object.item("2").operation_xyz == "-x,-y,-z"
    assert _object.item("2").operation_description == "inversion"
    assert _object.operation_xyz == ["x,y,z", "-x,-y,-z", "-x,1/2+y,1/2-z", "x,1/2-y,1/2+z"]
    assert _object.generator_xyz == []
    assert _object.item("1").generator_xyz == None
    assert _object.sg_id == []
    assert _object.item("2").sg_id == None


def test_from_cif_2():
    try: 
        _object = SpaceGroupSymop.from_cif(STR_FROM_CIF_2)
        flag = True
    except:
        flag = False
    assert flag
    assert _object.generator_xyz == ["x,1/2-y,1/2+z"]
    assert _object("1").generator_xyz == "x,1/2-y,1/2+z"
    assert _object.sg_id == []
    assert _object("1").sg_id == None
"""