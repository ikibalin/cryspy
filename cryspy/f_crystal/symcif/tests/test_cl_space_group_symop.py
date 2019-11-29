import pytest
from cryspy.f_crystal.symcif.cl_space_group_symop import SpaceGroupSymop, SpaceGroupSymopL

STR_FROM_CIF_1 = """
loop_
_space_group_symop.id
_space_group_symop.operation_xyz
_space_group_symop.operation_description
  1    x,y,z              'identity mapping'
  2    -x,-y,-z           inversion
  3    -x,1/2+y,1/2-z '2-fold screw rotation with axis in (0,y,1/4)'
  4    x,1/2-y,1/2+z 'c glide reflection through the plane (x,1/4,y)'
    """

STR_FROM_CIF_2 = """
loop_
_space_group_symop.id
_space_group_symop.generator_xyz
  1    x,1/2-y,1/2+z
    """

def test_init():
    try:
        _object = SpaceGroupSymopL()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif():
    try:
        _object = SpaceGroupSymopL()
        _str = _object.to_cif
        flag = True
    except:
        flag = False
    assert flag

def test_init_el():
    try:
        _object = SpaceGroupSymop()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif_el():
    try:
        _object = SpaceGroupSymop()
        _str = _object.to_cif
        flag = True
    except:
        flag = False
    assert flag

    
def test_init_2():
    i_1 = SpaceGroupSymop(id=1, operation_xyz="x,y,z")
    i_2 = SpaceGroupSymop(id="2", operation_xyz="-x,-y,-z")
    _object = SpaceGroupSymopL(item=[i_1, i_2])
    assert all([_1==_2 for _1, _2 in  zip(_object.id, ["1", "2"])])
    assert all([_1==_2 for _1, _2 in  zip(_object.operation_xyz, ["x,y,z", "-x,-y,-z"])])
    assert _object["1"].operation_xyz == "x,y,z"
    assert _object["2"].operation_xyz == "-x,-y,-z"
    assert _object["1"].id == "1"
    assert _object.prefix == "space_group_symop"


def test_from_cif():
    l_obj = SpaceGroupSymopL.from_cif(STR_FROM_CIF_1)
    assert len(l_obj) == 1
    _obj = l_obj[0]
    assert _obj.id == ["1", "2", "3", "4"]
    assert _obj.operation_xyz == ["x,y,z","-x,-y,-z","-x,1/2+y,1/2-z","x,1/2-y,1/2+z"]
    assert _obj["3"].operation_xyz == "-x,1/2+y,1/2-z"
    _str_1 = STR_FROM_CIF_1.replace(" ","").replace("\n","").lower()
    _str_2 = _obj.to_cif(separator=".").replace(" ","").replace("\n","").lower()
    print(_obj.to_cif(separator="."))
    assert _str_1 == _str_2 


