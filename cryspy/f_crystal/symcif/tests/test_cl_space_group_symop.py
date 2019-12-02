import pytest
import numpy
from cryspy.f_crystal.symcif.cl_space_group_symop import SpaceGroupSymop, SpaceGroupSymopL

STR_FROM_CIF_1 = """
loop_
_space_group_symop.operation_xyz
_space_group_symop.id
_space_group_symop.operation_description
      x,y,z          1    'identity mapping'
      -x,-y,-z       2    inversion
      -x,1/2+y,1/2-z 3 '2-fold screw rotation with axis in (0,y,1/4)'
      x,1/2-y,1/2+z 4 'c glide reflection through the plane (x,1/4,y)'
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


def test_get_symops_by_generator_xyz():
    _obj = SpaceGroupSymop(operation_xyz="x,y,z")
    symops = _obj.get_symops_by_generator_xyz("x, 1/4-y, 1/4+z")
    assert len(symops) == 4
    assert symops[0].operation_xyz == "x,y,z"
    assert symops[1].operation_xyz == "x,-y+1/4,z+1/4"
    assert symops[2].operation_xyz == "x,y,z+1/2"
    assert symops[3].operation_xyz == "x,-y+1/4,z+3/4"

def test_init_2():
    i_1 = SpaceGroupSymop(id=1, operation_xyz="x,y,z")
    i_2 = SpaceGroupSymop(id="2", operation_xyz="-x+1/2,-y,-z")
    assert float(i_2.r_11) == -1.
    assert float(i_2.r_12) ==  0.
    assert float(i_2.r_13) ==  0.
    assert float(i_2.r_21) ==  0.
    assert float(i_2.r_22) == -1.
    assert float(i_2.r_23) ==  0.
    assert float(i_2.r_31) ==  0.
    assert float(i_2.r_32) ==  0.
    assert float(i_2.r_33) == -1.
    assert float(i_2.b_1) == 0.5
    assert float(i_2.b_2) == 0.
    assert float(i_2.b_3) == 0.
    assert all(numpy.array(i_2.b, dtype=float) == numpy.array([0.5,0.,0.], dtype=float))
    assert numpy.array(i_2.r, dtype=float)[0,0] == -1.

    _object = SpaceGroupSymopL(item=[i_1, i_2])
    assert all([_1==_2 for _1, _2 in  zip(_object.id, ["1", "2"])])
    assert all([_1==_2 for _1, _2 in  zip(_object.operation_xyz, ["x,y,z", "-x+1/2,-y,-z"])])
    assert _object["1"].operation_xyz == "x,y,z"
    assert _object["2"].operation_xyz == "-x+1/2,-y,-z"
    assert _object["1"].id == "1"
    assert _object.prefix == "space_group_symop"

    assert all([float(_1)==_2 for _1, _2 in zip(_object.r_11, [1., -1.])])

def test_from_cif():
    l_obj = SpaceGroupSymopL.from_cif(STR_FROM_CIF_1)
    assert len(l_obj) == 1
    _obj = l_obj[0]
    assert _obj.id == ["1", "2", "3", "4"]
    assert _obj.operation_xyz == ["x,y,z","-x,-y,-z","-x,1/2+y,1/2-z","x,1/2-y,1/2+z"]
    assert _obj["3"].operation_xyz == "-x,1/2+y,1/2-z"
    _str_1 = STR_FROM_CIF_1.replace(" ","").replace("\n","").lower()
    _str_2 = _obj.to_cif(separator=".").replace(" ","").replace("\n","").lower()
    assert _str_1 == _str_2 


def test_create_by_generators_xyz():
    obj = SpaceGroupSymopL.create_by_generators_xyz(("-x,-y,-z", "x, 1/2+y, 1/2+z", "x, 1/4+y, 1/4+z"))
    assert len(obj.item) == 8
    assert obj["4"].operation_xyz == "x,y+3/4,z+3/4"
    assert obj["6"].operation_xyz == "-x,-y+1/4,-z+1/4"
