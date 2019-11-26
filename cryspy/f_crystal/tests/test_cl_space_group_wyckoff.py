import pytest
from cryspy.f_crystal.cl_space_group_wyckoff import SpaceGroupWyckoffEl, SpaceGroupWyckoff

STR_FROM_CIF_1 = """
loop_
_space_group_Wyckoff_id
_space_group_Wyckoff_coord_xyz
_space_group_Wyckoff_letter
_space_group_Wyckoff_multiplicity
_space_group_Wyckoff_site_symmetry
   1    x,y,z         h 192   1     
   2    1/4,y,-y      g  96   ..2   
   3    x,1/8,1/8     f  96   2..   
   4    1/4,1/4,1/4   b  32   .32   
    """


def test_init():
    _object = SpaceGroupWyckoff()
    try:
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

def test_init_el():
    try:
        _object = SpaceGroupWyckoffEl()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif_el():
    try:
        _object = SpaceGroupWyckoffEl()
        _str = _object.to_cif
        flag = True
    except:
        flag = False
    assert flag

    
def test_init_2():
    i_1 = SpaceGroupWyckoffEl(id=1, coord_xyz="x,y,z")
    i_2 = SpaceGroupWyckoffEl(id="2", coord_xyz="-x,-y,-z")
    _object = SpaceGroupWyckoff(item=[i_1, i_2])
    assert all([_1==_2 for _1, _2 in  zip(_object.id, ["1", "2"])])
    assert all([_1==_2 for _1, _2 in  zip(_object.coord_xyz, ["x,y,z", "-x,-y,-z"])])
    assert _object["1"].coord_xyz == "x,y,z"
    assert _object["2"].coord_xyz == "-x,-y,-z"
    assert _object["1"].id == "1"
    assert _object.prefix == "space_group_Wyckoff"

def test_from_cif():
    l_obj = SpaceGroupWyckoff.from_cif(STR_FROM_CIF_1)
    assert len(l_obj) == 1
    _obj = l_obj[0]
    assert _obj.id == ["1", "2", "3", "4"]
    assert _obj.multiplicity == [192, 96, 96, 32]
    assert _obj["3"].multiplicity == 96
    _str_1 = STR_FROM_CIF_1.replace(" ","").replace("\n","").lower()
    _str_2 = _obj.to_cif().replace(" ","").replace("\n","").lower()
    assert _str_1 == _str_2 