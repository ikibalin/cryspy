import pytest
from cryspy.symcif.cl_space_group import SpaceGroup




STR_FROM_CIF_1 = """
    _space_group.id                    1
    _space_group.name_HM_ref            'C 2/c'
    _space_group.name_Hall           '-C 2yc'
    _space_group.name_Schoenflies      C2h.6
    _space_group.IT_number             15
    _space_group.Laue_class            2/m
    _space_group.Patterson_name_HM  'C 2/m'
    _space_group.centring_type         C
    _space_group.Bravais_type          mS
    _space_group.crystal_system        monoclinic
    """  # TODO: temporary test, it should be: _space_group.Patterson_name_H-M-M_ref            'C 2/m'
    #                                          _space_group.name_H-M_ref            'C 2/c'

def test_init():
    _obj = SpaceGroup(name_schoenflies= "C2h.6")
    assert _obj.it_number == 15
    try:
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif():
    _obj = SpaceGroup(it_number=15)
    assert _obj.name_hm_ref == 'C 2/c'
    assert _obj.crystal_system == 'monoclinic'
    assert _obj.name_schoenflies == 'C2h.6'
    assert _obj.it_coordinate_system_code == "b1"    
    _str = _obj.to_cif(separator=".")
    try:
        flag = True
    except:
        flag = False
    assert flag


def test_from_cif():
    _obj = SpaceGroup.from_cif(STR_FROM_CIF_1)
    assert _obj is not None
    assert _obj.id == "1"
    assert _obj.it_number == 15
    assert _obj.name_hm_ref == 'C 2/c'
    assert _obj.crystal_system == 'monoclinic'

    assert _obj.it_coordinate_system_code == "b1"    
    #_str_1 = STR_FROM_CIF_1.replace(" ","").replace("\n","").lower()
    #_str_2 = _obj.to_cif(separator=".").replace(" ","").replace("\n","").lower()
    #print(_obj.to_cif(separator="."))
    #assert _str_1 == _str_2 

def test_space_group_wyckoff():
    _obj = SpaceGroup(it_number=64)
    assert _obj.space_group_wyckoff["5"].coord_xyz == "1/4,1/4,0"

    obj_2 = _obj.space_group_wyckoff
    assert obj_2.get_id_for_fract(0.75,0.75,0.5) == "5"
    assert obj_2.get_letter_for_fract(0.75,0.75,0.5) == "c"
    x_s, y_s, z_s, mult = _obj.calc_xyz_mult(0,0,0)
    assert x_s.size == 4
    assert ((x_s[0] == 0.) & (x_s[1] == 0.)  & (x_s[2] == 0.5) & (x_s[3] == 0.5) )
    assert ((y_s[0] == 0.) & (y_s[1] == 0.5) & (y_s[2] == 0.5) & (y_s[3] == 0. ))
    assert ((z_s[0] == 0.) & (z_s[1] == 0.5) & (z_s[2] == 0. ) & (z_s[3] == 0.5))
    assert mult == 4
