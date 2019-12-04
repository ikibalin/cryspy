import pytest
import math
import numpy
from cryspy.f_common.cl_fitable import Fitable

from cryspy.f_crystal.cl_crystal_new import Crystal, Cell

STR_FROM_CIF_1 = """
data_Fe3O4                                
_cell_angle_alpha 90.0                    
_cell_angle_beta 90.0
_cell_angle_gamma 90.0
_cell_length_a 8.56212()
_cell_length_b 8.56212
_cell_length_c 8.56212
_space_group_it_coordinate_system_code 2  
_space_group_name_H-M_alt "F d -3 m"
_space_group_IT_number    232

loop_                                     
_atom_site_adp_type
_atom_site_B_iso_or_equiv
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_label
_atom_site_occupancy
_atom_site_type_symbol
 uani 0.0 0.125 0.125 0.125 Fe3A 1.0 Fe3+
 uani 0.0 0.5 0.5 0.5 Fe3B 1.0 Fe3+
 uiso 0.0 0.25521 0.25521 0.25521 O1 1.0 O2-

loop_                                     
_atom_type_scat_length_neutron
_atom_type_symbol
  0.945 Fe3+
 0.5803 O2-

loop_
_atom_site_aniso_U_11
_atom_site_aniso_U_12
_atom_site_aniso_U_13
_atom_site_aniso_U_22
_atom_site_aniso_U_23
_atom_site_aniso_U_33
_atom_site_aniso_label
 0.0 0.0 0.0 0.0 0.0 0.0 Fe3A
 0.0 0.0 0.0 0.0 0.0 0.0 Fe3B

loop_
_atom_site_magnetism_label
_atom_site_magnetism_lande
_atom_site_magnetism_kappa
Fe3A 2.0 1.0()
Fe3B 2.0() 1.0

loop_     
_atom_site_magnetism_aniso_label
_atom_site_magnetism_aniso_chi_type
_atom_site_magnetism_aniso_chi_11
_atom_site_magnetism_aniso_chi_12
_atom_site_magnetism_aniso_chi_13
_atom_site_magnetism_aniso_chi_22
_atom_site_magnetism_aniso_chi_23
_atom_site_magnetism_aniso_chi_33
 Fe3A cani -3.468(74) 0.0 0.0 -3.468 0.0 -3.468
 Fe3B cani 3.041      0.0 0.0  3.041 0.0  3.041
    """

def test_init():
    try:
        _object = Crystal()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif():
    try:
        _object = Crystal()
        _str = _object.to_cif
        flag = True
    except:
        flag = False
    assert flag

def test_attribute_cell():
    cell = Cell(9,9,9,90,90,90)
    cc = Crystal()
    cc.cell = cell

    cell_2 = Cell(1,1.5,2,74,90,90)
    cc.cell = cell_2
    l_obj_1 = cc[Cell]
    assert len(l_obj_1) == 1
    assert cc.cell.length_a.value == 1.


#def test_from_cif():
#    _obj = Cell.from_cif(STR_FROM_CIF_1)
#    assert math.isclose(float(_obj.length_a), 5.959, rel_tol =0.001, abs_tol=0.001)
#    cell=Cell(9.)
#    cell.apply_constraint("cP")
#    assert math.isclose(float(cell.length_b), 9.000, rel_tol =0.001, abs_tol=0.001)
    
