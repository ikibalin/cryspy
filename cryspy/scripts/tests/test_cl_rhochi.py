import pytest
import math
import numpy
from cryspy.common.cl_fitable import Fitable

from cryspy.scripts.cl_rhochi import RhoChi

STR_FROM_CIF_1 = """
 global_

 data_Fe3O4             
 _cell_angle_alpha 90.0                    
 _cell_angle_beta 90.0
 _cell_angle_gamma 90.0
 _cell_length_a 8.56212()
 _cell_length_b 8.56212
 _cell_length_c 8.56212
 _space_group_it_coordinate_system_code 2  
 _space_group_IT_number    227
 
 loop_                                     
 _atom_site_adp_type
 _atom_site_B_iso_or_equiv
 _atom_site_fract_x
 _atom_site_fract_y
 _atom_site_fract_z
 _atom_site_label
 _atom_site_occupancy
 _atom_site_type_symbol
  Uani 0.0 0.125 0.125 0.125 Fe3A 1.0 Fe3+
  Uani 0.0 0.5 0.5 0.5 Fe3B 1.0 Fe3+
  Uiso 0.0 0.25521 0.25521 0.25521 O1 1.0 O2-
 
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
 _atom_site_scat_label
 _atom_site_scat_lande
 Fe3A 2.0 
 Fe3B 2.0 
 
 loop_     
 _atom_site_susceptibility_label
 _atom_site_susceptibility_chi_type
 _atom_site_susceptibility_chi_11
 _atom_site_susceptibility_chi_12
 _atom_site_susceptibility_chi_13
 _atom_site_susceptibility_chi_22
 _atom_site_susceptibility_chi_23
 _atom_site_susceptibility_chi_33
  Fe3A Cani -3.468(74) 0.0 0.0 -3.468 0.0 -3.468
  Fe3B Cani 3.041      0.0 0.0  3.041 0.0  3.041

 data_mono
 _setup_wavelength     0.840
 _setup_field          1.000
 
 _diffrn_radiation_polarization 1.0
 _diffrn_radiation_efficiency   1.0

 _extinction_mosaicity 100.0
 _extinction_radius    50.0
 _extinction_model     gauss

 _diffrn_orient_matrix_UB_11 6.59783
 _diffrn_orient_matrix_UB_12 -6.99807
 _diffrn_orient_matrix_UB_13 3.3663
 _diffrn_orient_matrix_UB_21 2.18396
 _diffrn_orient_matrix_UB_22 -2.60871
 _diffrn_orient_matrix_UB_23 -9.5302
 _diffrn_orient_matrix_UB_31 7.4657
 _diffrn_orient_matrix_UB_32 6.94702
 _diffrn_orient_matrix_UB_33 -0.18685

 _phase_label  Fe3O4

 loop_
 _diffrn_refln_index_h
 _diffrn_refln_index_k
 _diffrn_refln_index_l
 _diffrn_refln_fr
 _diffrn_refln_fr_sigma
 0 0 8 0.64545 0.01329 
 2 0 6 1.75682 0.04540 
 0 2 6 1.67974 0.03711 
    """

def test_init():
    try:
        _object = RhoChi()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif():
    try:
        _object = RhoChi()
        _str = _object.to_cif
        flag = True
    except:
        flag = False
    assert flag


def test_from_cif():
    rel_tol, abs_tol =0.001, 0.001
    _obj = RhoChi.from_cif(STR_FROM_CIF_1)
    assert _obj.crystals[0].space_group.it_number == 227
    assert len(_obj.get_variables()) == 2
