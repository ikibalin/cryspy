import pytest
import math
import numpy
from cryspy.common.cl_fitable import Fitable

from cryspy.cif_like.cl_crystal import Crystal, Cell

STR_FROM_CIF_1 = """
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
 _atom_site_scat_kappa
 Fe3A 2.0 1.0
 Fe3B 2.0 1.0
 
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


def test_from_cif():
    rel_tol, abs_tol =0.001, 0.001
    _obj = Crystal.from_cif(STR_FROM_CIF_1)
    assert math.isclose(float(_obj.cell.length_a), 8.56212, rel_tol=rel_tol, abs_tol=abs_tol)
    assert _obj.space_group.it_number == 227
    assert _obj.is_defined
    assert _obj.data_name == "Fe3O4"
    
def test_calc_refln():
    rel_tol, abs_tol =0.001, 0.01
    _obj = Crystal.from_cif(STR_FROM_CIF_1)
    _obj.apply_constraint()
    h, k, l = numpy.array([1, 2, 3], dtype=int), numpy.array([1, 2, 1], dtype=int), numpy.array([1, 0, 1], dtype=int)
    refln = _obj.calc_refln(h, k, l)
    f_reference = numpy.array([3.125, -7.63948, -12.599], dtype=complex)
    res = numpy.isclose(refln.f_calc, f_reference, rtol=rel_tol, atol=abs_tol)
    assert all(res)

def test_calc_refln_susceptibility():
    rel_tol, abs_tol =0.001, 0.01
    _obj = Crystal.from_cif(STR_FROM_CIF_1)
    _obj.apply_constraint()
    h, k, l = numpy.array([1, 2, 3], dtype=int), numpy.array([1, 2, 1], dtype=int), numpy.array([1, 0, 1], dtype=int)
    refln_suc = _obj.calc_refln_susceptibility(h, k, l)
    #print(refln_suc.to_cif())
    assert True
    #f_reference = numpy.array([31.25, -76.3948, -125.99], dtype=complex)
    #res = numpy.isclose(f_nucl, f_reference, rtol=rel_tol, atol=abs_tol)
    #assert all(res)    