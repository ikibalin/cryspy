import pytest
import math
import numpy
from cryspy.common.cl_fitable import Fitable

from cryspy.cif_like.cl_diffrn import Diffrn

STR_FROM_CIF_1 = """
 
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
        _object = Diffrn()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif():
    try:
        _object = Diffrn()
        _str = _object.to_cif
        flag = True
    except:
        flag = False
    assert flag

#def test_attribute_cell():
#    cell = Cell(9,9,9,90,90,90)
#    cc = Crystal()
#    cc.cell = cell
#
#    cell_2 = Cell(1,1.5,2,74,90,90)
#    cc.cell = cell_2
#    l_obj_1 = cc[Cell]
#    assert len(l_obj_1) == 1
#    assert cc.cell.length_a.value == 1.


def test_from_cif():
    rel_tol, abs_tol =0.001, 0.001
    _obj = Diffrn.from_cif(STR_FROM_CIF_1)
    assert _obj.is_defined
    assert _obj.extinction.mosaicity.value == 100.
    assert _obj.phase.label == "Fe3O4"

    print(_obj.to_cif(separator="."))
