import pytest
import math
import numpy
from cryspy.common.cl_fitable import Fitable

from cryspy.cif_like.cl_pd2d import Pd2d

STR_FROM_CIF_1 = """
 data_pnd
 _setup_wavelength     0.840
 _setup_field          1.000
 _setup_offset_2theta -0.385
 
 _chi2_sum  True
 _chi2_diff False
 _chi2_up   False
 _chi2_down False
 
 _range_2theta_min     4.000
 _range_2theta_max    80.000
 _range_phi_min       -2.000
 _range_phi_max       40.000
 
 _pd2d_background_2theta_phi_intensity
;
           2   4.50        40.00        80.00      
       -3.00 -350.0        -369(1)     -400.0      
       41.00 -350.0        -404(2)     -400.0      
;

 
 loop_
 _exclude_2theta_min
 _exclude_2theta_max
 _exclude_phi_min
 _exclude_phi_max
 0.0 1.0 0.0 1.0
 
 _pd2d_instr_reflex_asymmetry_p1 0.0
 _pd2d_instr_reflex_asymmetry_p2 0.0
 _pd2d_instr_reflex_asymmetry_p3 0.0
 _pd2d_instr_reflex_asymmetry_p4 0.0
 
 _diffrn_radiation_polarization 1.0
 _diffrn_radiation_efficiency   1.0
 
 _pd2d_instr_resolution_u 16.9776
 _pd2d_instr_resolution_v -2.8357
 _pd2d_instr_resolution_w  0.5763
 _pd2d_instr_resolution_x  0.0
 _pd2d_instr_resolution_y  0.0
 
 loop_
 _phase_label
 _phase_scale
 _phase_igsize
 Fe3O4 0.02381 0.0
 
 _pd2d_meas_2theta_phi_intensity_up
;
     2    4.5     40.0     80.0
-3.000 -350.0   -350.0   -400.0
41.000 -351.0   -350.0   -400.0
;

_pd2d_meas_2theta_phi_intensity_up_sigma
;
     2    4.5     40.0     80.0
-3.000 -352.0   -350.0   -400.0
41.000 -353.0   -350.0   -400.0
;

_pd2d_meas_2theta_phi_intensity_down
;
     2    4.5     40.0     80.0
-3.000 -354.0   -350.0   -400.0
41.000 -355.0   -350.0   -400.0
;

_pd2d_meas_2theta_phi_intensity_down_sigma
;
     2    4.5     40.0     80.0
-3.000 -356.0   -350.0   -400.0
41.000 -357.0   -350.0   -400.0
;
     """

def test_init():
    try:
        _object = Pd2d()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif():
    try:
        _object = Pd2d()
        _str = _object.to_cif
        flag = True
    except:
        flag = False
    assert flag



def test_from_cif():
    rel_tol, abs_tol =0.001, 0.001
    _obj = Pd2d.from_cif(STR_FROM_CIF_1)
    assert _obj.is_defined
    assert _obj.chi2.sum
    assert not(_obj.chi2.diff)
    print(_obj.to_cif(separator="."))
