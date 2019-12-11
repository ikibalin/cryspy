import pytest
import math
import numpy
from cryspy.common.cl_fitable import Fitable

from cryspy.cif_like.cl_pd import Pd

STR_FROM_CIF_1 = """
 data_pnd
 _setup_wavelength     0.840
 _setup_field          1.000
 _setup_offset_ttheta -0.385
 
 _chi2_sum  True
 _chi2_diff False
 _chi2_up   False
 _chi2_down False
 
 _range_ttheta_min     4.000
 _range_ttheta_max    80.000
 
 loop_
 _pd_background_ttheta
 _pd_background_intensity
  4.5 256.0
 40.0 158.0
 80.0  65.0
 
 loop_
 _pd_exclude_ttheta_min
 _pd_exclude_ttheta_max
 0.0 1.0
 
 _pd_instr_reflex_asymmetry_p1 0.0
 _pd_instr_reflex_asymmetry_p2 0.0
 _pd_instr_reflex_asymmetry_p3 0.0
 _pd_instr_reflex_asymmetry_p4 0.0
 
 _diffrn_radiation_polarization 0.0
 _diffrn_radiation_efficiency   1.0
 
 _pd_instr_resolution_u 16.9776
 _pd_instr_resolution_v -2.8357
 _pd_instr_resolution_w  0.5763
 _pd_instr_resolution_x  0.0
 _pd_instr_resolution_y  0.0
 
 loop_
 _pd_phase_label
 _pd_phase_scale
 _pd_phase_igsize
 Fe3O4 0.02381 0.0
 
 loop_
 _pd_meas_ttheta
 _pd_meas_intensity_up
 _pd_meas_intensity_up_sigma
 _pd_meas_intensity_down
 _pd_meas_intensity_down_sigma
 4.0 465.80 128.97 301.88 129.30
 4.2 323.78 118.22 206.06 120.00
 4.4 307.14 115.90 230.47 116.53
 
 
 loop_
 _pd_peak_index_h
 _pd_peak_index_k
 _pd_peak_index_l
 _pd_peak_mult
 _pd_peak_ttheta
 _pd_peak_intensity_up
 _pd_peak_intensity_down
 _pd_peak_width_2theta
 1 1 1  8  9.748 128.15576 128.15576 0.677
 2 0 0  6 11.260   0.00000   0.00000 0.680
 2 2 0 12 15.950  94.21107  94.21107 0.716
 
 
 loop_
 _pd_proc_ttheta
 _pd_proc_ttheta_corrected
 _pd_proc_intensity_up_net
 _pd_proc_intensity_down_net
 _pd_proc_intensity_up_total
 _pd_proc_intensity_down_total
 _pd_proc_intensity_bkg_calc
 _pd_proc_intensity_up
 _pd_proc_intensity_up_sigma
 _pd_proc_intensity_down
 _pd_proc_intensity_down_sigma
 4.000 4.385 0.00000 0.00000 256.00000 256.00000 256.00000 465.80000 128.97000 301.88000 129.30000
 4.200 4.585 0.00000 0.00000 256.00000 256.00000 256.00000 323.78000 118.22000 206.06000 120.00000
 4.400 4.785 0.00000 0.00000 256.00000 256.00000 256.00000 307.14000 115.90000 230.47000 116.53000
     """

def test_init():
    try:
        _object = Pd()
        flag = True
    except:
        flag = False
    assert flag

def test_to_cif():
    try:
        _object = Pd()
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
    _obj = Pd.from_cif(STR_FROM_CIF_1)
    assert _obj.is_defined
    assert _obj.proc.ttheta == [4.0, 4.2, 4.4]
    assert _obj.peak.index_h == [1, 2, 2]
    assert float(_obj.asymmetry.p1) == 0.0
    print(_obj.to_cif(separator="."))
