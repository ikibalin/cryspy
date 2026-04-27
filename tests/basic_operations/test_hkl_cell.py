import os
import math
import numpy

import cryspy
from cryspy.B_parent_classes.cl_3_hkl_cell import calc_d_sthovl_for_hkl

s_cont_l = """
 loop_
 _pd_peak_index_h
 _pd_peak_index_k
 _pd_peak_index_l
 _pd_peak_index_mult
 _pd_peak_2theta
 _pd_peak_intensity_plus
 _pd_peak_intensity_minus
 _pd_peak_width_2theta
  2  2  0  4 17.2 100.0  90.0  2.3
"""

s_cont_i = """
 _pd_peak_index_h 2
 _pd_peak_index_k 2
 _pd_peak_index_l 0
"""

def test_hkl_cell():
    obj_hkl_0 = cryspy.PdPeakL.from_cif(s_cont_l) 
    obj_hkl = cryspy.PdPeak.from_cif(s_cont_i)

    cell = cryspy.Cell(
        length_a=5, length_b=5, length_c=5, 
        angle_alpha=90, angle_beta=90, angle_gamma=90
    )
    calc_d_sthovl_for_hkl(cell, obj_hkl_0)
    calc_d_sthovl_for_hkl(cell, obj_hkl)

    assert math.isclose(obj_hkl.d_spacing, 0.5*3.53553, rel_tol=1e-4)
    assert math.isclose(obj_hkl.sintlambda, 0.28284, rel_tol=1e-4)