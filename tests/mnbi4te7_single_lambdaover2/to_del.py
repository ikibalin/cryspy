import os
import cryspy

def test_refine():
    dir = os.path.dirname(__file__)
    f_name = os.path.join(dir, "main.rcif")
    rhochi = cryspy.file_to_globaln(f_name)
    rhochi.refine()
    chi_sq, n_points = rhochi.calc_chi_sq()

def test_estimate_fm():
    dir = os.path.dirname(__file__)
    f_name = os.path.join(dir, "main.rcif")
    rhochi = cryspy.file_to_globaln(f_name)
    l_estimation = rhochi.estimate_f_mag_for_diffrn()
    estim = l_estimation[0]
    l_f_m = estim.f_m
    l_f_m_sigma = estim.f_m_sigma

test_refine()
test_estimate_fm()