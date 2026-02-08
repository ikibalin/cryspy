import os
import cryspy

def test_refine():
    dir = os.path.dirname(__file__)
    f_name = os.path.join(dir, "main.rcif")
    rhochi = cryspy.file_to_globaln(f_name)
    d_out = cryspy.rhochi_no_refinement(rhochi)
    chi_sq = d_out["chi_sq"]
    n_points = d_out["n_point"]
    assert chi_sq < 439
    assert int(n_points) == 142

def test_estimate_fm():
    dir = os.path.dirname(__file__)
    f_name = os.path.join(dir, "main.rcif")
    rhochi = cryspy.file_to_globaln(f_name)
    diffrn = rhochi.diffrn_exp_mnbite 
    crystal = rhochi.crystal_mnbite
    estim = diffrn.estimate_FM(crystal)
    l_f_m = estim.f_m
    l_f_m_sigma = estim.f_m_sigma
    assert abs(float(l_f_m[0])-(-0.7881)) < 0.001
    assert abs(float(l_f_m[1])-(-0.2003)) < 0.001
    assert abs(float(l_f_m[2])-(0.68611)) < 0.001
    assert abs(float(l_f_m[3])-(0.41953)) < 0.001

    assert abs(float(l_f_m_sigma[0]) - 0.16078) < 0.001
    assert abs(float(l_f_m_sigma[1]) - 0.14083) < 0.001
    assert abs(float(l_f_m_sigma[2]) - 0.47949) < 0.001
    assert abs(float(l_f_m_sigma[3]) - 0.34928) < 0.001

