import os
import cryspy

def test_UNPD_PbSo4():
    dir = os.path.dirname(__file__)
    f_name = os.path.join(dir, "model_pV.rcif")
    rhochi = cryspy.file_to_globaln(f_name)
    chi_sq, n_points = rhochi.calc_chi_sq()
    assert chi_sq < 41558
    assert int(n_points) == 3692
    f_name = os.path.join(dir, "model_gauss.rcif")
    rhochi = cryspy.file_to_globaln(f_name)
    chi_sq, n_points = rhochi.calc_chi_sq()
    assert chi_sq < 41558
    assert int(n_points) == 3692

