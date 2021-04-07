import os
import cryspy

def test_FeMn_zero():
    dir = os.path.dirname(__file__)
    f_name = os.path.join(dir, "main.rcif")
    rhochi = cryspy.file_to_globaln(f_name)
    chi_sq, n_points = rhochi.calc_chi_sq()
    assert chi_sq < 463.0
    assert int(n_points) == 1091
