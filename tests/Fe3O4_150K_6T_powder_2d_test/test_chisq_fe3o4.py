import os
import cryspy

def test_UNPD_PbSo4():
    dir = os.path.dirname(__file__)
    f_name = os.path.join(dir, "main.rcif")
    rhochi = cryspy.file_to_globaln(f_name)
    chi_sq, n_points = rhochi.calc_chi_sq()
    assert chi_sq < 55600.9
    assert int(n_points) == 25553
