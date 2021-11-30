import os
import cryspy


def test_UNPD_PbSo4():
    dir = os.path.dirname(__file__)
    f_name = os.path.join(dir, "main.rcif")
    rhochi = cryspy.load_file(f_name)
    d_out = cryspy.rhochi_no_refinement(rhochi)
    chi_sq = d_out["chi_sq"]
    n_points = d_out["n_points"]    
    assert chi_sq < 4427.6
    assert int(n_points) == 756
