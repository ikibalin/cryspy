import os
import cryspy

def test_UNPD_PbSo4():
    dir = os.path.dirname(__file__)
    f_name = os.path.join(dir, "main.rcif")
    rhochi = cryspy.file_to_globaln(f_name)
    cryspy.rhochi_run(rhochi)
    d_out = cryspy.rhochi_no_refinement(rhochi)
    chi_sq = d_out["chi_sq"]
    n_points = d_out["n_points"]    
    assert chi_sq < 1136 # 981.7
    assert int(n_points) == 94
