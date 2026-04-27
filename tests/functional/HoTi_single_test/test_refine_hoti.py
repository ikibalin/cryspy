import os
import cryspy

def test_UNPD_PbSo4():
    dir = os.path.dirname(__file__)
    f_name = os.path.join(dir, "main.rcif")
    rhochi = cryspy.file_to_globaln(f_name)
    cryspy.rhochi_rietveld_refinement( rhochi )
    res = cryspy.rhochi_no_refinement(rhochi)
    chi_sq = res["chi_sq"]
    n_points = res["n_point"]
    assert chi_sq < 1136
    assert int(n_points) == 94
