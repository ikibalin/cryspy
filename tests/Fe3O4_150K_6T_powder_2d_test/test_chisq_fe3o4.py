import os
import cryspy

def test_powder2D_Fe3O4():
    dir = os.path.dirname(__file__)
    f_name = os.path.join(dir, "Fe3O4_150K_6.0T.rcif")
    rhochi = cryspy.file_to_globaln(f_name)
    d_out = cryspy.rhochi_no_refinement(rhochi)
    #d_out = cryspy.rhochi_rietveld_refinement(rhochi)
    chi_sq = d_out["chi_sq"]
    n_point = d_out["n_point"]
    # chi_sq, n_p
    # oints = rhochi.calc_chi_sq()
    assert chi_sq < 299963
    assert int(n_point) == 153955