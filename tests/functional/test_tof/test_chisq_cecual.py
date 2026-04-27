import os
import cryspy

def test_UNPD_PbSo4_pV():
    dir = os.path.dirname(__file__)
    f_name = os.path.join(dir, "model_pV.rcif")
    rhochi = cryspy.file_to_globaln(f_name)
    res = cryspy.rhochi_no_refinement(rhochi)
    chi_sq, n_point = res["chi_sq"], res["n_point"]
    assert chi_sq < 41558
    assert int(n_point) == 3692

def test_UNPD_PbSo4_gauss():
    dir = os.path.dirname(__file__)
    f_name = os.path.join(dir, "model_gauss.rcif")
    rhochi = cryspy.file_to_globaln(f_name)
    res = cryspy.rhochi_no_refinement(rhochi)
    chi_sq, n_point = res["chi_sq"], res["n_point"]
    assert chi_sq < 41558
    assert int(n_point) == 3692

