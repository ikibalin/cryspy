import os 

import cryspy

def test_mem_magnetization_density_for_chi():
    f_name = os.path.join(os.path.dirname(__file__), "mempy_magnetization_density_Yb2Ti2O7_2K_1T.rcif")
    rcif_obj = cryspy.load_file(f_name)
    # cryspy.rhochi_rietveld_refinement(rcif_obj)
    chi_sq = cryspy.mempy_magnetization_density_reconstruction(rcif_obj)
    assert chi_sq  < 7.5
    # cryspy.mempy_spin_density_reconstruction(rcif_obj)
