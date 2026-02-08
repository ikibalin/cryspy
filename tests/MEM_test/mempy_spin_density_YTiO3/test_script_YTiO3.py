import os 

import cryspy

def test_mem_spin_density():
    f_name = os.path.join(os.path.dirname(__file__), "mempy_spin_density_YTiO3.rcif")
    rcif_obj = cryspy.load_file(f_name)
    # cryspy.rhochi_rietveld_refinement(rcif_obj)
    # cryspy.mempy_magnetization_density_reconstruction(rcif_obj)
    chi_sq = cryspy.mempy_spin_density_reconstruction(rcif_obj)
    assert chi_sq < 17.89
