import os 

import cryspy

def test_mem_magnetization_density_for_momemts():
    f_name = os.path.join(os.path.dirname(__file__), "main.rcif")
    rcif_obj = cryspy.load_file(f_name)
    # cryspy.rhochi_rietveld_refinement(rcif_obj)
    res =cryspy.mempy_magnetization_density_reconstruction(rcif_obj)
    # cryspy.mempy_spin_density_reconstruction(rcif_obj)
    assert res < 2.76