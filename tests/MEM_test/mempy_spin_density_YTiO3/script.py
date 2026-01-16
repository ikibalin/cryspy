import os 

import cryspy


f_name = os.path.join(os.path.dirname(__file__), "mempy_spin_density_YTiO3.rcif")
rcif_obj = cryspy.load_file(f_name)
# cryspy.rhochi_rietveld_refinement(rcif_obj)
# cryspy.mempy_magnetization_density_reconstruction(rcif_obj)
cryspy.mempy_spin_density_reconstruction(rcif_obj)
