import os 

import cryspy


f_name = os.path.join(os.path.dirname(__file__), "main.rcif")
rcif_obj = cryspy.load_file(f_name)

cryspy.mempy_magnetization_density_reconstruction(rcif_obj)