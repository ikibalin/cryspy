import os
import cryspy 



f_dir_examples = os.path.dirname(__file__)
f_name = os.path.join(f_dir_examples, "mcif_chi_even_hkl.rcif")
rcif_obj = cryspy.load_file(f_name)
cryspy.rhochi_no_refinement(rcif_obj)
pass