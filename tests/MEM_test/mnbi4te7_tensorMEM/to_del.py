import os
import cryspy

def test_calc_FR_tensorMEM():
    dir = os.path.dirname(__file__)
    f_name = os.path.join(dir, "main.rcif")
    mem = cryspy.file_to_globaln(f_name)
    mem.calc_fr()
    diffrn = mem.diffrn_exp_mnbite
    refine_ls = diffrn.refine_ls
    
    mem.create_prior_density()

test_calc_FR_tensorMEM()
print("Done.")