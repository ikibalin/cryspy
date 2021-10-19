import os
import cryspy

def test_calc_FR_2channel():
    dir = os.path.dirname(__file__)
    f_name = os.path.join(dir, "main.rcif")
    mem = cryspy.file_to_globaln(f_name)
    res = cryspy.mempy_density(mem)
    # mem.calc_fr()
    diffrn = mem.diffrn_exp_mnbite
    refine_ls = diffrn.refine_ls
    print(refine_ls)
    
test_calc_FR_2channel()