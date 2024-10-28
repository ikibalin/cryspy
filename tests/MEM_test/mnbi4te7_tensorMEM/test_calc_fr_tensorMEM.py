import os
import cryspy

def test_calc_FR_tensorMEM():
    dir = os.path.dirname(__file__)
    f_name = os.path.join(dir, "main.rcif")
    mem = cryspy.file_to_globaln(f_name)
    mem.calc_fr()
    diffrn = mem.diffrn_exp_mnbite
    refine_ls = diffrn.refine_ls
    
    assert refine_ls.goodness_of_fit_all < 4.99
    assert refine_ls.number_reflns == 142
    
    mem.create_prior_density()
