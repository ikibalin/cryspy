import os
import sys

from model import *
from variable import *
from read_rcif import *

from crystal import *

f_name = r"C:\Users\yurik\Documents\GitHub\rhochi\library_calc\TbMnO3\full.rcif"
def rhochi_refinement(f_name):
    """
    refinement,
    parameters are defined in given .rcif fiel
    """
    rcif = RCif()
    rcif.load_from_file(f_name)
    model = rcif.trans_to_model()
    model.get_variables()
    d_map = {}
    res = model.refine_model(d_map)
    print("\nRefined parameters:\n")
    for var in model.get_variables():
        print(var)
    model.save_exp_mod_data()
    print("Data are saved in the files")
    
    
if (__name__ == "__main__"):
    l_arg = sys.argv
    if len(l_arg) >= 2:
        f_name = l_arg[1]
        if os.path.isfile(f_name):
            rhochi_refinement(f_name)

