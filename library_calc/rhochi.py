import os
import sys

from model import *
from variable import *
from read_rcif import *

from crystal import *

def rhochi_refinement(f_name_in, f_name_out):
    """
    refinement,
    parameters are defined in given .rcif fiel
    """
    rcif = RCif()
    rcif.load_from_file(f_name_in)
    print("Parameters are taken from file '{:}'".format(f_name_in))
    model = rcif.trans_to_model()
    model.get_variables()
    d_map = {}
    res = model.refine_model(d_map)
    print("\nRefined parameters:\n")
    for var in model.get_variables():
        print(var)
    model.save_exp_mod_data()
    print("Data are saved in the files")

    if f_name_out is not None:
        rcif_2 = RCif()
        rcif_2.take_from_model(model)
        rcif_2.save_to_file(f_name_out)
        print("Parameters are saved in file '{:}'".format(f_name_out))
    
    
if (__name__ == "__main__"):
    l_arg = sys.argv
    if len(l_arg) >= 2:
        f_name_in = l_arg[1]
        
        f_name_out = None
        if len(l_arg) >= 3:
            f_name_out = l_arg[2]
        if os.path.isfile(f_name_in):
            rhochi_refinement(f_name_in, f_name_out)
