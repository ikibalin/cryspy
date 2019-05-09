__author__ = 'ikibalin'
__version__ = "2019_05_09"
import os
import sys

from model import *

from variable import *
from read_rcif import *

f_name = r"C:\Users\ikibalin\Documents\GitHub\rhochi\library_calc\Fe3O4_150K_6T_2d\full.rcif"
def rhochi_refinement(f_name):
    """
    refinement,
    parameters are defined in given .rcif fiel
    """
    rcif = RCif()
    rcif.load_from_file(f_name)
    model = rcif.trans_to_model()
    
    d_map = {}
    res = model.refine_model(d_map)
    print("\nRefined parameters:\n")
    for var in model.get_variables():
        print(var)
    
if (__name__ == "__main__"):
    l_arg = sys.argv
    
    if len(l_arg) >= 2:
        f_name = l_arg[1]
        if os.path.isfile(f_name):
            rhochi_refinement(f_name)
