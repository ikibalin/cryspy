import os
import sys
import matplotlib.pyplot
import numpy

from experiment import *
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
    model.apply_constraint()
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
    
    for experiment in model._list_experiment:
        if isinstance(experiment, ExperimentPowder1D):
            name = experiment.get_val("name")
            f_out = experiment.get_val("f_out")
        
            fid = open(f_out, 'r')
            l_cont = fid.readlines()
            fid.close
            l_name = l_cont[0].strip().split()
            dd = {}
            for hh in l_name:
                dd[hh] = []
            for line in l_cont[1:]:
                for hh_1, hh_2 in zip(l_name, line.strip().split()):
                    dd[hh_1].append(float(hh_2))
            for hh in l_name:
                dd[hh] = numpy.array(dd[hh], dtype=float)
            matplotlib.pyplot.plot(dd["ttheta"],dd["exp_up"]+dd["exp_down"],".",
                               dd["ttheta"],dd["mod_up"]+dd["mod_down"],"-",
                               dd["ttheta"],dd["exp_up"]+dd["exp_down"]-
                                            dd["mod_up"]-dd["mod_down"],"-")
            matplotlib.pyplot.xlabel("diffraction angle, degrees")
            matplotlib.pyplot.ylabel("intensity")
            matplotlib.pyplot.title("experiment: '{:}', up+down".format(name))
            matplotlib.pyplot.show()
            matplotlib.pyplot.plot(dd["ttheta"],dd["exp_up"]-dd["exp_down"],".",
                               dd["ttheta"],dd["mod_up"]-dd["mod_down"],"-",
                               dd["ttheta"],dd["exp_up"]-dd["exp_down"]-
                                            dd["mod_up"]+dd["mod_down"],"-")
            matplotlib.pyplot.xlabel("diffraction angle, degrees")
            matplotlib.pyplot.ylabel("intensity")
            matplotlib.pyplot.title("experiment: '{:}', up-down".format(name))
            matplotlib.pyplot.show()


if (__name__ == "__main__"):
    l_arg = sys.argv
    if len(l_arg) >= 2:
        f_name_in = l_arg[1]
        
        f_name_out = None
        if len(l_arg) >= 3:
            f_name_out = l_arg[2]
        if os.path.isfile(f_name_in):
            rhochi_refinement(f_name_in, f_name_out)
