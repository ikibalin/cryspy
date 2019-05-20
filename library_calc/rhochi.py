import os
import sys
import matplotlib.pyplot
import numpy

from setup_powder_1d import *
from setup_powder_2d import *
from observed_data import *
from calculated_data import *
from experiment import *

from model import *
from variable import *
from read_rcif import *
from crystal import *


def read_model(f_name_in):
    rcif = RCif()
    rcif.load_from_file(f_name_in)
    print("\nParameters are taken from file '{:}'.".format(f_name_in))
    model = rcif.trans_to_model()
    return model

def model_refinement(model):
    model.apply_constraint()
    model.get_variables()
    d_map = {}
    res = model.refine_model(d_map)
    print("\nRefined parameters:\n")
    for var in model.get_variables():
        print(var, "", var.sigma)
    model.save_exp_mod_data()
    print("\nData are saved in the files.")

    f_out = model.get_val("f_out")
    if f_out is not None:
        s_out = model.print_report()
        fid = open(f_out, "w")
        fid.write(s_out)
        fid.close()
        print("\nListing is saved in the file '{:}'.".format(f_out))
        
def write_to_rcif(model, f_name_out):
    if f_name_out is not None:
        rcif_2 = RCif()
        rcif_2.take_from_model(model)
        rcif_2.save_to_file(f_name_out)
        print("\nParameters are saved in the file '{:}'.".format(f_name_out))
    
def rhochi_refinement(f_name_in, f_name_out):
    """
    refinement,
    parameters are defined in given .rcif fiel
    """
    print(70*"*"+"\n"+"RhoChi program. Console version.".center(70)+"\n"+
          70*"*")
    model = read_model(f_name_in)
    model_refinement(model)
    write_to_rcif(model, f_name_out)
    plot_data(model)
    
    print(70*"*"+"\n"+70*"*")

def plot_data(model):
    for experiment in model._list_experiment:
        if isinstance(experiment, ExperimentSingle):
            name = experiment.get_val("name")
            l_crystal = model._list_crystal
            observed_data = experiment.get_val("observed_data")
            h, k = observed_data.get_val("h"), observed_data.get_val("k")
            l = observed_data.get_val("l")
            f_r_exp = observed_data.get_val("flip_ratio")
            sf_r_exp = observed_data.get_val("sflip_ratio")
            
            iint_u, iint_d, f_r_mod, d_info = experiment.calc_iint_u_d_flip_ratio(h, k, l, l_crystal)

            matplotlib.pyplot.errorbar(f_r_mod, f_r_exp, sf_r_exp, fmt=".")
            matplotlib.pyplot.plot(f_r_mod, f_r_mod, "-")
            matplotlib.pyplot.xlabel("model")
            matplotlib.pyplot.ylabel("experiment")
            matplotlib.pyplot.title("experiment: '{:}', flip ratio".format(name))
            matplotlib.pyplot.show()
            
        elif isinstance(experiment, ExperimentPowder1D):
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
            
        elif isinstance(experiment, ExperimentPowder2D):
            name = experiment.get_val("name")
            l_crystal = model._list_crystal
            observed_data = experiment.get_val("observed_data")
            tth = observed_data.get_val("tth")
            phi = observed_data.get_val("phi")
            int_u_exp = observed_data.get_val("int_u")
            int_d_exp = observed_data.get_val("int_d")
            
            int_u_mod, int_d_mod, d_info = experiment.calc_profile(tth, phi, 
                                                                   l_crystal)
            i_u_exp = numpy.where(numpy.isnan(int_u_exp), 0., int_u_exp).sum(axis=1)
            i_u_mod = numpy.where(numpy.isnan(int_u_exp), 0., int_u_mod).sum(axis=1)
            i_d_exp = numpy.where(numpy.isnan(int_d_exp), 0., int_d_exp).sum(axis=1)
            i_d_mod = numpy.where(numpy.isnan(int_d_exp), 0., int_d_mod).sum(axis=1)

            matplotlib.pyplot.plot(tth, i_u_exp+i_d_exp,".",
                                   tth, i_u_mod+i_d_mod,"-",
                                   tth, i_u_exp+i_d_exp-i_u_mod-i_d_mod,"-")
            matplotlib.pyplot.xlabel("diffraction angle, degrees")
            matplotlib.pyplot.ylabel("intensity")
            matplotlib.pyplot.title("experiment: '{:}', up+down, projection".format(name))
            matplotlib.pyplot.show()

            matplotlib.pyplot.plot(tth, i_u_exp-i_d_exp,".",
                                   tth, i_u_mod-i_d_mod,"-",
                                   tth, i_u_exp-i_d_exp-i_u_mod+i_d_mod,"-")
            matplotlib.pyplot.xlabel("diffraction angle, degrees")
            matplotlib.pyplot.ylabel("intensity")
            matplotlib.pyplot.title("experiment: '{:}', up-down, projection".format(name))
            matplotlib.pyplot.show()


def create_temporary(f_name_in):
    print("Master to create .rcif file in\n'{:}'\n".format(f_name_in))
    print("""You would like to work with:
 1. single diffraction data;
 2. 1D powder diffraction data;
 3. 2D powder diffraction data. """)
    s_help = input("")
    if type(s_help) is int:
        s_help = str(s_help)
    f_dir = os.path.dirname(f_name_in)
    model = Model()

    atom_type = AtomType()
    crystal_name="Phase1"
    crystal = Crystal(name=crystal_name)
    crystal.add_atom(atom_type)
    model.add_crystal(crystal)
    
    if "1" in s_help:
        f_out = os.path.join(f_dir, "full_sd.out")
        observed_data = ObservedDataSingle(file_dir=f_dir, file_name="full_sd.dat")
        observed_data.create_input_file()
        experiment = ExperimentSingle(name="exp_sd", observed_data=observed_data,
                                      f_out=f_out)

        calculated_data = CalculatedDataSingle(name=crystal_name)
        experiment.add_calculated_data(calculated_data)
        model.add_experiment(experiment)
    if "2" in s_help:
        f_out = os.path.join(f_dir, "full_pd.out")
        observed_data = ObservedDataPowder1D(file_dir=f_dir, file_name="full_pd.dat")
        observed_data.create_input_file()
        
        background = BackgroundPowder1D(file_dir=f_dir, file_name="full_pd.bkg")
        background.create_input_file()
        
        setup = SetupPowder1D(background=background)
        
        experiment = ExperimentPowder1D(name="exp_pd", setup=setup, f_out=f_out,
                                        observed_data=observed_data)

        calculated_data = CalculatedDataPowder1D(name=crystal_name)
        experiment.add_calculated_data(calculated_data)
        model.add_experiment(experiment)
    if "3" in s_help:
        f_out = os.path.join(f_dir, "full_2dpd.out")
        observed_data = ObservedDataPowder2D(file_dir=f_dir, file_name="full_2dpd.dat")
        observed_data.create_input_file()
        
        background = BackgroundPowder2D(file_dir=f_dir, file_name="full_2dpd.bkg")
        background.create_input_file()
        
        setup = SetupPowder2D(background=background)
        
        experiment = ExperimentPowder2D(name="exp_2dpd", setup=setup, f_out=f_out, 
                                        observed_data=observed_data)

        calculated_data = CalculatedDataPowder2D(name=crystal_name)
        experiment.add_calculated_data(calculated_data)
        model.add_experiment(experiment)
    write_to_rcif(model, f_name_in)

if (__name__ == "__main__"):
    l_arg = sys.argv
    if len(l_arg) >= 2:
        f_name_in = l_arg[1]
        
        f_name_out = None
        if len(l_arg) >= 3:
            f_name_out = l_arg[2]
        if os.path.isfile(f_name_in):
            rhochi_refinement(f_name_in, f_name_out)
        else:
            create_temporary(f_name_in)
