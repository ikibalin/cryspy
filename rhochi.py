import os
import sys
import matplotlib.pyplot
import numpy

import f_experiment.f_single.cl_setup_single 
import f_experiment.f_powder_1d.cl_setup_powder_1d 
import f_experiment.f_powder_2d.cl_setup_powder_2d 

import f_experiment.f_single.cl_observed_data_single 
import f_experiment.f_single_domain.cl_observed_data_single_domain
import f_experiment.f_powder_1d.cl_observed_data_powder_1d 
import f_experiment.f_powder_2d.cl_observed_data_powder_2d

import f_experiment.f_single.cl_calculated_data_single 
import f_experiment.f_powder_1d.cl_calculated_data_powder_1d 
import f_experiment.f_powder_2d.cl_calculated_data_powder_2d 

import f_experiment.f_single.cl_experiment_single 
import f_experiment.f_single_domain.cl_experiment_single_domain 
import f_experiment.f_powder_1d.cl_experiment_powder_1d 
import f_experiment.f_powder_2d.cl_experiment_powder_2d 
import f_experiment.f_powder_texture_2d.cl_experiment_powder_texture_2d 

import f_rhochi_model.cl_model 
import f_common.cl_variable 
import f_crystal.cl_crystal 

import f_rcif.cl_rcif 
import f_api_rcif.api_rcif_model


def read_model(f_name_in):
    rcif = f_rcif.cl_rcif.RCif()
    rcif.load_from_file(f_name_in)
    print("\nParameters are taken from file '{:}'.".format(f_name_in))
    model = f_api_rcif.api_rcif_model.conv_rcif_to_model(rcif)
    return model

def model_refinement(model):
    model.apply_constraint()
    model.get_variables()
    d_info_in = {}
    res = model.refine_model(d_info_in)
    print("\nRefined parameters:\n")
    for var in model.get_variables():
        print(var)
    model.save_exp_mod_data()
    print("\nData are saved in the files.")

    file_out = model.get_val("file_out")
    file_dir = model.get_val("file_dir")
    if ((file_out is not None)|(file_dir is not None)):
        f_name = os.path.join(file_dir, file_out)
        s_out = model.print_report()
        fid = open(f_name, "w")
        fid.write(s_out)
        fid.close()
        print("\nListing is saved in the file '{:}'.".format(f_name))
        
def write_to_rcif(model, f_name_out):
    if f_name_out is not None:
        rcif_2 = f_api_rcif.api_rcif_model.conv_model_to_rcif(model)
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
        if isinstance(experiment, f_experiment.f_single.cl_experiment_single.ExperimentSingle):
            name = experiment.get_val("name")
            l_crystal = model._list_crystal
            observed_data = experiment.get_val("observed_data")
            h, k = observed_data.get_val("h"), observed_data.get_val("k")
            l = observed_data.get_val("l")
            f_r_exp = observed_data.get_val("flip_ratio")
            sf_r_exp = observed_data.get_val("sflip_ratio")
            
            as_exp = (f_r_exp-1.)/(f_r_exp+1.)
            sas_exp = ((as_exp**2+1.)**0.5)*sf_r_exp/(f_r_exp+1.)
            
            iint_u, iint_d, f_r_mod, d_info = experiment.calc_iint_u_d_flip_ratio(h, k, l, l_crystal)

            as_mod = (f_r_mod-1.)/(f_r_mod+1.)

            matplotlib.pyplot.errorbar(as_mod, as_exp, sas_exp, fmt=".")
            matplotlib.pyplot.plot(as_mod, as_mod, "-")
            matplotlib.pyplot.xlabel("model")
            matplotlib.pyplot.ylabel("experiment")
            matplotlib.pyplot.title("experiment: '{:}', assymetry".format(name))
            matplotlib.pyplot.show()

        if isinstance(experiment, f_experiment.f_single_domain.cl_experiment_single_domain.ExperimentSingleDomain):
            name = experiment.get_val("name")
            l_crystal = model._list_crystal
            
            
            d_info_in = {}
            d_info_out = experiment.calc_chi_sq(l_crystal, d_info_in)[2]

            f_r_mod = d_info_out["flip_ratio"]
            observed_data = experiment.get_val("observed_data")
            f_r_exp = observed_data.get_val("flip_ratio")
            sf_r_exp = observed_data.get_val("sflip_ratio")

            as_exp = (f_r_exp-1.)/(f_r_exp+1.)
            sas_exp = ((as_exp**2+1.)**0.5)*sf_r_exp/(f_r_exp+1.)
            
            as_mod = (f_r_mod-1.)/(f_r_mod+1.)

            matplotlib.pyplot.errorbar(as_mod, as_exp, sas_exp, fmt=".")
            matplotlib.pyplot.plot(as_mod, as_mod, "-")
            matplotlib.pyplot.xlabel("model")
            matplotlib.pyplot.ylabel("experiment")
            matplotlib.pyplot.title("experiment: '{:}', assymetry".format(name))
            matplotlib.pyplot.show()

        elif isinstance(experiment, f_experiment.f_powder_1d.cl_experiment_powder_1d.ExperimentPowder1D):
            pass
            #name = experiment.get_val("name")
            #file_dir = experiment.get_val("file_dir")
            #file_out = experiment.get_val("file_out")
            #f_name = os.path.join(file_dir, file_out)
            #fid = open(f_name, 'r')
            #l_cont = fid.readlines()
            #fid.close
            #l_name = l_cont[0].strip().split()
            #dd = {}
            #for hh in l_name:
            #    dd[hh] = []
            #for line in l_cont[1:]:
            #    if line.strip()=="":
            #        break
            #    for hh_1, hh_2 in zip(l_name, line.strip().split()):
            #        dd[hh_1].append(float(hh_2))
            #for hh in l_name:
            #    dd[hh] = numpy.array(dd[hh], dtype=float)
            #    
            #matplotlib.pyplot.plot(dd["ttheta"],dd["exp_up"]+dd["exp_down"],".",
            #                   dd["ttheta"],dd["mod_up"]+dd["mod_down"],"-",
            #                   dd["ttheta"],dd["exp_up"]+dd["exp_down"]-
            #                                dd["mod_up"]-dd["mod_down"],"-")
            #matplotlib.pyplot.xlabel("diffraction angle, degrees")
            #matplotlib.pyplot.ylabel("intensity")
            #matplotlib.pyplot.title("experiment: '{:}', up+down".format(name))
            #matplotlib.pyplot.show()
            #matplotlib.pyplot.plot(dd["ttheta"],dd["exp_up"]-dd["exp_down"],".",
            #                   dd["ttheta"],dd["mod_up"]-dd["mod_down"],"-",
            #                   dd["ttheta"],dd["exp_up"]-dd["exp_down"]-
            #                                dd["mod_up"]+dd["mod_down"],"-")
            #matplotlib.pyplot.xlabel("diffraction angle, degrees")
            #matplotlib.pyplot.ylabel("intensity")
            #matplotlib.pyplot.title("experiment: '{:}', up-down".format(name))
            #matplotlib.pyplot.show()
            
        elif isinstance(experiment, f_experiment.f_powder_2d.cl_experiment_powder_2d.ExperimentPowder2D):
            name = experiment.get_val("name")
            l_crystal = model._list_crystal
            observed_data = experiment.get_val("observed_data")
            tth = observed_data.get_val("tth")
            phi = observed_data.get_val("phi")
            int_u_exp = observed_data.get_val("int_u")
            int_d_exp = observed_data.get_val("int_d")

            d_info_in = {}
            int_u_mod, int_d_mod, d_info_out = experiment.calc_profile(
                                         tth, phi, l_crystal, d_info_in)
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

        elif isinstance(experiment, f_experiment.f_powder_texture_2d.cl_experiment_powder_texture_2d.ExperimentPowderTexture2D):
            name = experiment.get_val("name")
            l_crystal = model._list_crystal
            observed_data = experiment.get_val("observed_data")
            tth = observed_data.get_val("tth")
            phi = observed_data.get_val("phi")
            int_u_exp = observed_data.get_val("int_u")
            int_d_exp = observed_data.get_val("int_d")

            d_info_in = {}
            int_u_mod, int_d_mod, d_info_out = experiment.calc_profile(
                                         tth, phi, l_crystal, d_info_in)
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
    model = f_rhochi_model.cl_model.Model()

    atom_type_1 = f_crystal.cl_crystal.AtomType(flag_m=True)
    crystal_name = "Phase1"
    crystal = f_crystal.cl_crystal.Crystal(name=crystal_name)
    crystal.add_atom(atom_type_1)
    model.add_crystal(crystal)
    
    if "1" in s_help:
        file_out = "full_sd.out"
        observed_data = f_experiment.f_single.cl_observed_data_single.ObservedDataSingle(file_dir=f_dir, file_name="full_sd.dat")
        observed_data.create_input_file()
        experiment = f_experiment.f_single.cl_experiment_single.ExperimentSingle(name="exp_sd", observed_data=observed_data,
                                      file_out=file_out, file_dir=f_dir)

        calculated_data = f_experiment.f_single.cl_calculated_data_single.CalculatedDataSingle(name=crystal_name)
        experiment.add_calculated_data(calculated_data)
        model.add_experiment(experiment)

    if "2" in s_help:
        file_out = "full_pd.out"
        observed_data = f_experiment.f_powder_1d.cl_observed_data_powder_1d.ObservedDataPowder1D(file_dir=f_dir, file_name="full_pd.dat")
        observed_data.create_input_file()
        
        background = f_experiment.f_powder_1d.cl_setup_powder_1d.BackgroundPowder1D(file_dir=f_dir, file_name="full_pd.bkg")
        background.create_input_file()
        
        setup = f_experiment.f_powder_1d.cl_setup_powder_1d.SetupPowder1D(background=background)
        
        experiment = f_experiment.f_powder_1d.cl_experiment_powder_1d.ExperimentPowder1D(name="exp_pd", setup=setup, file_dir=f_dir, 
                        file_out=file_out, observed_data=observed_data, excl_tth_min=[35], excl_tth_max=[40])

        calculated_data = f_experiment.f_powder_1d.cl_calculated_data_powder_1d.CalculatedDataPowder1D(name=crystal_name)
        experiment.add_calculated_data(calculated_data)
        model.add_experiment(experiment)

    if "3" in s_help:
        file_out = "full_2dpd.out"
        observed_data = f_experiment.f_powder_2d.cl_observed_data_powder_2d.ObservedDataPowder2D(file_dir=f_dir, file_name="full_2dpd.dat")
        observed_data.create_input_file()
        
        background = f_experiment.f_powder_2d.cl_setup_powder_2d.BackgroundPowder2D(file_dir=f_dir, file_name="full_2dpd.bkg")
        background.create_input_file()
        
        setup = f_experiment.f_powder_2d.cl_setup_powder_2d.SetupPowder2D(background=background)
        
        experiment = f_experiment.f_powder_2d.cl_experiment_powder_2d.ExperimentPowder2D(name="exp_2dpd", setup=setup, file_out=file_out,
                               file_dir=f_dir, observed_data=observed_data)

        calculated_data = f_experiment.f_powder_2d.cl_calculated_data_powder_2d.CalculatedDataPowder2D(name=crystal_name)
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
