import os
import sys
import numpy

from neupy import (
    BackgroundPowder1D,
    SetupPowder1D,
    ObservedDataPowder1D,
    CalculatedDataPowder1D,
    ExperimentPowder1D,
    BackgroundPowder2D,
    SetupPowder2D,
    ObservedDataPowder2D,
    CalculatedDataPowder2D,
    ExperimentPowder2D,
    ExperimentPowderTexture2D,
    ExperimentSingle,
    ObservedDataSingle,
    CalculatedDataSingle,
    ExperimentSingleDomain,
    Cell,
    SpaceGroup,
    AtomType,
    AtomSite,
    Crystal,
    Model,
    Variable,
    RCif,
    conv_rcif_to_model,
    conv_model_to_rcif
    )





def read_model(f_name_in):
    rcif = RCif()
    rcif.load_from_file(f_name_in)
    print("\nParameters are taken from file '{:}'.".format(f_name_in))
    model = conv_rcif_to_model(rcif)
    return model

def model_refinement(model):
    model.apply_constraint()
    model.get_variables()
    d_info_in = {}
    model.refine_model(d_info_in)
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
        rcif_2 = conv_model_to_rcif(model)
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
            pass

        if isinstance(experiment, ExperimentSingleDomain):
            pass

        elif isinstance(experiment, ExperimentPowder1D):
            pass
            
        elif isinstance(experiment, ExperimentPowder2D):
            pass

        elif isinstance(experiment, ExperimentPowderTexture2D):
            pass

def create_temporary(f_name_in, exp_type="2"):
    f_dir = os.path.dirname(f_name_in)
    model = Model()

    atom_type_1 = AtomType(flag_m=True)
    crystal_name = "Phase1"
    crystal = Crystal(name=crystal_name)
    crystal.add_atom(atom_type_1)
    model.add_crystal(crystal)
    
    if "1" in exp_type:
        file_out = "full_sd.out"
        observed_data = ObservedDataSingle(file_dir=f_dir, file_name="full_sd.dat")
        observed_data.create_input_file()
        experiment = ExperimentSingle(name="exp_sd", observed_data=observed_data,
                                      file_out=file_out, file_dir=f_dir)

        calculated_data = CalculatedDataSingle(name=crystal_name)
        experiment.add_calculated_data(calculated_data)
        model.add_experiment(experiment)

    if "2" in exp_type:
        file_out = "full_pd.out"
        observed_data = ObservedDataPowder1D(file_dir=f_dir, 
                        file_name="full_pd.dat", tth_min=2., tth_max=100)
        observed_data.create_input_file()
        
        background = BackgroundPowder1D(file_dir=f_dir, file_name="full_pd.bkg")
        background.create_input_file()
        
        setup = SetupPowder1D(background=background)
        
        experiment = ExperimentPowder1D(name="exp_pd", setup=setup, file_dir=f_dir, 
                        file_out=file_out, observed_data=observed_data, excl_tth_min=[35], excl_tth_max=[40])

        calculated_data = CalculatedDataPowder1D(name=crystal_name)
        experiment.add_calculated_data(calculated_data)
        model.add_experiment(experiment)

    if "3" in exp_type:
        file_out = "full_2dpd.out"
        observed_data = ObservedDataPowder2D(file_dir=f_dir, file_name="full_2dpd.dat",
        tth_min=2., tth_max=80, phi_min=0., phi_max=20.)
        observed_data.create_input_file()
        
        background = BackgroundPowder2D(file_dir=f_dir, file_name="full_2dpd.bkg")
        background.create_input_file()
        
        setup = SetupPowder2D(background=background)
        
        experiment = ExperimentPowder2D(name="exp_2dpd", setup=setup, file_out=file_out,
                               file_dir=f_dir, observed_data=observed_data)

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
            print("Master to create .rcif file in\n'{:}'\n".format(f_name_in))
            print("""You would like to work with:
 1. single diffraction data;
 2. 1D powder diffraction data;
 3. 2D powder diffraction data. """)
            exp_type = input("")
            if type(exp_type) is int:
                exp_type = str(exp_type)
            create_temporary(f_name_in, exp_type)
