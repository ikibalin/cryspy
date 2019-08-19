"""
transform inforamtion from cif file to class phase
"""
import sys
import os


from neupy.f_rcif.cl_rcif import RCif 

from neupy.f_crystal import (
    Crystal,
    Cell,
    SpaceGroup,
    AtomSite,
    AtomType
    )

from neupy.f_experiment import (
    BeamPolarization,
    CalculatedDataPowder1D,
    CalculatedDataPowder2D,
    CalculatedDataSingle,
    AsymmetryPowder1D,
    AsymmetryPowder2D,
    BackgroundPowder1D,
    BackgroundPowder2D,
    FactorLorentzPowder1D,
    FactorLorentzPowder2D,
    ResolutionPowder1D,
    ResolutionPowder2D,
    SetupPowder1D,
    SetupPowder2D,
    SetupSingle,
    ObservedDataPowder1D,
    ObservedDataPowder2D,
    ObservedDataSingle,
    ObservedDataSingleDomain,
    ExperimentPowder1D,
    ExperimentPowder2D,
    ExperimentPowderTexture2D,
    ExperimentSingle,
    ExperimentSingleDomain,
    Model,
    CellDensity,
    MemReconstruction,
    ObservedDataMem
    )


from neupy.f_common.cl_variable import Variable

from neupy.f_api_rcif.api_rcif_crystal import conv_data_to_crystal
from neupy.f_api_rcif.api_rcif_common import from_dict_to_obj


def conv_rcif_to_mem_reconstruction(rcif):
    p_glob = rcif.glob
    l_data = p_glob["data"]
    f_dir = rcif._p_file_dir

    l_crystal = conv_data_to_crystal(l_data)

    for data in l_data:
        l_key = data.keys()

        lab_mem = "_mem_"

        flag_mem = any([hh.startswith(lab_mem) for hh in l_key])
       
        if flag_mem:
            cell_density = conv_data_to_cell_density(data, f_dir)
            l_observed_data_mem = conv_data_to_observed_data_mem(data, f_dir)
            break

    name_p_glob = p_glob["name"]
    l_key_g = p_glob.keys()

    mem_reconstruction = MemReconstruction(name=name_p_glob, file_dir=f_dir, cell_density=cell_density)
    if "_file_name_output_listing" in l_key_g:
        file_out = p_glob["_file_name_output_listing"]
        mem_reconstruction.set_val(file_out=file_out)

    for crystal in l_crystal:
        mem_reconstruction.add_crystal(crystal)

    for observed_data_mem in l_observed_data_mem: 
        mem_reconstruction.add_observed_data_mem(observed_data_mem)
    
    return mem_reconstruction

def conv_data_to_cell_density(data, f_dir):
    cell_density = CellDensity(file_dir=f_dir)
    l_relation = data_cell_density_relation()
    from_dict_to_obj(data, l_relation, cell_density)
    return cell_density


def conv_data_to_observed_data_mem(data, f_dir):
    l_observed_data_mem = []

    l_key = data.keys()
    if (not ("loops" in l_key)):
        return l_observed_data_mem 

    l_relation = data_observed_data_mem_relation()
    for loop in data["loops"]:
        l_key = list(loop.keys())
        flag_name = "_mem_phase_name" in l_key


        if flag_name:
            n_crystal = len(loop["_mem_phase_name"])
            
            for i_crystal in range(n_crystal):
                dd = {}
                for key in l_key:
                    dd.update({key: loop[key][i_crystal]})
                
                observed_data_mem = ObservedDataMem(file_dir=f_dir) 
                l_relation = data_observed_data_mem_relation()
                from_dict_to_obj(dd, l_relation, observed_data_mem)
                l_observed_data_mem.append(observed_data_mem)     
    return l_observed_data_mem


def conv_mem_reconstruction_to_rcif(mem_reconstruction):

    def temp_func(ddata, relation, obj):
        lab_d, lab_o, type_val = relation[:3]
        key_d = lab_d
        
        if type_val == "logic":
            sval = "True" if obj.get_val(lab_o) else "False"
            ddata[key_d] = sval 
        elif type_val == "text":
            ddata[key_d] = obj.get_val(lab_o)
        elif type_val == "val":
            val = obj.get_val(lab_o)
            if isinstance(val, Variable):
                sval = "{:}".format(val.print_with_sigma())
            else:
                sval = "{:}".format(val)
            ddata[key_d] = sval 
        else:
            print(50*"ERROR ")
            print("type_val: ", type_val)
        return
        
    l_data = conv_crystal_to_data(mem_reconstruction._list_crystal)

    cell_density = mem_reconstruction.get_val("cell_density")
    dd = {"name": cell_density.get_val("name")}
    l_relation = data_cell_density_relation()
    for relation in l_relation:
        temp_func(dd, relation, cell_density)

    dd_2 = {}
    l_relation = data_observed_data_mem_relation()
    l_dhelp = []
    for observed_data_mem in mem_reconstruction._list_observed_data_mem:
        name = observed_data_mem.get_val("name")
        dhelp = {"_mem_phase_name": name}
        for relation in l_relation:
            temp_func(dhelp, relation, observed_data_mem)
        l_dhelp.append(dhelp)

    llab_dhelp = list(dhelp.keys())
    for lab_d in llab_dhelp:
        dd_2[lab_d] = [hh[lab_d] for hh in l_dhelp]

    dd.update({"loops":[dd_2]})
    l_data.append(dd)

    rcif = cl_rcif.RCif()
    rcif.glob["data"] = l_data
    rcif.glob["name"] = mem_reconstruction.get_val("name")

    file_out = mem_reconstruction.get_val("file_out")
    if file_out is None:
        file_out = "full.lis"
    rcif.glob["_file_name_output_listing"] = os.path.basename(file_out)
    return rcif

def conv_crystal_to_data(crystal):
    data = {}
    return data


def data_cell_density_relation():
    l_relation = [("name", "name", "text"),
        ("_mem_points_number_a", "points_number_a", "val"),
        ("_mem_points_number_b", "points_number_b", "val"),
        ("_mem_points_number_c", "points_number_c", "val"),
        ("_mem_file_name_output", "file_name", "text")]
    return l_relation



def data_observed_data_mem_relation():
    l_relation = [("_mem_beam_polarization_up", "beam_polarization_up", "val"),
                  ("_mem_beam_flipper_efficiency", "flipper_efficiency", "val"),
                  ("_mem_file_name_input", "file_name", "text"),
                  ("_mem_phase_extinction_mosaicity", "extinction_mosaicity", "val"),
                  ("_mem_phase_extinction_radius", "extinction_radius", "val"),
                  ("_mem_phase_name", "name", "text")]
    return l_relation




def main(larg):
    if len(larg) > 1:
        fname = "powder.rcif"
        fname = larg[1]
    else:
        fname = input("What is the name of the 'rcif' file?\n")
        
    
    

if __name__ == "__main__":
    main(sys.argv)
