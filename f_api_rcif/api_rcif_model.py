"""
transform inforamtion from cif file to class phase
"""
import sys
import os

import f_rcif.cl_rcif 

import f_crystal.cl_crystal

import f_experiment.f_single.cl_calculated_data_single 
import f_experiment.f_powder_1d.cl_calculated_data_powder_1d 
import f_experiment.f_powder_2d.cl_calculated_data_powder_2d 

import f_experiment.f_single.cl_setup_single 
import f_experiment.f_powder_1d.cl_setup_powder_1d 
import f_experiment.f_powder_2d.cl_setup_powder_2d 

import f_experiment.f_single.cl_observed_data_single 
import f_experiment.f_powder_1d.cl_observed_data_powder_1d 
import f_experiment.f_powder_2d.cl_observed_data_powder_2d 
import f_experiment.f_single_domain.cl_observed_data_single_domain 

import f_experiment.f_single.cl_experiment_single 
import f_experiment.f_powder_1d.cl_experiment_powder_1d 
import f_experiment.f_powder_2d.cl_experiment_powder_2d 
import f_experiment.f_single_domain.cl_experiment_single_domain 

import f_experiment.f_powder_texture_2d.cl_experiment_powder_texture_2d 

import f_rhochi_model.cl_model 
import f_common.cl_variable 


def conv_rcif_to_model(rcif):
    p_glob = rcif.glob
    l_data = p_glob["data"]
    f_dir = rcif._p_file_dir

    l_crystal = []
    for data in l_data:
        l_key = data.keys()

        lab_cry = "_cell_"

        flag_cry = any([hh.startswith(lab_cry) for hh in l_key])
        if flag_cry:
            crystal = conv_data_to_crystal(data)
            l_crystal.append(crystal)

    l_experiment = []
    for data in l_data:
        l_key = data.keys()

        lab_exp_pd_1d, lab_exp_pd_2d = "_pd_", "_2dpd_"
        lab_exp_pdt_2d = "_2dpdt_"
        lab_exp_sd, lab_exp_sdd = "_sd_", "_sdd_"

        flag_exp_pd = any([hh.startswith(lab_exp_pd_1d) for hh in l_key])
        flag_exp_pd_2d = any([hh.startswith(lab_exp_pd_2d) for hh in l_key])
        flag_exp_pdt_2d = any([hh.startswith(lab_exp_pdt_2d) for hh in l_key])
        flag_exp_sd = any([hh.startswith(lab_exp_sd) for hh in l_key])
        flag_exp_sdd = any([hh.startswith(lab_exp_sdd) for hh in l_key])
        
        if flag_exp_pd:
            experiment = conv_data_to_experiment_powder_1d(data, f_dir, l_crystal)
            l_experiment.append(experiment)
        elif flag_exp_pd_2d:
            experiment = conv_data_to_experiment_powder_2d(data, f_dir, l_crystal)
            l_experiment.append(experiment)
        elif flag_exp_pdt_2d:
            experiment = conv_data_to_experiment_powder_texture_2d(data, f_dir, l_crystal)
            l_experiment.append(experiment)
        elif flag_exp_sd:
            experiment = conv_data_to_experiment_single(data, f_dir, l_crystal)
            l_experiment.append(experiment)
        elif flag_exp_sdd:
            experiment = conv_data_to_experiment_single_domain(data, f_dir, l_crystal)
            l_experiment.append(experiment)


    name_p_glob = p_glob["name"]
    l_key_g = p_glob.keys()

    model = f_rhochi_model.cl_model.Model(name=name_p_glob, file_dir=f_dir)
    if "_file_name_output_listing" in l_key_g:
        file_out = p_glob["_file_name_output_listing"]
        model.set_val(file_out=file_out)

    for experiment in l_experiment:
        model.add_experiment(experiment)
    for crystal in l_crystal:
        model.add_crystal(crystal) 
            
    return model

def conv_data_to_crystal(data):
    """
    transform info from dictionary to Crystal  class
    """
    space_groupe = f_crystal.cl_crystal.SpaceGroupe()
    cell = f_crystal.cl_crystal.Cell()
    atom_site = f_crystal.cl_crystal.AtomSite()

    crystal = f_crystal.cl_crystal.Crystal(cell=cell, atom_site=atom_site, 
                      space_groupe=space_groupe)

    crystal.set_val(name=data["name"])
    
    l_relation = data_space_groupe_relation()
    from_dict_to_obj(data, l_relation, space_groupe)

    l_relation = data_cell_relation()
    from_dict_to_obj(data, l_relation, cell)
    
    cell.set_val(singony=space_groupe.get_val("singony"))

    l_key = data.keys()
    if (not ("loops" in l_key)):
        return crystal

    lab_xyz = "_atom_site_fract"
    lab_beta = "_atom_site_aniso_U"
    lab_chi = "_atom_site_susceptibility"
    llabs = [lab_xyz, lab_beta, lab_chi]
    lnumb = [[] for hh in llabs]

    for iloop, loop in enumerate(data["loops"]):
        l_key = loop.keys()
        for hh, numb in zip(llabs, lnumb):
            flag = any([True if hh1.startswith(hh) else False for hh1 in l_key])
            if flag:
                numb.append(iloop)
    
    l_atom_type = []
    if lnumb[0] != []:
        l_relation = data_atom_relation()
        for numb in lnumb[0]:
            l_key = list(data["loops"][numb].keys())
            n_atom = len(data["loops"][numb][l_key[0]])
            for i_atom in range(n_atom):
                dd = {}
                for key in l_key:
                    dd.update({key:data["loops"][numb][key][i_atom]})
                atom_type = f_crystal.cl_crystal.AtomType()
                from_dict_to_obj(dd, l_relation, atom_type)
                l_atom_type.append(atom_type)
    else:
        print("Fractional coordinates of atoms are not found")
        return crystal

    if lnumb[1] != []:
        l_relation = data_adp_relation()
        for numb in lnumb[1]:
            l_key = list(data["loops"][numb].keys())
            n_atom = len(data["loops"][numb][l_key[0]])
            for i_atom in range(n_atom):
                dd = {}
                for key in l_key:
                    dd.update({key:data["loops"][numb][key][i_atom]})
                for atom_type in l_atom_type:
                    if dd["_atom_site_aniso_label"] == atom_type.get_val("name"):
                        from_dict_to_obj(dd, l_relation, atom_type)
    
    if lnumb[2] != []:
        l_relation = data_chi_relation()
        for numb in lnumb[2]:
            l_key = list(data["loops"][numb].keys())
            n_atom = len(data["loops"][numb][l_key[0]])
            for i_atom in range(n_atom):
                dd = {}
                for key in l_key:
                    dd.update({key:data["loops"][numb][key][i_atom]})
                for atom_type in l_atom_type:
                    if dd["_atom_site_magnetism_aniso_label"] == atom_type.get_val("name"):
                        from_dict_to_obj(dd, l_relation, atom_type)
                        atom_type.set_val(flag_m=True)
            
    for atom_type in l_atom_type:
        atom_site.add_atom(atom_type)
    return crystal

def conv_data_to_experiment_single(data, f_dir, l_crystal):
    beam_polarization = f_experiment.f_powder_1d.cl_setup_powder_1d.BeamPolarization()
    setup = f_experiment.f_single.cl_setup_single.SetupSingle(beam_polarization=beam_polarization)
    observed_data = f_experiment.f_single.cl_observed_data_single.ObservedDataSingle(file_dir=f_dir)
    experiment = f_experiment.f_single.cl_experiment_single.ExperimentSingle(setup=setup, 
                                    observed_data=observed_data, file_dir=f_dir)

    l_relation = data_experiment_single_relation()
    from_dict_to_obj(data, l_relation, experiment)    

    l_relation = data_beam_polarization_single_relation()
    from_dict_to_obj(data, l_relation, beam_polarization)

    l_relation = data_observed_data_single_relation()
    from_dict_to_obj(data, l_relation, observed_data)
    
    observed_data.read_data()
    wave_length = observed_data.get_val("wave_length")
    field = observed_data.get_val("field")
    orientation = observed_data.get_val("orientation")

    setup.set_val(wave_length=wave_length)

    l_key = data.keys()
    if (not ("loops" in l_key)):
        return experiment 

    l_relation = data_extinction_single_relation()
    for loop in data["loops"]:
        l_key = list(loop.keys())
        flag_scale = "_sd_phase_extinction_radius" in l_key
        if flag_scale:
            n_crystal = len(loop["_sd_phase_name"])
            for i_crystal in range(n_crystal):
                dd = {}
                for key in l_key:
                    dd.update({key: loop[key][i_crystal]})
                for crystal in l_crystal:
                    if  dd["_sd_phase_name"] == crystal.get_val("name"):
                        extinction = crystal.get_val("extinction")
                        from_dict_to_obj(dd, l_relation, extinction)
                        name = crystal.get_val("name")
                        calculated_data = f_experiment.f_single.cl_calculated_data_single.CalculatedDataSingle(field=field,
                                orientation=orientation, name=name)
                        experiment.add_calculated_data(calculated_data)    
    return experiment

def conv_data_to_experiment_single_domain(data, f_dir, l_crystal):
    beam_polarization = f_experiment.f_powder_1d.cl_setup_powder_1d.BeamPolarization()
    setup = f_experiment.f_single.cl_setup_single.SetupSingle(beam_polarization=beam_polarization)
    observed_data = f_experiment.f_single_domain.cl_observed_data_single_domain.ObservedDataSingleDomain(file_dir=f_dir)
    experiment = f_experiment.f_single_domain.cl_experiment_single_domain.ExperimentSingleDomain(setup=setup, 
                                    observed_data=observed_data, file_dir=f_dir)

    l_relation = data_experiment_single_domain_relation()
    from_dict_to_obj(data, l_relation, experiment)    

    l_relation = data_beam_polarization_single_domain_relation()
    from_dict_to_obj(data, l_relation, beam_polarization)

    l_relation = data_observed_data_single_domain_relation()
    from_dict_to_obj(data, l_relation, observed_data)
    
    observed_data.read_data()
    wave_length = observed_data.get_val("wave_length")
    field = observed_data.get_val("field")
    #print(wave_length, field)
    setup.set_val(wave_length=wave_length)

    l_key = data.keys()
    if (not ("loops" in l_key)):
        return experiment 

    l_relation = data_extinction_single_domain_relation()
    for loop in data["loops"]:
        l_key = list(loop.keys())
        flag_scale = "_sdd_phase_extinction_radius" in l_key
        if flag_scale:
            n_crystal = len(loop["_sdd_phase_name"])
            for i_crystal in range(n_crystal):
                dd = {}
                for key in l_key:
                    dd.update({key: loop[key][i_crystal]})
                for crystal in l_crystal:
                    if  dd["_sdd_phase_name"] == crystal.get_val("name"):
                        extinction = crystal.get_val("extinction")
                        from_dict_to_obj(dd, l_relation, extinction)
                        name = crystal.get_val("name")
                        calculated_data = f_experiment.f_single.cl_calculated_data_single.CalculatedDataSingle(field=field,
                                name=name)
                        experiment.add_calculated_data(calculated_data) 
        flag_domain = "_sdd_scale_domain" in l_key   
        if flag_domain:
            l_relation = data_domain_single_domain_relation()
            dd = {}
            for relation in l_relation:
                lab_d = relation[0]
                val = [conv_str_to_text_float_logic(hh, "domain_{}".format(ihh+1)) for ihh, hh in enumerate(loop[lab_d])]
                dd.update({lab_d: val})
            from_dict_to_obj(dd, l_relation, experiment)

    #print(wave_length, field)
    return experiment

def conv_data_to_experiment_powder_1d(data, f_dir, l_crystal):
    beam_polarization = f_experiment.f_powder_1d.cl_setup_powder_1d.BeamPolarization()
    resolution = f_experiment.f_powder_1d.cl_setup_powder_1d.ResolutionPowder1D()
    factor_lorentz = f_experiment.f_powder_1d.cl_setup_powder_1d.FactorLorentzPowder1D()
    asymmetry = f_experiment.f_powder_1d.cl_setup_powder_1d.AsymmetryPowder1D()
    background = f_experiment.f_powder_1d.cl_setup_powder_1d.BackgroundPowder1D(file_dir=f_dir)
    setup = f_experiment.f_powder_1d.cl_setup_powder_1d.SetupPowder1D(resolution=resolution, 
            factor_lorentz=factor_lorentz, asymmetry=asymmetry, 
            beam_polarization=beam_polarization, background=background)
        
    observed_data = f_experiment.f_powder_1d.cl_observed_data_powder_1d.ObservedDataPowder1D(file_dir=f_dir)
    
    experiment = f_experiment.f_powder_1d.cl_experiment_powder_1d.ExperimentPowder1D(setup=setup, 
                                    observed_data=observed_data, file_dir=f_dir)

    l_relation = data_experiment_powder_1d_relation()
    from_dict_to_obj(data, l_relation, experiment)    
    
    l_relation = data_resolution_powder_1d_relation()
    from_dict_to_obj(data, l_relation, resolution)    

    l_relation = data_asymmetry_powder_1d_relation()
    from_dict_to_obj(data, l_relation, asymmetry)    

    l_relation = data_zero_shift_powder_1d_relation()
    from_dict_to_obj(data, l_relation, setup)    

    l_relation = data_beam_polarization_powder_1d_relation()
    from_dict_to_obj(data, l_relation, beam_polarization)

    l_relation = data_observed_data_powder_1d_relation()
    from_dict_to_obj(data, l_relation, observed_data)

    l_relation = data_background_powder_1d_relation()
    from_dict_to_obj(data, l_relation, background)
    
    background.read_data()

    observed_data.read_data()
    wave_length = observed_data.get_val("wave_length")
    field = observed_data.get_val("field")

    setup.set_val(wave_length=wave_length)

    l_key = data.keys()
    if (not ("loops" in l_key)):
        return experiment 

    for loop in data["loops"]:
        l_key = list(loop.keys())
        flag_scale = "_pd_phase_scale" in l_key
        flag_exclude = "_pd_exclude_tth_min" in l_key
        if flag_scale:
            n_crystal = len(loop["_pd_phase_name"])
            for i_crystal in range(n_crystal):
                dd = {}
                for key in l_key:
                    dd.update({key: loop[key][i_crystal]})
                for crystal in l_crystal:
                    if  dd["_pd_phase_name"] == crystal.get_val("name"):
                        scale = conv_str_to_text_float_logic(dd["_pd_phase_scale"], "scale")
                        i_g = conv_str_to_text_float_logic(dd["_pd_phase_igsize"], "i_g")
                        crystal.set_val(i_g=i_g)
                        name = crystal.get_val("name")
                        calculated_data = f_experiment.f_powder_1d.cl_calculated_data_powder_1d.CalculatedDataPowder1D(field=field,
                                scale=scale, name=name)
                        experiment.add_calculated_data(calculated_data)  
        if flag_exclude:
            l_excl_tth_min = [conv_str_to_text_float_logic(hh, "tth_min") for hh in loop["_pd_exclude_tth_min"]]
            l_excl_tth_max = [conv_str_to_text_float_logic(hh, "tth_max") for hh in loop["_pd_exclude_tth_max"]]
            experiment.set_val(excl_tth_min=l_excl_tth_min, excl_tth_max=l_excl_tth_max)
    return experiment

def conv_data_to_experiment_powder_2d(data, f_dir, l_crystal):
    beam_polarization = f_experiment.f_powder_1d.cl_setup_powder_1d.BeamPolarization()
    resolution = f_experiment.f_powder_2d.cl_setup_powder_2d.ResolutionPowder2D()
    factor_lorentz = f_experiment.f_powder_2d.cl_setup_powder_2d.FactorLorentzPowder2D()
    asymmetry = f_experiment.f_powder_2d.cl_setup_powder_2d.AsymmetryPowder2D()
    background = f_experiment.f_powder_2d.cl_setup_powder_2d.BackgroundPowder2D(file_dir=f_dir)
    setup = f_experiment.f_powder_2d.cl_setup_powder_2d.SetupPowder2D(resolution=resolution, 
            factor_lorentz=factor_lorentz, asymmetry=asymmetry, 
            beam_polarization=beam_polarization, background=background)
        
    observed_data = f_experiment.f_powder_2d.cl_observed_data_powder_2d.ObservedDataPowder2D(file_dir=f_dir)
    
    experiment = f_experiment.f_powder_2d.cl_experiment_powder_2d.ExperimentPowder2D(setup=setup, 
                                    observed_data=observed_data, file_dir=f_dir)

    l_relation = data_experiment_powder_2d_relation()
    from_dict_to_obj(data, l_relation, experiment)    

    
    l_relation = data_resolution_powder_2d_relation()
    from_dict_to_obj(data, l_relation, resolution)    

    l_relation = data_asymmetry_powder_2d_relation()
    from_dict_to_obj(data, l_relation, asymmetry)    

    l_relation = data_zero_shift_powder_2d_relation()
    from_dict_to_obj(data, l_relation, setup)    

    l_relation = data_beam_polarization_powder_2d_relation()
    from_dict_to_obj(data, l_relation, beam_polarization)

    l_relation = data_observed_data_powder_2d_relation()
    from_dict_to_obj(data, l_relation, observed_data)

    l_relation = data_background_powder_2d_relation()
    from_dict_to_obj(data, l_relation, background)
    
    background.read_data()

    observed_data.read_data()
    wave_length = observed_data.get_val("wave_length")
    field = observed_data.get_val("field")

    setup.set_val(wave_length=wave_length)

    l_key = data.keys()
    if (not ("loops" in l_key)):
        return experiment 
    
    for loop in data["loops"]:
        l_key = list(loop.keys())
        flag_scale = "_2dpd_phase_scale" in l_key
        if flag_scale:
            n_crystal = len(loop["_2dpd_phase_name"])
            for i_crystal in range(n_crystal):
                dd = {}
                for key in l_key:
                    dd.update({key: loop[key][i_crystal]})
                for crystal in l_crystal:
                    if  dd["_2dpd_phase_name"] == crystal.get_val("name"):
                        scale = conv_str_to_text_float_logic(dd["_2dpd_phase_scale"], "scale")
                        i_g = conv_str_to_text_float_logic(dd["_2dpd_phase_igsize"], "i_g")
                        crystal.set_val(i_g=i_g)
                        name = crystal.get_val("name")
                        calculated_data = f_experiment.f_powder_2d.cl_calculated_data_powder_2d.CalculatedDataPowder2D(field=field,
                                scale=scale, name=name)
                        experiment.add_calculated_data(calculated_data)   
    return experiment


def conv_data_to_experiment_powder_texture_2d(data, f_dir, l_crystal):
    beam_polarization = f_experiment.f_powder_1d.cl_setup_powder_1d.BeamPolarization()
    resolution = f_experiment.f_powder_2d.cl_setup_powder_2d.ResolutionPowder2D()
    factor_lorentz = f_experiment.f_powder_2d.cl_setup_powder_2d.FactorLorentzPowder2D()
    asymmetry = f_experiment.f_powder_2d.cl_setup_powder_2d.AsymmetryPowder2D()
    background = f_experiment.f_powder_2d.cl_setup_powder_2d.BackgroundPowder2D(file_dir=f_dir)
    setup = f_experiment.f_powder_2d.cl_setup_powder_2d.SetupPowder2D(resolution=resolution, 
            factor_lorentz=factor_lorentz, asymmetry=asymmetry, 
            beam_polarization=beam_polarization, background=background)
        
    observed_data = f_experiment.f_powder_2d.cl_observed_data_powder_2d.ObservedDataPowder2D(file_dir=f_dir)
    
    experiment = f_experiment.f_powder_texture_2d.cl_experiment_powder_texture_2d.ExperimentPowderTexture2D(setup=setup, 
                                    observed_data=observed_data, file_dir=f_dir)

    l_relation = data_experiment_powder_texture_2d_relation()
    from_dict_to_obj(data, l_relation, experiment)    

    
    l_relation = data_resolution_powder_texture_2d_relation()
    from_dict_to_obj(data, l_relation, resolution)    

    l_relation = data_asymmetry_powder_texture_2d_relation()
    from_dict_to_obj(data, l_relation, asymmetry)    

    l_relation = data_zero_shift_powder_texture_2d_relation()
    from_dict_to_obj(data, l_relation, setup)    

    l_relation = data_beam_polarization_powder_texture_2d_relation()
    from_dict_to_obj(data, l_relation, beam_polarization)

    l_relation = data_observed_data_powder_texture_2d_relation()
    from_dict_to_obj(data, l_relation, observed_data)

    l_relation = data_background_powder_texture_2d_relation()
    from_dict_to_obj(data, l_relation, background)
    
    background.read_data()

    observed_data.read_data()
    wave_length = observed_data.get_val("wave_length")
    field = observed_data.get_val("field")

    setup.set_val(wave_length=wave_length)

    l_key = data.keys()
    if (not ("loops" in l_key)):
        return experiment 
    
    for loop in data["loops"]:
        l_key = list(loop.keys())
        flag_scale = "_2dpdt_phase_scale" in l_key
        if flag_scale:
            n_crystal = len(loop["_2dpdt_phase_name"])
            for i_crystal in range(n_crystal):
                dd = {}
                for key in l_key:
                    dd.update({key: loop[key][i_crystal]})
                for crystal in l_crystal:
                    if  dd["_2dpdt_phase_name"] == crystal.get_val("name"):
                        scale = conv_str_to_text_float_logic(dd["_2dpdt_phase_scale"], "scale")
                        i_g = conv_str_to_text_float_logic(dd["_2dpdt_phase_igsize"], "i_g")
                        crystal.set_val(i_g=i_g)
                        name = crystal.get_val("name")
                        calculated_data = f_experiment.f_powder_2d.cl_calculated_data_powder_2d.CalculatedDataPowder2D(field=field,
                                scale=scale, name=name)
                        experiment.add_calculated_data(calculated_data)   
    return experiment


def conv_model_to_rcif(model):
    l_data = []
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
            if isinstance(val, f_common.cl_variable.Variable):
                sval = "{:}".format(val.print_with_sigma())
            else:
                sval = "{:}".format(val)
            ddata[key_d] = sval 
        else:
            print(50*"ERROR ")
            print("type_val: ", type_val)
        return
        
    
    for obj in model._list_crystal:
        dd = {"name": obj.get_val("name")}
        l_relation = data_space_groupe_relation()
        space_groupe = obj.get_val("space_groupe")

        for relation in l_relation:
            temp_func(dd, relation, space_groupe)

        l_relation = data_cell_relation()
        cell = obj.get_val("cell")

        for relation in l_relation:
            temp_func(dd, relation, cell)
            
        lloop_d = []
        
        d_at = {}
        l_relation = data_atom_relation()

        atom_site = obj.get_val("atom_site")
        l_atom_type = atom_site._list_atom_type

        l_lab_d_used = []
        l_dhelp = []
        for obj_2 in l_atom_type:
            dhelp = {}
            for relation in l_relation:
                lab_d = relation[0]
                temp_func(dhelp, relation, obj_2)
                if (not(lab_d in l_lab_d_used)): 
                    l_lab_d_used.append(lab_d)
            l_dhelp.append(dhelp) 
        for lab_d in l_lab_d_used:
            d_at[lab_d] = [hh[lab_d] for hh in l_dhelp]
        lloop_d.append(d_at)
        
        d_bscat = {}
        l_relation = data_bscat_relation()

        l_lab_d_used = []
        l_dhelp = []
        for obj_2 in l_atom_type:
            dhelp = {}
            for relation in l_relation:
                lab_d = relation[0]
                temp_func(dhelp, relation, obj_2)
                if (not(lab_d in l_lab_d_used)):
                    l_lab_d_used.append(lab_d)
            l_dhelp.append(dhelp) 
        lab_uniq = "_atom_site_type_symbol"
        l_val = [hh[lab_uniq] for hh in l_dhelp]
        s_uniq = set(l_val)
        l_ind_uniq = [l_val.index(hh) for hh in s_uniq]
        for lab_d in l_lab_d_used:
            d_bscat[lab_d] = [l_dhelp[ind][lab_d] for ind in l_ind_uniq]
        lloop_d.append(d_bscat)

        d_adp = {}
        l_relation = data_adp_relation()
        l_lab_d_used = []
        l_dhelp = []
        for obj_2 in l_atom_type:
            dhelp = {}
            flag_beta = obj_2.get_val("adp_type") == "uani"
            if flag_beta:
                for relation in l_relation:
                    lab_d = relation[0]
                    temp_func(dhelp, relation, obj_2)
                    if (not(lab_d in l_lab_d_used)): 
                        l_lab_d_used.append(lab_d)
                l_dhelp.append(dhelp) 
        for lab_d in l_lab_d_used:
            d_adp[lab_d] = [hh[lab_d] for hh in l_dhelp]
        lloop_d.append(d_adp)
        
        d_chi = {}
        l_relation = data_chi_relation()
        l_lab_d_used = []
        l_dhelp = []
        for obj_2 in l_atom_type:
            dhelp = {}
            flag_m = obj_2.get_val("flag_m")
            if flag_m:
                for relation in l_relation:
                    lab_d = relation[0]
                    temp_func(dhelp, relation, obj_2)
                    if (not(lab_d in l_lab_d_used)): 
                        l_lab_d_used.append(lab_d)
                l_dhelp.append(dhelp)
        for lab_d in l_lab_d_used:
            d_chi[lab_d] = [hh[lab_d] for hh in l_dhelp]
        lloop_d.append(d_chi)
        
        if lloop_d != []:
            dd["loops"] = lloop_d
            
        l_data.append(dd)

    for obj in model._list_experiment:
        dd = {"name": obj.get_val("name")}

        if isinstance(obj, f_experiment.f_single.cl_experiment_single.ExperimentSingle):
            l_relation = data_experiment_single_relation()
            for relation in l_relation:
                temp_func(dd, relation, obj)
            
            observed_data = obj.get_val("observed_data")
            l_relation = data_observed_data_single_relation()
            for relation in l_relation:
                temp_func(dd, relation, observed_data)
            
            setup = obj.get_val("setup")
            beam_polarization = setup.get_val("beam_polarization")

            l_relation = data_beam_polarization_single_relation()
            for relation in l_relation:
                temp_func(dd, relation, beam_polarization)

            
        if isinstance(obj, f_experiment.f_single_domain.cl_experiment_single_domain.ExperimentSingleDomain):
            l_relation = data_experiment_single_domain_relation()
            for relation in l_relation:
                temp_func(dd, relation, obj)
            
            observed_data = obj.get_val("observed_data")
            l_relation = data_observed_data_single_domain_relation()
            for relation in l_relation:
                temp_func(dd, relation, observed_data)
            
            setup = obj.get_val("setup")
            beam_polarization = setup.get_val("beam_polarization")

            l_relation = data_beam_polarization_single_domain_relation()
            for relation in l_relation:
                temp_func(dd, relation, beam_polarization)

        if isinstance(obj, f_experiment.f_powder_1d.cl_experiment_powder_1d.ExperimentPowder1D):
            l_relation = data_experiment_powder_1d_relation()
            for relation in l_relation:
                temp_func(dd, relation, obj)

            observed_data = obj.get_val("observed_data")
            l_relation = data_observed_data_powder_1d_relation()
            for relation in l_relation:
                temp_func(dd, relation, observed_data)

            setup = obj.get_val("setup")
            l_relation = data_zero_shift_powder_1d_relation()
            for relation in l_relation:
                temp_func(dd, relation, setup)

            background = setup.get_val("background")
            l_relation = data_background_powder_1d_relation()
            for relation in l_relation:
                temp_func(dd, relation, background)
            background.save_data()

            beam_polarization = setup.get_val("beam_polarization")
            l_relation = data_beam_polarization_powder_1d_relation()
            for relation in l_relation:
                temp_func(dd, relation, beam_polarization)
            
            resolution = setup.get_val("resolution")
            l_relation = data_resolution_powder_1d_relation()
            for relation in l_relation:
                temp_func(dd, relation, resolution)
            
            asymmetry = setup.get_val("asymmetry")
            l_relation = data_asymmetry_powder_1d_relation()
            for relation in l_relation:
                temp_func(dd, relation, asymmetry)

        if isinstance(obj, f_experiment.f_powder_2d.cl_experiment_powder_2d.ExperimentPowder2D):
            l_relation = data_experiment_powder_2d_relation()
            for relation in l_relation:
                temp_func(dd, relation, obj)

            observed_data = obj.get_val("observed_data")
            l_relation = data_observed_data_powder_2d_relation()
            for relation in l_relation:
                temp_func(dd, relation, observed_data)

            setup = obj.get_val("setup")
            l_relation = data_zero_shift_powder_2d_relation()
            for relation in l_relation:
                temp_func(dd, relation, setup)

            background = setup.get_val("background")
            l_relation = data_background_powder_2d_relation()
            for relation in l_relation:
                temp_func(dd, relation, background)
            background.save_data()
            
            beam_polarization = setup.get_val("beam_polarization")
            l_relation = data_beam_polarization_powder_2d_relation()
            for relation in l_relation:
                temp_func(dd, relation, beam_polarization)
            
            resolution = setup.get_val("resolution")
            l_relation = data_resolution_powder_2d_relation()
            for relation in l_relation:
                temp_func(dd, relation, resolution)
            
            asymmetry = setup.get_val("asymmetry")
            l_relation = data_asymmetry_powder_2d_relation()
            for relation in l_relation:
                temp_func(dd, relation, asymmetry)


        if isinstance(obj, f_experiment.f_powder_texture_2d.cl_experiment_powder_texture_2d.ExperimentPowderTexture2D):
            l_relation = data_experiment_powder_texture_2d_relation()
            for relation in l_relation:
                temp_func(dd, relation, obj)

            observed_data = obj.get_val("observed_data")
            l_relation = data_observed_data_powder_texture_2d_relation()
            for relation in l_relation:
                temp_func(dd, relation, observed_data)

            setup = obj.get_val("setup")
            l_relation = data_zero_shift_powder_texture_2d_relation()
            for relation in l_relation:
                temp_func(dd, relation, setup)

            background = setup.get_val("background")
            l_relation = data_background_powder_texture_2d_relation()
            for relation in l_relation:
                temp_func(dd, relation, background)
            background.save_data()
            
            beam_polarization = setup.get_val("beam_polarization")
            l_relation = data_beam_polarization_powder_texture_2d_relation()
            for relation in l_relation:
                temp_func(dd, relation, beam_polarization)
            
            resolution = setup.get_val("resolution")
            l_relation = data_resolution_powder_texture_2d_relation()
            for relation in l_relation:
                temp_func(dd, relation, resolution)
            
            asymmetry = setup.get_val("asymmetry")
            l_relation = data_asymmetry_powder_texture_2d_relation()
            for relation in l_relation:
                temp_func(dd, relation, asymmetry)

        d_eph = {}
        lloop_d = []
        l_calculated_data = obj._list_calculated_data

        if isinstance(obj, f_experiment.f_single.cl_experiment_single.ExperimentSingle):
            crystal = model._list_crystal[0]
            extinction = crystal.get_val("extinction")
            l_relation = data_extinction_single_relation()
            
            l_dhelp = []
            for calculated_data in l_calculated_data:
                name = calculated_data.get_val("name")
                dhelp = {"_sd_phase_name": name}
                for relation in l_relation:
                    temp_func(dhelp, relation, extinction)
                l_dhelp.append(dhelp) 
            llab_dhelp = list(dhelp.keys())
            for lab_d in llab_dhelp:
                d_eph[lab_d] = [hh[lab_d] for hh in l_dhelp]
            lloop_d.append(d_eph)

        if isinstance(obj, f_experiment.f_single_domain.cl_experiment_single_domain.ExperimentSingleDomain):
            l_val = ["{:}".format(hh) for hh in obj.get_val("scale_domain")]
            d_domain = {"_sdd_scale_domain": l_val}
            lloop_d.append(d_domain)

            crystal = model._list_crystal[0]
            extinction = crystal.get_val("extinction")
            l_relation = data_extinction_single_domain_relation()
            
            l_dhelp = []
            for calculated_data in l_calculated_data:
                name = calculated_data.get_val("name")
                dhelp = {"_sdd_phase_name": name}
                for relation in l_relation:
                    temp_func(dhelp, relation, extinction)
                l_dhelp.append(dhelp) 
            llab_dhelp = list(dhelp.keys())
            for lab_d in llab_dhelp:
                d_eph[lab_d] = [hh[lab_d] for hh in l_dhelp]
            lloop_d.append(d_eph)

        if isinstance(obj, f_experiment.f_powder_1d.cl_experiment_powder_1d.ExperimentPowder1D):
            l_relation = [("_pd_phase_name", "name", "text"),
                          ("_pd_phase_scale", "scale", "val")]

            l_dhelp = []
            for calculated_data in l_calculated_data:
                dhelp = {}
                for relation in l_relation:
                    temp_func(dhelp, relation, calculated_data)
                name = calculated_data.get_val("name")
                ind = 0
                for i_crystal, crystal in enumerate(model._list_crystal):
                    if crystal.get_val("name") == name:
                        ind = i_crystal
                        break
                crystal = model._list_crystal[ind]
                
                relation = ("_pd_phase_igsize", "i_g", "val")
                temp_func(dhelp, relation, crystal)
                l_dhelp.append(dhelp) 
            l_lab_d = list(dhelp.keys())

            for lab_d in l_lab_d:
                d_eph[lab_d] = [hh[lab_d] for hh in l_dhelp]
            lloop_d.append(d_eph)

            l_excl_tth_min = ["{:}".format(hh) for hh in obj._p_excl_tth_min]
            l_excl_tth_max = ["{:}".format(hh) for hh in obj._p_excl_tth_max]
            if len(l_excl_tth_min) != 0:
                d_eph_2 = {"_pd_exclude_tth_min":l_excl_tth_min, "_pd_exclude_tth_max":l_excl_tth_max}
                lloop_d.append(d_eph_2)

        if isinstance(obj, f_experiment.f_powder_2d.cl_experiment_powder_2d.ExperimentPowder2D):
            l_relation = [("_2dpd_phase_name", "name", "text"),
                          ("_2dpd_phase_scale", "scale", "val")]

            l_dhelp = []
            for calculated_data in l_calculated_data:
                dhelp = {}
                for relation in l_relation:
                    temp_func(dhelp, relation, calculated_data)
                name = calculated_data.get_val("name")
                ind = 0
                for i_crystal, crystal in enumerate(model._list_crystal):
                    if crystal.get_val("name") == name:
                        ind = i_crystal
                        break
                crystal = model._list_crystal[ind]
                relation = ("_2dpd_phase_igsize", "i_g", "val")
                temp_func(dhelp, relation, crystal)
                l_dhelp.append(dhelp) 

            l_lab_d = list(dhelp.keys())
            for lab_d in l_lab_d:
                d_eph[lab_d] = [hh[lab_d] for hh in l_dhelp]
            lloop_d.append(d_eph)
            

        if isinstance(obj, f_experiment.f_powder_texture_2d.cl_experiment_powder_texture_2d.ExperimentPowderTexture2D):
            l_relation = [("_2dpdt_phase_name", "name", "text"),
                          ("_2dpdt_phase_scale", "scale", "val")]

            l_dhelp = []
            for calculated_data in l_calculated_data:
                dhelp = {}
                for relation in l_relation:
                    temp_func(dhelp, relation, calculated_data)
                name = calculated_data.get_val("name")
                ind = 0
                for i_crystal, crystal in enumerate(model._list_crystal):
                    if crystal.get_val("name") == name:
                        ind = i_crystal
                        break
                crystal = model._list_crystal[ind]
                relation = ("_2dpdt_phase_igsize", "i_g", "val")
                temp_func(dhelp, relation, crystal)
                l_dhelp.append(dhelp) 

            l_lab_d = list(dhelp.keys())
            for lab_d in l_lab_d:
                d_eph[lab_d] = [hh[lab_d] for hh in l_dhelp]
            lloop_d.append(d_eph)
            
        if lloop_d != []:
            dd["loops"] = lloop_d
            
        l_data.append(dd)

    rcif = f_rcif.cl_rcif.RCif()
    rcif.glob["data"] = l_data
    rcif.glob["name"] = model.get_val("name")

    file_out = model.get_val("file_out")
    if file_out is None:
        file_out = "full.lis"
    rcif.glob["_file_name_output_listing"] = os.path.basename(file_out)
    return rcif

def conv_crystal_to_data(crystal):
    data = {}
    return data

def conv_experiment_single_to_data(experiment):
    data = {}
    return data

def conv_experiment_single_domain_to_data(experiment):
    data = {}
    return data

def conv_experiment_powder_1d_to_data(experiment):
    data = {}
    return data

def conv_experiment_powder_2d_to_data(experiment):
    data = {}
    return data



def from_dict_to_obj(dict_i, l_relation, obj):
    """
    l_relation is list of (lab_d, lab_o, val_type)
    """
    
    l_key = list(dict_i.keys())
    l_key_lower = [hh.lower() for hh in l_key]
    l_numb = [ihh for ihh, hh in enumerate(l_relation) if (hh[0].lower() in l_key_lower)]
    d_args = {}
    for numb in l_numb:
        lab_d, lab_o, val_type = l_relation[numb]
        key_d = l_key[l_key_lower.index(lab_d.lower())]
        if val_type == "val":
            #val = [dict_i[llab_d[numb]], False, ""]
            val = dict_i[key_d]
            val = conv_str_to_text_float_logic(val, lab_o)
        elif val_type == "text":
            val = dict_i[key_d]
        elif val_type == "list":
            val = dict_i[key_d]
        elif val_type == "logic":
            val = dict_i[key_d]
            val = conv_str_to_text_float_logic(val, lab_o)
        else:
            print("mistake in type variable of 'from_dict_to_obj' function")
            val = None
        d_args.update({lab_o: val})
    obj.set_val(**d_args)
    return


def data_space_groupe_relation():
    l_relation = [ ("_space_group_name_H-M_alt", "spgr_given_name", "text"),
        ("_space_group_it_coordinate_system_code", "spgr_choice", "text")]
    return l_relation

def data_cell_relation():
    l_relation = [("_cell_length_a", "a", "val"), ("_cell_length_b", "b", "val"),
        ("_cell_length_c", "c", "val"), ("_cell_angle_alpha", "alpha", "val"),
        ("_cell_angle_beta", "beta", "val"), ("_cell_angle_gamma", "gamma", "val")]
    return l_relation

def data_atom_relation():
    l_relation = [("_atom_site_label", "name", "text"),
        ("_atom_site_type_symbol", "type_n", "text"),
        ("_atom_site_fract_x", "x", "val"),
        ("_atom_site_fract_y", "y", "val"),
        ("_atom_site_fract_z", "z", "val"),
        ("_atom_site_occupancy", "occupation", "val"),
        ("_atom_site_b_iso_or_equiv", "b_iso", "val"),
        ("_atom_site_adp_type", "adp_type", "text")]
    return l_relation

def data_bscat_relation():
    l_relation = [("_atom_site_type_symbol", "type_n", "text"),
        ("_atom_site_bscat", "b_scat", "val")]
    return l_relation

def data_adp_relation():
    l_relation = [("_atom_site_aniso_label", "name", "text"),
        ("_atom_site_aniso_U_11", "u_11", "val"),
        ("_atom_site_aniso_U_12", "u_12", "val"),
        ("_atom_site_aniso_U_13", "u_13", "val"),
        ("_atom_site_aniso_U_22", "u_22", "val"),
        ("_atom_site_aniso_U_23", "u_23", "val"),
        ("_atom_site_aniso_U_33", "u_33", "val")]
    return l_relation


def data_chi_relation():
    l_relation = [("_atom_site_magnetism_aniso_label", "name", "text"),
        ("_atom_site_magnetism_type_symbol", "type_m", "text"),
        ("_atom_site_magnetism_kappa", "kappa", "val"),
        ("_atom_site_magnetism_factor_lande", "factor_lande", "val"),
        ("_atom_site_magnetism_type", "chi_type", "text"),
        ("_atom_site_susceptibility_aniso_chi_11", "chi_11", "val"),
        ("_atom_site_susceptibility_aniso_chi_12", "chi_12", "val"),
        ("_atom_site_susceptibility_aniso_chi_13", "chi_13", "val"),
        ("_atom_site_susceptibility_aniso_chi_22", "chi_22", "val"),
        ("_atom_site_susceptibility_aniso_chi_23", "chi_23", "val"),
        ("_atom_site_susceptibility_aniso_chi_33", "chi_33", "val")]
    return l_relation



def data_experiment_powder_1d_relation():
    l_relation = [("name", "name", "text"),
        ("_pd_file_name_output", "file_out", "text"),
        ("_pd_chi2_sum", "flag_chi2_sum", "logic"),
        ("_pd_chi2_diff", "flag_chi2_diff", "logic"),
        ("_pd_chi2_up", "flag_chi2_up", "logic"),
        ("_pd_chi2_down", "flag_chi2_down", "logic")]
    return l_relation


def data_background_powder_1d_relation():
    l_relation = [("_pd_file_name_bkgr", "file_name", "text")]
    return l_relation

def data_observed_data_powder_1d_relation():
    l_relation = [("_pd_file_name_input", "file_name", "text"),
        ("_pd_tth_min", "tth_min" , "val"),
        ("_pd_tth_max", "tth_max" , "val")]
    return l_relation

def data_beam_polarization_powder_1d_relation():
    l_relation = [("_pd_beam_polarization_up", "p_u", "val"),
        ("_pd_beam_flipper_efficiency", "flipper_efficiency", "val")]
    return l_relation

def data_resolution_powder_1d_relation():
    l_relation = [("_pd_resolution_u", "u", "val"),
        ("_pd_resolution_v", "v", "val"),
        ("_pd_resolution_w", "w", "val"),
        ("_pd_resolution_x", "x", "val"),
        ("_pd_resolution_y", "y", "val")]
    return l_relation

def data_asymmetry_powder_1d_relation():
    l_relation = [("_pd_reflex_asymmetry_p1", "p1", "val"),
        ("_pd_reflex_asymmetry_p2", "p2", "val"),
        ("_pd_reflex_asymmetry_p3", "p3", "val"),
        ("_pd_reflex_asymmetry_p4", "p4", "val")]
    return l_relation

def data_zero_shift_powder_1d_relation():
    l_relation = [("_pd_shift_const", "zero_shift", "val")]
    return l_relation

def data_exclude_powder_1d_relation():
    l_relation = [("_pd_exclude_tth_min", "excl_tth_min", "val"),
        ("_pd_exclude_tth_max", "excl_tth_max", "val")]
    return l_relation

def data_experiment_powder_2d_relation():
    l_relation = [("name", "name", "text"),
        ("_2dpd_file_name_output", "file_out", "text"),
        ("_2dpd_chi2_sum", "flag_chi2_sum", "logic"),
        ("_2dpd_chi2_diff", "flag_chi2_diff", "logic"),
        ("_2dpd_chi2_up", "flag_chi2_up", "logic"),
        ("_2dpd_chi2_down", "flag_chi2_down", "logic")]
    return l_relation

def data_background_powder_2d_relation():
    l_relation =[("_2dpd_file_name_bkgr", "file_name", "text")]
    return l_relation

def data_observed_data_powder_2d_relation():
    l_relation = [("_2dpd_file_name_input", "file_name", "text"),
        ("_2dpd_tth_min", "tth_min" , "val"),
        ("_2dpd_tth_max", "tth_max" , "val"),
        ("_2dpd_phi_min", "phi_min" , "val"),
        ("_2dpd_phi_max", "phi_max" , "val")]
    return l_relation

def data_beam_polarization_powder_2d_relation():
    l_relation = [("_2dpd_beam_polarization_up", "p_u", "val"),
        ("_2dpd_beam_flipper_efficiency", "flipper_efficiency", "val")]
    return l_relation

def data_resolution_powder_2d_relation():
    l_relation = [("_2dpd_resolution_u", "u", "val"),
        ("_2dpd_resolution_v", "v", "val"),
        ("_2dpd_resolution_w", "w", "val"),
        ("_2dpd_resolution_x", "x", "val"),
        ("_2dpd_resolution_y", "y", "val")]
    return l_relation

def data_asymmetry_powder_2d_relation():
    l_relation = [("_2dpd_reflex_asymmetry_p1", "p1", "val"),
        ("_2dpd_reflex_asymmetry_p2", "p2", "val"),
        ("_2dpd_reflex_asymmetry_p3", "p3", "val"),
        ("_2dpd_reflex_asymmetry_p4", "p4", "val")]
    return l_relation

def data_zero_shift_powder_2d_relation():
    l_relation =[("_2dpd_shift_const", "zero_shift", "val")]
    return l_relation





def data_experiment_powder_texture_2d_relation():
    l_relation = [("name", "name", "text"),
        ("_2dpdt_file_name_output", "file_out", "text"),
        ("_2dpdt_chi2_sum", "flag_chi2_sum", "logic"),
        ("_2dpdt_chi2_diff", "flag_chi2_diff", "logic"),
        ("_2dpdt_chi2_up", "flag_chi2_up", "logic"),
        ("_2dpdt_chi2_down", "flag_chi2_down", "logic"),
        ("_2dpdt_h_ax", "h_ax", "val"),
        ("_2dpdt_k_ax", "k_ax", "val"),
        ("_2dpdt_l_ax", "l_ax", "val"),
        ("_2dpdt_g_1", "g_1", "val"),
        ("_2dpdt_g_2", "g_2", "val"),
        ("_2dpdt_phi_0", "phi_0", "val")]
    return l_relation

def data_background_powder_texture_2d_relation():
    l_relation =[("_2dpdt_file_name_bkgr", "file_name", "text")]
    return l_relation

def data_observed_data_powder_texture_2d_relation():
    l_relation = [("_2dpdt_file_name_input", "file_name", "text"),
        ("_2dpdt_tth_min", "tth_min" , "val"),
        ("_2dpdt_tth_max", "tth_max" , "val"),
        ("_2dpdt_phi_min", "phi_min" , "val"),
        ("_2dpdt_phi_max", "phi_max" , "val")]
    return l_relation

def data_beam_polarization_powder_texture_2d_relation():
    l_relation = [("_2dpdt_beam_polarization_up", "p_u", "val"),
        ("_2dpdt_beam_flipper_efficiency", "flipper_efficiency", "val")]
    return l_relation

def data_resolution_powder_texture_2d_relation():
    l_relation = [("_2dpdt_resolution_u", "u", "val"),
        ("_2dpdt_resolution_v", "v", "val"),
        ("_2dpdt_resolution_w", "w", "val"),
        ("_2dpdt_resolution_x", "x", "val"),
        ("_2dpdt_resolution_y", "y", "val")]
    return l_relation

def data_asymmetry_powder_texture_2d_relation():
    l_relation = [("_2dpdt_reflex_asymmetry_p1", "p1", "val"),
        ("_2dpdt_reflex_asymmetry_p2", "p2", "val"),
        ("_2dpdt_reflex_asymmetry_p3", "p3", "val"),
        ("_2dpdt_reflex_asymmetry_p4", "p4", "val")]
    return l_relation

def data_zero_shift_powder_texture_2d_relation():
    l_relation =[("_2dpdt_shift_const", "zero_shift", "val")]
    return l_relation




def data_experiment_single_relation():
    l_relation = [("name", "name", "text"),
        ("_sd_file_name_output", "file_out", "text")]
    return l_relation

def data_observed_data_single_relation():
    l_relation = [("_sd_file_name_input", "file_name", "text")]
    return l_relation

def data_beam_polarization_single_relation():
    l_relation = [("_sd_beam_polarization_up", "p_u", "val"),
        ("_sd_beam_flipper_efficiency", "flipper_efficiency", "val")]
    return l_relation

def data_extinction_single_relation():
    l_relation = [("_sd_phase_extinction_radius", "domain_radius", "val"),
        ("_sd_phase_extinction_mosaicity", "mosaicity", "val")]
    return l_relation



def data_experiment_single_domain_relation():
    l_relation = [("name", "name", "text"),
        ("_sdd_file_name_output", "file_out", "text")]
    return l_relation

def data_observed_data_single_domain_relation():
    l_relation = [("_sdd_file_name_input", "file_name", "text")]
    return l_relation

def data_beam_polarization_single_domain_relation():
    l_relation = [("_sdd_beam_polarization_up", "p_u", "val"),
        ("_sdd_beam_flipper_efficiency", "flipper_efficiency", "val")]
    return l_relation

def data_extinction_single_domain_relation():
    l_relation = [("_sdd_phase_extinction_radius", "domain_radius", "val"),
        ("_sdd_phase_extinction_mosaicity", "mosaicity", "val")]
    return l_relation

def data_domain_single_domain_relation():
    l_relation = [("_sdd_scale_domain", "scale_domain", "list")]
    return l_relation


def conv_str_to_text_float_logic(sval, name=""):
    if (not (isinstance(sval,str))):
        return sval
    try:
        if len(sval.strip().split("("))>1:
            l_help = sval.split("(")
            val_1 = float(l_help[0])
            val = f_common.cl_variable.Variable(val_1, True, name)
            pass
        else:
            val = float(sval)
        return val
    except:
        val_1 = sval.strip()
        if len(val_1) == 0:
            val = None
        if (((val_1[0]=="\"")|(val_1[0]=="'"))&((val_1[-1]=="\"")|(val_1[-1]=="'"))):
            val_1=val_1[1:-1]
        if len(val_1) == 0:
            val = None
        elif val_1.lower() == "true":
            val = True
        elif val_1.lower() == "false":
            val = False
        else:
            val = val_1
    return val



def main(larg):
    if len(larg) > 1:
        fname = "powder.rcif"
        fname = larg[1]
    else:
        fname = input("What is the name of the 'rcif' file?\n")
        
    
    

if __name__ == "__main__":
    main(sys.argv)
