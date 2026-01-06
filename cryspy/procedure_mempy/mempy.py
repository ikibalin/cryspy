import cryspy

from .mempy_by_dictionary import mempy_reconstruction_by_dictionary, mempy_cycle_density_susceptibility

def mempy_spin_density_reconstruction(obj: cryspy.GlobalN):
    if not(obj.is_attribute("mem_parameters")):
        mem_parameters = cryspy.MEMParameters(
            points_a=48, points_b=48, points_c=48,
            channel_col=True, magnetization_plus=4, magnetization_minus=-1)
        obj.add_items([mem_parameters,])
    obj.mem_parameters.channel_col = True
    obj.mem_parameters.channel_ani = False

    res = mempy_reconstruction_with_parameters(obj)
    return res 


def mempy_magnetization_density_reconstruction(obj: cryspy.GlobalN):
    if not(obj.is_attribute("mem_parameters")):
        mem_parameters = cryspy.MEMParameters(
            points_a=48, points_b=48, points_c=48, channel_ani=True, only_magnetic_basins=True)
        obj.add_items([mem_parameters,])
    obj.mem_parameters.channel_col = False
    obj.mem_parameters.channel_ani = True

    res = mempy_reconstruction_with_parameters(obj)
    return res  

def mempy_reconstruction_with_parameters(obj: cryspy.GlobalN,
        parameter_lambda:float=1.e-5, iteration_max:int=1000, parameter_lambda_min:float=1.e-9, delta_density:float=1.e-5):
    l_diffrn = [item for item in obj.items if isinstance(item, cryspy.Diffrn)]

    dict_obj = obj.get_dictionary()
    l_dict_crystal = [value for name, value in dict_obj.items() if name.startswith("crystal_")]
    l_dict_mem_parameters = [value for name, value in dict_obj.items() if name.startswith("mem_parameters")]
    l_dict_diffrn = [value for name, value in dict_obj.items() if name.startswith("diffrn_")]

    if len(l_dict_crystal) == 0:
        raise KeyError("The crystal is not found")

    elif len(l_dict_crystal) > 1:
        raise UserWarning("The first crystal will be used for MEM reconstruction procedure")

    if len(l_dict_diffrn) == 0:
        raise KeyError("The single crystall polarized experiments are not found")

    dict_crystal = l_dict_crystal[0]
    dict_mem_parameters = l_dict_mem_parameters[0]
    dict_in_out = {}
    mempy_reconstruction_by_dictionary(dict_crystal, dict_mem_parameters, l_dict_diffrn, dict_in_out,
        parameter_lambda=parameter_lambda, iteration_max=iteration_max,
        parameter_lambda_min=parameter_lambda_min, delta_density=delta_density)

    dict_in_out_keys = dict_in_out.keys()

    for diffrn in l_diffrn:
        flip_ratio = dict_in_out["dict_in_out_"+diffrn.get_name()]["flip_ratio"]
        diffrn_refln = diffrn.diffrn_refln
        diffrn_refln.numpy_fr_calc = flip_ratio
        diffrn_refln.numpy_to_items()

    if (("symm_elem_channel_ani" in dict_in_out_keys) and
            ("density_channel_ani" in dict_in_out_keys) and
            (("susceptibility_channel_ani" in dict_in_out_keys) or ("moment_channel_ani" in dict_in_out_keys)) and
            ("atom_multiplicity_channel_ani" in dict_in_out_keys) and
            ("point_multiplicity_channel_ani" in dict_in_out_keys)):

        symm_elem_channel_ani = dict_in_out["symm_elem_channel_ani"]  
        density_channel_ani = dict_in_out["density_channel_ani"]
        
        atom_multiplicity_channel_ani = dict_in_out["atom_multiplicity_channel_ani"]  
        point_multiplicity_channel_ani = dict_in_out["point_multiplicity_channel_ani"]  

        channel_ani = cryspy.ChannelAniL(loop_name=dict_crystal["name"])

        channel_ani.numpy_numerator_x = symm_elem_channel_ani[0]
        channel_ani.numpy_numerator_y = symm_elem_channel_ani[1]
        channel_ani.numpy_numerator_z = symm_elem_channel_ani[2]
        channel_ani.numpy_denominator_xyz = symm_elem_channel_ani[3]
        channel_ani.numpy_density = density_channel_ani
        channel_ani.numpy_atom_multiplicity = atom_multiplicity_channel_ani
        channel_ani.numpy_point_multiplicity = point_multiplicity_channel_ani

        if "susceptibility_channel_ani" in dict_in_out_keys:
            susceptibility_channel_ani = dict_in_out["susceptibility_channel_ani"]
            channel_ani.numpy_chi_11 = susceptibility_channel_ani[0]
            channel_ani.numpy_chi_21 = susceptibility_channel_ani[1]
            channel_ani.numpy_chi_31 = susceptibility_channel_ani[2]
            channel_ani.numpy_chi_12 = susceptibility_channel_ani[3]
            channel_ani.numpy_chi_22 = susceptibility_channel_ani[4]
            channel_ani.numpy_chi_32 = susceptibility_channel_ani[5]
            channel_ani.numpy_chi_13 = susceptibility_channel_ani[6]
            channel_ani.numpy_chi_23 = susceptibility_channel_ani[7]
            channel_ani.numpy_chi_33 = susceptibility_channel_ani[8]

        if "moment_channel_ani" in dict_in_out_keys:
            moment_channel_ani = dict_in_out["moment_channel_ani"]
            channel_ani.numpy_m_x = moment_channel_ani[0]
            channel_ani.numpy_m_y = moment_channel_ani[1]
            channel_ani.numpy_m_z = moment_channel_ani[2]

        channel_ani.numpy_to_items()
        obj.add_items([channel_ani, ])

    if (("symm_elem_channel_col" in dict_in_out_keys) and
            ("magnetization_plus" in dict_in_out_keys) and
            ("magnetization_minus" in dict_in_out_keys) and
            ("density_channel_col" in dict_in_out_keys) and
            ("point_multiplicity_channel_col" in dict_in_out_keys)):

        symm_elem_channel_col = dict_in_out["symm_elem_channel_col"]
        density_channel_col = dict_in_out["density_channel_col"]
        magnetization_plus = dict_in_out["magnetization_plus"]
        magnetization_minus = dict_in_out["magnetization_minus"]
        point_multiplicity_channel_col = dict_in_out["point_multiplicity_channel_col"]


        channel_col = cryspy.ChannelColL(loop_name=dict_crystal["name"])

        channel_col.numpy_numerator_x = symm_elem_channel_col[0]
        channel_col.numpy_numerator_y = symm_elem_channel_col[1]
        channel_col.numpy_numerator_z = symm_elem_channel_col[2]
        channel_col.numpy_denominator_xyz = symm_elem_channel_col[3]
        channel_col.numpy_density_plus = density_channel_col[0]
        channel_col.numpy_density_minus = density_channel_col[1]
        channel_col.numpy_point_multiplicity = point_multiplicity_channel_col
        channel_col.numpy_to_items()
        obj.add_items([channel_col, ])
    return 



def mempy_cycle_with_parameters(obj: cryspy.GlobalN,
        parameter_lambda:float=1.e-5, iteration_max:int=1000, parameter_lambda_min:float=1.e-9, delta_density:float=1.e-5,
        n_cycle:int=10):

    l_diffrn = [item for item in obj.items if isinstance(item, cryspy.Diffrn)]

    dict_obj = obj.get_dictionary()
    l_dict_crystal = [value for name, value in dict_obj.items() if name.startswith("crystal_")]
    l_dict_mem_parameters = [value for name, value in dict_obj.items() if name.startswith("mem_parameters")]
    l_dict_diffrn = [value for name, value in dict_obj.items() if name.startswith("diffrn_")]

    if len(l_dict_crystal) == 0:
        raise KeyError("The crystal is not found")

    elif len(l_dict_crystal) > 1:
        raise UserWarning("The first crystal will be used for MEM reconstruction procedure")

    if len(l_dict_diffrn) == 0:
        raise KeyError("The single crystall polarized experiments are not found")

    dict_crystal = l_dict_crystal[0]
    dict_mem_parameters = l_dict_mem_parameters[0]
    

    dict_in_out = {}
    mempy_cycle_density_susceptibility(dict_crystal, dict_mem_parameters, l_dict_diffrn, dict_in_out,
        parameter_lambda=parameter_lambda, iteration_max=iteration_max, parameter_lambda_min=parameter_lambda_min,
        delta_density=delta_density, n_cycle=n_cycle)

    dict_in_out_den = dict_in_out["dict_in_out_den"]
    dict_in_out_susc = dict_in_out["dict_in_out_susc"]

    if "atom_para_susceptibility_sigma" in dict_in_out_susc.keys():
        atom_para_susceptibility_sigma = dict_in_out_susc["atom_para_susceptibility_sigma"]
        atom_para_susceptibility = dict_crystal["atom_para_susceptibility"]
        flags_atom_para_susceptibility = dict_crystal["flags_atom_para_susceptibility"]
        # obj.take_parameters_from_dictionary(dict_obj, l_parameter_name = parameter_name, l_sigma=sigma_p)
        #FIXME: temporary solution
        obj_crystal = getattr(obj, f"crystal_{dict_crystal['name']:}")
        atom_site_susceptibility = obj_crystal.atom_site_susceptibility
        l_attr = ["chi_11", "chi_22", "chi_33", "chi_12", "chi_13", "chi_23"]
        for i_item, item in enumerate(atom_site_susceptibility.items):
            for i_attr, attr in enumerate(l_attr):
                setattr(item, attr, atom_para_susceptibility[i_attr, i_item])
                if flags_atom_para_susceptibility[i_attr, i_item]:
                    setattr(item, f"{attr:}_sigma", atom_para_susceptibility_sigma[i_attr, i_item])


    dict_in_out_den_keys = dict_in_out_den.keys()

    for diffrn in l_diffrn:
        flip_ratio = dict_in_out_den["dict_in_out_"+diffrn.get_name()]["flip_ratio"]
        diffrn_refln = diffrn.diffrn_refln
        diffrn_refln.numpy_fr_calc = flip_ratio
        diffrn_refln.numpy_to_items()

    if (("symm_elem_channel_ani" in dict_in_out_den_keys) and
            ("density_channel_ani" in dict_in_out_den_keys) and
            (("susceptibility_channel_ani" in dict_in_out_den_keys) or ("moment_channel_ani" in dict_in_out_den_keys)) and
            ("atom_multiplicity_channel_ani" in dict_in_out_den_keys) and
            ("point_multiplicity_channel_ani" in dict_in_out_den_keys)):

        symm_elem_channel_ani = dict_in_out_den["symm_elem_channel_ani"]  
        density_channel_ani = dict_in_out_den["density_channel_ani"]
        susceptibility_channel_ani = dict_in_out_den["susceptibility_channel_ani"]
        atom_multiplicity_channel_ani = dict_in_out_den["atom_multiplicity_channel_ani"]  
        point_multiplicity_channel_ani = dict_in_out_den["point_multiplicity_channel_ani"]  

        channel_ani = cryspy.ChannelAniL(loop_name=dict_crystal["name"])

        channel_ani.numpy_numerator_x = symm_elem_channel_ani[0]
        channel_ani.numpy_numerator_y = symm_elem_channel_ani[1]
        channel_ani.numpy_numerator_z = symm_elem_channel_ani[2]
        channel_ani.numpy_denominator_xyz = symm_elem_channel_ani[3]
        channel_ani.numpy_density = density_channel_ani
        
        if "moment_channel_ani" in dict_in_out_den_keys:
            moment_channel_ani = dict_in_out_den["moment_channel_ani"]
            channel_ani.numpy_m_x = moment_channel_ani[0]
            channel_ani.numpy_m_y = moment_channel_ani[1]
            channel_ani.numpy_m_z = moment_channel_ani[2]
        
        if "susceptibility_channel_ani" in dict_in_out_den_keys:
            susceptibility_channel_ani = dict_in_out_den["susceptibility_channel_ani"]
            channel_ani.numpy_chi_11 = susceptibility_channel_ani[0]
            channel_ani.numpy_chi_21 = susceptibility_channel_ani[1]
            channel_ani.numpy_chi_31 = susceptibility_channel_ani[2]
            channel_ani.numpy_chi_12 = susceptibility_channel_ani[3]
            channel_ani.numpy_chi_22 = susceptibility_channel_ani[4]
            channel_ani.numpy_chi_32 = susceptibility_channel_ani[5]
            channel_ani.numpy_chi_13 = susceptibility_channel_ani[6]
            channel_ani.numpy_chi_23 = susceptibility_channel_ani[7]
            channel_ani.numpy_chi_33 = susceptibility_channel_ani[8]
        
        channel_ani.numpy_atom_multiplicity = atom_multiplicity_channel_ani
        channel_ani.numpy_point_multiplicity = point_multiplicity_channel_ani
        channel_ani.numpy_to_items()
        obj.add_items([channel_ani, ])

    if (("symm_elem_channel_col" in dict_in_out_den_keys) and
            ("magnetization_plus" in dict_in_out_den_keys) and
            ("magnetization_minus" in dict_in_out_den_keys) and
            ("density_channel_col" in dict_in_out_den_keys) and
            ("point_multiplicity_channel_col" in dict_in_out_den_keys)):

        symm_elem_channel_col = dict_in_out_den["symm_elem_channel_col"]
        density_channel_col = dict_in_out_den["density_channel_col"]
        magnetization_plus = dict_in_out_den["magnetization_plus"]
        magnetization_minus = dict_in_out_den["magnetization_minus"]
        point_multiplicity_channel_col = dict_in_out_den["point_multiplicity_channel_col"]


        channel_col = cryspy.ChannelColL(loop_name=dict_crystal["name"])

        channel_col.numpy_numerator_x = symm_elem_channel_col[0]
        channel_col.numpy_numerator_y = symm_elem_channel_col[1]
        channel_col.numpy_numerator_z = symm_elem_channel_col[2]
        channel_col.numpy_denominator_xyz = symm_elem_channel_col[3]
        channel_col.numpy_density_plus = density_channel_col[0]
        channel_col.numpy_density_minus = density_channel_col[1]
        channel_col.numpy_point_multiplicity = point_multiplicity_channel_col
        channel_col.numpy_to_items()
        obj.add_items([channel_col, ])
    print(70*"*")
    print("End of cycle procedure (MEMPy module)")
    print(70*"*")
