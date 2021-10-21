import cryspy

from .mempy_by_dictionary import mempy_reconstruction_by_dictionary

def mempy_spin_density_reconstruction(obj: cryspy.GlobalN):
    if not(obj.is_attribute("mem_parameters")):
        mem_parameters = cryspy.MEMParameters(
            points_a=48, points_b=48, points_c=48,
            channel_plus_minus=True, magnetization_plus=4, magnetization_minus=-1)
        obj.add_items([mem_parameters,])
    obj.mem_parameters.channel_plus_minus = True
    obj.mem_parameters.channel_chi = False

    res = mempy_reconstruction_with_parameters(obj)
    return res 


def mempy_magnetization_density_reconstruction(obj: cryspy.GlobalN):
    if not(obj.is_attribute("mem_parameters")):
        mem_parameters = cryspy.MEMParameters(
            points_a=48, points_b=48, points_c=48, channel_chi=True, only_magnetic_basins=True)
        obj.add_items([mem_parameters,])
    obj.mem_parameters.channel_plus_minus = False
    obj.mem_parameters.channel_chi = True

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

    if (("symm_elem_channel_chi" in dict_in_out_keys) and
            ("density_channel_chi" in dict_in_out_keys) and
            ("susceptibility_channel_chi" in dict_in_out_keys) and
            ("atom_multiplicity_channel_chi" in dict_in_out_keys) and
            ("point_multiplicity_channel_chi" in dict_in_out_keys)):

        symm_elem_channel_chi = dict_in_out["symm_elem_channel_chi"]  
        density_channel_chi = dict_in_out["density_channel_chi"]
        susceptibility_channel_chi = dict_in_out["susceptibility_channel_chi"]
        atom_multiplicity_channel_chi = dict_in_out["atom_multiplicity_channel_chi"]  
        point_multiplicity_channel_chi = dict_in_out["point_multiplicity_channel_chi"]  

        channel_chi = cryspy.ChannelChiL(loop_name=dict_crystal["name"])

        channel_chi.numpy_numerator_x = symm_elem_channel_chi[0]
        channel_chi.numpy_numerator_y = symm_elem_channel_chi[1]
        channel_chi.numpy_numerator_z = symm_elem_channel_chi[2]
        channel_chi.numpy_denominator_xyz = symm_elem_channel_chi[3]
        channel_chi.numpy_density = density_channel_chi
        channel_chi.numpy_chi_11 = susceptibility_channel_chi[0]
        channel_chi.numpy_chi_21 = susceptibility_channel_chi[1]
        channel_chi.numpy_chi_31 = susceptibility_channel_chi[2]
        channel_chi.numpy_chi_12 = susceptibility_channel_chi[3]
        channel_chi.numpy_chi_22 = susceptibility_channel_chi[4]
        channel_chi.numpy_chi_32 = susceptibility_channel_chi[5]
        channel_chi.numpy_chi_13 = susceptibility_channel_chi[6]
        channel_chi.numpy_chi_23 = susceptibility_channel_chi[7]
        channel_chi.numpy_chi_33 = susceptibility_channel_chi[8]
        channel_chi.numpy_atom_multiplicity = atom_multiplicity_channel_chi
        channel_chi.numpy_point_multiplicity = point_multiplicity_channel_chi
        channel_chi.numpy_to_items()
        obj.add_items([channel_chi, ])

    if (("symm_elem_channel_plus_minus" in dict_in_out_keys) and
            ("magnetization_plus" in dict_in_out_keys) and
            ("magnetization_minus" in dict_in_out_keys) and
            ("density_channel_plus_minus" in dict_in_out_keys) and
            ("point_multiplicity_channel_plus_minus" in dict_in_out_keys)):

        symm_elem_channel_plus_minus = dict_in_out["symm_elem_channel_plus_minus"]
        density_channel_plus_minus = dict_in_out["density_channel_plus_minus"]
        magnetization_plus = dict_in_out["magnetization_plus"]
        magnetization_minus = dict_in_out["magnetization_minus"]
        point_multiplicity_channel_plus_minus = dict_in_out["point_multiplicity_channel_plus_minus"]


        channel_plus_minus = cryspy.ChannelPlusMinusL(loop_name=dict_crystal["name"])

        channel_plus_minus.numpy_numerator_x = symm_elem_channel_plus_minus[0]
        channel_plus_minus.numpy_numerator_y = symm_elem_channel_plus_minus[1]
        channel_plus_minus.numpy_numerator_z = symm_elem_channel_plus_minus[2]
        channel_plus_minus.numpy_denominator_xyz = symm_elem_channel_plus_minus[3]
        channel_plus_minus.numpy_density_plus = density_channel_plus_minus[0]
        channel_plus_minus.numpy_density_minus = density_channel_plus_minus[1]
        channel_plus_minus.numpy_point_multiplicity = point_multiplicity_channel_plus_minus
        channel_plus_minus.numpy_to_items()
        obj.add_items([channel_plus_minus, ])
    return 
