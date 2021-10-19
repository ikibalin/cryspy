import cryspy

from .mempy_by_dictionary import mempy_reconstruction_by_dictionary

def mempy_spin_density_reconstruction(obj: cryspy.GlobalN, d_info:dict=None):
    if not(obj.is_attribute("mem_parameters")):
        mem_parameters = cryspy.MEMParameters(
            points_a=48, points_b=48, points_c=48,
            channel_plus_minus=True, magnetization_plus=4, magnetization_minus=-1)
        obj.add_items([mem_parameters,])
    obj.mem_parameters.channel_plus_minus = True
    obj.mem_parameters.channel_chi = False

    dict_obj = obj.get_dictionary()
    dict_obj_keys = dict_obj.keys()
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
    mempy_reconstruction_by_dictionary(dict_crystal, dict_mem_parameters, l_dict_diffrn, dict_in_out)
    return dict_in_out


def mempy_magnetization_density_reconstruction(obj: cryspy.GlobalN, d_info:dict=None):
    if not(obj.is_attribute("mem_parameters")):
        mem_parameters = cryspy.MEMParameters(
            points_a=48, points_b=48, points_c=48, channel_chi=True, only_magnetic_basins=True)
        obj.add_items([mem_parameters,])
    obj.mem_parameters.channel_plus_minus = False
    obj.mem_parameters.channel_chi = True


    dict_obj = obj.get_dictionary()
    dict_obj_keys = dict_obj.keys()
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
    mempy_reconstruction_by_dictionary(dict_crystal, dict_mem_parameters, l_dict_diffrn, dict_in_out)

    return dict_in_out

def mempy_reconstruction_with_parameters(obj: cryspy.GlobalN,
        parameter_lambda:float=1.e-5, iteration_max:int=1000, parameter_lambda_min:float=1.e-9, delta_density:float=1.e-5, d_info:dict=None):
    dict_obj = obj.get_dictionary()
    dict_obj_keys = dict_obj.keys()
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

    return dict_in_out
