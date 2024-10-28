from typing import List, Union
import numpy
import scipy
import scipy.optimize

from warnings import warn


from cryspy.A_functions_base.function_1_inversed_hessian import \
    estimate_inversed_hessian_matrix
from cryspy.A_functions_base.function_1_error_simplex import \
    error_estimation_simplex

from cryspy.B_parent_classes.cl_3_data import DataN
from cryspy.B_parent_classes.cl_4_global import GlobalN

from cryspy.C_item_loop_classes.cl_1_inversed_hessian import InversedHessian, inversed_hessian_to_correlation

from cryspy.E_data_classes.cl_1_crystal import Crystal
# from cryspy.E_data_classes.cl_1_mag_crystal import MagCrystal
from cryspy.E_data_classes.cl_2_diffrn import Diffrn
from cryspy.E_data_classes.cl_2_pd import Pd
from cryspy.E_data_classes.cl_2_pd2d import Pd2d
from cryspy.E_data_classes.cl_2_tof import TOF

from cryspy.procedure_rhochi.rhochi_by_dictionary import \
    rhochi_lsq_by_dictionary, rhochi_rietveld_refinement_by_dictionary,\
    rhochi_calc_chi_sq_by_dictionary

import cryspy

na = numpy.newaxis

def rhochi_check_items(cryspy_object: cryspy.GlobalN):
    """Check items for RhoChi refinement procedure.

    Output is list of experiments and list of crystals taken from cryspy_object
    """
    f_crystal, f_experiment = False, False

    if not(isinstance(cryspy_object, cryspy.GlobalN)):
        raise AttributeError("Incorrect type of object")

    cl_crysts = (cryspy.Crystal, ) #cryspy.MagCrystal
    cl_exps = (cryspy.Diffrn, cryspy.Pd, cryspy.Pd2d, cryspy.TOF)
    for item in cryspy_object.items:

        if isinstance(item, cl_crysts):
            f_crystal = True

        if isinstance(item, cl_exps):
            f_experiment = True

        if f_crystal & f_experiment:
            break

    if not(f_crystal):
        raise AttributeError("Crystal is not defined")

    if not(f_experiment):
        raise AttributeError("Experiment is not defined")



def rhochi_rietveld_refinement(cryspy_object: cryspy.GlobalN) -> dict:
    """Run refinement by RhoChi procedure.
    """
    # check object
    rhochi_check_items(cryspy_object)

    method = "BFGS"

    dict_out = rhochi_rietveld_refinement_with_parameters(
        cryspy_object,
        optimization_method=method) 

    return dict_out


def rhochi_rietveld_refinement_with_parameters(
        cryspy_object: cryspy.GlobalN,
        optimization_method: str = "BFGS") -> dict:
    """Run refinement by RhoChi procedure with non-default parameters.
    """

    # check object
    rhochi_check_items(cryspy_object)

    obj_dict = cryspy_object.get_dictionary()
    flag_scipy_refinements = True
    if flag_scipy_refinements:
        
        DICT_PARAMS["previous_arg"] = ()
        DICT_PARAMS["iteration"] = 0
        chi_sq, parameter_name, dict_in_out, res = rhochi_rietveld_refinement_by_dictionary(
            obj_dict, method=optimization_method, callback=_f_callback)
        dict_out = {"chi_sq": chi_sq, "parameter_name": parameter_name}
        if "hess_inv" in res.keys():
            hess_inv = res["hess_inv"]
            sigma_p = numpy.sqrt(numpy.abs(numpy.diag(hess_inv)))
            correlation_matrix = hess_inv/(sigma_p[:, na]*sigma_p[na, :])
            dict_out["correlation_matrix"] = correlation_matrix

            l_label = [hh[-1][0] for hh in parameter_name]
            if len(l_label) == hess_inv.shape[0]:
                inv_hessian = InversedHessian()
                inv_hessian.set_labels(l_label)
                inv_hessian.set_inversed_hessian(hess_inv)
                inv_hessian.form_inversed_hessian()
                inv_hessian.form_object()
                cryspy_object.add_items([inv_hessian, ])
        else:
            sigma_p = numpy.zeros((len(parameter_name),), dtype=float)
    else:
        chi_sq, delta_p, parameter_name, der_chi_sq, dder_chi_sq, dict_in_out = rhochi_lsq_by_dictionary(obj_dict)
        hessian = numpy.linalg.inv(dder_chi_sq)
        sigma_p = numpy.sqrt(numpy.diag(hessian))
        correlation_matrix = hessian/(sigma_p[:, na]*sigma_p[na, :])
        sigma_p[numpy.isnan(sigma_p)] = 0.
        correlation_matrix[numpy.isnan(correlation_matrix)] = 0.
        dict_out = {"chi_sq": chi_sq, "parameter_name": parameter_name,"der_chi_sq":der_chi_sq,
            "dder_chi_sq": dder_chi_sq, "correlation_matrix": correlation_matrix}
    cryspy_object.take_parameters_from_dictionary(obj_dict, l_parameter_name = parameter_name, l_sigma=sigma_p)
    cryspy_object.take_parameters_from_dictionary(dict_in_out, l_parameter_name = None, l_sigma=None)
    
    var_names = cryspy_object.get_variable_names()
    if len(var_names) > 0:
        print("Optimal parameters:")
    for name in var_names:
        value = cryspy_object.get_variable_by_name(name)
        print(f" - {name[-1][0]:} {value:.5f}")
    
    print("Errors are estimating by numerical derivatives.")
    rhochi_inversed_hessian(cryspy_object)
    return dict_out


def rhochi_no_refinement(cryspy_object: cryspy.GlobalN) -> dict:
    """Run calculations by RhoChi procedure.
    """
    flag_calc_analytical_derivatives = False
    flag_use_precalculated_data = False

    obj_dict = cryspy_object.get_dictionary()
    dict_in_out = {}

    chi_sq, n_point, der_chi_sq, dder_chi_sq, parameter_names = rhochi_calc_chi_sq_by_dictionary(
        obj_dict,
        dict_in_out=dict_in_out,
        flag_use_precalculated_data=flag_use_precalculated_data, flag_calc_analytical_derivatives=flag_calc_analytical_derivatives)
    
    cryspy_object.take_parameters_from_dictionary(obj_dict, l_parameter_name = None, l_sigma=None)
    cryspy_object.take_parameters_from_dictionary(dict_in_out, l_parameter_name = None, l_sigma=None)
    
    dict_out = {"chi_sq": chi_sq, "n_point": n_point}
    return dict_out

DICT_PARAMS = {"previous_arg": (), "iteration": 0}

def _f_callback(*arg, d_info: dict = None) -> bool:
    flag_out = False
    res_x = arg[0]
    if len(DICT_PARAMS["previous_arg"]) != len(res_x):
        DICT_PARAMS["previous_arg"] = res_x
    else:
        DICT_PARAMS["iteration"] += 1
        diff = numpy.array(res_x, dtype=float) - numpy.array(DICT_PARAMS["previous_arg"], dtype=float)
        shift = numpy.sqrt(numpy.square(diff).sum()) * 100
        print(f"Average shift of parameters is {shift:.5f} ({DICT_PARAMS['iteration']:}).                ", end="\r")
        DICT_PARAMS["previous_arg"] = res_x
    return flag_out
    # ls_out = ["{:12.5f}".format(_1) for _1 in res_x]
    # print(" ".join(ls_out), end="\r")
    # return flag_out


def rhochi_inversed_hessian(global_object: GlobalN):
    """Estimate inversed Hessian matrix."""
    if global_object.is_attribute("inversed_hessian"):
        global_object.items.remove(global_object.inversed_hessian)

    global_dict = global_object.get_dictionary()

    flag_calc_analytical_derivatives = False
    flag_use_precalculated_data = False

    obj_dict = global_object.get_dictionary()
    dict_in_out = {}

    chi_sq, n_point, der_chi_sq, dder_chi_sq, parameter_names = rhochi_calc_chi_sq_by_dictionary(
        obj_dict,
        dict_in_out=dict_in_out,
        flag_use_precalculated_data=flag_use_precalculated_data, flag_calc_analytical_derivatives=flag_calc_analytical_derivatives)
    
    param_0 = [global_dict[way[0]][way[1]][way[2]] for way  in parameter_names]

    flag_use_precalculated_data = True
    def tempfunc(l_param):
        for way, param in zip(parameter_names, l_param):
            global_dict[way[0]][way[1]][way[2]] = param
        chi_sq = rhochi_calc_chi_sq_by_dictionary(
            global_dict,
            dict_in_out=dict_in_out,
            flag_use_precalculated_data=flag_use_precalculated_data,
            flag_calc_analytical_derivatives=flag_calc_analytical_derivatives)[0]
        return chi_sq

    hess_inv, np_first_der = estimate_inversed_hessian_matrix(tempfunc, param_0)
    if (numpy.all(hess_inv == numpy.zeros_like(hess_inv)) or (hess_inv is None)):
        return 
    corr_matrix, sigmas = inversed_hessian_to_correlation(hess_inv) 
    global_object.take_parameters_from_dictionary(
        global_dict, l_parameter_name = parameter_names, l_sigma=sigmas)
    
    l_label = []
    for way in parameter_names:
        way_1, way_2, way_3 = way[0], way[1], way[2]
        l_h = way_1.split("_")
        if len(l_h) > 1:
            s_1 = "_".join(l_h[1:])
        else:
            s_1 = l_h[0]

        l_h = way_2.split("_")
        if len(l_h) > 1:
            s_2 = "_".join(l_h[1:])
        else:
            s_2 = l_h[0]

        s_3 = str(way_3).replace(" ", "")

        l_label.append(f"{s_1:},{s_2:},{s_3:}")


    # l_label = [f"{str(var_name).replace(' ',''):}" for i_param, var_name in enumerate(parameter_names)]

    inv_hessian = InversedHessian()
    inv_hessian.set_labels(l_label)
    inv_hessian.set_inversed_hessian(hess_inv)
    inv_hessian.form_inversed_hessian()
    inv_hessian.form_object()

    global_object.add_items([inv_hessian, ])

    return inv_hessian
