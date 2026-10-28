import logging
import numpy



from cryspy.A_functions_base.function_1_inversed_hessian import \
    estimate_inversed_hessian_matrix

from cryspy.B_parent_classes.cl_4_global import GlobalN
from cryspy.exceptions import RefinementError

from cryspy.C_item_loop_classes.cl_1_inversed_hessian import InversedHessian, inversed_hessian_to_correlation

# from cryspy.E_data_classes.cl_1_mag_crystal import MagCrystal

from cryspy.procedure_rhochi.rhochi_by_dictionary import \
    rhochi_lsq_by_dictionary, rhochi_rietveld_refinement_by_dictionary,\
    rhochi_calc_chi_sq_by_dictionary, split_parameters_by_constraints, \
    set_dependent_parameters, calc_punishement_function

import cryspy
from cryspy.procedure_assistant.parameters import (
    is_number, take_dict_paths, take_parameter_records)

from cryspy.A_functions_base.progress import PROGRESS

logger = logging.getLogger(__name__)

na = numpy.newaxis

def rhochi_check_items(cryspy_object: cryspy.GlobalN):
    """Check items for RhoChi refinement procedure.

    Output is list of experiments and list of crystals taken from cryspy_object
    """
    f_crystal, f_experiment = False, False

    if not(isinstance(cryspy_object, cryspy.GlobalN)):
        raise RefinementError("Incorrect type of object: GlobalN is expected")

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
        raise RefinementError("Crystal is not defined")

    if not(f_experiment):
        raise RefinementError("Experiment is not defined")



# `callback` is keyword-only and has no annotation on purpose: CrysPy Editor
# builds its menus and argument widgets from the positional arguments and
# their annotations, so the procedure stays a one-click (autorun) procedure.
def rhochi_rietveld_refinement(cryspy_object: cryspy.GlobalN, *, callback=None) -> dict:
    """Run refinement by RhoChi procedure.

    callback (keyword-only, optional) is called after every iteration of the
    minimizer as callback(info); info is a dict with "iteration", "chi_sq",
    "n_point", "chi_sq_per_n", "parameter_names" and "parameter_values".
    Returning False stops the refinement: the object keeps the parameters
    of the last completed iteration and the output dict has "stopped": True.
    """
    # check object
    rhochi_check_items(cryspy_object)

    method = "BFGS"

    dict_out = rhochi_rietveld_refinement_with_parameters(
        cryspy_object,
        optimization_method=method, callback=callback)

    return dict_out


def rhochi_rietveld_refinement_with_parameters(
        cryspy_object: cryspy.GlobalN,
        optimization_method: str = "BFGS", *, callback=None) -> dict:
    """Run refinement by RhoChi procedure with non-default parameters.

    callback (keyword-only, optional): see rhochi_rietveld_refinement.

    Parameters
    ----------
    cryspy_object : GlobalN
        The object with the crystals and the experiments.
    optimization_method : str, optional
        Minimization method: 'BFGS', 'basinhopping' or other methods of
        scipy.optimize.minimize
        (https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html).
    """

    # check object
    rhochi_check_items(cryspy_object)

    obj_dict = cryspy_object.get_dictionary()
    flag_scipy_refinements = True
    flag_stopped = False
    if flag_scipy_refinements:
        
        DICT_PARAMS["previous_arg"] = ()
        DICT_PARAMS["iteration"] = 0
        d_path = take_dict_paths(cryspy_object, obj_dict)
        chi_sq, parameter_name, dict_in_out, res = rhochi_rietveld_refinement_by_dictionary(
            obj_dict, method=optimization_method, callback=_f_callback,
            progress_callback=callback, parameter_paths=d_path)
        flag_stopped = bool(res.get("stopped", False))
        # paths as in the log of the refinement, in the order of parameter_name
        dict_out = {"chi_sq": chi_sq, "parameter_name": parameter_name,
                    "parameter_path": form_hessian_labels(parameter_name, d_path),
                    "stopped": flag_stopped}
        if "hess_inv" in res.keys():
            hess_inv = res["hess_inv"]
            sigma_p = numpy.sqrt(numpy.abs(numpy.diag(hess_inv)))
            correlation_matrix = hess_inv/(sigma_p[:, na]*sigma_p[na, :])
            dict_out["correlation_matrix"] = correlation_matrix

            l_label = dict_out["parameter_path"]
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
    if flag_stopped:
        # sigmas of a stopped refinement are unknown: keep the previous ones
        cryspy_object.take_parameters_from_dictionary(obj_dict, l_parameter_name = None, l_sigma=None)
    else:
        cryspy_object.take_parameters_from_dictionary(obj_dict, l_parameter_name = parameter_name, l_sigma=sigma_p)
    cryspy_object.take_parameters_from_dictionary(dict_in_out, l_parameter_name = None, l_sigma=None)
    
    if flag_stopped:
        log_refined_parameters(cryspy_object)
        return dict_out
    logger.info("Errors are estimating by numerical derivatives.")
    rhochi_inversed_hessian(cryspy_object)
    log_refined_parameters(cryspy_object)
    return dict_out


def log_refined_parameters(cryspy_object: GlobalN):
    """Log the refined parameters by their full paths
    ('crystal_fe.atom_site[Fe1].occupancy'), with sigmas and marks."""
    l_record = [record for record in take_parameter_records(cryspy_object)
                if record.flag_refinement]
    if len(l_record) == 0:
        return
    logger.info("Optimal parameters:")
    for record in l_record:
        mark = str(getattr(record.item, f"{record.attribute:}_mark", ""))
        s_mark = f" [{mark:}]" if mark != "" else ""
        if is_number(record.value):
            s_value = f"{record.value:.5f} ± {record.sigma:.5f}"
        else:
            s_value = str(record.value)
        logger.info(f" - {record.path:} {s_value:}{s_mark:}")


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
        logger.info(f"Average shift of parameters is {shift:.5f} ({DICT_PARAMS['iteration']:}).", extra=PROGRESS)
        DICT_PARAMS["previous_arg"] = res_x
    return flag_out
    # ls_out = ["{:12.5f}".format(_1) for _1 in res_x]
    # print(" ".join(ls_out), end="\r")
    # return flag_out


def form_hessian_labels(parameter_names, parameter_paths: dict = None) -> list:
    """Labels of the parameters of the inversed Hessian matrix.

    The label is the path of the parameter as in the log of the refinement
    ('crystal_fe.atom_site[Fe1].occupancy'). A parameter without a known
    path keeps the old label 'fe,atom_occupancy,(0,)'. Spaces are removed:
    the labels are separated by spaces in the rcif text of the matrix.
    """
    l_label = []
    for way in parameter_names:
        way_sh = (way[0], way[1], tuple(int(hh) for hh in way[2]))
        if parameter_paths is not None and way_sh in parameter_paths:
            label = parameter_paths[way_sh]
        else:
            l_h = str(way_sh[0]).split("_")
            s_1 = "_".join(l_h[1:]) if len(l_h) > 1 else l_h[0]
            l_h = str(way_sh[1]).split("_")
            s_2 = "_".join(l_h[1:]) if len(l_h) > 1 else l_h[0]
            label = f"{s_1:},{s_2:},{way_sh[2]:}"
        l_label.append("".join(str(label).split()))
    return l_label


def rhochi_inversed_hessian(global_object: GlobalN):
    """Estimate inversed Hessian matrix."""
    if global_object.is_attribute("inversed_hessian"):
        global_object.items.remove(global_object.inversed_hessian)

    global_dict = global_object.get_dictionary()

    flag_calc_analytical_derivatives = False
    flag_use_precalculated_data = False

    obj_dict = global_object.get_dictionary()
    dict_in_out = {}
    # paths are found before the values in global_dict are changed
    d_path = take_dict_paths(global_object, global_dict)

    chi_sq, n_point, der_chi_sq, dder_chi_sq, parameter_names = rhochi_calc_chi_sq_by_dictionary(
        obj_dict,
        dict_in_out=dict_in_out,
        flag_use_precalculated_data=flag_use_precalculated_data, flag_calc_analytical_derivatives=flag_calc_analytical_derivatives)
    
    # constrained parameters (equal marks, constraint loop) follow the free ones
    names_free, names_dep, matrix_m, vector_c = split_parameters_by_constraints(
        global_dict, parameter_names)
    param_0 = [global_dict[way[0]][way[1]][way[2]] for way  in names_free]

    flag_use_precalculated_data = True
    def tempfunc(l_param):
        for way, param in zip(names_free, l_param):
            global_dict[way[0]][way[1]][way[2]] = param
        set_dependent_parameters(global_dict, names_dep, matrix_m, vector_c, l_param)
        chi_sq = rhochi_calc_chi_sq_by_dictionary(
            global_dict,
            dict_in_out=dict_in_out,
            flag_use_precalculated_data=flag_use_precalculated_data,
            flag_calc_analytical_derivatives=flag_calc_analytical_derivatives)[0]
        # restraints narrow the errors as they narrow the minimum
        return chi_sq + calc_punishement_function(global_dict)

    hess_inv, np_first_der = estimate_inversed_hessian_matrix(tempfunc, param_0)
    if (hess_inv is None) or numpy.all(hess_inv == numpy.zeros_like(hess_inv)):
        return 
    if len(names_dep) > 0:
        # covariance of all parameters: p = J p_free
        jacobian = numpy.zeros((len(parameter_names), len(names_free)), dtype=float)
        for i_name, name in enumerate(parameter_names):
            if name in names_free:
                jacobian[i_name, names_free.index(name)] = 1.
            else:
                jacobian[i_name] = matrix_m[names_dep.index(name)]
        hess_inv = jacobian.dot(hess_inv).dot(jacobian.transpose())
    corr_matrix, sigmas = inversed_hessian_to_correlation(hess_inv) 
    global_object.take_parameters_from_dictionary(
        global_dict, l_parameter_name = parameter_names, l_sigma=sigmas)
    
    # the order of the log of the refinement ('Optimal parameters')
    l_label = form_hessian_labels(parameter_names, d_path)
    l_path = list(d_path.values())
    l_order = sorted(range(len(l_label)), key=lambda ind: (
        l_path.index(l_label[ind]) if l_label[ind] in l_path else len(l_path), ind))
    l_label = [l_label[ind] for ind in l_order]
    hess_inv = hess_inv[l_order, :][:, l_order]

    inv_hessian = InversedHessian()
    inv_hessian.set_labels(l_label)
    inv_hessian.set_inversed_hessian(hess_inv)
    inv_hessian.form_inversed_hessian()
    inv_hessian.form_object()

    global_object.add_items([inv_hessian, ])

    return inv_hessian
