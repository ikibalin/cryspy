import logging
import warnings
from typing import Callable
import re
import numpy
import scipy
import scipy.optimize

from .rhochi_diffrn import calc_chi_sq_for_diffrn_by_dictionary
from .rhochi_pd import calc_chi_sq_for_pd_by_dictionary
from .rhochi_pd2d import calc_chi_sq_for_pd2d_by_dictionary
from .rhochi_tof import calc_chi_sq_for_tof_by_dictionary

from cryspy.A_functions_base.linear_constraints import \
    calc_dependent_parameters, calc_soft_constraint_chi_sq, form_linear_system
from cryspy.exceptions import RefinementError

logger = logging.getLogger(__name__)

na = numpy.newaxis



#     if optimization_method == "basinhopping":
#         # basinhopping
#         res = scipy.optimize.basinhopping(
#             tempfunc, param_0, niter=10, T=10, stepsize=0.1, interval=20,
#             disp=disp)
#         _dict_out = {"flag": flag, "res": res}
#
#     elif optimization_method == "simplex":
#         # simplex
#         res = scipy.optimize.minimize(
#             tempfunc, param_0, method='Nelder-Mead',
#             callback=lambda x: _f_callback(
#                 disp, coeff_norm, x, param_name=l_var_name, d_info=d_info),
#             options={"fatol": 0.01*n})
#         m_error, dist_hh = error_estimation_simplex(
#             res["final_simplex"][0], res["final_simplex"][1], tempfunc)
#         l_sigma = []
#         for i, val_2 in zip(range(m_error.shape[0]), dist_hh):
#             # slightly change definition, instead of (n-k) here is n
#             error = (abs(m_error[i, i])*1./n)**0.5
#             if m_error[i, i] < 0.:
#                 warn("Negative diagonal elements of Hessian.", UserWarning)
#             if val_2 > error:
#                 warn("Minimum is not found.", UserWarning)
#             l_sigma.append(max(error, val_2))
#         obj_hm = cryspy.InversedHessian()
#         obj_hm.set_inversed_hessian(m_error*1./float(n))
#         l_label = [var_name[-1][0] for var_name in l_var_name]
#         obj_hm.set_labels(l_label)
#         obj_hm.form_inversed_hessian()
#         obj_hm.form_object()
#         cryspy_object.inversed_hessian = obj_hm
#         for var_name, sigma, coeff in \
#                 zip(l_var_name, l_sigma, coeff_norm):
#             hh = tuple((f"{var_name[-1][0]:}_sigma", ) + var_name[-1][1:])
#             var_name_sigma = tuple(var_name[:-1]+(hh, ))
#             cryspy_object.set_variable_by_name(var_name_sigma, sigma*coeff)
#         _dict_out = {"flag": flag, "res": res}

def calc_punishement_function(dict_obj):
    chi_sq_punishement = 0.
    if "punishment_function" in dict_obj.keys():
        expression, d_marked = dict_obj["punishment_function"]
        d_value = {}
        for mark in d_marked.keys():
            way = d_marked[mark]
            d_value[mark] = dict_obj[way[0]][way[1]][way[2]]
        chi_sq_punishement += eval(expression, {}, d_value)
    # soft constraints of the constraint loop: ((lhs - rhs)/sigma)**2
    for expression, sigma, d_mark_way in dict_obj.get("soft_constraints", []):
        d_value = {mark: sign * dict_obj[way[0]][way[1]][way[2]]
                   for mark, (sign, way) in d_mark_way.items()}
        chi_sq_punishement += calc_soft_constraint_chi_sq(expression, sigma, d_value)
    return chi_sq_punishement


def take_linear_equations(global_dict: dict, parameter_names: list) -> list:
    """Equations ``([(coeff, way), ...], const)`` between the refined
    parameters: equal marks (``a``, ``a-``) and the constraint loop
    (``x+y=1``)."""
    l_equation = []
    for (coeff_1, name_1), (coeff_2, name_2) in global_dict.get("linear_constraints", []):
        l_equation.append(([(coeff_1, name_1), (coeff_2, name_2)], 0.))
    l_equation.extend(global_dict.get("linear_equations", []))
    for l_term, const in l_equation:
        for coeff, name in l_term:
            if name not in parameter_names:
                raise RefinementError(
                    f"The constrained parameter {name:} is not refined.")
    return l_equation


def split_parameters_by_constraints(global_dict: dict, parameter_names: list):
    """Free and dependent refined parameters.

    Returns
    -------
    names_free, names_dep : list
        Ways of the free and of the dependent parameters.
    matrix_m : numpy.ndarray
        Shape (n_dep, n_free): ``p_dep = vector_c + matrix_m p_free``.
    vector_c : numpy.ndarray
        Shape (n_dep,).
    """
    l_equation = take_linear_equations(global_dict, parameter_names)
    if len(l_equation) == 0:
        return (list(parameter_names), [],
                numpy.zeros((0, len(parameter_names)), dtype=float),
                numpy.zeros((0,), dtype=float))
    matrix_a, vector_b = form_linear_system(l_equation, list(parameter_names))
    try:
        index_free, index_dep, matrix_m, vector_c = calc_dependent_parameters(matrix_a, vector_b)
    except ValueError as error:
        raise RefinementError(f"Constraints cannot be applied: {error:}.") from None
    names_free = [parameter_names[ind] for ind in index_free]
    names_dep = [parameter_names[ind] for ind in index_dep]
    return names_free, names_dep, matrix_m, vector_c


def set_dependent_parameters(global_dict: dict, names_dep: list, matrix_m, vector_c, param_free):
    """Write the dependent parameters calculated from the free ones."""
    if len(names_dep) == 0:
        return
    param_dep = vector_c + matrix_m.dot(numpy.asarray(param_free, dtype=float))
    for way, param in zip(names_dep, param_dep):
        global_dict[way[0]][way[1]][way[2]] = param


def calc_shift_p_by_constraints_hamilton(delta_p, dder_chi_sq, matrix_q):
    """Linear type Hamilton constraints.
    """
    q_t = matrix_q.transpose()
    q_dder_chi_sq = matrix_q.dot(dder_chi_sq)
    q_dder_chi_sq_qt = q_dder_chi_sq.dot(q_t)
    inv_q_dder_qt = numpy.linalg.inv(q_dder_chi_sq_qt)
    shift_delta_p = -delta_p.dot(q_t.dot(inv_q_dder_qt.dot(q_dder_chi_sq)))
    return shift_delta_p


def form_matrix_q(linear_constraints, parameter_name):
    matrix_q = numpy.zeros((len(linear_constraints), len(parameter_name)), dtype=float)
    for i_constraint, value_constraint in enumerate(linear_constraints):
        for coeff, p_name in value_constraint:
            ind_p_name = parameter_name.index(p_name)
            matrix_q[i_constraint, ind_p_name] = coeff
    return matrix_q


def rhochi_one_iteration_by_dictionary(
        global_dict, dict_in_out: dict = None, flag_use_precalculated_data: bool = False):
    chi_sq_sum, n_point_sum, der_chi_sq_sum, dder_chi_sq_sum, parameter_name_sum = \
        rhochi_calc_chi_sq_by_dictionary(
            global_dict, dict_in_out=dict_in_out, flag_use_precalculated_data=flag_use_precalculated_data,
            flag_calc_analytical_derivatives=True)

    delta_p = -1.* numpy.linalg.inv(dder_chi_sq_sum).dot(der_chi_sq_sum)
    if "linear_constraints" in global_dict.keys():
        linear_constraints = global_dict["linear_constraints"]
        matrix_q = form_matrix_q(linear_constraints, parameter_name_sum)
        shift_p = calc_shift_p_by_constraints_hamilton(delta_p, dder_chi_sq_sum, matrix_q)
        delta_p += shift_p
    for way, d_p in zip(parameter_name_sum, delta_p):
        global_dict[way[0]][way[1]][way[2]] += d_p
    return chi_sq_sum, n_point_sum, delta_p, parameter_name_sum, der_chi_sq_sum, dder_chi_sq_sum


class _RefinementStopped(Exception):
    """Raised inside the minimizer when the progress callback asks to stop."""


def form_parameter_name(way: tuple, parameter_paths: dict = None) -> str:
    """Path of the parameter for the log, or its way
    ('crystal_fe, atom_occupancy, (0,)') when the path is unknown."""
    way_sh = (way[0], way[1], tuple(int(hh) for hh in way[2]))
    if parameter_paths is not None and way_sh in parameter_paths:
        return parameter_paths[way_sh]
    return f"{way_sh[0]:}, {way_sh[1]:}, {way_sh[2]:}"


def rhochi_rietveld_refinement_by_dictionary(
        global_dict: dict, method: str = "BFGS", callback: Callable = None,
        *, progress_callback: Callable = None, parameter_paths: dict = None):
    """Refine the parameters of global_dict by minimization of chi_sq.

    Parameters
    ----------
    global_dict : dict
        Dictionary of the global object (``GlobalN.get_dictionary()``),
        refined in place.
    method : str
        Minimization method of ``scipy.optimize.minimize`` or
        ``"basinhopping [T=..., niter=...]"``.
    callback : callable, optional
        Passed to ``scipy.optimize.minimize`` unchanged, called as
        ``callback(xk)``.
    progress_callback : callable, optional
        Called after every iteration of the minimizer as
        ``progress_callback(info)``. ``info`` is a dict with keys
        ``"iteration"`` (int, from 1), ``"chi_sq"`` (chi_sq of the current
        parameters, punishment function included), ``"n_point"``,
        ``"chi_sq_per_n"``, ``"parameter_names"`` (dictionary ways of the
        free parameters) and ``"parameter_values"`` (numpy array).
        Returning ``False`` stops the refinement; the parameters of the
        last completed iteration are kept and ``res["stopped"]`` is True.
        Any other return value (``None`` too) continues the refinement.
    parameter_paths : dict, optional
        Readable names of the parameters in the log: way -> path
        (``cryspy.procedure_assistant.parameters.take_dict_paths``).

    Returns
    -------
    chi_sq, parameter_names, dict_in_out, res
        ``res`` is the result of the minimizer (``scipy.optimize.OptimizeResult``)
        or an empty dict when nothing is refined.
    """
    logger.info("*********************************************")
    logger.info("Rietveld refinement by CrysPy (module RhoChi)")
    logger.info("*********************************************\n")
    logger.info("Derivatives are calculated numerically.")
    dict_in_out = {}
    flag_use_precalculated_data = False
    flag_calc_analytical_derivatives = False
    logger.info("Preliminary calculations...")
    chi_sq, n_point, der_chi_sq, dder_chi_sq, parameter_names = rhochi_calc_chi_sq_by_dictionary(
        global_dict,
        dict_in_out=dict_in_out,
        flag_use_precalculated_data=flag_use_precalculated_data, flag_calc_analytical_derivatives=flag_calc_analytical_derivatives)


    # equal marks and the constraint loop: the dependent parameters are
    # calculated from the free ones
    parameter_names_free, parameter_names_fixed, matrix_m, vector_c = \
        split_parameters_by_constraints(global_dict, parameter_names)
    param_0 = [global_dict[way[0]][way[1]][way[2]] for way  in parameter_names_free]
    logger.info(f"Started chi_sq per number of points is {chi_sq/n_point:.2f}.")
    if len(param_0) == 0:
        res = {}
        warnings.warn(
            "Nothing is refined: mark the parameters to refine with parentheses "
            "after the value, for example 1.23()", stacklevel=2)
        return chi_sq, parameter_names, dict_in_out, res
    logger.info(f"Number of fitting parameters {len(param_0):}")
    for name, val in zip(parameter_names_free, param_0):
        logger.info(f" - {form_parameter_name(name, parameter_paths):} {val:.5f}")

    if len(parameter_names_fixed) >0:
        logger.info("Number of constrained parameters:")
    for name in parameter_names_fixed:
        logger.info(f" - {form_parameter_name(name, parameter_paths):}")


    flag_use_precalculated_data = True
    # chi_sq of the recent evaluations, to report the accepted point to
    # progress_callback without recalculating it
    d_recent_chi_sq = {}
    def tempfunc(l_param):
        for way, param in zip(parameter_names_free, l_param):
            global_dict[way[0]][way[1]][way[2]] = param
        set_dependent_parameters(global_dict, parameter_names_fixed, matrix_m, vector_c, l_param)

        chi_sq = rhochi_calc_chi_sq_by_dictionary(
            global_dict,
            dict_in_out=dict_in_out,
            flag_use_precalculated_data=flag_use_precalculated_data,
            flag_calc_analytical_derivatives=flag_calc_analytical_derivatives)[0]
        
        chi_sq_punishement = calc_punishement_function(global_dict)
        if progress_callback is not None:
            d_recent_chi_sq[tuple(numpy.asarray(l_param, dtype=float))] = chi_sq + chi_sq_punishement
        return chi_sq + chi_sq_punishement

    d_progress = {"iteration": 0, "x": numpy.array(param_0, dtype=float), "chi_sq": chi_sq}
    def callback_minimizer(xk, *arg):
        if callback is not None:
            callback(xk, *arg)
        if progress_callback is None:
            return
        x = numpy.array(xk, dtype=float)
        d_progress["iteration"] += 1
        d_progress["x"] = x
        d_progress["chi_sq"] = d_recent_chi_sq.get(tuple(x), d_progress["chi_sq"])
        d_recent_chi_sq.clear()
        info = {
            "iteration": d_progress["iteration"],
            "chi_sq": float(d_progress["chi_sq"]),
            "n_point": int(n_point),
            "chi_sq_per_n": float(d_progress["chi_sq"] / n_point),
            "parameter_names": parameter_names_free,
            "parameter_values": x.copy()}
        if progress_callback(info) is False:
            raise _RefinementStopped()


    logger.info("\nMinimization procedure of chi_sq is running... ")
    #if linear_constraint is not None:
    #    method="SLSQP"
    #    print(f"User constraints are used in the refinement by method '{method:}'.")
    #    res = scipy.optimize.minimize(tempfunc, param_0, method=method, constraints=linear_constraint)
    #    hess_inv, np_first_der = estimate_inversed_hessian_matrix(tempfunc, res.x)
    #    if hess_inv is not None:
    #        res["hess_inv"] = hess_inv
    #else:
    if method.lower().startswith("basinhopping"):
        logger.info("Basinhopping method is used for minimization.")
        niter=10
        T = 10.
        if len(method)>len('basinhopping '):
            for hh in re.split(r"[ ,)]", method[len('basinhopping '):]):
                if hh.lower().startswith("t="):
                    T = float(hh[len('t='):])
                elif hh.lower().startswith("niter="):
                    niter = int(hh[len('niter='):])
        try:
            res = scipy.optimize.basinhopping(
                tempfunc, param_0, niter=niter, T=T, stepsize=0.1, interval=20, disp=True,
                minimizer_kwargs={"callback": callback_minimizer})
            res["stopped"] = False
        except _RefinementStopped:
            res = _result_of_stopped_refinement(d_progress)
        tempfunc(res.x)
    else:
        logger.info(f"{method:} method is used for minimization.")
        try:
            res = scipy.optimize.minimize(tempfunc, param_0, method=method, callback=callback_minimizer)
            res["stopped"] = False
        except _RefinementStopped:
            res = _result_of_stopped_refinement(d_progress)
            # the last evaluation may be a finite-difference step
            tempfunc(res.x)
    if res["stopped"]:
        logger.info(f"Refinement is stopped by the callback after {d_progress['iteration']:} iterations.")
    else:
        logger.info("Optimization is done.")

    logger.info("Calculations for optimal parameters... ")
    dict_in_out = {}
    flag_use_precalculated_data = False
    chi_sq, n_point = rhochi_calc_chi_sq_by_dictionary(
        global_dict,
        dict_in_out=dict_in_out,
        flag_use_precalculated_data=flag_use_precalculated_data, flag_calc_analytical_derivatives=flag_calc_analytical_derivatives)[:2]
    logger.info(f"Optimal chi_sq per n is {chi_sq/n_point:.2f}")

    return chi_sq, parameter_names, dict_in_out, res

def _result_of_stopped_refinement(d_progress: dict):
    return scipy.optimize.OptimizeResult(
        x=d_progress["x"], fun=d_progress["chi_sq"], nit=d_progress["iteration"],
        success=False, stopped=True, message="Refinement is stopped by the callback.")


def rhochi_lsq_by_dictionary(global_dict):
    logger.info("*********************************************")
    logger.info("Rietveld refinement by CrysPy (module RhoChi)")
    logger.info("*********************************************\n")
    logger.info("Derivatives are calculated analytically.")
    logger.info("User constraints are working.")

    dict_in_out = {}
    dict_in_out_2 = {}
    flag_use_precalculated_data = False
    chi_sq_sum, n_point, delta_p, parameter_name_sum, der_chi_sq_sum, dder_chi_sq_sum = rhochi_one_iteration_by_dictionary(
        global_dict, dict_in_out=dict_in_out, flag_use_precalculated_data=flag_use_precalculated_data)
    flag = True
    flag_use_precalculated_data = True
    while flag:
        chi_sq_sum_2, n_point, delta_p_2, parameter_name_sum_2, der_chi_sq_sum_2, dder_chi_sq_sum_2 = rhochi_one_iteration_by_dictionary(
            global_dict, dict_in_out=dict_in_out_2, flag_use_precalculated_data=flag_use_precalculated_data)
        if chi_sq_sum_2 < chi_sq_sum:
            chi_sq_sum = chi_sq_sum_2
            delta_p = delta_p_2
            der_chi_sq_sum = der_chi_sq_sum_2
            dder_chi_sq_sum = dder_chi_sq_sum_2
            parameter_name_sum = parameter_name_sum_2
            dict_in_out = dict_in_out_2
        else:
            for way, d_p in zip(parameter_name_sum_2, delta_p_2):
                global_dict[way[0]][way[1]][way[2]] -= d_p
            for way, d_p in zip(parameter_name_sum, delta_p):
                global_dict[way[0]][way[1]][way[2]] -= d_p
            flag = False
    return chi_sq_sum, n_point, delta_p, parameter_name_sum, der_chi_sq_sum, dder_chi_sq_sum, dict_in_out


def rhochi_calc_chi_sq_by_dictionary(
        global_dict, dict_in_out:dict=None, flag_use_precalculated_data: bool=False,
        flag_calc_analytical_derivatives: bool = False):
    """Calculate chi_sq.
    """
    dict_in_out_keys = dict_in_out.keys()
    dict_keys = global_dict.keys()
    l_dict_crystal, l_dict_magcrystal = [], []
    l_dict_diffrn = []
    l_dict_pd, l_dict_pd2d, l_dict_tof = [], [], []
    for name_key in dict_keys:
        if name_key.startswith("crystal_"):
            l_dict_crystal.append((name_key, global_dict[name_key]))
        # if name_key.startswith("magcrystal_"):
        #     l_dict_magcrystal.append((name_key, global_dict[name_key]))
        if name_key.startswith("diffrn_"):
            l_dict_diffrn.append((name_key, global_dict[name_key]))
        if name_key.startswith("pd_"):
            l_dict_pd.append((name_key, global_dict[name_key]))
        if name_key.startswith("pd2d_"):
            l_dict_pd2d.append((name_key, global_dict[name_key]))
        if name_key.startswith("tof_"):
            l_dict_tof.append((name_key, global_dict[name_key]))

    l_chi_sq, l_der_chi_sq, l_dder_chi_sq = [], [], []
    l_n_point = []
    l_parameter_name, parameter_name_full = [], []
    dder = {}
    l_experiments, l_diff_chi = [], []
    for name_key_diffrn, dict_diffrn in l_dict_diffrn:
        if flag_use_precalculated_data and (name_key_diffrn in dict_in_out_keys):
            dict_in_out_diffrn = dict_in_out[name_key_diffrn]
        else:
            dict_in_out_diffrn = {}
            dict_in_out[name_key_diffrn] = dict_in_out_diffrn

        phase_label = dict_diffrn["phase_name"][0].lower()
        label_1 = f"crystal_{phase_label:}"
        label_2 = f"magcrystal_{phase_label:}"
        if label_1 in dict_keys:
            name_key_crystal = label_1
            dict_crystal = global_dict[label_1]
        elif label_2 in dict_keys:
            name_key_crystal = label_2
            dict_crystal = global_dict[label_2]
        elif (phase_label == "" and
                (len(l_dict_crystal)+len(l_dict_magcrystal)) > 0):
            name_key_crystal, dict_crystal = (l_dict_crystal+l_dict_magcrystal)[0]
        else:
            raise RefinementError(f"Phase {phase_label:} is not found")

        chi_sq, n_point, der_chi_sq, dder_chi_sq, parameter_name = calc_chi_sq_for_diffrn_by_dictionary(
            dict_diffrn, dict_crystal, dict_in_out=dict_in_out_diffrn,
            flag_use_precalculated_data=flag_use_precalculated_data,
            flag_calc_analytical_derivatives=flag_calc_analytical_derivatives)
        l_chi_sq.append(chi_sq)
        l_n_point.append(n_point)
        l_der_chi_sq.append(der_chi_sq)
        l_dder_chi_sq.append(dder_chi_sq)
        l_parameter_name.append(parameter_name)
        parameter_name_full.extend(parameter_name)
        l_experiments.append(name_key_diffrn)

    dict_crystals = [hh[1] for hh in l_dict_crystal]
    for name_key_exp, dict_exp in l_dict_pd + l_dict_pd2d + l_dict_tof:
        if flag_use_precalculated_data and (name_key_exp in dict_in_out_keys):
            dict_in_out_diffrn = dict_in_out[name_key_exp]
        else:
            dict_in_out_diffrn = {}
            dict_in_out[name_key_exp] = dict_in_out_diffrn

        if name_key_exp.startswith("pd_"):
            chi_sq, n_point, der_chi_sq, dder_chi_sq, parameter_name = \
                calc_chi_sq_for_pd_by_dictionary(
                    dict_exp, dict_crystals,
                    dict_in_out=dict_in_out_diffrn,
                    flag_use_precalculated_data=flag_use_precalculated_data,
                    flag_calc_analytical_derivatives=flag_calc_analytical_derivatives)
            l_chi_sq.append(chi_sq)
            l_n_point.append(n_point)
            l_der_chi_sq.append(der_chi_sq)
            l_dder_chi_sq.append(dder_chi_sq)
            l_parameter_name.append(parameter_name)
            parameter_name_full.extend(parameter_name)
            l_experiments.append(name_key_exp)

        elif name_key_exp.startswith("tof_"):
            chi_sq, n_point, der_chi_sq, dder_chi_sq, parameter_name = \
                calc_chi_sq_for_tof_by_dictionary(
                    dict_exp, dict_crystals,
                    dict_in_out=dict_in_out_diffrn,
                    flag_use_precalculated_data=flag_use_precalculated_data,
                    flag_calc_analytical_derivatives=flag_calc_analytical_derivatives)
            l_chi_sq.append(chi_sq)
            l_n_point.append(n_point)
            l_der_chi_sq.append(der_chi_sq)
            l_dder_chi_sq.append(dder_chi_sq)
            l_parameter_name.append(parameter_name)
            parameter_name_full.extend(parameter_name)
            l_experiments.append(name_key_exp)

        elif name_key_exp.startswith("pd2d_"):
            chi_sq, n_point, der_chi_sq, dder_chi_sq, parameter_name = \
                calc_chi_sq_for_pd2d_by_dictionary(
                    dict_exp, dict_crystals,
                    dict_in_out=dict_in_out_diffrn,
                    flag_use_precalculated_data=flag_use_precalculated_data,
                    flag_calc_analytical_derivatives=flag_calc_analytical_derivatives)

            l_chi_sq.append(chi_sq)
            l_n_point.append(n_point)
            l_der_chi_sq.append(der_chi_sq)
            l_dder_chi_sq.append(dder_chi_sq)
            l_parameter_name.append(parameter_name)
            parameter_name_full.extend(parameter_name)
            l_experiments.append(name_key_exp)


    chi_sq_sum = sum(l_chi_sq) #Unity weighting scheme
    n_point_sum = sum(l_n_point)

    # dict keeps the first-seen order; set order depends on PYTHONHASHSEED
    # and made the refinement path differ from one Python run to the next
    parameter_name_sum = list(dict.fromkeys(parameter_name_full))
    der_chi_sq_sum = numpy.zeros((len(parameter_name_sum),), dtype=float)
    dder_chi_sq_sum = numpy.zeros((len(parameter_name_sum), len(parameter_name_sum)), dtype=float)
    if flag_calc_analytical_derivatives:
        for l_pn, dc, ddc in zip(l_parameter_name, l_der_chi_sq, l_dder_chi_sq):
            for i_pn_1, pn_1 in enumerate(l_pn):
                ind_1 = parameter_name_sum.index(pn_1)
                der_chi_sq_sum[ind_1] += dc[i_pn_1]
                dder_chi_sq_sum[ind_1, ind_1] += ddc[i_pn_1, i_pn_1]
                if len(l_pn) > (i_pn_1+1):
                    for i_pn_2, pn_2 in enumerate(l_pn[i_pn_1+1:]):
                        ind_2 = parameter_name_sum.index(pn_2)
                        dder_chi_sq_sum[ind_1, ind_2] += ddc[i_pn_1, i_pn_2+i_pn_1+1]
                        dder_chi_sq_sum[ind_2, ind_1] += ddc[i_pn_2+i_pn_1+1, i_pn_1]
    return chi_sq_sum, n_point_sum, der_chi_sq_sum, dder_chi_sq_sum, parameter_name_sum
