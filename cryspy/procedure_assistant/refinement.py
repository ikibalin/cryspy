"""Refinements for AI assistants: correlations at the minimum and the
selection of parameters one by one.

A typical request of a user is "refine the parameters X one by one; keep a
parameter refined when chi_sq gets better and no high correlation appears,
otherwise fix it and go to the next one". A small language model cannot
run such a loop itself, so refine_one_by_one does it:

    >>> result = refine_one_by_one(cryspy_object, "atom_site_moment_fourier")
    >>> [(step["path"], step["decision"]) for step in result["steps"]]
"""
import math
from typing import Callable, List, Optional

import numpy

from cryspy.A_functions_base.significance_test import calc_f_test_added_parameters
from cryspy.C_item_loop_classes.cl_1_inversed_hessian import inversed_hessian_to_correlation
from cryspy.procedure_assistant.parameters import (
    ParameterRecord, find_parameter_records, take_parameter_records)

SIGNIFICANCE = 0.95  # confidence level of the F-test a kept parameter must pass
CORRELATION_LIMIT = 0.9  # |r| of a new pair of refined parameters that rejects a parameter


def take_numerical_correlation(cryspy_object, result: dict) -> dict:
    """The result with the correlations of the inversed Hessian estimated by
    numerical derivatives at the minimum (the sigmas of the refined values
    are from it too); the matrix of BFGS is accumulated along the path of
    the minimizer and may show |r| = 1 for independent parameters."""
    if not cryspy_object.is_attribute("inversed_hessian"):
        return result
    matrix = getattr(cryspy_object.inversed_hessian, "matrix", None)
    if matrix is None or len(matrix) != len(result.get("parameter_name", [])):
        return result
    matrix = numpy.array(matrix, dtype=float)
    # the matrix and the result may list the parameters in other orders
    l_label = [str(label) for label in getattr(cryspy_object.inversed_hessian, "label", [])]
    l_path = list(result.get("parameter_path", []))
    if len(l_path) == len(l_label) and sorted(l_path) == sorted(l_label):
        l_index = [l_label.index(path) for path in l_path]
        matrix = matrix[l_index, :][:, l_index]
    correlation_matrix, _ = inversed_hessian_to_correlation(matrix)
    return dict(result, correlation_matrix=correlation_matrix)


def find_high_correlations(result: dict, limit: float) -> dict:
    """{(path_1, path_2): r} of the pairs of refined parameters with |r| >= limit."""
    l_path = [str(path) for path in result.get("parameter_path", [])]
    matrix = result.get("correlation_matrix", None)
    if matrix is None or len(matrix) != len(l_path):
        return {}
    d_pair = {}
    for ind_1 in range(len(l_path)):
        for ind_2 in range(ind_1 + 1, len(l_path)):
            value = float(matrix[ind_1][ind_2])
            if math.isfinite(value) and abs(value) >= limit:
                d_pair[tuple(sorted((l_path[ind_1], l_path[ind_2])))] = value
    return d_pair


def take_values(cryspy_object) -> list:
    return [(record, record.value, record.sigma) for record in take_parameter_records(cryspy_object)]


def put_values(l_value: list):
    """Return the values and sigmas kept by take_values."""
    for record, value, sigma in l_value:
        if not record.flag_constraint:
            record.set_value(value)
        setattr(record.item, f"{record.attribute:}_sigma", sigma)


def calc_chi_sq_per_n(chi_sq: float, n_point: int) -> float:
    if n_point <= 0 or not math.isfinite(float(chi_sq)):
        return math.inf
    return float(chi_sq) / n_point


def refine_one_by_one(cryspy_object, pattern: str, significance: float = SIGNIFICANCE,
                      correlation_limit: float = CORRELATION_LIMIT,
                      start_value: Optional[float] = None,
                      rejected_value: Optional[float] = None, *,
                      progress: Optional[Callable] = None) -> dict:
    """Refine the parameters matching the pattern one by one and keep only
    those that make the fit better (forward selection).

    The parameters matching the pattern (not constrained by symmetry) are
    fixed first (at rejected_value if it is given); the other parameters
    marked for refinement stay refined in every step. Then, in the order of
    the paths, every candidate is refined together with the parameters kept
    so far. It is kept refined when the decrease of chi_sq is significant
    (Hamilton's F-test, W. C. Hamilton, Acta Cryst. 18 (1965) 502: the
    p-value of the decrease is below 1 - significance) and no new pair of
    refined parameters gets |r| >= correlation_limit; otherwise all values return to the state
    before the step and the candidate is fixed. At the end the kept
    parameters are refined once more, so that sigmas and correlations
    belong to the final model.

    Parameters
    ----------
    cryspy_object : GlobalN
        The object with the crystals and the experiments.
    pattern : str
        Pattern of the paths of the candidates, e.g. 'atom_site_moment_fourier'.
    significance : float
        Confidence level of the F-test that the decrease of chi_sq is not by
        chance (0.95; 0.99 is stricter; 0 keeps every decrease).
    correlation_limit : float
        |r| between two refined parameters that rejects the candidate.
    start_value : float, optional
        Value given before the try to a candidate that is exactly 0 (at 0
        the gradient of a magnetic moment in unpolarized data vanishes).
    rejected_value : float, optional
        Value of a parameter that is not refined, e.g. 0 for magnetic moments:
        all candidates get it first, so that each one is compared with the
        model without it, and the rejected ones keep it. By default the
        candidates keep their values (a value at its minimum then gives no
        decrease of chi_sq).
    progress : callable, keyword-only
        progress(info) after every iteration of the minimizer; info as in
        rhochi_rietveld_refinement with "step", "n_step" and "path" of the
        candidate. Returning False stops the selection.

    Returns
    -------
    dict
        "steps": a dict for every tried candidate with "path", "decision"
        ('kept' or 'fixed'), "reason", "chi_sq_per_n_before",
        "chi_sq_per_n_after", "f_value", "p_value", "value", "sigma", "correlations" ({other
        path: r} with |r| >= correlation_limit); "chi_sq_per_n_start",
        "chi_sq_per_n", "n_point", "result" (the last refinement or
        calculation of the kept model), "stopped".
    """
    from cryspy.procedure_rhochi.rhochi import rhochi_no_refinement, rhochi_rietveld_refinement
    l_candidate: List[ParameterRecord] = [
        record for record in find_parameter_records(cryspy_object, pattern)
        if not record.flag_constraint]
    if len(l_candidate) == 0:
        raise ValueError(f"no free parameters match '{pattern:}'")
    for record in l_candidate:
        record.set_refinement(False)
        if rejected_value is not None:
            # every candidate is compared with the model without it
            record.set_value(rejected_value)
            setattr(record.item, f"{record.attribute:}_sigma", 0.)
    d_stop = {"stopped": False}

    def refine(d_step: dict) -> dict:
        def callback(info: dict):
            if progress is None:
                return None
            flag = progress(dict(info, **d_step))
            if flag is False:
                d_stop["stopped"] = True
            return flag
        result = rhochi_rietveld_refinement(cryspy_object, callback=callback)
        return take_numerical_correlation(cryspy_object, result)

    def refine_or_calculate(i_step: int) -> dict:
        if len(cryspy_object.get_variable_names()) == 0:
            return rhochi_no_refinement(cryspy_object)
        return refine({"step": i_step, "n_step": len(l_candidate), "path": ""})

    n_point = int(rhochi_no_refinement(cryspy_object).get("n_point", 0))
    # the other refined parameters at their minimum: the candidates are compared with it
    result = refine_or_calculate(0)
    chi_sq_per_n_start = calc_chi_sq_per_n(result["chi_sq"], n_point)
    chi_sq_per_n = chi_sq_per_n_start
    d_high = find_high_correlations(result, correlation_limit)
    l_step = []
    flag_last_rejected = False
    for i_step, record in enumerate(l_candidate, start=1):
        if d_stop["stopped"]:
            break
        l_value = take_values(cryspy_object)
        if start_value is not None and record.value == 0.:
            record.set_value(start_value)
        record.set_refinement(True)
        d_out = {"path": record.path, "chi_sq_per_n_before": chi_sq_per_n}
        try:
            result_step = refine({"step": i_step, "n_step": len(l_candidate), "path": record.path})
            chi_sq_per_n_step = calc_chi_sq_per_n(result_step["chi_sq"], n_point)
            d_new = {pair: value for pair, value in
                     find_high_correlations(result_step, correlation_limit).items() if pair not in d_high}
            s_error = ""
        except Exception as error:  # a diverged step rejects the candidate only
            result_step, chi_sq_per_n_step, d_new = {}, math.inf, {}
            s_error = f"{type(error).__name__:}: {error:}"
        if d_stop["stopped"]:
            put_values(l_value)
            record.set_refinement(False)
            break
        d_out["chi_sq_per_n_after"] = chi_sq_per_n_step
        d_out["correlations"] = {
            (pair[1] if pair[0] == record.path else pair[0]): value for pair, value in d_new.items()}
        f_value, p_value = calc_f_test_added_parameters(
            chi_sq_per_n * n_point, chi_sq_per_n_step * n_point, n_point,
            len(result_step.get("parameter_name", [])))
        d_out["f_value"], d_out["p_value"] = f_value, p_value
        s_test = (f"chi_sq/n {chi_sq_per_n:.5g} -> {chi_sq_per_n_step:.5g}, F = {f_value:.3g}, "
                  f"p = {p_value:.2g}")
        flag_keep = False
        if s_error != "":
            d_out["reason"] = f"the refinement failed ({s_error:})"
        elif not math.isfinite(chi_sq_per_n_step):
            d_out["reason"] = "the refinement diverged"
        elif chi_sq_per_n_step >= chi_sq_per_n:
            d_out["reason"] = f"chi_sq/n is not reduced ({s_test:})"
        elif math.isnan(p_value):
            d_out["reason"] = f"no degrees of freedom left: {n_point:} points ({s_test:})"
        elif p_value >= 1. - significance:
            d_out["reason"] = (f"the decrease of chi_sq is not significant at {100 * significance:g} % "
                               f"({s_test:})")
        elif len(d_new) > 0:
            s_pair = ", ".join(f"{path_1:} - {path_2:}: {value:.2f}"
                               for (path_1, path_2), value in d_new.items())
            d_out["reason"] = f"high correlation |r| >= {correlation_limit:g} ({s_pair:})"
        else:
            d_out["reason"] = f"significant decrease of chi_sq ({s_test:})"
            flag_keep = True
        if flag_keep:
            d_out["decision"] = "kept"
            chi_sq_per_n = chi_sq_per_n_step
            d_high.update(find_high_correlations(result_step, correlation_limit))
            result = result_step
            flag_last_rejected = False
        else:
            d_out["decision"] = "fixed"
            put_values(l_value)
            record.set_refinement(False)
            flag_last_rejected = True
        d_out["value"], d_out["sigma"] = record.value, record.sigma
        l_step.append(d_out)
    if flag_last_rejected and not d_stop["stopped"]:
        # the inversed Hessian of the object is of the rejected step
        result = refine_or_calculate(len(l_candidate))
        chi_sq_per_n = calc_chi_sq_per_n(result["chi_sq"], n_point)
    return {"steps": l_step, "chi_sq_per_n_start": chi_sq_per_n_start,
            "chi_sq_per_n": chi_sq_per_n, "n_point": n_point,
            "result": dict(result, n_point=n_point), "stopped": d_stop["stopped"]}
