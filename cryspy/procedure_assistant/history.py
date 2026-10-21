"""Results of the refinements of an object kept for AI assistants.

An assistant (a language model) refines step by step and decides itself
whether a step is good. CrysPy gives it the facts: the change of
chi-square per point (the number of points may change between
calculations) and, for single crystals, of R_F2 since the previous
calculation, warnings about unphysical values, and the way back. Before every refinement the parameters are kept,
so that undo_refinement returns them:

    >>> keep_snapshot(cryspy_object)
    >>> result = rhochi_rietveld_refinement(cryspy_object)
    >>> restore_snapshot(cryspy_object)  # the refinement made it worse

The history of an object lives while the object lives (weak references),
nothing is written to the object itself.
"""
import math
import weakref
from typing import Dict, List, Optional

from cryspy.B_parent_classes.cl_3_data import DataN
from cryspy.procedure_assistant.parameters import take_parameter_records

MAX_SNAPSHOTS = 10
CHI_SQ_TOLERANCE = 1e-3  # relative increase of chi_sq/n or R_F2 taken as "not worse"
CORRELATION_WARNING = 0.95
OCCUPANCY_MAX = 1.05  # partial occupancies may slightly exceed 1 within sigma

D_HISTORY = weakref.WeakKeyDictionary()


class History:
    """Metrics of the last calculation, the refinement flags at that
    moment (checkpoint) and the snapshots before the refinements."""

    def __init__(self):
        self.d_metrics: Optional[dict] = None
        self.l_flag_checkpoint: Optional[list] = None
        self.l_snapshot: List[dict] = []


def take_history(cryspy_object) -> History:
    if cryspy_object not in D_HISTORY:
        D_HISTORY[cryspy_object] = History()
    return D_HISTORY[cryspy_object]


def take_flags(cryspy_object) -> list:
    return [(record.path, record.flag_refinement)
            for record in take_parameter_records(cryspy_object)]


def keep_checkpoint(cryspy_object, d_metrics: dict):
    """Remember the metrics and the refinement flags after a calculation."""
    history = take_history(cryspy_object)
    history.d_metrics = d_metrics
    history.l_flag_checkpoint = take_flags(cryspy_object)


def keep_snapshot(cryspy_object):
    """Keep the values and sigmas before a refinement, with the refinement
    flags of the last checkpoint: undo returns the parameters marked after
    it (by refine_parameters) to their previous state too."""
    history = take_history(cryspy_object)
    l_value = [(record.path, record.value, record.sigma)
               for record in take_parameter_records(cryspy_object)]
    l_flag = history.l_flag_checkpoint
    if l_flag is None:
        l_flag = take_flags(cryspy_object)
    history.l_snapshot.append(
        {"values": l_value, "flags": l_flag, "metrics": history.d_metrics})
    del history.l_snapshot[:-MAX_SNAPSHOTS]


def restore_snapshot(cryspy_object) -> Optional[dict]:
    """Return the parameters to the last snapshot; None if there is none.

    Constrained parameters are restored too, so the object stays
    consistent without applying the constraints again."""
    history = take_history(cryspy_object)
    if len(history.l_snapshot) == 0:
        return None
    snapshot = history.l_snapshot.pop()
    d_record = {record.path: record for record in take_parameter_records(cryspy_object)}
    for path, value, sigma in snapshot["values"]:
        record = d_record.get(path)
        if record is not None:
            record.set_value(value)
            setattr(record.item, f"{record.attribute:}_sigma", sigma)
    for path, flag in snapshot["flags"]:
        record = d_record.get(path)
        if record is not None and not record.flag_constraint:
            record.set_refinement(flag)
    history.d_metrics = snapshot["metrics"]
    history.l_flag_checkpoint = snapshot["flags"]
    return snapshot


def take_max_correlation(result: dict) -> Optional[float]:
    matrix = result.get("correlation_matrix", None)
    if matrix is None or len(matrix) < 2:
        return None
    n_param = len(matrix)
    return max(abs(float(matrix[ind_1][ind_2]))
               for ind_1 in range(n_param) for ind_2 in range(ind_1 + 1, n_param))


def calc_r_f2(diffrn_refln) -> Optional[float]:
    """R_F2 = sum |I_obs - I_calc| / sum I_obs of the measured reflections
    which are not excluded (unpolarized single crystal: I = F^2 with
    extinction and scale); None without measured intensities."""
    numerator, denominator = 0., 0.
    for item in diffrn_refln.items:
        if not (item.is_attribute("intensity") and item.is_attribute("intensity_calc")):
            continue
        if item.is_attribute("excluded") and bool(item.excluded):
            continue
        numerator += abs(float(item.intensity) - float(item.intensity_calc))
        denominator += abs(float(item.intensity))
    if denominator == 0.:
        return None
    return numerator / denominator


def take_r_f2(cryspy_object) -> Dict[str, float]:
    """R_F2 of every single crystal experiment with intensities, by the
    name of its data block."""
    d_r_f2 = {}
    for data in getattr(cryspy_object, "items", []):
        if isinstance(data, DataN) and data.is_attribute("diffrn_refln"):
            r_f2 = calc_r_f2(data.diffrn_refln)
            if r_f2 is not None:
                d_r_f2[str(data.data_name)] = r_f2
    return d_r_f2


def calc_change(value: float, previous: float) -> Optional[float]:
    if previous > 0. and math.isfinite(value):
        return (value - previous) / previous
    return None


def form_metrics(result: dict, cryspy_object, n_point: Optional[int] = None,
                 d_previous: Optional[dict] = None) -> dict:
    """Numbers of a calculation or refinement for the assistant and for
    the checks of a recipe; d_previous: the metrics of the previous
    calculation, compared by chi_sq/n and R_F2."""
    d_metrics = {}
    d_previous = d_previous or {}
    chi_sq = float(result["chi_sq"]) if "chi_sq" in result else None
    if chi_sq is not None:
        d_metrics["chi_sq"] = chi_sq
    if n_point is None and "n_point" in result:
        n_point = int(result["n_point"])
    if n_point is None:  # a refinement without iterations: the same data as before
        n_point = d_previous.get("n_point", None)
    if n_point is not None and n_point > 0:
        d_metrics["n_point"] = n_point
        if chi_sq is not None:
            d_metrics["chi_sq_per_n"] = chi_sq / n_point
    if "chi_sq_per_n" in d_metrics and "chi_sq_per_n" in d_previous:
        d_metrics["chi_sq_per_n_previous"] = d_previous["chi_sq_per_n"]
        change = calc_change(d_metrics["chi_sq_per_n"], d_previous["chi_sq_per_n"])
        if change is not None:
            d_metrics["chi_sq_per_n_change"] = change
    d_r_f2 = take_r_f2(cryspy_object)
    if len(d_r_f2) > 0:
        d_metrics["r_f2"] = d_r_f2
        d_r_f2_previous = {name: value for name, value in d_previous.get("r_f2", {}).items()
                           if name in d_r_f2}
        if len(d_r_f2_previous) > 0:
            d_metrics["r_f2_previous"] = d_r_f2_previous
    if "parameter_name" in result:
        d_metrics["n_refined"] = len(result["parameter_name"])
    max_correlation = take_max_correlation(result)
    if max_correlation is not None:
        d_metrics["max_correlation"] = max_correlation
    if result.get("stopped", False):
        d_metrics["stopped"] = True
    d_metrics["warnings"] = find_warnings(d_metrics, cryspy_object)
    return d_metrics


def find_warnings(d_metrics: dict, cryspy_object) -> List[str]:
    """Signs that the result is wrong or unphysical, as text for the
    assistant; it decides whether to undo the refinement."""
    l_warning = []
    chi_sq = d_metrics.get("chi_sq", None)
    if chi_sq is not None and not math.isfinite(chi_sq):
        l_warning.append("chi_sq is not a finite number: the refinement diverged.")
    change = d_metrics.get("chi_sq_per_n_change", None)
    if change is not None and change > CHI_SQ_TOLERANCE:
        l_warning.append(f"chi_sq/n got worse by {100 * change:.2g} % since the previous calculation.")
    for name, previous in d_metrics.get("r_f2_previous", {}).items():
        change = calc_change(d_metrics["r_f2"][name], previous)
        if change is not None and change > CHI_SQ_TOLERANCE:
            l_warning.append(f"R_F2 of {name:} got worse by {100 * change:.2g} % since the previous calculation.")
    # only the refined parameters: a fixed value is not a result of this refinement
    l_refined = [record for record in take_parameter_records(cryspy_object) if record.flag_refinement]
    l_adp = [record.path for record in l_refined
             if "iso_or_equiv" in record.attribute and record.value < 0.]
    if len(l_adp) > 0:
        l_warning.append("Negative isotropic thermal parameters (unphysical): "
                         + ", ".join(l_adp[:10]) + ".")
    l_occupancy = [record.path for record in l_refined
                   if record.attribute == "occupancy"
                   and not 0. <= record.value <= OCCUPANCY_MAX]
    if len(l_occupancy) > 0:
        l_warning.append("Occupancies outside [0, 1]: " + ", ".join(l_occupancy[:10]) + ".")
    max_correlation = d_metrics.get("max_correlation", None)
    if max_correlation is not None and max_correlation >= CORRELATION_WARNING:
        l_warning.append(f"Strong correlation of refined parameters (|r| = {max_correlation:.3f}).")
    if d_metrics.get("stopped", False):
        l_warning.append("The refinement was stopped before convergence.")
    return l_warning


def format_change(change: Optional[float]) -> str:
    if change is None:
        return ""
    if change < -CHI_SQ_TOLERANCE:
        s_verdict = "better"
    elif change > CHI_SQ_TOLERANCE:
        s_verdict = "worse"
    else:
        s_verdict = "unchanged"
    return f" ({100 * change:+.3g} %, {s_verdict:})"


def format_metrics(d_metrics: dict, flag_undo_hint: bool = True) -> str:
    """Comparison with the previous calculation (chi_sq/n, R_F2 of single
    crystals) and the warnings as text."""
    ls_out = []
    if "chi_sq_per_n_previous" in d_metrics:
        ls_out.append(f"chi_sq/n since the previous calculation: "
                      f"{d_metrics['chi_sq_per_n_previous']:.5g} -> {d_metrics['chi_sq_per_n']:.5g}"
                      + format_change(d_metrics.get("chi_sq_per_n_change", None)))
    d_r_f2_previous = d_metrics.get("r_f2_previous", {})
    for name, r_f2 in d_metrics.get("r_f2", {}).items():
        if name in d_r_f2_previous:
            ls_out.append(f"R_F2 of {name:}: {d_r_f2_previous[name]:.4f} -> {r_f2:.4f}"
                          + format_change(calc_change(r_f2, d_r_f2_previous[name])))
        else:
            ls_out.append(f"R_F2 of {name:} = {r_f2:.4f}")
    for warning in d_metrics.get("warnings", []):
        ls_out.append(f"Warning: {warning:}")
    if flag_undo_hint and len(d_metrics.get("warnings", [])) > 0:
        ls_out.append("undo_refinement returns the parameters before this refinement.")
    return "\n".join(ls_out)
