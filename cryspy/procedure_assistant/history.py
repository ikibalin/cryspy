"""Results of the refinements of an object kept for AI assistants.

An assistant (a language model) refines step by step and decides itself
whether a step is good. CrysPy gives it the facts: the change of
chi-square since the previous calculation, warnings about unphysical
values, and the way back. Before every refinement the parameters are kept,
so that undo_refinement returns them:

    >>> keep_snapshot(cryspy_object)
    >>> result = rhochi_rietveld_refinement(cryspy_object)
    >>> restore_snapshot(cryspy_object)  # the refinement made it worse

The history of an object lives while the object lives (weak references),
nothing is written to the object itself.
"""
import math
import weakref
from typing import List, Optional

from cryspy.procedure_assistant.parameters import take_parameter_records

MAX_SNAPSHOTS = 10
CHI_SQ_TOLERANCE = 1e-3  # relative increase of chi-square taken as "not worse"
CORRELATION_WARNING = 0.95
OCCUPANCY_MAX = 1.05  # partial occupancies may slightly exceed 1 within sigma

D_HISTORY = weakref.WeakKeyDictionary()


class History:
    """chi-square of the last calculation, the refinement flags at that
    moment (checkpoint) and the snapshots before the refinements."""

    def __init__(self):
        self.chi_sq: Optional[float] = None
        self.l_flag_checkpoint: Optional[list] = None
        self.l_snapshot: List[dict] = []


def take_history(cryspy_object) -> History:
    if cryspy_object not in D_HISTORY:
        D_HISTORY[cryspy_object] = History()
    return D_HISTORY[cryspy_object]


def take_flags(cryspy_object) -> list:
    return [(record.path, record.flag_refinement)
            for record in take_parameter_records(cryspy_object)]


def keep_checkpoint(cryspy_object, chi_sq: Optional[float]):
    """Remember chi-square and the refinement flags after a calculation."""
    history = take_history(cryspy_object)
    history.chi_sq = chi_sq
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
        {"values": l_value, "flags": l_flag, "chi_sq": history.chi_sq})
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
    history.chi_sq = snapshot["chi_sq"]
    history.l_flag_checkpoint = snapshot["flags"]
    return snapshot


def take_max_correlation(result: dict) -> Optional[float]:
    matrix = result.get("correlation_matrix", None)
    if matrix is None or len(matrix) < 2:
        return None
    n_param = len(matrix)
    return max(abs(float(matrix[ind_1][ind_2]))
               for ind_1 in range(n_param) for ind_2 in range(ind_1 + 1, n_param))


def form_metrics(result: dict, cryspy_object, n_point: Optional[int] = None,
                 chi_sq_previous: Optional[float] = None) -> dict:
    """Numbers of a calculation or refinement for the assistant and for
    the checks of a recipe."""
    d_metrics = {}
    chi_sq = float(result["chi_sq"]) if "chi_sq" in result else None
    if chi_sq is not None:
        d_metrics["chi_sq"] = chi_sq
    if n_point is None and "n_point" in result:
        n_point = int(result["n_point"])
    if n_point is not None and n_point > 0:
        d_metrics["n_point"] = n_point
        if chi_sq is not None:
            d_metrics["chi_sq_per_n"] = chi_sq / n_point
    if chi_sq_previous is not None and chi_sq is not None:
        d_metrics["chi_sq_previous"] = chi_sq_previous
        if chi_sq_previous > 0:
            d_metrics["chi_sq_change"] = (chi_sq - chi_sq_previous) / chi_sq_previous
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
    change = d_metrics.get("chi_sq_change", None)
    if change is not None and change > CHI_SQ_TOLERANCE:
        l_warning.append(f"chi_sq got worse by {100 * change:.2g} % since the previous calculation.")
    l_adp = [record.path for record in take_parameter_records(cryspy_object)
             if "iso_or_equiv" in record.attribute and record.value < 0.]
    if len(l_adp) > 0:
        l_warning.append("Negative isotropic thermal parameters (unphysical): "
                         + ", ".join(l_adp[:10]) + ".")
    l_occupancy = [record.path for record in take_parameter_records(cryspy_object)
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


def format_metrics(d_metrics: dict) -> str:
    """Comparison with the previous calculation and the warnings as text."""
    ls_out = []
    if "chi_sq_previous" in d_metrics:
        s_change = ""
        if "chi_sq_change" in d_metrics:
            change = d_metrics["chi_sq_change"]
            if change < -CHI_SQ_TOLERANCE:
                s_verdict = "better"
            elif change > CHI_SQ_TOLERANCE:
                s_verdict = "worse"
            else:
                s_verdict = "unchanged"
            s_change = f" ({100 * change:+.3g} %, {s_verdict:})"
        ls_out.append(f"chi_sq since the previous calculation: "
                      f"{d_metrics['chi_sq_previous']:.5g} -> {d_metrics['chi_sq']:.5g}{s_change:}")
    for warning in d_metrics.get("warnings", []):
        ls_out.append(f"Warning: {warning:}")
    if len(d_metrics.get("warnings", [])) > 0:
        ls_out.append("undo_refinement returns the parameters before this refinement.")
    return "\n".join(ls_out)
