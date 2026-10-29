"""Report of the powder experiments (1d, 2d, TOF) shown in the View tab."""
import numpy

from cryspy.A_functions_base.function_1_agreement_factors import \
    report_powder_agreement
from cryspy.B_parent_classes.cl_3_data import DataN

# attributes of PdProcL, TOFProcL and Pd2dProc used in the report
SIGNAL_NAMES = (
    "intensity_plus_net", "intensity_minus_net", "intensity_bkg_calc",
    "intensity", "intensity_sigma", "intensity_plus", "intensity_plus_sigma",
    "intensity_minus", "intensity_minus_sigma")
EXCLUDED_NAMES = ("excluded", "excluded_points")


def report_powder_experiment(experiment: DataN) -> str:
    """Report (markdown) of a powder experiment.

    Agreement of the processed profile (``pd_proc``, ``tof_proc`` or
    ``pd2d_proc``) with the model for the signals fitted according to
    ``chi2``, the phases and the refinement summary (``refine_ls``) when
    they are given.
    """
    ls_out = []
    flag_sum, flag_diff = True, None
    if experiment.is_attribute("chi2"):
        chi2 = experiment.chi2
        if chi2.is_attribute("sum"):
            flag_sum = bool(chi2.sum)
        if chi2.is_attribute("diff"):
            flag_diff = bool(chi2.diff)
    for name in ("pd_proc", "tof_proc", "pd2d_proc"):
        if not experiment.is_attribute(name):
            continue
        proc = getattr(experiment, name)
        dict_signal = {}
        for signal_name in SIGNAL_NAMES:
            if proc.is_attribute(signal_name):
                dict_signal[signal_name] = numpy.array(
                    getattr(proc, signal_name), dtype=float)
        for excluded_name in EXCLUDED_NAMES:
            if proc.is_attribute(excluded_name):
                dict_signal["excluded"] = numpy.array(
                    getattr(proc, excluded_name), dtype=bool)
        s_report = report_powder_agreement(
            dict_signal, flag_sum=flag_sum, flag_diff=flag_diff)
        if s_report != "":
            ls_out.append(s_report)
    for name in ("phase", "refine_ls"):
        if experiment.is_attribute(name):
            s_report = getattr(experiment, name).report()
            if s_report != "":
                ls_out.append(s_report)
    return "\n\n".join(ls_out)
