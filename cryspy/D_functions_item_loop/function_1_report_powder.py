"""Report of the powder experiments (1d, 2d, TOF) shown in the View tab."""
from cryspy.B_parent_classes.cl_3_data import DataN


def report_powder_experiment(experiment: DataN) -> str:
    """Report (markdown) of a powder experiment.

    Agreement of the processed profile with the model (proc loop with the
    method ``report_agreement``, fitted signals from ``chi2``), the phases
    and the refinement summary (``refine_ls``) when they are given.
    """
    ls_out = []
    flag_sum, flag_diff = True, None
    if experiment.is_attribute("chi2"):
        chi2 = experiment.chi2
        if chi2.is_attribute("sum"):
            flag_sum = bool(chi2.sum)
        if chi2.is_attribute("diff"):
            flag_diff = bool(chi2.diff)
    for item in experiment.items:
        if hasattr(type(item), "report_agreement"):
            s_report = item.report_agreement(
                flag_sum=flag_sum, flag_diff=flag_diff)
            if s_report != "":
                ls_out.append(s_report)
    for name in ("phase", "refine_ls"):
        if experiment.is_attribute(name):
            s_report = getattr(experiment, name).report()
            if s_report != "":
                ls_out.append(s_report)
    return "\n\n".join(ls_out)
