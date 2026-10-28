"""Agreement factors between experimental and calculated profiles.

Used by the reports of the powder experiments (1d, 2d, TOF).
"""
import numpy


def calc_agreement_factors(
        signal_exp: numpy.ndarray, signal_sigma: numpy.ndarray,
        signal_calc: numpy.ndarray, flag_in: numpy.ndarray = None) -> dict:
    """Calculate chi_sq, R-factors of a profile.

    Parameters
    ----------
    signal_exp, signal_sigma, signal_calc : numpy.ndarray
        Experimental signal, its sigma and calculated signal (any shape,
        the same for all three, arbitrary units). Points with not finite
        values are skipped (not measured points of 2d patterns).
    flag_in : numpy.ndarray, optional
        True for the points taken into account (not excluded).

    Returns
    -------
    dict
        "n_point", "chi_sq" (sum of ((exp - calc)/sigma)^2),
        "chi_sq_per_n", "r_p" (sum|exp - calc|/sum|exp|),
        "r_wp" (sqrt(sum w (exp - calc)^2 / sum w exp^2), w = 1/sigma^2).
    """
    signal_exp = numpy.asarray(signal_exp, dtype=float).flatten()
    signal_sigma = numpy.asarray(signal_sigma, dtype=float).flatten()
    signal_calc = numpy.asarray(signal_calc, dtype=float).flatten()
    flag = (numpy.isfinite(signal_exp) & numpy.isfinite(signal_sigma) &
            numpy.isfinite(signal_calc) & (signal_sigma > 0.))
    if flag_in is not None:
        flag = flag & numpy.asarray(flag_in, dtype=bool).flatten()
    n_point = int(flag.sum())
    res = {"n_point": n_point, "chi_sq": numpy.nan, "chi_sq_per_n": numpy.nan,
           "r_p": numpy.nan, "r_wp": numpy.nan}
    if n_point == 0:
        return res
    y_exp, y_sigma, y_calc = signal_exp[flag], signal_sigma[flag], signal_calc[flag]
    weight = 1./numpy.square(y_sigma)
    chi_sq = float(numpy.sum(numpy.square((y_exp - y_calc)/y_sigma)))
    res["chi_sq"] = chi_sq
    res["chi_sq_per_n"] = chi_sq/n_point
    denom_p = numpy.sum(numpy.abs(y_exp))
    if denom_p > 0.:
        res["r_p"] = float(numpy.sum(numpy.abs(y_exp - y_calc))/denom_p)
    denom_wp = numpy.sum(weight*numpy.square(y_exp))
    if denom_wp > 0.:
        res["r_wp"] = float(numpy.sqrt(
            numpy.sum(weight*numpy.square(y_exp - y_calc))/denom_wp))
    return res


def report_powder_agreement(
        dict_signal: dict, flag_sum: bool = True, flag_diff: bool = None) -> str:
    """Report (markdown) on agreement of powder profiles with the model.

    Parameters
    ----------
    dict_signal : dict
        Arrays of the processed profile named as attributes of PdProcL,
        TOFProcL, Pd2dProc: "intensity_plus_net", "intensity_minus_net",
        "intensity_bkg_calc" (model), "intensity", "intensity_sigma"
        (unpolarized experiment) or "intensity_plus", "intensity_plus_sigma",
        "intensity_minus", "intensity_minus_sigma" (polarized experiment),
        "excluded" (excluded points, optional).
    flag_sum : bool
        Sum of up and down signals is fitted.
    flag_diff : bool
        Difference of up and down signals is fitted (polarized experiment);
        None means "if the experiment is polarized".

    The sum signal is compared to plus_net + minus_net + bkg_calc without
    the excluded points; the difference to plus_net - minus_net at all
    points, as chi_sq is calculated in the refinement.
    """
    keys = dict_signal.keys()
    if not (("intensity_plus_net" in keys) and ("intensity_minus_net" in keys)
            and ("intensity_bkg_calc" in keys)):
        return ""
    plus_net = numpy.asarray(dict_signal["intensity_plus_net"], dtype=float)
    minus_net = numpy.asarray(dict_signal["intensity_minus_net"], dtype=float)
    bkg = numpy.asarray(dict_signal["intensity_bkg_calc"], dtype=float)
    flag_polarized = all([key in keys for key in (
        "intensity_plus", "intensity_plus_sigma", "intensity_minus",
        "intensity_minus_sigma")])
    if flag_polarized:
        up = numpy.asarray(dict_signal["intensity_plus"], dtype=float)
        s_up = numpy.asarray(dict_signal["intensity_plus_sigma"], dtype=float)
        down = numpy.asarray(dict_signal["intensity_minus"], dtype=float)
        s_down = numpy.asarray(dict_signal["intensity_minus_sigma"], dtype=float)
        signal_sum = up + down
        sigma_sum = numpy.sqrt(numpy.square(s_up) + numpy.square(s_down))
    elif ("intensity" in keys) and ("intensity_sigma" in keys):
        signal_sum = numpy.asarray(dict_signal["intensity"], dtype=float)
        sigma_sum = numpy.asarray(dict_signal["intensity_sigma"], dtype=float)
    else:
        return ""
    if flag_diff is None:
        flag_diff = flag_polarized
    flag_diff = flag_diff and flag_polarized

    if "excluded" in keys:
        flag_excluded = numpy.asarray(dict_signal["excluded"], dtype=bool)
    else:
        flag_excluded = numpy.zeros(signal_sum.shape, dtype=bool)
    flag_measured = numpy.isfinite(signal_sum)
    n_measured = int(flag_measured.sum())
    n_excluded = int((flag_excluded & flag_measured).sum())

    ls_out = ["# Agreement with the model", ""]
    ls_out.append(f"| Measured points | {n_measured:} |")
    ls_out.append(f"| Excluded points | {n_excluded:} |")
    l_signal = []
    res = calc_agreement_factors(
        signal_sum, sigma_sum, plus_net + minus_net + bkg,
        flag_in=numpy.logical_not(flag_excluded))
    s_name = "up + down" if flag_polarized else "intensity"
    s_fitted = "yes" if (flag_sum or not flag_polarized) else "no"
    l_signal.append((s_name, s_fitted, res))
    if flag_polarized:
        res = calc_agreement_factors(
            up - down, sigma_sum, plus_net - minus_net)
        s_fitted = "yes" if flag_diff else "no"
        l_signal.append(("up − down", s_fitted, res))
    ls_out.append("")
    ls_out.append("| Signal | Fitted | Points | χ²/n | R_p | R_wp |")
    ls_out.append("|:---|:---:|---:|---:|---:|---:|")
    for s_name, s_fitted, res in l_signal:
        ls_out.append(
            f"| {s_name:} | {s_fitted:} | {res['n_point']:} | "
            f"{res['chi_sq_per_n']:.2f} | {100*res['r_p']:.2f}% | "
            f"{100*res['r_wp']:.2f}% |")
    return "\n".join(ls_out)
