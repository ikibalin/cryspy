"""Signatures of the misfit of a model to diffraction data.

The residuals of a fit tell what is wrong with the model: calculated peaks
shifted (zero shift, unit cell), too narrow or too broad (profile),
intensities off by a trend with sin(theta)/lambda (displacement
parameters), by the size of the reflection (extinction) or one by one
(structure, texture), a background too low or too high, peaks that the
model has not got at all (an impurity, a wrong cell, magnetic
reflections). The functions find these signatures by linear least
squares on numpy arrays, so that an assistant gets them as facts, each a
line starting with a keyword.

Powder profiles are treated window by window (one window per group of
overlapping calculated peaks): inside a window the residual
``r = y_obs - y_calc`` is fitted by

    r = a n + b dn/dx + c d2n/dx2 + e,

where n is the calculated net intensity (without background) and x the
position (2theta in degrees or time of flight in μs). For a peak
``y_obs = (1 + a) n(x - delta)`` of the width ``s + ds`` one gets
``b = -(1 + a) delta`` and, for a peak of fixed area,
``c = (1 + a) s ds`` (``dn/ds = s d2n/dx2`` for a Gaussian), so the
windows give the shift, the relative change of the width and the
intensity ratio of every peak group. The approach follows the reading of
difference curves of B. H. Toby, Powder Diffr. 21 (2006) 67 and
J. Appl. Cryst. 57 (2024) 175, and L. B. McCusker et al., J. Appl. Cryst.
32 (1999) 36.
"""
import numpy

na = numpy.newaxis

MIN_WINDOW_POINTS = 5
PEAK_LEVEL = 0.1  # a window is where the net intensity exceeds this part of its peak
MIN_PEAK_SIGMA = 3.  # weaker calculated peaks (in sigma of the data) are not analysed
UNEXPLAINED_SIGMA = 5.  # a residual peak above this (in sigma) is reported
SCATTER_LIMIT = 0.15  # weighted rms of log(I_obs/I_calc) of the peak groups
BACKGROUND_SIGMA = 2.  # mean residual (in sigma) of a background region
MIN_SHIFT = 0.05  # smaller shifts (in peak widths) are not reported


def find_peak_windows(net, sigma, level: float = PEAK_LEVEL,
                      min_peak_sigma: float = MIN_PEAK_SIGMA):
    """Windows (start, stop) of the groups of calculated peaks.

    Parameters
    ----------
    net : numpy.ndarray
        (n,) calculated net intensity (without background).
    sigma : numpy.ndarray
        (n,) sigma of the measured intensity.
    level : float
        A window is a run of points with the net intensity above this part
        of the highest point of the run.
    min_peak_sigma : float
        Runs whose highest net intensity is below this many sigma are left out.

    Returns
    -------
    list of (int, int)
        Index ranges ``[start, stop)``.
    """
    flag = net > max(min_peak_sigma * float(numpy.median(sigma)) * level, 0.)
    l_window, start = [], None
    for index, value in enumerate(numpy.append(flag, False)):
        if value and start is None:
            start = index
        elif not value and start is not None:
            l_window.append((start, index))
            start = None
    l_split = []
    for start, stop in l_window:
        l_split.extend(split_window_at_valleys(net, start, stop))
    l_window = l_split
    out = []
    for start, stop in l_window:
        i_max = start + int(numpy.argmax(net[start:stop]))
        if net[i_max] < min_peak_sigma * sigma[i_max] or stop - start < MIN_WINDOW_POINTS:
            continue
        # shrink to the part above level of the highest point, keep a margin
        above = numpy.flatnonzero(net[start:stop] >= level * net[i_max])
        margin = max(2, (above[-1] - above[0]) // 4)
        out.append((max(start, start + above[0] - margin),
                    min(stop, start + above[-1] + margin + 1)))
    return out


VALLEY_LEVEL = 0.5  # a valley deeper than this part of the lower neighbour peak splits a window


def split_window_at_valleys(net, start: int, stop: int, level: float = VALLEY_LEVEL):
    """Split a run of overlapping peaks at the deep valleys of the net
    intensity, so that every window holds one peak or a few close ones."""
    segment = net[start:stop]
    if segment.size < 3:
        return [(start, stop)]
    inner = segment[1:-1]
    i_max = numpy.flatnonzero((inner >= segment[:-2]) & (inner > segment[2:])) + 1
    i_min = numpy.flatnonzero((inner <= segment[:-2]) & (inner < segment[2:])) + 1
    l_cut = []
    for index in i_min:
        left, right = i_max[i_max < index], i_max[i_max > index]
        if left.size == 0 or right.size == 0:
            continue
        lower = min(segment[left[-1]], segment[right[0]])
        if segment[index] < level * lower:
            l_cut.append(index)
    bounds = [0] + l_cut + [segment.size]
    return [(start + lo, start + hi) for lo, hi in zip(bounds[:-1], bounds[1:])]


def calc_window_signature(x, residual, sigma, net):
    """Shift, relative width change and intensity ratio of one peak group.

    Parameters
    ----------
    x : numpy.ndarray
        (m,) positions (2theta in degrees or time of flight).
    residual : numpy.ndarray
        (m,) ``y_obs - y_calc``.
    sigma : numpy.ndarray
        (m,) sigma of y_obs.
    net : numpy.ndarray
        (m,) calculated net intensity.

    Returns
    -------
    dict or None
        ``position`` (centre of mass of net), ``ratio`` (1 + a),
        ``shift`` (delta, units of x: positive when the observed peak is at
        larger x), ``width`` (ds/s), ``width_s`` (s, units of x), ``chi_sq``
        of the window before and ``chi_sq_left`` after the fit, the
        z-values ``z_shift`` and ``z_width``; None when the fit fails.
    """
    weight = 1. / numpy.square(sigma)
    derivative_1 = numpy.gradient(net, x)
    derivative_2 = numpy.gradient(derivative_1, x)
    basis = numpy.stack([net, derivative_1, derivative_2, numpy.ones_like(net)], axis=1)
    normal = basis.transpose() @ (basis * weight[:, na])
    try:
        covariance = numpy.linalg.inv(normal)
    except numpy.linalg.LinAlgError:
        return None
    coefficient = covariance @ (basis.transpose() @ (weight * residual))
    a_value, b_value, c_value, _ = coefficient
    ratio = 1. + a_value
    if ratio <= 0.:
        return None
    total = float(net.sum())
    position = float((net * x).sum() / total)
    width_s = float(numpy.sqrt(max((net * numpy.square(x - position)).sum() / total, 1e-30)))
    error = numpy.sqrt(numpy.diag(covariance))
    chi_sq_left = float((weight * numpy.square(residual - basis @ coefficient)).sum())
    return {"position": position, "ratio": float(ratio),
            "shift": float(-b_value / ratio), "width": float(c_value / (ratio * width_s ** 2)),
            "width_s": width_s, "chi_sq": float((weight * numpy.square(residual)).sum()),
            "chi_sq_left": chi_sq_left,
            "z_shift": float(abs(b_value) / error[1]) if error[1] > 0 else 0.,
            "z_width": float(abs(c_value) / error[2]) if error[2] > 0 else 0.}


def _weighted_line(x, y, weight):
    """Slope and intercept of y = p0 + p1 x and the sigma of p1."""
    basis = numpy.stack([numpy.ones_like(x), x], axis=1)
    normal = basis.transpose() @ (basis * weight[:, na])
    covariance = numpy.linalg.pinv(normal)
    coefficient = covariance @ (basis.transpose() @ (weight * y))
    return coefficient[1], coefficient[0], float(numpy.sqrt(max(covariance[1, 1], 0.)))


def find_unexplained_peaks(x, residual, sigma, net, max_peaks: int = 3,
                           limit: float = UNEXPLAINED_SIGMA):
    """Peaks of the residual where the model has no peak.

    The residual in sigma is smoothed over 5 points; its maxima above
    ``limit`` at points where the calculated net intensity is below 3
    sigma are returned, the strongest first.

    Returns
    -------
    list of (float, float)
        (position, height in sigma).
    """
    z_value = residual / sigma
    kernel = numpy.ones(5) / 5.
    smooth = numpy.convolve(z_value, kernel, mode="same") * numpy.sqrt(5.)
    flag = (smooth > limit) & (net < 3. * sigma)
    l_peak = []
    for index in numpy.flatnonzero(flag):
        lo, hi = max(index - 3, 0), min(index + 4, smooth.size)
        if smooth[index] == smooth[lo:hi].max():
            l_peak.append((float(x[index]), float(smooth[index])))
    l_peak.sort(key=lambda peak: -peak[1])
    out = []
    for position, height in l_peak:
        if all(abs(position - other) > 3. * abs(x[1] - x[0]) for other, _ in out):
            out.append((position, height))
        if len(out) == max_peaks:
            break
    return out


def diagnose_powder_profile(x, y_obs, sigma, y_calc, net, sthovl=None,
                            flag_time: bool = False):
    """Signatures of the misfit of a powder profile.

    Parameters
    ----------
    x : numpy.ndarray
        (n,) 2theta in degrees, or time of flight in μs when ``flag_time``.
    y_obs, sigma, y_calc, net : numpy.ndarray
        (n,) measured intensity, its sigma, calculated intensity with the
        background, calculated net intensity (all peaks, no background).
    sthovl : numpy.ndarray, optional
        (n,) sin(theta)/lambda in Å⁻¹ of the points, for the trend of the
        intensities (2theta or 1/time is used otherwise).
    flag_time : bool
        x is the time of flight.

    Returns
    -------
    list of (str, str)
        (keyword, text): SHIFT, WIDTH, INTENSITY_TREND, INTENSITY_SCATTER,
        SCALE, BACKGROUND, UNEXPLAINED_PEAK, FIT; only the signatures found.
    """
    x, y_obs, sigma = (numpy.asarray(v, dtype=float) for v in (x, y_obs, sigma))
    y_calc, net = numpy.asarray(y_calc, dtype=float), numpy.asarray(net, dtype=float)
    sigma = numpy.where(sigma > 0., sigma, numpy.inf)
    residual = y_obs - y_calc
    unit = "μs" if flag_time else "°"
    l_out = []
    chi_sq_total = float(numpy.square(residual / sigma).sum())
    l_sign, flag_window = [], numpy.zeros(x.size, dtype=bool)
    for start, stop in find_peak_windows(net, sigma):
        flag_window[start:stop] = True
        signature = calc_window_signature(
            x[start:stop], residual[start:stop], sigma[start:stop], net[start:stop])
        if signature is not None:
            l_sign.append(signature)
    chi_sq_peaks = sum(sign["chi_sq"] for sign in l_sign)
    if len(l_sign) >= 2:
        position = numpy.array([sign["position"] for sign in l_sign])
        weight = numpy.array([sign["chi_sq"] - sign["chi_sq_left"] + 1. for sign in l_sign])
        shift = numpy.array([sign["shift"] for sign in l_sign])
        width_s = numpy.array([sign["width_s"] for sign in l_sign])
        z_shift = numpy.array([sign["z_shift"] for sign in l_sign])
        # significant and not negligible against the peak width
        if numpy.median(z_shift) > 3. and numpy.median(numpy.abs(shift) / width_s) > MIN_SHIFT:
            # 2theta: a zero shift is constant, a wrong cell grows as tan(theta);
            # time of flight: zero is constant, a wrong cell (or dtt1) grows as t
            variable = position if flag_time else numpy.tan(numpy.radians(position) / 2.)
            slope, intercept, _ = _weighted_line(variable, shift, weight)
            spread_const = numpy.average(numpy.square(shift - numpy.average(shift, weights=weight)),
                                         weights=weight)
            spread_line = numpy.average(numpy.square(shift - intercept - slope * variable),
                                        weights=weight)
            mean_shift = float(numpy.average(shift, weights=weight))
            kind = ("grows with the position (unit cell" + (" or dtt1)" if flag_time else ")")
                    if spread_line < 0.5 * spread_const else "nearly constant (zero shift)")
            s_far = ""
            if numpy.median(numpy.abs(shift) / width_s) > 1.:
                s_far = (" The shift is larger than the peak width: refine may fail, use "
                         "scan_parameter for the zero shift or the cell first.")
            l_out.append(("SHIFT", f"observed peaks are shifted by {mean_shift:+.3g} {unit:} on "
                          f"average from the calculated ones, {kind:}.{s_far:}"))
        z_width = numpy.array([sign["z_width"] for sign in l_sign])
        width = numpy.array([sign["width"] for sign in l_sign])
        if numpy.median(z_width) > 3. and abs(numpy.average(width, weights=weight)) > 0.03:
            mean_width = float(numpy.average(width, weights=weight))
            half = position < numpy.median(position)
            s_where = ""
            if numpy.any(half) and numpy.any(~half):
                low = numpy.average(width[half], weights=weight[half])
                high = numpy.average(width[~half], weights=weight[~half])
                s_where = (f" ({100 * low:+.0f} % at small, {100 * high:+.0f} % at large "
                           f"{'time' if flag_time else '2theta'})")
            l_out.append(("WIDTH", f"observed peaks are {'broader' if mean_width > 0 else 'narrower'} "
                          f"than calculated by about {100 * abs(mean_width):.0f} %{s_where:}: "
                          "profile parameters."))
        ratio = numpy.array([sign["ratio"] for sign in l_sign])
        log_ratio = numpy.log(ratio)
        mean_log = float(numpy.average(log_ratio, weights=weight))
        if abs(mean_log) > 0.05:
            l_out.append(("SCALE", f"observed peaks are {numpy.exp(mean_log):.3g} times the "
                          "calculated ones on average: scale factor."))
        if sthovl is not None:
            variable = numpy.interp(position, x, sthovl)
            s_variable = "sin(theta)/lambda"
        else:
            variable = 1. / position if flag_time else position
            s_variable = "1/time" if flag_time else "2theta"
        slope, _, slope_sigma = _weighted_line(variable, log_ratio, weight)
        span = float(variable.max() - variable.min())
        if slope_sigma > 0. and abs(slope) > 3. * slope_sigma and abs(slope * span) > 0.1:
            direction = "falls" if slope < 0 else "grows"
            hint = ("calculated intensities fall off too slowly: larger displacement parameters"
                    if slope < 0 else "calculated intensities fall off too fast: smaller "
                    "displacement parameters (or absorption)")
            l_out.append(("INTENSITY_TREND", f"I_obs/I_calc {direction} with {s_variable:} "
                          f"(by a factor {numpy.exp(abs(slope) * span):.2g} over the range): {hint:}."))
        scatter = float(numpy.sqrt(numpy.average(
            numpy.square(log_ratio - mean_log), weights=weight)))
        if scatter > SCATTER_LIMIT:
            l_worst = numpy.argsort(-numpy.abs(log_ratio - mean_log))[:3]
            s_worst = ", ".join(f"{position[i]:.4g}{unit:} ({ratio[i]:.2g})" for i in l_worst)
            l_out.append(("INTENSITY_SCATTER", f"peak groups misfit one by one (rms of "
                          f"log(I_obs/I_calc) = {scatter:.2f}; worst: {s_worst:}): atomic "
                          "positions, occupancies, missing atoms, texture or a wrong magnetic model."))
    flag_bkg = ~flag_window & numpy.isfinite(sigma)
    if numpy.count_nonzero(flag_bkg) >= 10:
        index = numpy.flatnonzero(flag_bkg)
        l_part = numpy.array_split(index, min(5, index.size // 5))
        l_text = []
        for part in l_part:
            mean = float(numpy.mean(residual[part] / sigma[part]) * numpy.sqrt(part.size))
            mean_z = float(numpy.mean(residual[part] / sigma[part]))
            if abs(mean) > BACKGROUND_SIGMA * 3. and abs(mean_z) > 0.5:
                l_text.append(f"{'too low' if mean_z > 0 else 'too high'} at {x[part[0]]:.4g}-"
                              f"{x[part[-1]]:.4g} {unit:} ({mean_z:+.1f} sigma)")
        if l_text:
            l_out.append(("BACKGROUND", "the calculated background is " + "; ".join(l_text)
                          + ": background points or coefficients."))
    l_peak = find_unexplained_peaks(x, residual, sigma, net)
    if l_peak:
        s_peak = ", ".join(f"{position:.4g}{unit:} ({height:.0f} sigma)" for position, height in l_peak)
        l_out.append(("UNEXPLAINED_PEAK", f"measured peaks without calculated ones at {s_peak:}: "
                      "an impurity phase, a wrong cell or space group, or magnetic reflections."))
    share = chi_sq_peaks / chi_sq_total if chi_sq_total > 0 else 0.
    l_out.append(("FIT", f"chi_sq/n = {chi_sq_total / max(x.size, 1):.4g}; {100 * share:.0f} % of chi_sq "
                  f"is in {len(l_sign):} peak groups, the rest in the background regions."))
    return l_out


def diagnose_single_crystal(index_hkl, i_obs, i_sigma, i_calc, sthovl,
                            n_parameter: int = 0, max_worst: int = 5):
    """Signatures of the misfit of single-crystal intensities.

    Parameters
    ----------
    index_hkl : numpy.ndarray
        (3, n) Miller indices.
    i_obs, i_sigma, i_calc : numpy.ndarray
        (n,) measured intensities, their sigmas, calculated intensities.
    sthovl : numpy.ndarray
        (n,) sin(theta)/lambda in Å⁻¹.
    n_parameter : int
        Number of refined parameters (0: not reported).
    max_worst : int
        Number of the worst reflections listed.

    Returns
    -------
    list of (str, str)
        EXTINCTION, INTENSITY_TREND, SCALE, WORST, DATA, FIT.
    """
    i_obs, i_sigma = numpy.asarray(i_obs, dtype=float), numpy.asarray(i_sigma, dtype=float)
    i_calc, sthovl = numpy.asarray(i_calc, dtype=float), numpy.asarray(sthovl, dtype=float)
    flag = (i_sigma > 0.) & (i_calc > 0.) & (i_obs > 2. * i_sigma)
    l_out = []
    z_value = (i_obs - i_calc) / numpy.where(i_sigma > 0., i_sigma, numpy.inf)
    if numpy.count_nonzero(flag) >= 10:
        log_ratio = numpy.log(i_obs[flag] / i_calc[flag])
        # error of log(I_obs) is sigma/I_obs
        weight = numpy.square(i_obs[flag] / i_sigma[flag])
        weight = weight / weight.max()
        mean_log = float(numpy.average(log_ratio, weights=weight))
        order = numpy.argsort(-i_calc[flag])
        n_strong = max(3, order.size // 10)
        strong, rest = order[:n_strong], order[n_strong:]
        ratio_strong = float(numpy.exp(numpy.median(log_ratio[strong])))
        ratio_rest = float(numpy.exp(numpy.median(log_ratio[rest])))
        if ratio_strong < 0.9 * ratio_rest:
            l_out.append(("EXTINCTION", f"the strongest 10 % of reflections have I_obs/I_calc = "
                          f"{ratio_strong:.2f} against {ratio_rest:.2f} for the others: extinction "
                          "(give a start radius or mosaicity before refining it)."))
        elif ratio_strong > 1.1 * ratio_rest:
            l_out.append(("EXTINCTION", f"the strongest reflections are stronger than calculated "
                          f"(I_obs/I_calc = {ratio_strong:.2f} against {ratio_rest:.2f}): the "
                          "extinction correction is too large."))
        if abs(mean_log) > 0.05:
            l_out.append(("SCALE", f"I_obs/I_calc = {numpy.exp(mean_log):.3g} on average: scale factor."))
        slope, _, slope_sigma = _weighted_line(sthovl[flag], log_ratio, weight)
        span = float(sthovl[flag].max() - sthovl[flag].min())
        if slope_sigma > 0. and abs(slope) > 3. * slope_sigma and abs(slope * span) > 0.1:
            hint = ("larger displacement parameters" if slope < 0 else
                    "smaller displacement parameters")
            l_out.append(("INTENSITY_TREND", f"I_obs/I_calc {'falls' if slope < 0 else 'grows'} "
                          f"with sin(theta)/lambda (factor {numpy.exp(abs(slope) * span):.2g} over "
                          f"the range): {hint:}."))
    l_worst = numpy.argsort(-numpy.abs(z_value))[:max_worst]
    s_worst = ", ".join(
        f"({int(index_hkl[0, i]):d} {int(index_hkl[1, i]):d} {int(index_hkl[2, i]):d}) "
        f"{z_value[i]:+.1f}" for i in l_worst if numpy.isfinite(z_value[i]))
    if s_worst:
        l_out.append(("WORST", f"reflections with the largest (I_obs - I_calc)/sigma: {s_worst:}."))
    if n_parameter > 0:
        ratio = i_obs.size / n_parameter
        s_low = " (low: refine fewer parameters)" if ratio < 10 else ""
        l_out.append(("DATA", f"{i_obs.size:} reflections for {n_parameter:} refined parameters "
                      f"({ratio:.0f} per parameter){s_low:}."))
    finite = numpy.isfinite(z_value)
    l_out.append(("FIT", f"chi_sq/n = {float(numpy.mean(numpy.square(z_value[finite]))):.4g}."))
    return l_out


def diagnose_flip_ratios(index_hkl, fr_obs, fr_sigma, fr_calc, sthovl,
                         n_parameter: int = 0, max_worst: int = 5):
    """Signatures of the misfit of flip ratios.

    Parameters
    ----------
    index_hkl : numpy.ndarray
        (3, n) Miller indices.
    fr_obs, fr_sigma, fr_calc : numpy.ndarray
        (n,) measured flip ratios, their sigmas and the calculated ones.
    sthovl : numpy.ndarray
        (n,) sin(theta)/lambda in Å⁻¹.

    Returns
    -------
    list of (str, str)
        NO_MAGNETIC_SIGNAL, TREND, WORST, DATA, FIT.
    """
    fr_obs, fr_sigma = numpy.asarray(fr_obs, dtype=float), numpy.asarray(fr_sigma, dtype=float)
    fr_calc, sthovl = numpy.asarray(fr_calc, dtype=float), numpy.asarray(sthovl, dtype=float)
    z_value = (fr_obs - fr_calc) / numpy.where(fr_sigma > 0., fr_sigma, numpy.inf)
    finite = numpy.isfinite(z_value)
    l_out = []
    if numpy.allclose(fr_calc, 1., atol=1e-6):
        l_out.append(("NO_MAGNETIC_SIGNAL", "all calculated flip ratios are 1: the model has no "
                      "magnetic scattering (field, susceptibility or moments missing; "
                      "check_flip_ratio_model)."))
    if numpy.count_nonzero(finite) >= 10:
        # the magnetic signal is |FR - 1|: a model too weak at large sin(theta)/lambda
        # points to the form factor (kappa)
        signal = numpy.sign(fr_calc[finite] - 1.) * z_value[finite]
        slope, _, slope_sigma = _weighted_line(sthovl[finite], signal, numpy.ones(signal.size))
        if slope_sigma > 0. and abs(slope) > 3. * slope_sigma:
            l_out.append(("TREND", f"the calculated magnetic signal |FR-1| is too "
                          f"{'weak' if slope > 0 else 'strong'} at large sin(theta)/lambda: "
                          "the magnetic form factor (kappa) or the anisotropy of the susceptibility."))
    l_worst = numpy.argsort(-numpy.where(finite, numpy.abs(z_value), -1.))[:max_worst]
    s_worst = ", ".join(
        f"({int(index_hkl[0, i]):d} {int(index_hkl[1, i]):d} {int(index_hkl[2, i]):d}) "
        f"{z_value[i]:+.1f}" for i in l_worst if finite[i])
    if s_worst:
        l_out.append(("WORST", f"flip ratios with the largest (FR_obs - FR_calc)/sigma: {s_worst:}."))
    if n_parameter > 0:
        l_out.append(("DATA", f"{fr_obs.size:} flip ratios for {n_parameter:} refined parameters."))
    l_out.append(("FIT", f"chi_sq/n = {float(numpy.mean(numpy.square(z_value[finite]))):.4g}."))
    return l_out
