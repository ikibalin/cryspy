"""Asymmetry of powder peaks from the axial divergence (Finger-Cox-Jephcoat).

The Debye-Scherrer cone of a reflection at 2theta is cut by a detector
(receiving slit) of half height H at the distance L from a sample of half
height S. Rays out of the scattering plane reach the detector at the angle
2phi with cos(2phi) = cos(2theta) sqrt(1 + (h/L)^2), where h is the height
of the point of the cone on the detector: below 90 degrees the peak gets a
tail to low angles, above 90 degrees to high angles, none at 90 degrees.
The observed peak is the symmetric profile averaged over these shifts with
the weight of Finger, Cox & Jephcoat, J. Appl. Cryst. 27 (1994) 892, eq. 4-7:

    W(2phi) ~ A(h) / (h |cos 2phi|),
    A = 2 min(H, S)    for h < |H - S|,
    A = H + S - h      for |H - S| < h < H + S.

In the variable t = h/L the weight is W d(2phi) ~ A(t) dt / ((1 + t^2)
sin 2phi): no singularity at 2phi = 2theta, so a Gauss-Legendre quadrature
in t on [0, |H-S|/L] and [|H-S|/L, (H+S)/L] gives the shifts and weights.
Only the ratios S/L and H/L (FullProf S_L and D_L, GSAS-II S/L and H/L)
enter the formulas.
"""
import numpy

na = numpy.newaxis

N_NODE_SEGMENT = 8


def calc_axial_divergence_shifts(
        ttheta_hkl, s_l: float, d_l: float, n_node_segment: int = N_NODE_SEGMENT):
    """Shifts of the peak positions and their weights for the axial divergence.

    Parameters
    ----------
    ttheta_hkl : numpy.ndarray, shape (n_hkl,)
        Bragg angles 2theta of the reflections, radians.
    s_l : float
        Half height of the sample over the distance sample-detector, S/L
        (dimensionless).
    d_l : float
        Half height of the detector or receiving slit over the same
        distance, H/L (dimensionless).
    n_node_segment : int
        Nodes of the Gauss-Legendre quadrature on each of the two parts of
        the weight.

    Returns
    -------
    shift : numpy.ndarray, shape (n_node, n_hkl)
        2phi - 2theta, radians (negative below 90 degrees).
    weight : numpy.ndarray, shape (n_node, n_hkl)
        Weights summing to 1 for every reflection (a zero shift with the
        weight 1 when s_l or d_l is 0).
    """
    ttheta_hkl = numpy.atleast_1d(numpy.asarray(ttheta_hkl, dtype=float))
    s_l, d_l = abs(float(s_l)), abs(float(d_l))
    n_hkl = ttheta_hkl.size
    if s_l == 0. or d_l == 0.:
        return numpy.zeros((1, n_hkl)), numpy.ones((1, n_hkl))
    x_gl, w_gl = numpy.polynomial.legendre.leggauss(n_node_segment)
    t_inner = abs(s_l - d_l)
    t_outer = s_l + d_l
    l_t, l_w = [], []
    for t_1, t_2 in ((0., t_inner), (t_inner, t_outer)):
        if t_2 > t_1:
            l_t.append(0.5 * (t_2 - t_1) * x_gl + 0.5 * (t_2 + t_1))
            l_w.append(0.5 * (t_2 - t_1) * w_gl)
    t_node = numpy.concatenate(l_t)
    w_node = numpy.concatenate(l_w)
    a_node = numpy.where(t_node < t_inner, 2. * min(s_l, d_l), s_l + d_l - t_node)

    cos_2phi = numpy.cos(ttheta_hkl)[na, :] * numpy.sqrt(1. + numpy.square(t_node))[:, na]
    # a point of the cone beyond 0 or 180 degrees does not reach the detector
    flag_reach = numpy.abs(cos_2phi) < 1.
    ttheta_phi = numpy.arccos(numpy.clip(cos_2phi, -1., 1.))
    sin_2phi = numpy.sin(ttheta_phi)
    with numpy.errstate(divide="ignore", invalid="ignore"):
        weight = numpy.where(
            flag_reach,
            (w_node * a_node / (1. + numpy.square(t_node)))[:, na] / sin_2phi,
            0.)
    weight_sum = weight.sum(axis=0)
    flag_none = weight_sum <= 0.
    weight[:, flag_none] = 0.
    weight[0, flag_none] = 1.
    weight = weight / weight.sum(axis=0)[na, :]
    shift = numpy.where(weight > 0., ttheta_phi - ttheta_hkl[na, :], 0.)
    return shift, weight


def calc_profile_with_axial_divergence(
        profile_function, ttheta_hkl, s_l: float, d_l: float,
        ttheta_max: float = numpy.pi, n_node_segment: int = N_NODE_SEGMENT):
    """Profiles of the reflections averaged over the axial-divergence shifts.

    Parameters
    ----------
    profile_function : callable
        ``profile_function(ttheta_hkl, flags_hkl)`` gives the symmetric
        profiles (n_point, n_hkl_selected) of the reflections selected by
        the boolean array ``flags_hkl`` placed at ``ttheta_hkl`` (radians,
        shape (n_hkl_selected,)).
    ttheta_hkl : numpy.ndarray, shape (n_hkl,)
        Bragg angles 2theta, radians.
    s_l, d_l : float
        S/L and H/L, see :func:`calc_axial_divergence_shifts`.
    ttheta_max : float
        Reflections with 2theta above it (radians) keep the symmetric
        profile (FullProf AsymLim): the effect is small there and every
        node costs one more profile.
    n_node_segment : int
        Nodes of the quadrature on each part of the weight.

    Returns
    -------
    numpy.ndarray, shape (n_point, n_hkl)
        Normalized profiles (the same units as ``profile_function``).
    """
    ttheta_hkl = numpy.atleast_1d(numpy.asarray(ttheta_hkl, dtype=float))
    flags_all = numpy.ones(ttheta_hkl.size, dtype=bool)
    flags_asym = ttheta_hkl <= ttheta_max
    if s_l == 0. or d_l == 0. or not numpy.any(flags_asym):
        return profile_function(ttheta_hkl, flags_all)
    res = numpy.array(profile_function(ttheta_hkl, flags_all), dtype=float, copy=True)
    shift, weight = calc_axial_divergence_shifts(
        ttheta_hkl[flags_asym], s_l, d_l, n_node_segment=n_node_segment)
    profile_asym = 0.
    for shift_node, weight_node in zip(shift, weight):
        profile_asym = profile_asym + weight_node[na, :] * profile_function(
            ttheta_hkl[flags_asym] + shift_node, flags_asym)
    res[:, flags_asym] = profile_asym
    return res


def take_axial_divergence_by_dictionary(dict_pd: dict):
    """(S/L, H/L, 2theta_max in radians) of the experiment or None.

    From the keys ``axial_divergence_parameters`` and
    ``axial_divergence_ttheta_max`` of the dictionary of a powder
    experiment (item pd_instr_axial_divergence).
    """
    if "axial_divergence_parameters" not in dict_pd:
        return None
    s_l, d_l = dict_pd["axial_divergence_parameters"]
    ttheta_max = float(dict_pd.get("axial_divergence_ttheta_max", numpy.pi))
    return (float(s_l), float(d_l), ttheta_max)
