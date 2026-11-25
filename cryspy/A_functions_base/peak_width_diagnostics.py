"""Which powder reflections are broader than the instrument (peak width diagnosis).

The net pattern is described by three peaks for every group of reflections
at the same position: one with the width of the instrument (the profile of
the phase) and two broader Lorentzians (FWHM 2.5 and 5 times the
instrumental one), all with free non-negative intensities, plus a smooth
correction of the background (hat functions on 12 nodes). The part of the intensity of a group in the
broad peaks, f, tells whether the group is sharp (f ~ 0) or broad (f ~ 1);
the overlaps are handled because all the groups are fitted together.

If only some reflections are broad, the rule that selects them is looked
for among the shifts t of random stacking faults: a reflection H sees the
faults when H.t is not an integer (t = b/2: k odd; t = (a + b)/2: h + k
odd; t = (a/3, 2a/3): h - k not 3n, ...). The rule with the largest
difference of the mean f of the selected and the other groups is
reported. The same diagnosis is done by eye from a plot by experienced
users (e.g. in layered crystals: hk0 and 00l sharp, hkl broad).
"""
import itertools

import numpy
import scipy.optimize

na = numpy.newaxis

BROAD_FACTORS = (2.5, 5.)  # FWHM of the broad peaks in instrumental FWHM
BACKGROUND_NODES = 12  # nodes of the smooth correction of the background
MIN_SIGNIFICANCE = 5.  # groups weaker than this (in sigma of the intensity) are not used
MIN_SEPARATION = 0.3  # difference of the mean broad parts of the two classes of a rule
MAX_SHARP_PART = 0.35  # mean broad part of the class that stays sharp


def calc_lorentzian_columns(x, position, fwhm):
    """Normalized Lorentzians (1/unit of x) at the points x, shape (n_x, n_peak)."""
    half = 0.5 * fwhm[na, :]
    return half / (numpy.pi * (numpy.square(x[:, na] - position[na, :]) + numpy.square(half)))


def group_reflections(position, fwhm, tolerance: float = 0.3):
    """Groups of reflections closer than tolerance * FWHM (they cannot be separated).

    Parameters
    ----------
    position, fwhm : array (n_hkl,)
        Positions and instrumental FWHM, the same units.

    Returns
    -------
    array (n_hkl,) of int
        Index of the group of every reflection (groups ordered by position).
    """
    order = numpy.argsort(position)
    group = numpy.zeros(position.size, dtype=int)
    i_group = 0
    for i_prev, i_next in zip(order[:-1], order[1:]):
        if position[i_next] - position[i_prev] > tolerance * fwhm[i_prev]:
            i_group += 1
        group[i_next] = i_group
    return group


def calc_background_correction(x, n_node: int = BACKGROUND_NODES):
    """Hat functions on n_node equidistant nodes over x, shape (n_x, n_node).

    A smooth correction of the background: broad peaks must not take the
    errors of the background points, and the correction must not take the
    peaks (the nodes are several peak widths apart).
    """
    node = numpy.linspace(x.min(), x.max(), max(int(n_node), 2))
    step = node[1] - node[0] if node.size > 1 else 1.
    return numpy.maximum(1. - numpy.abs(x[:, na] - node[na, :]) / step, 0.)


def fit_sharp_and_broad(y_net, sigma, column_sharp, column_broad, x):
    """Non-negative intensities of the sharp and broad peaks of every group.

    Parameters
    ----------
    y_net : array (n_x,)
        Measured intensity minus the calculated background.
    sigma : array (n_x,)
        Its sigma (> 0).
    column_sharp : array (n_x, n_group)
        Peaks of the groups with the instrumental profile (unit intensity).
    column_broad : list of arrays (n_x, n_group)
        Broad peaks of the groups (unit intensity).
    x : array (n_x,)
        Positions of the points (for the background correction).

    Returns
    -------
    intensity_sharp : array (n_group,)
    intensity_broad : array (n_group,)
        Sum over the broad peaks.
    sigma_intensity : array (n_group,)
        Rough sigma of the intensity of a group (from the norms of its
        sharp and broad columns, mixed by its broad part).
    chi_sq : float
    """
    n_group = column_sharp.shape[1]
    background = calc_background_correction(x)
    background = numpy.concatenate([background, -background], axis=1)
    matrix = numpy.concatenate([column_sharp] + list(column_broad) + [background], axis=1)
    matrix = matrix / sigma[:, na]
    norm = numpy.sqrt(numpy.square(matrix).sum(axis=0))
    norm = numpy.where(norm > 0., norm, 1.)
    solution, residual = scipy.optimize.nnls(matrix / norm[na, :], y_net / sigma, maxiter=50 * matrix.shape[1])
    solution = solution / norm
    intensity_sharp = solution[:n_group]
    intensity_broad = solution[n_group:n_group * (1 + len(column_broad))].reshape(
        len(column_broad), n_group).sum(axis=0)
    # sigma of a sharp and of a broad intensity, mixed by the broad part
    sigma_sharp = 1. / norm[:n_group]
    sigma_broad = 1. / norm[n_group:n_group * (1 + len(column_broad))].reshape(
        max(len(column_broad), 1), n_group).min(axis=0) if column_broad else sigma_sharp
    total = intensity_sharp + intensity_broad
    part = numpy.where(total > 0., intensity_broad / numpy.where(total > 0., total, 1.), 0.)
    sigma_intensity = (1. - part) * sigma_sharp + part * sigma_broad
    return intensity_sharp, intensity_broad, sigma_intensity, float(residual ** 2)


def form_candidate_shifts():
    """Fault vectors t to test: halves of the cell and thirds in one plane."""
    l_shift = []
    for shift in itertools.product((0., 0.5), repeat=3):
        if any(shift):
            l_shift.append(numpy.array(shift))
    for i_1, i_2 in ((0, 1), (1, 2), (0, 2)):
        shift = numpy.zeros(3)
        shift[i_1], shift[i_2] = 1. / 3., 2. / 3.
        l_shift.append(shift)
    return l_shift


def describe_shift(shift) -> str:
    """The rule of a fault vector in words (k odd, h + k odd, h - k not 3n)."""
    names = "hkl"
    nonzero = [i for i in range(3) if abs(shift[i]) > 1e-9]
    if all(abs(shift[i] - 0.5) < 1e-9 for i in nonzero):
        return " + ".join(names[i] for i in nonzero) + " odd"
    if len(nonzero) == 2:
        i_1, i_2 = nonzero
        sign = "-" if abs(shift[i_1] + shift[i_2] - 1.) < 1e-9 else "+"
        return f"{names[i_1]:} {sign:} {names[i_2]:} not 3n"
    return "H.t not an integer"


def find_fault_rule(l_index_group, broad_part, weight, l_shift=None):
    """The fault vector whose reflections are broad and the others sharp.

    Parameters
    ----------
    l_index_group : list of arrays (3, n_member)
        Reflections of every significant group.
    broad_part : array (n_group,)
        Part of the intensity of the groups in the broad peaks (0 ... 1).
    weight : array (n_group,)
        Weights of the groups.
    l_shift : list of arrays (3,), optional
        Candidates (form_candidate_shifts by default).

    Returns
    -------
    list of dict
        For every candidate with both classes present, sorted by the
        separation (best first): "shift", "rule", "broad_faulted",
        "broad_other" (weighted mean broad parts), "n_faulted", "n_other",
        "separation".
    """
    if l_shift is None:
        l_shift = form_candidate_shifts()
    l_result = []
    for shift in l_shift:
        l_class = []
        for index_group in l_index_group:
            product = shift @ index_group
            flag = numpy.abs(product - numpy.round(product)) > 1e-6
            l_class.append(1 if numpy.all(flag) else (0 if not numpy.any(flag) else -1))
        group_class = numpy.array(l_class)
        flag_f, flag_o = group_class == 1, group_class == 0
        if numpy.count_nonzero(flag_f) < 1 or numpy.count_nonzero(flag_o) < 2:
            continue
        broad_f = float(numpy.average(broad_part[flag_f], weights=weight[flag_f]))
        broad_o = float(numpy.average(broad_part[flag_o], weights=weight[flag_o]))
        l_result.append({
            "shift": shift, "rule": describe_shift(shift), "broad_faulted": broad_f,
            "broad_other": broad_o, "n_faulted": int(numpy.count_nonzero(flag_f)),
            "n_other": int(numpy.count_nonzero(flag_o)), "separation": broad_f - broad_o})
    l_result.sort(key=lambda item: (-round(item["separation"], 3),
                                    numpy.count_nonzero(item["shift"])))
    return l_result
