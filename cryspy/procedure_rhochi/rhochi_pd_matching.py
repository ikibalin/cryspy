"""Matching mode of a 1D constant-wavelength powder pattern.

In the matching mode the integrated intensities of the reflections are not
calculated from the crystal structure: the pattern is the sum of the peaks of
all reflections of every phase with free integrated intensities

    y(2θ) = bkg(2θ) + Σ_hkl I_hkl P(2θ - 2θ_hkl) L(2θ) / L(2θ_hkl),

where P is the pseudo-Voigt profile of the phase, L the Lorentz factor and
I_hkl the integrated intensity written to the loop pd_peak (the same value
as ``intensity_plus + intensity_minus`` of a Rietveld calculation, so a
pattern calculated from pd_peak is the Rietveld pattern).

- ``refine_intensities``: for the current peak positions and shapes the
  model is linear in I_hkl, so the intensities minimizing chi_sq with
  I_hkl >= 0 are found by the non-negative linear least squares (Pawley fit
  with the positivity constraint; Lawson & Hanson NNLS, scipy.optimize.nnls)
  at every calculation. The cell, the zero shift, the resolution and the
  other parameters of the peak positions and shapes are refined with the
  optimal intensities (variable projection). Polarized data: the intensities
  I+ and I- (up and down) are fitted to the sum and the difference of the
  pattern, each >= 0.
- ``fixed_intensities``: I_hkl are taken from pd_peak and kept fixed.

Reflections at the same position in one phase (2θ equal within 1e-9 rad)
cannot be separated: they get one intensity, shared in proportion to their
multiplicity (Le Bail equipartition). Reflections without measured points
get zero intensity.
"""
import numpy
import scipy.optimize

from cryspy.A_functions_base.unit_cell import calc_sthovl_by_unit_cell_parameters
from cryspy.A_functions_base.powder_diffraction_const_wavelength import (
    calc_profile_pseudo_voight,
    calc_lorentz_factor,
)

na = numpy.newaxis

MATCHING_MODES = ("off", "refine_intensities", "fixed_intensities")


def take_matching_mode(dict_pd: dict) -> str:
    """Matching mode of the experiment ("off" when not given)."""
    matching_mode = dict_pd.get("matching_mode", "off")
    if matching_mode not in MATCHING_MODES:
        raise ValueError(
            f"Unknown matching mode '{matching_mode:}': use one of "
            f"{', '.join(MATCHING_MODES):}.")
    return matching_mode


def take_fixed_peaks(dict_pd: dict, p_name: str) -> dict:
    """Reflections and fixed intensities of a phase from the loop pd_peak."""
    dict_peaks = dict_pd.get("matching_peaks", {}).get(p_name)
    if dict_peaks is None:
        raise ValueError(
            f"Matching mode 'fixed_intensities': the experiment '{dict_pd['name']:}' "
            f"has no loop pd_peak with intensities for the phase '{p_name:}'. "
            "Run the matching mode 'refine_intensities' (or a calculation) first, "
            "it writes the loop.")
    return dict_peaks


def calc_matching_phase(
        matching_mode: str, dict_pd: dict, dict_in_out_phase: dict, p_name: str,
        index_hkl: numpy.ndarray, multiplicity_hkl: numpy.ndarray,
        unit_cell_parameters: numpy.ndarray, wavelength: float,
        ttheta_zs: numpy.ndarray, lorentz_factor: numpy.ndarray,
        resolution_parameters: numpy.ndarray, p_ig: float,
        asymmetry_parameters: tuple, k: float, cthm: float) -> dict:
    """Peaks of one phase in the matching mode.

    Parameters
    ----------
    index_hkl, multiplicity_hkl : numpy.ndarray
        Reflections of the phase, shapes (3, n_hkl) and (n_hkl,); ignored for
        "fixed_intensities", where the reflections of pd_peak are used.
    unit_cell_parameters : numpy.ndarray
        a, b, c (Å), alpha, beta, gamma (radians), shape (6,).
    wavelength : float
        Wavelength (Å).
    ttheta_zs : numpy.ndarray
        2θ of the points corrected for the offsets (radians), shape (n_point,).
    lorentz_factor : numpy.ndarray
        Lorentz factor at ttheta_zs, shape (n_point,).
    resolution_parameters : numpy.ndarray
        U, V, W, X, Y of the phase (instrument + phase), shape (5,).
    asymmetry_parameters : tuple
        p1, p2, p3, p4.

    Returns
    -------
    dict
        "column": peaks of unit integrated intensity, shape (n_point, n_hkl);
        "multiplicity", "ttheta_hkl" (radians, without the offsets),
        "dict_in_out_phase"; for "fixed_intensities" also "iint_plus",
        "iint_minus" (n_hkl,).
    """
    # the dictionary of the phase keeps only the matching results (no
    # structure factors of an earlier Rietveld calculation)
    dict_in_out_phase.clear()
    res = {}
    if matching_mode == "fixed_intensities":
        dict_peaks = take_fixed_peaks(dict_pd, p_name)
        index_hkl = dict_peaks["index_hkl"]
        multiplicity_hkl = dict_peaks["multiplicity"]
        res["iint_plus"] = dict_peaks["intensity_plus"]
        res["iint_minus"] = dict_peaks["intensity_minus"]

    sthovl_hkl = calc_sthovl_by_unit_cell_parameters(index_hkl, unit_cell_parameters)[0]
    ttheta_hkl = 2 * numpy.arcsin(numpy.clip(sthovl_hkl * wavelength, 0., 1.))
    offset_ttheta = dict_pd["offset_ttheta"]
    offset_sycos = dict_pd.get("offset_sycos", 0.0)
    offset_sysin = dict_pd.get("offset_sysin", 0.0)
    ttheta_hkl_shift = (
        offset_ttheta
        + numpy.radians(offset_sycos) * numpy.cos(ttheta_hkl)
        + numpy.radians(offset_sysin) * numpy.sin(ttheta_hkl))

    u, v, w, x, y = resolution_parameters
    p_1, p_2, p_3, p_4 = asymmetry_parameters
    profile_pv = calc_profile_pseudo_voight(
        ttheta_zs, ttheta_hkl, u, v, w, p_ig, x, y, p_1, p_2, p_3, p_4,
        cutoff_fwhm=dict_pd.get("profile_cutoff_fwhm", 0.))[0]
    lf_hkl = calc_lorentz_factor(ttheta_hkl, k=k, cthm=cthm, flag_ttheta=None)[0]
    # 2θ beyond the reflection limit (sinθ/λ > 1/λ) has no peak
    flag_reached = sthovl_hkl * wavelength < 1.
    column = numpy.where(
        flag_reached[na, :],
        lorentz_factor[:, na] * profile_pv / numpy.where(flag_reached, lf_hkl, 1.)[na, :],
        0.)

    dict_in_out_phase["index_hkl"] = index_hkl
    dict_in_out_phase["multiplicity_hkl"] = multiplicity_hkl
    dict_in_out_phase["ttheta_hkl"] = ttheta_hkl + ttheta_hkl_shift
    dict_in_out_phase["profile_pv"] = profile_pv
    res.update({
        "column": column, "multiplicity": multiplicity_hkl, "ttheta_hkl": ttheta_hkl,
        "dict_in_out_phase": dict_in_out_phase})
    return res


def group_equal_positions(l_phase: list):
    """Columns of the reflections that cannot be separated merged.

    Returns the merged columns (n_point, n_group), the group of every
    reflection (concatenated over the phases) and the share of the group
    intensity of every reflection (in proportion to its multiplicity).
    """
    l_key = []
    for i_phase, phase in enumerate(l_phase):
        ttheta_hkl = numpy.round(phase["ttheta_hkl"] / 1e-9)
        l_key.append(numpy.stack([numpy.full(ttheta_hkl.shape, i_phase), ttheta_hkl], axis=1))
    key = numpy.concatenate(l_key, axis=0)
    column = numpy.concatenate([phase["column"] for phase in l_phase], axis=1)
    multiplicity = numpy.concatenate(
        [numpy.asarray(phase["multiplicity"], dtype=float) for phase in l_phase])
    multiplicity = numpy.where(multiplicity > 0, multiplicity, 1.)
    if key.shape[0] == 0:
        return column, numpy.zeros((0,), dtype=int), multiplicity
    key_unique, index_first, index_group = numpy.unique(
        key, axis=0, return_index=True, return_inverse=True)
    index_group = numpy.ravel(index_group)
    multiplicity_group = numpy.bincount(index_group, weights=multiplicity,
                                        minlength=key_unique.shape[0])
    share = multiplicity / multiplicity_group[index_group]
    return column[:, index_first], index_group, share


def solve_nonnegative(a: numpy.ndarray, b: numpy.ndarray) -> numpy.ndarray:
    """x >= 0 minimizing |a x - b|^2 (columns of a without data get 0)."""
    x = numpy.zeros((a.shape[1],), dtype=float)
    norm = numpy.sqrt(numpy.square(a).sum(axis=0))
    if norm.size == 0 or a.shape[0] == 0:
        return x
    active = norm > 1e-12 * norm.max()
    if not numpy.any(active):
        return x
    # columns scaled to the unit norm: the bounded solver converges on
    # intensities differing by orders of magnitude
    a_active = a[:, active] / norm[na, active]
    x_active = scipy.optimize.nnls(a_active, b)[0]
    # the intensities of the passive set solved exactly: the iterations of
    # NNLS stop at a tolerance, which would make chi_sq noisy for the
    # numerical derivatives of the refined parameters
    passive = x_active > 0.
    if numpy.any(passive):
        x_passive = numpy.linalg.lstsq(a_active[:, passive], b, rcond=None)[0]
        if numpy.all(x_passive >= 0.):
            x_active[passive] = x_passive
    x[active] = x_active / norm[active]
    return x


def fit_matching_intensities(
        l_phase: list, signal_sum: numpy.ndarray, sigma_sum: numpy.ndarray,
        in_points: numpy.ndarray, signal_diff: numpy.ndarray = None,
        sigma_diff: numpy.ndarray = None):
    """Integrated intensities >= 0 of the reflections of all phases.

    Parameters
    ----------
    l_phase : list
        Results of calc_matching_phase for every phase.
    signal_sum, sigma_sum : numpy.ndarray
        Measured (up + down) intensity minus the background and its sigma,
        shape (n_point,); None when the sum is not fitted (polarized data
        with chi2 sum False).
    in_points : numpy.ndarray
        Points of the sum taken into account (not excluded), shape (n_point,).
    signal_diff, sigma_diff : numpy.ndarray, optional
        Measured up - down intensity and its sigma for polarized data
        (chi2 diff True); then I+ and I- are fitted separately.

    Writes "iint_plus", "iint_minus" (n_hkl,) to every phase of l_phase.
    """
    column, index_group, share = group_equal_positions(l_phase)
    l_a, l_b = [], []
    n_group = column.shape[1]
    flag_diff = signal_diff is not None
    if signal_sum is not None:
        rows = numpy.logical_and(in_points, sigma_sum > 0)
        a_sum = column[rows] / sigma_sum[rows, na]
        l_a.append(numpy.concatenate([a_sum, a_sum], axis=1) if flag_diff else a_sum)
        l_b.append(signal_sum[rows] / sigma_sum[rows])
    elif not flag_diff:
        raise ValueError("Matching mode: nothing to fit (chi2 sum and diff are False).")
    if flag_diff:
        if signal_sum is None:
            raise ValueError(
                "Matching mode 'refine_intensities' needs the sum of the polarized "
                "pattern (chi2 sum True): the difference alone does not define "
                "the intensities.")
        rows = sigma_diff > 0
        a_diff = column[rows] / sigma_diff[rows, na]
        l_a.append(numpy.concatenate([a_diff, -a_diff], axis=1))
        l_b.append(signal_diff[rows] / sigma_diff[rows])
    x = solve_nonnegative(numpy.concatenate(l_a, axis=0), numpy.concatenate(l_b, axis=0))
    if flag_diff:
        iint_plus = x[:n_group][index_group] * share
        iint_minus = x[n_group:][index_group] * share
    else:
        # unpolarized (or the sum only): up and down are equal halves
        iint_plus = 0.5 * x[index_group] * share
        iint_minus = iint_plus
    i_start = 0
    for phase in l_phase:
        n_hkl = phase["column"].shape[1]
        phase["iint_plus"] = iint_plus[i_start:i_start + n_hkl]
        phase["iint_minus"] = iint_minus[i_start:i_start + n_hkl]
        i_start += n_hkl


def calc_matching_signal(phase: dict):
    """Up and down signals of a phase from its intensities; the intensities
    are written to dict_in_out_phase as in a Rietveld calculation (pd_peak)."""
    column = phase["column"]
    signal_plus = column @ phase["iint_plus"]
    signal_minus = column @ phase["iint_minus"]
    dict_in_out_phase = phase["dict_in_out_phase"]
    dict_in_out_phase["iint_plus_with_factors"] = phase["iint_plus"]
    dict_in_out_phase["iint_minus_with_factors"] = phase["iint_minus"]
    dict_in_out_phase["signal_plus"] = signal_plus
    dict_in_out_phase["signal_minus"] = signal_minus
    return signal_plus, signal_minus
