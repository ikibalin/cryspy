"""Absorption correction of single crystal experiments: is it needed, and
its set up (categories exptl_absorpt and exptl_crystal of Diffrn).

The estimate compares the model of each experiment with and without the
transmission factors of the crystal (sphere of the given size or the
ellipsoid of exptl_crystal): for intensities after the best scale factor,
since the scale is refined anyway; for flip ratios the transmission
cancels exactly and only the lambda/2 contamination, absorbed differently,
changes them. The change of the model in units of the experimental sigmas,
dchi_sq/n = sum ((s M_A - M) / sigma)^2 / n, says how much chi_sq/n grows
when absorption is ignored. For unpolarized intensities with symmetry
equivalent reflections a test of the data is added (as in PyColl5): does
ln(I / <I>) of the equivalents depend on the directions of the beams as
absorption by a non-spherical crystal would make it (F-test)?
"""
import copy
import logging

import numpy
from scipy import stats

from cryspy.A_functions_base.absorption import (
    calc_absorption_coefficients_by_composition, calc_beam_directions_by_angles,
    calc_beam_directions_normal_beam, calc_transmission_sphere_zero_angle,
    take_element_symbol, take_resonant_elements)
from cryspy.A_functions_base.database import DATABASE
from cryspy.A_functions_base.unit_cell import (
    calc_eq_ccs_by_unit_cell_parameters, calc_sthovl_by_unit_cell_parameters)
from cryspy.B_parent_classes.cl_4_global import GlobalN
from cryspy.C_item_loop_classes.cl_1_exptl_absorpt import ExptlAbsorpt, ExptlCrystal
from cryspy.E_data_classes.cl_1_crystal import Crystal
from cryspy.E_data_classes.cl_2_diffrn import Diffrn
from cryspy.procedure_rhochi.rhochi_diffrn import calc_chi_sq_for_diffrn_by_dictionary

logger = logging.getLogger(__name__)

na = numpy.newaxis

AVOGADRO = 6.02214076e23
DCHI_SQ_NEEDED = 0.5  # dchi_sq/n above which the correction is needed
DCHI_SQ_SMALL = 0.05  # below it absorption is negligible
P_VALUE_SIGNIFICANT = 0.01
N_ANISOTROPY = 5  # parameters of the symmetric traceless matrix Q
HKL_DECIMALS = 3


def take_crystal_of_diffrn(rcif_object: GlobalN, diffrn: Diffrn):
    """Crystal of the phase of the experiment (the first crystal when the
    phase has no label); None when there is none."""
    l_crystal = [item for item in rcif_object.items if isinstance(item, Crystal)]
    label = ""
    if diffrn.is_attribute("phase"):
        label = str(diffrn.phase.label).strip().lower()
    for crystal in l_crystal:
        if crystal.data_name.lower() == label:
            return crystal
    return l_crystal[0] if len(l_crystal) > 0 else None


def take_wavelength(diffrn: Diffrn) -> float:
    return float(diffrn.setup.wavelength)


def calc_absorption_coefficient_of_crystal(crystal: Crystal, wavelength: float) -> dict:
    """mu (1/mm) of the crystal at the wavelength (Ang) and at half of it,
    from the atoms of atom_site (multiplicity x occupancy) and the unit cell;
    the density (g/cm^3) and the resonant absorbers (1/v law is wrong)."""
    # multiplicities are calculated from the space group by the dictionary
    dict_crystal = crystal.get_dictionary()
    atom_site = crystal.atom_site
    d_symbol = {str(item.label): str(item.type_symbol) for item in atom_site.items}
    l_symbol = [d_symbol[str(label)] for label in dict_crystal["atom_label"]]
    numbers = numpy.asarray(dict_crystal["atom_multiplicity"], dtype=float) * \
        numpy.asarray(dict_crystal["atom_occupancy"], dtype=float)
    volume = float(crystal.cell.volume)
    mu_abs, mu_scat = calc_absorption_coefficients_by_composition(
        l_symbol, numbers, volume, wavelength)
    mass = 0.
    for symbol, number in zip(l_symbol, numbers):
        element = take_element_symbol(symbol)
        element_sh = element.lstrip("0123456789")
        weight = DATABASE["Elements"].get(("Atomic weight", "H" if element_sh in ("D", "T") else element_sh))
        if element.startswith(("2H", "D")):
            weight = 2.014
        mass += number * float(weight)
    density = mass / (AVOGADRO * volume * 1e-24)
    return {"mu": mu_abs + mu_scat, "mu_abs": mu_abs, "mu_scat": mu_scat,
            "mu_half": 0.5 * mu_abs + mu_scat, "density": density,
            "resonant": take_resonant_elements(l_symbol)}


def calc_absorption_coefficient(rcif_object: GlobalN, experiment: str = ""):
    """Calculate the linear absorption coefficient mu of single crystal
    experiments from the crystal structure and write it to exptl_absorpt.

    mu = n (sigma_abs + sigma_scat): the absorption cross sections at
    1.798 A scaled as 1/v to the wavelength of the experiment and the bound
    scattering cross sections (Sears 1992), n from the atoms of atom_site
    (multiplicity x occupancy) and the unit cell. The density is written to
    exptl_crystal.density_diffrn.

    Parameters
    ----------
    rcif_object : GlobalN
        The object with the crystals and the single crystal experiments.
    experiment : str
        Name of the experiment; empty for all single crystal experiments.

    Returns
    -------
    str
        mu, its parts and the density of every experiment.
    """
    ls_out = []
    for diffrn in take_diffrn_experiments(rcif_object, experiment):
        crystal = take_crystal_of_diffrn(rcif_object, diffrn)
        if crystal is None:
            ls_out.append(f"{diffrn.data_name:}: no crystal structure, mu is not calculated.")
            continue
        d_mu = calc_absorption_coefficient_of_crystal(crystal, take_wavelength(diffrn))
        set_absorption_coefficient(diffrn, d_mu)
        ls_out.append(format_absorption_coefficient(diffrn, d_mu))
    return "\n".join(ls_out)


def set_absorption_coefficient(diffrn: Diffrn, d_mu: dict):
    if not diffrn.is_attribute("exptl_absorpt"):
        diffrn.add_items([ExptlAbsorpt(correction_type="none")])
    exptl_absorpt = diffrn.exptl_absorpt
    exptl_absorpt.coefficient_mu = round(d_mu["mu"], 6)
    exptl_absorpt.coefficient_mu_half = round(d_mu["mu_half"], 6)
    s_details = ("mu calculated by CrysPy from the composition: sigma_abs (1/v law) + "
                 "bound scattering cross sections, Sears (1992)")
    if d_mu["resonant"]:
        s_details += "; 1/v is approximate for " + ", ".join(d_mu["resonant"])
    exptl_absorpt.process_details = s_details
    if not diffrn.is_attribute("exptl_crystal"):
        diffrn.add_items([ExptlCrystal()])
    diffrn.exptl_crystal.density_diffrn = round(d_mu["density"], 4)


def format_absorption_coefficient(diffrn: Diffrn, d_mu: dict) -> str:
    s_out = (f"{diffrn.data_name:}: mu = {d_mu['mu']:.4g} 1/mm at {take_wavelength(diffrn):.4g} A "
             f"(absorption {d_mu['mu_abs']:.4g}, scattering {d_mu['mu_scat']:.4g} 1/mm; "
             f"{d_mu['mu_half']:.4g} 1/mm at lambda/2), density {d_mu['density']:.4g} g/cm^3.")
    if d_mu["resonant"]:
        s_out += (" Warning: the absorption of " + ", ".join(d_mu["resonant"]) +
                  " does not follow the 1/v law; give the measured or tabulated mu at this "
                  "wavelength in _exptl_absorpt_coefficient_mu.")
    return s_out


def take_diffrn_experiments(rcif_object: GlobalN, experiment: str = "") -> list:
    """Single crystal experiments, or the one with the given name."""
    name = experiment.strip().lower()
    return [item for item in rcif_object.items if isinstance(item, Diffrn)
            and (name == "" or item.data_name.lower() == name)]


def calc_beam_directions_of_dict(dict_diffrn: dict, unit_cell_parameters):
    """Incident and diffracted beams (3, n) of the measured reflections, as
    the absorption correction uses them."""
    index_hkl = dict_diffrn["index_hkl"]
    wavelength = dict_diffrn.get("wavelength_hkl", dict_diffrn["wavelength"])
    sthovl = calc_sthovl_by_unit_cell_parameters(index_hkl, unit_cell_parameters)[0]
    eq_ccs = calc_eq_ccs_by_unit_cell_parameters(index_hkl, unit_cell_parameters)[0]
    q_ccs = 2. * sthovl[na, :] * eq_ccs
    if "angle_omega_chi_phi" in dict_diffrn:
        omega, chi, phi = dict_diffrn["angle_omega_chi_phi"]
        e_in, e_out, misfit = calc_beam_directions_by_angles(
            q_ccs, dict_diffrn["matrix_u"], wavelength, omega, chi, phi)
        return e_in, e_out, misfit < 0.05
    return calc_beam_directions_normal_beam(
        q_ccs, dict_diffrn["matrix_u"][6:9], wavelength,
        str(dict_diffrn.get("absorption_diffracted_beam", "left")))


def calc_model_of_dict(dict_diffrn: dict, dict_crystal: dict):
    """(model values, experimental values, sigmas, included, is_flip_ratio)
    of the experiment dictionary."""
    dict_in_out = {}
    calc_chi_sq_for_diffrn_by_dictionary(dict_diffrn, dict_crystal, dict_in_out=dict_in_out)
    included = numpy.logical_not(dict_diffrn["flip_ratio_excluded"])
    if "flip_ratio_es" in dict_diffrn:
        model = dict_in_out["iint_plus"] / dict_in_out["iint_minus"]
        return model, dict_diffrn["flip_ratio_es"][0], dict_diffrn["flip_ratio_es"][1], included, True
    return (dict_in_out["intensity_calc"], dict_diffrn["intensity_es"][0],
            dict_diffrn["intensity_es"][1], included, False)


def calc_laue_rotations(dict_crystal: dict):
    """Rotations (n, 3, 3) of the Laue class acting on hkl as h' = R^T h."""
    # magnetic space groups (mcif) keep the same rows r_11 ... r_33
    symm_elems = dict_crystal["full_symm_elems"] if "full_symm_elems" in dict_crystal \
        else dict_crystal["full_mcif_elems"]
    r_ij = numpy.asarray(symm_elems[4:13], dtype=float)
    m_r = numpy.moveaxis(r_ij.reshape(3, 3, -1), 2, 0)
    return numpy.unique(numpy.concatenate([m_r, -m_r], axis=0), axis=0)


def group_equivalents(index_hkl, rotations):
    """Index of the group of symmetry equivalent reflections of each hkl."""
    equiv = numpy.round(numpy.einsum("rji,jn->nri", rotations, numpy.asarray(index_hkl, float)),
                        HKL_DECIMALS) + 0.
    d_group, l_index = {}, []
    for set_hkl in equiv:
        key = max(tuple(row) for row in set_hkl)
        l_index.append(d_group.setdefault(key, len(d_group)))
    return numpy.array(l_index, dtype=int)


def _group_means(values, weights, groups):
    n_group = int(groups.max()) + 1 if groups.size else 0
    sum_w = numpy.bincount(groups, weights, n_group)
    sum_wv = numpy.bincount(groups, weights * values, n_group)
    return (sum_wv / sum_w)[groups]


def _direction_products(u):
    x, y, z = u
    return numpy.stack([x * x - z * z, y * y - z * z, 2. * x * y, 2. * x * z, 2. * y * z], axis=1)


def _r_int(intensity, groups):
    mean = _group_means(intensity, numpy.ones_like(intensity), groups)
    return float(numpy.abs(intensity - mean).sum() / intensity.sum())


def test_absorption_anisotropy(intensity, sigma, groups, e_in, e_out) -> dict:
    """F-test: does ln(I / <I>) of equivalent reflections depend on the
    directions of the beams as ln A = -u^T Q u (u_in and u_out, Q symmetric
    traceless) of a non-spherical crystal? Groups of one reflection and
    I <= 0 are left out. Returns {'p_value', 'r_int', 'r_int_corrected',
    'n_used', 'n_groups', 'longest'} or {'error': text}."""
    flag = (intensity > 0) & (sigma > 0)
    count = numpy.bincount(groups[flag], minlength=int(groups.max()) + 1 if groups.size else 0)
    flag &= count[groups] >= 2
    l_group = numpy.unique(groups[flag])
    n_used = int(flag.sum())
    n_free = n_used - len(l_group) - N_ANISOTROPY
    if n_free < 3:
        return {"error": f"not enough equivalent reflections ({n_used} in {len(l_group)} groups "
                         "of two or more) for the test of the data"}
    groups_used = numpy.searchsorted(l_group, groups[flag])
    y = numpy.log(intensity[flag])
    w = numpy.square(intensity[flag] / sigma[flag])
    x = _direction_products(e_in[:, flag]) + _direction_products(e_out[:, flag])
    y_d = y - _group_means(y, w, groups_used)
    x_d = x - numpy.stack([_group_means(col, w, groups_used) for col in x.T], axis=1)
    sqrt_w = numpy.sqrt(w)
    coefficients = numpy.linalg.lstsq(x_d * sqrt_w[:, na], y_d * sqrt_w, rcond=None)[0]
    rss_0 = float((w * y_d**2).sum())
    rss_1 = float((w * (y_d - x_d @ coefficients)**2).sum())
    if rss_1 > 0:
        f_value = (rss_0 - rss_1) / N_ANISOTROPY / (rss_1 / n_free)
        p_value = float(stats.f.sf(f_value, N_ANISOTROPY, n_free))
    else:
        p_value = 0.
    qxx, qyy, qxy, qxz, qyz = -coefficients
    m_q = numpy.array([[qxx, qxy, qxz], [qxy, qyy, qyz], [qxz, qyz, -qxx - qyy]])
    eigenvalues, eigenvectors = numpy.linalg.eigh(m_q)
    corrected = intensity[flag] * numpy.exp(-(x @ coefficients))
    return {"p_value": p_value, "r_int": _r_int(intensity[flag], groups_used),
            "r_int_corrected": _r_int(corrected, groups_used), "n_used": n_used,
            "n_groups": len(l_group), "longest": eigenvectors[:, numpy.argmax(eigenvalues)]}


def form_correction_dict(dict_diffrn: dict, correction_type: str, mu: float, mu_half: float,
                         radius: float = None, sizes=None, euler=(0., 0., 0.)) -> dict:
    dict_out = copy.copy(dict_diffrn)
    for key in list(dict_out.keys()):
        if key.startswith("flags_absorption_"):
            dict_out[key] = numpy.zeros_like(dict_out[key], dtype=bool)
    dict_out["absorption_correction_type"] = correction_type
    dict_out["absorption_mu"] = numpy.array([mu], dtype=float)
    dict_out["absorption_mu_ratio_half"] = mu_half / mu if mu > 0 else 1.
    if radius is not None:
        dict_out["absorption_radius"] = numpy.array([radius], dtype=float)
    if sizes is not None:
        dict_out["absorption_size"] = numpy.asarray(sizes, dtype=float)
        dict_out["absorption_euler"] = numpy.asarray(euler, dtype=float)
    return dict_out


def estimate_absorption_of_diffrn(rcif_object: GlobalN, diffrn: Diffrn,
                                  crystal_size: float = 0., dict_global: dict = None) -> dict:
    """Estimate of the absorption correction of one experiment; see the
    module docs. crystal_size: mean size (diameter of the sphere) in mm, 0 to
    take exptl_crystal. The keys of the result: 'name', 'mu' (dict or None),
    'radius', 'mu_r', 'type' ('flip ratio' or 'intensity'), 'c_lambda2',
    'sphere' and 'ellipsoid' ({'t_min', 't_max', 'dchi_sq'}), 'anisotropy',
    'verdict' ('needed', 'small', 'not needed', 'unknown'), 'lines'."""
    d_res = {"name": diffrn.data_name, "lines": []}
    ls_line = d_res["lines"]
    crystal = take_crystal_of_diffrn(rcif_object, diffrn)
    if dict_global is None:
        dict_global = rcif_object.get_dictionary()
    dict_diffrn = dict_global[diffrn.get_name()]
    d_mu = None
    if crystal is not None:
        try:
            d_mu = calc_absorption_coefficient_of_crystal(crystal, take_wavelength(diffrn))
        except (ValueError, KeyError, TypeError, AttributeError) as error:
            ls_line.append(f"mu is not calculated from the structure: {error:}")
    mu, mu_half = None, None
    if diffrn.is_attribute("exptl_absorpt") and diffrn.exptl_absorpt.is_attribute("coefficient_mu"):
        mu = float(diffrn.exptl_absorpt.coefficient_mu)
        mu_half = float(diffrn.exptl_absorpt.coefficient_mu_half) \
            if diffrn.exptl_absorpt.is_attribute("coefficient_mu_half") else mu
        ls_line.append(f"mu = {mu:.4g} 1/mm (given in exptl_absorpt).")
    elif d_mu is not None:
        mu, mu_half = d_mu["mu"], d_mu["mu_half"]
        ls_line.append(format_absorption_coefficient(diffrn, d_mu).split(": ", 1)[1])
    d_res["mu"] = d_mu
    if d_mu is not None and d_mu["resonant"] and mu == d_mu["mu"]:
        ls_line.append("The 1/v law is approximate for " + ", ".join(d_mu["resonant"]) +
                       ": the estimate is rough.")
    if mu is None:
        d_res["verdict"] = "unknown"
        ls_line.append("No mu: give _exptl_absorpt_coefficient_mu or the crystal structure.")
        return d_res

    exptl_crystal = diffrn.exptl_crystal if diffrn.is_attribute("exptl_crystal") else None
    sizes = exptl_crystal.get_sizes() if exptl_crystal is not None else None
    if crystal_size is not None and crystal_size > 0.:
        radius = 0.5 * float(crystal_size)
    elif exptl_crystal is not None and exptl_crystal.get_radius() is not None:
        radius = exptl_crystal.get_radius()
    else:
        radius = None
    if radius is None:
        d_res["verdict"] = "unknown"
        ls_line.append("The size of the crystal is not given (exptl_crystal.size_rad or "
                       "size_max, size_mid, size_min). Transmission of a sphere at 2theta = 0: " +
                       ", ".join(f"D = {diameter:g} mm: {float(calc_transmission_sphere_zero_angle(mu * diameter / 2.)):.3f}"
                                 for diameter in (0.5, 1., 2., 5.)) + ".")
        return d_res
    d_res["radius"], d_res["mu_r"] = radius, mu * radius
    ls_line.append(f"Sphere of the radius {radius:.3g} mm: mu R = {mu * radius:.3g}, transmission at "
                   f"2theta = 0 is {float(calc_transmission_sphere_zero_angle(mu * radius)):.3f}.")

    unit_cell_parameters = dict_global[crystal.get_name()]["unit_cell_parameters"] \
        if crystal is not None else None
    dict_crystal = dict_global[crystal.get_name()] if crystal is not None else None
    dict_none = copy.copy(dict_diffrn)
    dict_none["absorption_correction_type"] = "none"
    model_0, exp_value, exp_sigma, included, flag_fr = calc_model_of_dict(dict_none, dict_crystal)
    d_res["type"] = "flip ratio" if flag_fr else "intensity"
    sigma_0 = numpy.where(exp_sigma > 0, exp_sigma, numpy.inf)
    d_res["chi_sq_per_n"] = float(numpy.square(
        (exp_value - model_0)[included] / sigma_0[included]).mean())
    s_data = "flip ratios" if flag_fr else "intensities"
    ls_line.append(f"Data: {s_data:}, {int(included.sum()):} reflections; chi_sq/n of the "
                   f"model without the correction is {d_res['chi_sq_per_n']:.4g}.")
    c_lambda2 = float(numpy.atleast_1d(dict_diffrn.get("c_lambda2", numpy.nan))[0])
    d_res["c_lambda2"] = c_lambda2
    flag_c_lambda2 = not numpy.isnan(c_lambda2) and c_lambda2 != 0.

    l_model = [("sphere", form_correction_dict(dict_diffrn, "sphere", mu, mu_half, radius=radius))]
    if sizes is not None and not numpy.allclose(sizes, sizes[0]):
        l_model.append(("ellipsoid", form_correction_dict(
            dict_diffrn, "gaussian", mu, mu_half, sizes=sizes,
            euler=exptl_crystal.get_euler_angles())))
    sigma = numpy.where(exp_sigma > 0, exp_sigma, numpy.inf)
    l_dchi = []
    for name, dict_model in l_model:
        model_a = calc_model_of_dict(dict_model, dict_crystal)[0]
        transmission = calc_chi_sq_transmission(dict_model, unit_cell_parameters)
        if flag_fr:
            scale = 1.
        else:
            w = 1. / numpy.square(sigma[included])
            scale = float((w * model_a[included] * model_0[included]).sum() /
                          (w * numpy.square(model_a[included])).sum())
        dchi_sq = float(numpy.square((scale * model_a - model_0)[included] / sigma[included]).mean())
        d_res[name] = {"t_min": float(transmission.min()), "t_max": float(transmission.max()),
                       "dchi_sq": dchi_sq}
        l_dchi.append(dchi_sq)
        s_shape = "sphere" if name == "sphere" else \
            "ellipsoid " + " x ".join(f"{size:g}" for size in sizes) + " mm"
        ls_line.append(f"{s_shape:}: transmission {transmission.min():.3f} .. {transmission.max():.3f} "
                       f"(T_max/T_min = {transmission.max() / transmission.min():.3f}); ignoring it "
                       f"changes chi_sq/n by {dchi_sq:.3g}.")
    if flag_fr:
        if flag_c_lambda2:
            ls_line.append("Flip ratios: the transmission cancels in I+/I-; only the lambda/2 "
                           f"contamination (c_lambda2 = {c_lambda2:.3g}) is absorbed differently.")
        else:
            ls_line.append("Flip ratios: the transmission is the same for I+ and I- and cancels "
                           "exactly (single scattering, no lambda/2 contamination).")
    else:
        ls_line.append("Intensities: the scale factor is fitted again with the correction; a sphere "
                       "changes intensities with 2theta only and correlates with B_iso.")
        if crystal is not None:
            e_in, e_out, flag_reachable = calc_beam_directions_of_dict(dict_diffrn, unit_cell_parameters)
            groups = group_equivalents(dict_diffrn["index_hkl"], calc_laue_rotations(dict_crystal))
            use = included & flag_reachable
            d_test = test_absorption_anisotropy(
                exp_value[use], exp_sigma[use], groups[use], e_in[:, use], e_out[:, use])
            d_res["anisotropy"] = d_test
            if "error" in d_test:
                ls_line.append("Test of the equivalents: " + d_test["error"] + ".")
            else:
                s_sign = "significant" if d_test["p_value"] < P_VALUE_SIGNIFICANT else "not significant"
                ls_line.append(
                    f"Test of the equivalents ({d_test['n_used']:} reflections in {d_test['n_groups']:} "
                    f"groups): the dependence on the beam directions is {s_sign:} "
                    f"(p = {d_test['p_value']:.2g}); R_int {100. * d_test['r_int']:.2f} % -> "
                    f"{100. * d_test['r_int_corrected']:.2f} % with an empirical correction.")
                if d_test["p_value"] < P_VALUE_SIGNIFICANT and \
                        d_test["r_int_corrected"] < 0.9 * d_test["r_int"]:
                    longest = d_test["longest"]
                    ls_line.append("The crystal is not spherical: use correction_type 'gaussian' with "
                                   "the ellipsoid; the longest path is along (" +
                                   " ".join(f"{x:.2f}" for x in longest) +
                                   ") of the Cartesian system x || a*, z || c.")
    dchi_max = max(l_dchi)
    if dchi_max > DCHI_SQ_NEEDED or (
            "anisotropy" in d_res and d_res["anisotropy"].get("p_value", 1.) < P_VALUE_SIGNIFICANT
            and d_res["anisotropy"]["r_int_corrected"] < 0.9 * d_res["anisotropy"]["r_int"]):
        d_res["verdict"] = "needed"
    elif dchi_max > DCHI_SQ_SMALL:
        d_res["verdict"] = "small"
    else:
        d_res["verdict"] = "not needed"
    return d_res


def calc_chi_sq_transmission(dict_model: dict, unit_cell_parameters):
    from cryspy.A_functions_base.absorption import calc_transmission_by_dictionary
    return calc_transmission_by_dictionary(
        dict_model, dict_model["index_hkl"], unit_cell_parameters,
        dict_model.get("wavelength_hkl", dict_model["wavelength"]))[0]


D_VERDICT = {
    "needed": "absorption correction is NEEDED",
    "small": "absorption is small but visible: the correction is advisable",
    "not needed": "absorption correction is NOT needed",
    "unknown": "cannot be estimated",
}


def estimate_absorption_correction(rcif_object: GlobalN, crystal_size: float = 0.):
    """Estimate whether single crystal data need an absorption correction.

    For every single crystal experiment (unpolarized intensities, flip
    ratios, magnetic satellites, data for MEM) mu is taken from
    exptl_absorpt or calculated from the crystal structure; the model is
    calculated with and without the transmission factors of a sphere (and
    of the ellipsoid of exptl_crystal) and the change is compared with the
    experimental sigmas. Equivalent reflections of unpolarized intensities
    are tested for an anisotropic absorption (non-spherical crystal).
    Nothing in the object is changed.

    Parameters
    ----------
    rcif_object : GlobalN
        The object with the crystals and the single crystal experiments.
    crystal_size : float
        Mean size of the crystal (diameter of the sphere), mm; 0 to take it
        from exptl_crystal.

    Returns
    -------
    str
        Report with the verdict for every experiment.
    """
    dict_global = rcif_object.get_dictionary()
    ls_out = []
    for diffrn in take_diffrn_experiments(rcif_object):
        d_res = estimate_absorption_of_diffrn(rcif_object, diffrn, crystal_size, dict_global)
        ls_out.append(f"# {diffrn.data_name:}: {D_VERDICT[d_res['verdict']]:}")
        ls_out.extend(f"- {line:}" for line in d_res["lines"])
        ls_out.append("")
    if len(ls_out) == 0:
        return "No single crystal experiments."
    return "\n".join(ls_out).strip()


def apply_absorption_correction(rcif_object: GlobalN, correction_type: str = "gaussian",
                                experiment: str = ""):
    """Switch on the absorption correction of single crystal experiments.

    The transmission factors of the crystal multiply the calculated
    intensities of the refinement (RhoChi) and of MEM. mu is calculated from
    the crystal structure when exptl_absorpt has none. The sizes come from
    exptl_crystal: size_rad for a sphere, size_max, size_mid, size_min and
    the Euler angles for an ellipsoid. T_min and T_max are reported.

    Parameters
    ----------
    rcif_object : GlobalN
        The object with the crystals and the single crystal experiments.
    correction_type : str
        'gaussian' (ellipsoid, numerical integration), 'sphere' or 'none'
        (switch off).
    experiment : str
        Name of the experiment; empty for all single crystal experiments.

    Returns
    -------
    str
        What was set and the range of the transmission factors.
    """
    correction_type = correction_type.strip().lower()
    if correction_type not in ("gaussian", "sphere", "none"):
        raise ValueError("The correction type should be 'gaussian', 'sphere' or 'none'.")
    ls_out = []
    l_diffrn = take_diffrn_experiments(rcif_object, experiment)
    if len(l_diffrn) == 0:
        return "No single crystal experiments" + (f" named '{experiment:}'." if experiment else ".")
    for diffrn in l_diffrn:
        if not diffrn.is_attribute("exptl_absorpt"):
            diffrn.add_items([ExptlAbsorpt(correction_type="none")])
        exptl_absorpt = diffrn.exptl_absorpt
        if correction_type == "none":
            exptl_absorpt.correction_type = "none"
            ls_out.append(f"{diffrn.data_name:}: absorption correction is switched off.")
            continue
        exptl_crystal = diffrn.exptl_crystal if diffrn.is_attribute("exptl_crystal") else None
        if exptl_crystal is None or exptl_crystal.get_sizes() is None:
            raise ValueError(f"{diffrn.data_name:}: give the size of the crystal in exptl_crystal "
                             "(size_rad, or size_max, size_mid, size_min in mm).")
        if correction_type == "sphere" and not exptl_crystal.is_attribute("size_rad"):
            exptl_crystal.size_rad = round(exptl_crystal.get_radius(), 4)
        if not exptl_absorpt.is_attribute("coefficient_mu"):
            crystal = take_crystal_of_diffrn(rcif_object, diffrn)
            if crystal is None:
                raise ValueError(f"{diffrn.data_name:}: give _exptl_absorpt_coefficient_mu "
                                 "(no crystal structure to calculate it).")
            d_mu = calc_absorption_coefficient_of_crystal(crystal, take_wavelength(diffrn))
            set_absorption_coefficient(diffrn, d_mu)
            ls_out.append(format_absorption_coefficient(diffrn, d_mu))
        exptl_absorpt.correction_type = correction_type
        dict_global = rcif_object.get_dictionary()
        crystal = take_crystal_of_diffrn(rcif_object, diffrn)
        if crystal is None:
            ls_out.append(f"{diffrn.data_name:}: correction '{correction_type:}' is set.")
            continue
        unit_cell_parameters = dict_global[crystal.get_name()]["unit_cell_parameters"]
        transmission = calc_chi_sq_transmission(dict_global[diffrn.get_name()], unit_cell_parameters)
        exptl_absorpt.correction_t_min = float(transmission.min())
        exptl_absorpt.correction_t_max = float(transmission.max())
        ls_out.append(f"{diffrn.data_name:}: correction '{correction_type:}', mu = "
                      f"{float(exptl_absorpt.coefficient_mu):.4g} 1/mm, transmission "
                      f"T_min = {transmission.min():.4f}, T_max = {transmission.max():.4f}.")
    return "\n".join(ls_out)
