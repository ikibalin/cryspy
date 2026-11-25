"""Random stacking faults in 1D powder experiments (loop stacking_fault).

A phase whose layers normal to N = (h k l) get at random, with the
probability alpha, an extra in-plane shift t has two kinds of reflections:
those with H.t an integer stay sharp Bragg peaks, the others are spread
along the normal of the layers and give broad or asymmetric lines (a sharp
front and a tail to high angles for the reflections in the plane of the
layers). The loop stacking_fault of a pd experiment gives, for a phase,
the layers, the fault vector t and the refinable probability alpha.
The spread reflections keep their integrated intensity (Rietveld:
structure factors; profile matching: free intensities), their line is the
powder average of the Lorentzian along the normal convolved with the
instrumental profile (resolution and asymmetry p1-p4 of the experiment).

The procedures add the rows and report them; diagnose_peak_widths finds
the groups of reflections broader than the instrument and the fault
vector and the layers that explain them. They are shown in the menus of
CrysPy Editor and used by AI assistants.
"""
import itertools
import logging

import numpy

from cryspy.A_functions_base.powder_diffraction_const_wavelength import (
    calc_h_g, calc_h_l, calc_h_pv, calc_profile_pseudo_voight)
from cryspy.A_functions_base.stacking_faults import (
    calc_fault_beta,
    calc_flag_faulted,
    calc_rod_geometry,
    calc_rod_mass,
    calc_rod_width,
)
from cryspy.A_functions_base.symmetry_elements import calc_equivalent_reflections
from cryspy.B_parent_classes.cl_4_global import GlobalN
from cryspy.C_item_loop_classes.cl_1_stacking_fault import StackingFault, StackingFaultL
from cryspy.E_data_classes.cl_1_crystal import Crystal
from cryspy.E_data_classes.cl_2_pd import Pd
from cryspy.E_data_classes.cl_2_pd2d import Pd2d
from cryspy.E_data_classes.cl_2_tof import TOF

logger = logging.getLogger(__name__)
na = numpy.newaxis

# sampling of the spread reflections before the convolution: points per
# FWHM of the instrument
POINTS_PER_FWHM = 5
SUB_BINS = 4
MAX_GRID_POINTS = 20000


def take_stacking_fault_of_phase(dict_pd: dict, phase_name: str):
    """Row of the phase in the loop stacking_fault of the experiment.

    Returns
    -------
    dict or None
        None when the phase has no row, otherwise "index", "normal_index"
        (3,), "shift" (3,), "probability" (clipped to 0 ... 1).
    """
    if "stacking_fault_phase_name" not in dict_pd:
        return None
    l_index = numpy.flatnonzero(dict_pd["stacking_fault_phase_name"] == phase_name)
    if l_index.size == 0:
        return None
    index = int(l_index[0])
    return {
        "index": index,
        "normal_index": numpy.array(dict_pd["stacking_fault_normal_index"][:, index], dtype=int),
        "shift": numpy.array(dict_pd["stacking_fault_shift"][:, index], dtype=float),
        "probability": float(numpy.clip(dict_pd["stacking_fault_probability"][index], 0., 1.)),
    }


def form_virtual_grid(ttheta_zs, u, v, w, p_ig, x, y):
    """Points (radians) where the spread reflections are sampled.

    The step is a fifth of the narrowest instrumental FWHM, never larger
    than the step of the data; the grid covers the data with a margin of
    5 FWHM.
    """
    h_pv = calc_h_pv(calc_h_g(u, v, w, p_ig, ttheta_zs)[0], calc_h_l(x, y, ttheta_zs)[0])[0]
    h_pv = numpy.abs(h_pv[numpy.isfinite(h_pv)])
    if h_pv.size == 0:
        h_pv = numpy.array([0.1])
    h_min = numpy.radians(max(float(numpy.min(h_pv)), 1e-3))
    step = h_min / POINTS_PER_FWHM
    if ttheta_zs.size > 1:
        step = min(step, float(numpy.median(numpy.abs(numpy.diff(ttheta_zs)))))
    margin = 5. * numpy.radians(float(numpy.max(h_pv)))
    t_min = max(float(ttheta_zs.min()) - margin, step)
    t_max = min(float(ttheta_zs.max()) + margin, numpy.pi - step)
    n_point = min(int(numpy.ceil((t_max - t_min) / step)) + 1, MAX_GRID_POINTS)
    return numpy.linspace(t_min, t_max, n_point)


def take_symmetry_of_crystal(dict_crystal: dict):
    """Reduced symmetry elements and centrosymmetry for the equivalent reflections."""
    if "reduced_symm_elems" in dict_crystal:
        return dict_crystal["reduced_symm_elems"], bool(dict_crystal["centrosymmetry"])
    if "full_mcif_elems" in dict_crystal:
        return dict_crystal["full_mcif_elems"][:13], False
    return numpy.array([[0], [0], [0], [1], [1], [0], [0], [0], [1], [0], [0], [0], [1]]), False


def calc_stacking_fault_mass(
        index_hkl, normal_index, shift, probability: float, unit_cell_parameters,
        reduced_symm_elems, centrosymmetry: bool, s_edges):
    """Spread of the reflections over the bins of s (averaged over the equivalents).

    Parameters
    ----------
    index_hkl : array (3, n_hkl)
        Reflections of the phase (one of every set of equivalents).
    normal_index, shift, probability
        The row of stacking_fault.
    unit_cell_parameters : array (6,)
    reduced_symm_elems, centrosymmetry
        Symmetry of the crystal (calc_equivalent_reflections): the faults
        break it, every equivalent reflection is spread on its own rod.
    s_edges : array (n_bin + 1,)
        Edges of the bins, 1/angstrom.

    Returns
    -------
    flag_faulted : array (n_hkl,) of bool
        Reflections spread by the faults (H.t not an integer for at least
        one equivalent).
    mass : array (n_bin, n_faulted)
        Part of the intensity of each spread reflection in the bins.
    """
    index_eq = calc_equivalent_reflections(index_hkl, reduced_symm_elems, centrosymmetry)
    index_eq = numpy.concatenate([index_eq, -index_eq], axis=2)  # Friedel pairs
    n_hkl, n_eq = index_eq.shape[1], index_eq.shape[2]
    flat = index_eq.reshape(3, n_hkl * n_eq)
    flag_eq = calc_flag_faulted(flat, shift).reshape(n_hkl, n_eq)
    flag_faulted = numpy.any(flag_eq, axis=1)
    n_bin = s_edges.size - 1
    if not numpy.any(flag_faulted):
        return flag_faulted, numpy.zeros((n_bin, 0))
    index_f = index_eq[:, flag_faulted, :].reshape(3, -1)
    g, z, d_layer = calc_rod_geometry(index_f, normal_index, unit_cell_parameters)
    gamma, delta = calc_rod_width(calc_fault_beta(index_f, shift, probability), d_layer)
    # an equivalent with H.t integer is a Bragg peak: its Lorentzian is
    # infinitely narrow, the same formula gives it in the bin of s_H
    z_centre = z + delta
    # equal rods are calculated once
    key = numpy.round(numpy.stack([g, z_centre, gamma], axis=0) * 1e9).astype(numpy.int64)
    key_unique, index_first, index_inverse = numpy.unique(
        key, axis=1, return_index=True, return_inverse=True)
    mass_unique = calc_rod_mass(s_edges, g[index_first], z_centre[index_first], gamma[index_first])
    mass = mass_unique[:, numpy.ravel(index_inverse)].reshape(n_bin, -1, n_eq).mean(axis=2)
    return flag_faulted, mass


def calc_stacking_fault_profile(
        row: dict, dict_crystal: dict, index_hkl, unit_cell_parameters, wavelength,
        ttheta_zs, u, v, w, p_ig, x, y, p_1, p_2, p_3, p_4):
    """Profiles of the spread reflections of a phase at the points.

    The spread reflection is sampled on a fine grid of 2theta (the parts of
    its intensity in the bins of the grid) and convolved with the
    instrumental profile; a reflection without faults would get the
    profile of the phase.

    Returns
    -------
    flag_faulted : array (n_hkl,) of bool
    profile : array (n_point, n_faulted)
        In the units of calc_profile_pseudo_voight (1/degree): it replaces
        the columns of the flagged reflections of the profile of the phase.
    """
    wavelength = float(numpy.atleast_1d(wavelength)[0])
    ttheta_grid = form_virtual_grid(ttheta_zs, u, v, w, p_ig, x, y)
    # the parts of the intensity in SUB_BINS sub-bins of every step are
    # shared linearly between the two nearest points of the grid: a narrow
    # line keeps its position (not moved to the nearest point)
    n_fine = (ttheta_grid.size - 1) * SUB_BINS
    tth_edges = numpy.linspace(ttheta_grid[0], ttheta_grid[-1], n_fine + 1)
    s_edges = 2. * numpy.sin(0.5 * numpy.clip(tth_edges, 0., numpy.pi)) / wavelength
    reduced_symm_elems, centrosymmetry = take_symmetry_of_crystal(dict_crystal)
    flag_faulted, mass_fine = calc_stacking_fault_mass(
        index_hkl, row["normal_index"], row["shift"], row["probability"],
        unit_cell_parameters, reduced_symm_elems, centrosymmetry, s_edges)
    if mass_fine.shape[1] == 0:
        return flag_faulted, numpy.zeros((ttheta_zs.size, 0))
    position = (numpy.arange(n_fine) + 0.5) / SUB_BINS
    index_low = numpy.floor(position).astype(int)
    fraction = (position - index_low)[:, na]
    mass = numpy.zeros((ttheta_grid.size, mass_fine.shape[1]))
    numpy.add.at(mass, index_low, mass_fine * (1. - fraction))
    numpy.add.at(mass, index_low + 1, mass_fine * fraction)
    profile_grid = calc_profile_pseudo_voight(
        ttheta_zs, ttheta_grid, u, v, w, p_ig, x, y, p_1, p_2, p_3, p_4)[0]
    return flag_faulted, profile_grid @ mass


# ---------------------------------------------------------------------------
# Procedures

def take_powder_1d_experiments(rcif_object: GlobalN, experiment: str = "") -> list:
    """Powder experiments at constant wavelength (pd), or the one with the given name."""
    name = experiment.strip().lower()
    return [item for item in rcif_object.items if isinstance(item, Pd)
            and (name == "" or item.data_name.lower() == name)]


def take_crystal_by_label(rcif_object: GlobalN, phase_label: str = ""):
    """Crystal with the data name phase_label (the first crystal for "")."""
    l_crystal = [item for item in rcif_object.items if isinstance(item, Crystal)]
    if phase_label.strip() == "":
        return l_crystal[0] if l_crystal else None
    for item in l_crystal:
        if item.data_name.lower() == phase_label.strip().lower():
            return item
    return None


def parse_triple(text: str, name: str, dtype=float):
    """Three numbers from "a b c" (fractions like 1/2 allowed)."""
    words = text.replace(",", " ").split()
    if len(words) != 3:
        raise ValueError(f"{name:} must be three numbers, e.g. '0 0 1'; got '{text:}'.")
    values = []
    for word in words:
        if "/" in word:
            numerator, denominator = word.split("/")
            values.append(float(numerator) / float(denominator))
        else:
            values.append(float(word))
    if dtype is int:
        if any(abs(value - round(value)) > 1e-9 for value in values):
            raise ValueError(f"{name:} must be integers; got '{text:}'.")
        return numpy.array([int(round(value)) for value in values], dtype=int)
    return numpy.array(values, dtype=float)


def add_stacking_fault(
        rcif_object: GlobalN, normal: str = "0 0 1", shift: str = "0 1/2 0",
        probability: float = 0.1, phase_label: str = "", experiment: str = ""):
    """Add random stacking faults to a powder phase: only the reflections that see the fault shift become broad or asymmetric.

    The layers normal to the reciprocal vector N = (h k l) (distance 1/|N|
    between them; "0 0 2" are two layers per c) get, between two
    neighbour layers with the probability alpha, an extra shift t of all
    the layers above. The reflections with H.t an integer stay sharp, the
    others are spread along the normal of the layers: a sharp front and a
    long tail to high angles for the reflections in the plane of the
    layers, broad peaks for the others. The integrated intensities are
    kept (Rietveld and profile matching). The probability is refined.
    Only the 1D powder experiments at constant wavelength (pd).

    Parameters
    ----------
    rcif_object : GlobalN
        The object with the crystals and the powder experiments.
    normal : str
        Miller indices of the layer planes as "h k l".
    shift : str
        Fault vector t in fractions of the cell as "x y z" (e.g. "0 1/2 0").
    probability : float
        Start value of the fault probability between two layers (0 ... 1).
    phase_label : str
        Label of the phase; empty for every phase of the experiments.
    experiment : str
        Name of the powder experiment; empty for all pd experiments.

    Returns
    -------
    str
        The phases with the faults and the reflections they spread.
    """
    l_experiment = take_powder_1d_experiments(rcif_object, experiment)
    if len(l_experiment) == 0:
        if any(isinstance(item, (Pd2d, TOF)) for item in rcif_object.items):
            raise ValueError("The stacking faults are implemented only for the 1D "
                             "powder experiments at constant wavelength (pd).")
        return "No powder experiments" + (f" named '{experiment:}'." if experiment else ".")
    normal_index = parse_triple(normal, "normal", int)
    if not numpy.any(normal_index):
        raise ValueError("The normal of the layers is (0 0 0).")
    shift_vector = parse_triple(shift, "shift")
    probability = float(probability)
    if not 0. <= probability <= 1.:
        raise ValueError("The probability must be between 0 and 1.")
    if abs(float(normal_index @ shift_vector) - round(float(normal_index @ shift_vector))) > 1e-6:
        logger.warning("The fault shift %s is not in the plane of the layers %s: "
                       "N.t is not an integer.", shift, normal)
    ls_out = []
    for pd in l_experiment:
        l_label = [label for label in pd.phase.label
                   if phase_label.strip() == "" or label.lower() == phase_label.strip().lower()]
        if len(l_label) == 0:
            raise ValueError(f"{pd.data_name:}: no phase '{phase_label:}'.")
        for label in l_label:
            if not pd.is_attribute("stacking_fault"):
                pd.add_items([StackingFaultL()])
            loop = pd.stacking_fault
            row = loop.take_row(label)
            if row is not None:
                loop.items.remove(row)
            row = StackingFault(
                phase_label=label, normal_h=int(normal_index[0]), normal_k=int(normal_index[1]),
                normal_l=int(normal_index[2]), shift_x=float(shift_vector[0]),
                shift_y=float(shift_vector[1]), shift_z=float(shift_vector[2]),
                probability=probability)
            row.probability_refinement = True
            loop.items.append(row)
            ls_out.append(f"{pd.data_name:}, phase {label:}: layers ({normal.strip():}), fault "
                          f"shift ({shift.strip():}), probability {probability:.3f} (refined). "
                          f"Spread reflections: {describe_rule(shift_vector):}.")
    return "\n".join(ls_out)


def describe_rule(shift_vector) -> str:
    """Words for the reflections H with H.t not an integer (k odd, h - k not 3n)."""
    from cryspy.A_functions_base.peak_width_diagnostics import describe_shift
    return describe_shift(numpy.asarray(shift_vector, dtype=float))


def report_stacking_fault(rcif_object: GlobalN, experiment: str = ""):
    """Report the random stacking faults of the powder phases.

    For every row of the loop stacking_fault: the layers, the fault
    vector, the probability with its error, the width of the spread
    reflections along the normal and the correlation length of the
    stacking, and the first spread reflections.

    Parameters
    ----------
    rcif_object : GlobalN
        The object with the crystals and the powder experiments.
    experiment : str
        Name of the powder experiment; empty for all pd experiments.

    Returns
    -------
    str
        Report in markdown.
    """
    ls_out = []
    for pd in take_powder_1d_experiments(rcif_object, experiment):
        if not pd.is_attribute("stacking_fault"):
            continue
        wavelength = float(pd.setup.wavelength)
        ls_out.append(f"# {pd.data_name:}")
        for row in pd.stacking_fault.items:
            crystal = take_crystal_by_label(rcif_object, row.phase_label)
            if crystal is None:
                ls_out.append(f"- {row.phase_label:}: no crystal with this label.")
                continue
            ls_out.append(form_report_of_row(row, crystal, wavelength))
        ls_out.append("")
    if len(ls_out) == 0:
        return "No pd experiment has stacking faults (stacking_fault)."
    return "\n".join(ls_out).strip()


def form_report_of_row(row, crystal, wavelength: float) -> str:
    """Lines of the report for one phase."""
    from cryspy.procedure_rhochi.rhochi_pd import calc_index_hkl_multiplicity_for_phase
    from cryspy.A_functions_base.structure_factor import \
        calc_extinction_rule_by_symmetry_elements
    normal_index, shift_vector = row.get_normal_index(), row.get_shift()
    dict_crystal = crystal.get_dictionary()
    ucp = numpy.array(dict_crystal["unit_cell_parameters"], dtype=float)
    s_sigma = f" ± {row.probability_sigma:.4f}" if row.probability_sigma else ""
    index_hkl = calc_index_hkl_multiplicity_for_phase(
        dict_crystal, 1e-4, 0.999 / wavelength, ucp)[0]
    if "full_symm_elems" in dict_crystal:
        flag_absent = calc_extinction_rule_by_symmetry_elements(
            index_hkl, dict_crystal["full_symm_elems"], dict_crystal["translation_elems"], False)
        index_hkl = index_hkl[:, ~flag_absent]
    reduced_symm_elems, centrosymmetry = take_symmetry_of_crystal(dict_crystal)
    index_eq = calc_equivalent_reflections(index_hkl, reduced_symm_elems, centrosymmetry)
    product = numpy.einsum("i,ihe->he", shift_vector, index_eq)
    flag = numpy.any(numpy.abs(product - numpy.round(product)) > 1e-6, axis=1)
    index_f = index_hkl[:, flag]
    g, z, d_layer = calc_rod_geometry(index_f, normal_index, ucp)
    gamma = calc_rod_width(calc_fault_beta(index_f, shift_vector, float(row.probability)),
                           d_layer)[0]
    s_line = (f"- {row.phase_label:}: layers ({' '.join(str(i) for i in normal_index):}) "
              f"{d_layer:.3f} Å apart, fault shift ({' '.join(f'{v:g}' for v in shift_vector):}), "
              f"probability {row.probability:.4f}{s_sigma}. Spread reflections: "
              f"{describe_rule(shift_vector):}")
    gamma_max = float(numpy.max(gamma)) if gamma.size else 0.
    if gamma_max > 0.:
        s_line += (f"; HWHM along the normal {gamma_max:.4f} 1/Å, correlation length "
                   f"of the stacking about {1. / (2. * numpy.pi * gamma_max):.1f} Å "
                   f"({1. / (2. * numpy.pi * gamma_max * d_layer):.1f} layers)")
    if index_f.shape[1] > 0:
        s_line += ". First spread reflections: " + ", ".join(
            "".join(str(int(v)) for v in index_f[:, i]) for i in range(min(8, index_f.shape[1])))
    return s_line + "."


# ---------------------------------------------------------------------------
# Diagnosis of the peak widths

ALPHA_GRID = (0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45)
MAX_TABLE_ROWS = 30


def take_signal_exp(dict_pd: dict):
    """Measured intensity (sum of up and down) and its sigma at the points."""
    if "signal_exp_plus" in dict_pd and "signal_exp_minus" in dict_pd:
        plus, minus = dict_pd["signal_exp_plus"], dict_pd["signal_exp_minus"]
        return plus[0] + minus[0], numpy.sqrt(numpy.square(plus[1]) + numpy.square(minus[1]))
    return dict_pd["signal_exp"][0], dict_pd["signal_exp"][1]


def form_candidate_normals(shift) -> list:
    """Primitive normals with indices -1 ... 1 (one of +-N) in which the shift lies."""
    l_normal = []
    for normal in itertools.product((-1, 0, 1), repeat=3):
        if normal == (0, 0, 0) or next(v for v in normal if v != 0) < 0:
            continue
        product = float(numpy.array(normal) @ shift)
        if abs(product - round(product)) < 1e-6:
            l_normal.append(numpy.array(normal, dtype=int))
    return l_normal


def merge_columns(column, group, multiplicity):
    """Columns of the groups: the columns of the members weighted by multiplicity."""
    n_group = int(group.max()) + 1 if group.size else 0
    weight = numpy.zeros((group.size, n_group))
    weight[numpy.arange(group.size), group] = multiplicity
    weight = weight / numpy.maximum(weight.sum(axis=0), 1e-30)[na, :]
    return column @ weight


def calc_peak_width_diagnosis(
        dict_pd: dict, dict_crystals: list, phase_name: str, flag_scan_layers: bool = True):
    """Sharp and broad reflection groups of a phase, the fault rule and the layers.

    Returns
    -------
    dict
        "groups": list of (2theta in degrees, hkl list, intensity, broad part,
        significance); "rules": find_fault_rule; "chi_sq_n": chi_sq/n of the
        model with the instrumental width only, with broad peaks and (when
        scanned) with the best stacking faults; "faults": the best
        (normal, shift, probability) or None.
    """
    from cryspy.procedure_rhochi.rhochi_pd import (
        calc_chi_sq_for_pd_by_dictionary, calc_index_hkl_multiplicity_for_phase)
    from cryspy.A_functions_base.peak_width_diagnostics import (
        BROAD_FACTORS, MIN_SEPARATION, MAX_SHARP_PART, MIN_SIGNIFICANCE,
        calc_lorentzian_columns, find_fault_rule, fit_sharp_and_broad, group_reflections)
    from cryspy.A_functions_base.powder_diffraction_const_wavelength import calc_lorentz_factor

    dict_in_out = {}
    calc_chi_sq_for_pd_by_dictionary(dict_pd, dict_crystals, dict_in_out=dict_in_out)
    ttheta_zs = dict_in_out["ttheta_corrected"]
    y_exp, sigma = take_signal_exp(dict_pd)
    y_net = y_exp - dict_in_out["signal_background"]
    # the other phases are kept as calculated
    for name in dict_pd["phase_name"]:
        if name.lower() != phase_name:
            d_phase = dict_in_out.get(f"dict_in_out_{name.lower():}", {})
            if "signal_plus" in d_phase:
                y_net = y_net - d_phase["signal_plus"] - d_phase["signal_minus"]
    flag_point = ~numpy.asarray(dict_pd["excluded_points"], dtype=bool) & (sigma > 0)

    dict_crystal = [d for d in dict_crystals if d["name"].lower() == phase_name][0]
    ind_phase = [name.lower() for name in dict_pd["phase_name"]].index(phase_name)
    ucp = dict_crystal["unit_cell_parameters"]
    wavelength = float(numpy.atleast_1d(dict_pd["wavelength"])[0])
    sthovl_min = max(numpy.sin(0.5 * ttheta_zs.min()) / wavelength, 1e-4)
    sthovl_max = numpy.sin(min(0.5 * ttheta_zs.max(), 0.5 * numpy.pi - 1e-6)) / wavelength
    index_hkl, multiplicity = calc_index_hkl_multiplicity_for_phase(
        dict_crystal, sthovl_min, sthovl_max, ucp)
    if "full_symm_elems" in dict_crystal:
        # free intensities of systematically absent reflections would take
        # any diffuse intensity
        from cryspy.A_functions_base.structure_factor import \
            calc_extinction_rule_by_symmetry_elements
        flag_absent = calc_extinction_rule_by_symmetry_elements(
            index_hkl, dict_crystal["full_symm_elems"], dict_crystal["translation_elems"],
            False)
        index_hkl, multiplicity = index_hkl[:, ~flag_absent], multiplicity[~flag_absent]
    from cryspy.A_functions_base.unit_cell import calc_sthovl_by_unit_cell_parameters
    sthovl_hkl = calc_sthovl_by_unit_cell_parameters(index_hkl, ucp)[0]
    ttheta_hkl = 2. * numpy.arcsin(numpy.clip(sthovl_hkl * wavelength, 0., 1.))
    u, v, w, x, y = dict_pd["resolution_parameters"] + dict_pd["phase_resolution_parameters"][:, ind_phase]
    p_ig = dict_pd["phase_ig"][ind_phase]
    p_1, p_2, p_3, p_4 = dict_pd.get("asymmetry_parameters", numpy.zeros(4))
    k, cthm = dict_pd["k"], dict_pd["cthm"]

    x_deg = numpy.degrees(ttheta_zs)
    lf_point = calc_lorentz_factor(ttheta_zs, k=k, cthm=cthm, flag_ttheta=None)[0]
    lf_hkl = calc_lorentz_factor(ttheta_hkl, k=k, cthm=cthm, flag_ttheta=None)[0]
    factor = lf_point[:, na] / lf_hkl[na, :]
    fwhm = calc_h_pv(calc_h_g(u, v, w, p_ig, ttheta_hkl)[0], calc_h_l(x, y, ttheta_hkl)[0])[0]
    column = calc_profile_pseudo_voight(
        ttheta_zs, ttheta_hkl, u, v, w, p_ig, x, y, p_1, p_2, p_3, p_4)[0] * factor
    group = group_reflections(numpy.degrees(ttheta_hkl), fwhm)
    column_sharp = merge_columns(column, group, multiplicity)
    column_broad = [merge_columns(calc_lorentzian_columns(
        x_deg, numpy.degrees(ttheta_hkl), factor_b * fwhm) * factor, group, multiplicity)
        for factor_b in BROAD_FACTORS]
    fp = flag_point
    i_sharp, i_broad, s_int, chi_sq_broad = fit_sharp_and_broad(
        y_net[fp], sigma[fp], column_sharp[fp], [c[fp] for c in column_broad], x_deg[fp])
    _, _, _, chi_sq_sharp = fit_sharp_and_broad(
        y_net[fp], sigma[fp], column_sharp[fp], [], x_deg[fp])
    n_point = max(int(numpy.count_nonzero(fp)), 1)
    total = i_sharp + i_broad
    significance = total / s_int
    broad_part = numpy.where(total > 0., i_broad / numpy.where(total > 0., total, 1.), 0.)

    l_group, l_index_group = [], []
    for i_group in range(int(group.max()) + 1 if group.size else 0):
        members = numpy.flatnonzero(group == i_group)
        l_index_group.append(index_hkl[:, members])
        l_group.append((float(numpy.degrees(ttheta_hkl[members]).mean()),
                        [tuple(int(v) for v in index_hkl[:, m]) for m in members],
                        float(total[i_group]), float(broad_part[i_group]),
                        float(significance[i_group])))
    flag_sig = significance > MIN_SIGNIFICANCE
    l_rule = find_fault_rule([l_index_group[i] for i in numpy.flatnonzero(flag_sig)],
                             broad_part[flag_sig], numpy.sqrt(significance[flag_sig])) \
        if numpy.count_nonzero(flag_sig) >= 3 else []
    result = {"groups": l_group, "rules": l_rule,
              "chi_sq_n": {"instrumental width": chi_sq_sharp / n_point,
                           "free broad peaks": chi_sq_broad / n_point},
              "faults": None, "all_broad": None}
    if numpy.any(flag_sig):
        result["all_broad"] = float(numpy.average(
            broad_part[flag_sig], weights=numpy.sqrt(significance[flag_sig])))
    if not (l_rule and l_rule[0]["separation"] > MIN_SEPARATION
            and l_rule[0]["broad_other"] < MAX_SHARP_PART and flag_scan_layers):
        return result

    # the layers: the normal and the probability with the best chi_sq
    shift = l_rule[0]["shift"]
    best = None
    for normal in form_candidate_normals(shift):
        for alpha in ALPHA_GRID:
            row = {"normal_index": normal, "shift": shift, "probability": alpha}
            flag_faulted, profile_faulted = calc_stacking_fault_profile(
                row, dict_crystal, index_hkl, ucp, wavelength, ttheta_zs,
                u, v, w, p_ig, x, y, p_1, p_2, p_3, p_4)
            column_f = numpy.array(column, copy=True)
            column_f[:, flag_faulted] = profile_faulted * factor[:, flag_faulted]
            chi_sq = fit_sharp_and_broad(
                y_net[fp], sigma[fp], merge_columns(column_f, group, multiplicity)[fp], [],
                x_deg[fp])[3]
            if best is None or chi_sq < best[0]:
                best = (chi_sq, normal, alpha)
    if best is not None:
        result["faults"] = {"normal": best[1], "shift": shift, "probability": best[2]}
        result["chi_sq_n"]["stacking faults"] = best[0] / n_point
    return result


def form_diagnosis_lines(result: dict, phase_label: str, flag_table: bool = True) -> list:
    """Lines of the diagnosis (keywords first, as diagnose_fit)."""
    from cryspy.A_functions_base.peak_width_diagnostics import (
        MIN_SEPARATION, MAX_SHARP_PART, MIN_SIGNIFICANCE)
    ls_out = []
    l_rule = result["rules"]
    flag_rule = bool(l_rule) and l_rule[0]["separation"] > MIN_SEPARATION and \
        l_rule[0]["broad_other"] < MAX_SHARP_PART
    if flag_rule:
        rule = l_rule[0]
        s_t = " ".join(f"{v:g}" for v in rule["shift"])
        ls_out.append(
            f"SELECTIVE_WIDTH ({phase_label:}): reflections with {rule['rule']:} are broad "
            f"(mean broad part {rule['broad_faulted']:.2f} of {rule['n_faulted']} groups), the "
            f"others are sharp ({rule['broad_other']:.2f} of {rule['n_other']} groups): random "
            f"stacking faults with the shift t = ({s_t:}) (H.t not an integer), not the "
            "profile parameters or the background.")
    elif result["all_broad"] is not None and result["all_broad"] > 0.5:
        ls_out.append(f"ALL_BROAD ({phase_label:}): all reflections are broader than the "
                      f"instrument (mean broad part {result['all_broad']:.2f}): size or strain "
                      "broadening of the phase (resolution parameters of the phase, IG, Y).")
    else:
        ls_out.append(f"WIDTHS ({phase_label:}): no group of reflections is clearly broader "
                      "than the others.")
    chi = result["chi_sq_n"]
    ls_out.append("CHI_SQ_N (" + phase_label + ", free intensities): " + "; ".join(
        f"{key:} {value:.4g}" for key, value in chi.items()) + ".")
    faults = result["faults"]
    if faults is not None:
        s_n = " ".join(str(int(v)) for v in faults["normal"])
        s_t = " ".join(f"{v:g}" for v in faults["shift"])
        ls_out.append(
            f"STACKING_FAULTS ({phase_label:}): best layers ({s_n:}), shift ({s_t:}), "
            f"probability about {faults['probability']:.2f} per layer of the normal ({s_n:}). "
            f"Use add_stacking_fault(normal=\"{s_n:}\", shift=\"{s_t:}\", probability="
            f"{faults['probability']:.2f}, phase_label=\"{phase_label:}\") and refine; if the "
            "cell has n layers along the normal, give n times the normal (e.g. \"0 0 2\"): "
            "the probability is then per layer.")
    if flag_table:
        ls_out.append("")
        ls_out.append("| 2theta | hkl | intensity | broad part | I/sigma |")
        ls_out.append("|---|---|---|---|---|")
        l_group = [item for item in result["groups"] if item[4] > MIN_SIGNIFICANCE]
        for position, l_hkl, intensity, broad, significance in l_group[:MAX_TABLE_ROWS]:
            s_hkl = ", ".join("".join(str(v) for v in hkl) for hkl in l_hkl[:4])
            if len(l_hkl) > 4:
                s_hkl += ", ..."
            ls_out.append(f"| {position:.2f} | {s_hkl:} | {intensity:.4g} | {broad:.2f} | "
                          f"{significance:.0f} |")
        if len(l_group) > MAX_TABLE_ROWS:
            ls_out.append(f"| ... | {len(l_group) - MAX_TABLE_ROWS} more groups | | | |")
    return ls_out


def diagnose_peak_widths(rcif_object: GlobalN, experiment: str = "", phase_label: str = ""):
    """Diagnose which powder reflections are broader than the instrument and why (stacking faults, size).

    Every group of reflections of the phase is fitted with a sharp peak
    (the profile of the phase) and two broad ones, with free intensities
    and a linear correction of the background (the structure is not
    used). The part of every group in the broad peaks shows which
    reflections are broad. If they follow a rule of random stacking faults
    (k odd, h + k odd, h - k not 3n, ...), the layers and the probability
    that describe the lines best are found and the call of
    add_stacking_fault is suggested; if all reflections are broad, the
    size or strain of the phase is suggested. A table of the groups
    (2theta, hkl, intensity, broad part, significance) follows. Only the
    1D powder experiments at constant wavelength (pd).

    Parameters
    ----------
    rcif_object : GlobalN
        The object with the crystals and the powder experiments.
    experiment : str
        Name of the powder experiment; empty for all pd experiments.
    phase_label : str
        Label of the phase; empty for every phase of the experiments.

    Returns
    -------
    str
        Diagnosis in markdown.
    """
    obj_dict = rcif_object.get_dictionary()
    dict_crystals = [value for key, value in obj_dict.items()
                     if isinstance(value, dict) and value.get("type_name", "").startswith("crystal")]
    ls_out = []
    for pd in take_powder_1d_experiments(rcif_object, experiment):
        dict_pd = obj_dict[f"pd_{pd.data_name.lower():}"]
        ls_out.append(f"# {pd.data_name:}")
        for name in dict_pd["phase_name"]:
            name = name.lower()
            if phase_label.strip() and name != phase_label.strip().lower():
                continue
            result = calc_peak_width_diagnosis(dict_pd, dict_crystals, name)
            ls_out.extend(form_diagnosis_lines(result, name))
            ls_out.append("")
    if len(ls_out) == 0:
        return "No powder experiments (pd)" + (f" named '{experiment:}'." if experiment else ".")
    return "\n".join(ls_out).strip()
