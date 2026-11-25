"""Stacking faults in 1D powder experiments (loop stacking_fault).

A phase built of layers normal to N = (h k l) can have stacking faults:
the next layer is placed, with the probability alpha after a normal
contact (beta after a fault), with an extra in-plane shift t. The ideal
stacking operation is found in the space group, the layer is cut from the
structure at the widest gap (or at layer_cut). The reflections that see
the faults are spread along rods normal to the layers with the intensity of
the Hendricks-Teller model (structure factor of the layer along the rod,
cryspy/A_functions_base/stacking_faults.py), the others stay sharp Bragg
peaks. Rietveld: the intensity of the rods comes from the structure (the
integrated intensity of a rod is kept); profile matching: the shapes come
from the model, the intensities are free. The spread lines are convolved
with the instrumental profile (resolution and asymmetry p1-p4).

The procedures add the rows and report them; diagnose_peak_widths finds
the groups of reflections broader than the instrument and the fault
vector and the layers that explain them. They are shown in the menus of
CrysPy Editor and used by AI assistants.
"""
import itertools
import logging

import numpy

from cryspy.A_functions_base.matrix_operations import calc_m_v
from cryspy.A_functions_base.powder_diffraction_const_wavelength import (
    calc_h_g, calc_h_l, calc_h_pv, calc_lorentz_factor, calc_profile_pseudo_voight)
from cryspy.A_functions_base.stacking_faults import (
    calc_atoms_of_cell,
    calc_fault_variants,
    calc_flag_faulted,
    calc_layer_cut,
    calc_layer_normal,
    calc_layer_structure_factor,
    calc_orientations,
    calc_reciprocal_metric,
    calc_rod_intensity,
    calc_rod_width,
    calc_symmetry_operations,
    PROBABILITY_MIN,
    find_stacking_operation,
    select_layer,
    split_normal,
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
MAX_GRID_POINTS = 20000
# points along the rods: at least 2 per step of the grid in |Q|, refined
# near the Bragg positions (geometric steps down to REFINE_MIN)
MAX_ROD_POINTS = 1500000
REFINE_FACTOR = 0.6
REFINE_MIN = 1e-7
CHUNK = 100000


def take_stacking_fault_of_phase(dict_pd: dict, phase_name: str):
    """Row of the phase in the loop stacking_fault of the experiment.

    Returns
    -------
    dict or None
        None when the phase has no row, otherwise "index", "normal_index"
        (3,), "shift" (3,), "probability" alpha, "probability_repeat" beta
        (alpha when not given) and "layer_cut" (None: automatic).
    """
    if "stacking_fault_phase_name" not in dict_pd:
        return None
    l_index = numpy.flatnonzero(dict_pd["stacking_fault_phase_name"] == phase_name)
    if l_index.size == 0:
        return None
    index = int(l_index[0])
    alpha = float(numpy.clip(dict_pd["stacking_fault_probability"][index], 0., 1.))
    beta = dict_pd.get("stacking_fault_probability_repeat", numpy.full(index + 1, numpy.nan))[index]
    layer_cut = dict_pd.get("stacking_fault_layer_cut", numpy.full(index + 1, numpy.nan))[index]
    return {
        "index": index,
        "normal_index": numpy.array(dict_pd["stacking_fault_normal_index"][:, index], dtype=int),
        "shift": numpy.array(dict_pd["stacking_fault_shift"][:, index], dtype=float),
        "probability": alpha,
        "probability_repeat": alpha if not numpy.isfinite(beta) else float(numpy.clip(beta, 0., 1.)),
        "layer_cut": None if not numpy.isfinite(layer_cut) else float(layer_cut),
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


def take_virtual_grid(ttheta_zs, u, v, w, p_ig, x, y, dict_cache: dict = None):
    """The grid of form_virtual_grid, kept in dict_cache while it still fits.

    The grid is made again only when the data are outside of it or the
    narrowest instrumental line is shorter than 3 steps: small changes of
    the zero shift and of the resolution during a refinement do not move
    it (the chi_sq stays a smooth function of the parameters).
    """
    grid = None if dict_cache is None else dict_cache.get("stacking_fault_grid")
    if grid is not None and grid.size > 1:
        h_pv = calc_h_pv(calc_h_g(u, v, w, p_ig, ttheta_zs)[0], calc_h_l(x, y, ttheta_zs)[0])[0]
        h_pv = numpy.abs(h_pv[numpy.isfinite(h_pv)])
        step = grid[1] - grid[0]
        if (h_pv.size and grid[0] <= ttheta_zs.min() and grid[-1] >= ttheta_zs.max()
                and numpy.radians(float(h_pv.min())) >= 3. * step):
            return grid
    grid = form_virtual_grid(ttheta_zs, u, v, w, p_ig, x, y)
    if dict_cache is not None:
        dict_cache["stacking_fault_grid"] = grid
    return grid


def take_symmetry_of_crystal(dict_crystal: dict):
    """Reduced symmetry elements and centrosymmetry for the equivalent reflections."""
    if "reduced_symm_elems" in dict_crystal:
        return dict_crystal["reduced_symm_elems"], bool(dict_crystal["centrosymmetry"])
    if "full_mcif_elems" in dict_crystal:
        return dict_crystal["full_mcif_elems"][:13], False
    return numpy.array([[0], [0], [0], [1], [1], [0], [0], [0], [1], [0], [0], [0], [1]]), False


def take_operations_of_crystal(dict_crystal: dict):
    """All symmetry operations (rotations (n, 3, 3), translations (n, 3)) of the crystal."""
    p1 = numpy.array([[0], [0], [0], [1]], dtype=int)
    if "full_symm_elems" in dict_crystal:
        return calc_symmetry_operations(dict_crystal["full_symm_elems"],
                                        dict_crystal.get("translation_elems", p1))
    if "full_mcif_elems" in dict_crystal:
        return calc_symmetry_operations(dict_crystal["full_mcif_elems"][:13], p1)
    reduced_symm_elems, centrosymmetry = take_symmetry_of_crystal(dict_crystal)
    from cryspy.A_functions_base.symmetry_elements import calc_full_symm_elems_by_reduced
    full = calc_full_symm_elems_by_reduced(
        reduced_symm_elems, centrosymmetry, dict_crystal.get("centrosymmetry_position"),
        dict_crystal.get("translation_elems", p1))
    return calc_symmetry_operations(full, p1)


def form_atoms_of_crystal(dict_crystal: dict, unit_cell_parameters, flag_xrays: bool = False):
    """Atoms of the cell with their scattering amplitudes and displacements.

    Returns
    -------
    dict
        "xyz" (n, 3), "amplitude" (n,) (neutrons) or "table_sthovl" and
        "table_amplitude" (n, n_t) (X-rays), "b_iso" (n,), "beta" (n, 3, 3),
        "metric" (3, 3).
    """
    rotations, translations = take_operations_of_crystal(dict_crystal)
    atom_fract_xyz = numpy.asarray(dict_crystal["atom_fract_xyz"], dtype=float)
    if "atom_site_sc_fract" in dict_crystal:
        atom_fract_xyz = calc_m_v(dict_crystal["atom_site_sc_fract"], numpy.mod(atom_fract_xyz, 1.),
                                  flag_m=False, flag_v=False)[0] + dict_crystal["atom_site_sc_b"]
    xyz, index_atom, index_operation = calc_atoms_of_cell(atom_fract_xyz, rotations, translations)
    occupancy = numpy.asarray(dict_crystal["atom_occupancy"], dtype=float)[index_atom]
    beta = numpy.asarray(dict_crystal["atom_beta"], dtype=float)
    if "atom_site_aniso_sc_beta" in dict_crystal:
        atom_sc_beta = numpy.zeros((6,) + beta.shape, dtype=float)
        atom_sc_beta[:, :, dict_crystal["atom_site_aniso_index"]] = dict_crystal["atom_site_aniso_sc_beta"]
        beta = (atom_sc_beta * beta[na, :, :]).sum(axis=1)
    b_11, b_22, b_33, b_12, b_13, b_23 = beta
    beta_matrix = numpy.array([[b_11, b_12, b_13], [b_12, b_22, b_23], [b_13, b_23, b_33]])
    beta_matrix = beta_matrix.transpose(2, 0, 1)[index_atom]
    rotation_copy = rotations[index_operation].astype(float)
    beta_copy = numpy.einsum("nij,njk,nlk->nil", rotation_copy, beta_matrix, rotation_copy)
    atoms = {
        "xyz": xyz,
        "b_iso": numpy.asarray(dict_crystal["atom_b_iso"], dtype=float)[index_atom],
        "beta": beta_copy,
        "metric": calc_reciprocal_metric(unit_cell_parameters)[0],
    }
    if flag_xrays and "table_atom_scattering_amplitude" in dict_crystal:
        atoms["table_sthovl"] = numpy.asarray(dict_crystal["table_sthovl"], dtype=float)
        atoms["table_amplitude"] = numpy.real(numpy.asarray(
            dict_crystal["table_atom_scattering_amplitude"]))[index_atom] * occupancy[:, na]
    else:
        atoms["amplitude"] = numpy.asarray(
            dict_crystal["atom_scat_length_neutron"])[index_atom] * occupancy
    return atoms


def take_part(atoms: dict, xyz, index_copy) -> dict:
    """The atoms index_copy at the positions xyz (a layer)."""
    part = {"xyz": xyz, "b_iso": atoms["b_iso"][index_copy], "beta": atoms["beta"][index_copy],
            "metric": atoms["metric"]}
    if "table_amplitude" in atoms:
        part["table_sthovl"] = atoms["table_sthovl"]
        part["table_amplitude"] = atoms["table_amplitude"][index_copy]
    else:
        part["amplitude"] = atoms["amplitude"][index_copy]
    return part


def form_stacking_model(row: dict, dict_crystal: dict, unit_cell_parameters,
                        flag_xrays: bool = False) -> dict:
    """Layer, stacking operation, fault variants and orientations of a phase.

    Raises ValueError when the shift is not in the plane of the layers or
    when the layers are not related by a symmetry operation.
    """
    normal_index = numpy.asarray(row["normal_index"], dtype=int)
    shift = numpy.asarray(row["shift"], dtype=float)
    primitive, g = split_normal(normal_index)
    if abs(float(normal_index @ shift)) > 1e-6:
        raise ValueError(
            f"The fault shift {tuple(shift)} is not in the plane of the layers "
            f"{tuple(normal_index)} (N.t = {float(normal_index @ shift):g}, must be 0).")
    rotations, translations = take_operations_of_crystal(dict_crystal)
    operation = find_stacking_operation(normal_index, rotations, translations)
    if operation is None:
        raise ValueError(
            f"The layers normal to {tuple(normal_index)} are not related by a symmetry "
            f"operation of the crystal (no R with R^T N = N and N.t = 1): give the normal "
            f"{tuple(primitive)} or another multiple of it.")
    atoms = form_atoms_of_crystal(dict_crystal, unit_cell_parameters, flag_xrays)
    layer_cut = row.get("layer_cut")
    if layer_cut is None:
        layer_cut = calc_layer_cut(atoms["xyz"], normal_index)
    xyz_layer, index_copy = select_layer(atoms["xyz"], normal_index, layer_cut)
    rotation, translation = operation
    variants = calc_fault_variants(shift, normal_index, rotations, translations)
    if not numpy.any(numpy.abs(variants - numpy.round(variants)) > 1e-6):
        variants = numpy.zeros((0, 3))
    return {
        "normal_index": normal_index, "primitive": primitive, "g": g,
        "rotation": rotation, "translation": translation, "variants": variants,
        "orientations": calc_orientations(rotation), "layer_cut": layer_cut,
        "atoms": atoms, "layer": take_part(atoms, xyz_layer, index_copy),
        "n_atom_layer": int(index_copy.size),
    }


def calc_rod_samples(base, primitive, metric, s_max: float, d_tau: float):
    """Points tau of a rod Q = base + tau P with |Q| < s_max, refined near the integers."""
    pp = float(primitive @ metric @ primitive)
    bp = float(base @ metric @ primitive)
    bb = float(base @ metric @ base)
    discriminant = bp * bp - pp * (bb - s_max * s_max)
    if discriminant <= 0.:
        return numpy.zeros(0)
    root = numpy.sqrt(discriminant)
    tau_min, tau_max = (-bp - root) / pp, (-bp + root) / pp
    n_step = max(int(numpy.ceil((tau_max - tau_min) / d_tau)), 2)
    tau = numpy.linspace(tau_min, tau_max, n_step + 1)
    l_extra = []
    n_refine = int(numpy.ceil(numpy.log(REFINE_MIN / d_tau) / numpy.log(REFINE_FACTOR)))
    offset = d_tau * numpy.power(REFINE_FACTOR, numpy.arange(1, max(n_refine, 1) + 1))
    for integer in range(int(numpy.ceil(tau_min)), int(numpy.floor(tau_max)) + 1):
        l_extra.append(integer + offset)
        l_extra.append(integer - offset)
        l_extra.append(numpy.array([float(integer)]))
    if l_extra:
        tau = numpy.concatenate([tau] + l_extra)
    tau = numpy.unique(tau[(tau >= tau_min) & (tau <= tau_max)])
    return tau


def form_stacking_fault_setup(model: dict, index_hkl, unit_cell_parameters, wavelength: float,
                              reduced_symm_elems, centrosymmetry: bool, s_max: float,
                              d_s: float) -> dict:
    """Points of the rods and the terms of the model that do not depend on alpha and beta.

    Parameters
    ----------
    model : dict
        form_stacking_model.
    index_hkl : array (3, n_hkl)
        Reflections of the phase (one of every set of equivalents).
    s_max : float
        Largest |Q| = 2 sin(theta) / lambda (1/angstrom) of the rods.
    d_s : float
        Wanted step in |Q| along the rods (1/angstrom).
    """
    index_hkl = numpy.asarray(index_hkl, dtype=int)
    n_hkl = index_hkl.shape[1]
    orientations, g = model["orientations"], model["g"]
    primitive = model["primitive"].astype(float)
    metric = model["atoms"]["metric"]
    flag_faulted = calc_flag_faulted(index_hkl, orientations, model["variants"]) \
        if model["variants"].shape[0] else numpy.zeros(n_hkl, dtype=bool)
    # Bragg |F|^2 of every reflection (perfect crystal, the same model)
    f2_bragg = numpy.square(numpy.abs(calc_layer_structure_factor(
        index_hkl.transpose().astype(float), model["atoms"])))
    setup = {"flag_faulted": flag_faulted, "f2_bragg": f2_bragg, "g": g}
    index_f = numpy.flatnonzero(flag_faulted)
    if index_f.size == 0:
        return setup
    index_eq = calc_equivalent_reflections(index_hkl[:, index_f], reduced_symm_elems, centrosymmetry)
    index_eq = numpy.concatenate([index_eq, -index_eq], axis=2)
    n_eq = index_eq.shape[2]
    # points of the rods: column of the reflection and its weight
    d_point = {}
    for i_column in range(index_f.size):
        for i_eq in range(n_eq):
            key = tuple(int(v) for v in index_eq[:, i_column, i_eq])
            d_point.setdefault(key, [i_column, 0.])
            d_point[key][1] += 1. / n_eq
    points = numpy.array(list(d_point.keys()), dtype=float)  # (n_p, 3)
    point_column = numpy.array([value[0] for value in d_point.values()], dtype=int)
    point_weight = numpy.array([value[1] for value in d_point.values()], dtype=float)
    pp = float(primitive @ metric @ primitive)
    tau_point = (points @ metric @ primitive) / pp
    m_point = numpy.floor(tau_point + 0.5 - 1e-9)
    base = points - m_point[:, na] * primitive[na, :]
    key_base = numpy.round(base * 1000).astype(numpy.int64)
    _, rod_of_point = numpy.unique(key_base, axis=0, return_inverse=True)
    rod_of_point = numpy.ravel(rod_of_point)
    d_tau = d_s / numpy.sqrt(pp) / 2.
    l_q, l_cell = [], []
    n_total = 0
    for i_rod in range(int(rod_of_point.max()) + 1):
        members = numpy.flatnonzero(rod_of_point == i_rod)
        base_rod = base[members[0]]
        tau = calc_rod_samples(base_rod, primitive, metric, s_max, d_tau)
        if tau.size < 2:
            continue
        order = numpy.argsort(m_point[members])
        m_present = m_point[members][order]
        member_sorted = members[order]
        tau_mid = 0.5 * (tau[1:] + tau[:-1])
        boundaries = 0.5 * (m_present[1:] + m_present[:-1])
        nearest = member_sorted[numpy.searchsorted(boundaries, tau_mid)]
        l_q.append(base_rod[na, :] + tau[:, na] * primitive[na, :])
        l_cell.append((n_total + numpy.arange(tau.size - 1), numpy.diff(tau), tau_mid,
                       numpy.full(tau.size - 1, i_rod), nearest, base_rod))
        n_total += tau.size
    if n_total > MAX_ROD_POINTS:
        logger.warning("Stacking faults: %d points along the rods (more than %d), the "
                       "calculation is slow.", n_total, MAX_ROD_POINTS)
    q = numpy.concatenate(l_q, axis=0)
    i_first = numpy.concatenate([item[0] for item in l_cell])
    d_tau_cell = numpy.concatenate([item[1] for item in l_cell])
    nearest = numpy.concatenate([item[4] for item in l_cell])
    q_mid = numpy.concatenate([item[5][na, :] + item[2][:, na] * primitive[na, :] for item in l_cell])
    s_mid = numpy.sqrt(numpy.maximum(numpy.einsum("ni,ij,nj->n", q_mid, metric, q_mid), 0.))
    ttheta_mid = 2. * numpy.arcsin(numpy.clip(0.5 * s_mid * wavelength, 0., 1.))
    # terms of the model at the points
    shifts = [model["translation"]] + [model["translation"] + v for v in model["variants"]]
    n_g = orientations.shape[0]
    f_layer = numpy.zeros((n_g, q.shape[0]), dtype=complex)
    phase = numpy.zeros((n_g, len(shifts), q.shape[0]), dtype=complex)
    for k, rotation in enumerate(orientations):
        q_rot = q @ rotation  # R_k^T Q
        for start in range(0, q.shape[0], CHUNK):
            f_layer[k, start:start + CHUNK] = calc_layer_structure_factor(
                q_rot[start:start + CHUNK], model["layer"])
        for j, shift in enumerate(shifts):
            phase[k, j] = numpy.exp(2j * numpy.pi * (q_rot @ shift))
    sthovl_hkl = 0.5 * numpy.sqrt(numpy.einsum(
        "in,ij,jn->n", index_hkl[:, index_f].astype(float), metric, index_hkl[:, index_f].astype(float)))
    ttheta_hkl = 2. * numpy.arcsin(numpy.clip(sthovl_hkl * wavelength, 0., 1.))
    setup.update({
        "index_faulted": index_f, "f_layer": f_layer, "phase": phase,
        "cell_first": i_first, "cell_d_tau": d_tau_cell, "cell_ttheta": ttheta_mid,
        "cell_column": point_column[nearest], "cell_weight": point_weight[nearest],
        "ttheta_hkl_faulted": ttheta_hkl, "n_point": q.shape[0],
    })
    return setup


def calc_stacking_fault_intensity(setup: dict, alpha: float, beta: float):
    """Intensity of the cells of the rods (units of |F_cell|^2, per equivalent of the reflection)."""
    n_point = setup["n_point"]
    intensity = numpy.zeros(n_point)
    for start in range(0, n_point, CHUNK):
        stop = min(start + CHUNK, n_point)
        intensity[start:stop] = calc_rod_intensity(
            setup["f_layer"][:, start:stop], setup["phase"][:, :, start:stop], alpha, beta,
            setup["g"])
    first = setup["cell_first"]
    return 0.5 * (intensity[first] + intensity[first + 1]) * setup["cell_d_tau"] * \
        setup["cell_weight"]


def form_setup_key(row: dict, unit_cell_parameters, dict_crystal: dict, index_hkl,
                   wavelength: float, s_max: float, flag_xrays: bool):
    """Parameters the setup depends on (all but alpha and beta)."""
    l_part = [numpy.asarray(row["normal_index"], dtype=float), numpy.asarray(row["shift"], dtype=float),
              numpy.array([numpy.nan if row.get("layer_cut") is None else row["layer_cut"],
                           wavelength, round(s_max, 4), float(flag_xrays)]),
              numpy.asarray(unit_cell_parameters, dtype=float),
              numpy.asarray(index_hkl, dtype=float).ravel()]
    for key in ("atom_fract_xyz", "atom_occupancy", "atom_b_iso", "atom_beta"):
        if key in dict_crystal:
            l_part.append(numpy.real(numpy.asarray(dict_crystal[key], dtype=complex)).ravel())
    return numpy.round(numpy.concatenate(l_part), 9).tobytes()


def calc_stacking_fault_profile(
        row: dict, dict_crystal: dict, index_hkl, unit_cell_parameters, wavelength,
        ttheta_zs, u, v, w, p_ig, x, y, p_1, p_2, p_3, p_4, dict_cache: dict = None,
        flag_xrays: bool = False, k: float = 0., cthm: float = 0.91):
    """Profiles of the spread reflections of a phase at the points (Hendricks-Teller model).

    The intensity along the rods is put on a fine grid of 2theta and
    convolved with the instrumental profile; a reflection without faults
    keeps the profile of the phase.

    Parameters
    ----------
    dict_cache : dict, optional
        The setup of the rods is kept there ("stacking_fault_setup") and
        used again while only alpha, beta or the profile change.

    Returns
    -------
    flag_faulted : array (n_hkl,) of bool
    profile : array (n_point, n_faulted)
        Normalized lines (the units of calc_profile_pseudo_voight,
        1/degree), the Lorentz factor along the line relative to the one of
        the reflection included: they replace the columns of the flagged
        reflections of the profile of the phase.
    intensity : array (n_faulted,)
        Intensity of the spread reflections (per equivalent, units of
        |F|^2 of the model).
    f2_bragg : array (n_hkl,)
        |F|^2 of the reflections of the perfect crystal (the same units).
    """
    wavelength = float(numpy.atleast_1d(wavelength)[0])
    ttheta_grid = take_virtual_grid(ttheta_zs, u, v, w, p_ig, x, y, dict_cache)
    s_max = 2. * numpy.sin(0.5 * min(float(ttheta_grid[-1]), numpy.pi)) / wavelength
    step = float(ttheta_grid[1] - ttheta_grid[0]) if ttheta_grid.size > 1 else 1e-3
    d_s = step * max(numpy.cos(0.5 * float(ttheta_grid[-1])), 0.2) / wavelength
    key = form_setup_key(row, unit_cell_parameters, dict_crystal, index_hkl, wavelength, s_max,
                         flag_xrays)
    setup = None if dict_cache is None else dict_cache.get("stacking_fault_setup")
    if setup is None or setup.get("key") != key:
        reduced_symm_elems, centrosymmetry = take_symmetry_of_crystal(dict_crystal)
        model = form_stacking_model(row, dict_crystal, unit_cell_parameters, flag_xrays)
        setup = form_stacking_fault_setup(
            model, index_hkl, unit_cell_parameters, wavelength, reduced_symm_elems,
            centrosymmetry, s_max, d_s)
        setup["key"] = key
        if dict_cache is not None:
            dict_cache["stacking_fault_setup"] = setup
    flag_faulted = setup["flag_faulted"]
    if row["probability"] <= PROBABILITY_MIN:
        # no faults: the Bragg reflections
        flag_faulted = numpy.zeros_like(flag_faulted)
    if not numpy.any(flag_faulted):
        return flag_faulted, numpy.zeros((ttheta_zs.size, 0)), numpy.zeros(0), setup["f2_bragg"]
    cell_mass = calc_stacking_fault_intensity(setup, row["probability"], row["probability_repeat"])
    column = setup["cell_column"]
    n_faulted = int(numpy.count_nonzero(flag_faulted))
    intensity = numpy.bincount(column, weights=cell_mass, minlength=n_faulted)
    ttheta_cell = setup["cell_ttheta"]
    lf_ratio = calc_lorentz_factor(ttheta_cell, k=k, cthm=cthm)[0] / \
        calc_lorentz_factor(setup["ttheta_hkl_faulted"], k=k, cthm=cthm)[0][column]
    position = (ttheta_cell - ttheta_grid[0]) / step
    flag = (position >= 0.) & (position < ttheta_grid.size - 1)
    index_low = numpy.floor(position[flag]).astype(int)
    fraction = position[flag] - index_low
    weight = (cell_mass * lf_ratio)[flag] / numpy.where(intensity > 0., intensity, 1.)[column[flag]]
    mass = numpy.zeros((ttheta_grid.size, n_faulted))
    numpy.add.at(mass, (index_low, column[flag]), weight * (1. - fraction))
    numpy.add.at(mass, (index_low + 1, column[flag]), weight * fraction)
    profile_grid = calc_profile_pseudo_voight(
        ttheta_zs, ttheta_grid, u, v, w, p_ig, x, y, p_1, p_2, p_3, p_4)[0]
    return flag_faulted, profile_grid @ mass, intensity, setup["f2_bragg"]


def calc_stacking_fault_scale(f2_structure, f2_model, multiplicity) -> float:
    """Factor from the |F|^2 of the model of the rods to the |F|^2 of the structure factors of the phase."""
    weight = numpy.asarray(multiplicity, dtype=float)
    model = float((numpy.asarray(f2_model) * weight).sum())
    if model <= 0.:
        return 0.
    return float((numpy.asarray(f2_structure) * weight).sum()) / model


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
        probability: float = 0.1, phase_label: str = "", experiment: str = "",
        probability_repeat: float = -1., layer_cut: float = -1.):
    """Add stacking faults to a powder phase: only the reflections that see the fault shift become broad or asymmetric.

    The layers normal to the reciprocal vector N = (h k l) (distance 1/|N|
    between them; "0 0 2" are two layers per c) are stacked by a symmetry
    operation of the crystal (found in the space group) and, at a fault,
    with an extra in-plane shift t. The faults follow a normal contact
    with the probability alpha and a fault with the probability beta
    (probability_repeat; equal to alpha when not given: independent
    faults). The reflections that see the shift are spread along rods
    normal to the layers with the intensity of the Hendricks-Teller model
    (the structure factor of the layer changes along the rod: asymmetric
    lines, a sharp front and a tail to high angles for the reflections in
    the plane of the layers); the others stay sharp. Rietveld: the
    intensities of the rods come from the structure; profile matching:
    the shapes from the model, the intensities free. The probabilities are
    refined. Only the 1D powder experiments at constant wavelength (pd).

    Parameters
    ----------
    rcif_object : GlobalN
        The object with the crystals and the powder experiments.
    normal : str
        Miller indices of the layer planes as "h k l".
    shift : str
        Fault vector t in the plane of the layers, fractions of the cell, as "x y z" (e.g. "0 1/2 0").
    probability : float
        Start value of the probability of a fault after a normal contact (0 ... 1).
    phase_label : str
        Label of the phase; empty for every phase of the experiments.
    experiment : str
        Name of the powder experiment; empty for all pd experiments.
    probability_repeat : float
        Start value of the probability of a fault after a fault (0 ... 1); negative: not used (independent faults).
    layer_cut : float
        N.r of the cut between two layers; negative: the widest gap between the atoms.

    Returns
    -------
    str
        The phases with the faults, the stacking operation, the layer and the reflections they spread.
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
    if abs(float(normal_index @ shift_vector)) > 1e-6:
        raise ValueError(f"The fault shift ({shift:}) is not in the plane of the layers "
                         f"({normal:}): N.t must be 0.")
    probability_repeat = float(probability_repeat)
    if probability_repeat > 1.:
        raise ValueError("probability_repeat must be between 0 and 1 (negative: not used).")
    layer_cut = float(layer_cut)
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
            crystal = take_crystal_by_label(rcif_object, label)
            if crystal is not None:
                # the layers must be related by a symmetry operation
                form_stacking_model({"normal_index": normal_index, "shift": shift_vector,
                                     "layer_cut": None}, crystal.get_dictionary(),
                                    crystal.get_dictionary()["unit_cell_parameters"])
            row = StackingFault(
                phase_label=label, normal_h=int(normal_index[0]), normal_k=int(normal_index[1]),
                normal_l=int(normal_index[2]), shift_x=float(shift_vector[0]),
                shift_y=float(shift_vector[1]), shift_z=float(shift_vector[2]),
                probability=probability)
            row.probability_refinement = True
            if probability_repeat >= 0.:
                row.probability_repeat = probability_repeat
                row.probability_repeat_refinement = True
            if layer_cut >= 0.:
                row.layer_cut = layer_cut
            loop.items.append(row)
            s_line = f"{pd.data_name:}, phase {label:}: "
            if crystal is not None:
                s_line += form_report_of_row(row, crystal, float(pd.setup.wavelength))[2:]
            else:
                s_line += (f"layers ({normal.strip():}), fault shift ({shift.strip():}), "
                           f"probability {probability:.3f}.")
            ls_out.append(s_line + " The probabilities are refined.")
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


def describe_operation(rotation, translation) -> str:
    """The operation (R, t) as x, y, z (e.g. "x, -y, z+1/2")."""
    from fractions import Fraction
    names = "xyz"
    l_part = []
    for i in range(3):
        text = ""
        for j in range(3):
            value = int(rotation[i, j])
            if value != 0:
                text += ("+" if value > 0 and text else "") + ("-" if value < 0 else "") + names[j]
        fraction = Fraction(float(translation[i])).limit_denominator(12)
        if fraction != 0:
            text += ("+" if fraction > 0 else "-") + str(abs(fraction))
        l_part.append(text)
    return ", ".join(l_part)


def form_report_of_row(row, crystal, wavelength: float) -> str:
    """Lines of the report for one phase."""
    from cryspy.procedure_rhochi.rhochi_pd import calc_index_hkl_multiplicity_for_phase
    from cryspy.A_functions_base.structure_factor import \
        calc_extinction_rule_by_symmetry_elements
    normal_index, shift_vector = row.get_normal_index(), row.get_shift()
    dict_crystal = crystal.get_dictionary()
    ucp = numpy.array(dict_crystal["unit_cell_parameters"], dtype=float)
    alpha = float(row.probability)
    beta = row.get_probability_repeat()
    d_row = {"normal_index": normal_index, "shift": shift_vector, "probability": alpha,
             "probability_repeat": alpha if not numpy.isfinite(beta) else beta,
             "layer_cut": None if not numpy.isfinite(row.get_layer_cut()) else row.get_layer_cut()}
    model = form_stacking_model(d_row, dict_crystal, ucp)
    s_sigma = f" ± {row.probability_sigma:.4f}" if row.probability_sigma else ""
    d_layer = calc_layer_normal(normal_index, ucp)[1]
    s_line = (f"- {row.phase_label:}: layers ({' '.join(str(i) for i in normal_index):}) "
              f"{d_layer:.3f} Å apart, cut at N.r = {model['layer_cut']:.3f} "
              f"({model['n_atom_layer']} atoms per layer); next layer by "
              f"({describe_operation(model['rotation'], model['translation'])}); fault shift "
              f"({' '.join(f'{v:g}' for v in shift_vector):})")
    if model["variants"].shape[0] > 1:
        s_line += f" and {model['variants'].shape[0] - 1} equivalent shift(s)"
    s_line += f"; probability after a normal contact {alpha:.4f}{s_sigma}"
    if numpy.isfinite(beta):
        s_sigma_b = f" ± {row.probability_repeat_sigma:.4f}" if row.probability_repeat_sigma else ""
        s_line += f", after a fault {beta:.4f}{s_sigma_b}"
        if beta > alpha + 0.05:
            s_line += " (runs of faults: domains of the other stacking)"
        elif beta < alpha - 0.05:
            s_line += " (faults avoid each other)"
    else:
        s_line += " (independent faults)"
    if model["variants"].shape[0] == 0:
        return s_line + ". The shift is a lattice vector: no reflection sees the faults."
    index_hkl = calc_index_hkl_multiplicity_for_phase(
        dict_crystal, 1e-4, 0.999 / wavelength, ucp)[0]
    if "full_symm_elems" in dict_crystal:
        flag_absent = calc_extinction_rule_by_symmetry_elements(
            index_hkl, dict_crystal["full_symm_elems"], dict_crystal["translation_elems"], False)
        index_hkl = index_hkl[:, ~flag_absent]
    flag = calc_flag_faulted(index_hkl, model["orientations"], model["variants"])
    index_f = index_hkl[:, flag]
    s_line += f". Spread reflections: {describe_rule(shift_vector):}"
    if index_f.shape[1] > 0:
        # width of the rod at the first spread reflection
        q = index_f[:, :1].transpose().astype(float)
        shifts = [model["translation"]] + [model["translation"] + v for v in model["variants"]]
        phase = numpy.array([[numpy.exp(2j * numpy.pi * ((q @ rotation) @ shift))
                              for shift in shifts] for rotation in model["orientations"]])
        width = float(calc_rod_width(phase, alpha, d_row["probability_repeat"], model["g"])[0])
        primitive_length = 1. / (d_layer * model["g"])
        hwhm = width * primitive_length
        if hwhm > 0.:
            s_line += (f"; HWHM of the rods about {hwhm:.4f} 1/Å, correlation length of the "
                       f"stacking about {1. / (2. * numpy.pi * hwhm):.1f} Å "
                       f"({1. / (2. * numpy.pi * hwhm * d_layer):.1f} layers)")
        s_line += ". First spread reflections: " + ", ".join(
            "".join(str(int(v)) for v in index_f[:, i]) for i in range(min(8, index_f.shape[1])))
    return s_line + "."


# ---------------------------------------------------------------------------
# Diagnosis of the peak widths

ALPHA_GRID = (0.03, 0.07, 0.12, 0.18, 0.25, 0.35)
MAX_TABLE_ROWS = 30


def take_signal_exp(dict_pd: dict):
    """Measured intensity (sum of up and down) and its sigma at the points."""
    if "signal_exp_plus" in dict_pd and "signal_exp_minus" in dict_pd:
        plus, minus = dict_pd["signal_exp_plus"], dict_pd["signal_exp_minus"]
        return plus[0] + minus[0], numpy.sqrt(numpy.square(plus[1]) + numpy.square(minus[1]))
    return dict_pd["signal_exp"][0], dict_pd["signal_exp"][1]


LAYERS_PER_CELL = (1, 2, 3, 4, 6)


def form_candidate_normals(shift, dict_crystal: dict = None) -> list:
    """Normals g P (P primitive with indices -1 ... 1, one of +-P) with the shift in the plane.

    With dict_crystal, only the multiples g P whose layers are related by
    a symmetry operation of the crystal (the layers of the model).
    """
    rotations = translations = None
    if dict_crystal is not None:
        rotations, translations = take_operations_of_crystal(dict_crystal)
    l_normal = []
    for normal in itertools.product((-1, 0, 1), repeat=3):
        if normal == (0, 0, 0) or next(v for v in normal if v != 0) < 0:
            continue
        if abs(float(numpy.array(normal) @ shift)) > 1e-6:
            continue
        for g in LAYERS_PER_CELL:
            candidate = g * numpy.array(normal, dtype=int)
            if rotations is None:
                if g == 1:
                    l_normal.append(candidate)
                continue
            if find_stacking_operation(candidate, rotations, translations) is not None:
                l_normal.append(candidate)
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
    flag_xrays = str(numpy.atleast_1d(dict_pd.get("radiation", ["neutrons"]))[0]).startswith("X-rays")
    for normal in form_candidate_normals(shift, dict_crystal):
        dict_cache = {}
        for alpha in ALPHA_GRID:
            row = {"normal_index": normal, "shift": shift, "probability": alpha,
                   "probability_repeat": alpha, "layer_cut": None}
            try:
                flag_faulted, profile_faulted = calc_stacking_fault_profile(
                    row, dict_crystal, index_hkl, ucp, wavelength, ttheta_zs,
                    u, v, w, p_ig, x, y, p_1, p_2, p_3, p_4, dict_cache=dict_cache,
                    flag_xrays=flag_xrays, k=k, cthm=cthm)[:2]
            except ValueError:
                break
            column_f = numpy.array(column, copy=True)
            column_f[:, flag_faulted] = profile_faulted
            # the spread and the sharp members of a group have their own
            # intensities (their lines differ)
            group_f = numpy.unique(2 * group + flag_faulted, return_inverse=True)[1].ravel()
            chi_sq = fit_sharp_and_broad(
                y_net[fp], sigma[fp], merge_columns(column_f, group_f, multiplicity)[fp], [],
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
            f"probability about {faults['probability']:.2f} per layer. "
            f"Use add_stacking_fault(normal=\"{s_n:}\", shift=\"{s_t:}\", probability="
            f"{faults['probability']:.2f}, phase_label=\"{phase_label:}\") and refine the "
            "probability; then try probability_repeat (correlated faults).")
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
