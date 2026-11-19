"""Layers without correlations between them in powder experiments (Warren line shape).

A phase whose order (magnetic, nuclear or both) exists only inside the
layers normal to N = (h k l) scatters along rods instead of Bragg
reflections: its powder pattern has the asymmetric saw-tooth lines of
B. E. Warren, Phys. Rev. 59 (1941) 693. The loop layer_correlation of a
pd experiment gives, for a phase, the layer planes, which order is
two-dimensional, the in-plane correlation length (refined) and the scale
of the rods (refined). The two-dimensional part is taken away from the
Bragg reflections of the phase and calculated as rods from the atoms and
moments of the crystal; the rods are convolved with the instrumental
profile of the experiment (resolution and asymmetry p1-p4 included).

The procedures suggest the layers from the crystal structure (the
direction with the widest empty gap and the longest bonds across it),
add the rows and report them. They are shown in the menus of CrysPy
Editor and used by AI assistants.
"""
import itertools
import logging

import numpy

from cryspy.A_functions_base.layer_correlations import (
    calc_layer_cuts,
    calc_layer_index,
    calc_layer_normal,
    calc_layer_powder_density,
    calc_layer_rod_s_z,
    calc_layers_of_atoms,
    calc_rod_intensity,
    calc_rods,
)
from cryspy.A_functions_base.magnetic_form_factor import calc_form_factor
from cryspy.A_functions_base.matrix_operations import calc_m_v
from cryspy.A_functions_base.powder_diffraction_const_wavelength import (
    calc_h_g, calc_h_l, calc_h_pv, calc_profile_pseudo_voight)
from cryspy.A_functions_base.structure_scene import calc_atoms_in_unit_cell
from cryspy.B_parent_classes.cl_4_global import GlobalN
from cryspy.C_item_loop_classes.cl_1_layer_correlation import (
    CONTRIBUTIONS, LayerCorrelation, LayerCorrelationL)
from cryspy.E_data_classes.cl_1_crystal import Crystal
from cryspy.E_data_classes.cl_2_pd import Pd
from cryspy.E_data_classes.cl_2_pd2d import Pd2d
from cryspy.E_data_classes.cl_2_tof import TOF

logger = logging.getLogger(__name__)
na = numpy.newaxis

# sampling of the rods before the convolution: points per FWHM of the
# instrument and sub-points per point for the sharp front of the line
POINTS_PER_FWHM = 5
SUB_POINTS = 4


def take_layer_correlation_of_phase(dict_pd: dict, phase_name: str):
    """Row of the phase in the loop layer_correlation of the experiment.

    Returns
    -------
    dict or None
        None when the phase has no row, otherwise "index", "normal_index"
        (3,), "contribution", "length" (angstrom), "scale".
    """
    if "layer_correlation_phase_name" not in dict_pd:
        return None
    l_index = numpy.flatnonzero(dict_pd["layer_correlation_phase_name"] == phase_name)
    if l_index.size == 0:
        return None
    index = int(l_index[0])
    contribution = str(dict_pd["layer_correlation_contribution"][index]).strip().lower()
    if contribution not in CONTRIBUTIONS:
        raise ValueError(f"layer_correlation of '{phase_name:}': contribution "
                         f"'{contribution:}' is not one of {', '.join(CONTRIBUTIONS):}.")
    return {
        "index": index,
        "normal_index": numpy.array(dict_pd["layer_correlation_normal_index"][:, index], dtype=int),
        "contribution": contribution,
        "length": float(dict_pd["layer_correlation_length"][index]),
        "scale": float(dict_pd["layer_correlation_scale"][index]),
    }


def form_atoms_of_unit_cell(dict_crystal: dict, unit_cell_parameters,
                            flag_magnetic: bool):
    """All atoms of the unit cell from the dictionary of a crystal.

    Returns
    -------
    dict
        "fract_xyz" (3, n), "atom_index" (n,) in atom_site, "moment_ccs"
        (3, n) Cartesian moments in Bohr magnetons (0 for the atoms
        without a moment; None when flag_magnetic is False).
    """
    keys = dict_crystal.keys()
    if "full_mcif_elems" in keys:
        symm_elems = dict_crystal["full_mcif_elems"]
    else:
        symm_elems = dict_crystal["full_symm_elems"]
    fract = calc_m_v(dict_crystal["atom_site_sc_fract"],
                     numpy.mod(dict_crystal["atom_fract_xyz"], 1.),
                     flag_m=False, flag_v=False)[0] + dict_crystal["atom_site_sc_b"]
    n_atom = fract.shape[1]
    m_m = calc_layer_normal((0, 0, 1), unit_cell_parameters)[2]
    moment_fract = None
    if flag_magnetic:
        if "atom_ordered_index" not in keys:
            if "atom_fourier_index" in keys:
                raise ValueError(
                    "layer_correlation: the moments of the crystal are given by Fourier "
                    "components (k != 0); the magnetic layers need the moments of a "
                    "commensurate magnetic cell (atom_site_moment in a magnetic space group).")
            raise ValueError(
                "layer_correlation with a magnetic contribution needs ordered moments "
                "of the crystal (atom_site_moment in a magnetic space group).")
        moment_dn = numpy.zeros((3, n_atom), dtype=float)
        moment_dn[:, dict_crystal["atom_ordered_index"]] = \
            dict_crystal["atom_ordered_moment_crystalaxis_dn"]
        moment_fract = moment_dn / unit_cell_parameters[:3, na]
    fract_cell, atom_index, moment_cell, _ = calc_atoms_in_unit_cell(
        symm_elems, fract, moment_fract)
    moment_ccs = None if moment_cell is None else m_m @ moment_cell
    return {"fract_xyz": fract_cell, "atom_index": atom_index, "moment_ccs": moment_ccs}


def calc_magnetic_form_factors(dict_crystal: dict, atom_index, s):
    """Magnetic form factors (n_s, n) of the atoms; 0 for atoms without moments."""
    res = numpy.zeros((s.size, atom_index.size), dtype=float)
    if "atom_ordered_index" not in dict_crystal or "mag_atom_ordered_index" not in dict_crystal:
        return res
    atom_ordered_index = dict_crystal["atom_ordered_index"]
    mag_index = dict_crystal["mag_atom_ordered_index"]
    form_factor = calc_form_factor(
        0.5 * s[:, na],
        dict_crystal["mag_atom_lande_factor"][mag_index][na, :],
        dict_crystal["mag_atom_kappa"][mag_index][na, :],
        dict_crystal["mag_atom_j0_parameters"][:, mag_index][:, na, :],
        dict_crystal["mag_atom_j2_parameters"][:, mag_index][:, na, :])[0]
    for i_ordered, i_atom in enumerate(atom_ordered_index):
        res[:, atom_index == i_atom] = form_factor[:, i_ordered][:, na]
    return res


def form_virtual_grid(ttheta_zs, u, v, w, p_ig, x, y):
    """Points (radians) at which the rods are calculated before the convolution.

    The step is a fifth of the narrowest instrumental FWHM, never larger
    than the step of the data; the grid covers the data with a margin of
    5 FWHM.
    """
    h_pv = calc_h_pv(calc_h_g(u, v, w, p_ig, ttheta_zs)[0],
                     calc_h_l(x, y, ttheta_zs)[0])[0]
    h_min = numpy.radians(max(float(numpy.min(h_pv)), 1e-3))
    step = h_min / POINTS_PER_FWHM
    if ttheta_zs.size > 1:
        step = min(step, float(numpy.median(numpy.abs(numpy.diff(ttheta_zs)))))
    margin = 5. * numpy.radians(float(numpy.max(h_pv)))
    t_min = max(float(ttheta_zs.min()) - margin, step)
    t_max = min(float(ttheta_zs.max()) + margin, numpy.pi - step)
    n_point = int(numpy.ceil((t_max - t_min) / step)) + 1
    if n_point > 20000:
        n_point = 20000
    return numpy.linspace(t_min, t_max, n_point)


def calc_layer_density_on_grid(
        ttheta_grid, wavelength, dict_crystal, unit_cell_parameters, row,
        flag_nuclear: bool, flag_magnetic: bool):
    """Powder-averaged density of the rods averaged over every grid step."""
    wavelength = float(numpy.atleast_1d(wavelength)[0])
    normal_index = row["normal_index"]
    atoms = form_atoms_of_unit_cell(dict_crystal, unit_cell_parameters, flag_magnetic)
    fract_layer, layer_index, n_layer = calc_layers_of_atoms(
        atoms["fract_xyz"], normal_index, unit_cell_parameters)
    atom_index = atoms["atom_index"]
    s = 2. * numpy.sin(0.5 * ttheta_grid) / wavelength
    index_rod, g_rod, g_ccs = calc_rods(normal_index, unit_cell_parameters, float(s.max()))
    if g_rod.size == 0:
        return numpy.zeros_like(ttheta_grid)

    occupancy = dict_crystal["atom_occupancy"][atom_index]
    scattering_length = None
    if flag_nuclear:
        scattering_length = dict_crystal["atom_scat_length_neutron"][atom_index] * occupancy
    moment_ccs, form_factor = None, None
    if flag_magnetic:
        moment_ccs = atoms["moment_ccs"] * occupancy[na, :]
        form_factor = calc_magnetic_form_factors(dict_crystal, atom_index, s)
    b_iso = dict_crystal["atom_b_iso"][atom_index]

    # the structure factors change slowly: on the grid; the Warren factor
    # with the sharp front: averaged over sub-points of every step
    s_z = calc_layer_rod_s_z(s, g_rod)
    args = (normal_index, unit_cell_parameters, fract_layer, layer_index, n_layer)
    kwargs = {"scattering_length": scattering_length, "moment_ccs": moment_ccs,
              "magnetic_form_factor": form_factor, "b_iso": b_iso}
    h_plus = calc_rod_intensity(s, g_ccs, index_rod, s_z, *args, **kwargs)
    h_minus = calc_rod_intensity(s, g_ccs, index_rod, -s_z, *args, **kwargs)
    step = ttheta_grid[1] - ttheta_grid[0] if ttheta_grid.size > 1 else 0.
    density = numpy.zeros_like(ttheta_grid)
    for shift in (numpy.arange(SUB_POINTS) + 0.5) / SUB_POINTS - 0.5:
        s_sub = 2. * numpy.sin(0.5 * (ttheta_grid + shift * step)) / wavelength
        density += calc_layer_powder_density(s_sub, g_rod, h_plus, h_minus, row["length"])
    return density / SUB_POINTS


def calc_signal_layer_correlation(
        dict_pd: dict, dict_crystal: dict, dict_in_out_phase: dict, phase_name: str,
        unit_cell_parameters, wavelength, ttheta_zs, u, v, w, p_ig, x, y,
        p_1, p_2, p_3, p_4, k, cthm, radiation: str):
    """Signal of the rods of a phase (half of it: the part of one neutron spin).

    signal(2theta) = 0.5 scale P(2theta) (pi/180) (2/lambda^3) D(s) convolved
    with the instrumental profile; D(s) is the powder-averaged density of
    the rods (layer_correlations.py) normalized as m |F|^2 delta(s - s_hkl)
    / s^2, so that a Bragg reflection with the same normalization gives
    0.5 scale L(theta) m |F|^2 (CrysPy). P is the polarization factor of
    the Lorentz factor (1 for neutrons).

    Returns
    -------
    signal : array (n_point,) or None
        None when the phase has no row in layer_correlation.
    row : dict or None
        The row (take_layer_correlation_of_phase).
    """
    row = take_layer_correlation_of_phase(dict_pd, phase_name)
    if row is None:
        return None, None
    if not radiation.startswith("neutrons"):
        raise ValueError("layer_correlation is implemented for neutrons only.")
    flag_nuclear = row["contribution"] in ("nuclear", "nuclear and magnetic")
    flag_magnetic = row["contribution"] in ("magnetic", "nuclear and magnetic")
    ttheta_grid = form_virtual_grid(ttheta_zs, u, v, w, p_ig, x, y)
    density = calc_layer_density_on_grid(
        ttheta_grid, wavelength, dict_crystal, unit_cell_parameters, row,
        flag_nuclear, flag_magnetic)
    wavelength = float(numpy.atleast_1d(wavelength)[0])
    polarization = 1. - k + k * cthm * numpy.square(numpy.cos(ttheta_grid))
    intensity = 0.5 * row["scale"] * polarization * numpy.pi / 180. * 2. / \
        wavelength ** 3 * density
    dict_in_out_phase["layer_correlation_ttheta"] = ttheta_grid
    dict_in_out_phase["layer_correlation_intensity"] = intensity

    # the instrumental profile does not depend on the rods: kept while the
    # resolution and the grid stay the same
    key = (tuple(numpy.round(ttheta_zs[[0, -1]], 12)), ttheta_zs.size,
           tuple(numpy.round(ttheta_grid[[0, -1]], 12)), ttheta_grid.size,
           float(u), float(v), float(w), float(p_ig), float(x), float(y),
           float(p_1), float(p_2), float(p_3), float(p_4))
    if dict_in_out_phase.get("layer_correlation_profile_key") == key:
        profile = dict_in_out_phase["layer_correlation_profile"]
    else:
        profile = calc_profile_pseudo_voight(
            ttheta_zs, ttheta_grid, u, v, w, p_ig, x, y, p_1, p_2, p_3, p_4)[0]
        dict_in_out_phase["layer_correlation_profile"] = profile
        dict_in_out_phase["layer_correlation_profile_key"] = key
    step_deg = numpy.degrees(ttheta_grid[1] - ttheta_grid[0]) if ttheta_grid.size > 1 else 1.
    signal = profile @ (intensity * step_deg)
    dict_in_out_phase["signal_layer_correlation"] = 2. * signal
    return signal, row


def calc_in_plane_mask(index_hkl, normal_index, unit_cell_parameters):
    """True for the reflections off the line 00l of the stacking (they become rods)."""
    normal_ccs, d_n, m_m = calc_layer_normal(normal_index, unit_cell_parameters)
    h_ccs = numpy.linalg.inv(m_m).transpose() @ index_hkl
    g_ccs = h_ccs - normal_ccs[:, na] * (normal_ccs @ h_ccs)[na, :]
    return numpy.linalg.norm(g_ccs, axis=0) > 1e-8


# ---------------------------------------------------------------------------
# Layers from the crystal structure

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


def form_candidate_normals(max_index: int = 2) -> list:
    """Primitive Miller indices (h k l) up to max_index, one of +-N."""
    l_normal = []
    rng = range(-max_index, max_index + 1)
    for normal in itertools.product(rng, rng, rng):
        if normal == (0, 0, 0) or numpy.gcd.reduce(numpy.abs(normal)) != 1:
            continue
        first = next(value for value in normal if value != 0)
        if first > 0:
            l_normal.append(normal)
    return l_normal


def calc_cross_distances(fract_xyz, ucp, normal_index, cuts, flag_atoms=None):
    """Shortest distances between atoms of the same layer and of different layers.

    Parameters
    ----------
    fract_xyz : array (3, n)
        All atoms of the unit cell.
    ucp : array (6,)
        Unit cell parameters (angstrom, radians).
    normal_index : array (3,)
        Miller indices of the layer planes.
    cuts : array
        Cutting planes (calc_layer_cuts).
    flag_atoms : array (n,) of bool, optional
        Only these atoms (the magnetic ones, for example).

    Returns
    -------
    d_intra, d_inter : float
        Shortest distance inside a layer and across a gap, angstrom
        (inf when there is no such pair).
    """
    if flag_atoms is not None:
        fract_xyz = fract_xyz[:, flag_atoms]
    m_m = calc_layer_normal(normal_index, ucp)[2]
    shifts = numpy.array(list(itertools.product((-1, 0, 1), repeat=3)), dtype=float).transpose()
    images = (fract_xyz[:, :, na] + shifts[:, na, :]).reshape(3, -1)
    heights_0 = numpy.asarray(normal_index) @ fract_xyz
    heights_1 = numpy.asarray(normal_index) @ images
    layer_0 = calc_layer_index(heights_0, cuts)
    layer_1 = calc_layer_index(heights_1, cuts)
    diff = m_m @ (images[:, na, :] - fract_xyz[:, :, na]).reshape(3, -1)
    dist = numpy.linalg.norm(diff, axis=0).reshape(fract_xyz.shape[1], -1)
    same = layer_0[:, na] == layer_1[na, :]
    dist = numpy.where(dist < 1e-4, numpy.inf, dist)
    d_intra = float(numpy.min(numpy.where(same, dist, numpy.inf))) if dist.size else numpy.inf
    d_inter = float(numpy.min(numpy.where(same, numpy.inf, dist))) if dist.size else numpy.inf
    return d_intra, d_inter


def calc_layer_candidates(crystal: Crystal, max_index: int = 2, magnetic_labels=None) -> list:
    """Directions of the crystal with weak bonds between layers, best first.

    For every primitive normal N up to max_index the atoms of the unit
    cell are cut into layers at the widest empty gap along N (all gaps of
    the same width are cut: equivalent layers). A layered direction has a
    wide gap and its shortest bond across the gap is clearly longer than
    the shortest bond inside a layer.

    Returns
    -------
    list of dict
        "normal" (h, k, l), "d_n" (angstrom), "n_layer" (layers in the cell),
        "gap" (empty width, angstrom), "d_intra", "d_inter" (shortest
        distances inside and across, angstrom), "ratio" d_inter/d_intra,
        and for the magnetic atoms "m_intra", "m_inter".
    """
    atoms = crystal.calc_atoms_in_unit_cell()
    fract_xyz = atoms["fract_xyz"]
    ucp = numpy.array(crystal.get_dictionary()["unit_cell_parameters"], dtype=float)
    flag_magnetic = None
    if magnetic_labels:
        labels = numpy.array(crystal.atom_site.label, dtype=str)[atoms["atom_index"]]
        flag_magnetic = numpy.isin(labels, list(magnetic_labels))
        if not numpy.any(flag_magnetic):
            flag_magnetic = None
    l_res = []
    for normal in form_candidate_normals(max_index):
        d_n = calc_layer_normal(normal, ucp)[1]
        heights = numpy.asarray(normal) @ fract_xyz
        cuts, gaps = calc_layer_cuts(heights, d_n)
        d_intra, d_inter = calc_cross_distances(fract_xyz, ucp, normal, cuts)
        d_res = {"normal": normal, "d_n": d_n, "n_layer": len(cuts), "gap": float(gaps[0]),
                 "d_intra": d_intra, "d_inter": d_inter,
                 "ratio": d_inter / d_intra if d_intra > 0. else 0.}
        if flag_magnetic is not None:
            d_res["m_intra"], d_res["m_inter"] = calc_cross_distances(
                fract_xyz, ucp, normal, cuts, flag_magnetic)
        l_res.append(d_res)
    # weak bonds across first, then the wider gap, then the lower index
    l_res.sort(key=lambda d: (-round(d["d_inter"], 2), -round(d["gap"], 2),
                              sum(abs(i) for i in d["normal"])))
    return l_res


def take_magnetic_labels(crystal: Crystal) -> list:
    """Labels of the magnetic atoms of the crystal (moments, atom_site_scat or magnetic elements)."""
    from cryspy.D_functions_item_loop.function_1_magnetic_ions import guess_magnetic_atoms
    if crystal.is_attribute("atom_site_moment"):
        return list(crystal.atom_site_moment.label)
    return list(guess_magnetic_atoms(crystal))


def suggest_layer_direction(rcif_object: GlobalN, phase_label: str = "",
                            max_index: int = 2):
    """Suggest the planes of the layers with weak bonds between them from the crystal structure.

    All atoms of the unit cell are generated by symmetry. For every
    direction N = (h k l) (primitive, |h|, |k|, |l| up to max_index) the
    cell is cut into layers at the widest empty gap between the atomic
    planes normal to N. The best direction has the longest shortest
    distance between atoms of neighbouring layers (weak, van der Waals
    like bonds) and a wide gap; a ratio d_inter / d_intra clearly above 1
    means that the bonds across the gap are weaker than inside the
    layers. For the magnetic atoms the shortest magnetic distances inside
    and between the layers are given too: a large ratio means
    two-dimensional magnetic exchange. These are geometric hints, the
    bonding itself (charges, orbitals) is not calculated.

    Parameters
    ----------
    rcif_object : GlobalN
        The object with the crystals.
    phase_label : str
        Data name of the crystal; empty for the first crystal.
    max_index : int
        Largest Miller index of the candidate directions.

    Returns
    -------
    str
        Table of the best directions in markdown and the suggestion.
    """
    crystal = take_crystal_by_label(rcif_object, phase_label)
    if crystal is None:
        return "No crystal" + (f" '{phase_label:}'." if phase_label.strip() else ".")
    magnetic_labels = take_magnetic_labels(crystal)
    l_res = calc_layer_candidates(crystal, int(max_index), magnetic_labels)
    if len(l_res) == 0:
        return "No candidate directions."
    flag_mag = "m_intra" in l_res[0]
    ls_out = [f"# Layers of {crystal.data_name:}", "",
              "| planes (h k l) | layers in cell | gap, Å | shortest bond inside, Å | "
              "shortest bond across, Å | across / inside |"
              + (" magnetic inside, Å | magnetic across, Å |" if flag_mag else ""),
              "|---|---|---|---|---|---|" + ("---|---|" if flag_mag else "")]
    for d_res in l_res[:6]:
        line = (f"| ({' '.join(str(i) for i in d_res['normal'])}) | {d_res['n_layer']} | "
                f"{d_res['gap']:.2f} | {d_res['d_intra']:.3f} | {d_res['d_inter']:.3f} | "
                f"{d_res['ratio']:.2f} |")
        if flag_mag:
            line += f" {d_res['m_intra']:.3f} | {d_res['m_inter']:.3f} |"
        ls_out.append(line)
    best = l_res[0]
    s_normal = " ".join(str(i) for i in best["normal"])
    ls_out.append("")
    if best["ratio"] < 1.1:
        ls_out.append(
            f"No layered direction: the shortest bonds across every gap are as short as "
            f"inside the layers (best ({s_normal:}), ratio {best['ratio']:.2f}). The "
            "structure is three-dimensional for the bonds.")
    else:
        ls_out.append(
            f"Suggested layers: planes ({s_normal:}), normal along the reciprocal vector "
            f"{s_normal:}, {best['n_layer']} layer(s) per cell, gap {best['gap']:.2f} Å; the "
            f"shortest bond across the gap {best['d_inter']:.3f} Å is {best['ratio']:.2f} "
            f"times the shortest bond inside ({best['d_intra']:.3f} Å).")
    if flag_mag and numpy.isfinite(best["m_intra"]):
        ratio_m = best["m_inter"] / best["m_intra"] if best["m_intra"] > 0 else 0.
        ls_out.append(
            f"Magnetic atoms ({', '.join(magnetic_labels)}): shortest distance inside a "
            f"layer {best['m_intra']:.3f} Å, between layers {best['m_inter']:.3f} Å "
            f"(ratio {ratio_m:.2f}" + ("; two-dimensional magnetic correlations are "
                                       "likely)." if ratio_m > 1.3 else ")."))
    return "\n".join(ls_out)


def take_suggested_normal(crystal: Crystal) -> tuple:
    """The normal of the best layers of the crystal (suggest_layer_direction)."""
    return tuple(calc_layer_candidates(crystal)[0]["normal"])


def parse_normal(normal: str):
    """'0 1 0' or '(0,1,0)' -> (0, 1, 0); '' -> None."""
    text = str(normal).replace("(", " ").replace(")", " ").replace(",", " ").split()
    if len(text) == 0:
        return None
    if len(text) != 3:
        raise ValueError(f"The normal '{normal:}' must be three integers h k l.")
    values = tuple(int(round(float(value))) for value in text)
    if values == (0, 0, 0):
        return None
    if numpy.gcd.reduce(numpy.abs(values)) != 1:
        divisor = int(numpy.gcd.reduce(numpy.abs(values)))
        values = tuple(value // divisor for value in values)
    return values


def add_layer_correlation(
        rcif_object: GlobalN, length: float = 100., normal: str = "",
        contribution: str = "magnetic", phase_label: str = "",
        experiment: str = ""):
    """Make the order of a powder phase two-dimensional: layers without correlations between them (Warren line shape).

    The order of the phase given by contribution exists only inside the
    layers normal to the reciprocal vector N = (h k l): its Bragg
    reflections are replaced by rods (hk), whose powder average is the
    asymmetric saw-tooth line of Warren (sharp front at the 2theta of the
    two-dimensional reflection, long tail to high angles). For "magnetic"
    the magnetic Bragg intensity of the phase is taken away and the
    magnetic rods are calculated from the ordered moments of the crystal
    (they must be given, refinable as usual); for "nuclear" the nuclear
    reflections off the 00l line are replaced by nuclear rods. The
    in-plane correlation length and the scale of the rods (start: the
    scale of the phase) are refined. The rods are convolved with the
    instrumental profile, asymmetry p1-p4 included. Only the 1D powder
    experiments at constant wavelength (pd) with neutrons are supported.

    Parameters
    ----------
    rcif_object : GlobalN
        The object with the crystals and the powder experiments.
    length : float
        Start value of the correlation length in the plane of the layers, angstrom.
    normal : str
        Miller indices of the layer planes as "h k l"; empty for the
        planes suggested from the crystal structure.
    contribution : str
        The two-dimensional order: magnetic, nuclear, or nuclear and magnetic.
    phase_label : str
        Label of the phase; empty for every phase of the experiments.
    experiment : str
        Name of the powder experiment; empty for all pd experiments.

    Returns
    -------
    str
        The phases with the layers and the parameters to refine.
    """
    contribution = contribution.strip().lower()
    if contribution not in CONTRIBUTIONS:
        raise ValueError(f"contribution must be one of {', '.join(CONTRIBUTIONS):}.")
    l_experiment = take_powder_1d_experiments(rcif_object, experiment)
    if len(l_experiment) == 0:
        if any(isinstance(item, (Pd2d, TOF)) for item in rcif_object.items):
            raise ValueError("The layer correlations are implemented only for the 1D "
                             "powder experiments at constant wavelength (pd).")
        return "No powder experiments" + (f" named '{experiment:}'." if experiment else ".")
    normal_index = parse_normal(normal)
    ls_out = []
    for pd in l_experiment:
        l_label = [label for label in pd.phase.label
                   if phase_label.strip() == "" or label.lower() == phase_label.strip().lower()]
        if len(l_label) == 0:
            raise ValueError(f"{pd.data_name:}: no phase '{phase_label:}'.")
        for label in l_label:
            crystal = take_crystal_by_label(rcif_object, label)
            if crystal is None:
                raise ValueError(f"{pd.data_name:}: no crystal '{label:}'.")
            if contribution != "nuclear" and not crystal.is_attribute("atom_site_moment"):
                raise ValueError(
                    f"Crystal '{label:}' has no ordered moments (atom_site_moment): give the "
                    "magnetic structure of the layers first, or use contribution nuclear.")
            normal_phase = normal_index or take_suggested_normal(crystal)
            s_line = set_layer_correlation_of_phase(
                pd, label, normal_phase, contribution, float(length))
            if normal_index is None:
                s_line += " The planes are suggested from the crystal structure " \
                    "(suggest_layer_direction)."
            ls_out.append(s_line)
    return "\n".join(ls_out)


def set_layer_correlation_of_phase(pd, label: str, normal_index, contribution: str,
                                   length: float) -> str:
    """Write the row of the phase to the loop layer_correlation of pd."""
    if not pd.is_attribute("layer_correlation"):
        pd.add_items([LayerCorrelationL()])
    loop = pd.layer_correlation
    row = loop.take_row(label)
    if row is not None:
        loop.items.remove(row)
    scale = 1.
    for item in pd.phase.items:
        if item.label.lower() == label.lower() and item.is_attribute("scale"):
            scale = float(item.scale)
    row = LayerCorrelation(
        phase_label=label, normal_h=int(normal_index[0]), normal_k=int(normal_index[1]),
        normal_l=int(normal_index[2]), contribution=contribution, length=length, scale=scale)
    row.length_refinement = True
    row.scale_refinement = True
    loop.items.append(row)
    s_normal = " ".join(str(int(i)) for i in normal_index)
    return (f"{pd.data_name:}, phase {label:}: {contribution:} order in the layers ({s_normal:}), "
            f"correlation length {length:.1f} Å and scale {scale:.5f} (refined).")


def report_layer_correlation(rcif_object: GlobalN, experiment: str = ""):
    """Report the layers without correlations between them (Warren line shape) of the powder phases.

    For every row of the loop layer_correlation: the layer planes, the
    two-dimensional order, the in-plane correlation length with its error,
    the scale of the rods and its ratio to the scale of the phase (the
    part of the phase in the layers when the whole phase has the same
    order), the layers of the unit cell and the positions of the first
    rods (2theta of the two-dimensional reflections).

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
        if not pd.is_attribute("layer_correlation"):
            continue
        wavelength = float(pd.setup.wavelength)
        ls_out.append(f"# {pd.data_name:}")
        for row in pd.layer_correlation.items:
            crystal = take_crystal_by_label(rcif_object, row.phase_label)
            if crystal is None:
                ls_out.append(f"- {row.phase_label:}: no crystal with this label.")
                continue
            ls_out.append(form_report_of_row(pd, row, crystal, wavelength))
        ls_out.append("")
    if len(ls_out) == 0:
        return "No pd experiment has layers without correlations (layer_correlation)."
    return "\n".join(ls_out).strip()


def form_report_of_row(pd, row, crystal, wavelength: float) -> str:
    """Lines of the report for one phase."""
    normal_index = row.get_normal_index()
    s_normal = " ".join(str(i) for i in normal_index)
    ucp = numpy.array(crystal.get_dictionary()["unit_cell_parameters"], dtype=float)
    d_n = calc_layer_normal(normal_index, ucp)[1]
    n_layer = calc_layers_of_atoms(crystal.calc_atoms_in_unit_cell()["fract_xyz"],
                                   normal_index, ucp)[2]
    s_sigma = f" ± {row.length_sigma:.1f}" if row.length_sigma else ""
    s_scale_sigma = f" ± {row.scale_sigma:.5f}" if row.scale_sigma else ""
    s_line = (f"- {row.phase_label:}: {row.contribution:} order in the layers ({s_normal:}) "
              f"(distance between the planes {d_n:.3f} Å, {n_layer} layer(s) per cell); "
              f"in-plane correlation length L = {row.length:.1f}{s_sigma} Å; scale "
              f"{row.scale:.5f}{s_scale_sigma}")
    for item in pd.phase.items:
        if item.label.lower() == row.phase_label.lower() and item.is_attribute("scale") \
                and float(item.scale) != 0.:
            s_line += f" ({row.scale / float(item.scale):.3f} of the scale of the phase)"
    s_max = 2. / wavelength
    g_rod = numpy.unique(numpy.round(calc_rods(normal_index, ucp, s_max)[1], 6))
    l_tth = [f"{numpy.degrees(2. * numpy.arcsin(0.5 * g * wavelength)):.2f}"
             for g in g_rod[:6]]
    if l_tth:
        s_line += f". First rods (fronts of the lines) at 2theta {', '.join(l_tth):} deg."
    return s_line
