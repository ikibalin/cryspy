"""Spherical harmonics texture of the phases of powder experiments.

The texture correction multiplies the integrated intensities of a phase in
the pd, pd2d and tof experiments (GSAS-II model, see
A_functions_base/texture_spherical_harmonics.py). The procedures add the
texture to the experiments and report it; they are shown in the menus of
CrysPy Editor and used by AI assistants.
"""
import logging

import numpy

from cryspy.A_functions_base.texture_spherical_harmonics import (
    SAMPLE_SYMMETRIES,
    calc_crystal_harmonics,
    calc_laue_class_by_rotations,
    calc_sample_harmonics,
    calc_sample_angles,
    calc_texture_index,
    calc_texture_sh_correction_by_harmonics,
    calc_texture_sh_harmonics,
    form_texture_sh_lmn,
    is_allowed_harmonic,
)
from cryspy.B_parent_classes.cl_4_global import GlobalN
from cryspy.C_item_loop_classes.cl_1_texture_sh import (
    DEFAULT_GONIOMETER_ANGLES,
    TextureSh,
    TextureShL,
    TextureShCoefficient,
    TextureShCoefficientL,
)
from cryspy.E_data_classes.cl_1_crystal import Crystal
from cryspy.E_data_classes.cl_2_pd import Pd
from cryspy.E_data_classes.cl_2_pd2d import Pd2d
from cryspy.E_data_classes.cl_2_tof import TOF
from cryspy.procedure_rhochi.rhochi_absorption import calc_laue_rotations

logger = logging.getLogger(__name__)
na = numpy.newaxis

N_GRID_POLAR = 37
N_GRID_AZIMUTH = 72


def take_laue_class_of_crystal(dict_crystal: dict, dict_in_out_phase: dict = None):
    """Laue class and unique axis of the phase for the texture harmonics.

    The result is kept in dict_in_out_phase (the symmetry does not change
    during a refinement).
    """
    if dict_in_out_phase is not None and "texture_sh_laue_class" in dict_in_out_phase:
        return dict_in_out_phase["texture_sh_laue_class"]
    rotations = calc_laue_rotations(dict_crystal)
    result = calc_laue_class_by_rotations(rotations, dict_crystal["unit_cell_parameters"])
    if dict_in_out_phase is not None:
        dict_in_out_phase["texture_sh_laue_class"] = result
    return result


def take_texture_sh_of_phase(dict_pd: dict, phase_name: str):
    """Texture of the phase from the dictionary of the experiment.

    Returns
    -------
    dict or None
        None when the phase has no spherical harmonics texture, otherwise
        the keys "index" (row of texture_sh), "sample_symmetry",
        "order_max", "sample_angles" (3,), "flag_sample_angles",
        "index_coefficient" (rows of texture_sh_coefficient), "lmn" (3, n),
        "coefficient" (n,), "flag_coefficient".
    """
    if "texture_sh_phase_name" not in dict_pd:
        return None
    l_index = numpy.flatnonzero(dict_pd["texture_sh_phase_name"] == phase_name)
    if l_index.size == 0:
        return None
    index = int(l_index[0])
    if "texture_sh_coefficient_phase_name" in dict_pd:
        index_coefficient = numpy.flatnonzero(
            dict_pd["texture_sh_coefficient_phase_name"] == phase_name)
    else:
        index_coefficient = numpy.zeros((0,), dtype=int)
    lmn = dict_pd["texture_sh_lmn"][:, index_coefficient] if index_coefficient.size > 0 \
        else numpy.zeros((3, 0), dtype=int)
    coefficient = dict_pd["texture_sh_coefficient"][index_coefficient] \
        if index_coefficient.size > 0 else numpy.zeros((0,), dtype=float)
    flag_coefficient = numpy.any(dict_pd["flags_texture_sh_coefficient"][index_coefficient]) \
        if index_coefficient.size > 0 else False
    return {
        "index": index,
        "sample_symmetry": str(dict_pd["texture_sh_sample_symmetry"][index]),
        "order_max": int(dict_pd["texture_sh_order_max"][index]),
        "sample_angles": dict_pd["texture_sh_sample_angles"][:, index],
        "flag_sample_angles": bool(numpy.any(
            dict_pd["flags_texture_sh_sample_angles"][:, index])),
        "index_coefficient": index_coefficient,
        "lmn": lmn,
        "coefficient": coefficient,
        "flag_coefficient": bool(flag_coefficient),
    }


def check_texture_sh_lmn(lmn, laue_class: str, sample_symmetry: str, order_max: int,
                         phase_name: str):
    """Raise ValueError for a coefficient not allowed by the symmetry."""
    for l, m, n in lmn.transpose():
        if l > order_max:
            raise ValueError(
                f"Texture of '{phase_name:}': the coefficient ({l:}, {m:}, {n:}) has l larger "
                f"than order_max = {order_max:}.")
        if not (is_allowed_harmonic(sample_symmetry, l, m) and
                is_allowed_harmonic(laue_class, l, n)):
            raise ValueError(
                f"Texture of '{phase_name:}': the coefficient ({l:}, {m:}, {n:}) is not allowed "
                f"by the Laue class {laue_class:} and the sample symmetry {sample_symmetry:}.")


def calc_texture_sh_for_phase(
        dict_pd: dict, dict_crystal: dict, dict_in_out_phase: dict, phase_name: str,
        index_hkl, unit_cell_parameters, theta, azimuth,
        flag_use_precalculated_data: bool = False):
    """Texture correction of the reflections of a phase in an experiment.

    Parameters
    ----------
    dict_pd : dict
        Dictionary of the experiment (pd, pd2d or tof).
    dict_crystal : dict
        Dictionary of the crystal of the phase.
    dict_in_out_phase : dict
        Dictionary of the calculated values of the phase in the experiment.
    phase_name : str
        Label of the phase (lower case).
    index_hkl : numpy.ndarray
        Miller indexes, shape (3, n_hkl).
    unit_cell_parameters : numpy.ndarray
        a, b, c (Å), alpha, beta, gamma (radians).
    theta : numpy.ndarray
        Bragg angle of the reflections, radians (shape (n_hkl,) or a
        scalar for a time-of-flight bank).
    azimuth : numpy.ndarray
        Azimuth of the detector, radians: a scalar, or (..., 1) for the
        points of a 2D detector.
    flag_use_precalculated_data : bool
        Keep the harmonics of the previous call when their arguments are
        the same.

    Returns
    -------
    numpy.ndarray or None
        Correction of the shape of the broadcast of theta and azimuth with
        the reflections on the last axis; None when the phase has no
        spherical harmonics texture.
    """
    texture = take_texture_sh_of_phase(dict_pd, phase_name)
    if texture is None:
        return None
    goniometer = numpy.asarray(dict_pd.get(
        "goniometer_angles", numpy.radians(numpy.array(DEFAULT_GONIOMETER_ANGLES))), dtype=float)
    flag_coupled = str(dict_pd.get("goniometer_geometry", "debye_scherrer")) == "bragg_brentano"
    theta = numpy.asarray(theta, dtype=float)
    if theta.ndim == 0:
        theta = numpy.full((index_hkl.shape[1],), float(theta))
    azimuth = numpy.asarray(azimuth, dtype=float)
    # the correction is linear in the coefficients: the harmonics are kept
    # while their arguments stay the same (any refined angle, cell
    # parameter or offset changes them and the harmonics are recalculated)
    l_argument = [numpy.asarray(index_hkl), numpy.asarray(unit_cell_parameters, dtype=float),
                  theta, azimuth, goniometer, numpy.asarray(texture["sample_angles"], dtype=float),
                  texture["lmn"]]
    cache = dict_in_out_phase.get("texture_sh_harmonics")
    flag_keep = (flag_use_precalculated_data and cache is not None
                 and cache["flag_coupled"] == flag_coupled
                 and cache["sample_symmetry"] == texture["sample_symmetry"]
                 and all(hh_1.shape == hh_2.shape and numpy.array_equal(hh_1, hh_2)
                         for hh_1, hh_2 in zip(cache["arguments"], l_argument)))
    if not flag_keep:
        laue_class, unique_axis = take_laue_class_of_crystal(dict_crystal, dict_in_out_phase)
        check_texture_sh_lmn(texture["lmn"], laue_class, texture["sample_symmetry"],
                             texture["order_max"], phase_name)
        sample_omega, sample_chi, sample_phi = texture["sample_angles"]
        psi, gamma = calc_sample_angles(
            theta, goniometer[0], goniometer[1], goniometer[2], goniometer[3] + azimuth,
            sample_omega, sample_chi, sample_phi, flag_coupled=flag_coupled)
        harmonics = calc_texture_sh_harmonics(
            index_hkl, unit_cell_parameters, laue_class, unique_axis,
            texture["sample_symmetry"], texture["lmn"], psi, gamma)
        cache = {"arguments": [numpy.copy(hh) for hh in l_argument],
                 "flag_coupled": flag_coupled, "sample_symmetry": texture["sample_symmetry"],
                 "harmonics": harmonics}
        dict_in_out_phase["texture_sh_harmonics"] = cache
    correction = calc_texture_sh_correction_by_harmonics(
        *cache["harmonics"], texture["coefficient"])
    dict_in_out_phase["texture_sh_correction"] = correction
    return correction


def take_powder_experiments(rcif_object: GlobalN, experiment: str = "") -> list:
    """Powder experiments (pd, pd2d, tof), or the one with the given name."""
    name = experiment.strip().lower()
    return [item for item in rcif_object.items if isinstance(item, (Pd, Pd2d, TOF))
            and (name == "" or item.data_name.lower() == name)]


def take_crystal_by_label(rcif_object: GlobalN, phase_label: str):
    """Crystal with the data name phase_label."""
    for item in rcif_object.items:
        if isinstance(item, Crystal) and \
                item.data_name.lower() == phase_label.lower():
            return item
    return None


def add_texture_spherical_harmonics(
        rcif_object: GlobalN, phase_label: str = "", order_max: int = 4,
        sample_symmetry: str = "cylindrical", experiment: str = ""):
    """Describe the texture of a phase in powder experiments by spherical harmonics.

    The intensities of the phase are multiplied by the texture correction
    of GSAS-II, 1 + sum 4pi/(2l+1) C_l^mn k_l^m(sample) k_l^n(crystal). The
    coefficients C_l^mn allowed by the Laue class and the sample symmetry
    up to order_max are added to texture_sh_coefficient with the value 0
    and the refinement flag; the existing coefficients keep their values,
    the ones above order_max are removed (and the ones with m != 0 when the
    sample symmetry changes). A phase in several experiments gets one set of
    coefficients shared by them (equal constraint marks): banks, sectors of
    a 2D detector or orientations of the same sample. The Bragg-Brentano geometry, the
    goniometer angles and the azimuth of the detector are set in
    sample_goniometer (by default Debye-Scherrer with the cylinder axis
    vertical, normal to the horizontal scattering plane).

    Parameters
    ----------
    rcif_object : GlobalN
        The object with the crystals and the powder experiments.
    phase_label : str
        Label of the textured phase; empty for every phase of the experiments.
    order_max : int
        Maximal even order of the harmonics (2 to 34; 4 to 8 for a weak texture).
    sample_symmetry : str
        Sample symmetry: 'cylindrical' (fibre texture, one pattern), '2/m',
        'mmm' (rolling) or '-1'.
    experiment : str
        Name of the powder experiment; empty for all powder experiments.

    Returns
    -------
    str
        The phases with the number of coefficients to refine.

    Adds
    ----
    texture_sh : phase_label = phase_label
        The order and the sample symmetry of the texture of the phase.
    texture_sh_coefficient : phase_label = phase_label
        The coefficients C_l^mn of the phase (refined).
    sample_goniometer
        The geometry of the sample (Debye-Scherrer by default).
    """
    sample_symmetry = sample_symmetry.strip()
    if sample_symmetry not in SAMPLE_SYMMETRIES:
        raise ValueError(f"The sample symmetry should be one of {', '.join(SAMPLE_SYMMETRIES)}.")
    order_max = int(order_max)
    if order_max < 2 or order_max > 34:
        raise ValueError("order_max should be between 2 and 34.")
    order_max -= order_max % 2
    l_experiment = take_powder_experiments(rcif_object, experiment)
    if len(l_experiment) == 0:
        return "No powder experiments" + (f" named '{experiment:}'." if experiment else ".")
    ls_out = []
    d_phase_experiments = {}
    for pd in l_experiment:
        l_label = [label for label in pd.phase.label
                   if phase_label.strip() == "" or label.lower() == phase_label.strip().lower()]
        if len(l_label) == 0:
            raise ValueError(f"{pd.data_name:}: no phase '{phase_label:}'.")
        for label in l_label:
            crystal = take_crystal_by_label(rcif_object, label)
            if crystal is None:
                raise ValueError(f"{pd.data_name:}: no crystal '{label:}'.")
            dict_crystal = crystal.get_dictionary()
            laue_class, unique_axis = take_laue_class_of_crystal(dict_crystal)
            if pd.is_attribute("texture") and any(
                    item.label.lower() == label.lower() for item in pd.texture.items
                    if item.is_attribute("label")):
                raise ValueError(f"{pd.data_name:}: phase '{label:}' has the March-Dollase "
                                 "texture; remove it before adding the spherical harmonics.")
            ls_out.append(set_texture_sh_of_phase(
                pd, label, laue_class, unique_axis, order_max, sample_symmetry))
            d_phase_experiments.setdefault(label.lower(), []).append(pd)
    for label, l_pd in d_phase_experiments.items():
        if len(l_pd) > 1:
            ls_out.append(tie_texture_sh_coefficients(l_pd, label))
    return "\n".join(ls_out)


def tie_texture_sh_coefficients(l_pd: list, label: str) -> str:
    """Give the coefficients of a phase the same values and constraint
    marks in all the experiments: the texture belongs to the sample (as the
    phase texture of GSAS-II), the experiments are banks, detector sectors
    or orientations of it."""
    d_value = {}
    for pd in l_pd:
        for item in pd.texture_sh_coefficient.take_coefficients(label):
            key = (int(item.l), int(item.m), int(item.n))
            if key not in d_value:
                d_value[key] = (item.value, f"{label:}_c{key[0]:}_{key[1]:}_{key[2]:}".replace(
                    "-", "m"))
            item.value, item.value_mark = d_value[key]
    return (f"Phase '{label:}': the texture coefficients are shared by the experiments "
            f"{', '.join(pd.data_name for pd in l_pd):} (equal marks such as "
            f"{next(iter(d_value.values()))[1] if d_value else '':}); remove the marks to "
            "refine them separately.")


def set_texture_sh_of_phase(pd, label: str, laue_class: str, unique_axis: str,
                            order_max: int, sample_symmetry: str) -> str:
    """Add or change the texture_sh row and the coefficients of one phase."""
    if not pd.is_attribute("texture_sh"):
        pd.add_items([TextureShL()])
    texture_sh = pd.texture_sh
    l_row = [item for item in texture_sh.items if item.phase_label.lower() == label.lower()]
    flag_new_symmetry = False
    if len(l_row) == 0:
        row = TextureSh(phase_label=label, order_max=order_max,
                        sample_symmetry=sample_symmetry)
        texture_sh.items.append(row)
    else:
        row = l_row[0]
        flag_new_symmetry = row.sample_symmetry != sample_symmetry
        row.sample_symmetry = sample_symmetry
        row.order_max = order_max

    if not pd.is_attribute("texture_sh_coefficient"):
        pd.add_items([TextureShCoefficientL()])
    coefficients = pd.texture_sh_coefficient
    # k_l^0 of the sample is the same for every sample symmetry; the
    # harmonics of m != 0 change their meaning with the symmetry
    d_old = {(int(item.l), int(item.m), int(item.n)): item
             for item in coefficients.take_coefficients(label)
             if not (flag_new_symmetry and int(item.m) != 0)}
    lmn = form_texture_sh_lmn(laue_class, sample_symmetry, order_max)
    l_new, n_added = [], 0
    for l, m, n in lmn.transpose():
        key = (int(l), int(m), int(n))
        if key in d_old:
            l_new.append(d_old[key])
        else:
            n_added += 1
            l_new.append(TextureShCoefficient(
                phase_label=label, l=key[0], m=key[1], n=key[2], value=0.,
                value_refinement=True))
    l_other = [item for item in coefficients.items if item.phase_label.lower() != label.lower()]
    coefficients.items.clear()
    coefficients.items.extend(l_other + l_new)
    s_axis = f", unique axis {unique_axis:}" if unique_axis else ""
    return (f"{pd.data_name:}: phase '{label:}' (Laue class {laue_class:}{s_axis:}), sample "
            f"symmetry {sample_symmetry:}, order_max {order_max:}: {len(l_new):} coefficients "
            f"C_l^mn to refine ({n_added:} new).")


def calc_inverse_pole_figure_range(laue_class: str, sample_symmetry: str, lmn, coefficient):
    """Minimum and maximum of the texture correction over the crystal
    directions for the cylinder axis (psi = 0) of the sample, in multiples
    of the random distribution."""
    polar = numpy.linspace(0., numpy.pi, N_GRID_POLAR)
    azimuth = numpy.linspace(0., 2. * numpy.pi, N_GRID_AZIMUTH, endpoint=False)
    polar, azimuth = [hh.flatten() for hh in numpy.meshgrid(polar, azimuth, indexing="ij")]
    lmn = numpy.asarray(lmn, dtype=int).reshape(3, -1)
    k_crystal = calc_crystal_harmonics(lmn[[0, 2]], laue_class, polar, azimuth)
    k_sample = calc_sample_harmonics(lmn[[0, 1]], sample_symmetry, 0., 0.)
    l_norm = 4. * numpy.pi / (2. * lmn[0] + 1.)
    value = 1. + numpy.sum((coefficient * l_norm * k_sample)[:, na] * k_crystal, axis=0)
    return float(value.min()), float(value.max())


def report_texture_spherical_harmonics(rcif_object: GlobalN, experiment: str = ""):
    """Report the spherical harmonics texture of the powder experiments.

    For every textured phase: the Laue class, the sample symmetry, the
    number of coefficients, the texture index J (1 for a random powder) and
    the range of the correction over the crystal directions along the
    cylinder axis of the sample (inverse pole figure, multiples of the
    random distribution). A negative minimum means an unphysical texture:
    lower order_max or refine fewer coefficients.

    Parameters
    ----------
    rcif_object : GlobalN
        The object with the crystals and the powder experiments.
    experiment : str
        Name of the powder experiment; empty for all powder experiments.

    Returns
    -------
    str
        Report in markdown.
    """
    ls_out = []
    for pd in take_powder_experiments(rcif_object, experiment):
        if not pd.is_attribute("texture_sh"):
            continue
        dict_pd = pd.get_dictionary()
        ls_out.append(f"# {pd.data_name:}")
        for row in pd.texture_sh.items:
            label = row.phase_label.lower()
            texture = take_texture_sh_of_phase(dict_pd, label)
            crystal = take_crystal_by_label(rcif_object, label)
            if crystal is None:
                ls_out.append(f"- {row.phase_label:}: no crystal with this label.")
                continue
            dict_crystal = crystal.get_dictionary()
            laue_class, unique_axis = take_laue_class_of_crystal(dict_crystal)
            j_index = calc_texture_index(texture["lmn"], texture["coefficient"])
            v_min, v_max = calc_inverse_pole_figure_range(
                laue_class, texture["sample_symmetry"], texture["lmn"], texture["coefficient"])
            s_line = (f"- {row.phase_label:}: Laue class {laue_class:}, sample symmetry "
                      f"{texture['sample_symmetry']:}, order_max {texture['order_max']:}, "
                      f"{texture['lmn'].shape[1]:} coefficients, texture index J = "
                      f"{j_index:.3f}, inverse pole figure of the sample axis from "
                      f"{v_min:.2f} to {v_max:.2f} MRD")
            if v_min < 0.:
                s_line += " (negative: unphysical, lower order_max)"
            ls_out.append(s_line + ".")
        ls_out.append("")
    if len(ls_out) == 0:
        return "No powder experiment has a spherical harmonics texture (texture_sh)."
    return "\n".join(ls_out).strip()

