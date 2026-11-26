"""Extra broadening of the magnetic reflections of the phases of powder experiments.

The magnetic intensity and the nuclear-magnetic interference of a phase
are convolved with an extra Lorentzian, whose half width in Q is the
parameter width (1/xi_mag - 1/xi_nucl, isotropic or depending on the
direction of hkl); the nuclear intensity keeps the profile of the phase
(FullProf: a "magnetic phase" with its own size). width = 0 gives the
magnetic reflections the profile of the nuclear ones. The numerical part
is in A_functions_base/magnetic_broadening.py; the procedures add the
broadening to the experiments and report it. They are shown in the menus
of CrysPy Editor and used by AI assistants.
"""
import logging

import numpy

from cryspy.A_functions_base.magnetic_broadening import (
    calc_aniso_by_isotropic_width,
    calc_aniso_from_independent,
    calc_independent_aniso_components,
    calc_principal_widths,
    calc_width_hkl,
    calc_y_mag_by_width,
)
from cryspy.B_parent_classes.cl_4_global import GlobalN
from cryspy.C_item_loop_classes.cl_1_magnetic_broadening import (
    ANISO_NAMES,
    MagneticBroadening,
    MagneticBroadeningL,
)
from cryspy.E_data_classes.cl_1_crystal import Crystal
from cryspy.E_data_classes.cl_2_pd import Pd
from cryspy.E_data_classes.cl_2_pd2d import Pd2d
from cryspy.E_data_classes.cl_2_tof import TOF
from cryspy.procedure_rhochi.rhochi_absorption import calc_laue_rotations

logger = logging.getLogger(__name__)
na = numpy.newaxis

AXIS_NAMES = ("a*", "b*", "c*")


def take_magnetic_broadening_of_phase(dict_pd: dict, phase_name: str):
    """Row of the phase in the loop magnetic_broadening of the experiment.

    Returns
    -------
    dict or None
        None when the phase has no extra magnetic broadening (no row, or
        all its widths are 0), otherwise the keys "index" (row),
        "is_aniso", "width" (inverse angstrom) and "aniso" (6,) (inverse
        angstrom).
    """
    if "magnetic_broadening_phase_name" not in dict_pd:
        return None
    l_index = numpy.flatnonzero(dict_pd["magnetic_broadening_phase_name"] == phase_name)
    if l_index.size == 0:
        return None
    index = int(l_index[0])
    is_aniso = bool(dict_pd["magnetic_broadening_is_aniso"][index])
    width = float(dict_pd["magnetic_broadening_width"][index])
    aniso = numpy.array(dict_pd["magnetic_broadening_aniso"][:, index], dtype=float)
    flag_refined = bool(numpy.any(dict_pd["flags_magnetic_broadening_aniso"][:, index])) \
        if is_aniso else bool(dict_pd["flags_magnetic_broadening_width"][index])
    # zero width: the nuclear profile, the same calculation as without the row
    if not flag_refined and ((is_aniso and not numpy.any(aniso != 0.))
                             or (not is_aniso and width == 0.)):
        return None
    return {"index": index, "is_aniso": is_aniso, "width": width, "aniso": aniso}


def take_laue_rotations(dict_crystal: dict, dict_in_out_phase: dict = None):
    """Rotations of the Laue class of the phase (kept in dict_in_out_phase)."""
    if dict_in_out_phase is not None and "magnetic_broadening_rotations" in dict_in_out_phase:
        return dict_in_out_phase["magnetic_broadening_rotations"]
    rotations = calc_laue_rotations(dict_crystal)
    if dict_in_out_phase is not None:
        dict_in_out_phase["magnetic_broadening_rotations"] = rotations
    return rotations


def form_y_mag_by_hkl(dict_pd: dict, dict_crystal: dict, dict_in_out: dict,
                      dict_in_out_phase: dict, phase_name: str,
                      unit_cell_parameters, wavelength):
    """Function giving the extra Lorentzian size term of the magnetic reflections.

    Returns
    -------
    function or None
        None when the phase has no extra magnetic broadening; otherwise
        y_mag_by_hkl(index_hkl) -> array (n_hkl,) of Y_mag in degrees (the
        FWHM of the extra Lorentzian in 2theta is Y_mag / cos(theta)). The
        anisotropic tensor is symmetrized by the Laue class of the phase and
        written to dict_in_out["magnetic_broadening_aniso_symmetrized"], so
        the dependent components follow the refined ones.
    """
    row = take_magnetic_broadening_of_phase(dict_pd, phase_name)
    if row is None:
        return None
    wavelength = float(numpy.atleast_1d(wavelength)[0])
    if not row["is_aniso"]:
        y_mag = float(calc_y_mag_by_width(row["width"], wavelength))
        dict_in_out_phase["magnetic_broadening_y"] = y_mag
        return lambda index_hkl: numpy.full(numpy.shape(index_hkl)[1], y_mag, dtype=float)

    rotations = take_laue_rotations(dict_crystal, dict_in_out_phase)
    projection, independent = calc_independent_aniso_components(
        rotations, unit_cell_parameters)
    aniso = calc_aniso_from_independent(row["aniso"], projection, independent)
    key = "magnetic_broadening_aniso_symmetrized"
    if key not in dict_in_out:
        dict_in_out[key] = numpy.array(dict_pd["magnetic_broadening_aniso"], dtype=float)
    dict_in_out[key][:, row["index"]] = aniso

    def y_mag_by_hkl(index_hkl):
        return calc_y_mag_by_width(
            calc_width_hkl(index_hkl, unit_cell_parameters, aniso), wavelength)
    return y_mag_by_hkl


def take_powder_1d_experiments(rcif_object: GlobalN, experiment: str = "") -> list:
    """Powder experiments at constant wavelength (pd), or the one with the given name."""
    name = experiment.strip().lower()
    return [item for item in rcif_object.items if isinstance(item, Pd)
            and (name == "" or item.data_name.lower() == name)]


def take_crystal_by_label(rcif_object: GlobalN, phase_label: str):
    """Crystal with the data name phase_label."""
    for item in rcif_object.items:
        if isinstance(item, Crystal) and \
                item.data_name.lower() == phase_label.lower():
            return item
    return None


def add_magnetic_broadening(
        rcif_object: GlobalN, width: float = 0., phase_label: str = "",
        anisotropic: bool = False, experiment: str = ""):
    """Give the magnetic reflections of a phase their own extra width (magnetic domain size).

    The magnetic intensity |F_M|^2 and the nuclear-magnetic interference
    of the phase are convolved with an extra Lorentzian, whose half width
    in Q is width = 1/xi_mag - 1/xi_nucl (inverse angstrom, xi are the
    magnetic and nuclear correlation lengths); the nuclear intensity
    |F_N|^2 keeps the profile of the phase. width = 0 gives the magnetic
    reflections the nuclear profile (as without it), width > 0 means
    shorter magnetic correlations (smaller magnetic domains), width < 0
    longer ones. In 2theta the full width of the Lorentzian is
    180 lambda width / (pi^2 cos(theta)) degrees, the size term Y of
    FullProf. With anisotropic the width depends on the direction of hkl
    (width_ij in inverse angstrom, convention of U_ij) restricted by the
    symmetry of the crystal class: only the independent components are
    refined. The refinement flags are switched on. Only the 1D powder
    experiments at constant wavelength (pd) are supported.

    Parameters
    ----------
    rcif_object : GlobalN
        The object with the crystals and the powder experiments.
    width : float
        Start value of the extra half width of the magnetic reflections in
        Q, inverse angstrom (0 for the nuclear width).
    phase_label : str
        Label of the magnetic phase; empty for every phase of the experiments.
    anisotropic : bool
        Width depending on the direction of hkl (as for magnetic
        correlations shorter along c).
    experiment : str
        Name of the powder experiment; empty for all pd experiments.

    Returns
    -------
    str
        The phases with the parameters to refine.

    Needs
    -----
    phase : label = phase_label
        The magnetic phase in the pd experiments.

    Adds
    ----
    magnetic_broadening : phase_label = phase_label
        The extra width of the magnetic reflections of the phase (refined).
    """
    width = float(width)
    l_experiment = take_powder_1d_experiments(rcif_object, experiment)
    if len(l_experiment) == 0:
        if any(isinstance(item, (Pd2d, TOF)) for item in rcif_object.items):
            raise ValueError("The magnetic broadening is implemented only for the 1D "
                             "powder experiments at constant wavelength (pd).")
        return "No powder experiments" + (f" named '{experiment:}'." if experiment else ".")
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
            ls_out.append(set_magnetic_broadening_of_phase(
                pd, crystal, label, width, bool(anisotropic)))
    return "\n".join(ls_out)


def set_magnetic_broadening_of_phase(pd, crystal, label: str, width: float,
                                     anisotropic: bool) -> str:
    """Write the row of the phase to the loop magnetic_broadening of pd."""
    if not pd.is_attribute("magnetic_broadening"):
        pd.add_items([MagneticBroadeningL()])
    loop = pd.magnetic_broadening
    row = loop.take_row(label)
    if row is not None:
        loop.items.remove(row)
    row = MagneticBroadening(phase_label=label)
    if not anisotropic:
        row.width = width
        row.width_refinement = True
        loop.items.append(row)
        return (f"{pd.data_name:}, phase {label:}: extra magnetic width "
                f"{width:.5f} 1/A (refined).")
    dict_crystal = crystal.get_dictionary()
    unit_cell_parameters = numpy.array(dict_crystal["unit_cell_parameters"], dtype=float)
    rotations = calc_laue_rotations(dict_crystal)
    projection, independent = calc_independent_aniso_components(rotations, unit_cell_parameters)
    aniso = calc_aniso_by_isotropic_width(width, unit_cell_parameters)
    aniso = calc_aniso_from_independent(aniso, projection, independent)
    l_refined = []
    for i_comp, name in enumerate(ANISO_NAMES):
        value = float(aniso[i_comp])
        setattr(row, name, 0. if abs(value) < 1e-12 else value)
        setattr(row, f"{name:}_refinement", bool(independent[i_comp]))
        setattr(row, f"{name:}_constraint", not bool(independent[i_comp]))
        if independent[i_comp]:
            l_refined.append(name)
    loop.items.append(row)
    return (f"{pd.data_name:}, phase {label:}: extra magnetic width depending on hkl, "
            f"start {width:.5f} 1/A; refined {', '.join(l_refined):} (the other "
            "components follow the symmetry).")


def report_magnetic_broadening(rcif_object: GlobalN, experiment: str = ""):
    """Report the extra broadening of the magnetic reflections (magnetic domain size) of the powder phases.

    For every phase of the loop magnetic_broadening: the extra half width
    in Q with its error, the Lorentzian size term Y_mag (degrees), the
    Lorentzian size term of the nuclear profile (instrument and phase) and
    the limit of a negative width, the equivalent extra Scherrer size
    2/width for a positive width, and for the width depending on hkl the
    widths along a*, b*, c* and the principal widths with their axes
    (Cartesian x || a*, z || c).

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
        if not pd.is_attribute("magnetic_broadening"):
            continue
        wavelength = float(pd.setup.wavelength)
        ls_out.append(f"# {pd.data_name:}")
        for row in pd.magnetic_broadening.items:
            crystal = take_crystal_by_label(rcif_object, row.phase_label)
            if crystal is None:
                ls_out.append(f"- {row.phase_label:}: no crystal with this label.")
                continue
            y_nucl = take_nuclear_lorentz_y(pd, row.phase_label)
            ls_out.append(form_report_of_row(row, crystal, wavelength, y_nucl))
        ls_out.append("")
    if len(ls_out) == 0:
        return "No pd experiment has an extra magnetic broadening (magnetic_broadening)."
    return "\n".join(ls_out).strip()


def take_nuclear_lorentz_y(pd, phase_label: str) -> float:
    """Lorentzian size term Y (degrees) of the nuclear profile: instrument plus phase."""
    y_nucl = float(pd.pd_instr_resolution.y) if pd.pd_instr_resolution.is_attribute("y") else 0.
    for item in pd.phase.items:
        if item.label.lower() == phase_label.lower() and item.is_attribute("y"):
            y_nucl += float(item.y)
    return y_nucl


def form_width_text(width: float, sigma: float, wavelength: float) -> str:
    """Width with its error, Y_mag and the extra Scherrer size."""
    s_sigma = f" ± {sigma:.5f}" if sigma else ""
    y_mag = float(calc_y_mag_by_width(width, wavelength))
    s_out = f"{width:.5f}{s_sigma:} 1/A (Y_mag = {y_mag:.4f} deg"
    if width > 0.:
        s_out += f", extra Scherrer size 2/width = {2. / width:.0f} A"
    return s_out + ")"


def form_report_of_row(row, crystal, wavelength: float, y_nucl: float) -> str:
    """Lines of the report for one phase."""
    width_min = -numpy.square(numpy.pi) * y_nucl / (180. * wavelength)
    s_limit = (f" The Lorentzian term of the nuclear profile is Y = {y_nucl:.4f} deg: a "
               f"width below {width_min:.5f} 1/A cannot make the magnetic peaks narrower.")
    if not row.is_anisotropic():
        width = float(row.width) if row.is_attribute("width") else 0.
        s_line = f"- {row.phase_label:}: extra magnetic half width in Q "
        s_line += form_width_text(width, float(row.width_sigma), wavelength) + "."
        if width == 0.:
            s_line += " The magnetic reflections have the nuclear profile."
        elif width > 0.:
            s_line += " The magnetic correlations are shorter than the nuclear ones."
        else:
            s_line += " The magnetic correlations are longer than the nuclear ones." + s_limit
        return s_line
    unit_cell_parameters = numpy.array(crystal.get_dictionary()["unit_cell_parameters"],
                                       dtype=float)
    aniso = numpy.array(row.get_aniso(), dtype=float)
    l_axis = [f"{axis_name:} {aniso[i_axis]:.5f}" for i_axis, axis_name in enumerate(AXIS_NAMES)]
    widths, axes = calc_principal_widths(aniso, unit_cell_parameters)
    l_principal = [f"{width:.5f} along ({axes[0, i]:.3f}, {axes[1, i]:.3f}, {axes[2, i]:.3f})"
                   for i, width in enumerate(widths)]
    s_line = (f"- {row.phase_label:}: extra magnetic half width in Q (1/A) along "
              f"{', '.join(l_axis):}; principal widths {'; '.join(l_principal):} "
              "(Cartesian x || a*, z || c). A larger width means shorter magnetic "
              "correlations along this direction.")
    if widths[0] < 0.:
        s_line += s_limit
    return s_line
