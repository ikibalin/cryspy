"""Finite magnetic correlation length of the phases of powder experiments.

The magnetic intensity and the nuclear-magnetic interference of a phase
are convolved with an extra Lorentzian, the half width of which in Q is
1/xi (isotropic) or 1/xi(h) for an ellipsoid of correlation lengths; the
nuclear intensity keeps the profile of the phase (FullProf: a "magnetic
phase" with its own size). The numerical part is in
A_functions_base/magnetic_correlation.py; the procedures add the
correlation length to the experiments and report it. They are shown in the
menus of CrysPy Editor and used by AI assistants.
"""
import logging

import numpy

from cryspy.A_functions_base.magnetic_correlation import (
    calc_aniso_by_isotropic_length,
    calc_aniso_from_independent,
    calc_correlation_length_hkl,
    calc_independent_aniso_components,
    calc_principal_correlation_lengths,
    calc_y_mag_by_correlation_length,
)
from cryspy.B_parent_classes.cl_4_global import GlobalN
from cryspy.C_item_loop_classes.cl_1_magnetic_correlation import (
    ANISO_NAMES,
    MagneticCorrelation,
    MagneticCorrelationL,
)
from cryspy.E_data_classes.cl_1_crystal import Crystal
from cryspy.E_data_classes.cl_2_pd import Pd
from cryspy.E_data_classes.cl_2_pd2d import Pd2d
from cryspy.E_data_classes.cl_2_tof import TOF
from cryspy.procedure_rhochi.rhochi_absorption import calc_laue_rotations

logger = logging.getLogger(__name__)
na = numpy.newaxis

AXIS_NAMES = ("a*", "b*", "c*")


def take_magnetic_correlation_of_phase(dict_pd: dict, phase_name: str):
    """Row of the phase in the loop magnetic_correlation of the experiment.

    Returns
    -------
    dict or None
        None when the phase has no magnetic correlation length (or it is
        infinite), otherwise the keys "index" (row), "is_aniso", "length"
        (angstrom) and "aniso" (6,) (angstrom^2).
    """
    if "magnetic_correlation_phase_name" not in dict_pd:
        return None
    l_index = numpy.flatnonzero(dict_pd["magnetic_correlation_phase_name"] == phase_name)
    if l_index.size == 0:
        return None
    index = int(l_index[0])
    is_aniso = bool(dict_pd["magnetic_correlation_is_aniso"][index])
    length = float(dict_pd["magnetic_correlation_length"][index])
    aniso = numpy.array(dict_pd["magnetic_correlation_aniso"][:, index], dtype=float)
    if not is_aniso and length <= 0.:
        return None
    return {"index": index, "is_aniso": is_aniso, "length": length, "aniso": aniso}


def take_laue_rotations(dict_crystal: dict, dict_in_out_phase: dict = None):
    """Rotations of the Laue class of the phase (kept in dict_in_out_phase)."""
    if dict_in_out_phase is not None and "magnetic_correlation_rotations" in dict_in_out_phase:
        return dict_in_out_phase["magnetic_correlation_rotations"]
    rotations = calc_laue_rotations(dict_crystal)
    if dict_in_out_phase is not None:
        dict_in_out_phase["magnetic_correlation_rotations"] = rotations
    return rotations


def form_y_mag_by_hkl(dict_pd: dict, dict_crystal: dict, dict_in_out: dict,
                      dict_in_out_phase: dict, phase_name: str,
                      unit_cell_parameters, wavelength):
    """Function giving the extra Lorentzian size term of the magnetic reflections.

    Returns
    -------
    function or None
        None when the phase has no magnetic correlation length; otherwise
        y_mag_by_hkl(index_hkl) -> array (n_hkl,) of Y_mag in degrees
        (FWHM of the Lorentzian in 2theta is Y_mag / cos(theta)). For the
        ellipsoid the tensor is symmetrized by the Laue class of the phase
        and written to dict_in_out["magnetic_correlation_aniso_symmetrized"],
        so the dependent components follow the refined ones.
    """
    row = take_magnetic_correlation_of_phase(dict_pd, phase_name)
    if row is None:
        return None
    wavelength = float(numpy.atleast_1d(wavelength)[0])
    if not row["is_aniso"]:
        y_mag = float(calc_y_mag_by_correlation_length(row["length"], wavelength))
        dict_in_out_phase["magnetic_correlation_y"] = y_mag
        return lambda index_hkl: numpy.full(numpy.shape(index_hkl)[1], y_mag, dtype=float)

    rotations = take_laue_rotations(dict_crystal, dict_in_out_phase)
    projection, independent = calc_independent_aniso_components(
        rotations, unit_cell_parameters)
    aniso = calc_aniso_from_independent(row["aniso"], projection, independent)
    key = "magnetic_correlation_aniso_symmetrized"
    if key not in dict_in_out:
        dict_in_out[key] = numpy.array(dict_pd["magnetic_correlation_aniso"], dtype=float)
    dict_in_out[key][:, row["index"]] = aniso

    def y_mag_by_hkl(index_hkl):
        xi_hkl = calc_correlation_length_hkl(index_hkl, unit_cell_parameters, aniso)
        return calc_y_mag_by_correlation_length(xi_hkl, wavelength)
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


def take_unit_cell_parameters(crystal) -> numpy.ndarray:
    """a, b, c (angstrom), alpha, beta, gamma (radians) of the crystal."""
    return numpy.array(crystal.get_dictionary()["unit_cell_parameters"], dtype=float)


def add_magnetic_correlation_length(
        rcif_object: GlobalN, length: float = 100., phase_label: str = "",
        anisotropic: bool = False, experiment: str = ""):
    """Give the magnetic reflections of a phase their own Lorentzian width (magnetic domain size).

    The magnetic intensity |F_M|^2 and the nuclear-magnetic interference
    of the phase are convolved with an extra Lorentzian, the half width of
    which in Q is 1/xi (xi is the magnetic correlation length); the nuclear
    intensity |F_N|^2 keeps the profile of the phase. In 2theta the full
    width of the Lorentzian is 180 lambda / (pi^2 xi cos(theta)) degrees, as
    the size term Y of FullProf for the Scherrer size 2 xi. With anisotropic
    the correlation lengths form an ellipsoid (aniso_ij in angstrom^2,
    convention of U_ij) restricted by the symmetry of the crystal class:
    only the independent components are refined. The values are set to the
    given length and their refinement flags are switched on. Only the 1D
    powder experiments at constant wavelength (pd) are supported.

    Parameters
    ----------
    rcif_object : GlobalN
        The object with the crystals and the powder experiments.
    length : float
        Start value of the magnetic correlation length xi, angstrom (about
        half of the magnetic domain size).
    phase_label : str
        Label of the magnetic phase; empty for every phase of the experiments.
    anisotropic : bool
        Ellipsoid of correlation lengths instead of one length (when the
        width depends on the direction of hkl, as for correlations shorter
        along c).
    experiment : str
        Name of the powder experiment; empty for all pd experiments.

    Returns
    -------
    str
        The phases with the parameters to refine.
    """
    length = float(length)
    if length <= 0.:
        raise ValueError("The magnetic correlation length should be positive (angstrom).")
    l_experiment = take_powder_1d_experiments(rcif_object, experiment)
    if len(l_experiment) == 0:
        if any(isinstance(item, (Pd2d, TOF)) for item in rcif_object.items):
            raise ValueError("The magnetic correlation length is implemented only for "
                             "the 1D powder experiments at constant wavelength (pd).")
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
            ls_out.append(set_magnetic_correlation_of_phase(
                pd, crystal, label, length, bool(anisotropic)))
    return "\n".join(ls_out)


def set_magnetic_correlation_of_phase(pd, crystal, label: str, length: float,
                                      anisotropic: bool) -> str:
    """Write the row of the phase to the loop magnetic_correlation of pd."""
    if not pd.is_attribute("magnetic_correlation"):
        pd.add_items([MagneticCorrelationL()])
    loop = pd.magnetic_correlation
    row = loop.take_row(label)
    if row is not None:
        loop.items.remove(row)
    row = MagneticCorrelation(phase_label=label)
    if not anisotropic:
        row.length = length
        row.length_refinement = True
        loop.items.append(row)
        return (f"{pd.data_name:}, phase {label:}: magnetic correlation length "
                f"{length:.1f} A (refined).")
    dict_crystal = crystal.get_dictionary()
    unit_cell_parameters = numpy.array(dict_crystal["unit_cell_parameters"], dtype=float)
    rotations = calc_laue_rotations(dict_crystal)
    projection, independent = calc_independent_aniso_components(rotations, unit_cell_parameters)
    aniso = calc_aniso_by_isotropic_length(length, unit_cell_parameters)
    aniso = calc_aniso_from_independent(aniso, projection, independent)
    l_refined = []
    for i_comp, name in enumerate(ANISO_NAMES):
        value = float(aniso[i_comp])
        setattr(row, name, 0. if abs(value) < 1e-8 else value)
        setattr(row, f"{name:}_refinement", bool(independent[i_comp]))
        setattr(row, f"{name:}_constraint", not bool(independent[i_comp]))
        if independent[i_comp]:
            l_refined.append(name)
    loop.items.append(row)
    return (f"{pd.data_name:}, phase {label:}: ellipsoid of magnetic correlation "
            f"lengths starting from {length:.1f} A; refined {', '.join(l_refined):} "
            "(the other components follow the symmetry).")


def report_magnetic_correlation(rcif_object: GlobalN, experiment: str = ""):
    """Report the magnetic correlation lengths (magnetic domain sizes) of the powder phases.

    For every phase with the loop magnetic_correlation: the correlation
    length xi (angstrom) with its error, the Scherrer size 2 xi (the
    apparent size of FullProf), the Lorentzian size term Y_mag (degrees)
    and, for the ellipsoid, the lengths along a*, b*, c* and the principal
    lengths with their axes (Cartesian x || a*, z || c).

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
        if not pd.is_attribute("magnetic_correlation"):
            continue
        wavelength = float(pd.setup.wavelength)
        ls_out.append(f"# {pd.data_name:}")
        for row in pd.magnetic_correlation.items:
            crystal = take_crystal_by_label(rcif_object, row.phase_label)
            if crystal is None:
                ls_out.append(f"- {row.phase_label:}: no crystal with this label.")
                continue
            ls_out.append(form_report_of_row(row, crystal, wavelength))
        ls_out.append("")
    if len(ls_out) == 0:
        return "No pd experiment has a magnetic correlation length (magnetic_correlation)."
    return "\n".join(ls_out).strip()


def form_report_of_row(row, crystal, wavelength: float) -> str:
    """Lines of the report for one phase."""
    if not row.is_anisotropic():
        length = float(row.length) if row.is_attribute("length") else 0.
        if length <= 0.:
            return f"- {row.phase_label:}: infinite magnetic correlation length (no broadening)."
        s_sigma = f"({float(row.length_sigma):.1f})" if row.length_sigma else ""
        y_mag = float(calc_y_mag_by_correlation_length(length, wavelength))
        return (f"- {row.phase_label:}: xi = {length:.1f}{s_sigma:} A, Scherrer size "
                f"2 xi = {2. * length:.1f} A, half width in Q 1/xi = {1. / length:.5f} 1/A, "
                f"Y_mag = {y_mag:.4f} deg.")
    unit_cell_parameters = take_unit_cell_parameters(crystal)
    aniso = numpy.array(row.get_aniso(), dtype=float)
    l_axis = []
    for i_axis, axis_name in enumerate(AXIS_NAMES):
        value = aniso[i_axis]
        xi_axis = numpy.sqrt(value) if value > 0. else 0.
        l_axis.append(f"{axis_name:} {xi_axis:.1f} A")
    lengths, axes = calc_principal_correlation_lengths(aniso, unit_cell_parameters)
    l_principal = [f"{length:.1f} A along ({axes[0, i]:.3f}, {axes[1, i]:.3f}, {axes[2, i]:.3f})"
                   for i, length in enumerate(lengths)]
    s_line = (f"- {row.phase_label:}: ellipsoid of correlation lengths, along "
              f"{', '.join(l_axis):}; principal lengths {'; '.join(l_principal):} "
              "(Cartesian x || a*, z || c). Scherrer sizes are 2 xi.")
    if numpy.any(lengths <= 0.):
        s_line += " Not positive definite: the ellipsoid is unphysical, refine the isotropic length."
    elif lengths[-1] > 1.5 * lengths[0]:
        s_line += (f" The shortest principal length is {lengths[-1] / lengths[0]:.1f} times "
                   "shorter than the longest one.")
    return s_line
