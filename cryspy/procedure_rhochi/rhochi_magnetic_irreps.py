"""Choice of the magnetic representation of a propagation vector by the data.

For every irreducible representation of the group of k (and every mode of
the multidimensional ones) a copy of the object gets the propagation vector
(Crystal.add_propagation_vector), its Fourier components are refined against
the experiments of the object and the chi_sq values are compared, as the
representation analysis of SARAh, BasIreps or ISODISTORT followed by a
refinement of every model (FullProf, Jana2020).
"""
import logging
import warnings

import numpy

from cryspy.A_functions_base.magnetic_k_symmetry import snap_k_vector
from cryspy.B_parent_classes.cl_4_global import GlobalN
from cryspy.D_functions_item_loop.function_1_magnetic_k import (
    classify_moment_modulation,
    form_irrep_name,
    take_irrep_candidates,
)
from cryspy.E_data_classes.cl_1_crystal import Crystal
from cryspy.exceptions import CrysPyError
from cryspy.H_functions_global.function_1_cryspy_objects import str_to_globaln
from cryspy.procedure_rhochi.rhochi import rhochi_no_refinement
from cryspy.procedure_rhochi.rhochi_by_dictionary import (
    rhochi_rietveld_refinement_by_dictionary,
)

logger = logging.getLogger(__name__)


def take_crystal(rcif_object: GlobalN, crystal_name: str = "") -> Crystal:
    """The crystal named ``crystal_name`` (data name without 'crystal_'),
    or the only crystal of the object."""
    l_crystal = [item for item in rcif_object.items if isinstance(item, Crystal)]
    if crystal_name:
        l_crystal = [item for item in l_crystal
                     if item.data_name.lower() == crystal_name.strip().lower()]
    if len(l_crystal) != 1:
        names = [item.data_name for item in rcif_object.items if isinstance(item, Crystal)]
        raise CrysPyError(
            f"Give the name of the crystal ({', '.join(names) or 'no crystals'}).")
    return l_crystal[0]


def take_moment_ranges(crystal: Crystal) -> dict:
    """{label: (smallest, largest |M| over the cells and the atoms of the
    orbit in μB, kind of modulation)} of the magnetic atoms."""
    from cryspy.A_functions_base.structure_factor_magnetic_k import (
        calc_moment_modulation_by_dictionary)
    with warnings.catch_warnings():
        # symmform and gauge warnings were given by the refinement already
        warnings.simplefilter("ignore")
        dict_crystal = crystal.get_dictionary()
    modulation = calc_moment_modulation_by_dictionary(dict_crystal)
    k_entry = numpy.asarray(modulation["k_vector"], dtype=float)
    d_range = {}
    for i_atom, label in enumerate(modulation["label"]):
        k_vector = k_entry if k_entry.ndim == 1 else k_entry[:, i_atom]
        moment_min, moment_max, kind = classify_moment_modulation(
            modulation["amplitude"][:, i_atom], k_vector,
            dict_crystal["unit_cell_parameters"])
        if label in d_range:
            old_min, old_max, old_kind = d_range[label]
            moment_min, moment_max = min(old_min, moment_min), max(old_max, moment_max)
            kind = old_kind if old_kind == kind else "mixed"
        d_range[label] = (moment_min, moment_max, kind)
    return d_range


def refine_representation(rcif_object: GlobalN, crystal_name: str, k_vector,
                          magnetic_atoms: str, number: int, mode, moment: float,
                          use_minus_k_operations: bool, magnetic_ion: str = "",
                          lande: float = 0.) -> dict:
    """Copy of the object with the representation Γ<number> (mode) of k,
    its Fourier components refined; chi_sq, n_point, the number of free
    parameters and the moment ranges."""
    # a copy through the rcif text: GlobalN cannot be deep-copied
    obj = str_to_globaln(rcif_object.to_cif())
    crystal = take_crystal(obj, crystal_name)
    crystal.add_propagation_vector(
        float(k_vector[0]), float(k_vector[1]), float(k_vector[2]), magnetic_atoms,
        str(number), moment, use_minus_k_operations, mode or 1, magnetic_ion, lande)
    for item in crystal.atom_site_moment_fourier.items:
        item.cos_refinement, item.sin_refinement = True, True
    with warnings.catch_warnings():
        # the Fourier components the data cannot define are fixed with a
        # warning per representation; the table gives the free ones
        warnings.simplefilter("ignore")
        obj_dict = obj.get_dictionary()
        _, parameter_name, _, _ = rhochi_rietveld_refinement_by_dictionary(obj_dict)
        obj.take_parameters_from_dictionary(obj_dict, l_parameter_name=parameter_name,
                                            l_sigma=None)
        d_out = rhochi_no_refinement(obj)
    n_fourier = sum(1 for name in parameter_name if "atom_fourier" in str(name))
    return {"object": obj, "chi_sq": float(d_out["chi_sq"]),
            "n_point": int(d_out["n_point"]), "n_fourier": n_fourier,
            "n_parameter": len(parameter_name),
            "moments": take_moment_ranges(crystal)}


def rhochi_rank_magnetic_representations(
        rcif_object: GlobalN, k_1: float, k_2: float, k_3: float,
        magnetic_atoms: str, moment: float = 1.0, crystal_name: str = "",
        keep_best: bool = False, use_minus_k_operations: bool = True,
        magnetic_ion: str = "", lande: float = 0.0) -> str:
    """Refine the magnetic structure of every representation of k and compare chi_sq.

    The paramagnetic crystal (space_group, fitted to the nuclear data) gets
    the propagation vector k with every one-dimensional irreducible
    representation of the group of k and with every mode (isotropy
    subgroup) of the multidimensional ones, as add_propagation_vector does,
    in a copy of the object. The Fourier components (cos and sin) of the
    magnetic atoms allowed by the representation are refined against all
    experiments of the object (powder patterns, single-crystal intensities
    with fractional indices of the satellites, flip ratios); the other
    parameters keep their refinement flags (fix them to compare the
    representations faster); single-crystal experiments are calculated for
    one magnetic domain unless they have jana_twin. The table lists chi_sq/n, the number of the
    refined Fourier components and the moment ranges, the best first. The
    object is not changed unless keep_best is true. Superpositions of
    representations (irrep '1,2', cycloids and helices of two modes) are
    not tried: add them with add_propagation_vector.

    Parameters
    ----------
    rcif_object : GlobalN
        The object with the crystal and the experiments.
    k_1 : float
        First component of the propagation vector (reciprocal lattice units).
    k_2 : float
        Second component of the propagation vector (reciprocal lattice units).
    k_3 : float
        Third component of the propagation vector (reciprocal lattice units).
    magnetic_atoms : str
        Labels of the magnetic atoms separated by commas, as Fe1, Fe2.
    moment : float
        Start moment of the magnetic atoms, μB (not zero: the intensities
        of the satellites do not change at zero moment to first order).
    crystal_name : str
        Name of the crystal; empty when the object has one crystal.
    keep_best : bool
        Replace the crystal of the object by the refined crystal of the
        best representation.
    use_minus_k_operations : bool
        Add the operations mapping k to -k when the representation allows it.
    magnetic_ion : str
        Magnetic ion of the form factor (as Tb3+) for the magnetic atoms
        without atom_site_scat; empty takes the type_symbol of the atom
        when it has a charge.
    lande : float
        Lande factor of these atoms; 0 takes g_J of Hund's rules for a 4f
        ion and 2 otherwise.

    Returns
    -------
    str
        Markdown table of the representations, the best first.
    """
    if abs(moment) < 1e-6:
        raise CrysPyError("Give a start moment (μB) that is not zero.")
    crystal = take_crystal(rcif_object, crystal_name)
    crystal_name = crystal.data_name
    k_vector, l_note = snap_k_vector([k_1, k_2, k_3])
    l_candidate = take_irrep_candidates(crystal.take_group_of_k(k_vector), k_vector)
    l_row, l_fail = [], []
    for number, mode in l_candidate:
        name = form_irrep_name([number], mode)
        logger.info(f"Refinement of {name}.")
        try:
            d_res = refine_representation(
                rcif_object, crystal_name, k_vector, magnetic_atoms, number, mode,
                moment, use_minus_k_operations, magnetic_ion, lande)
        except CrysPyError as error:
            l_fail.append(f"- {name}: {str(error).splitlines()[0]}")
            continue
        d_res["name"], d_res["number"], d_res["mode"] = name, number, mode
        l_row.append(d_res)
    if not l_row:
        raise CrysPyError("No representation could be refined:\n" + "\n".join(l_fail))
    l_row.sort(key=lambda row: row["chi_sq"])

    k_text = ", ".join(f"{float(value):g}" for value in k_vector)
    ls_out = [f"# Representations of k = ({k_text}) refined against the data", ""]
    if l_note:
        ls_out += ["k is taken as commensurate: " + ", ".join(l_note) + ".", ""]
    ls_out += ["| Rank | Representation | chi_sq/n | Fourier parameters | Moments (μB) |",
               "|---|---|---|---|---|"]
    for i_row, row in enumerate(l_row):
        text_moment = "; ".join(
            f"{label} {low:.2f}–{high:.2f} {kind}"
            for label, (low, high, kind) in row["moments"].items())
        ls_out.append(
            f"| {i_row + 1} | {row['name']} | {row['chi_sq'] / max(row['n_point'], 1):.3f} "
            f"| {row['n_fourier']} | {text_moment} |")
    ls_out += ["", f"n = {l_row[0]['n_point']} points. A representation with more Fourier "
               "parameters fits better by itself: prefer it only when chi_sq/n is clearly "
               "lower. Every model starts from one start moment; a model far from the "
               "minimum can be tried again with add_propagation_vector and a refinement."]
    if l_fail:
        ls_out += ["", "Not refined:"] + l_fail
    best = l_row[0]
    if keep_best:
        i_item = rcif_object.items.index(crystal)
        rcif_object.items[i_item] = take_crystal(best["object"], crystal_name)
        ls_out += ["", f"The crystal now has k with {best['name']}, refined."]
    else:
        mode_text = f", mode {best['mode']}" if best["mode"] is not None else ""
        ls_out += ["", f"The object is not changed: add_propagation_vector with irrep "
                   f"{best['number']}{mode_text} gives the best model."]
    return "\n".join(ls_out)
