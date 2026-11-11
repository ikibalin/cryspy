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


def add_domains_to_experiments(rcif_object: GlobalN, crystal: Crystal) -> int:
    """Magnetic domains of the crystal (equal fractions, not refined) in
    every single-crystal experiment of its phase; the number of domains."""
    from cryspy.E_data_classes.cl_2_diffrn import Diffrn
    from cryspy.D_functions_item_loop.function_1_magnetic_k import (
        calc_magnetic_domain_matrices)
    n_domain = int(calc_magnetic_domain_matrices(crystal).shape[1])
    if n_domain < 2:
        return n_domain
    for item in rcif_object.items:
        if isinstance(item, Diffrn) and item.is_attribute("phase") and \
                str(item.phase.label).lower() == crystal.data_name.lower():
            item.define_magnetic_domains(crystal, flag_refine=False)
    return n_domain


def is_satellite_of_k(index_hkl, k_vector, tolerance: float = 1e-4):
    """Reflections (3, n) that are satellites H = h ± k of k (h integer)."""
    index_hkl = numpy.asarray(index_hkl, dtype=float)
    k_vector = numpy.asarray(k_vector, dtype=float)[:, numpy.newaxis]
    flag_fractional = numpy.any(numpy.abs(index_hkl - numpy.round(index_hkl)) > tolerance, axis=0)
    flag = numpy.zeros(index_hkl.shape[1], dtype=bool)
    for sign in (1., -1.):
        shifted = index_hkl - sign * k_vector
        flag |= numpy.all(numpy.abs(shifted - numpy.round(shifted)) < tolerance, axis=0)
    return flag & flag_fractional


def calc_satellite_agreement(rcif_object: GlobalN, k_vector):
    """R_F2 = sum |I_obs - I_calc| / sum I_obs and the number of the
    satellites of k measured in the single-crystal intensity experiments
    (excluded reflections are left out); (None, 0) without them."""
    from cryspy.E_data_classes.cl_2_diffrn import Diffrn
    sum_diff, sum_obs, n_satellite = 0., 0., 0
    for item in rcif_object.items:
        if not (isinstance(item, Diffrn) and item.is_attribute("diffrn_refln")):
            continue
        refln = item.diffrn_refln
        if not (refln.is_attribute("intensity") and refln.is_attribute("intensity_calc")):
            continue
        index_hkl = numpy.array([refln.index_h, refln.index_k, refln.index_l], dtype=float)
        flag = is_satellite_of_k(index_hkl, k_vector)
        if refln.is_attribute("excluded"):
            flag &= numpy.logical_not(numpy.array(refln.excluded, dtype=bool))
        intensity = numpy.array(refln.intensity, dtype=float)[flag]
        intensity_calc = numpy.array(refln.intensity_calc, dtype=float)[flag]
        sum_diff += float(numpy.abs(intensity - intensity_calc).sum())
        sum_obs += float(numpy.abs(intensity).sum())
        n_satellite += int(flag.sum())
    if n_satellite == 0 or sum_obs == 0.:
        return None, 0
    return sum_diff / sum_obs, n_satellite


def refine_representation(rcif_object: GlobalN, crystal_name: str, k_vector,
                          magnetic_atoms: str, number, mode, moment: float,
                          use_minus_k_operations: bool, magnetic_ion: str = "",
                          lande: float = 0., with_domains: bool = True) -> dict:
    """Copy of the object with the representation Γ<number> (mode) of k
    (number may be '1,2' for a superposition of one-dimensional ones), its
    Fourier components refined; chi_sq, n_point, the number of free
    parameters, the moment ranges, R_F2 of the satellites of k and the
    number of magnetic domains given to the single-crystal experiments."""
    # a copy through the rcif text: GlobalN cannot be deep-copied
    obj = str_to_globaln(rcif_object.to_cif())
    crystal = take_crystal(obj, crystal_name)
    crystal.add_propagation_vector(
        float(k_vector[0]), float(k_vector[1]), float(k_vector[2]), magnetic_atoms,
        str(number), moment, use_minus_k_operations, mode or 1, magnetic_ion, lande)
    for item in crystal.atom_site_moment_fourier.items:
        item.cos_refinement, item.sin_refinement = True, True
    n_domain = add_domains_to_experiments(obj, crystal) if with_domains else 1
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
    r_satellite, n_satellite = calc_satellite_agreement(obj, k_vector)
    return {"object": obj, "chi_sq": float(d_out["chi_sq"]),
            "n_point": int(d_out["n_point"]), "n_fourier": n_fourier,
            "n_parameter": len(parameter_name),
            "moments": take_moment_ranges(crystal),
            "r_satellite": r_satellite, "n_satellite": n_satellite, "n_domain": n_domain}


# models fit equally well when chi_sq/n differs by less than 5 % of the best
# one or by less than 0.05 (data fitted within the sigmas)
CHI_SQ_SIMILAR = 0.05


def take_recommended_row(l_row: list) -> dict:
    """The row with the fewest Fourier parameters among those whose chi_sq/n
    is within CHI_SQ_SIMILAR (relative, or absolute near a perfect fit) of
    the best one (rows sorted by chi_sq): more parameters always fit a
    little better by themselves."""
    def chi_sq_n(row):
        return row["chi_sq"] / max(row["n_point"], 1)
    limit = chi_sq_n(l_row[0]) * (1. + CHI_SQ_SIMILAR) + CHI_SQ_SIMILAR
    l_similar = [row for row in l_row if chi_sq_n(row) <= limit]
    return min(l_similar, key=lambda row: (row["n_fourier"], row["chi_sq"]))


def take_rank_candidates(k_group: dict, k_vector, try_superpositions: bool) -> list:
    """(irrep text, mode) of every candidate: the one-dimensional
    representations, the modes of the multidimensional ones and, with
    try_superpositions, the sums of two one-dimensional ones ('1,2')."""
    l_candidate = [(str(number), mode)
                   for number, mode in take_irrep_candidates(k_group, k_vector)]
    if try_superpositions:
        l_one = [number for number, mode in l_candidate if mode is None]
        l_candidate += [(f"{first},{second}", None) for i_first, first in enumerate(l_one)
                        for second in l_one[i_first + 1:]]
    return l_candidate


def rhochi_rank_magnetic_representations(
        rcif_object: GlobalN, k_1: float, k_2: float, k_3: float,
        magnetic_atoms: str, moment: float = 1.0, crystal_name: str = "",
        keep_best: bool = False, use_minus_k_operations: bool = True,
        magnetic_ion: str = "", lande: float = 0.0, with_domains: bool = True,
        try_superpositions: bool = True) -> str:
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
    representations faster). Single-crystal experiments of the crystal get
    its magnetic domains with equal fractions (with_domains), as
    Diffrn.define_magnetic_domains does: the arms of the star of k and the
    orientation domains change the satellites seen by the data. With
    try_superpositions the sums of two one-dimensional representations
    (irrep '1,2': cycloids, helices, canted structures) are refined too.
    The table lists chi_sq/n of all experiments, R_F2 of the satellites of
    this k only (independent of the satellites of other propagation
    vectors), the number of the refined Fourier components and the moment
    ranges, the best first. A crystal with a propagation vector already
    gets k as a second one (multi-k, see add_propagation_vector). The
    object is not changed unless keep_best is true.

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
        recommended representation (the fewest Fourier parameters among
        those with chi_sq/n within 5 % of the best).
    use_minus_k_operations : bool
        Add the operations mapping k to -k when the representation allows it.
    magnetic_ion : str
        Magnetic ion of the form factor (as Tb3+) for the magnetic atoms
        without atom_site_scat; empty takes the type_symbol of the atom
        when it has a charge.
    lande : float
        Lande factor of these atoms; 0 takes g_J of Hund's rules for a 4f
        ion and 2 otherwise.
    with_domains : bool
        Magnetic domains of every model in the single-crystal experiments
        (equal fractions).
    try_superpositions : bool
        Sums of two one-dimensional representations refined as well.

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
    l_candidate = take_rank_candidates(
        crystal.take_group_of_k(k_vector), k_vector, try_superpositions)
    l_row, l_fail = [], []
    for number, mode in l_candidate:
        name = form_irrep_name(number.split(","), mode)
        logger.info(f"Refinement of {name}.")
        try:
            d_res = refine_representation(
                rcif_object, crystal_name, k_vector, magnetic_atoms, number, mode,
                moment, use_minus_k_operations, magnetic_ion, lande, with_domains)
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
    ls_out += ["| Rank | Representation | chi_sq/n | R_F2 satellites of k | Fourier parameters "
               "| Domains | Moments (μB) |",
               "|---|---|---|---|---|---|---|"]
    for i_row, row in enumerate(l_row):
        text_moment = "; ".join(
            f"{label} {low:.2f}–{high:.2f} {kind}"
            for label, (low, high, kind) in row["moments"].items())
        text_r = "–" if row["r_satellite"] is None else f"{100. * row['r_satellite']:.1f} %"
        ls_out.append(
            f"| {i_row + 1} | {row['name']} | {row['chi_sq'] / max(row['n_point'], 1):.3f} "
            f"| {text_r} | {row['n_fourier']} | {row['n_domain']} | {text_moment} |")
    n_satellite = l_row[0]["n_satellite"]
    ls_out += ["", f"n = {l_row[0]['n_point']} points, {n_satellite} satellites of k in the "
               "single-crystal experiments. A representation with more Fourier "
               "parameters fits better by itself: prefer it only when chi_sq/n is clearly "
               "lower. The moments are the range of one wave; with several propagation "
               "vectors report_magnetic_moments gives the total moment, which should stay "
               "below g·J of the free ion. Every model starts from one start moment; a model "
               "far from the minimum can be tried again with add_propagation_vector and a "
               "refinement."]
    if l_fail:
        ls_out += ["", "Not refined:"] + l_fail
    best = take_recommended_row(l_row)
    if best is not l_row[0]:
        ls_out += ["", f"Recommended: {best['name']} ({best['n_fourier']} Fourier parameters "
                   f"instead of {l_row[0]['n_fourier']}): its chi_sq/n is within "
                   f"{100. * CHI_SQ_SIMILAR:.0f} % of the best, so the simpler model is "
                   "preferred."]
    else:
        ls_out += ["", f"Recommended: {best['name']}, the best chi_sq/n."]
    mode_text = f", mode {best['mode']}" if best["mode"] is not None else ""
    if keep_best:
        i_item = rcif_object.items.index(crystal)
        rcif_object.items[i_item] = take_crystal(best["object"], crystal_name)
        text_domain = ""
        if best["n_domain"] > 1 and with_domains:
            # the experiments keep the magnetic domains the model was refined with
            from cryspy.E_data_classes.cl_2_diffrn import Diffrn
            d_best = {item.data_name: item for item in best["object"].items
                      if isinstance(item, Diffrn)}
            for i_item, item in enumerate(rcif_object.items):
                if isinstance(item, Diffrn) and item.data_name in d_best:
                    rcif_object.items[i_item] = d_best[item.data_name]
            text_domain = f" The single-crystal experiments have its {best['n_domain']} " \
                "magnetic domains (equal fractions)."
        ls_out += ["", f"The crystal now has k with {best['name']}, refined.{text_domain}"]
    else:
        text_domain = ""
        if best["n_domain"] > 1 and with_domains:
            text_domain = " Then define_magnetic_domains gives the single-crystal experiments " \
                f"its {best['n_domain']} magnetic domains, as in the ranking."
        ls_out += ["", f"The object is not changed: add_propagation_vector with irrep "
                   f"{best['number']}{mode_text} gives the recommended model.{text_domain}"]
    return "\n".join(ls_out)


def take_canonical_k(fraction, tolerance: float = 1e-4):
    """The propagation vector of a satellite as k or -k reduced to [0, 1):
    the one of the two that is lexicographically smaller, so that h + k
    and h - k give the same k."""
    l_k = []
    for sign in (1., -1.):
        k_vector = numpy.mod(sign * numpy.asarray(fraction, dtype=float), 1.)
        k_vector[numpy.abs(k_vector - 1.) < tolerance] = 0.
        k_vector[numpy.abs(k_vector) < tolerance] = 0.
        l_k.append(tuple(float(value) for value in numpy.round(k_vector, 4)))
    return min(l_k)


def find_propagation_vectors(rcif_object: GlobalN) -> str:
    """Find the propagation vectors of the magnetic satellites in the data.

    The measured reflections of the single-crystal experiments with
    fractional indices H = h ± k (the satellites) give k as the fractional
    part of H; k and -k, and k shifted by a reciprocal lattice vector,
    are the same vector. Commensurate components are written as fractions
    with a denominator up to 12. The vectors are listed with the number of
    their satellites and the summed measured intensity in every
    experiment, so the magnetic structure can be built for each of them
    (add_propagation_vector, rhochi_rank_magnetic_representations). Nothing
    in the object is changed.

    Parameters
    ----------
    rcif_object : GlobalN
        The object with the single-crystal experiments.

    Returns
    -------
    str
        Markdown table of the propagation vectors.
    """
    from cryspy.E_data_classes.cl_2_diffrn import Diffrn
    d_k, l_experiment, n_integer = {}, [], 0
    for item in rcif_object.items:
        if not (isinstance(item, Diffrn) and item.is_attribute("diffrn_refln")):
            continue
        refln = item.diffrn_refln
        index_hkl = numpy.array([refln.index_h, refln.index_k, refln.index_l], dtype=float)
        if refln.is_attribute("intensity"):
            intensity = numpy.array(refln.intensity, dtype=float)
        else:
            intensity = numpy.full(index_hkl.shape[1], numpy.nan)
        flag_fractional = numpy.any(numpy.abs(index_hkl - numpy.round(index_hkl)) > 1e-4, axis=0)
        n_integer += int((~flag_fractional).sum())
        l_experiment.append(item.data_name)
        for hkl, value in zip(index_hkl[:, flag_fractional].transpose(), intensity[flag_fractional]):
            key = take_canonical_k(hkl - numpy.round(hkl))
            d_exp = d_k.setdefault(key, {})
            count, total = d_exp.get(item.data_name, (0, 0.))
            d_exp[item.data_name] = (count + 1, total + (0. if numpy.isnan(value) else value))
    if not l_experiment:
        raise CrysPyError("The object has no single-crystal experiments with reflections.")
    ls_out = ["# Propagation vectors of the measured satellites", ""]
    if not d_k:
        ls_out.append(f"No satellites: all {n_integer} reflections have integer indices "
                      "(k = 0 or a paramagnetic state).")
        return "\n".join(ls_out)
    ls_out += ["| k | Satellites | Experiments (satellites, summed intensity) |", "|---|---|---|"]
    l_key = sorted(d_k, key=lambda key: -sum(count for count, _ in d_k[key].values()))
    for key in l_key:
        k_vector, l_note = snap_k_vector(list(key))
        k_text = ", ".join(f"{float(value):g}" for value in k_vector)
        text_exp = "; ".join(f"{name} ({count}, {total:.4g})"
                             for name, (count, total) in d_k[key].items())
        n_total = sum(count for count, _ in d_k[key].values())
        ls_out.append(f"| ({k_text}) | {n_total} | {text_exp} |")
    ls_out += ["", f"{n_integer} reflections with integer indices (nuclear and k = 0)."]
    if len(l_key) > 1:
        ls_out.append(
            f"{len(l_key)} propagation vectors: add them one by one (the first k with "
            "rhochi_rank_magnetic_representations and add_propagation_vector, then the "
            "second one on the crystal that has the first: a multi-k structure). Vectors of "
            "one star (related by the paramagnetic group) are one k with magnetic domains.")
    return "\n".join(ls_out)
