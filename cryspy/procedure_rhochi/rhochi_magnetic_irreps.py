"""Choice of the magnetic representation of a propagation vector by the data.

For every irreducible representation of the group of k (and every mode of
the multidimensional ones) a copy of the object gets the propagation vector
(Crystal.add_propagation_vector), its Fourier components are refined against
the experiments of the object and the chi_sq values are compared, as the
representation analysis of SARAh, BasIreps or ISODISTORT followed by a
refinement of every model (FullProf, Jana2020).
"""
import copy
import logging
from fractions import Fraction
import warnings

import numpy

from cryspy.A_functions_base.magnetic_k_symmetry import snap_k_vector
from cryspy.B_parent_classes.cl_4_global import GlobalN
from cryspy.D_functions_item_loop.function_1_magnetic_k import (
    add_magnetic_scattering,
    take_k_group,
    classify_moment_modulation,
    form_irrep_name,
    split_labels,
    take_model_candidates,
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


def refine_representation(rcif_object: GlobalN, crystal_name: str, k_vector,  # noqa: PLR0913
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
    """(model number, irrep text, mode) of every candidate, numbered as in
    Crystal.report_magnetic_models (take_model_candidates)."""
    return [(i_model + 1, irrep, mode) for i_model, (irrep, mode) in enumerate(
        take_model_candidates(k_group, k_vector, try_superpositions))]


def rhochi_rank_magnetic_representations(  # noqa: C901, PLR0912, PLR0913, PLR0915
        rcif_object: GlobalN, k_1: float, k_2: float, k_3: float,
        magnetic_atoms: str, moment: float = 1.0, crystal_name: str = "",
        keep_rank: int = 0, use_minus_k_operations: bool = True,
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
    object is not changed unless keep_rank is the number of a row of the
    table: then its refined model replaces the crystal.

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
    keep_rank : int
        Rank (row number of the table) of the representation whose refined
        crystal replaces the crystal of the object: 1 is the best chi_sq/n,
        2 the second and so on; the report names the rank of the
        recommended one (the fewest Fourier parameters among those with
        chi_sq/n within 5 % of the best). 0 (default) or a number that is
        not in the table leaves the object unchanged.
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
    keep_rank = int(keep_rank)
    crystal = take_crystal(rcif_object, crystal_name)
    crystal_name = crystal.data_name
    k_vector, l_note = snap_k_vector([k_1, k_2, k_3])
    # the form factor every model gets: the Lande factor taken for lande 0 is shown
    l_scat_added = add_magnetic_scattering(
        copy.deepcopy(crystal), split_labels(magnetic_atoms), magnetic_ion, lande)
    l_candidate = take_rank_candidates(
        crystal.take_group_of_k(k_vector), k_vector, try_superpositions)
    l_row, l_fail = [], []
    for model, number, mode in l_candidate:
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
        d_res["model"] = model
        l_row.append(d_res)
    if not l_row:
        raise CrysPyError("No representation could be refined:\n" + "\n".join(l_fail))
    l_row.sort(key=lambda row: row["chi_sq"])

    k_text = ", ".join(f"{float(value):g}" for value in k_vector)
    ls_out = [f"# Representations of k = ({k_text}) refined against the data", ""]
    if l_note:
        ls_out += ["k is taken as commensurate: " + ", ".join(l_note) + ".", ""]
    if l_scat_added:
        ls_out += ["Magnetic form factor of every model (atom_site_scat): "
                   + "; ".join(l_scat_added) + ".", ""]
    ls_out += ["| Rank | Model | Representation | chi_sq/n | R_F2 satellites of k "
               "| Fourier parameters | Domains | Moments (μB) |",
               "|---|---|---|---|---|---|---|---|"]
    for i_row, row in enumerate(l_row):
        text_moment = "; ".join(
            f"{label} {low:.2f}–{high:.2f} {kind}"
            for label, (low, high, kind) in row["moments"].items())
        text_r = "–" if row["r_satellite"] is None else f"{100. * row['r_satellite']:.1f} %"
        ls_out.append(
            f"| {i_row + 1} | {row['model']} | {row['name']} | {row['chi_sq'] / max(row['n_point'], 1):.3f} "
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
    if any(row["mode"] is not None for row in l_row):
        ls_out += ["", "Γn mode m: representation Γn has dimension 2 or more, so its "
                   "Fourier components are a combination of several basis functions; mode m "
                   "is one isotropy subgroup (direction of the order parameter in the space "
                   "of Γn, numbered as in report_irreps_of_k, mode 1 keeps the most "
                   "symmetry). The other directions of the same mode are its magnetic "
                   "domains."]
    if l_fail:
        ls_out += ["", "Not refined:"] + l_fail
    best = take_recommended_row(l_row)
    rank_best = l_row.index(best) + 1
    if best is not l_row[0]:
        ls_out += ["", f"Recommended: {best['name']} (rank {rank_best}, {best['n_fourier']} "
                   f"Fourier parameters instead of {l_row[0]['n_fourier']}): its chi_sq/n is "
                   f"within {100. * CHI_SQ_SIMILAR:.0f} % of the best, so the simpler model "
                   f"is preferred (keep_rank {rank_best} keeps it)."]
    else:
        ls_out += ["", f"Recommended: {best['name']} (rank 1), the best chi_sq/n."]
    if 1 <= keep_rank <= len(l_row):
        kept = l_row[keep_rank - 1]
        i_item = rcif_object.items.index(crystal)
        rcif_object.items[i_item] = take_crystal(kept["object"], crystal_name)
        text_domain = ""
        if kept["n_domain"] > 1 and with_domains:
            # the experiments keep the magnetic domains the model was refined with
            from cryspy.E_data_classes.cl_2_diffrn import Diffrn
            d_kept = {item.data_name: item for item in kept["object"].items
                      if isinstance(item, Diffrn)}
            for i_item, item in enumerate(rcif_object.items):
                if isinstance(item, Diffrn) and item.data_name in d_kept:
                    rcif_object.items[i_item] = d_kept[item.data_name]
            text_domain = f" The single-crystal experiments have its {kept['n_domain']} " \
                "magnetic domains (equal fractions)."
        ls_out += ["", f"The crystal now has k with {kept['name']} (rank {keep_rank}), "
                   f"refined.{text_domain} Another row is kept by a new ranking with its "
                   "keep_rank, or by add_propagation_vector on the paramagnetic crystal."]
    else:
        text_domain = ""
        if best["n_domain"] > 1 and with_domains:
            text_domain = " Then define_magnetic_domains gives the single-crystal experiments " \
                f"its {best['n_domain']} magnetic domains, as in the ranking."
        mode_text = f", mode {best['mode']}" if best["mode"] is not None else ""
        text_rank = "" if keep_rank == 0 else \
            f"keep_rank {keep_rank} is not a row of the table (1 to {len(l_row)}). "
        ls_out += ["", f"{text_rank}The object is not changed: add_propagation_vector with "
                   f"irrep {best['number']}{mode_text} gives the recommended model (or a "
                   f"ranking with keep_rank {rank_best}).{text_domain}"]
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


def take_propagation_vectors(rcif_object: GlobalN) -> dict:
    """Propagation vectors of the measured satellites (see
    find_propagation_vectors) as data.

    Returns
    -------
    dict
        "vectors": [{"k" (tuple of 3 floats, snapped to fractions),
        "n_satellite", "experiments" ({name: (count, summed intensity)}),
        "point" (label of the special point of the Brillouin zone of the
        only crystal, see name_special_k_point; '' if none)}],
        the vector with more satellites first; "n_integer" (reflections
        with integer indices); "experiments" (names of the single-crystal
        experiments with reflections).
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
    l_crystal = [item for item in rcif_object.items if isinstance(item, Crystal)]
    crystal = l_crystal[0] if len(l_crystal) == 1 else None
    l_key = sorted(d_k, key=lambda key: -sum(count for count, _ in d_k[key].values()))
    l_vector = []
    for key in l_key:
        k_vector = snap_k_vector(list(key))[0]
        l_vector.append({"k": tuple(float(value) for value in k_vector),
                         "n_satellite": sum(count for count, _ in d_k[key].values()),
                         "experiments": d_k[key],
                         "point": "" if crystal is None else name_special_k_point(crystal, k_vector)})
    return {"vectors": l_vector, "n_integer": n_integer, "experiments": l_experiment}


def form_k_fraction_text(k_vector) -> str:
    """'(1/2, 1/2, 0)' of a propagation vector (fractions up to 1/12)."""
    l_text = []
    for value in k_vector:
        fraction = Fraction(float(value)).limit_denominator(12)
        text = str(fraction) if abs(float(fraction) - float(value)) < 1e-4 else f"{float(value):g}"
        l_text.append(text)
    return "(" + ", ".join(l_text) + ")"


def name_special_k_point(crystal: Crystal, k_vector) -> str:
    """Label of the special point of the Brillouin zone (knowledge base,
    special_k_points.md) equal to k up to the operations of the
    paramagnetic group and a reciprocal lattice vector of the centred cell
    ('S' for (1/2, 1/2, 0) of C m c m); '' when it is not one of them or the
    lattice is not in the table."""
    from cryspy.A_functions_base.knowledge_tables import take_special_k_points
    from cryspy.D_functions_item_loop.function_1_magnetic_k import take_parent_space_group
    try:
        space_group = take_parent_space_group(crystal)
        lattice = f"{space_group.crystal_system} {str(space_group.name_hm_alt).strip()[0]}"
        k_group = take_k_group(space_group, numpy.asarray(k_vector, dtype=float))
    except (AttributeError, IndexError, CrysPyError):
        return ""
    l_point = take_special_k_points().get(lattice, [])
    rotation = numpy.asarray(k_group["parent_rotation"], dtype=float).reshape(3, 3, -1)
    centring = numpy.asarray(k_group["centring"], dtype=float)
    k_vector = numpy.asarray(k_vector, dtype=float)
    for label, k_point in l_point:
        for i_s in range(rotation.shape[2]):
            for sign in (1., -1.):
                # a reciprocal vector transforms with the transposed rotation
                difference = sign * rotation[:, :, i_s].transpose() @ k_vector - numpy.array(k_point)
                if not numpy.allclose(difference, numpy.round(difference), atol=1e-4):
                    continue
                phase = numpy.round(difference) @ centring
                if numpy.allclose(phase, numpy.round(phase), atol=1e-4):
                    return label
    return ""


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
    d_found = take_propagation_vectors(rcif_object)
    n_integer = d_found["n_integer"]
    ls_out = ["# Propagation vectors of the measured satellites", ""]
    if not d_found["vectors"]:
        ls_out.append(f"No satellites: all {n_integer} reflections have integer indices "
                      "(k = 0 or a paramagnetic state).")
        return "\n".join(ls_out)
    ls_out += ["| k | Point | Satellites | Experiments (satellites, summed intensity) |",
               "|---|---|---|---|"]
    for d_vector in d_found["vectors"]:
        k_text = ", ".join(f"{float(value):g}" for value in d_vector["k"])
        text_exp = "; ".join(f"{name} ({count}, {total:.4g})"
                             for name, (count, total) in d_vector["experiments"].items())
        ls_out.append(f"| ({k_text}) | {d_vector['point'] or '–'} | {d_vector['n_satellite']} "
                      f"| {text_exp} |")
    ls_out += ["", f"{n_integer} reflections with integer indices (nuclear and k = 0). Point: "
               "the special point of the Brillouin zone of the paramagnetic group "
               "(search_knowledge 'special points')."]
    if len(d_found["vectors"]) > 1:
        ls_out.append(
            f"{len(d_found['vectors'])} propagation vectors: add them one by one (the first k with "
            "rhochi_rank_magnetic_representations and add_propagation_vector, then the "
            "second one on the crystal that has the first: a multi-k structure). Vectors of "
            "one star (related by the paramagnetic group) are one k with magnetic domains.")
    return "\n".join(ls_out)
