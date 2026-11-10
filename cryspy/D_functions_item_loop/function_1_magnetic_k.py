"""From a paramagnetic crystal to a magnetic structure with a propagation vector.

The functions work on the items and loops of a crystal data block and are
called by the methods ``Crystal.report_irreps_of_k``,
``Crystal.add_propagation_vector`` and ``Crystal.report_magnetic_domains``.
The numerical part is in ``A_functions_base.magnetic_k_symmetry``.
"""
import copy
import re
from fractions import Fraction
from math import gcd

import numpy

from cryspy.exceptions import CrysPyError
from cryspy.A_functions_base.magnetic_k_symmetry import (
    calc_axial_moment_matrix, calc_domain_matrices, calc_fourier_site_projector,
    calc_k_group_flags, calc_minus_k_flags, calc_one_dimensional_irreps,
    form_fourier_symmform, form_operation_xyz, form_superspace_group,
    split_atom_orbit)
from cryspy.A_functions_base.structure_factor_magnetic_k import (
    calc_superspace_operation_arrays, is_k_equivalent_to_minus_k,
    parse_operation_xyz, parse_ssg_operation)
from cryspy.A_functions_base.unit_cell import (
    calc_reciprocal_by_unit_cell_parameters)

na = numpy.newaxis


def take_space_group_operations(space_group):
    """Rotations (9, n), translations (3, n) without the centring and the
    centring translations (3, n_c) (the first one is zero) of a SpaceGroup;
    the identity is the first operation."""
    l_rotation, l_translation = [], []
    for operation_xyz in space_group.full_space_group_symop.operation_xyz:
        rotation, translation = parse_operation_xyz(operation_xyz)
        if any(numpy.allclose(rotation, other) for other in l_rotation):
            continue  # the same operation with a centring translation
        l_rotation.append(rotation)
        l_translation.append(translation)
    centring = numpy.array(
        [[float(value) for value in shift] for shift in space_group.shift],
        dtype=float).transpose()
    return (numpy.stack(l_rotation, axis=1), numpy.stack(l_translation, axis=1),
            centring)


def take_k_group(space_group, k_vector):
    """Operations of the group of k of a SpaceGroup.

    Returns
    -------
    dict
        ``parent_rotation``, ``parent_translation`` (all operations),
        ``rotation``, ``translation`` (group of k), ``centring``,
        ``irreps`` (list of complex characters), ``moment`` (9, n_s) axial
        matrices.
    """
    k_vector = numpy.asarray(k_vector, dtype=float)
    parent_rotation, parent_translation, centring = take_space_group_operations(
        space_group)
    flags = calc_k_group_flags(parent_rotation, k_vector)
    rotation, translation = parent_rotation[:, flags], parent_translation[:, flags]
    irreps = calc_one_dimensional_irreps(rotation, translation, centring, k_vector)
    return {
        "parent_rotation": parent_rotation, "parent_translation": parent_translation,
        "rotation": rotation, "translation": translation, "centring": centring,
        "irreps": irreps, "moment": calc_axial_moment_matrix(rotation)}


def form_phase_by_character(chi):
    """Moment phases (units of 2π, in [0, 1)) of the operations from the
    characters: ``chi = exp(-2 pi i phase)``."""
    phase = numpy.mod(-numpy.angle(chi) / (2. * numpy.pi), 1.)
    phase[numpy.isclose(phase, 1.)] = 0.
    return numpy.round(phase, 8)


def form_vector_text(vector) -> str:
    """'(0, 0, 1)' or '(1, 0.5i, 0)' for a complex vector."""
    l_text = []
    for value in vector:
        real, imag = float(numpy.round(value.real, 4)), float(numpy.round(value.imag, 4))
        if abs(imag) < 1e-4:
            l_text.append(f"{real:g}")
        elif abs(real) < 1e-4:
            l_text.append(f"{imag:g}i")
        else:
            l_text.append(f"{real:g}{imag:+g}i")
    return "(" + ", ".join(l_text) + ")"


def take_irrep_numbers(irrep, n_irrep: int) -> list:
    """Numbers of the representations from 1, '2' or '2,3' (Γ2 ⊕ Γ3).

    Raises CrysPyError when a number is out of range or the group of k has
    no one-dimensional representations.
    """
    if n_irrep == 0:
        raise CrysPyError(
            "The group of k has no one-dimensional representations (all "
            "small representations are multidimensional); this is not "
            "supported yet: give the magnetic superspace group (mCIF, e.g. "
            "from ISODISTORT or MAGNDATA) instead.")
    text = str(irrep).replace("Γ", "").replace("+", ",").replace("⊕", ",")
    try:
        numbers = sorted({int(part) for part in text.split(",") if part.strip()})
    except ValueError:
        raise CrysPyError(f"irrep '{irrep}' is not a number or a list such as 2,3.") from None
    if not numbers or not all(1 <= number <= n_irrep for number in numbers):
        raise CrysPyError(f"irrep must be from 1 to {n_irrep} (or a list such as 1,2).")
    return numbers


def form_irrep_name(numbers: list) -> str:
    """'Γ2' or 'Γ2⊕Γ3'."""
    return "⊕".join(f"Γ{number}" for number in numbers)


def form_superspace_by_irrep(k_group: dict, k_vector, irrep,
                             use_minus_k_operations: bool = True) -> dict:
    """Magnetic superspace group of a one-dimensional representation.

    For several representations (``irrep`` = '2,3') the group keeps the
    operations of the group of k on which their characters are equal: the
    superposition of the modes (cycloid, helix) is invariant there, up to
    a shift along x4.

    Parameters
    ----------
    k_group : dict
        Result of :func:`take_k_group`.
    k_vector : numpy.ndarray
        (3,) propagation vector.
    irrep : int or str
        Number of the representation (1 for Γ1) or numbers separated by
        commas ('2,3' for Γ2 ⊕ Γ3).
    use_minus_k_operations : bool
        Add an operation of the space group mapping k to -k when the
        representation allows it (the highest superspace group).

    Returns
    -------
    dict
        ``operations``, ``centrings`` (algebraic strings), ``arrays`` (see
        ``calc_superspace_operation_arrays``) and ``minus_k`` (True when the
        operations mapping k to -k are in the group).
    """
    k_vector = numpy.asarray(k_vector, dtype=float)
    numbers = take_irrep_numbers(irrep, len(k_group["irreps"]))
    l_chi = [k_group["irreps"][number - 1] for number in numbers]
    flags = numpy.all([numpy.isclose(chi, l_chi[0], atol=1e-6) for chi in l_chi], axis=0)
    phase = form_phase_by_character(l_chi[0][flags])
    arguments = (k_group["rotation"][:, flags], k_group["translation"][:, flags],
                 phase, k_group["centring"], k_vector)
    ssg, minus_k = None, False
    flags = calc_minus_k_flags(k_group["parent_rotation"], k_vector)
    if (use_minus_k_operations and numpy.any(flags)
            and not is_k_equivalent_to_minus_k(k_vector)):
        i_minus = int(numpy.flatnonzero(flags)[0])
        ssg = form_superspace_group(
            *arguments, minus_rotation=k_group["parent_rotation"][:, i_minus],
            minus_translation=k_group["parent_translation"][:, i_minus])
        minus_k = ssg is not None
    if ssg is None:
        ssg = form_superspace_group(*arguments)
    ssg["arrays"] = calc_superspace_operation_arrays(
        ssg["operations"], ssg["centrings"], k_vector)
    ssg["minus_k"] = minus_k
    return ssg


def calc_site_symmform(position, k_vector, arrays: dict) -> list:
    """Symmforms (6,) of cos x, y, z and sin x, y, z of an atom at a position
    in the superspace group given by its operation arrays."""
    projector = calc_fourier_site_projector(
        position, k_vector, arrays["rotation"], arrays["translation"],
        arrays["moment"], arrays["phase"], arrays["conjugate"])
    return form_fourier_symmform(projector), projector


def format_symmform(symmform: list) -> str:
    """'cos (0, myc1, 0), sin (mxs1, 0, mzs1)' or 'zero'."""
    if all(text == "0" for text in symmform):
        return "zero"
    return f"cos ({', '.join(symmform[:3])}), sin ({', '.join(symmform[3:])})"


def report_irreps(k_group: dict, k_vector, atoms: dict,
                  use_minus_k_operations: bool = True) -> str:
    """Markdown report: the operations of the group of k, their
    one-dimensional representations, the superspace group of every
    representation and the Fourier components allowed for every magnetic
    atom (symmforms as in MAGNDATA: ``M = cos cos(2 pi x4) + sin sin(2 pi
    x4)``, equal symbols are equal parameters).

    ``atoms`` is {label: fractional position (3,)}.
    """
    rotation, translation = k_group["rotation"], k_group["translation"]
    k_text = ", ".join(f"{float(value):g}" for value in k_vector)
    ls_out = [f"# Group of k = ({k_text})", "",
              f"{rotation.shape[1]} of {k_group['parent_rotation'].shape[1]} "
              "operations of the space group keep k (modulo the reciprocal lattice).",
              ""]
    header = "| Operation | " + " | ".join(
        f"Γ{i + 1}" for i in range(len(k_group["irreps"]))) + " |"
    ls_out += ["Characters of the one-dimensional representations:", "",
               header, "|---" * (len(k_group["irreps"]) + 1) + "|"]
    for s in range(rotation.shape[1]):
        cells = [form_vector_text([chi[s]])[1:-1] for chi in k_group["irreps"]]
        ls_out.append(f"| {form_operation_xyz(rotation[:, s], translation[:, s])} | "
                      + " | ".join(cells) + " |")
    if len(k_group["irreps"]) == 0:
        ls_out.append("No one-dimensional representation: write the "
                      "superspace operations (space_group_symop_magn_ssg_operation) "
                      "by hand.")
    for i_irrep in range(len(k_group["irreps"])):
        ssg = form_superspace_by_irrep(
            k_group, k_vector, i_irrep + 1, use_minus_k_operations)
        ls_out += ["", f"## Γ{i_irrep + 1}", "",
                   f"Superspace group: {len(ssg['operations'])} operations, "
                   f"{len(ssg['centrings'])} centrings"
                   + (" (with the operations mapping k to -k)." if ssg["minus_k"]
                      else " (group of k only).")]
        for label, xyz in atoms.items():
            symmform = calc_site_symmform(xyz, k_vector, ssg["arrays"])[0]
            ls_out.append(f"- {label}: {format_symmform(symmform)}")
    if len(k_group["irreps"]) > 1:
        ls_out += ["", "A superposition of representations (a cycloid or a helix "
                   "from two modes) is given as irrep = '1,2': the superspace "
                   "group keeps the operations with equal characters."]
    if len(k_group["irreps"]) < rotation.shape[1]:
        ls_out += ["", "Representations of higher dimension are not listed."]
    return "\n".join(ls_out)


def split_labels(text: str) -> list:
    return [label for label in re.split(r"[\s,;]+", text.strip()) if label]


def _rotate_aniso(item, rotation, unit_cell_parameters):
    """U' = D^-1 R D U D R^T D^-1 with D = diag(a*, b*, c*) (CIF U_ij refer
    to the reciprocal axes); the same for B_ij."""
    reciprocal = calc_reciprocal_by_unit_cell_parameters(unit_cell_parameters)[0][:3]
    matrix_d = numpy.diag(reciprocal)
    transform = numpy.linalg.inv(matrix_d) @ rotation @ matrix_d
    for prefix in ("u", "b"):
        names = [f"{prefix}_{i}{j}" for i, j in ("11", "22", "33", "12", "13", "23")]
        if not all(item.is_attribute(name) for name in names):
            continue
        values = [getattr(item, name) for name in names]
        tensor = numpy.array([[values[0], values[3], values[4]],
                              [values[3], values[1], values[5]],
                              [values[4], values[5], values[2]]], dtype=float)
        tensor = transform @ tensor @ transform.transpose()
        for name, (i, j) in zip(names, ((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2))):
            setattr(item, name, float(numpy.round(tensor[i, j], 8)))


def _reset_constraints(item):
    for name in getattr(item, "ATTR_CONSTR_FLAG", ()):
        setattr(item, name, False)


UNSUPPORTED_SPLIT = ("atom_site_moment", "atom_site_susceptibility",
                     "atom_site_exchange", "atom_rho_multipole", "atom_local_axes")


def split_atoms(crystal, k_group: dict) -> dict:
    """New atom_site, atom_site_aniso and atom_site_scat loops: every atom
    whose orbit in the space group falls into several orbits of the group
    of k gives several atoms ('O1', 'O1_2', ...).

    Returns
    -------
    dict
        {loop name: new loop}, and "labels" {old label: [new labels]}.
    """
    unit_cell_parameters = crystal.cell.get_unit_cell_parameters()
    parent_rotation = k_group["parent_rotation"].reshape(3, 3, -1)
    atom_site = crystal.atom_site
    d_new = {"atom_site": [], "atom_site_aniso": [], "atom_site_scat": []}
    d_labels = {}
    for item in atom_site.items:
        xyz = numpy.array([item.fract_x, item.fract_y, item.fract_z], dtype=float)
        representatives = split_atom_orbit(
            xyz, k_group["parent_rotation"], k_group["parent_translation"],
            k_group["rotation"], k_group["translation"], k_group["centring"])
        d_labels[item.label] = []
        for i_new, (position, i_op) in enumerate(representatives):
            label = item.label if i_new == 0 else f"{item.label}_{i_new + 1}"
            d_labels[item.label].append(label)
            new_item = copy.deepcopy(item)
            _reset_constraints(new_item)
            for name in ("multiplicity", "wyckoff_symbol", "space_group_wyckoff",
                         "constr_number"):
                new_item.__dict__.pop(name, None)  # of the old space group
            new_item.label = label
            if i_new > 0:
                for name, value in zip(("fract_x", "fract_y", "fract_z"), position):
                    setattr(new_item, name, float(numpy.round(value, 8)))
            d_new["atom_site"].append(new_item)
            for loop_name in ("atom_site_aniso", "atom_site_scat"):
                if not crystal.is_attribute(loop_name):
                    continue
                for other in getattr(crystal, loop_name).items:
                    if other.label != item.label:
                        continue
                    new_other = copy.deepcopy(other)
                    _reset_constraints(new_other)
                    new_other.label = label
                    if loop_name == "atom_site_aniso" and i_new > 0:
                        _rotate_aniso(new_other, parent_rotation[:, :, i_op],
                                      unit_cell_parameters)
                    d_new[loop_name].append(new_other)
        if len(representatives) > 1:
            for loop_name in UNSUPPORTED_SPLIT:
                if crystal.is_attribute(loop_name) and item.label in [
                        other.label for other in getattr(crystal, loop_name).items
                        if hasattr(other, "label")]:
                    raise CrysPyError(
                        f"Atom '{item.label}' splits into "
                        f"{len(representatives)} atoms in the group of k, and "
                        f"its {loop_name} cannot be converted. Remove "
                        f"{loop_name} of this atom first.")
    out = {"labels": d_labels}
    for loop_name, l_item in d_new.items():
        if l_item:
            loop = copy.deepcopy(getattr(crystal, loop_name))
            loop.items = l_item
            out[loop_name] = loop
    return out


def take_parent_space_group(crystal):
    """SpaceGroup of the crystal or of its parent_space_group, or an error."""
    if crystal.is_attribute("space_group"):
        return crystal.space_group
    if crystal.is_attribute("parent_space_group"):
        if not crystal.parent_space_group.is_cell_of_structure():
            raise CrysPyError(
                "parent_space_group is given in another cell "
                "(transform_Pp_abc or child_transform_Pp_abc is not a,b,c): "
                "its operations cannot be used in the cell of the structure.")
        return crystal.parent_space_group.take_space_group()
    raise CrysPyError(
        "The crystal has neither space_group nor parent_space_group: the "
        "paramagnetic space group is needed.")


def calc_magnetic_domain_matrices(crystal):
    """Matrices (9, n_domain) on hkl of the magnetic domains of a crystal
    described by a magnetic superspace group and its parent space group.

    The domains are the cosets of the point group of the superspace
    operations in the point group of the parent space group; Friedel pairs
    are merged (unpolarized intensities)."""
    if not crystal.is_attribute("space_group_symop_magn_ssg_operation"):
        raise CrysPyError("The crystal has no propagation vector (superspace "
                          "group); use add_propagation_vector.")
    space_group = take_parent_space_group(crystal)
    parent_rotation = take_space_group_operations(space_group)[0]
    sub_rotation = numpy.stack([
        parse_ssg_operation(item.algebraic)["rotation"]
        for item in crystal.space_group_symop_magn_ssg_operation.items], axis=1)
    return calc_domain_matrices(parent_rotation, sub_rotation)


def form_start_fourier(projector, moment: float, cell_angles=(90., 90., 90.)):
    """Start values (6,) of cos x, y, z, sin x, y, z allowed by the site
    symmetry; zeros for ``moment`` = 0.

    The components are along the unit vectors of a, b, c (μB); the largest
    moment over x4 is ``moment`` (μB). A cos-only start (amplitude-modulated
    moments) is taken when the site symmetry allows it.
    ``cell_angles`` are alpha, beta, gamma in degrees (the lengths depend on
    them in a non-orthogonal cell).
    """
    values = numpy.zeros(6)
    if moment == 0. or numpy.allclose(projector, 0.):
        return values
    # a cos-only start (amplitude-modulated) when the symmetry allows it
    values = projector @ numpy.array([1., 1., 1., 0., 0., 0.])
    if numpy.linalg.norm(values) < 1e-8:
        values = projector @ numpy.ones(6)
    if numpy.linalg.norm(values) < 1e-8:
        i_col = int(numpy.argmax(numpy.linalg.norm(projector, axis=0)))
        values = projector[:, i_col]
    cos_alpha, cos_beta, cos_gamma = numpy.cos(numpy.radians(cell_angles))
    # metric of the unit vectors along a, b, c
    metric = numpy.array([[1., cos_gamma, cos_beta],
                          [cos_gamma, 1., cos_alpha],
                          [cos_beta, cos_alpha, 1.]])
    # largest moment over x4: M = cos cos(t) + sin sin(t) is an ellipse
    gram = numpy.array([[values[:3] @ metric @ values[:3], values[:3] @ metric @ values[3:]],
                        [values[:3] @ metric @ values[3:], values[3:] @ metric @ values[3:]]])
    amplitude = numpy.sqrt(numpy.linalg.eigvalsh(gram).max())
    return numpy.round(moment * values / amplitude, 6)


def form_superspace_items(k_vector, ssg: dict) -> list:
    """cell_wave_vector, superspace operations and centrings, and
    atom_site_Fourier_wave_vector loops of a superspace group."""
    from cryspy.C_item_loop_classes.cl_1_cell_wave_vector import (
        CellWaveVector, CellWaveVectorL)
    from cryspy.C_item_loop_classes.cl_1_space_group_symop_magn_ssg_operation import (
        SpaceGroupSymopMagnSsgOperation, SpaceGroupSymopMagnSsgOperationL)
    from cryspy.C_item_loop_classes.cl_1_space_group_symop_magn_ssg_centering import (
        SpaceGroupSymopMagnSsgCentering, SpaceGroupSymopMagnSsgCenteringL)
    from cryspy.C_item_loop_classes.cl_1_atom_site_fourier_wave_vector import (
        AtomSiteFourierWaveVector, AtomSiteFourierWaveVectorL)
    loop_k = CellWaveVectorL()
    loop_k.items = [CellWaveVector(
        seq_id=1, x=float(k_vector[0]), y=float(k_vector[1]), z=float(k_vector[2]))]
    loop_op = SpaceGroupSymopMagnSsgOperationL()
    loop_op.items = [SpaceGroupSymopMagnSsgOperation(id=str(i + 1), algebraic=text)
                     for i, text in enumerate(ssg["operations"])]
    loop_centring = SpaceGroupSymopMagnSsgCenteringL()
    loop_centring.items = [SpaceGroupSymopMagnSsgCentering(id=str(i + 1), algebraic=text)
                           for i, text in enumerate(ssg["centrings"])]
    loop_harmonic = AtomSiteFourierWaveVectorL()
    loop_harmonic.items = [AtomSiteFourierWaveVector(seq_id=1, q1_coeff=1)]
    return [loop_k, loop_op, loop_centring, loop_harmonic]


def form_fourier_items(rows: list):
    """atom_site_moment_Fourier loop from rows (label, values (6,),
    symmforms (6,) and optionally refinement flags (6,)): one row per atom
    and axis."""
    from cryspy.C_item_loop_classes.cl_1_atom_site_moment_fourier import (
        AtomSiteMomentFourier, AtomSiteMomentFourierL)
    l_item = []
    for row in rows:
        label, values, symmform = row[:3]
        flags = row[3] if len(row) > 3 else [False] * 6
        for i_axis, axis in enumerate("xyz"):
            l_item.append(AtomSiteMomentFourier(
                atom_site_label=label, axis=axis, wave_vector_seq_id=1,
                cos=float(values[i_axis]), sin=float(values[3 + i_axis]),
                cos_symmform=symmform[i_axis], sin_symmform=symmform[3 + i_axis],
                cos_refinement=bool(flags[i_axis]),
                sin_refinement=bool(flags[3 + i_axis])))
    loop = AtomSiteMomentFourierL()
    loop.items = l_item
    return loop


def form_parent_space_group(space_group):
    """ParentSpaceGroup item of a SpaceGroup."""
    from cryspy.C_item_loop_classes.cl_1_parent_space_group import ParentSpaceGroup
    parent = ParentSpaceGroup(name_hm_alt=space_group.name_hm_alt)
    if space_group.is_attribute("it_number"):
        parent.it_number = space_group.it_number
    if space_group.is_attribute("it_coordinate_system_code"):
        parent.it_coordinate_system_code = space_group.it_coordinate_system_code
    return parent


def convert_magnetic_space_group(crystal) -> tuple:
    """Superspace description (k = 0) of a crystal with a magnetic space group.

    The operations ``x,y,z,theta`` of space_group_symop_magn_operation and
    _centering become ``x1,x2,x3,x4,theta`` (k = 0, epsilon = +1, tau = 0)
    and the moments of atom_site_moment become the cos components of
    atom_site_moment_Fourier (sin = 0): for k = 0, ``x4 = k.r = 0`` and
    ``M = Mcos``. Symmforms are taken from atom_site_moment.symmform or
    calculated from the site symmetry.

    Returns
    -------
    tuple
        (list of new items and loops, report in markdown).
    """
    if not crystal.is_attribute("space_group_symop_magn_operation"):
        raise CrysPyError("The crystal has no magnetic space group "
                          "(space_group_symop_magn_operation).")
    if crystal.is_attribute("space_group_symop_magn_ssg_operation"):
        raise CrysPyError("The crystal already has a superspace group.")

    def to_algebraic(xyz):
        parts = [part.strip() for part in xyz.split(",")]
        if len(parts) != 4:
            raise CrysPyError(f"Operation '{xyz}' must be 'x,y,z,+1'.")
        rotation, translation = parse_operation_xyz(",".join(parts[:3]))
        rows = form_operation_xyz(rotation, translation)
        rows = rows.replace("x", "x1").replace("y", "x2").replace("z", "x3")
        theta = "+1" if float(parts[3]) > 0 else "-1"
        return f"{rows},x4,{theta}"

    l_operation = [to_algebraic(item.xyz)
                   for item in crystal.space_group_symop_magn_operation.items]
    l_centring = [to_algebraic(item.xyz)
                  for item in crystal.space_group_symop_magn_centering.items] \
        if crystal.is_attribute("space_group_symop_magn_centering") else []
    if not l_centring:
        l_centring = ["x1,x2,x3,x4,+1"]
    k_vector = numpy.zeros(3)
    ssg = {"operations": l_operation, "centrings": l_centring}
    arrays = calc_superspace_operation_arrays(l_operation, l_centring, k_vector)
    rows = []
    if crystal.is_attribute("atom_site_moment"):
        for item in crystal.atom_site_moment.items:
            site = crystal.atom_site[item.label]
            position = numpy.array([site.fract_x, site.fract_y, site.fract_z])
            values = numpy.zeros(6)
            values[:3] = [getattr(item, f"crystalaxis_{axis}")
                          if item.is_attribute(f"crystalaxis_{axis}") else 0.
                          for axis in "xyz"]
            if item.is_attribute("symmform") and len(item.symmform.split(",")) == 3:
                symmform = [part.strip() for part in item.symmform.split(",")] + ["0"] * 3
            else:
                symmform = calc_site_symmform(position, k_vector, arrays)[0]
                symmform = symmform[:3] + ["0"] * 3
            flags = [bool(getattr(item, f"crystalaxis_{axis}_refinement"))
                     for axis in "xyz"] + [False] * 3
            rows.append((item.label, values, symmform, flags))
    l_new = form_superspace_items(k_vector, ssg)
    if rows:
        l_new.append(form_fourier_items(rows))
    report = (f"Magnetic space group ({len(l_operation)} operations, "
              f"{len(l_centring)} centrings) written as the superspace group "
              f"with k = 0; {len(rows)} moments became cos components of "
              "atom_site_moment_Fourier.")
    return l_new, report


def _take_commensurate_denominator(k_vector, max_denominator: int = 12):
    """q with k.L in (1/q) Z for every lattice vector L, or None."""
    denominator = 1
    for value in k_vector:
        fraction = Fraction(float(value)).limit_denominator(max_denominator)
        if abs(float(fraction) - value) > 1e-8:
            return None
        denominator = denominator * fraction.denominator // gcd(denominator, fraction.denominator)
    return denominator if denominator <= max_denominator else None


def _calc_unit_axes(unit_cell_parameters):
    """Cartesian unit vectors (3, 3) along a, b, c (columns) and the direct
    basis (columns a, b, c in Å); angles in radians."""
    a, b, c, alpha, beta, gamma = unit_cell_parameters
    c_x = c * numpy.cos(beta)
    c_y = c * (numpy.cos(alpha) - numpy.cos(beta) * numpy.cos(gamma)) / numpy.sin(gamma)
    basis = numpy.array([[a, b * numpy.cos(gamma), c_x],
                         [0., b * numpy.sin(gamma), c_y],
                         [0., 0., numpy.sqrt(max(c * c - c_x * c_x - c_y * c_y, 0.))]])
    return basis / numpy.linalg.norm(basis, axis=0)[na, :], basis


def classify_moment_modulation(amplitude, k_vector, unit_cell_parameters):
    """Smallest and largest moment (μB) and the kind of the modulation.

    The moment of the atom in the cell L is ``Re(A exp(-2 pi i k.L)) = a
    cos(t) + b sin(t)`` with ``a = Re A``, ``b = Im A`` and ``t = 2 pi k.L``:
    an ellipse (or a segment) in real space. For a commensurate k only the
    values ``t = 2 pi n / q`` occur.

    Parameters
    ----------
    amplitude : numpy.ndarray
        (3,) complex A along the unit vectors of a, b, c (μB).
    k_vector : numpy.ndarray
        (3,) propagation vector (reciprocal lattice units).
    unit_cell_parameters : numpy.ndarray
        a, b, c (Å), alpha, beta, gamma (radians).

    Returns
    -------
    tuple
        (moment_min, moment_max, kind) with kind such as 'helix',
        'cycloid', 'amplitude-modulated'.
    """
    unit_axes, basis = _calc_unit_axes(unit_cell_parameters)
    vector_a = unit_axes @ amplitude.real
    vector_b = unit_axes @ amplitude.imag
    denominator = _take_commensurate_denominator(k_vector)
    if denominator is not None:
        angle = 2. * numpy.pi * numpy.arange(denominator) / denominator
    else:
        angle = numpy.linspace(0., 2. * numpy.pi, 721)
    moment = numpy.linalg.norm(
        vector_a[:, na] * numpy.cos(angle)[na, :] + vector_b[:, na] * numpy.sin(angle)[na, :],
        axis=0)
    moment_min, moment_max = float(moment.min()), float(moment.max())
    if moment_max < 1e-6:
        return 0., 0., "zero"
    normal = numpy.cross(vector_a, vector_b)
    flag_collinear = numpy.linalg.norm(normal) < 1e-6 * moment_max ** 2
    if denominator is not None and denominator <= 2:
        kind = "collinear, sign alternates" if denominator == 2 else "collinear"
    elif flag_collinear:
        kind = "collinear, amplitude-modulated"
    else:
        k_cart = 2. * numpy.pi * numpy.linalg.inv(basis).transpose() @ numpy.asarray(k_vector)
        cos_angle = abs(normal @ k_cart) / (numpy.linalg.norm(normal) * numpy.linalg.norm(k_cart))
        shape = "circular" if moment_max - moment_min < 1e-3 * moment_max else "elliptic"
        if cos_angle > 0.99:
            kind = f"{shape} helix (rotation plane ⊥ k)"
        elif cos_angle < 0.01:
            kind = f"{shape} cycloid (k in the rotation plane)"
        else:
            kind = (f"{shape} spiral, rotation plane at "
                    f"{numpy.degrees(numpy.arcsin(cos_angle)):.0f}° to k")
    if denominator is not None and denominator > 2:
        kind += f", commensurate (period {denominator} cells)"
    return moment_min, moment_max, kind


def report_moment_modulation(modulation: dict, unit_cell_parameters,
                             max_atoms: int = 60) -> str:
    """Markdown table of the moments of the magnetic atoms of the cell
    (:func:`calc_moment_modulation_by_dictionary`)."""
    k_vector = modulation["k_vector"]
    k_text = ", ".join(f"{float(value):g}" for value in k_vector)
    n_atom = len(modulation["label"])
    ls_out = [f"Magnetic moments for k = ({k_text}): {n_atom} magnetic atoms in the cell.",
              "M(r) = Re(A exp(-2 pi i k.L)) in the cell L; |M| is the range over the cells.",
              "", "| atom | from | x | y | z | A_a, A_b, A_c (μB) | |M| min–max (μB) | modulation |",
              "|---|---|---|---|---|---|---|---|"]
    for i_atom in range(min(n_atom, max_atoms)):
        amplitude = modulation["amplitude"][:, i_atom]
        moment_min, moment_max, kind = classify_moment_modulation(
            amplitude, k_vector, unit_cell_parameters)
        xyz = modulation["fract_xyz"][:, i_atom]
        text_a = ", ".join(form_vector_text([value])[1:-1] for value in amplitude)
        ls_out.append(
            f"| {i_atom + 1} | {modulation['label'][i_atom]} | {xyz[0]:.4f} | "
            f"{xyz[1]:.4f} | {xyz[2]:.4f} | {text_a} | {moment_min:.3f}–{moment_max:.3f} | {kind} |")
    if n_atom > max_atoms:
        ls_out.append(f"[... {n_atom - max_atoms} more atoms are cut ...]")
    return "\n".join(ls_out)
