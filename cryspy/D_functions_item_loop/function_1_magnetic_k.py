"""From a paramagnetic crystal to a magnetic structure with a propagation vector.

The functions work on the items and loops of a crystal data block and are
called by the methods ``Crystal.report_irreps_of_k``,
``Crystal.add_propagation_vector`` and ``Crystal.report_magnetic_domains``.
The numerical part is in ``A_functions_base.magnetic_k_symmetry``.
"""
import copy
import re

import numpy

from cryspy.exceptions import CrysPyError
from cryspy.A_functions_base.magnetic_k_symmetry import (
    calc_axial_moment_matrix, calc_domain_matrices, calc_fourier_basis,
    calc_fourier_projector, calc_k_group_flags, calc_one_dimensional_irreps,
    form_operation_xyz, split_atom_orbit)
from cryspy.A_functions_base.structure_factor_magnetic_k import (
    parse_operation_xyz)
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


def report_irreps(k_group: dict, k_vector, atoms: dict) -> str:
    """Markdown report: the operations of the group of k, their
    one-dimensional representations and the Fourier components allowed for
    every magnetic atom.

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
    ls_out += ["Characters of the one-dimensional representations "
               "(the moment phase of the operation is -arg χ / 2π):", "",
               header, "|---" * (len(k_group["irreps"]) + 1) + "|"]
    for s in range(rotation.shape[1]):
        cells = [form_vector_text([chi[s]])[1:-1] for chi in k_group["irreps"]]
        ls_out.append(f"| {form_operation_xyz(rotation[:, s], translation[:, s])} | "
                      + " | ".join(cells) + " |")
    if len(k_group["irreps"]) == 0:
        ls_out.append("No one-dimensional representation: edit moment_uvw and "
                      "phase of space_group_symop_magn_k by hand.")
    if atoms:
        ls_out += ["", "Allowed Fourier components S (μB along a, b, c) of the "
                   "magnetic atoms:", ""]
        for i_irrep, chi in enumerate(k_group["irreps"]):
            phase = form_phase_by_character(chi)
            l_text = []
            for label, xyz in atoms.items():
                basis = calc_fourier_basis(calc_fourier_projector(
                    xyz, rotation, translation, k_group["moment"], phase,
                    k_group["centring"], k_vector))
                if basis.shape[1] == 0:
                    l_text.append(f"{label}: zero")
                elif basis.shape[1] == 3:
                    l_text.append(f"{label}: any direction")
                else:
                    l_text.append(f"{label}: " + " + ".join(
                        form_vector_text(vector) for vector in basis.transpose()))
            ls_out.append(f"- Γ{i_irrep + 1}: " + "; ".join(l_text))
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
        return crystal.parent_space_group.take_space_group()
    raise CrysPyError(
        "The crystal has neither space_group nor parent_space_group: the "
        "paramagnetic space group is needed.")


def calc_magnetic_domain_matrices(crystal):
    """Matrices (9, n_domain) on hkl of the magnetic domains of a crystal
    with a propagation vector and its parent space group."""
    if not crystal.is_attribute("magnetic_structure_k"):
        raise CrysPyError("The crystal has no propagation vector "
                          "(magnetic_structure_k); use add_propagation_vector.")
    space_group = take_parent_space_group(crystal)
    parent_rotation = take_space_group_operations(space_group)[0]
    l_flag = []
    for item in crystal.magnetic_structure_k.items:
        k_vector = numpy.array([item.vector_1, item.vector_2, item.vector_3], dtype=float)
        l_flag.append(calc_k_group_flags(parent_rotation, k_vector))
    flags = numpy.all(l_flag, axis=0)
    return calc_domain_matrices(parent_rotation, parent_rotation[:, flags])
