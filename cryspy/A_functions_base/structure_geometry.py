"""Interatomic distances, bond angles and torsion angles of a crystal.

The geometric quantities of the CIF categories ``_geom_bond``,
``_geom_angle`` and ``_geom_torsion``: every atom is an atom of the
asymmetric unit (``atom_site``) moved by a symmetry element and a lattice
translation. The equivalent position is written as a symmetry operation
'-x+1,y+1/2,-z' (translation included).

Coordinates are fractional; the Cartesian system is x || a*, z || c, in Å.
Unit cell parameters are (a, b, c in Å, alpha, beta, gamma in radians).

Standard uncertainties come from the standard uncertainties of the atomic
coordinates and of the unit cell parameters, taken as independent
(no correlations, as in the programs that have no covariance matrix):
sigma_q^2 = sum_p (dq/dp sigma_p)^2, with numerical derivatives.
A coordinate that two atoms of the quantity share (equivalent atoms of one
atom_site atom) is one parameter, so its correlation is taken into account.

Torsion angle sign: IUPAC convention, positive when, looking from atom 2
to atom 3, the bond 1-2 has to be turned clockwise to eclipse the bond 3-4
(International Tables for Crystallography, Vol. C, 9.4; Klyne & Prelog,
Experientia 16, 521 (1960)).
"""
from fractions import Fraction
from typing import List

import numpy

from cryspy.A_functions_base.unit_cell import (
    calc_m_m_by_unit_cell_parameters, calc_m_inv_m_by_unit_cell_parameters)
from cryspy.A_functions_base.function_1_strings import transform_r_b_to_string

na = numpy.newaxis

# atoms closer than this (Å) are one position (special positions)
TOLERANCE_DISTANCE = 1e-4


def calc_m_m(unit_cell_parameters) -> numpy.ndarray:
    """Matrix (3, 3) from fractional coordinates to Cartesian ones (Å)."""
    m_m, _ = calc_m_m_by_unit_cell_parameters(
        numpy.asarray(unit_cell_parameters, dtype=float))
    return m_m.reshape(3, 3)


def calc_images_near_position(
        unit_cell_parameters, full_symm_elems, fract_atom, fract_center,
        max_distance: float) -> dict:
    """All equivalent positions of an atom near a point.

    Parameters
    ----------
    unit_cell_parameters : numpy.ndarray, shape (6,)
        a, b, c (Å), alpha, beta, gamma (radians).
    full_symm_elems : numpy.ndarray, shape (13, n_symm) or (14, n_symm)
        Symmetry elements with centring: numerators of the translation (3),
        common denominator, rotation matrix (9) [, time reversal].
    fract_atom : numpy.ndarray, shape (3,)
        Fractional coordinates of the atom of the asymmetric unit.
    fract_center : numpy.ndarray, shape (3,)
        Fractional coordinates of the point.
    max_distance : float
        Largest distance from the point, Å.

    Returns
    -------
    dict
        ``fract`` (3, n) positions, ``distance`` (n,) Å, ``symm_index`` (n,)
        symmetry element, ``translation`` (3, n) integer lattice
        translation added after it; sorted by distance, positions at
        zero distance from the point (the point itself) left out.
    """
    symm = numpy.asarray(full_symm_elems, dtype=float)
    b = symm[0:3] / symm[3:4]
    r = symm[4:13].reshape(3, 3, -1)
    fract_atom = numpy.asarray(fract_atom, dtype=float)
    fract_center = numpy.asarray(fract_center, dtype=float)
    m_m = calc_m_m(unit_cell_parameters)
    m_inv_m, _ = calc_m_inv_m_by_unit_cell_parameters(
        numpy.asarray(unit_cell_parameters, dtype=float))
    m_inv_m = m_inv_m.reshape(3, 3)

    # (3, n_symm): the symmetry copies moved next to the point
    xyz = numpy.einsum("ijs,j->is", r, fract_atom) + b
    shift = -numpy.round(xyz - fract_center[:, na])
    # a vector of length max_distance changes the fractional coordinate i
    # by at most max_distance * |row i of M^-1|
    n_cell = numpy.ceil(max_distance * numpy.linalg.norm(m_inv_m, axis=1) + 0.5).astype(int)
    grid = numpy.stack(numpy.meshgrid(
        *[numpy.arange(-n, n + 1) for n in n_cell], indexing="ij"), axis=0).reshape(3, -1)
    # (3, n_symm, n_grid)
    translation = shift[:, :, na] + grid[:, na, :]
    fract = xyz[:, :, na] + translation
    distance = numpy.linalg.norm(
        numpy.einsum("ij,jsg->isg", m_m, fract - fract_center[:, na, na]), axis=0)
    symm_index = numpy.broadcast_to(numpy.arange(xyz.shape[1])[:, na], distance.shape)
    flag = (distance <= max_distance) & (distance > TOLERANCE_DISTANCE)
    fract, translation = fract[:, flag], translation[:, flag]
    distance, symm_index = distance[flag], symm_index[flag]

    order = numpy.lexsort((symm_index, distance))
    l_kept = []
    for i_pos in order:  # one position made by several symmetry elements
        if not any(numpy.linalg.norm(m_m @ (fract[:, i_pos] - fract[:, i_kept]))
                   < TOLERANCE_DISTANCE for i_kept in l_kept):
            l_kept.append(i_pos)
    l_kept = numpy.array(l_kept, dtype=int)
    return {"fract": fract[:, l_kept], "distance": distance[l_kept],
            "symm_index": symm_index[l_kept],
            "translation": numpy.round(translation[:, l_kept]).astype(int)}


def calc_nearest_image(unit_cell_parameters, full_symm_elems, fract_atom,
                       fract_center, l_fract_excluded=()) -> dict:
    """The equivalent position of an atom closest to a point, the point
    itself and the positions ``l_fract_excluded`` (fractional (3,)) left
    out; keys of :func:`calc_images_near_position` for one position."""
    m_m = calc_m_m(unit_cell_parameters)
    # every position has a lattice translation within (a+b+c)/2 of the point;
    # a small sphere is enough in most cases and much faster
    max_distance_all = float(numpy.sum(unit_cell_parameters[:3])) + 0.1
    max_distance = min(4., max_distance_all)
    while True:
        images = calc_images_near_position(
            unit_cell_parameters, full_symm_elems, fract_atom, fract_center, max_distance)
        for i_pos in range(images["distance"].size):
            if not any(numpy.linalg.norm(m_m @ (images["fract"][:, i_pos] - numpy.asarray(excluded)))
                       < TOLERANCE_DISTANCE for excluded in l_fract_excluded):
                return {key: value[..., i_pos] for key, value in images.items()}
        if max_distance >= max_distance_all:
            break
        max_distance = min(2. * max_distance, max_distance_all)
    raise ValueError("the atom has no other equivalent position near the point.")


def form_symmetry_code(full_symm_elems, symm_index: int, translation) -> str:
    """Symmetry operation of an equivalent position with its translation,
    '-x+1,y+1/2,-z'; 'x,y,z' for the atom itself."""
    symm = numpy.asarray(full_symm_elems)
    elem = symm[:, symm_index]
    r = numpy.round(elem[4:13].astype(float)).astype(int).reshape(3, 3)
    b = [Fraction(int(round(elem[i])), int(round(elem[3]))) + int(translation[i])
         for i in range(3)]
    return transform_r_b_to_string(r, b)


def calc_positions_xyz(unit_cell_parameters, atom_fract_xyz, l_position) -> numpy.ndarray:
    """Cartesian positions (3, n_position), Å, of atoms given by
    (atom index, rotation (3, 3), translation (3,) with symmetry element
    translation included): R x + b."""
    fract = numpy.stack([
        rotation @ atom_fract_xyz[:, i_atom] + offset
        for i_atom, rotation, offset in l_position], axis=1)
    return calc_m_m(unit_cell_parameters) @ fract


def calc_distance_xyz(p_1, p_2) -> float:
    """Distance between two Cartesian points, Å."""
    return float(numpy.linalg.norm(p_2 - p_1))


def calc_angle_xyz(p_1, p_2, p_3) -> float:
    """Angle p_1-p_2-p_3 at p_2, degrees."""
    v_1, v_2 = p_1 - p_2, p_3 - p_2
    cos_angle = numpy.dot(v_1, v_2) / (numpy.linalg.norm(v_1) * numpy.linalg.norm(v_2))
    return float(numpy.degrees(numpy.arccos(numpy.clip(cos_angle, -1., 1.))))


def calc_torsion_xyz(p_1, p_2, p_3, p_4) -> float:
    """Torsion angle 1-2-3-4, degrees in (-180, 180], IUPAC sign."""
    b_1, b_2, b_3 = p_2 - p_1, p_3 - p_2, p_4 - p_3
    n_1, n_2 = numpy.cross(b_1, b_2), numpy.cross(b_2, b_3)
    m_1 = numpy.cross(n_1, b_2 / numpy.linalg.norm(b_2))
    return float(numpy.degrees(numpy.arctan2(numpy.dot(m_1, n_2), numpy.dot(n_1, n_2))))


GEOMETRY_FUNCTIONS = {2: calc_distance_xyz, 3: calc_angle_xyz, 4: calc_torsion_xyz}


def calc_geometry_value(unit_cell_parameters, atom_fract_xyz, l_position) -> float:
    """Distance (2 positions, Å), angle (3, degrees, at the second one) or
    torsion angle (4, degrees) of the positions of :func:`calc_positions_xyz`."""
    p_xyz = calc_positions_xyz(unit_cell_parameters, atom_fract_xyz, l_position)
    return GEOMETRY_FUNCTIONS[len(l_position)](*p_xyz.T)


def calc_geometry_sigma(unit_cell_parameters, unit_cell_parameters_sigma,
                        atom_fract_xyz, atom_fract_xyz_sigma, l_position) -> float:
    """Standard uncertainty of :func:`calc_geometry_value` from independent
    uncertainties of the cell parameters (6,) (angles in radians) and of
    the coordinates (3, n_atom), by central differences."""
    ucp = numpy.asarray(unit_cell_parameters, dtype=float)
    fract = numpy.asarray(atom_fract_xyz, dtype=float)
    sigma_sq = 0.
    for i_par in range(6):
        sigma = float(unit_cell_parameters_sigma[i_par])
        if sigma > 0.:
            step = numpy.zeros(6)
            step[i_par] = 1e-6
            derivative = (calc_geometry_value(ucp + step, fract, l_position)
                          - calc_geometry_value(ucp - step, fract, l_position)) / 2e-6
            sigma_sq += (derivative * sigma) ** 2
    for i_atom in sorted({position[0] for position in l_position}):
        for i_xyz in range(3):
            sigma = float(atom_fract_xyz_sigma[i_xyz, i_atom])
            if sigma > 0.:
                step = numpy.zeros_like(fract)
                step[i_xyz, i_atom] = 1e-7
                derivative = (calc_geometry_value(ucp, fract + step, l_position)
                              - calc_geometry_value(ucp, fract - step, l_position)) / 2e-7
                sigma_sq += (derivative * sigma) ** 2
    return float(numpy.sqrt(sigma_sq))


def form_position(full_symm_elems, i_atom: int, symm_index: int, translation) -> tuple:
    """(atom index, rotation (3, 3), offset (3,)) of an equivalent position
    for :func:`calc_positions_xyz`."""
    symm = numpy.asarray(full_symm_elems, dtype=float)
    rotation = symm[4:13, symm_index].reshape(3, 3)
    offset = symm[0:3, symm_index] / symm[3, symm_index] + numpy.asarray(translation, dtype=float)
    return (i_atom, rotation, offset)


def identity_index(full_symm_elems) -> int:
    """Index of the symmetry element x,y,z (no translation)."""
    symm = numpy.asarray(full_symm_elems, dtype=float)
    flag = (numpy.all(numpy.isclose(symm[4:13], numpy.eye(3).reshape(9, 1)), axis=0)
            & numpy.all(numpy.isclose(symm[0:3], 0.), axis=0))
    return int(numpy.flatnonzero(flag)[0])


def find_chain_positions(unit_cell_parameters, full_symm_elems, atom_fract_xyz,
                         l_atom_index: List[int]) -> List[tuple]:
    """Positions (atom index, symmetry index, translation) of a chain of
    2–4 atoms: for a distance or a torsion angle the first atom stays where
    atom_site gives it and every next atom is the equivalent position
    closest to the previous one; for an angle the middle atom stays and the
    two others are the closest to it. A position is used once, so
    'O1 Fe1 O1' gives the angle between two different O1 atoms."""
    fract = numpy.asarray(atom_fract_xyz, dtype=float)
    i_identity = identity_index(full_symm_elems)
    zero = numpy.zeros(3, dtype=int)
    l_used = []

    def nearest(i_atom, center):
        image = calc_nearest_image(unit_cell_parameters, full_symm_elems,
                                   fract[:, i_atom], center, l_used)
        l_used.append(image["fract"])
        return (i_atom, int(image["symm_index"]), image["translation"]), image["fract"]

    if len(l_atom_index) == 3:
        center = fract[:, l_atom_index[1]]
        l_used.append(center)
        position_1, _ = nearest(l_atom_index[0], center)
        position_3, _ = nearest(l_atom_index[2], center)
        return [position_1, (l_atom_index[1], i_identity, zero), position_3]
    l_out = [(l_atom_index[0], i_identity, zero)]
    center = fract[:, l_atom_index[0]]
    l_used.append(center)
    for i_atom in l_atom_index[1:]:
        position, center = nearest(i_atom, center)
        l_out.append(position)
    return l_out


def calc_closest_translation(unit_cell_parameters, fract, fract_reference) -> numpy.ndarray:
    """Integer lattice translation (3,) that puts the fractional position
    closest to the reference position (Cartesian distance)."""
    m_m = calc_m_m(unit_cell_parameters)
    t_0 = -numpy.round(fract - fract_reference)
    grid = numpy.stack(numpy.meshgrid(*[numpy.arange(-1, 2)] * 3, indexing="ij"),
                       axis=0).reshape(3, -1)
    translation = t_0[:, na] + grid
    distance = numpy.linalg.norm(
        m_m @ (fract[:, na] + translation - fract_reference[:, na]), axis=0)
    return translation[:, numpy.argmin(distance)]


def calc_restraint_value(unit_cell_parameters, atom_fract_xyz, l_atom_index,
                         l_rotation, l_offset) -> float:
    """Distance (2 atoms, Å), angle (3, degrees) or torsion angle (4,
    degrees) of a restraint.

    The atom i of the chain is R_i x + b_i of the atom l_atom_index[i]
    moved by the lattice translation that puts it closest to the previous
    atom of the chain (to the middle atom for an angle), so that the value
    does not jump when a coordinate passes 0 or 1 in a refinement.

    Parameters
    ----------
    unit_cell_parameters : numpy.ndarray, shape (6,)
        a, b, c (Å), alpha, beta, gamma (radians).
    atom_fract_xyz : numpy.ndarray, shape (3, n_atom)
        Fractional coordinates of the atoms of atom_site.
    l_atom_index : list of int
    l_rotation : list of numpy.ndarray, shape (3, 3)
    l_offset : list of numpy.ndarray, shape (3,)
        Fractional translation of the symmetry element.
    """
    fract = [rotation @ atom_fract_xyz[:, i_atom] + offset
             for i_atom, rotation, offset in zip(l_atom_index, l_rotation, l_offset)]
    if len(fract) == 3:
        for i_pos in (0, 2):
            fract[i_pos] = fract[i_pos] + calc_closest_translation(
                unit_cell_parameters, fract[i_pos], fract[1])
    else:
        for i_pos in range(1, len(fract)):
            fract[i_pos] = fract[i_pos] + calc_closest_translation(
                unit_cell_parameters, fract[i_pos], fract[i_pos - 1])
    p_xyz = calc_m_m(unit_cell_parameters) @ numpy.stack(fract, axis=1)
    return GEOMETRY_FUNCTIONS[len(fract)](*p_xyz.T)


def calc_restraint_difference(value: float, target: float, n_atom: int) -> float:
    """value - target; for a torsion angle (4 atoms) in (-180, 180]
    degrees, so that -179° and 179° differ by 2°."""
    diff = value - target
    if n_atom == 4:
        diff = (diff + 180.) % 360. - 180.
    return diff
