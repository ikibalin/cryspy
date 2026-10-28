"""Atoms, magnetic moments and unit cell edges for a 3D view of a crystal.

Pure numpy: the positions are computed here and any 3D library draws them
(see ``cryspy.D_functions_item_loop.function_1_view_structure``).

Cartesian coordinate system: x || a*, z || c (as everywhere in CrysPy),
in Å. Magnetic moments are in μB.
"""
import itertools

import numpy

from cryspy.A_functions_base.unit_cell import calc_m_m_by_unit_cell_parameters

na = numpy.newaxis


def calc_atoms_in_unit_cell(
        full_symm_elems, atom_fract_xyz, atom_moment_fract=None,
        tolerance: float = 1e-4):
    """Generate all atoms of the unit cell from the asymmetric unit.

    Parameters
    ----------
    full_symm_elems : numpy.ndarray, shape (13, n_symm) or (14, n_symm)
        Symmetry elements including centring: numerators of the
        translation (3), common denominator, rotation matrix (9) and, for
        magnetic groups, time reversal (+1 or -1).
    atom_fract_xyz : numpy.ndarray, shape (3, n_atom)
        Fractional coordinates of the asymmetric unit.
    atom_moment_fract : numpy.ndarray, shape (3, n_atom), optional
        Magnetic moments in the basis (a, b, c) (μB / Å along each axis);
        NaN for atoms without a moment.
    tolerance : float
        Atoms closer than this (fractional units) are one atom.

    Returns
    -------
    fract_xyz : numpy.ndarray, shape (3, n), in [0, 1)
    index_atom : numpy.ndarray, shape (n,)
        Index of the atom of the asymmetric unit.
    moment_fract : numpy.ndarray, shape (3, n) or None
        Moments transformed as axial vectors: theta * det(R) * R m.
    """
    symm = numpy.asarray(full_symm_elems, dtype=float)
    fract = numpy.asarray(atom_fract_xyz, dtype=float)
    b = symm[0:3] / symm[3:4]
    r = symm[4:13].reshape(3, 3, -1)  # r[i, j, s]
    theta = symm[13] if symm.shape[0] > 13 else numpy.ones(symm.shape[1])

    # (3, n_symm, n_atom)
    xyz = numpy.einsum("ijs,ja->isa", r, fract) + b[:, :, na]
    xyz = numpy.mod(xyz, 1.)
    xyz = numpy.where(numpy.isclose(xyz, 1., atol=tolerance), 0., xyz)

    if atom_moment_fract is not None:
        moment = numpy.asarray(atom_moment_fract, dtype=float)
        det_r = numpy.linalg.det(numpy.moveaxis(r, 2, 0))
        moment_symm = (numpy.einsum("ijs,ja->isa", r, moment) *
                       (theta * det_r)[na, :, na])

    l_xyz, l_index, l_moment = [], [], []
    for i_atom in range(fract.shape[1]):
        xyz_a = xyz[:, :, i_atom]
        kept = []
        for i_s in range(xyz_a.shape[1]):
            diff = xyz_a[:, kept] - xyz_a[:, i_s:i_s+1]
            diff = diff - numpy.round(diff)
            if not numpy.any(numpy.all(numpy.abs(diff) < tolerance, axis=0)):
                kept.append(i_s)
        l_xyz.append(xyz_a[:, kept])
        l_index.append(numpy.full(len(kept), i_atom, dtype=int))
        if atom_moment_fract is not None:
            l_moment.append(moment_symm[:, kept, i_atom])
    fract_xyz = numpy.concatenate(l_xyz, axis=1)
    index_atom = numpy.concatenate(l_index)
    moment_fract = None
    if atom_moment_fract is not None:
        moment_fract = numpy.concatenate(l_moment, axis=1)
    return fract_xyz, index_atom, moment_fract


def calc_atoms_in_cells(
        fract_xyz, n_cells=(1, 1, 1), flag_boundary: bool = True,
        tolerance: float = 1e-4):
    """Translate the atoms of one unit cell to a block of n_a x n_b x n_c cells.

    Parameters
    ----------
    fract_xyz : numpy.ndarray, shape (3, n)
        Fractional coordinates in [0, 1).
    n_cells : tuple of 3 int
        Number of unit cells along a, b and c.
    flag_boundary : bool
        Also show atoms lying on the far faces, edges and corners of the
        block (x = n_a, ...), so that the block looks complete.
    tolerance : float
        Fractional distance to a face to count an atom as lying on it.

    Returns
    -------
    fract_block : numpy.ndarray, shape (3, n_out)
        Fractional coordinates in the block, from 0 to n_cells.
    index : numpy.ndarray, shape (n_out,)
        Index of the atom in ``fract_xyz``.
    """
    fract = numpy.asarray(fract_xyz, dtype=float)
    n_cells = numpy.array(n_cells, dtype=int)
    shifts = numpy.array(list(itertools.product(
        *[range(n_c + 1) for n_c in n_cells])), dtype=float).transpose()
    xyz = fract[:, :, na] + shifts[:, na, :]  # (3, n, n_shift)
    upper = numpy.isclose(xyz, n_cells[:, na, na], atol=tolerance)
    inside = numpy.all((xyz < n_cells[:, na, na] - tolerance) |
                       (upper & flag_boundary), axis=0)
    index, i_shift = numpy.nonzero(inside)
    return xyz[:, index, i_shift], index


def calc_cell_edges(n_cells=(1, 1, 1), flag_inner: bool = False):
    """Edges of the block of unit cells.

    Parameters
    ----------
    n_cells : tuple of 3 int
        Number of unit cells along a, b and c.
    flag_inner : bool
        Also draw the edges of every unit cell inside the block.

    Returns
    -------
    edges : numpy.ndarray, shape (2, 3, n_edge)
        Start and end of every edge in fractional coordinates.
    """
    n_cells = [int(n_c) for n_c in n_cells]
    l_start, l_end = [], []
    for axis in range(3):
        other = [i for i in range(3) if i != axis]
        if flag_inner:
            l_1, l_2 = range(n_cells[other[0]] + 1), range(n_cells[other[1]] + 1)
        else:
            l_1, l_2 = (0, n_cells[other[0]]), (0, n_cells[other[1]])
        for p_1, p_2 in itertools.product(l_1, l_2):
            start = [0., 0., 0.]
            start[other[0]], start[other[1]] = p_1, p_2
            end = list(start)
            end[axis] = n_cells[axis]
            l_start.append(start)
            l_end.append(end)
    return numpy.array([l_start, l_end], dtype=float).transpose(0, 2, 1)


def calc_xyz_by_fract(unit_cell_parameters, fract_xyz):
    """Cartesian coordinates (Å, x || a*, z || c) from fractional ones.

    Parameters
    ----------
    unit_cell_parameters : numpy.ndarray, shape (6,)
        a, b, c (Å), alpha, beta, gamma (radians).
    fract_xyz : numpy.ndarray, shape (3, ...)
    """
    m_m = calc_m_m_by_unit_cell_parameters(
        numpy.asarray(unit_cell_parameters, dtype=float))[0].reshape(3, 3)
    fract = numpy.asarray(fract_xyz, dtype=float)
    return numpy.tensordot(m_m, fract, axes=(1, 0))


def calc_structure_scene(
        unit_cell_parameters, full_symm_elems, atom_fract_xyz,
        atom_moment_crystalaxis=None, n_cells=(1, 1, 1),
        flag_boundary: bool = True, flag_inner_edges: bool = False,
        tolerance: float = 1e-4) -> dict:
    """Everything a 3D view of a crystal structure draws.

    Parameters
    ----------
    unit_cell_parameters : numpy.ndarray, shape (6,)
        a, b, c (Å), alpha, beta, gamma (radians).
    full_symm_elems : numpy.ndarray, shape (13, n_symm) or (14, n_symm)
        Symmetry elements with centring (14 rows: magnetic group with
        time reversal).
    atom_fract_xyz : numpy.ndarray, shape (3, n_atom)
        Fractional coordinates of the asymmetric unit.
    atom_moment_crystalaxis : numpy.ndarray, shape (3, n_atom), optional
        Magnetic moments (μB) along a/|a|, b/|b|, c/|c| (mcif
        ``_atom_site_moment.crystalaxis_x``); NaN for non-magnetic atoms.
    n_cells : tuple of 3 int
        Number of unit cells along a, b and c.
    flag_boundary : bool
        Show atoms on the far faces of the block.
    flag_inner_edges : bool
        Draw the edges of every unit cell, not only of the block.
    tolerance : float
        Fractional distance under which two positions are the same.

    Returns
    -------
    dict
        ``atom_index`` (n,) index in the asymmetric unit,
        ``atom_fract_xyz`` (3, n), ``atom_xyz`` (3, n) Cartesian in Å,
        ``moment_xyz`` (3, n) Cartesian in μB or None,
        ``edges_xyz`` (2, 3, n_edge) Cartesian in Å,
        ``axes_xyz`` (3, 3) Cartesian vectors a, b, c as columns.
    """
    ucp = numpy.asarray(unit_cell_parameters, dtype=float)
    moment_fract = None
    if atom_moment_crystalaxis is not None:
        moment_fract = (numpy.asarray(atom_moment_crystalaxis, dtype=float) /
                        ucp[:3, na])
    fract_cell, index_cell, moment_cell = calc_atoms_in_unit_cell(
        full_symm_elems, atom_fract_xyz, moment_fract, tolerance=tolerance)
    fract_block, index_block = calc_atoms_in_cells(
        fract_cell, n_cells, flag_boundary=flag_boundary, tolerance=tolerance)

    moment_xyz = None
    if moment_cell is not None:
        moment_xyz = calc_xyz_by_fract(ucp, moment_cell[:, index_block])
    edges = calc_cell_edges(n_cells, flag_inner=flag_inner_edges)
    return {
        "atom_index": index_cell[index_block],
        "atom_fract_xyz": fract_block,
        "atom_xyz": calc_xyz_by_fract(ucp, fract_block),
        "moment_xyz": moment_xyz,
        "edges_xyz": calc_xyz_by_fract(ucp, edges.transpose(1, 0, 2)
                                       ).transpose(1, 0, 2),
        "axes_xyz": calc_xyz_by_fract(ucp, numpy.eye(3)),
    }
