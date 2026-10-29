"""Atoms, magnetic moments and unit cell edges for a 3D view of a crystal.

Pure numpy: the positions are computed here and any 3D library draws them
(see ``cryspy.D_functions_item_loop.function_1_view_structure``).

Cartesian coordinate system: x || a*, z || c (as everywhere in CrysPy),
in Å. Magnetic moments are in μB.

Atoms can also carry a symmetric tensor in this Cartesian system (atomic
displacement parameters U, Å², or the susceptibility χ, μB/T); it is
rotated to every symmetry-equivalent atom and drawn as an ellipsoid.
"""
import itertools
import math

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
    symm_index : numpy.ndarray, shape (n,)
        Index of the symmetry element that made the atom.
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

    l_xyz, l_index, l_moment, l_symm = [], [], [], []
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
        l_symm.append(numpy.array(kept, dtype=int))
        if atom_moment_fract is not None:
            l_moment.append(moment_symm[:, kept, i_atom])
    fract_xyz = numpy.concatenate(l_xyz, axis=1)
    index_atom = numpy.concatenate(l_index)
    moment_fract = None
    if atom_moment_fract is not None:
        moment_fract = numpy.concatenate(l_moment, axis=1)
    symm_index = numpy.concatenate(l_symm)
    return fract_xyz, index_atom, moment_fract, symm_index


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


def calc_tensors_xyz_by_symmetry(
        unit_cell_parameters, full_symm_elems, symm_index, tensor_xyz):
    """Rotate Cartesian tensors of atoms by their symmetry elements.

    T' = R T R^T with the Cartesian rotation R = M R_fract M^-1. Both the
    displacement tensor U (polar) and the susceptibility χ (axial to
    axial, time even) transform so.

    Parameters
    ----------
    unit_cell_parameters : numpy.ndarray, shape (6,)
        a, b, c (Å), alpha, beta, gamma (radians).
    full_symm_elems : numpy.ndarray, shape (13, n_symm) or (14, n_symm)
    symm_index : numpy.ndarray, shape (n,)
        Symmetry element of every atom.
    tensor_xyz : numpy.ndarray, shape (3, 3, n)
        Tensor of every atom before the symmetry element.

    Returns
    -------
    numpy.ndarray, shape (3, 3, n)
    """
    m_m = calc_m_m_by_unit_cell_parameters(
        numpy.asarray(unit_cell_parameters, dtype=float))[0].reshape(3, 3)
    symm = numpy.asarray(full_symm_elems, dtype=float)
    r_fract = symm[4:13].reshape(3, 3, -1)[:, :, symm_index]
    r_xyz = numpy.einsum("ij,jks,kl->ils", m_m, r_fract, numpy.linalg.inv(m_m))
    return numpy.einsum("ijs,jks,lks->ils", r_xyz, tensor_xyz, r_xyz)


def calc_adp_probability_factor(probability: float = 0.5) -> float:
    """Scale of displacement ellipsoids at a probability level (ORTEP).

    The ellipsoid with semi-axes c * sqrt(U_i) holds the atom with the
    given probability: P(c) = erf(c/√2) - sqrt(2/π) c exp(-c²/2), the
    chi distribution with 3 degrees of freedom; c = 1.5382 for 50 %.
    """
    probability = min(max(float(probability), 1e-6), 1. - 1e-9)
    low, high = 0., 10.
    for _ in range(80):
        c = 0.5 * (low + high)
        p = math.erf(c / math.sqrt(2.)) - math.sqrt(2. / math.pi) * c * math.exp(-0.5 * c * c)
        if p < probability:
            low = c
        else:
            high = c
    return 0.5 * (low + high)


def calc_ellipsoid_matrix_by_adp(u_xyz, probability: float = 0.5):
    """Matrices that turn a unit sphere into displacement ellipsoids.

    Parameters
    ----------
    u_xyz : numpy.ndarray, shape (3, 3, n)
        Cartesian displacement tensors U, Å².
    probability : float
        Probability level of the ellipsoids.

    Returns
    -------
    numpy.ndarray, shape (3, 3, n)
        c * U^(1/2) (negative eigenvalues of non-positive U are set to 0).
    """
    factor = calc_adp_probability_factor(probability)
    u = numpy.moveaxis(numpy.asarray(u_xyz, dtype=float), 2, 0)
    val, vec = numpy.linalg.eigh(0.5 * (u + numpy.swapaxes(u, 1, 2)))
    root = numpy.sqrt(numpy.clip(val, 0., None))
    res = numpy.einsum("nij,nj,nkj->nik", vec, root, vec) * factor
    return numpy.moveaxis(res, 0, 2)


def calc_ellipsoid_matrix_by_susceptibility(chi_xyz, field: float = 1.):
    """Matrices that turn a unit sphere into magnetization ellipsoids.

    The magnetization ellipsoid is the set of moments M = χ B induced by a
    field B of the given magnitude in every direction; its semi-axes are
    |χ_i| B along the principal axes of χ (``local_susceptibility``,
    ``calc_magnetization_ellipsoid_axes``).

    Parameters
    ----------
    chi_xyz : numpy.ndarray, shape (3, 3, n)
        Cartesian susceptibility tensors, μB/T.
    field : float
        Magnitude of the field, T.

    Returns
    -------
    numpy.ndarray, shape (3, 3, n), μB
        (χ^T χ)^(1/2) B: the same surface as χ applied to the sphere.
    """
    chi = numpy.moveaxis(numpy.asarray(chi_xyz, dtype=float), 2, 0)
    q = numpy.einsum("nji,njk->nik", chi, chi)
    val, vec = numpy.linalg.eigh(q)
    root = numpy.sqrt(numpy.clip(val, 0., None))
    res = numpy.einsum("nij,nj,nkj->nik", vec, root, vec) * float(field)
    return numpy.moveaxis(res, 0, 2)


def calc_ellipsoid_mesh(center, matrix, n_theta: int = 12, n_phi: int = 24):
    """Triangle mesh of the ellipsoid center + matrix @ (unit sphere).

    Parameters
    ----------
    center : numpy.ndarray, shape (3,)
    matrix : numpy.ndarray, shape (3, 3)
    n_theta, n_phi : int
        Number of divisions along the polar and the azimuthal angles.

    Returns
    -------
    vertices : numpy.ndarray, shape (3, n_vertex)
    normals : numpy.ndarray, shape (3, n_vertex), unit vectors
    faces : numpy.ndarray, shape (n_face, 3), vertex indices
    """
    theta = numpy.linspace(0., numpy.pi, n_theta + 1)
    phi = numpy.linspace(0., 2. * numpy.pi, n_phi, endpoint=False)
    t, p = numpy.meshgrid(theta, phi, indexing="ij")
    sphere = numpy.stack([numpy.sin(t) * numpy.cos(p),
                          numpy.sin(t) * numpy.sin(p),
                          numpy.cos(t)], axis=0).reshape(3, -1)
    m = numpy.asarray(matrix, dtype=float)
    vertices = numpy.asarray(center, dtype=float)[:, na] + m @ sphere
    # normals of an affine image of a sphere: M^-T n (pinv: flat ellipsoids)
    normals = numpy.linalg.pinv(m).transpose() @ sphere
    norm = numpy.sqrt(numpy.square(normals).sum(axis=0))
    normals = normals / numpy.where(norm > 0., norm, 1.)[na, :]
    faces = []
    for i_t in range(n_theta):
        for i_p in range(n_phi):
            i_1 = i_t * n_phi + i_p
            i_2 = i_t * n_phi + (i_p + 1) % n_phi
            i_3, i_4 = i_1 + n_phi, i_2 + n_phi
            if i_t > 0:
                faces.append((i_1, i_3, i_2))
            if i_t < n_theta - 1:
                faces.append((i_2, i_3, i_4))
    return vertices, normals, numpy.array(faces, dtype=int)


def calc_structure_scene(
        unit_cell_parameters, full_symm_elems, atom_fract_xyz,
        atom_moment_crystalaxis=None, atom_tensor_xyz=None, n_cells=(1, 1, 1),
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
    atom_tensor_xyz : numpy.ndarray, shape (3, 3, n_atom), optional
        Symmetric Cartesian tensor of every atom (U in Å² or χ in μB/T);
        NaN for atoms without one.
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
        ``tensor_xyz`` (3, 3, n) Cartesian tensors or None,
        ``edges_xyz`` (2, 3, n_edge) Cartesian in Å,
        ``axes_xyz`` (3, 3) Cartesian vectors a, b, c as columns.
    """
    ucp = numpy.asarray(unit_cell_parameters, dtype=float)
    moment_fract = None
    if atom_moment_crystalaxis is not None:
        moment_fract = (numpy.asarray(atom_moment_crystalaxis, dtype=float) /
                        ucp[:3, na])
    fract_cell, index_cell, moment_cell, symm_cell = calc_atoms_in_unit_cell(
        full_symm_elems, atom_fract_xyz, moment_fract, tolerance=tolerance)
    fract_block, index_block = calc_atoms_in_cells(
        fract_cell, n_cells, flag_boundary=flag_boundary, tolerance=tolerance)

    moment_xyz = None
    if moment_cell is not None:
        moment_xyz = calc_xyz_by_fract(ucp, moment_cell[:, index_block])
    tensor_xyz = None
    if atom_tensor_xyz is not None:
        tensor = numpy.asarray(atom_tensor_xyz, dtype=float)
        tensor_xyz = calc_tensors_xyz_by_symmetry(
            ucp, full_symm_elems, symm_cell[index_block],
            tensor[:, :, index_cell[index_block]])
    edges = calc_cell_edges(n_cells, flag_inner=flag_inner_edges)
    return {
        "atom_index": index_cell[index_block],
        "atom_fract_xyz": fract_block,
        "atom_xyz": calc_xyz_by_fract(ucp, fract_block),
        "moment_xyz": moment_xyz,
        "tensor_xyz": tensor_xyz,
        "edges_xyz": calc_xyz_by_fract(ucp, edges.transpose(1, 0, 2)
                                       ).transpose(1, 0, 2),
        "axes_xyz": calc_xyz_by_fract(ucp, numpy.eye(3)),
    }
