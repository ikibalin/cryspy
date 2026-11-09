"""Symmetry of magnetic structures with a propagation vector k.

Helpers to go from the paramagnetic space group to the description used by
CrysPy (``space_group_symop_magn_k``): the group of k, its one-dimensional
irreducible representations, the symmetry-allowed Fourier components of an
atom, the splitting of the atoms into the orbits of the group of k and the
magnetic domains.

Conventions (as in ``structure_factor_magnetic_k``): operations ``{R|t}``
act on fractional coordinates, ``r' = R r + t``; the rotation is stored as
9 components along axis 0 (row-major). The Fourier component of the atom
``R_s r + t_s`` is ``S_s = M_s S exp(-2 pi i phase_s)``; for a lattice
(or centring) translation L, ``S(r + L) = S(r) exp(-2 pi i k.L)``. Moments
are in μB along a/|a|, b/|b|, c/|c|; the symmetry operations only mix axes
of equal length, so an axial vector transforms with ``det(R) R``.

Reference: Yu. A. Izyumov, V. E. Naish, R. P. Ozerov, "Neutron Diffraction
of Magnetic Materials" (1991), ch. 2 (group of k, small representations);
J. Rodriguez-Carvajal, FullProf manual (MSYM, Jbt = 10).
"""
from fractions import Fraction
from itertools import product

import numpy

na = numpy.newaxis

TOLERANCE = 1e-5


def _is_integer(vector, tolerance: float = TOLERANCE) -> bool:
    vector = numpy.asarray(vector, dtype=float)
    return bool(numpy.all(numpy.abs(vector - numpy.round(vector)) < tolerance))


def calc_axial_moment_matrix(op_rotation):
    """``det(R) R`` for every rotation: the action on an axial vector.

    Parameters
    ----------
    op_rotation : numpy.ndarray
        (9, n_s) rotations in fractional coordinates.

    Returns
    -------
    numpy.ndarray
        (9, n_s) matrices acting on the moment components.
    """
    rotation = numpy.asarray(op_rotation, dtype=float).reshape(3, 3, -1)
    det = numpy.linalg.det(rotation.transpose(2, 0, 1))
    return (rotation * numpy.round(det)[na, na, :]).reshape(9, -1)


def _form_linear_expression(row, names, translation: float = 0.) -> str:
    terms = []
    for coefficient, name in zip(row, names):
        coefficient = int(round(coefficient))
        if coefficient == 0:
            continue
        sign = "-" if coefficient < 0 else "+"
        value = "" if abs(coefficient) == 1 else str(abs(coefficient))
        terms.append(f"{sign}{value}{name}")
    fraction = Fraction(float(translation)).limit_denominator(24)
    if fraction != 0:
        sign = "-" if fraction < 0 else "+"
        terms.append(f"{sign}{abs(fraction)}")
    text = "".join(terms) if terms else "0"
    return text[1:] if text.startswith("+") else text


def form_operation_xyz(rotation, translation) -> str:
    """'-y,x-y,z+1/2' for a rotation (9,) and a translation (3,)."""
    rotation = numpy.asarray(rotation, dtype=float).reshape(3, 3)
    return ",".join(_form_linear_expression(rotation[i], "xyz", translation[i])
                    for i in range(3))


def form_moment_uvw(matrix) -> str:
    """'-v,u,w' for a matrix (9,) acting on the moment components."""
    matrix = numpy.asarray(matrix, dtype=float).reshape(3, 3)
    return ",".join(_form_linear_expression(matrix[i], "uvw") for i in range(3))


def calc_k_group_flags(op_rotation, k_vector, tolerance: float = TOLERANCE):
    """Operations of the group of k: ``R^T k ≡ k`` modulo the reciprocal
    lattice (k transforms as a reciprocal vector, ``H' = R^T H``... with
    ``H.(R r) = (R^T H).r``).

    Parameters
    ----------
    op_rotation : numpy.ndarray
        (9, n_s) rotations.
    k_vector : numpy.ndarray
        (3,) propagation vector in the reciprocal basis.

    Returns
    -------
    numpy.ndarray
        (n_s,) bool.
    """
    rotation = numpy.asarray(op_rotation, dtype=float).reshape(3, 3, -1)
    k_vector = numpy.asarray(k_vector, dtype=float)
    k_rotated = (rotation * k_vector[:, na, na]).sum(axis=0)  # R^T k
    difference = k_rotated - k_vector[:, na]
    return numpy.all(numpy.abs(difference - numpy.round(difference)) < tolerance,
                     axis=0)


def calc_multiplication_table(op_rotation, op_translation, centring_translations):
    """Product of operations modulo lattice and centring translations.

    ``g_i g_j = {I|L_ij} g_m`` with ``m = index[i, j]`` and the lattice
    translation ``L_ij`` (integer plus centring), which gives the phase
    ``exp(-2 pi i k.L_ij)`` of the Fourier components.

    Returns
    -------
    index : numpy.ndarray
        (n_s, n_s) int.
    lattice : numpy.ndarray
        (3, n_s, n_s) translations L_ij.
    """
    rotation = numpy.asarray(op_rotation, dtype=float).reshape(3, 3, -1)
    translation = numpy.asarray(op_translation, dtype=float)
    centring = numpy.asarray(centring_translations, dtype=float)
    n_s = rotation.shape[2]
    index = numpy.zeros((n_s, n_s), dtype=int)
    lattice = numpy.zeros((3, n_s, n_s), dtype=float)
    for i in range(n_s):
        for j in range(n_s):
            r_ij = rotation[:, :, i] @ rotation[:, :, j]
            t_ij = rotation[:, :, i] @ translation[:, j] + translation[:, i]
            for m in range(n_s):
                if not numpy.allclose(r_ij, rotation[:, :, m]):
                    continue
                shift = t_ij - translation[:, m]
                if any(_is_integer(shift - centring[:, c])
                       for c in range(centring.shape[1])):
                    index[i, j], lattice[:, i, j] = m, shift
                    break
            else:
                raise ValueError(
                    "The operations do not form a group: the product of "
                    f"operations {i + 1} and {j + 1} is not among them.")
    return index, lattice


def _find_generators(index):
    n_s = index.shape[0]
    generators, reached = [], {0}
    while len(reached) < n_s:
        generator = min(set(range(n_s)) - reached)
        generators.append(generator)
        reached = {0}
        changed = True
        while changed:
            new = {int(index[a, g]) for a in reached for g in generators} | reached
            changed = new != reached
            reached = new
    return generators


def calc_one_dimensional_irreps(
        op_rotation, op_translation, centring_translations, k_vector):
    """One-dimensional irreducible (small) representations of the group of k.

    A representation gives ``chi_s`` for every operation such that the
    Fourier components ``S_s = chi_s det(R_s) R_s S`` are consistent:
    ``chi_i chi_j = chi_m exp(-2 pi i k.L_ij)`` (see
    :func:`calc_multiplication_table`). For a generator of order n,
    ``chi^n`` is fixed by the translations, so its n roots are the
    candidates; every combination of the generators is checked on the
    whole table. Representations of higher dimension are not found.

    Parameters
    ----------
    op_rotation, op_translation : numpy.ndarray
        (9, n_s), (3, n_s) operations of the group of k; the first one is
        the identity.
    centring_translations : numpy.ndarray
        (3, n_c) centring translations including zero.
    k_vector : numpy.ndarray
        (3,) propagation vector.

    Returns
    -------
    list of numpy.ndarray
        complex (n_s,) characters ``chi_s`` of every representation; the
        moment phase of the operation is ``-arg(chi_s) / (2 pi)``.
    """
    k_vector = numpy.asarray(k_vector, dtype=float)
    index, lattice = calc_multiplication_table(
        op_rotation, op_translation, centring_translations)
    omega = numpy.exp(-2j * numpy.pi * (k_vector[:, na, na] * lattice).sum(axis=0))
    n_s = index.shape[0]
    generators = _find_generators(index)
    l_candidates = []
    for generator in generators:
        power, factor, order = generator, 1. + 0j, 1
        while power != 0:
            factor = factor * omega[power, generator]
            power = int(index[power, generator])
            order += 1
        # chi_g^order = product of omega along the powers of g
        root = factor ** (1. / order)
        l_candidates.append([root * numpy.exp(2j * numpy.pi * m / order)
                             for m in range(order)])
    l_irrep = []
    for values in product(*l_candidates):
        chi = numpy.full(n_s, numpy.nan, dtype=complex)
        chi[0] = 1.
        changed = True
        while changed:
            changed = False
            for a in numpy.flatnonzero(~numpy.isnan(chi)):
                for generator, value in zip(generators, values):
                    m = int(index[a, generator])
                    if numpy.isnan(chi[m]):
                        chi[m] = chi[a] * value / omega[a, generator]
                        changed = True
        if numpy.any(numpy.isnan(chi)):
            continue
        product_chi = chi[:, na] * chi[na, :]
        if not numpy.allclose(product_chi, chi[index] * omega, atol=1e-6):
            continue
        if not any(numpy.allclose(chi, other, atol=1e-6) for other in l_irrep):
            l_irrep.append(chi)
    return l_irrep


def calc_fourier_projector(
        atom_fract_xyz, op_rotation, op_translation, op_moment, op_phase,
        centring_translations, k_vector):
    """Projector on the Fourier components allowed by the site symmetry.

    An operation which keeps the atom in place up to a lattice translation
    L (``R r + t = r + L``) requires
    ``S = M S exp(-2 pi i phase) exp(2 pi i k.L)``; the projector is the
    average of these matrices over the stabilizer.

    Parameters
    ----------
    atom_fract_xyz : numpy.ndarray
        (3,) fractional position.
    op_rotation, op_translation, op_moment, op_phase, centring_translations :
        See ``structure_factor_magnetic_k.calc_f_m_k``.
    k_vector : numpy.ndarray
        (3,) propagation vector.

    Returns
    -------
    numpy.ndarray
        (3, 3) complex projector acting on S (μB, crystal axes).
    """
    xyz = numpy.asarray(atom_fract_xyz, dtype=float)
    k_vector = numpy.asarray(k_vector, dtype=float)
    rotation = numpy.asarray(op_rotation, dtype=float).reshape(3, 3, -1)
    moment = numpy.asarray(op_moment, dtype=float).reshape(3, 3, -1)
    centring = numpy.asarray(centring_translations, dtype=float)
    l_matrix = []
    for s in range(rotation.shape[2]):
        shift = rotation[:, :, s] @ xyz + op_translation[:, s] - xyz
        for c in range(centring.shape[1]):
            if _is_integer(shift - centring[:, c], 1e-4):
                l_matrix.append(moment[:, :, s] * numpy.exp(
                    -2j * numpy.pi * op_phase[s]
                    + 2j * numpy.pi * numpy.dot(k_vector, shift)))
                break
    return numpy.mean(l_matrix, axis=0)


def calc_fourier_basis(projector, tolerance: float = 1e-6):
    """Orthonormal basis (3, n_free) of the allowed Fourier components."""
    u, singular, _ = numpy.linalg.svd(projector)
    basis = u[:, singular > tolerance]
    # make the vectors real when possible and the largest component positive
    for i in range(basis.shape[1]):
        vector = basis[:, i]
        j = int(numpy.argmax(numpy.abs(vector)))
        vector = vector * numpy.exp(-1j * numpy.angle(vector[j]))
        vector[numpy.abs(vector) < 1e-12] = 0.
        basis[:, i] = vector
    return basis


def split_atom_orbit(
        atom_fract_xyz, parent_rotation, parent_translation,
        sub_rotation, sub_translation, centring_translations):
    """Representatives of the orbits of the subgroup in the parent orbit.

    Parameters
    ----------
    atom_fract_xyz : numpy.ndarray
        (3,) position of the atom.
    parent_rotation, parent_translation : numpy.ndarray
        (9, n_p), (3, n_p) operations of the paramagnetic space group.
    sub_rotation, sub_translation : numpy.ndarray
        (9, n_s), (3, n_s) operations of the group of k.
    centring_translations : numpy.ndarray
        (3, n_c) centring translations including zero (common to both).

    Returns
    -------
    list of (numpy.ndarray, int)
        Position in [0, 1) of every representative and the index of the
        parent operation which makes it from the atom (0, the identity,
        for the first one).
    """
    def reduce(position):
        return numpy.mod(numpy.round(position, 8), 1.)

    def orbit(position, rotation, translation):
        rotation = numpy.asarray(rotation, dtype=float).reshape(3, 3, -1)
        l_out = []
        for s in range(rotation.shape[2]):
            for c in range(centring_translations.shape[1]):
                l_out.append(reduce(rotation[:, :, s] @ position + translation[:, s]
                                    + centring_translations[:, c]))
        return l_out

    def contains(l_position, position):
        return any(_is_integer(position - other, 1e-4) for other in l_position)

    xyz = numpy.asarray(atom_fract_xyz, dtype=float)
    rotation = numpy.asarray(parent_rotation, dtype=float).reshape(3, 3, -1)
    l_covered, l_out = [], []
    for p in range(rotation.shape[2]):
        position = reduce(rotation[:, :, p] @ xyz + parent_translation[:, p])
        if contains(l_covered, position):
            continue
        l_out.append((position, p))
        l_covered.extend(orbit(position, sub_rotation, sub_translation))
    return l_out


def calc_domain_matrices(parent_rotation, sub_rotation):
    """Matrices on hkl of the magnetic domains (cosets of the group of k).

    Every left coset ``g G_k`` of the point group of the paramagnetic group
    gives a domain whose intensity at H is the intensity of the first domain
    at ``R_g^T H``. Unpolarized intensities obey Friedel's law,
    ``I(-H) = I(H)``, so cosets related by the inversion are merged.

    Parameters
    ----------
    parent_rotation : numpy.ndarray
        (9, n_p) rotations of the paramagnetic group.
    sub_rotation : numpy.ndarray
        (9, n_s) rotations of the group of k.

    Returns
    -------
    numpy.ndarray
        (9, n_domain) integer matrices ``R_g^T``; the first one is the
        identity.
    """
    parent = numpy.asarray(parent_rotation, dtype=float).reshape(3, 3, -1)
    sub = numpy.asarray(sub_rotation, dtype=float).reshape(3, 3, -1)
    l_sub = [sub[:, :, s] for s in range(sub.shape[2])]
    l_sub += [-matrix for matrix in l_sub]
    l_covered, l_matrix = [], []
    order = numpy.argsort([not numpy.allclose(parent[:, :, p], numpy.eye(3))
                           for p in range(parent.shape[2])], kind="stable")
    for p in order:
        rotation = parent[:, :, p]
        if any(numpy.allclose(rotation, other) for other in l_covered):
            continue
        l_matrix.append(numpy.round(rotation.transpose()).flatten())
        l_covered.extend(rotation @ matrix for matrix in l_sub)
    return numpy.array(l_matrix, dtype=int).transpose()
