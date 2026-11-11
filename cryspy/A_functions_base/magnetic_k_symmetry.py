"""Symmetry of magnetic structures with a propagation vector k.

Helpers to go from the paramagnetic space group to the magnetic superspace
group used by CrysPy (msCIF, ``space_group_symop_magn_ssg_operation``): the
group of k, its one-dimensional irreducible representations, the
superspace group of a representation (with the operations mapping k to
-k), the symmetry-allowed Fourier components of an atom (symmform), the
splitting of the atoms into orbits and the magnetic domains.

Conventions (as in ``structure_factor_magnetic_k``): operations ``{R|t}``
act on fractional coordinates, ``r' = R r + t``; the rotation is stored as
9 components along axis 0 (row-major). A character ``chi_s`` of a
representation of the group of k gives the moment phase ``phi_s =
-arg(chi_s) / 2 pi`` and the superspace operation ``x4' = h.x + x4 +
tau`` with ``tau = k.t - phi``. Moments are in μB along a/|a|, b/|b|,
c/|c|; the symmetry operations only mix axes of equal length, so an axial
vector transforms with ``det(R) R``.

References: Yu. A. Izyumov, V. E. Naish, R. P. Ozerov, "Neutron Diffraction
of Magnetic Materials" (1991), ch. 2 (group of k, small representations);
J. M. Perez-Mato et al., J. Phys.: Condens. Matter 24 (2012) 163201
(magnetic superspace groups, msCIF); Bilbao MAGNDATA (symmform).
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


def calc_k_group_flags(op_rotation, k_vector, tolerance: float = TOLERANCE,
                       centring_translations=None):
    """Operations of the group of k: ``R^T k ≡ k`` modulo the reciprocal
    lattice (k transforms as a reciprocal vector, ``H' = R^T H``... with
    ``H.(R r) = (R^T H).r``).

    Parameters
    ----------
    op_rotation : numpy.ndarray
        (9, n_s) rotations.
    k_vector : numpy.ndarray
        (3,) propagation vector in the reciprocal basis.
    centring_translations : numpy.ndarray, optional
        (3, n_c) centring translations: ``R^T k - k`` must be a vector of
        the reciprocal lattice of the centred cell, not only integer.

    Returns
    -------
    numpy.ndarray
        (n_s,) bool.
    """
    rotation = numpy.asarray(op_rotation, dtype=float).reshape(3, 3, -1)
    k_vector = numpy.asarray(k_vector, dtype=float)
    k_rotated = (rotation * k_vector[:, na, na]).sum(axis=0)  # R^T k
    return _is_reciprocal(k_rotated - k_vector[:, na], centring_translations, tolerance)


def _is_reciprocal(vector, centring_translations, tolerance):
    """Integer vectors (3, n) with integer products with the centrings."""
    flag = numpy.all(numpy.abs(vector - numpy.round(vector)) < tolerance, axis=0)
    if centring_translations is not None:
        centring = numpy.asarray(centring_translations, dtype=float).reshape(3, -1)
        product_t = centring.transpose() @ vector
        flag = flag & numpy.all(numpy.abs(product_t - numpy.round(product_t)) < tolerance, axis=0)
    return flag


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


SYMBOL_TERM = r"([+-]?)([0-9.]+(?:/[0-9]+)?)?\*?([A-Za-z_][A-Za-z0-9_]*)?"


def parse_symmform(text: str) -> dict:
    """{symbol: coefficient} of a symmform such as '-mxs1', '2mys1+mzs1' or
    '0' (empty dict). Constant terms other than zero are not allowed."""
    import re
    from fractions import Fraction as _Fraction
    pattern = re.compile(SYMBOL_TERM)
    text = str(text).replace(" ", "")
    d_out, position = {}, 0
    while position < len(text):
        match = pattern.match(text, position)
        if match is None or match.end() == position:
            raise ValueError(f"Cannot parse symmform '{text}'.")
        sign = -1. if match.group(1) == "-" else 1.
        value = float(_Fraction(match.group(2))) if match.group(2) else 1.
        if match.group(3):
            d_out[match.group(3)] = d_out.get(match.group(3), 0.) + sign * value
        elif value * sign != 0.:
            raise ValueError(f"Symmform '{text}' has a constant term.")
        position = match.end()
    return d_out


def form_symmform_constraint(l_symmform, scope):
    """Constraint matrix of parameters given by symmforms.

    Every parameter with a symmform of a single symbol (``mxs1`` or
    ``2mxs1``) seen first defines the symbol (it is free, the symbol is its
    value divided by the coefficient); the other parameters are linear
    combinations of the symbols; ``0`` is zero; a parameter without symmform
    (None or '.') is free. Symbols are local to a scope (the atom).

    Parameters
    ----------
    l_symmform : list of str or None
        (n,) symmforms.
    scope : numpy.ndarray
        (n,) scope (atom index) of every parameter.

    Returns
    -------
    constraint : numpy.ndarray
        (n, n) matrix: the effective values are ``constraint @ values``.
    flag_free : numpy.ndarray
        (n,) bool, the parameters that can be refined.
    """
    n_param = len(l_symmform)
    constraint = numpy.zeros((n_param, n_param), dtype=float)
    flag_free = numpy.zeros(n_param, dtype=bool)
    l_parsed = [None if text in (None, ".", "?") else parse_symmform(text)
                for text in l_symmform]
    d_symbol = {}
    for i_param, parsed in enumerate(l_parsed):
        if parsed is not None and len(parsed) == 1:
            symbol, coefficient = list(parsed.items())[0]
            key = (int(scope[i_param]), symbol)
            if key not in d_symbol and coefficient != 0.:
                d_symbol[key] = (i_param, coefficient)
    for i_param, parsed in enumerate(l_parsed):
        if parsed is None:
            constraint[i_param, i_param] = 1.
            flag_free[i_param] = True
            continue
        for symbol, coefficient in parsed.items():
            key = (int(scope[i_param]), symbol)
            if key not in d_symbol:
                raise ValueError(
                    f"Symbol '{symbol}' of a symmform is not defined by a "
                    "parameter of its own.")
            j_param, coefficient_def = d_symbol[key]
            constraint[i_param, j_param] += coefficient / coefficient_def
            if j_param == i_param:
                flag_free[i_param] = True
    return constraint, flag_free


FOURIER_NAMES = ("mxc", "myc", "mzc", "mxs", "mys", "mzs")


def calc_fourier_site_projector(
        atom_fract_xyz, k_vector, op_rotation, op_translation, op_moment,
        op_phase, op_conjugate, tolerance: float = 1e-4):
    """Projector on the Fourier components allowed by the site symmetry.

    The components are ``(Mcos_x, Mcos_y, Mcos_z, Msin_x, Msin_y, Msin_z)``
    of ``M(x4) = Mcos cos(2 pi x4) + Msin sin(2 pi x4)``, ``x4 = k.r``,
    along the unit vectors of a, b, c (msCIF). An operation of the
    superspace group (stored as in ``structure_factor_magnetic_k``) that
    maps the atom at p to ``p + L`` (L a lattice vector) requires, with
    ``A = Mcos + i Msin``,

        A = M A exp(-2 pi i (phi - k.L))               (epsilon = +1),
        A = M A* exp(-2 pi i (phi - 2 k.p - k.L))      (epsilon = -1),

    where M is ``theta det(R) R`` and ``phi = k.t - tau``. The result does
    not depend on the lattice translation chosen for p.

    Parameters
    ----------
    atom_fract_xyz : numpy.ndarray
        (3,) fractional coordinates of the atom.
    k_vector : numpy.ndarray
        (3,) propagation vector (reciprocal lattice units).
    op_rotation, op_translation, op_moment : numpy.ndarray
        (9, n_s), (3, n_s), (9, n_s) operations (centrings included).
    op_phase : numpy.ndarray
        (n_s,) phases ``phi`` in units of 2π.
    op_conjugate : numpy.ndarray
        (n_s,) bool, epsilon = -1.
    tolerance : float
        Tolerance to decide that ``R p + t - p`` is a lattice vector.

    Returns
    -------
    numpy.ndarray
        (6, 6) real projector: allowed components are ``projector @ v``.
    """
    position = numpy.asarray(atom_fract_xyz, dtype=float)
    k_vector = numpy.asarray(k_vector, dtype=float)
    rotation = op_rotation.reshape(3, 3, -1)
    moment = op_moment.reshape(3, 3, -1)
    projector = numpy.zeros((6, 6), dtype=float)
    n_stabilizer = 0
    for i_s in range(rotation.shape[2]):
        lattice = rotation[:, :, i_s] @ position + op_translation[:, i_s] - position
        if not _is_integer(lattice, tolerance):
            continue
        lattice = numpy.round(lattice)
        conjugate = bool(op_conjugate[i_s])
        psi = op_phase[i_s] - numpy.dot(k_vector, lattice)
        if conjugate:
            psi -= 2. * numpy.dot(k_vector, position)
        a_value, b_value = numpy.cos(2. * numpy.pi * psi), numpy.sin(2. * numpy.pi * psi)
        sign = -1. if conjugate else 1.
        # (c, s) -> M (a c + sign b s, -b c + sign a s)
        block = numpy.array([[a_value, sign * b_value], [-b_value, sign * a_value]])
        projector += numpy.kron(block, moment[:, :, i_s])
        n_stabilizer += 1
    return projector / n_stabilizer


def _format_coefficient(value: float) -> str:
    from fractions import Fraction as _Fraction
    fraction = _Fraction(float(value)).limit_denominator(12)
    if abs(float(fraction) - value) < 1e-6:
        return str(fraction)
    return f"{value:.6f}"


def form_fourier_symmform(projector, tolerance: float = 1e-6):
    """Symmforms of the six Fourier components allowed by a projector.

    The allowed space (column space of the projector) is written in reduced
    row echelon form; every pivot component defines a symbol (``mxc1``,
    ``mys1``, ... as in MAGNDATA) and the other components are linear
    combinations of the symbols.

    Parameters
    ----------
    projector : numpy.ndarray
        (6, 6) projector of :func:`calc_fourier_site_projector`.
    tolerance : float
        Numbers smaller than it are zero.

    Returns
    -------
    list of str
        (6,) symmforms of ``cos_x, cos_y, cos_z, sin_x, sin_y, sin_z``.
    """
    matrix = numpy.array(projector, dtype=float).transpose().copy()
    pivots, i_row = [], 0
    for i_col in range(6):
        if i_row >= 6:
            break
        i_max = i_row + int(numpy.argmax(numpy.abs(matrix[i_row:, i_col])))
        if abs(matrix[i_max, i_col]) < tolerance:
            continue
        matrix[[i_row, i_max]] = matrix[[i_max, i_row]]
        matrix[i_row] = matrix[i_row] / matrix[i_row, i_col]
        for j_row in range(6):
            if j_row != i_row:
                matrix[j_row] = matrix[j_row] - matrix[j_row, i_col] * matrix[i_row]
        pivots.append(i_col)
        i_row += 1
    l_symmform = []
    for i_comp in range(6):
        terms = []
        for i_basis, i_pivot in enumerate(pivots):
            coefficient = matrix[i_basis, i_comp]
            if abs(coefficient) < tolerance:
                continue
            symbol = f"{FOURIER_NAMES[i_pivot]}1"
            if abs(coefficient - 1.) < tolerance:
                text = symbol
            elif abs(coefficient + 1.) < tolerance:
                text = f"-{symbol}"
            else:
                text = f"{_format_coefficient(coefficient)}{symbol}"
            if terms and not text.startswith("-"):
                text = "+" + text
            terms.append(text)
        l_symmform.append("".join(terms) if terms else "0")
    return l_symmform


def _normalize_superspace_element(matrix, translation, theta, tolerance):
    """Translation in [0, 1); tau in [0, 1/2) with the time reversal centring
    (x4 + 1/2, theta = -1), which every first-harmonic structure has."""
    translation = numpy.mod(translation, 1.)
    translation[numpy.abs(translation - 1.) < tolerance] = 0.
    translation[numpy.abs(translation) < tolerance] = 0.
    if translation[3] >= 0.5 - tolerance:
        translation[3] = translation[3] - 0.5
        theta = -theta
        if abs(translation[3]) < tolerance:
            translation[3] = 0.
    return matrix, translation, theta


def _key_superspace_element(matrix, translation):
    return tuple(numpy.round(matrix[:3, :3]).astype(int).flatten()) + tuple(
        numpy.round(translation[:3], 5) + 0.)


def form_superspace_group(
        op_rotation, op_translation, op_phase, centring_translations,
        k_vector, minus_rotation=None, minus_translation=None,
        tolerance: float = 1e-5):
    """Magnetic superspace group of a one-dimensional representation.

    The operations ``{R_s|t_s}`` of the group of k with the moment phases
    ``phi_s`` of a one-dimensional representation (``exp(-2 pi i phi_s)`` is
    the character) become ``x4' = h.x + x4 + tau`` with
    ``tau = k.t - phi`` and ``h = R^T k - k`` (msCIF, see
    ``structure_factor_magnetic_k``); the centrings get ``tau = 0`` and the
    time reversal centring ``(x4 + 1/2, -1)`` is added. An operation of the
    space group mapping k to -k (``minus_rotation``, ``minus_translation``)
    is added with epsilon = -1, theta = +1 and tau = 0 (the origin of x4);
    the group is closed and checked. When the closure is not consistent (the
    representation is not equivalent to its conjugate under that
    operation), None is returned.

    Parameters
    ----------
    op_rotation, op_translation : numpy.ndarray
        (9, n_s), (3, n_s) operations of the group of k.
    op_phase : numpy.ndarray
        (n_s,) moment phases (units of 2π) of the representation.
    centring_translations : numpy.ndarray
        (3, n_c) centring translations including zero.
    k_vector : numpy.ndarray
        (3,) propagation vector.
    minus_rotation, minus_translation : numpy.ndarray, optional
        (9,), (3,) an operation with ``R^T k = -k + h``.

    Returns
    -------
    dict or None
        ``operations`` and ``centrings``: lists of algebraic strings such as
        ``'-x1,-x2,x3,x4+1/4,+1'`` (``_space_group_symop_magn_ssg_*``).
    """
    k_vector = numpy.asarray(k_vector, dtype=float)
    rotation = numpy.asarray(op_rotation, dtype=float).reshape(3, 3, -1)

    def element(rot, trans, epsilon, tau, theta):
        matrix = numpy.zeros((4, 4))
        matrix[:3, :3] = rot
        matrix[3, :3] = numpy.round(rot.transpose() @ k_vector - epsilon * k_vector)
        matrix[3, 3] = epsilon
        return _normalize_superspace_element(
            matrix, numpy.concatenate([trans, [tau]]).astype(float), theta, tolerance)

    generators = []
    for i_s in range(rotation.shape[2]):
        tau = numpy.dot(k_vector, op_translation[:, i_s]) - op_phase[i_s]
        generators.append(element(rotation[:, :, i_s], op_translation[:, i_s], 1, tau, 1))
    for i_c in range(centring_translations.shape[1]):
        generators.append(element(numpy.eye(3), centring_translations[:, i_c], 1, 0., 1))
    if minus_rotation is not None:
        generators.append(element(
            numpy.asarray(minus_rotation, dtype=float).reshape(3, 3),
            numpy.asarray(minus_translation, dtype=float), -1, 0., 1))

    group = {}
    queue = list(generators)
    while queue:
        matrix, translation, theta = queue.pop()
        key = _key_superspace_element(matrix, translation)
        if key in group:
            other_matrix, other_translation, other_theta = group[key]
            if (not numpy.allclose(other_matrix, matrix, atol=1e-6)
                    or abs(other_translation[3] - translation[3]) > 1e-4
                    or other_theta != theta):
                return None
            continue
        group[key] = (matrix, translation, theta)
        for g_matrix, g_translation, g_theta in generators:
            queue.append(_normalize_superspace_element(
                g_matrix @ matrix, g_matrix @ translation + g_translation,
                g_theta * theta, tolerance))
            queue.append(_normalize_superspace_element(
                matrix @ g_matrix, matrix @ g_translation + translation,
                g_theta * theta, tolerance))

    identity = numpy.eye(3)
    l_centring = [value for value in group.values()
                  if numpy.allclose(value[0][:3, :3], identity)]
    l_centring.sort(key=lambda value: tuple(value[1]))
    covered, l_operation = set(), []
    for value in sorted(group.values(), key=lambda value: (
            not numpy.allclose(value[0][:3, :3], identity),
            value[0][3, 3] < 0, tuple(value[1]))):
        key = _key_superspace_element(value[0], value[1])
        if key in covered:
            continue
        l_operation.append(value)
        for c_matrix, c_translation, _ in l_centring:
            covered.add(_key_superspace_element(
                *_normalize_superspace_element(
                    c_matrix @ value[0], c_matrix @ value[1] + c_translation,
                    1, tolerance)[:2]))
    names = ("x1", "x2", "x3", "x4")

    def algebraic(matrix, translation, theta):
        rows = [_form_linear_expression(matrix[i], names, translation[i])
                for i in range(4)]
        return ",".join(rows) + (",+1" if theta > 0 else ",-1")

    l_text_centring = []
    for matrix, translation, theta in l_centring:
        l_text_centring.append(algebraic(matrix, translation, theta))
        shifted = translation.copy()
        shifted[3] += 0.5
        l_text_centring.append(algebraic(matrix, shifted, -theta))
    return {"operations": [algebraic(*value) for value in l_operation],
            "centrings": l_text_centring}


def calc_minus_k_flags(op_rotation, k_vector, tolerance: float = TOLERANCE,
                       centring_translations=None):
    """Operations mapping k to -k: ``R^T k ≡ -k`` (and not to k).

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
    k_rotated = (rotation * k_vector[:, na, na]).sum(axis=0)
    flags = _is_reciprocal(k_rotated + k_vector[:, na], centring_translations, tolerance)
    return flags & ~calc_k_group_flags(op_rotation, k_vector, tolerance,
                                       centring_translations)


def snap_k_vector(k_vector, max_denominator: int = 12, tolerance: float = 1e-3):
    """Replace rounded components of k by the fractions they stand for.

    A component closer than ``tolerance`` to a fraction with a denominator
    up to ``max_denominator`` (0.3333 → 1/3) is replaced by it: the group
    of k is decided with a tolerance of 1e-5, so a rounded commensurate k
    would otherwise keep only part of the symmetry.

    Parameters
    ----------
    k_vector : array_like
        (3,) propagation vector in reciprocal lattice units.
    max_denominator : int
        Largest denominator of the fractions.
    tolerance : float
        Largest change of a component (reciprocal lattice units).

    Returns
    -------
    k_snapped : numpy.ndarray
        (3,) propagation vector.
    l_note : list of str
        One line per replaced component, as '0.3333 → 1/3'.
    """
    k_snapped = numpy.array(k_vector, dtype=float)
    l_note = []
    for i_component, value in enumerate(k_snapped):
        fraction = Fraction(float(value)).limit_denominator(max_denominator)
        delta = abs(float(fraction) - value)
        if 1e-12 < delta < tolerance:
            k_snapped[i_component] = float(fraction)
            l_note.append(f"{value:g} → {fraction}")
    return k_snapped, l_note



def calc_projective_regular_representation(op_rotation, op_translation,
                                           centring_translations, k_vector):
    """Projective regular representation of the group of k.

    ``D(g_i) |e_j> = omega_ij |e_m>`` with ``g_i g_j = {I|L_ij} g_m`` and
    ``omega_ij = exp(-2 pi i k.L_ij)`` (:func:`calc_multiplication_table`), so
    that ``D(g_i) D(g_j) = omega_ij D(g_m)``: the factor system of the small
    representations of the space group at k, the same as that of the
    one-dimensional ones (:func:`calc_one_dimensional_irreps`).

    Returns
    -------
    numpy.ndarray
        Complex (n_s, n_s, n_s): ``D[i]`` is the matrix of operation i.
    """
    k_vector = numpy.asarray(k_vector, dtype=float)
    index, lattice = calc_multiplication_table(
        op_rotation, op_translation, centring_translations)
    omega = numpy.exp(-2j * numpy.pi * (k_vector[:, na, na] * lattice).sum(axis=0))
    n_s = index.shape[0]
    matrices = numpy.zeros((n_s, n_s, n_s), dtype=complex)
    for i in range(n_s):
        matrices[i, index[i], numpy.arange(n_s)] = omega[i]
    return matrices


def calc_small_irreps(op_rotation, op_translation, centring_translations,
                      k_vector, seed: int = 0, tolerance: float = 1e-6):
    """All irreducible small representations of the group of k, as matrices.

    The projective regular representation is split by the eigenvectors of a
    random Hermitian matrix of its commutant: an irreducible representation
    of dimension d occurs d times in it, and every eigenvalue of such a
    matrix belongs to one copy, so its eigenspace (of dimension d) is
    irreducible (J. D. Dixon, Numer. Math. 10 (1967) 446). Equivalent
    copies are recognized by their characters; ``sum d^2 = n_s`` is checked.

    Parameters
    ----------
    op_rotation, op_translation : numpy.ndarray
        (9, n_s), (3, n_s) operations of the group of k (the first one is the identity).
    centring_translations : numpy.ndarray
        (3, n_c) centring translations including zero.
    k_vector : numpy.ndarray
        (3,) propagation vector.
    seed : int
        Seed of the random matrix (the result does not depend on it).

    Returns
    -------
    list of numpy.ndarray
        Complex (n_s, d, d) unitary matrices of every representation,
        ``D(g_i) D(g_j) = exp(-2 pi i k.L_ij) D(g_m)``; a one-dimensional
        one is ``chi_s`` of :func:`calc_one_dimensional_irreps`.
    """
    regular = calc_projective_regular_representation(
        op_rotation, op_translation, centring_translations, k_vector)
    n_s = regular.shape[0]
    rng = numpy.random.default_rng(seed)
    random_matrix = rng.normal(size=(n_s, n_s)) + 1j * rng.normal(size=(n_s, n_s))
    commutant = numpy.einsum("sij,jk,slk->il", regular, random_matrix, regular.conj())
    commutant = commutant + commutant.conj().transpose()
    eigenvalue, eigenvector = numpy.linalg.eigh(commutant)
    l_irrep, l_character = [], []
    start = 0
    scale = max(1., float(numpy.abs(eigenvalue).max()))
    while start < n_s:
        stop = start + 1
        while stop < n_s and eigenvalue[stop] - eigenvalue[start] < 1e-7 * scale:
            stop += 1
        basis = eigenvector[:, start:stop]
        matrices = numpy.einsum("ia,sij,jb->sab", basis.conj(), regular, basis)
        character = numpy.trace(matrices, axis1=1, axis2=2)
        if not any(numpy.allclose(character, other, atol=1e-5) for other in l_character):
            l_character.append(character)
            l_irrep.append(matrices)
        start = stop
    if sum(irrep.shape[1] ** 2 for irrep in l_irrep) != n_s:
        raise ValueError("The small representations of the group of k were not "
                         "separated (sum of d^2 differs from the order of the group).")
    l_irrep.sort(key=lambda irrep: irrep.shape[1])
    return l_irrep


def _generate_subgroup(index, generators):
    subgroup, changed = {0}, True
    while changed:
        new = subgroup | {int(index[a, g]) for a in subgroup for g in generators}
        changed = new != subgroup
        subgroup = new
    return frozenset(subgroup)


def calc_stabilizer(matrices, order_parameter, tolerance: float = 1e-6):
    """Operations that keep the order parameter up to a phase.

    Parameters
    ----------
    matrices : numpy.ndarray
        (n_s, d, d) matrices of a representation.
    order_parameter : numpy.ndarray
        (d,) complex order-parameter direction (OPD).

    Returns
    -------
    flags : numpy.ndarray
        (n_s,) bool, ``D(g) eta = lambda_g eta``.
    eigen : numpy.ndarray
        (n_s,) complex lambda_g (NaN where the flag is False).
    """
    eta = numpy.asarray(order_parameter, dtype=complex)
    eta = eta / numpy.linalg.norm(eta)
    image = matrices @ eta  # (n_s, d)
    eigen = image @ eta.conj()
    flags = numpy.linalg.norm(image - eigen[:, na] * eta[na, :], axis=1) < tolerance
    return flags, numpy.where(flags, eigen, numpy.nan)


def calc_isotropy_subgroups(op_rotation, op_translation, centring_translations,
                            k_vector, matrices, seed: int = 0):
    """Isotropy subgroups of a small representation, one per class of domains.

    For every subgroup H of the group of k (generated by one or two
    operations) and every one-dimensional representation lambda of H with
    the factor system of k, the order-parameter directions (OPD) with
    ``D(h) eta = lambda_h eta`` are found; the full stabilizer of a generic
    eta of that space is an isotropy subgroup. Isotropy subgroups whose
    OPDs are related by an operation of the group of k (domains) are
    listed once. The magnetic structure of an OPD is invariant under the
    operations of its isotropy subgroup with the moment phases
    ``-arg(lambda_h) / 2 pi`` (H. T. Stokes and D. M. Hatch, "Isotropy
    Subgroups of the 230 Crystallographic Space Groups" (1988); ISODISTORT).

    Parameters
    ----------
    op_rotation, op_translation, centring_translations, k_vector :
        As in :func:`calc_small_irreps`.
    matrices : numpy.ndarray
        (n_s, d, d) matrices of the representation.

    Returns
    -------
    list of dict
        ``flags`` (n_s,) bool operations of the subgroup, ``eigen`` (n_s,)
        complex lambda, ``order_parameter`` (d,) a generic OPD of it and
        ``dimension`` of its OPD space; the largest subgroups first.
    """
    index, _ = calc_multiplication_table(op_rotation, op_translation, centring_translations)
    n_s, dimension = matrices.shape[0], matrices.shape[1]
    rotation = numpy.asarray(op_rotation, dtype=float)
    translation = numpy.asarray(op_translation, dtype=float)
    l_subgroup = {_generate_subgroup(index, [g]) for g in range(n_s)}
    l_subgroup |= {_generate_subgroup(index, [g, h])
                   for g in range(1, n_s) for h in range(g + 1, n_s)}
    rng = numpy.random.default_rng(seed)
    l_found = []
    for subgroup in sorted(l_subgroup, key=len, reverse=True):
        elements = sorted(subgroup)
        for chi in calc_one_dimensional_irreps(
                rotation[:, elements], translation[:, elements], centring_translations,
                k_vector):
            # OPDs with D(h) eta = chi_h eta: null space of the stacked D(h) - chi_h
            stacked = numpy.concatenate(
                [matrices[h] - chi[i] * numpy.eye(dimension) for i, h in enumerate(elements)])
            _, singular, right = numpy.linalg.svd(stacked)
            n_null = dimension - int(numpy.sum(singular > 1e-6))
            if n_null == 0:
                continue
            space = right[dimension - n_null:].conj().transpose()  # (d, n_null)
            eta = space @ (rng.normal(size=n_null) + 1j * rng.normal(size=n_null))
            flags, eigen = calc_stabilizer(matrices, eta)
            if any(numpy.array_equal(flags, other["flags"])
                   and numpy.allclose(eigen[flags], other["eigen"][flags], atol=1e-5)
                   for other in l_found):
                continue
            l_found.append({"flags": flags, "eigen": eigen, "order_parameter": eta,
                            "dimension": n_null})
    # one subgroup per class of domains: eta and D(g) eta are domains of each other
    l_class = []
    for found in sorted(l_found, key=lambda d: (-int(d["flags"].sum()), d["dimension"])):
        flag_new = True
        for other in l_class:
            if int(other["flags"].sum()) != int(found["flags"].sum()):
                continue
            for g in range(n_s):
                flags, eigen = calc_stabilizer(matrices, matrices[g] @ other["order_parameter"])
                if numpy.array_equal(flags, found["flags"]) and numpy.allclose(
                        eigen[flags], found["eigen"][flags], atol=1e-5):
                    flag_new = False
                    break
            if not flag_new:
                break
        if flag_new:
            l_class.append(found)
    return l_class


def _expand_superspace_elements(l_operation, l_centring, tolerance: float = 1e-5):
    """All elements (products of the centrings and the operations) of a
    superspace group as (matrix (3+d)x(3+d), translation (3+d,) in [0, 1),
    theta)."""
    from cryspy.A_functions_base.structure_factor_magnetic_k import parse_ssg_operation

    def affine(text):
        operation = parse_ssg_operation(text)
        n_wave = operation["epsilon_all"].size
        matrix = numpy.zeros((3 + n_wave, 3 + n_wave))
        matrix[:3, :3] = operation["rotation"].reshape(3, 3)
        matrix[3:, :3] = operation["h_all"]
        matrix[3:, 3:] = numpy.diag(operation["epsilon_all"])
        translation = numpy.concatenate([operation["translation"], operation["tau_all"]])
        return matrix, translation, operation["theta"]

    l_op = [affine(text) for text in l_operation]
    n_dim = l_op[0][0].shape[0]
    l_cen = [affine(text) for text in l_centring] or [
        (numpy.eye(n_dim), numpy.zeros(n_dim), 1)]
    l_element = []
    for c_matrix, c_translation, c_theta in l_cen:
        for o_matrix, o_translation, o_theta in l_op:
            translation = numpy.mod(c_matrix @ o_translation + c_translation, 1.)
            translation[numpy.abs(translation - 1.) < tolerance] = 0.
            l_element.append((c_matrix @ o_matrix, translation, c_theta * o_theta))
    return l_element


def combine_superspace_groups(l_operation_1, l_centring_1, l_operation_2, l_centring_2,
                              tolerance: float = 1e-5) -> dict:
    """Superspace group of two independent sets of modulation waves.

    The group 1 is a (3+d) group of the waves k1..kd, the group 2 a (3+1)
    group of a new wave k_new built on operations of the group 1 (their
    x1, x2, x3 parts). An element {R|t} of both groups acts on every wave
    by its own x(3+n) part, so the elements with the same R and t (modulo
    the lattice) are combined into one (3+d+1) element. The time reversal
    of the two parts must agree: for a first harmonic, theta = -1 is the
    shift of x(3+n) by 1/2, so the x4 part of the group 2 is shifted by 1/2
    where its theta differs from that of the group 1. The elements of the
    group 1 that are not in the group 2 (they do not keep k_new or break
    its representation) are left out.

    Parameters
    ----------
    l_operation_1, l_centring_1 : list of str
        (3+d) superspace operations and centrings of the present waves.
    l_operation_2, l_centring_2 : list of str
        (3+1) superspace operations and centrings of the new wave.

    Returns
    -------
    dict
        ``operations`` and ``centrings``: (3+d+1) algebraic strings
        (x1, ..., x(4+d), time reversal).
    """
    l_element_1 = _expand_superspace_elements(l_operation_1, l_centring_1, tolerance)
    l_element_2 = _expand_superspace_elements(l_operation_2, l_centring_2, tolerance)
    n_wave_1 = l_element_1[0][0].shape[0] - 3
    n_dim = 4 + n_wave_1

    def key_3d(matrix, translation):
        return tuple(numpy.round(matrix[:3, :3]).astype(int).flatten()) + tuple(
            numpy.round(numpy.mod(translation[:3] + tolerance, 1.) - tolerance, 5) + 0.)

    d_element_1 = {}
    for element in l_element_1:
        d_element_1.setdefault(key_3d(element[0], element[1]), []).append(element)
    d_group, l_order = {}, []
    for matrix_2, translation_2, theta_2 in l_element_2:
        for matrix_1, translation_1, theta_1 in d_element_1.get(
                key_3d(matrix_2, translation_2), []):
            matrix = numpy.zeros((n_dim, n_dim))
            matrix[:3 + n_wave_1, :3 + n_wave_1] = matrix_1
            matrix[-1, :3] = matrix_2[3, :3]
            matrix[-1, -1] = matrix_2[3, 3]
            tau_new = translation_2[3] + (0.5 if theta_1 != theta_2 else 0.)
            translation = numpy.mod(numpy.concatenate([translation_1, [tau_new]]), 1.)
            translation[numpy.abs(translation - 1.) < tolerance] = 0.
            key = tuple(numpy.round(matrix).astype(int).flatten()) + tuple(
                numpy.round(translation, 5) + 0.) + (theta_1,)
            if key not in d_group:
                d_group[key] = (matrix, translation, theta_1)
                l_order.append(key)
    l_all = [d_group[key] for key in l_order]
    identity = numpy.eye(n_dim)
    l_centring = [element for element in l_all if numpy.allclose(element[0], identity)]
    l_centring.sort(key=lambda element: (element[2] < 0, tuple(element[1])))

    def key_full(matrix, translation, theta):
        translation = numpy.mod(translation, 1.)
        translation[numpy.abs(translation - 1.) < tolerance] = 0.
        return tuple(numpy.round(matrix).astype(int).flatten()) + tuple(
            numpy.round(translation, 5) + 0.) + (theta,)

    covered, l_operation = set(), []
    for element in sorted(l_all, key=lambda element: (
            not numpy.allclose(element[0][:3, :3], numpy.eye(3)),
            tuple(numpy.diag(element[0])[3:] < 0), element[2] < 0, tuple(element[1]))):
        if key_full(*element) in covered:
            continue
        l_operation.append(element)
        for c_matrix, c_translation, c_theta in l_centring:
            covered.add(key_full(c_matrix @ element[0], c_matrix @ element[1] + c_translation,
                                 c_theta * element[2]))
    names = tuple(f"x{i + 1}" for i in range(n_dim))

    def algebraic(matrix, translation, theta):
        rows = [_form_linear_expression(matrix[i], names, translation[i]) for i in range(n_dim)]
        return ",".join(rows) + (",+1" if theta > 0 else ",-1")

    return {"operations": [algebraic(*element) for element in l_operation],
            "centrings": [algebraic(*element) for element in l_centring]}


def calc_supercell_magnetic_space_group(
        op_rotation, op_translation, n_cell, fract_xyz, atom_kind, moment,
        cell_lengths, tolerance: float = 1e-4) -> dict:
    """Magnetic space group of a commensurate structure in its supercell.

    The candidates are the operations ``{R|t}`` of the (superspace) group
    combined with the lattice translations L of the parent cell inside the
    supercell: in the supercell axes ``R_s = N^-1 R N`` and ``t_s = N^-1
    (t + L)`` (N = diag(n_cell)). A candidate is an operation of the
    supercell group with the time reversal theta when it maps every atom on
    an atom of the same kind and ``M(g r) = theta det(R) R M(r)`` for all
    moments; with theta = -1 only when +1 does not fit. This is the
    t = t0 section of the superspace group (the supercell structure itself
    decides, so any phase of the modulation is treated alike).

    Parameters
    ----------
    op_rotation, op_translation : numpy.ndarray
        (9, n_s) and (3, n_s) operations of the parent cell (fractional),
        centrings included.
    n_cell : numpy.ndarray
        (3,) integer multipliers of a, b, c.
    fract_xyz : numpy.ndarray
        (3, n) fractional positions of all atoms of the supercell
        (supercell axes).
    atom_kind : numpy.ndarray
        (n,) int, atoms with different kinds are never mapped on each other.
    moment : numpy.ndarray
        (3, n) moments along a/|a|, b/|b|, c/|c| (μB); NaN for non-magnetic atoms.
    cell_lengths : numpy.ndarray
        (3,) a, b, c of the supercell (Å).
    tolerance : float
        Tolerance of the positions (fractions of the supercell) and of the
        moments (μB).

    Returns
    -------
    dict
        ``rotation`` (9, n_g) int, ``translation`` (3, n_g) in [0, 1),
        ``theta`` (n_g,) int, and ``representative`` (n,) bool: one atom
        of every orbit of the group (the first one in the given order).
    """
    n_cell = numpy.asarray(n_cell, dtype=float)
    fract_xyz = numpy.asarray(fract_xyz, dtype=float)
    atom_kind = numpy.asarray(atom_kind)
    moment = numpy.asarray(moment, dtype=float)
    flag_magnetic = ~numpy.any(numpy.isnan(moment), axis=0)
    moment_fract = numpy.where(flag_magnetic[na, :], moment, 0.) / numpy.asarray(
        cell_lengths, dtype=float)[:, na]
    lattice = numpy.stack(numpy.meshgrid(
        *[numpy.arange(int(n)) for n in n_cell], indexing="ij"), axis=0).reshape(3, -1)
    l_rotation, l_translation, l_theta, l_image = [], [], [], []
    for i_op in range(op_rotation.shape[1]):
        rotation = op_rotation[:, i_op].reshape(3, 3)
        rotation_s = rotation * n_cell[na, :] / n_cell[:, na]
        if not numpy.allclose(rotation_s, numpy.round(rotation_s), atol=1e-8):
            continue  # R does not keep the supercell lattice
        rotation_s = numpy.round(rotation_s)
        for shift in lattice.transpose():
            translation_s = numpy.mod((op_translation[:, i_op] + shift) / n_cell, 1.)
            translation_s[numpy.isclose(translation_s, 1., atol=1e-8)] = 0.
            xyz_new = rotation_s @ fract_xyz + translation_s[:, na]
            difference = xyz_new[:, :, na] - fract_xyz[:, na, :]
            flag_match = numpy.all(
                numpy.abs(difference - numpy.round(difference)) < tolerance, axis=0)
            flag_match &= atom_kind[:, na] == atom_kind[na, :]
            if not numpy.all(flag_match.sum(axis=1) == 1):
                continue
            image = numpy.argmax(flag_match, axis=1)
            m_new = numpy.linalg.det(rotation_s) * rotation_s @ moment_fract
            m_image = moment_fract[:, image]
            l_theta_ok = [theta for theta in (1, -1) if numpy.allclose(
                theta * m_new[:, flag_magnetic] * numpy.asarray(cell_lengths)[:, na],
                m_image[:, flag_magnetic] * numpy.asarray(cell_lengths)[:, na],
                atol=tolerance)]
            if len(l_theta_ok) == 0:
                continue
            key = (tuple(rotation_s.astype(int).flatten()),
                   tuple(numpy.round(translation_s, 6)), l_theta_ok[0])
            if key in [(tuple(r), tuple(numpy.round(t, 6)), th) for r, t, th in
                       zip(l_rotation, l_translation, l_theta)]:
                continue
            l_rotation.append(rotation_s.astype(int).flatten())
            l_translation.append(translation_s)
            l_theta.append(l_theta_ok[0])
            l_image.append(image)
    representative = numpy.zeros(fract_xyz.shape[1], dtype=bool)
    flag_done = numpy.zeros(fract_xyz.shape[1], dtype=bool)
    for i_atom in range(fract_xyz.shape[1]):
        if flag_done[i_atom]:
            continue
        representative[i_atom] = True
        for image in l_image:
            flag_done[image[i_atom]] = True
    return {"rotation": numpy.stack(l_rotation, axis=1),
            "translation": numpy.stack(l_translation, axis=1),
            "theta": numpy.array(l_theta, dtype=int),
            "representative": representative}
