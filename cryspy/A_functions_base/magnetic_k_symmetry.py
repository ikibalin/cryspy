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


def calc_minus_k_flags(op_rotation, k_vector, tolerance: float = TOLERANCE):
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
    total = k_rotated + k_vector[:, na]
    flags = numpy.all(numpy.abs(total - numpy.round(total)) < tolerance, axis=0)
    return flags & ~calc_k_group_flags(op_rotation, k_vector, tolerance)


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

