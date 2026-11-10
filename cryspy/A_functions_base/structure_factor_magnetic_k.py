"""Magnetic structure factor of a structure with a propagation vector k.

The structure is described as in magnetic superspace groups (mCIF / msCIF,
J. M. Perez-Mato et al., J. Phys.: Condens. Matter 24 (2012) 163201; the
convention of the Bilbao MAGNDATA database and Jana2006). The moment of the
atom at the average position r (fractional, any lattice translation
included) is

    M(x4) = Mcos cos(2 pi x4) + Msin sin(2 pi x4),   x4 = k . r,

with the Fourier components ``Mcos``, ``Msin`` (``_atom_site_moment_Fourier``)
along the unit vectors of a, b, c in μB. Written as a sum over ``±k``,

    M = S exp(-2 pi i k . r) + c.c.,   S = (Mcos + i Msin) / 2,

so the term ``+k`` scatters at ``H = h + k`` and the term ``-k`` at
``H = h - k``. For a general k only one of ``H - k`` and ``H + k`` is
integer; when ``2k`` is integer (k equivalent to -k) both terms scatter at
the same H and add to ``S = Re[(Mcos + i Msin) exp(-2 pi i k.r)]``.

Internally an atom at the position p carries ``S(p) = (Mcos + i Msin)
exp(-2 pi i k.p) / 2`` (this holds for every lattice translation of p).
A superspace operation ``{R, epsilon | t, tau} theta`` (``x' = R x + t``,
``x4' = h.x + epsilon x4 + tau``, ``R^T k = epsilon k + h``) gives the atom
at ``R p + t`` the component

    S' = theta det(R) R S^(epsilon) exp(-2 pi i phi),   phi = k.t - tau,

where ``S^(+1) = S`` and ``S^(-1) = S*`` (the operation maps k to -k). This
follows from ``M'(epsilon x4 + tau + h.p) = theta det(R) R M(x4)``. So an
operation is stored as the rotation and translation acting on positions,
the matrix ``theta det(R) R`` acting on moments, the phase ``phi`` and the
flag of complex conjugation. Centring operations are operations too.

Magnetic reflections appear at ``H = h ± k``. The magnetic structure factor
of the term +k is

    F_M(H) = p sum_j O_j f_j(Q) T_j / n_s sum_s S_js exp(2 pi i H . (R_s r_j + t_s)),

with ``p = 0.2695`` (10^-12 cm per μB); every distinct position is counted
once (``O_j`` is the true occupancy).

All functions work with plain numpy arrays and know nothing about CIF
objects.

Coordinate systems: fractional coordinates for atoms, Miller indices (non
integer for satellites) for reflections, moment components along the unit
vectors of the crystal axes a, b, c (as ``_atom_site_moment.crystalaxis_x``),
and the Cartesian system x || a*, z || c for the output structure factor.
Angles of ``unit_cell_parameters`` are in radians.
"""
import warnings
from fractions import Fraction

import numpy

from cryspy.A_functions_base.unit_cell import (
    calc_eq_ccs_by_unit_cell_parameters,
    calc_m_m_norm_by_unit_cell_parameters,
    calc_sthovl_by_unit_cell_parameters,
)
from cryspy.A_functions_base.magnetic_form_factor import calc_form_factor

na = numpy.newaxis

MAGNETIC_SCATTERING_LENGTH = 0.2695  # 10^-12 cm per μB, gamma r_0 / 2


def parse_operation_xyz(operation_xyz: str):
    """Parse a symmetry operation such as ``"-x+1/2,y,-z"``.

    Parameters
    ----------
    operation_xyz : str
        Three comma separated expressions of ``x``, ``y``, ``z`` with
        rational translations.

    Returns
    -------
    rotation : numpy.ndarray
        Rotation matrix stored as 9 components ``(r11, r12, r13, r21, ...)``.
    translation : numpy.ndarray
        Translation (3,) in fractions of the cell.
    """
    return _parse_linear_expression(operation_xyz, ("x", "y", "z"))


def parse_ssg_operation(algebraic: str) -> dict:
    """Parse a (3+1)-dimensional magnetic superspace operation.

    Parameters
    ----------
    algebraic : str
        mCIF ``_space_group_symop_magn_ssg_operation.algebraic`` such as
        ``"-x1+1/2,x2,-x3,-x4+1/2,-1"``: x1', x2', x3', x4' and the time
        reversal +1 or -1.

    Returns
    -------
    dict
        ``rotation`` (9,), ``translation`` (3,) acting on x1, x2, x3; ``h``
        (3,) coefficients of x1, x2, x3 in x4'; ``epsilon`` (+1 or -1,
        coefficient of x4 in x4'); ``tau`` (translation of x4');
        ``theta`` (+1 or -1).
    """
    parts = [part.strip() for part in algebraic.replace(" ", "").split(",")]
    if len(parts) != 5:
        raise ValueError(
            f"Superspace operation '{algebraic}' must have five components "
            "(x1, x2, x3, x4, time reversal); only one propagation vector "
            "(3+1 dimensions) is supported.")
    names = ("x1", "x2", "x3", "x4")
    matrix, translation = _parse_affine(parts[:4], names)
    epsilon = matrix[3, 3]
    if (not numpy.allclose(matrix[:3, 3], 0.) or abs(abs(epsilon) - 1.) > 1e-8):
        raise ValueError(f"Superspace operation '{algebraic}' is not of the "
                         "form x4' = h.x + (±1) x4 + tau.")
    theta = int(round(float(Fraction(parts[4]))))
    if theta not in (1, -1):
        raise ValueError(f"Time reversal of '{algebraic}' must be +1 or -1.")
    return {"rotation": matrix[:3, :3].flatten(), "translation": translation[:3],
            "h": matrix[3, :3].copy(), "epsilon": int(round(epsilon)),
            "tau": float(translation[3]), "theta": theta}


TERM_PATTERN = r"([+-]?)([0-9.]+(?:/[0-9]+)?)?\*?({names})?"


def calc_superspace_operation_arrays(l_operation, l_centring, k_vector,
                                     tolerance: float = 1e-4) -> dict:
    """Arrays of the operations of a magnetic superspace group.

    Every product of a centring and an operation (algebraic strings of
    ``_space_group_symop_magn_ssg_operation`` and ``_centering``) becomes
    one operation ``{R|t}`` with the matrix ``theta det(R) R`` acting on the
    moments, the phase ``phi = k.t - tau`` and the flag ``epsilon = -1``.

    Parameters
    ----------
    l_operation, l_centring : list of str
        Algebraic strings such as ``'-x1,x2+1/2,-x3,-x4+1/4,+1'``.
    k_vector : numpy.ndarray
        (3,) propagation vector (reciprocal lattice units).
    tolerance : float
        Tolerance of the check ``R^T k = epsilon k + h``.

    Returns
    -------
    dict
        ``rotation`` (9, n_s), ``translation`` (3, n_s), ``moment``
        (9, n_s), ``phase`` (n_s,), ``conjugate`` (n_s,) bool.
    """
    k_vector = numpy.asarray(k_vector, dtype=float)
    l_op = [parse_ssg_operation(text) for text in l_operation]
    l_cen = [parse_ssg_operation(text) for text in l_centring] or [
        parse_ssg_operation("x1,x2,x3,x4,+1")]
    l_rotation, l_translation, l_moment, l_phase, l_conjugate = [], [], [], [], []
    for centring in l_cen:
        if (not numpy.allclose(centring["rotation"], numpy.eye(3).flatten())
                or centring["epsilon"] != 1 or numpy.any(centring["h"] != 0.)):
            raise ValueError("A centring of the superspace group must be a "
                             "pure translation.")
        for operation in l_op:
            rotation = operation["rotation"].reshape(3, 3)
            translation = operation["translation"] + centring["translation"]
            tau = operation["tau"] + centring["tau"]
            theta = operation["theta"] * centring["theta"]
            epsilon = operation["epsilon"]
            if not numpy.allclose(rotation.transpose() @ k_vector,
                                  epsilon * k_vector + operation["h"],
                                  atol=tolerance):
                raise ValueError(
                    "The superspace operations do not agree with the wave "
                    f"vector {k_vector}: R^T k must be epsilon k + h.")
            l_rotation.append(rotation.flatten())
            l_translation.append(translation)
            l_moment.append(theta * numpy.round(numpy.linalg.det(rotation))
                            * rotation.flatten())
            l_phase.append(float(numpy.dot(k_vector, translation) - tau))
            l_conjugate.append(epsilon == -1)
    return {"rotation": numpy.stack(l_rotation, axis=1),
            "translation": numpy.stack(l_translation, axis=1),
            "moment": numpy.stack(l_moment, axis=1),
            "phase": numpy.array(l_phase, dtype=float),
            "conjugate": numpy.array(l_conjugate, dtype=bool)}


def _parse_affine(parts, names):
    """Matrix (n_part, n_name) and translation (n_part,) of linear
    expressions such as '-x1+x4+1/2' or '2x-y'."""
    import re
    pattern = re.compile(TERM_PATTERN.format(
        names="|".join(sorted(names, key=len, reverse=True))))
    matrix = numpy.zeros((len(parts), len(names)), dtype=float)
    translation = numpy.zeros((len(parts),), dtype=float)
    for i_row, part in enumerate(parts):
        position = 0
        if part == "":
            raise ValueError("Empty component of a symmetry operation.")
        while position < len(part):
            match = pattern.match(part, position)
            if match is None or match.end() == position:
                raise ValueError(f"Cannot parse '{part}' of a symmetry operation.")
            sign = -1. if match.group(1) == "-" else 1.
            value = float(Fraction(match.group(2))) if match.group(2) else 1.
            if match.group(3):
                matrix[i_row, names.index(match.group(3))] += sign * value
            elif match.group(2):
                translation[i_row] += sign * value
            else:
                raise ValueError(f"Cannot parse '{part}' of a symmetry operation.")
            position = match.end()
    return matrix, translation


def _parse_linear_expression(expression: str, names: tuple):
    parts = [part.strip() for part in expression.replace(" ", "").split(",")]
    if len(parts) != 3:
        raise ValueError(f"Operation '{expression}' must have three components.")
    matrix, translation = _parse_affine(parts, names)
    return matrix.flatten(), translation


def is_k_equivalent_to_minus_k(k_vector, tolerance: float = 1e-6) -> bool:
    """Return True when ``2k`` is a reciprocal lattice vector."""
    two_k = 2. * numpy.asarray(k_vector, dtype=float)
    return bool(numpy.all(numpy.abs(two_k - numpy.round(two_k)) < tolerance))


def calc_satellite_index_hkl(
        sthovl_min: float, sthovl_max: float, unit_cell_parameters,
        k_vector, centring_translations=None):
    """List magnetic satellites ``H = h + k`` (and ``h - k``) in a range.

    No symmetry reduction is done: every satellite is returned once with
    multiplicity 1, which is what a powder pattern of a single k domain
    needs (all arms of the star give the same powder pattern).

    Parameters
    ----------
    sthovl_min, sthovl_max : float
        Range of sin(theta)/lambda in Å⁻¹.
    unit_cell_parameters : numpy.ndarray
        a, b, c (Å), alpha, beta, gamma (radians).
    k_vector : numpy.ndarray
        Propagation vector (3,) in reciprocal lattice units.
    centring_translations : numpy.ndarray, optional
        Centring translations (3, n_c) of the nuclear cell. Satellites with
        ``h`` forbidden by the centring are removed.

    Returns
    -------
    index_hkl : numpy.ndarray
        (3, n_hkl) float indices of the satellites.
    sign_k : numpy.ndarray
        (n_hkl,) +1 for ``h + k`` and -1 for ``h - k``.
    """
    k_vector = numpy.asarray(k_vector, dtype=float)
    lengths = unit_cell_parameters[:3]
    h_max = numpy.ceil(2. * sthovl_max * lengths).astype(int) + 1
    grid = numpy.meshgrid(
        *[numpy.arange(-h_m, h_m + 1) for h_m in h_max], indexing="ij")
    index_h = numpy.stack([g.flatten() for g in grid], axis=0).astype(float)

    if centring_translations is not None:
        phase = 2. * numpy.pi * (
            index_h[:, :, na] * centring_translations[:, na, :]).sum(axis=0)
        allowed = numpy.abs(numpy.exp(1j * phase).sum(axis=1)) > 1e-6
        index_h = index_h[:, allowed]

    signs = (1,) if is_k_equivalent_to_minus_k(k_vector) else (1, -1)
    l_index, l_sign = [], []
    for sign in signs:
        index_hkl = index_h + sign * k_vector[:, na]
        sthovl = calc_sthovl_by_unit_cell_parameters(
            index_hkl, unit_cell_parameters)[0]
        flag = numpy.logical_and(sthovl >= sthovl_min, sthovl <= sthovl_max)
        l_index.append(index_hkl[:, flag])
        l_sign.append(numpy.full(flag.sum(), sign, dtype=int))
    return numpy.concatenate(l_index, axis=1), numpy.concatenate(l_sign)


def calc_generated_multiplicity(
        atom_fract_xyz, op_rotation, op_translation, centring_translations,
        tolerance: float = 1e-5):
    """Number of distinct positions generated by the operations and centring.

    Parameters
    ----------
    atom_fract_xyz : numpy.ndarray
        (3, n_atom) fractional coordinates.
    op_rotation : numpy.ndarray
        (9, n_s) rotation parts.
    op_translation : numpy.ndarray
        (3, n_s) translation parts.
    centring_translations : numpy.ndarray
        (3, n_c) centring translations (including zero).

    Returns
    -------
    numpy.ndarray
        (n_atom,) integer multiplicities.
    """
    positions = _calc_generated_positions(
        atom_fract_xyz, op_rotation, op_translation, centring_translations)
    multiplicity = []
    for i_atom in range(positions.shape[-1]):
        xyz = numpy.mod(positions[:, :, :, i_atom].reshape(3, -1), 1.)
        xyz = numpy.where(numpy.abs(xyz - 1.) < tolerance, 0., xyz)
        keys = {tuple(numpy.round(xyz[:, i] / tolerance).astype(int))
                for i in range(xyz.shape[1])}
        multiplicity.append(len(keys))
    return numpy.array(multiplicity, dtype=int)


def _calc_generated_positions(
        atom_fract_xyz, op_rotation, op_translation, centring_translations):
    """Positions (3, n_s, n_c, n_atom) without reduction to the first cell."""
    rotation = op_rotation.reshape(3, 3, -1)
    xyz = (rotation[:, :, :, na] * atom_fract_xyz[na, :, na, :]).sum(axis=1)
    xyz = xyz + op_translation[:, :, na]
    return xyz[:, :, na, :] + centring_translations[:, na, :, na]


def calc_f_m_k(
        index_hkl, k_vector, unit_cell_parameters,
        op_rotation, op_translation, op_moment, op_phase,
        centring_translations,
        atom_fract_xyz, atom_occupancy, atom_b_iso, atom_form_factor,
        atom_fourier_crystalaxis, atom_multiplicity=None, op_conjugate=None):
    """Magnetic structure factor F_M(H) of a k structure.

    Parameters
    ----------
    index_hkl : numpy.ndarray
        (3, n_hkl) satellite indices ``H = h + k`` (reciprocal lattice
        units of the nuclear cell).
    k_vector : numpy.ndarray
        (3,) propagation vector. Pass ``-k`` (and conjugated Fourier
        components with ``-op_phase``) for the satellites ``h - k``.
    unit_cell_parameters : numpy.ndarray
        a, b, c (Å), alpha, beta, gamma (radians).
    op_rotation, op_translation : numpy.ndarray
        (9, n_s) and (3, n_s) operations ``{R_s|t_s}`` acting on positions.
    op_moment : numpy.ndarray
        (9, n_s) matrices ``M_s`` acting on the moment components along a,
        b, c.
    op_phase : numpy.ndarray
        (n_s,) phases ``phi_s`` in units of 2π.
    centring_translations : numpy.ndarray
        (3, n_c) centring translations of the nuclear cell including zero.
    atom_fract_xyz : numpy.ndarray
        (3, n_atom) fractional coordinates of the magnetic atoms.
    atom_occupancy, atom_b_iso : numpy.ndarray
        (n_atom,) occupancy and isotropic B (Å²).
    atom_form_factor : numpy.ndarray
        (n_hkl, n_atom) magnetic form factor.
    atom_fourier_crystalaxis : numpy.ndarray
        (3, n_atom) complex Fourier components ``S_kj`` along the unit
        vectors of a, b, c in μB.
    atom_multiplicity : numpy.ndarray, optional
        (n_atom,) number of distinct generated positions; calculated when
        not given.
    op_conjugate : numpy.ndarray, optional
        (n_s,) bool: the operation maps k to -k (epsilon = -1), so it acts
        on the complex conjugate of ``S_kj``.

    Returns
    -------
    numpy.ndarray
        (3, n_hkl) complex magnetic structure factor in the Cartesian
        system x || a*, z || c, in 10^-12 cm.
    """
    k_vector = numpy.asarray(k_vector, dtype=float)
    if atom_multiplicity is None:
        atom_multiplicity = calc_generated_multiplicity(
            atom_fract_xyz, op_rotation, op_translation,
            centring_translations)
    n_s, n_c = op_rotation.shape[1], centring_translations.shape[1]

    # positions (3, n_s, n_c, n_atom)
    positions = _calc_generated_positions(
        atom_fract_xyz, op_rotation, op_translation, centring_translations)
    # phase of the position (n_hkl, n_s, n_c, n_atom)
    phase_position = 2. * numpy.pi * (
        index_hkl[:, :, na, na, na] * positions[:, na, :, :, :]).sum(axis=0)
    # phase of the moment: exp(-2 pi i (phi_s + k.t_c))
    phase_moment = -2. * numpy.pi * (
        op_phase[:, na] + (k_vector[:, na] * centring_translations).sum(axis=0)[na, :])
    exp_total = numpy.exp(1j * (phase_position + phase_moment[na, :, :, na]))
    # sum over centring (n_hkl, n_s, n_atom)
    exp_s = exp_total.sum(axis=2)

    # moments of generated atoms in crystal axes (3, n_s, n_atom)
    moment = op_moment.reshape(3, 3, -1)
    fourier_s = numpy.repeat(
        numpy.asarray(atom_fourier_crystalaxis, dtype=complex)[:, na, :], n_s, axis=1)
    if op_conjugate is not None:
        fourier_s[:, op_conjugate, :] = numpy.conjugate(fourier_s[:, op_conjugate, :])
    s_kjs = (moment[:, :, :, na] * fourier_s[na, :, :, :]).sum(axis=1)

    sthovl = calc_sthovl_by_unit_cell_parameters(index_hkl, unit_cell_parameters)[0]
    debye_waller = numpy.exp(-atom_b_iso[na, :] * numpy.square(sthovl)[:, na])
    weight = (atom_multiplicity * atom_occupancy)[na, :] * atom_form_factor * \
        debye_waller / (n_s * n_c)  # (n_hkl, n_atom)

    # F in crystal axes (3, n_hkl)
    f_m_crystalaxis = (
        s_kjs[:, na, :, :] * exp_s[na, :, :, :] * weight[na, :, na, :]
        ).sum(axis=(2, 3))

    m_norm = calc_m_m_norm_by_unit_cell_parameters(unit_cell_parameters)[0]
    m_norm = m_norm.reshape(3, 3)
    f_m_ccs = (m_norm[:, :, na] * f_m_crystalaxis[na, :, :]).sum(axis=1)
    return MAGNETIC_SCATTERING_LENGTH * f_m_ccs


def calc_f_m_perp_k(index_hkl, unit_cell_parameters, f_m_ccs):
    """Component of F_M perpendicular to the scattering vector.

    Parameters
    ----------
    index_hkl : numpy.ndarray
        (3, n_hkl) satellite indices.
    unit_cell_parameters : numpy.ndarray
        a, b, c (Å), alpha, beta, gamma (radians).
    f_m_ccs : numpy.ndarray
        (3, n_hkl) magnetic structure factor, Cartesian x || a*, z || c.

    Returns
    -------
    numpy.ndarray
        (3, n_hkl) F_M_perp = F_M - e (e . F_M).
    """
    eq_ccs = calc_eq_ccs_by_unit_cell_parameters(index_hkl, unit_cell_parameters)[0]
    return f_m_ccs - eq_ccs * (eq_ccs * f_m_ccs).sum(axis=0)[na, :]


def calc_powder_iint_magnetic_k(
        sthovl_min: float, sthovl_max: float, unit_cell_parameters,
        k_vector, op_rotation, op_translation, op_moment, op_phase,
        centring_translations,
        atom_fract_xyz, atom_occupancy, atom_b_iso, form_factor_function,
        atom_fourier_crystalaxis, op_conjugate=None):
    """Satellites and their unpolarized powder intensities |F_M_perp|².

    Parameters
    ----------
    sthovl_min, sthovl_max : float
        Range of sin(theta)/lambda (Å⁻¹).
    unit_cell_parameters : numpy.ndarray
        a, b, c (Å), alpha, beta, gamma (radians).
    k_vector : numpy.ndarray
        (3,) propagation vector.
    op_rotation, op_translation, op_moment, op_phase : numpy.ndarray
        Operations of the group of k, see :func:`calc_f_m_k`.
    centring_translations : numpy.ndarray
        (3, n_c) centring translations including zero.
    atom_fract_xyz, atom_occupancy, atom_b_iso : numpy.ndarray
        Magnetic atoms, see :func:`calc_f_m_k`.
    form_factor_function : callable
        Takes sin(theta)/lambda (n_hkl,) and returns (n_hkl, n_atom).
    atom_fourier_crystalaxis : numpy.ndarray
        (3, n_atom) complex Fourier components ``S_kj`` (μB).

    Returns
    -------
    index_hkl : numpy.ndarray
        (3, n_hkl) satellite indices.
    iint : numpy.ndarray
        (n_hkl,) |F_M_perp|² in 10^-24 cm², multiplicity 1.
    f_m_perp : numpy.ndarray
        (3, n_hkl) complex F_M_perp.
    """
    index_hkl, _ = calc_satellite_index_hkl(
        sthovl_min, sthovl_max, unit_cell_parameters, k_vector,
        centring_translations)
    f_m_perp = calc_f_m_perp_k_at_index_hkl(
        index_hkl, unit_cell_parameters, k_vector, op_rotation,
        op_translation, op_moment, op_phase, centring_translations,
        atom_fract_xyz, atom_occupancy, atom_b_iso, form_factor_function,
        atom_fourier_crystalaxis, op_conjugate=op_conjugate)
    iint = numpy.square(numpy.abs(f_m_perp)).sum(axis=0)
    return index_hkl, iint, f_m_perp


def calc_f_m_perp_k_at_index_hkl(
        index_hkl, unit_cell_parameters, k_vector, op_rotation,
        op_translation, op_moment, op_phase, centring_translations,
        atom_fract_xyz, atom_occupancy, atom_b_iso, form_factor_function,
        atom_fourier_crystalaxis, op_conjugate=None, tolerance: float = 1e-4):
    """F_M_perp of one propagation vector at arbitrary (fractional) indices.

    The term ``+k`` of the moments scatters at ``H = h + k`` and the term
    ``-k`` (conjugated components) at ``H = h - k`` (h integer). For a
    general k only one of ``H - k`` and ``H + k`` can be integer, so every
    fractional H belongs to one term without ambiguity. When k is
    equivalent to -k (2k integer) both terms scatter at the same H and are
    added. Indices which are not satellites of k get zero.

    Parameters
    ----------
    index_hkl : numpy.ndarray
        (3, n_hkl) float Miller indices H in the reciprocal basis of the
        crystal cell.
    unit_cell_parameters : numpy.ndarray
        a, b, c (Å), alpha, beta, gamma (radians).
    k_vector : numpy.ndarray
        (3,) propagation vector in the reciprocal basis.
    op_rotation, op_translation, op_moment, op_phase, centring_translations,
    atom_fract_xyz, atom_occupancy, atom_b_iso, op_conjugate :
        See :func:`calc_f_m_k`.
    form_factor_function : callable
        Takes sin(theta)/lambda (n,) (Å⁻¹) and returns (n, n_atom).
    atom_fourier_crystalaxis : numpy.ndarray
        (3, n_atom) complex Fourier components ``S`` of the ``+k`` term at
        the positions ``atom_fract_xyz`` (μB, crystal axes).
    tolerance : float
        Tolerance to decide that ``H ∓ k`` is integer.

    Returns
    -------
    f_m_perp : numpy.ndarray
        (3, n_hkl) complex F_M_perp (10^-12 cm) in the Cartesian system
        x || a*, z || c.
    """
    index_hkl = numpy.asarray(index_hkl, dtype=float)
    k_vector = numpy.asarray(k_vector, dtype=float)
    multiplicity = calc_generated_multiplicity(
        atom_fract_xyz, op_rotation, op_translation, centring_translations)
    f_m_ccs = numpy.zeros(index_hkl.shape, dtype=complex)
    flag_any = numpy.zeros(index_hkl.shape[1], dtype=bool)
    for sign in (1, -1):
        index_h = index_hkl - sign * k_vector[:, na]
        flag = numpy.all(
            numpy.abs(index_h - numpy.round(index_h)) < tolerance, axis=0)
        if not numpy.any(flag):
            continue
        flag_any = flag_any | flag
        index_sign = index_hkl[:, flag]
        sthovl = calc_sthovl_by_unit_cell_parameters(
            index_sign, unit_cell_parameters)[0]
        if sign == 1:
            fourier, phase = atom_fourier_crystalaxis, op_phase
        else:
            fourier, phase = numpy.conjugate(atom_fourier_crystalaxis), -op_phase
        f_m_ccs[:, flag] += calc_f_m_k(
            index_sign, sign * k_vector,
            unit_cell_parameters, op_rotation, op_translation, op_moment,
            phase, centring_translations, atom_fract_xyz, atom_occupancy,
            atom_b_iso, form_factor_function(sthovl), fourier,
            atom_multiplicity=multiplicity, op_conjugate=op_conjugate)
    f_m_perp = numpy.zeros(index_hkl.shape, dtype=complex)
    if numpy.any(flag_any):
        f_m_perp[:, flag_any] = calc_f_m_perp_k(
            index_hkl[:, flag_any], unit_cell_parameters, f_m_ccs[:, flag_any])
    return f_m_perp


def calc_full_symm_elems_by_k_operations(
        op_rotation, op_translation, centring_translations):
    """Integer symmetry elements of the crystal from the operations of k.

    The result has the layout of ``full_mcif_elems`` (rows: translation
    numerators (3), common denominator (1), rotation (9), time reversal
    (1, always +1 here)), so the nuclear structure factor, the reflection
    list, the atom multiplicities and the symmetry constraints use the same
    symmetry as the magnetic structure.

    Parameters
    ----------
    op_rotation, op_translation : numpy.ndarray
        (9, n_s) and (3, n_s) operations ``{R_s|t_s}``.
    centring_translations : numpy.ndarray
        (3, n_c) centring translations including zero.

    Returns
    -------
    numpy.ndarray
        (14, n_s * n_c) integer elements.
    """
    l_elems = []
    for i_c in range(centring_translations.shape[1]):
        for i_s in range(op_rotation.shape[1]):
            fractions = [
                Fraction(float(hh)).limit_denominator(24) % 1 for hh in
                op_translation[:, i_s] + centring_translations[:, i_c]]
            denominator = int(numpy.lcm.reduce(
                [fr.denominator for fr in fractions]))
            numerators = [int(fr * denominator) for fr in fractions]
            rotation = [int(round(hh)) for hh in op_rotation[:, i_s]]
            elem = numerators + [denominator] + rotation + [1]
            if elem not in l_elems:  # superspace operations differing only in x4
                l_elems.append(elem)
    return numpy.array(l_elems, dtype=int).transpose()


def is_magnetic_k_refined(dict_crystal) -> bool:
    """True when a parameter of the magnetic k-structure part of a crystal
    dict is refined (Fourier components, positions, occupancies, B_iso,
    form factor), so that cached F_M_perp of satellites must be recalculated.
    """
    for key in (
            "flags_atom_fourier_cos", "flags_atom_fourier_sin",
            "flags_atom_fract_xyz",
            "flags_atom_occupancy", "flags_atom_b_iso",
            "flags_mag_atom_lande_factor", "flags_mag_atom_kappa"):
        if key in dict_crystal and numpy.any(dict_crystal[key]):
            return True
    return False


def calc_fourier_modulation_by_dictionary(dict_crystal):
    """Complex amplitudes ``Mcos + i Msin`` (3, n_atom) of the magnetic atoms.

    The parameters are the rows of ``_atom_site_moment_Fourier``
    (``atom_fourier_cos``, ``atom_fourier_sin``, one value per row); the
    symmetry constraints (symmform) are applied by the matrices
    ``atom_fourier_constraint_cos``/``_sin`` (2 n_row, 2 n_row): the
    effective values are ``C @ [cos, sin]``. Row ``i`` is the component
    ``atom_fourier_row_axis[i]`` of atom ``atom_fourier_row_atom[i]``.
    """
    cos = dict_crystal["atom_fourier_cos"]
    sin = dict_crystal["atom_fourier_sin"]
    values = dict_crystal["atom_fourier_constraint"] @ numpy.concatenate([cos, sin])
    n_row = cos.size
    n_atom = len(dict_crystal["atom_fourier_label"])
    modulation = numpy.zeros((3, n_atom), dtype=complex)
    numpy.add.at(modulation, (dict_crystal["atom_fourier_row_axis"],
                              dict_crystal["atom_fourier_row_atom"]),
                 values[:n_row] + 1j * values[n_row:])
    return modulation


def _take_magnetic_k_atoms(dict_crystal):
    """Arguments of :func:`calc_f_m_k` per propagation vector of a crystal dict.

    Yields ``(i_k, arguments)`` where ``arguments`` is the dict of
    keyword arguments of :func:`calc_f_m_perp_k_at_index_hkl` except
    ``index_hkl`` and ``unit_cell_parameters``. Atom positions are reduced to
    the first cell [0, 1), as in the rest of the crystal dictionary; the
    Fourier component ``S(p)`` is calculated for that position.
    """
    atom_index = dict_crystal["atom_fourier_index"]
    mag_index = dict_crystal["mag_atom_fourier_index"]
    atom_fract_xyz = (
        (dict_crystal["atom_site_sc_fract"].reshape(3, 3, -1)
         * dict_crystal["atom_fract_xyz"][na, :, :]).sum(axis=1)
        + dict_crystal["atom_site_sc_b"])[:, atom_index]
    atom_occupancy = dict_crystal["atom_occupancy"][atom_index]
    atom_b_iso = dict_crystal["atom_b_iso"][atom_index]
    lande = dict_crystal["mag_atom_lande_factor"][mag_index]
    kappa = dict_crystal["mag_atom_kappa"][mag_index]
    j0 = dict_crystal["mag_atom_j0_parameters"][:, mag_index]
    j2 = dict_crystal["mag_atom_j2_parameters"][:, mag_index]
    centring = dict_crystal["magnetic_k_centring"]
    # Mcos + i Msin along a, b, c of every magnetic atom (symmform constraints)
    modulation = calc_fourier_modulation_by_dictionary(dict_crystal)
    k_all = dict_crystal["magnetic_k_vector"]
    k_of_atom = k_all[:, dict_crystal["atom_fourier_k_index"]]
    # S(p) = (Mcos + i Msin) exp(-2 pi i k.p) / 2 at the position p used below
    fourier = 0.5 * modulation * numpy.exp(
        -2j * numpy.pi * (k_of_atom * atom_fract_xyz).sum(axis=0))[na, :]
    op_conjugate = dict_crystal.get("magnetic_k_op_conjugate")

    for i_k, k_vector in enumerate(dict_crystal["magnetic_k_vector"].transpose()):
        flag = dict_crystal["atom_fourier_k_index"] == i_k
        if not numpy.any(flag):
            continue

        def form_factor_function(sthovl, flag=flag):
            return calc_form_factor(
                sthovl[:, na], lande[na, flag], kappa[na, flag],
                j0[:, na, flag], j2[:, na, flag])[0]

        yield i_k, dict(
            k_vector=k_vector,
            op_rotation=dict_crystal["magnetic_k_op_rotation"],
            op_translation=dict_crystal["magnetic_k_op_translation"],
            op_moment=dict_crystal["magnetic_k_op_moment"],
            op_phase=dict_crystal["magnetic_k_op_phase"],
            centring_translations=centring,
            atom_fract_xyz=atom_fract_xyz[:, flag],
            atom_occupancy=atom_occupancy[flag],
            atom_b_iso=atom_b_iso[flag],
            form_factor_function=form_factor_function,
            atom_fourier_crystalaxis=fourier[:, flag],
            op_conjugate=op_conjugate)


def calc_f_m_perp_k_by_dictionary(
        dict_crystal, index_hkl, unit_cell_parameters):
    """F_M_perp of all propagation vectors of a crystal dict at given indices.

    Used for single crystals, where the measured reflections are given by
    fractional indices H. The contributions of all propagation vectors
    scattering at the same H are added coherently.

    Parameters
    ----------
    dict_crystal : dict
        Crystal dictionary with the ``magnetic_k_*`` and ``atom_fourier_*``
        keys (see ``form_magnetic_k_dictionary``).
    index_hkl : numpy.ndarray
        (3, n_hkl) float Miller indices in the reciprocal basis of the cell.
    unit_cell_parameters : numpy.ndarray
        a, b, c (Å), alpha, beta, gamma (radians).

    Returns
    -------
    numpy.ndarray
        (3, n_hkl) complex F_M_perp (10^-12 cm), Cartesian x || a*, z || c;
        zero at indices which are not satellites.
    """
    f_m_perp = numpy.zeros(numpy.shape(index_hkl), dtype=complex)
    for _, arguments in _take_magnetic_k_atoms(dict_crystal):
        f_m_perp += calc_f_m_perp_k_at_index_hkl(
            index_hkl, unit_cell_parameters, **arguments)
    return f_m_perp


def calc_powder_iint_magnetic_k_by_dictionary(
        dict_crystal, sthovl_min: float, sthovl_max: float,
        unit_cell_parameters):
    """Magnetic satellites of all propagation vectors of a crystal dict.

    The dict keys are written by ``form_magnetic_k_dictionary`` (see
    ``cryspy.E_data_classes.cl_1_crystal_dictionary``). Atom positions are
    reduced to the first cell [0, 1), as in the rest of the crystal
    dictionary; the moment phase refers to that cell.

    Parameters
    ----------
    dict_crystal : dict
        Crystal dictionary with the ``magnetic_k_*`` and ``atom_fourier_*``
        keys.
    sthovl_min, sthovl_max : float
        Range of sin(theta)/lambda (Å⁻¹).
    unit_cell_parameters : numpy.ndarray
        a, b, c (Å), alpha, beta, gamma (radians).

    Returns
    -------
    index_hkl : numpy.ndarray
        (3, n_hkl) satellite indices of all propagation vectors.
    iint : numpy.ndarray
        (n_hkl,) |F_M_perp|² (10^-24 cm²), multiplicity 1.
    k_index : numpy.ndarray
        (n_hkl,) index of the propagation vector of each satellite.
    """
    l_index, l_iint, l_k = [], [], []
    for i_k, arguments in _take_magnetic_k_atoms(dict_crystal):
        index_hkl, iint, _ = calc_powder_iint_magnetic_k(
            sthovl_min, sthovl_max, unit_cell_parameters, **arguments)
        l_index.append(index_hkl)
        l_iint.append(iint)
        l_k.append(numpy.full(iint.shape, i_k, dtype=int))
    if len(l_index) == 0:
        return (numpy.zeros((3, 0), dtype=float), numpy.zeros((0,), dtype=float),
                numpy.zeros((0,), dtype=int))
    return (numpy.concatenate(l_index, axis=1), numpy.concatenate(l_iint),
            numpy.concatenate(l_k))


def calc_unique_satellite_index_hkl_by_dictionary(
        dict_crystal, sthovl_min: float, sthovl_max: float,
        unit_cell_parameters):
    """Satellites of all propagation vectors, one per symmetry orbit.

    The satellites ``h ± k`` in the range are reduced by the rotations of
    the group of k and by Friedel's law (``H`` and ``-H`` give the same
    unpolarized intensity), so that a simulated single-crystal measurement
    has no equivalent reflections.

    Parameters
    ----------
    dict_crystal : dict
        Crystal dictionary with the ``magnetic_k_*`` keys.
    sthovl_min, sthovl_max : float
        Range of sin(theta)/lambda (Å⁻¹).
    unit_cell_parameters : numpy.ndarray
        a, b, c (Å), alpha, beta, gamma (radians).

    Returns
    -------
    numpy.ndarray
        (3, n_hkl) float indices of the satellites.
    """
    rotation = dict_crystal["magnetic_k_op_rotation"].reshape(3, 3, -1)
    l_matrix = [sign * rotation[:, :, s].transpose()
                for s in range(rotation.shape[2]) for sign in (1., -1.)]
    centring = dict_crystal["magnetic_k_centring"]
    l_index, keys = [], set()
    for k_vector in dict_crystal["magnetic_k_vector"].transpose():
        index_hkl, _ = calc_satellite_index_hkl(
            sthovl_min, sthovl_max, unit_cell_parameters, k_vector, centring)
        for index in index_hkl.transpose():
            key = tuple(numpy.round(index, 4) + 0.)
            if key in keys:
                continue
            l_index.append(index)
            for matrix in l_matrix:
                keys.add(tuple(numpy.round(matrix @ index, 4) + 0.))
    if len(l_index) == 0:
        return numpy.zeros((3, 0), dtype=float)
    return numpy.stack(l_index, axis=1)


def calc_moment_modulation_by_dictionary(dict_crystal):
    """Modulation of the moments of all magnetic atoms of the cell.

    Every magnetic atom of ``atom_site_moment_Fourier`` is copied by the
    superspace operations (with the centrings) into the cell [0, 1); a copy
    at ``p`` gets the complex amplitude ``A``, so that its moment in the
    cell shifted by the lattice vector L is ``M(p + L) = Re(A exp(-2 pi i
    k.L))``. For the atom of atom_site_moment_Fourier itself, at r = p + L,
    this is the msCIF rule ``Mcos cos(2 pi k.r) + Msin sin(2 pi k.r)``
    written for the cell: ``A = (Mcos + i Msin) exp(-2 pi i k.p)``.

    Parameters
    ----------
    dict_crystal : dict
        Crystal dictionary with the ``magnetic_k_*`` and ``atom_fourier_*``
        keys.

    Returns
    -------
    dict
        ``label`` (list of str, the atom of atom_site_moment_Fourier),
        ``fract_xyz`` (3, n) fractional positions in [0, 1), ``amplitude``
        (3, n) complex A along the unit vectors of a, b, c (μB),
        ``k_vector`` (3,).
    """
    atom_index = dict_crystal["atom_fourier_index"]
    position = (
        (dict_crystal["atom_site_sc_fract"].reshape(3, 3, -1)
         * dict_crystal["atom_fract_xyz"][na, :, :]).sum(axis=1)
        + dict_crystal["atom_site_sc_b"])[:, atom_index]
    k_vector = dict_crystal["magnetic_k_vector"][:, 0]
    modulation = calc_fourier_modulation_by_dictionary(dict_crystal)
    fourier = 0.5 * modulation * numpy.exp(-2j * numpy.pi * (k_vector @ position))[na, :]
    rotation = dict_crystal["magnetic_k_op_rotation"]
    translation = dict_crystal["magnetic_k_op_translation"]
    moment = dict_crystal["magnetic_k_op_moment"]
    phase = dict_crystal["magnetic_k_op_phase"]
    conjugate = dict_crystal.get(
        "magnetic_k_op_conjugate", numpy.zeros(phase.shape, dtype=bool))
    l_label, l_xyz, l_amplitude = [], [], []
    for i_atom, label in enumerate(dict_crystal["atom_fourier_label"]):
        l_key = []
        for i_op in range(phase.size):
            xyz = rotation[:, i_op].reshape(3, 3) @ position[:, i_atom] + translation[:, i_op]
            xyz_cell = numpy.mod(xyz, 1.)
            xyz_cell[numpy.isclose(xyz_cell, 1., atol=1e-6)] = 0.
            key = tuple(numpy.round(xyz_cell, 5))
            if key in l_key:
                continue
            l_key.append(key)
            fourier_op = fourier[:, i_atom]
            if conjugate[i_op]:
                fourier_op = numpy.conjugate(fourier_op)
            fourier_new = moment[:, i_op].reshape(3, 3) @ fourier_op * numpy.exp(
                -2j * numpy.pi * phase[i_op])
            l_label.append(str(label))
            l_xyz.append(xyz_cell)
            # S' refers to the unreduced position xyz: shift to the cell [0, 1)
            l_amplitude.append(2. * fourier_new * numpy.exp(
                2j * numpy.pi * (k_vector @ (xyz - xyz_cell))))
    return {"label": l_label, "fract_xyz": numpy.stack(l_xyz, axis=1),
            "amplitude": numpy.stack(l_amplitude, axis=1), "k_vector": k_vector}


def warn_not_satellites(index_hkl, k_vectors, tolerance: float = 1e-4):
    """Warn about measured fractional indices that are no satellites.

    A satellite ``h ± k`` given with rounded indices (0.282 for k = 0.2817)
    is not recognised by :func:`calc_f_m_perp_k_at_index_hkl` and gets a
    zero magnetic structure factor; call this for the measured indices
    (not for indices of other domains, which may legitimately belong to
    other arms of k).

    Parameters
    ----------
    index_hkl : numpy.ndarray
        (3, n_hkl) measured Miller indices.
    k_vectors : numpy.ndarray
        (3, n_k) propagation vectors (reciprocal lattice units).
    tolerance : float
        Tolerance on the indices.
    """
    index_hkl = numpy.asarray(index_hkl, dtype=float)
    flag_known = numpy.all(numpy.abs(index_hkl - numpy.round(index_hkl)) < tolerance, axis=0)
    for k_vector in numpy.asarray(k_vectors, dtype=float).transpose():
        for sign in (1., -1.):
            index_h = index_hkl - sign * k_vector[:, na]
            flag_known |= numpy.all(numpy.abs(index_h - numpy.round(index_h)) < tolerance, axis=0)
    if numpy.all(flag_known):
        return
    example = ", ".join(f"{value:g}" for value in index_hkl[:, ~flag_known][:, 0])
    warnings.warn(
        f"{int((~flag_known).sum())} reflections with fractional indices (e.g. "
        f"({example})) are not h ± k of the propagation vector within "
        f"{tolerance:g}: their magnetic structure factor is zero. Give the "
        "indices with enough digits or check k.", UserWarning, stacklevel=2)
