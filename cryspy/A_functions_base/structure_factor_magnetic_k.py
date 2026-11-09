"""Magnetic structure factor of a structure with a propagation vector k.

The magnetic moment of atom ``j`` in the lattice cell ``L`` is written as a
Fourier series (FullProf convention, J. Rodriguez-Carvajal, Physica B 192
(1993) 55; FullProf manual, "Magnetic structures with propagation vectors")

    m_Lj = sum_k S_kj exp(-2 pi i k . R_L),

where the sum runs over ``k`` and ``-k`` (``S_-kj = S_kj^*``). When ``k`` is
equivalent to ``-k`` (``2k`` is a reciprocal lattice vector, for example
``k = (1/2, 0, 0)``) only one term is present and ``S_kj`` is real up to a
global phase, so its components are the moment components in μB.

The atom ``j`` is generated from an independent atom by the symmetry
operation ``{R_s|t_s}`` of the group of ``k``. Its Fourier component is

    S_kjs = M_s S_kj exp(-2 pi i phi_s),

where ``M_s`` acts on the moment components (FullProf ``MSYM u,v,w,phase``).
Lattice centring translations ``t_c`` of the nuclear cell are lattice
translations, so they bring the phase factor ``exp(-2 pi i k . t_c)``.

Magnetic reflections appear at ``H = h + k``, where ``h`` is a reciprocal
lattice vector of the nuclear cell. The magnetic structure factor is

    F_M(H) = p sum_j O_j f_j(Q) T_j sum_{s,c} S_kjs exp(-2 pi i k . t_c)
             exp(2 pi i H . (R_s r_j + t_s + t_c)),

with ``p = 0.2695`` (10^-12 cm per μB). Here the sum over the generated
positions is normalized as in the rest of CrysPy: every distinct position is
counted once (``O_j`` is the true occupancy).

All functions work with plain numpy arrays and know nothing about CIF
objects.

Coordinate systems: fractional coordinates for atoms, Miller indices (non
integer for satellites) for reflections, moment components along the unit
vectors of the crystal axes a, b, c (as ``_atom_site_moment.crystalaxis_x``
and FullProf ``Rx, Ry, Rz``), and the Cartesian system x || a*, z || c for
the output structure factor. Angles of ``unit_cell_parameters`` are in
radians.
"""
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


def parse_moment_uvw(moment_uvw: str):
    """Parse a moment operation such as ``"-u,v,-w"`` (FullProf MSYM).

    Parameters
    ----------
    moment_uvw : str
        Three comma separated expressions of ``u``, ``v``, ``w`` (moment
        components along a, b, c). ``mx, my, mz`` and ``x, y, z`` are
        accepted as synonyms of ``u, v, w``.

    Returns
    -------
    numpy.ndarray
        Matrix stored as 9 components acting on the moment components.
    """
    text = moment_uvw.lower()
    for old, new in (("mx", "u"), ("my", "v"), ("mz", "w"),
                     ("x", "u"), ("y", "v"), ("z", "w")):
        text = text.replace(old, new)
    matrix, translation = _parse_linear_expression(text, ("u", "v", "w"))
    if numpy.any(translation != 0.):
        raise ValueError(f"Moment operation '{moment_uvw}' has a constant term.")
    return matrix


def _parse_linear_expression(expression: str, names: tuple):
    parts = [part.strip() for part in expression.replace(" ", "").split(",")]
    if len(parts) != 3:
        raise ValueError(f"Operation '{expression}' must have three components.")
    matrix = numpy.zeros((3, 3), dtype=float)
    translation = numpy.zeros((3,), dtype=float)
    for i_row, part in enumerate(parts):
        terms = part.replace("-", "+-").split("+")
        for term in terms:
            if term in ("", "+"):
                continue
            sign = 1.
            if term.startswith("-"):
                sign, term = -1., term[1:]
            name = term[-1]
            if name in names:
                coeff = term[:-1].rstrip("*")
                value = float(Fraction(coeff)) if coeff else 1.
                matrix[i_row, names.index(name)] += sign * value
            else:
                translation[i_row] += sign * float(Fraction(term))
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
        atom_fourier_crystalaxis, atom_multiplicity=None):
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
    s_kjs = (moment[:, :, :, na] * atom_fourier_crystalaxis[na, :, na, :]).sum(axis=1)

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
        atom_fourier_crystalaxis):
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
        atom_fourier_crystalaxis)
    iint = numpy.square(numpy.abs(f_m_perp)).sum(axis=0)
    return index_hkl, iint, f_m_perp


def calc_f_m_perp_k_at_index_hkl(
        index_hkl, unit_cell_parameters, k_vector, op_rotation,
        op_translation, op_moment, op_phase, centring_translations,
        atom_fract_xyz, atom_occupancy, atom_b_iso, form_factor_function,
        atom_fourier_crystalaxis, tolerance: float = 1e-4):
    """F_M_perp of one propagation vector at arbitrary (fractional) indices.

    The moments are ``m_Lj = sum_{±k} S_kj exp(-2 pi i k R_L)`` with
    ``S_-k = S_k*`` (FullProf convention), so the term ``+k`` scatters at
    ``H = h + k`` and the term ``-k`` at ``H = h - k`` (h integer). For a
    general k only one of ``H - k`` and ``H + k`` can be integer, so every
    fractional H belongs to one term without ambiguity. When k is
    equivalent to -k (2k integer) the two terms are one and the same, the
    sum has a single term and ``S_kj`` is the real moment. Indices which
    are not satellites of k get zero.

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
    atom_fract_xyz, atom_occupancy, atom_b_iso :
        See :func:`calc_f_m_k`.
    form_factor_function : callable
        Takes sin(theta)/lambda (n,) (Å⁻¹) and returns (n, n_atom).
    atom_fourier_crystalaxis : numpy.ndarray
        (3, n_atom) complex Fourier components ``S_kj`` of the ``+k`` term
        (μB, crystal axes).
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
    f_m_perp = numpy.zeros(index_hkl.shape, dtype=complex)
    signs = (1,) if is_k_equivalent_to_minus_k(k_vector) else (1, -1)
    for sign in signs:
        index_h = index_hkl - sign * k_vector[:, na]
        flag = numpy.all(
            numpy.abs(index_h - numpy.round(index_h)) < tolerance, axis=0)
        if not numpy.any(flag):
            continue
        index_sign = index_hkl[:, flag]
        sthovl = calc_sthovl_by_unit_cell_parameters(
            index_sign, unit_cell_parameters)[0]
        if sign == 1:
            fourier, phase = atom_fourier_crystalaxis, op_phase
        else:
            fourier, phase = numpy.conjugate(atom_fourier_crystalaxis), -op_phase
        f_m_ccs = calc_f_m_k(
            index_sign, sign * k_vector,
            unit_cell_parameters, op_rotation, op_translation, op_moment,
            phase, centring_translations, atom_fract_xyz, atom_occupancy,
            atom_b_iso, form_factor_function(sthovl), fourier,
            atom_multiplicity=multiplicity)
        f_m_perp[:, flag] = calc_f_m_perp_k(
            index_sign, unit_cell_parameters, f_m_ccs)
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
            l_elems.append(numerators + [denominator] + rotation + [1])
    return numpy.array(l_elems, dtype=int).transpose()


def is_magnetic_k_refined(dict_crystal) -> bool:
    """True when a parameter of the magnetic k-structure part of a crystal
    dict is refined (Fourier components, positions, occupancies, B_iso,
    form factor), so that cached F_M_perp of satellites must be recalculated.
    """
    for key in (
            "flags_atom_fourier_real", "flags_atom_fourier_imag",
            "flags_atom_fourier_phase", "flags_atom_fract_xyz",
            "flags_atom_occupancy", "flags_atom_b_iso",
            "flags_mag_atom_lande_factor", "flags_mag_atom_kappa"):
        if key in dict_crystal and numpy.any(dict_crystal[key]):
            return True
    return False


def _take_magnetic_k_atoms(dict_crystal):
    """Arguments of :func:`calc_f_m_k` per propagation vector of a crystal dict.

    Yields ``(i_k, arguments)`` where ``arguments`` is the dict of
    keyword arguments of :func:`calc_f_m_perp_k_at_index_hkl` except
    ``index_hkl`` and ``unit_cell_parameters``. Atom positions are reduced to
    the first cell [0, 1), as in the rest of the crystal dictionary; the
    moment phase refers to that cell.
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
    fourier = (
        dict_crystal["atom_fourier_real"]
        + 1j * dict_crystal["atom_fourier_imag"]) * numpy.exp(
            -2j * numpy.pi * dict_crystal["atom_fourier_phase"])[na, :]

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
            atom_fourier_crystalaxis=fourier[:, flag])


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
