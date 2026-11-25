"""Powder lines of crystals with stacking faults (Hendricks-Teller model).

The crystal is built of layers normal to the reciprocal vector
N = (h_n k_n l_n) = g P (P primitive, g layers per cell along the normal):
the layer is the slab c_0 <= N.r < c_0 + 1 of the atoms of the cell (the
cut c_0 is put in the widest gap between the atoms when it is not given).
In the perfect crystal the next layer is placed by a symmetry operation
S_0 = (R_0, t_0) of the crystal with R_0^T N = N and N.t_0 = 1 (a
translation, a screw axis or a glide plane); at a fault it is placed by
S_v = (R_0, t_0 + t_v) instead, t_v being the fault vector t (in the plane
of the layers, N.t = 0, coordinates of the cell of the layer below) or one
of its images by the symmetry of the layer (the variants are equally
probable). The operations form a Markov chain: a fault follows a normal
contact with the probability alpha and a fault with the probability beta
(beta = alpha: independent, random faults; beta > alpha: the faults come in
runs, domains of the other stacking; beta < alpha: they avoid each other).

Along every rod Q = H + tau P (H integer, tau real) the intensity per cell
is (S. B. Hendricks and E. Teller, J. Chem. Phys. 10 (1942) 147;
M. M. J. Treacy, J. M. Newsam and M. W. Deem, Proc. R. Soc. Lond. A 433
(1991) 499, DIFFaX)

    I(Q) = g / n_G sum_k { |F_k|^2 + 2 Re[ F_k^* v . T (1 - T)^-1 psi_k ] },

F_k = F_L(R_k^T Q) the structure factor of the layer in its k-th
orientation (R_k = R_0^k, n_G orientations), T the transfer matrix of the
states (orientation, last operation) with the phases
exp(2 pi i (R_k^T Q).t_j), psi_k the stationary probabilities of the
operations for a reference layer in the orientation k and v the structure
factors of the states. The structure factor of the layer changes along the
rod, so the lines are asymmetric where |F_L| is. The reflections whose
phases do not depend on the operation stay Bragg peaks. In the limit
alpha -> 0 the rods become the Bragg reflections of the crystal. The faults
move intensity between the reflections of a rod (the layer structure factor
is not constant along it); the total over many reflections is kept.
"""
import itertools

import numpy

from cryspy.A_functions_base.unit_cell import calc_m_m_by_unit_cell_parameters

na = numpy.newaxis

PROBABILITY_MIN = 1e-6
PROBABILITY_MAX = 0.999
TOLERANCE = 1e-6


def calc_reciprocal_metric(unit_cell_parameters):
    """Metric tensor G* of the reciprocal lattice (1/angstrom^2) and B (H -> Cartesian)."""
    m_m = calc_m_m_by_unit_cell_parameters(unit_cell_parameters)[0].reshape(3, 3)
    m_b = numpy.linalg.inv(m_m).transpose()
    return m_b.transpose() @ m_b, m_b


def calc_layer_normal(normal_index, unit_cell_parameters):
    """Unit normal of the layers (Cartesian) and the distance between them 1/|N| (angstrom)."""
    m_b = calc_reciprocal_metric(unit_cell_parameters)[1]
    vector = m_b @ numpy.asarray(normal_index, dtype=float)
    length = numpy.linalg.norm(vector)
    if length == 0.:
        raise ValueError("The normal of the layers is (0 0 0).")
    return vector / length, 1. / length


def split_normal(normal_index):
    """Primitive normal P and the number of layers g per cell, N = g P."""
    normal_index = numpy.asarray(normal_index, dtype=int)
    g = int(numpy.gcd.reduce(numpy.abs(normal_index)))
    if g == 0:
        raise ValueError("The normal of the layers is (0 0 0).")
    return normal_index // g, g


def calc_unit_step(primitive):
    """Shortest integer vector e with P.e = 1 (lattice translation to the next set of layers)."""
    best = None
    for vector in itertools.product(range(-3, 4), repeat=3):
        if int(numpy.dot(primitive, vector)) == 1:
            norm = sum(abs(v) for v in vector)
            if best is None or norm < best[0]:
                best = (norm, numpy.array(vector, dtype=int))
    if best is None:
        raise ValueError(f"No lattice translation for the normal {tuple(primitive)}.")
    return best[1]


def calc_symmetry_operations(full_symm_elems, translation_elems):
    """All operations (R, t) of the space group: rotations (n, 3, 3), translations (n, 3)."""
    symm = numpy.asarray(full_symm_elems)
    rotation = symm[4:13].transpose().reshape(-1, 3, 3).astype(int)
    translation = (symm[0:3] / symm[3][na, :]).transpose().astype(float)
    centring = numpy.asarray(translation_elems)
    centring = (centring[0:3] / centring[3][na, :]).transpose().astype(float)
    l_rotation, l_translation = [], []
    for vector in centring:
        l_rotation.append(rotation)
        l_translation.append(translation + vector[na, :])
    return numpy.concatenate(l_rotation, axis=0), numpy.concatenate(l_translation, axis=0)


def calc_atoms_of_cell(atom_fract_xyz, rotations, translations):
    """Atoms of the cell (the copies of every atom of the asymmetric unit).

    Returns
    -------
    xyz : array (n, 3)
        Fractional coordinates, 0 <= x < 1.
    index_atom : array (n,) of int
        Atom of the asymmetric unit of every copy.
    index_operation : array (n,) of int
        Operation that made the copy (for the anisotropic displacements).
    """
    l_xyz, l_atom, l_operation = [], [], []
    xyz_atoms = numpy.asarray(atom_fract_xyz, dtype=float).reshape(3, -1).transpose()
    for i_atom, position in enumerate(xyz_atoms):
        copies = numpy.mod(numpy.einsum("nij,j->ni", rotations, position) + translations, 1.)
        l_kept = []
        for i_operation, copy in enumerate(copies):
            if any(numpy.all(numpy.abs(numpy.mod(copy - other + 0.5, 1.) - 0.5) < 1e-4)
                   for other in l_kept):
                continue
            l_kept.append(copy)
            l_xyz.append(numpy.mod(copy, 1.) % 1.)
            l_atom.append(i_atom)
            l_operation.append(i_operation)
    return (numpy.array(l_xyz).reshape(-1, 3), numpy.array(l_atom, dtype=int),
            numpy.array(l_operation, dtype=int))


def calc_layer_cut(xyz, normal_index):
    """Position c_0 of the cut in the widest gap of N.r (mod 1) between the atoms."""
    values = numpy.sort(numpy.mod(xyz @ numpy.asarray(normal_index, dtype=float), 1.))
    if values.size == 0:
        return 0.
    gaps = numpy.diff(numpy.concatenate([values, values[:1] + 1.]))
    i_gap = int(numpy.argmax(gaps))
    return float(numpy.mod(values[i_gap] + 0.5 * gaps[i_gap], 1.))


def select_layer(xyz, normal_index, layer_cut: float):
    """Atoms of the layer c_0 <= N.r < c_0 + 1, moved by lattice vectors into it.

    Returns
    -------
    xyz_layer : array (m, 3)
    index_copy : array (m,) of int
        Index of the copy of calc_atoms_of_cell.
    """
    normal_index = numpy.asarray(normal_index, dtype=int)
    primitive, g = split_normal(normal_index)
    step = calc_unit_step(primitive)
    value = xyz @ normal_index.astype(float)
    n_step = numpy.floor((value - layer_cut) / g + 1e-9)
    moved = xyz - n_step[:, na] * step[na, :]
    value = moved @ normal_index.astype(float)
    flag = (value - layer_cut) < 1. - 1e-9
    return moved[flag], numpy.flatnonzero(flag)


def calc_rotation_order(rotation) -> int:
    """Order n of R (R^n = 1)."""
    power = numpy.eye(3, dtype=int)
    for order in range(1, 13):
        power = power @ numpy.asarray(rotation, dtype=int)
        if numpy.array_equal(power, numpy.eye(3, dtype=int)):
            return order
    raise ValueError("The rotation has no finite order.")


def find_stacking_operation(normal_index, rotations, translations):
    """Symmetry operation S_0 = (R_0, t_0) that puts a layer onto the next one.

    R_0^T N = N and N.t_0 = 1 (a lattice translation along the normal is
    added to t_0). The operation with the rotation of the lowest order is
    taken (a translation if there is one); the others give the same next
    layer.

    Returns
    -------
    (rotation, translation) or None
        None when the layers are not related by a symmetry operation.
    """
    normal_index = numpy.asarray(normal_index, dtype=int)
    primitive, g = split_normal(normal_index)
    step = calc_unit_step(primitive)
    best = None
    for rotation, translation in zip(rotations, translations):
        if not numpy.array_equal(rotation.transpose() @ normal_index, normal_index):
            continue
        remainder = (1. - float(normal_index @ translation)) / g
        if abs(remainder - round(remainder)) > TOLERANCE:
            continue
        shift = translation + round(remainder) * step
        key = (calc_rotation_order(rotation), float(numpy.abs(shift).sum()))
        if best is None or key < best[0]:
            best = (key, rotation.astype(int), shift)
    if best is None:
        return None
    return best[1], best[2]


def calc_fault_variants(shift, normal_index, rotations, translations):
    """The fault vector and its images by the operations of the layer (different up to a lattice vector).

    The operations of the layer keep the normal (R^T N = N) and the layer
    (N.t = 0 up to a lattice vector).
    """
    normal_index = numpy.asarray(normal_index, dtype=int)
    g = split_normal(normal_index)[1]
    shift = numpy.asarray(shift, dtype=float)
    l_variant = []
    for rotation, translation in zip(rotations, translations):
        if not numpy.array_equal(rotation.transpose() @ normal_index, normal_index):
            continue
        remainder = float(normal_index @ translation) / g
        if abs(remainder - round(remainder)) > TOLERANCE:
            continue
        image = rotation @ shift
        if any(numpy.all(numpy.abs(image - other - numpy.round(image - other)) < TOLERANCE)
               for other in l_variant):
            continue
        l_variant.append(image)
    return numpy.array(l_variant).reshape(-1, 3)


def calc_orientations(rotation):
    """R_0^k, k = 0 ... n - 1 (the orientations of the layers)."""
    order = calc_rotation_order(rotation)
    l_power = [numpy.eye(3, dtype=int)]
    for _ in range(order - 1):
        l_power.append(l_power[-1] @ numpy.asarray(rotation, dtype=int))
    return numpy.array(l_power)


def calc_flag_faulted(index_hkl, orientations, variants):
    """True for the reflections that see the faults: (R_k^T H).t_v not an integer."""
    index_hkl = numpy.asarray(index_hkl, dtype=float).reshape(3, -1)
    flag = numpy.zeros(index_hkl.shape[1], dtype=bool)
    for rotation in orientations:
        rotated = rotation.transpose() @ index_hkl
        for variant in variants:
            product = variant @ rotated
            flag |= numpy.abs(product - numpy.round(product)) > TOLERANCE
    return flag


def calc_layer_structure_factor(q, layer: dict):
    """Structure factor of the layer at the points q (n, 3) of the reciprocal space (real indices).

    layer: "xyz" (m, 3), "amplitude" (m,) complex or "table_sthovl" (n_t,)
    with "table_amplitude" (m, n_t), "b_iso" (m,), "beta" (m, 3, 3) (the
    anisotropic displacements of every copy, exp(-h^T beta h)), "metric" G*.
    """
    q = numpy.atleast_2d(q)
    sthovl = 0.5 * numpy.sqrt(numpy.maximum(numpy.einsum("ni,ij,nj->n", q, layer["metric"], q), 0.))
    power = layer["b_iso"][na, :] * numpy.square(sthovl)[:, na] + numpy.einsum(
        "ni,mij,nj->nm", q, layer["beta"], q)
    if "table_amplitude" in layer:
        amplitude = numpy.stack([numpy.interp(sthovl, layer["table_sthovl"], row)
                                 for row in layer["table_amplitude"]], axis=1)
    else:
        amplitude = layer["amplitude"][na, :]
    phase = numpy.exp(2j * numpy.pi * (q @ layer["xyz"].transpose()))
    return (amplitude * numpy.exp(-power) * phase).sum(axis=1)


def calc_markov_probabilities(alpha: float, beta: float, n_variant: int):
    """Transition matrix p[j, j'] of the operations (0: normal, 1 ... n_v: faults) and the stationary pi."""
    alpha = float(numpy.clip(alpha, PROBABILITY_MIN, PROBABILITY_MAX))
    beta = float(numpy.clip(beta, 0., PROBABILITY_MAX))
    n_op = 1 + n_variant
    p = numpy.zeros((n_op, n_op))
    p[0, 0], p[0, 1:] = 1. - alpha, alpha / n_variant
    p[1:, 0], p[1:, 1:] = 1. - beta, beta / n_variant
    part_fault = alpha / (1. - beta + alpha)
    pi = numpy.concatenate([[1. - part_fault], numpy.full(n_variant, part_fault / n_variant)])
    return p, pi


def form_transfer_matrix(phase, p):
    """T[n, (k + 1, j'), (k, j)] = p[j, j'] exp(2 pi i (R_k^T Q).t_j')."""
    n_g, n_op, n_point = phase.shape
    n_state = n_g * n_op
    t_matrix = numpy.zeros((n_point, n_state, n_state), dtype=complex)
    for k in range(n_g):
        k_next = (k + 1) % n_g
        for j in range(n_op):
            for j_next in range(n_op):
                t_matrix[:, k_next * n_op + j_next, k * n_op + j] = p[j, j_next] * phase[k, j_next]
    return t_matrix


def calc_rod_intensity(f_layer, phase, alpha: float, beta: float, g: int):
    """Intensity per cell at points of the rods (Markov chain of the stacking operations).

    Parameters
    ----------
    f_layer : array (n_g, n) complex
        F_L(R_k^T Q) for the orientations k of the layers.
    phase : array (n_g, n_op, n) complex
        exp(2 pi i (R_k^T Q).t_j) of the operation j for a layer in the
        orientation k (j = 0 the normal stacking, then the faults).
    alpha, beta : float
        P(fault | normal contact) and P(fault | fault).
    g : int
        Number of layers per cell along the normal.

    Returns
    -------
    array (n,)
        In the units of |F_cell|^2 per unit of tau (Q = H + tau P).
    """
    n_g, n_op, n_point = phase.shape
    p, pi = calc_markov_probabilities(alpha, beta, n_op - 1)
    n_state = n_g * n_op
    t_matrix = form_transfer_matrix(phase, p)
    psi = numpy.zeros((n_point, n_state, n_g), dtype=complex)
    for k in range(n_g):
        psi[:, k * n_op:(k + 1) * n_op, k] = pi[na, :]
    solution = numpy.linalg.solve(numpy.eye(n_state)[na, :, :] - t_matrix, t_matrix @ psi)
    v = numpy.repeat(f_layer.transpose(), n_op, axis=1)  # (n, n_state)
    total = numpy.zeros(n_point)
    for k in range(n_g):
        corr = numpy.conj(f_layer[k]) * (v * solution[:, :, k]).sum(axis=1)
        total += numpy.square(numpy.abs(f_layer[k])) + 2. * corr.real
    return g * total / n_g


def calc_rod_width(phase, alpha: float, beta: float, g: int):
    """HWHM along the rods (units of tau) from the largest eigenvalue of the transfer matrix.

    Near a maximum of a rod the intensity is a Lorentzian of the HWHM
    -g ln|lambda| / (2 pi) in tau (Q = H + tau P).
    """
    p = calc_markov_probabilities(alpha, beta, phase.shape[1] - 1)[0]
    modulus = numpy.abs(numpy.linalg.eigvals(form_transfer_matrix(phase, p))).max(axis=1)
    return -g * numpy.log(numpy.clip(modulus, 1e-12, 1.)) / (2. * numpy.pi)
