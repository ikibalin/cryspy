"""Powder diffraction of layers without correlations between them (Warren line shape).

A stack of layers ordered in the plane but not correlated from one layer
to the next scatters along rods (hk) perpendicular to the layers instead of
Bragg points. The powder average of a rod gives the asymmetric saw-tooth
line of B. E. Warren, Phys. Rev. 59 (1941) 693: a sharp front at the
2theta of the two-dimensional reflection and a long tail to high angles.

Notation (no 2 pi in the reciprocal space, s = 1/d = 2 sin(theta)/lambda):

- n: unit vector normal to the layers (Cartesian, x || a*, z || c), the
  direction of the reciprocal vector N = (h_n, k_n, l_n) of the layer
  planes; d_N the distance between the planes of N;
- a rod passes through the in-plane vector G (|G| = g) along n; s_z is
  the coordinate along the rod;
- the layers have a Gaussian in-plane correlation, the radial profile of
  a rod is p(v - g) = L exp(-pi L^2 (v - g)^2) (L, angstrom, is the mean
  layer dimension of Warren, the correlation length in the plane).

The powder-averaged intensity of one rod at s (the same normalization as
m |F|^2 delta(s - s_hkl) / s^2 of a Bragg reflection) is

    D(s) = 1/s [h(+u) + h(-u)] J(s, g, L),   u = sqrt(max(s^2 - g^2, 0)),
    J(s, g, L) = sqrt(4 L / (sqrt(pi) (s + g))) F(a),
    a = sqrt(pi) L (s - g),
    F(a) = int_0^inf exp(-(x^2 - a)^2) dx,

where h(s_z) is the intensity of the rod per unit s_z (d_N times the sum
of |F|^2 of the single layers of the unit cell). J is the integral
int_0^s p(v - g) / sqrt(s^2 - v^2) dv with sqrt(s + v) ~ sqrt(s + g):
for s = g it is Warren's formula, far in the tail it gives the exact
1 / sqrt(s^2 - g^2) of an infinitely sharp rod. The intensity of a rod
along s_z is the sum of |F|^2 of the single layers of the cell: it keeps
the total intensity (the mean of |F_cell|^2 over the Bragg points of the
rod), so the scale factor has the meaning of the scale of the phase.
"""
import numpy

from cryspy.A_functions_base.unit_cell import calc_m_m_by_unit_cell_parameters

na = numpy.newaxis

# F(a) is tabulated on [A_MIN, A_TABLE_MAX]; below A_MIN it is below 1e-25,
# above A_TABLE_MAX the Gauss-Hermite rule is used
A_MIN = -8.
A_TABLE_MAX = 30.
A_TABLE_STEP = 0.002
_WARREN_TABLE = {}


def _calc_warren_function_direct(a):
    """F(a) by Gauss-Legendre quadrature over x (accurate for a <= 40)."""
    a = numpy.asarray(a, dtype=float)
    nodes, weights = numpy.polynomial.legendre.leggauss(600)
    x_max = numpy.sqrt(numpy.maximum(a, 0.) + 7.)
    x = 0.5 * (nodes[na, :] + 1.) * x_max[:, na]
    integrand = numpy.exp(-numpy.square(numpy.square(x) - a[:, na]))
    return 0.5 * x_max * (integrand * weights[na, :]).sum(axis=1)


def _calc_warren_function_large(a):
    """F(a) = int exp(-y^2) / (2 sqrt(y + a)) dy by Gauss-Hermite (a >= 20)."""
    nodes, weights = numpy.polynomial.hermite.hermgauss(40)
    return (weights[na, :] / (2. * numpy.sqrt(
        nodes[na, :] + numpy.asarray(a, dtype=float)[:, na]))).sum(axis=1)


def calc_warren_function(a):
    """Warren function F(a) = int_0^inf exp(-(x^2 - a)^2) dx.

    F(0) = Gamma(5/4) = 0.9064; F(a) -> 0 fast for a < 0 (the front of
    the line), F(a) ~ sqrt(pi) / (2 sqrt(a)) for large a (the tail).
    Reference: B. E. Warren, Phys. Rev. 59 (1941) 693.

    Parameters
    ----------
    a : array
        Argument (dimensionless).

    Returns
    -------
    array of the shape of a
    """
    a = numpy.asarray(a, dtype=float)
    if "a" not in _WARREN_TABLE:
        grid = numpy.arange(A_MIN, A_TABLE_MAX + A_TABLE_STEP, A_TABLE_STEP)
        _WARREN_TABLE["a"] = grid
        _WARREN_TABLE["f"] = _calc_warren_function_direct(grid)
    res = numpy.interp(a, _WARREN_TABLE["a"], _WARREN_TABLE["f"], left=0.)
    large = a > A_TABLE_MAX
    if numpy.any(large):
        res[large] = _calc_warren_function_large(a[large])
    return res


def calc_rod_powder_factor(s, g, length):
    """Powder average J(s, g, L) of a rod with Gaussian in-plane correlations.

    Parameters
    ----------
    s : array
        Length of the scattering vector, 1/angstrom (s = 1/d).
    g : array
        Distance of the rod from the origin, 1/angstrom; broadcast with s.
    length : float
        Correlation length L in the plane (Warren's mean layer dimension),
        angstrom.

    Returns
    -------
    array
        J = sqrt(4 L / (sqrt(pi) (s + g))) F(sqrt(pi) L (s - g)),
        angstrom.
    """
    length = abs(float(length))
    a = numpy.sqrt(numpy.pi) * length * (s - g)
    return numpy.sqrt(4. * length / (numpy.sqrt(numpy.pi) * (s + g))) * \
        calc_warren_function(a)


def calc_layer_normal(normal_index, unit_cell_parameters):
    """Unit normal of the layer planes and their distance.

    Parameters
    ----------
    normal_index : array (3,)
        Miller indices (h, k, l) of the layer planes (integers).
    unit_cell_parameters : array (6,)
        a, b, c (angstrom), alpha, beta, gamma (radians).

    Returns
    -------
    normal_ccs : array (3,)
        Unit vector along the reciprocal vector N, Cartesian (x || a*, z || c).
    d_n : float
        Distance between the planes of N, angstrom.
    m_m : array (3, 3)
        Cell vectors a, b, c (columns), Cartesian.
    """
    m_m = calc_m_m_by_unit_cell_parameters(unit_cell_parameters)[0].reshape(3, 3)
    normal_index = numpy.asarray(normal_index, dtype=float)
    vector = numpy.linalg.inv(m_m).transpose() @ normal_index
    length = numpy.linalg.norm(vector)
    if length == 0.:
        raise ValueError("The normal of the layers is (0 0 0).")
    return vector / length, 1. / length, m_m


def calc_layer_cuts(heights, d_n: float, tolerance: float = 0.05):
    """Planes cutting a period of the stacking into layers at its widest gaps.

    Parameters
    ----------
    heights : array (n_atom,)
        Heights u = N.r of the atoms (fractional, period 1).
    d_n : float
        Distance between the planes of N, angstrom.
    tolerance : float
        Gaps narrower than the widest one by less than this (angstrom)
        are cut too (equivalent layers of the cell).

    Returns
    -------
    cuts : array (n_cut,)
        Heights of the cutting planes in [0, 1), in the middle of the gaps.
    gaps : array (n_atom,)
        Widths of all gaps between the atomic planes (angstrom), largest first.
    """
    u = numpy.sort(numpy.mod(numpy.asarray(heights, dtype=float), 1.))
    if u.size == 0:
        return numpy.zeros(0), numpy.zeros(0)
    u_next = numpy.concatenate([u[1:], u[:1] + 1.])
    gaps = (u_next - u) * d_n
    gap_max = gaps.max()
    flag_cut = gaps >= gap_max - tolerance
    cuts = numpy.mod(0.5 * (u + u_next)[flag_cut], 1.)
    return numpy.sort(cuts), numpy.sort(gaps)[::-1]


def calc_layer_index(heights, cuts):
    """Index of the layer of every height (counted over the whole crystal).

    Parameters
    ----------
    heights : array
        Heights u = N.r (any real values).
    cuts : array (n_cut,)
        Cutting planes in [0, 1) (calc_layer_cuts).

    Returns
    -------
    array of int of the shape of heights
    """
    heights = numpy.asarray(heights, dtype=float)
    n_cut = len(cuts)
    floor_u = numpy.floor(heights)
    return (n_cut * floor_u).astype(int) + numpy.searchsorted(
        cuts, heights - floor_u, side="right")


def calc_lattice_vector_across(normal_index):
    """Lattice vector T (integers) with N.T = 1 (N primitive)."""
    normal_index = numpy.asarray(normal_index, dtype=int)
    for size in range(1, 4):
        rng = range(-size, size + 1)
        l_t = [numpy.array((i, j, k)) for i in rng for j in rng for k in rng
               if i * normal_index[0] + j * normal_index[1] + k * normal_index[2] == 1]
        if l_t:
            return min(l_t, key=lambda t: numpy.abs(t).sum())
    raise ValueError(f"The normal {tuple(normal_index)} of the layers must be a "
                     "primitive Miller index (no common divisor).")


def calc_layers_of_atoms(fract_xyz, normal_index, unit_cell_parameters,
                         tolerance: float = 0.05):
    """Split the atoms of the unit cell into layers along N.

    Every atom gets the index of its layer in the cell (0 ... n_layer - 1)
    and is moved by a lattice translation so that the atoms of a layer
    have contiguous heights (needed for the structure factor of a layer
    at non-integer reciprocal coordinates).

    Parameters
    ----------
    fract_xyz : array (3, n_atom)
        Fractional coordinates of all atoms of the unit cell.
    normal_index : array (3,)
        Miller indices of the layer planes (primitive integers).
    unit_cell_parameters : array (6,)
        a, b, c (angstrom), alpha, beta, gamma (radians).
    tolerance : float
        See calc_layer_cuts, angstrom.

    Returns
    -------
    fract_layer : array (3, n_atom)
        Fractional coordinates with the atoms of each layer together.
    layer_index : array (n_atom,)
    n_layer : int
        Number of layers in the unit cell.
    """
    fract_xyz = numpy.asarray(fract_xyz, dtype=float)
    normal_index = numpy.asarray(normal_index, dtype=int)
    d_n = calc_layer_normal(normal_index, unit_cell_parameters)[1]
    heights = normal_index @ fract_xyz
    cuts = calc_layer_cuts(heights, d_n, tolerance)[0]
    layer = calc_layer_index(heights, cuts)
    n_layer = max(len(cuts), 1)
    # layer j of the crystal lies in the cell period floor(j / n_layer)
    shift = -numpy.floor_divide(layer, n_layer)
    t_vector = calc_lattice_vector_across(normal_index)
    fract_layer = fract_xyz + t_vector[:, na] * shift[na, :]
    return fract_layer, numpy.mod(layer, n_layer), n_layer


def calc_rods(normal_index, unit_cell_parameters, s_max: float):
    """Rods of the layers inside the sphere s <= s_max.

    The rods pass through the projections of the reciprocal lattice
    vectors on the plane of the layers; the rod through the origin
    (the 00l line of the stacking) is not included.

    Parameters
    ----------
    normal_index : array (3,)
        Miller indices of the layer planes.
    unit_cell_parameters : array (6,)
        a, b, c (angstrom), alpha, beta, gamma (radians).
    s_max : float
        Maximal 1/d, 1/angstrom.

    Returns
    -------
    index_rod : array (3, n_rod) of int
        A reciprocal lattice vector on every rod.
    g_rod : array (n_rod,)
        Distance of the rods from the origin, 1/angstrom.
    g_ccs : array (3, n_rod)
        In-plane vectors of the rods, Cartesian, 1/angstrom.
    """
    normal_ccs, d_n, m_m = calc_layer_normal(normal_index, unit_cell_parameters)
    inv_m_t = numpy.linalg.inv(m_m).transpose()
    h_max = numpy.sqrt(s_max ** 2 + 0.25 / d_n ** 2)
    limits = numpy.ceil(h_max * numpy.linalg.norm(m_m, axis=0)).astype(int)
    grid = numpy.meshgrid(*[numpy.arange(-lim, lim + 1) for lim in limits],
                          indexing="ij")
    index_hkl = numpy.stack([hh.ravel() for hh in grid], axis=0)
    h_ccs = inv_m_t @ index_hkl
    g_ccs = h_ccs - normal_ccs[:, na] * (normal_ccs @ h_ccs)[na, :]
    g_rod = numpy.linalg.norm(g_ccs, axis=0)
    flag = (g_rod > 1e-8) & (g_rod <= s_max)
    index_hkl, g_ccs, g_rod = index_hkl[:, flag], g_ccs[:, flag], g_rod[flag]
    key = numpy.round(g_ccs * 1e6).astype(numpy.int64)
    unique_index = numpy.unique(key, axis=1, return_index=True)[1]
    unique_index = numpy.sort(unique_index)
    return index_hkl[:, unique_index], g_rod[unique_index], g_ccs[:, unique_index]


def calc_rod_intensity(
        s, g_ccs, index_rod, s_z, normal_index, unit_cell_parameters,
        fract_layer, layer_index, n_layer, scattering_length=None,
        moment_ccs=None, magnetic_form_factor=None, b_iso=None):
    """Intensity of the rods per unit length along the normal, h(s_z).

    h(s_z) = d_N sum_layers |F_layer(Q)|^2 with Q = G + s_z n. The nuclear
    part uses the scattering lengths, the magnetic one the component of
    the magnetic structure factor perpendicular to Q (0.2695 x 10^-12 cm
    per Bohr magneton, as everywhere in CrysPy).

    Parameters
    ----------
    s : array (n_s,)
        1/d of the points (for the form factors and the Debye-Waller factor).
    g_ccs : array (3, n_rod)
        In-plane vectors of the rods (calc_rods), 1/angstrom.
    index_rod : array (3, n_rod)
        Reciprocal lattice vectors on the rods (calc_rods).
    s_z : array (n_s, n_rod)
        Coordinate along the rods (sign included), 1/angstrom.
    normal_index, unit_cell_parameters
        See calc_rods.
    fract_layer, layer_index, n_layer
        Atoms of the cell split into layers (calc_layers_of_atoms).
    scattering_length : array (n_atom,), complex, optional
        Scattering length times occupancy (10^-12 cm); None for no nuclear part.
    moment_ccs : array (3, n_atom), optional
        Ordered moments times occupancy, Cartesian, Bohr magnetons; None
        for no magnetic part.
    magnetic_form_factor : array (n_s, n_atom), optional
        Magnetic form factors at s.
    b_iso : array (n_atom,), optional
        Isotropic displacement parameters B, angstrom^2.

    Returns
    -------
    array (n_s, n_rod)
        h(s_z), 10^-24 cm^2 angstrom.
    """
    normal_ccs, d_n, m_m = calc_layer_normal(normal_index, unit_cell_parameters)
    normal_index = numpy.asarray(normal_index, dtype=float)
    # Q in reciprocal fractional coordinates: H_rod + (s_z - H_rod.n) d_N N
    h_rod_z = normal_ccs @ (numpy.linalg.inv(m_m).transpose() @ index_rod)
    tau = (s_z - h_rod_z[na, :]) * d_n
    q_fract = index_rod[:, na, :] + tau[na, :, :] * normal_index[:, na, na]
    # (n_s, n_rod, n_atom)
    phase = numpy.exp(2j * numpy.pi * numpy.einsum("isr,ia->sra", q_fract, fract_layer))
    n_atom = fract_layer.shape[1]
    weight = numpy.ones((s.size, n_atom), dtype=float)
    if b_iso is not None:
        weight = weight * numpy.exp(-numpy.asarray(b_iso)[na, :] * numpy.square(0.5 * s)[:, na])
    res = numpy.zeros(s_z.shape, dtype=float)
    q_ccs = g_ccs[:, na, :] + normal_ccs[:, na, na] * s_z[na, :, :]
    q_norm = numpy.linalg.norm(q_ccs, axis=0)
    q_unit = q_ccs / numpy.where(q_norm == 0., 1., q_norm)[na, :, :]
    for i_layer in range(n_layer):
        flag = layer_index == i_layer
        if not numpy.any(flag):
            continue
        phase_l = phase[:, :, flag]
        if scattering_length is not None:
            f_n = (phase_l * (weight[:, flag] * numpy.asarray(scattering_length)[na, flag])[:, na, :]
                   ).sum(axis=2)
            res += numpy.square(numpy.abs(f_n))
        if moment_ccs is not None:
            ff = weight[:, flag] * magnetic_form_factor[:, flag]
            # (3, n_s, n_rod)
            f_m = 0.2695 * numpy.einsum("sra,sa,ia->isr", phase_l, ff, moment_ccs[:, flag])
            f_m_perp = f_m - q_unit * (q_unit * f_m).sum(axis=0)[na, :, :]
            res += numpy.square(numpy.abs(f_m_perp)).sum(axis=0)
    return d_n * res


def calc_layer_powder_density(s, g_rod, h_plus, h_minus, length):
    """Powder-averaged density sum_rods D(s) of the rods.

    Parameters
    ----------
    s : array (n_s,)
        1/d, 1/angstrom.
    g_rod : array (n_rod,)
        Distances of the rods from the origin.
    h_plus, h_minus : array (n_s, n_rod)
        Intensity of the rods at s_z = +u and -u, u = sqrt(max(s^2 - g^2, 0)).
    length : float
        Correlation length in the plane, angstrom.

    Returns
    -------
    array (n_s,)
        sum_rods (h_plus + h_minus) J / s, normalized as
        m |F|^2 delta(s - s_hkl) / s^2 of a Bragg reflection.
    """
    j_rod = calc_rod_powder_factor(s[:, na], g_rod[na, :], length)
    return ((h_plus + h_minus) * j_rod).sum(axis=1) / s


def calc_layer_rod_s_z(s, g_rod):
    """u = sqrt(max(s^2 - g^2, 0)), the height on the rod at the sphere s."""
    return numpy.sqrt(numpy.maximum(numpy.square(s[:, na]) - numpy.square(g_rod[na, :]), 0.))
