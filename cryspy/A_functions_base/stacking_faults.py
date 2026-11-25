"""Powder lines of crystals with random stacking faults (displacement faults).

A crystal built of layers normal to the reciprocal vector N = (h_n k_n l_n)
(distance d_L = 1/|N| between the layers) has at every layer, with the
probability alpha, an extra in-plane shift t (fractional coordinates of the
cell) of all the layers above it. The reflections H with H.t integer do not
see the faults and stay Bragg reflections. The others are spread along the
normal n of the layers: the correlation of two layers m apart is
multiplied by beta^m,

    beta = 1 - alpha + alpha exp(2 pi i H.t),

and every Bragg point H becomes a Lorentzian along n (B. E. Warren,
X-ray Diffraction (1969), ch. 13; S. B. Hendricks and E. Teller,
J. Chem. Phys. 10 (1942) 147):

    P(z) = gamma / (pi ((z - z_H - delta)^2 + gamma^2)),
    gamma = -ln|beta| / (2 pi d_L),    delta = -arg(beta) / (2 pi d_L),

z being the coordinate along n (no 2 pi: s = 1/d units, 1/angstrom). For a
shift t with exp(2 pi i H.t) = -1 (e.g. t = b/2 and k odd) beta is real,
delta = 0 and gamma = -ln|1 - 2 alpha| / (2 pi d_L). The Lorentzian form is
the sum of the Lorentzians of all the Bragg points of a rod: it holds as
long as gamma is smaller than the distance between them.

The powder average keeps the intensity of every reflection: the part of
the Lorentzian at |Q| between s_1 and s_2 (u_i = sqrt(s_i^2 - g^2), g the
in-plane distance of H from the origin) is

    [C(u_2) - C(u_1)] + [C(-u_1) - C(-u_2)],
    C(z) = 1/2 + arctan((z - z_H - delta) / gamma) / pi,

exact for the 1/u singularity at the front of the line (the 110 line of a
rod normal to it gets a sharp front at s = g and a tail to high angles,
the hkl lines with l != 0 a broad, nearly symmetric peak).
"""
import numpy

from cryspy.A_functions_base.unit_cell import calc_m_m_by_unit_cell_parameters

na = numpy.newaxis

# |beta| is limited from below (alpha = 1/2 with exp(i phi) = -1 gives 0)
BETA_MIN = 1e-6


def calc_layer_normal(normal_index, unit_cell_parameters):
    """Unit normal of the layers and the distance between them.

    Parameters
    ----------
    normal_index : array (3,)
        Miller indices (h, k, l) of the layer planes; (0 0 2) are layers
        at c/2 (two layers per cell).
    unit_cell_parameters : array (6,)
        a, b, c (angstrom), alpha, beta, gamma (radians).

    Returns
    -------
    normal_ccs : array (3,)
        Unit vector along N = h a* + k b* + l c*, Cartesian (x || a*, z || c).
    d_layer : float
        Distance between the layers 1/|N|, angstrom.
    """
    m_m = calc_m_m_by_unit_cell_parameters(unit_cell_parameters)[0].reshape(3, 3)
    vector = numpy.linalg.inv(m_m).transpose() @ numpy.asarray(normal_index, dtype=float)
    length = numpy.linalg.norm(vector)
    if length == 0.:
        raise ValueError("The normal of the layers is (0 0 0).")
    return vector / length, 1. / length


def calc_fault_beta(index_hkl, shift, probability: float):
    """Factor beta = 1 - alpha + alpha exp(2 pi i H.t) of the reflections.

    Parameters
    ----------
    index_hkl : array (3, n_hkl)
        Miller indices.
    shift : array (3,)
        Fault vector t, fractional coordinates of the cell.
    probability : float
        Probability alpha of a fault between two neighbour layers (0 ... 1).

    Returns
    -------
    array (n_hkl,), complex
    """
    phase = 2. * numpy.pi * (numpy.asarray(shift, dtype=float) @ numpy.asarray(index_hkl, dtype=float))
    return 1. - probability + probability * numpy.exp(1j * phase)


def calc_flag_faulted(index_hkl, shift, tolerance: float = 1e-6):
    """True for the reflections that see the faults (H.t not an integer)."""
    product = numpy.asarray(shift, dtype=float) @ numpy.asarray(index_hkl, dtype=float)
    return numpy.abs(product - numpy.round(product)) > tolerance


def calc_rod_width(beta, d_layer: float):
    """HWHM gamma and centre shift delta of the Lorentzian along the normal.

    Parameters
    ----------
    beta : array, complex
        calc_fault_beta.
    d_layer : float
        Distance between the layers, angstrom.

    Returns
    -------
    gamma, delta : arrays of the shape of beta
        1/angstrom (s = 1/d units).
    """
    modulus = numpy.maximum(numpy.abs(beta), BETA_MIN)
    gamma = -numpy.log(numpy.minimum(modulus, 1.)) / (2. * numpy.pi * d_layer)
    delta = -numpy.angle(beta) / (2. * numpy.pi * d_layer)
    return gamma, delta


def calc_probability_by_width(gamma: float, d_layer: float, phase: float = numpy.pi):
    """Fault probability alpha (0 ... 1/2 side) giving the HWHM gamma.

    Inverse of calc_rod_width for one phase 2 pi H.t (radians):
    |1 - alpha + alpha exp(i phase)| = exp(-2 pi d_L gamma).
    """
    target = numpy.exp(-2. * numpy.pi * d_layer * gamma)
    grid = numpy.linspace(0., 0.5, 5001)
    modulus = numpy.abs(1. - grid + grid * numpy.exp(1j * phase))
    return float(grid[numpy.argmin(numpy.abs(modulus - target))])


def calc_rod_geometry(index_hkl, normal_index, unit_cell_parameters):
    """In-plane distance g and height z of the reflections on their rods.

    Parameters
    ----------
    index_hkl : array (3, n_hkl)
    normal_index, unit_cell_parameters
        See calc_layer_normal.

    Returns
    -------
    g : array (n_hkl,)
        Distance of the rod from the origin, 1/angstrom.
    z : array (n_hkl,)
        Coordinate of H along the normal, 1/angstrom.
    d_layer : float
        Distance between the layers, angstrom.
    """
    normal_ccs, d_layer = calc_layer_normal(normal_index, unit_cell_parameters)
    m_m = calc_m_m_by_unit_cell_parameters(unit_cell_parameters)[0].reshape(3, 3)
    h_ccs = numpy.linalg.inv(m_m).transpose() @ numpy.asarray(index_hkl, dtype=float)
    z = normal_ccs @ h_ccs
    s = numpy.linalg.norm(h_ccs, axis=0)
    g = numpy.sqrt(numpy.maximum(numpy.square(s) - numpy.square(z), 0.))
    return g, z, d_layer


def calc_lorentz_cdf(z, centre, gamma):
    """Cumulative distribution of a normalized Lorentzian."""
    return 0.5 + numpy.arctan((z - centre) / gamma) / numpy.pi


def calc_rod_mass(s_edges, g, z_centre, gamma):
    """Part of the intensity of every reflection between the edges in s.

    Parameters
    ----------
    s_edges : array (n_bin + 1,)
        Increasing edges of the bins, s = 2 sin(theta) / lambda, 1/angstrom.
    g : array (n_hkl,)
        In-plane distance of the rods, 1/angstrom.
    z_centre : array (n_hkl,)
        Centre of the Lorentzian along the normal (z_H + delta), 1/angstrom.
    gamma : array (n_hkl,)
        HWHM along the normal, 1/angstrom (> 0).

    Returns
    -------
    array (n_bin, n_hkl)
        The sum over the bins is the part of the reflection inside
        [s_edges[0], s_edges[-1]] (1 for a full range).
    """
    gamma = numpy.maximum(numpy.asarray(gamma, dtype=float), 1e-12)
    u = numpy.sqrt(numpy.maximum(numpy.square(s_edges[:, na]) - numpy.square(g[na, :]), 0.))
    cdf_plus = calc_lorentz_cdf(u, z_centre[na, :], gamma[na, :])
    cdf_minus = calc_lorentz_cdf(-u, z_centre[na, :], gamma[na, :])
    return (cdf_plus[1:] - cdf_plus[:-1]) + (cdf_minus[:-1] - cdf_minus[1:])
