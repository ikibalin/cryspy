"""Real spherical harmonics shared by the multipole densities and the texture.

The functions work on numpy arrays of any shape and for any order l. Two
conventions are used in CrysPy:

- orthonormal real spherical harmonics ``Y_l^m`` (this module), with the
  normalized associated Legendre functions of Bunge without the
  Condon-Shortley phase;
- density normalized harmonics of the multipole model (``multipoles.py``),
  which are the same functions multiplied by constants.

References
----------
H.-J. Bunge, Texture Analysis in Materials Science, Butterworths (1982),
chapter 14 (normalized associated Legendre functions).
International Tables for Crystallography, Vol. B (2010), chapter 1.2
(real spherical harmonics of the multipole model).
"""
import numpy

INV_SQRT_PI = 1.0 / numpy.sqrt(numpy.pi)
INV_SQRT_2PI = 1.0 / numpy.sqrt(2.0 * numpy.pi)


def cartesian_to_spherical(x, y, z):
    """Spherical coordinates of Cartesian vectors.

    Parameters
    ----------
    x, y, z : array_like
        Cartesian components (any common shape).

    Returns
    -------
    r : numpy.ndarray
        Length of the vector.
    theta : numpy.ndarray
        Polar angle from z in radians, [0, pi]; 0 for the zero vector.
    phi : numpy.ndarray
        Azimuthal angle from x towards y in radians, (-pi, pi].
    """
    x = numpy.asarray(x, dtype=float)
    y = numpy.asarray(y, dtype=float)
    z = numpy.asarray(z, dtype=float)

    r = numpy.sqrt(x**2 + y**2 + z**2)

    theta = numpy.zeros_like(r)
    phi = numpy.zeros_like(r)

    mask = r > 0
    theta[mask] = numpy.arccos(numpy.clip(z[mask] / r[mask], -1.0, 1.0))
    phi[mask] = numpy.arctan2(y[mask], x[mask])
    return r, theta, phi


def calc_normalized_associated_legendre(l_max: int, cos_theta):
    """Normalized associated Legendre functions for 0 <= m <= l <= l_max.

    P_l^m(x) = sqrt((2l+1)/2 (l-m)!/(l+m)!) (1-x^2)^(m/2) d^m P_l(x)/dx^m,
    so that the integral of P_l^m(x)^2 over x in [-1, 1] is 1. The
    Condon-Shortley phase (-1)^m is not included (Bunge, 1982, eq. 14.1;
    GSAS-II ``pyplmpsi``). The stable recurrence over l at fixed m is used,
    so l_max may be large (texture analysis uses l up to 34).

    Parameters
    ----------
    l_max : int
        Maximal order l.
    cos_theta : array_like
        cos of the polar angle (any shape).

    Returns
    -------
    numpy.ndarray
        Array of the shape (l_max+1, l_max+1) + cos_theta.shape with
        P_l^m at index [l, m]; zero for m > l.
    """
    x = numpy.asarray(cos_theta, dtype=float)
    p_lm = numpy.zeros((l_max + 1, l_max + 1) + x.shape, dtype=float)
    for m, l, p_l_m in iterate_normalized_associated_legendre(
            {m: l_max for m in range(l_max + 1)}, x):
        p_lm[l, m] = p_l_m
    return p_lm


def iterate_normalized_associated_legendre(d_l_max: dict, cos_theta):
    """Normalized associated Legendre functions one by one.

    Yields (m, l, P_l^m) for every m of d_l_max and m <= l <= d_l_max[m]
    (see calc_normalized_associated_legendre). Only two orders are kept in
    memory at a time, which matters for large grids of directions.

    Parameters
    ----------
    d_l_max : dict
        Maximal order l for every index m >= 0 needed.
    cos_theta : array_like
        cos of the polar angle (any shape).
    """
    x = numpy.clip(numpy.asarray(cos_theta, dtype=float), -1.0, 1.0)
    sin_theta = numpy.sqrt(1.0 - numpy.square(x))
    p_mm = numpy.full(x.shape, numpy.sqrt(0.5))
    m_max = max(d_l_max) if len(d_l_max) > 0 else -1
    for m in range(m_max + 1):
        if m > 0:
            p_mm = numpy.sqrt((2.0 * m + 1.0) / (2.0 * m)) * sin_theta * p_mm
        if m not in d_l_max or d_l_max[m] < m:
            continue
        l_max = d_l_max[m]
        yield m, m, p_mm
        if m + 1 > l_max:
            continue
        p_l_2, p_l_1 = p_mm, numpy.sqrt(2.0 * m + 3.0) * x * p_mm
        yield m, m + 1, p_l_1
        for l in range(m + 2, l_max + 1):
            a_lm = numpy.sqrt((4.0 * l * l - 1.0) / (l * l - m * m))
            b_lm = numpy.sqrt(((l - 1.0)**2 - m * m) / (4.0 * (l - 1.0)**2 - 1.0))
            p_l_2, p_l_1 = p_l_1, a_lm * (x * p_l_1 - b_lm * p_l_2)
            yield m, l, p_l_1


def calc_real_spherical_harmonics_by_lm(index_lm, theta, phi):
    """Orthonormal real spherical harmonics Y_l^m for a list of (l, m).

    Y_l^0 = P_l^0(cos theta) / sqrt(2 pi),
    Y_l^m = P_l^m(cos theta) cos(m phi) / sqrt(pi) for m > 0,
    Y_l^m = P_l^|m|(cos theta) sin(|m| phi) / sqrt(pi) for m < 0,
    with the normalized associated Legendre functions P_l^m of
    ``calc_normalized_associated_legendre`` (no Condon-Shortley phase), so
    that the integral of Y_l^m Y_l'^m' over the unit sphere is
    delta_ll' delta_mm'. Y_1^1, Y_1^-1, Y_1^0 are proportional to x, y, z.

    Parameters
    ----------
    index_lm : array_like of int
        Orders and indexes, shape (2, n_lm): l >= 0, -l <= m <= l.
    theta : array_like
        Polar angle from z, radians.
    phi : array_like
        Azimuthal angle from x towards y, radians (same shape as theta, or
        broadcastable to it).

    Returns
    -------
    numpy.ndarray
        Y_l^m of the shape (n_lm,) + broadcast shape of theta and phi.
    """
    index_lm = numpy.asarray(index_lm, dtype=int).reshape(2, -1)
    theta, phi = numpy.broadcast_arrays(
        numpy.asarray(theta, dtype=float), numpy.asarray(phi, dtype=float))
    l_index, m_index = index_lm[0], index_lm[1]
    if numpy.any(l_index < 0) or numpy.any(numpy.abs(m_index) > l_index):
        raise ValueError("Spherical harmonics need l >= 0 and -l <= m <= l.")
    d_l_max = {}
    for l, m in zip(l_index, m_index):
        d_l_max[abs(int(m))] = max(int(l), d_l_max.get(abs(int(m)), 0))
    d_index = {}
    for i_lm, (l, m) in enumerate(zip(l_index, m_index)):
        d_index.setdefault((int(l), abs(int(m))), []).append((i_lm, int(m)))
    res = numpy.empty((l_index.size,) + theta.shape, dtype=float)
    for m_abs, l, p_l_m in iterate_normalized_associated_legendre(d_l_max, numpy.cos(theta)):
        for i_lm, m in d_index.get((l, m_abs), []):
            if m == 0:
                res[i_lm] = INV_SQRT_2PI * p_l_m
            elif m > 0:
                res[i_lm] = INV_SQRT_PI * p_l_m * numpy.cos(m * phi)
            else:
                res[i_lm] = INV_SQRT_PI * p_l_m * numpy.sin(-m * phi)
    return res


def calc_real_spherical_harmonic(l: int, m: int, theta, phi):
    """Orthonormal real spherical harmonic Y_l^m (see
    ``calc_real_spherical_harmonics_by_lm``).

    Parameters
    ----------
    l : int
        Order, l >= 0.
    m : int
        Index, -l <= m <= l.
    theta : array_like
        Polar angle from z, radians.
    phi : array_like
        Azimuthal angle from x towards y, radians.

    Returns
    -------
    numpy.ndarray
        Y_l^m of the broadcast shape of theta and phi.
    """
    return calc_real_spherical_harmonics_by_lm([[l], [m]], theta, phi)[0]
