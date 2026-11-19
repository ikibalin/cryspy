"""Shared real spherical harmonics (texture, multipoles)."""
import math

import numpy
import pytest

from cryspy.A_functions_base.spherical_harmonics import (
    calc_normalized_associated_legendre,
    calc_real_spherical_harmonics_by_lm,
)
from cryspy.A_functions_base.multipoles import realsphharm2df


def gsas2_plmpsi(l, m, psi_deg):
    """GSAS-II PLMPSI (texturesubs/plmpsi.for, aplms.for, qlmn.for): the
    normalized associated Legendre function as a Fourier series in psi."""
    def factln(k):
        return math.lgamma(k + 1)

    def qlmn(l_, mm, nn):
        m_, n_ = abs(mm), abs(nn)
        if n_ > m_:
            m_, n_ = n_, m_
        temp = 0.5 * (factln(l_ + n_) + factln(l_ + m_) + factln(l_ - m_) + factln(l_ - n_))
        total = 0.
        for i in range(0, l_ - n_ + 1):
            if l_ - m_ - i >= 0 and m_ + n_ + i >= 0:
                term = math.exp(temp - factln(i) - factln(l_ - m_ - i) - factln(l_ - n_ - i)
                                - factln(m_ + n_ + i))
                total += -term if i % 2 else term
        q = total / 2.**l_
        if (l_ - m_ - n_) % 2:
            q = -q
        if mm < 0 and (l_ + nn) % 2:
            q = -q
        if nn < 0 and (l_ + mm) % 2:
            q = -q
        return q

    def aplms(l_, m_, s_):
        a = qlmn(l_, m_, s_) * qlmn(l_, 0, s_) * math.sqrt((2. * l_ + 1.) / 2.)
        if abs(m_) % 2 == 0:
            ap = a if s_ == 0 else 2. * a
            return -ap if abs(m_) % 4 == 2 else ap
        ap = 0. if s_ == 0 else 2. * a
        return -ap if abs(m_) % 4 == 1 else ap

    psi = numpy.radians(psi_deg)
    if abs(m) % 2 == 0:
        return sum(aplms(l, m, s) * numpy.cos(s * psi) for s in range(0, l + 1, 2))
    return sum(aplms(l, m, s) * numpy.sin(s * psi) for s in range(2, l + 1, 2))


def test_legendre_as_gsas2():
    psi = numpy.array([5., 33., 71., 90., 127., 166.])
    p_lm = calc_normalized_associated_legendre(16, numpy.cos(numpy.radians(psi)))
    for l in range(0, 17, 2):
        for m in range(0, l + 1):
            assert numpy.allclose(p_lm[l, m], gsas2_plmpsi(l, m, psi), atol=1e-9), (l, m)


def test_real_spherical_harmonics_orthonormal():
    cos_theta, w_theta = numpy.polynomial.legendre.leggauss(40)
    phi = numpy.arange(80) * 2. * numpy.pi / 80
    theta, phi = numpy.meshgrid(numpy.arccos(cos_theta), phi, indexing="ij")
    weight = (w_theta[:, None] * numpy.full(phi.shape, 2. * numpy.pi / 80))
    index_lm = numpy.array([(l, m) for l in range(0, 13) for m in range(-l, l + 1)]).T
    y_lm = calc_real_spherical_harmonics_by_lm(index_lm, theta, phi)
    overlap = numpy.einsum("aij,bij,ij->ab", y_lm, y_lm, weight)
    assert numpy.allclose(overlap, numpy.eye(index_lm.shape[1]), atol=1e-10)


def test_real_spherical_harmonics_l1_are_xyz():
    theta, phi = numpy.array([0.3, 1.2, 2.5]), numpy.array([0.1, -2.0, 4.0])
    y = calc_real_spherical_harmonics_by_lm([[1, 1, 1], [1, -1, 0]], theta, phi)
    coeff = numpy.sqrt(3. / (4. * numpy.pi))
    xyz = numpy.stack([numpy.sin(theta) * numpy.cos(phi), numpy.sin(theta) * numpy.sin(phi),
                       numpy.cos(theta)])
    assert numpy.allclose(y, coeff * xyz)


def test_wrong_lm():
    with pytest.raises(ValueError):
        calc_real_spherical_harmonics_by_lm([[2], [3]], 0.1, 0.2)


# Implementation of realsphharm2df before it used the shared harmonics.
def old_realsphharm2df(n, m, theta, phi):
    """Real spherical harmonics taken from [International tables] normalized for density functions\n
    n is order or real spherical harmonics up to 4\n
    m is in range [-n:n]\n
    theta  -- polar angle in range [0:180]\n
    phi -- azimuthal angle in range [0:360]"""
    one = numpy.ones_like(theta)
    x = numpy.sin(theta) * numpy.cos(phi)
    y = numpy.sin(theta) * numpy.sin(phi)
    z = numpy.cos(theta)
    if (n, m) == (0, 0):
        res = one * (1.0 / (4.0 * numpy.pi))
    elif (n, m) == (1, 1):
        res = x * (1.0 / (numpy.pi))
    elif (n, m) == (1, -1):
        res = y * (1.0 / (numpy.pi))
    elif (n, m) == (1, 0):
        res = z * (1.0 / (numpy.pi))
    elif (n, m) == (2, 0):
        res = (3.0 * z**2 - 1.0) * ((3.0 * 3.0**0.5) / (8.0 * numpy.pi))
    elif (n, m) == (2, 1):
        res = (z * x) * (3.0 / 4.0)
    elif (n, m) == (2, -1):
        res = (z * y) * (3.0 / 4.0)
    elif (n, m) == (2, 2):
        res = 0.5 * (x**2 - y**2) * (3.0 / 4.0)
    elif (n, m) == (2, -2):
        res = (x * y) * (3.0 / 4.0)
    elif (n, m) == (3, 0):
        res = (5.0 * z**3 - 3.0 * z) * (10.0 / (13.0 * numpy.pi))
    elif (n, m) == (3, 1):
        res = x * (5.0 * z**2 - 1.0) * (0.32033)
    elif (n, m) == (3, -1):
        res = y * (5.0 * z**2 - 1.0) * (0.32033)
    elif (n, m) == (3, 2):
        res = z * (x**2 - y**2) * (1.0)
    elif (n, m) == (3, -2):
        res = 2 * x * y * z * (1.0)
    elif (n, m) == (3, 3):
        res = (x**3 - 3 * x * y**2) * (4.0 / (3.0 * numpy.pi))
    elif (n, m) == (3, -3):
        res = (-(y**3) + 3 * y * x**2) * (4.0 / (3.0 * numpy.pi))
    elif (n, m) == (4, 0):
        res = (35.0 * z**4 - 30.0 * z**2 + 3.0) * (0.06942)
    elif (n, m) == (4, 1):
        res = x * (7.0 * z**3 - 3.0 * z) * (735.0 / (512.0 * 7.0**0.5 + 196.0))
    elif (n, m) == (4, -1):
        res = y * (7.0 * z**3 - 3.0 * z) * (735.0 / (512.0 * 7.0**0.5 + 196.0))
    elif (n, m) == (4, 2):
        res = (
            (x**2 - y**2)
            * (7.0 * z**2 - 1.0)
            * ((105.0 * 7.0**0.5) / (4.0 * (136.0 + 28.0 * 7.0**0.5)))
        )
    elif (n, m) == (4, -2):
        res = (
            2
            * x
            * y
            * (7.0 * z**2 - 1.0)
            * ((105.0 * 7.0**0.5) / (4.0 * (136.0 + 28.0 * 7.0**0.5)))
        )
    elif (n, m) == (4, 3):
        res = z * (x**3 - 3 * x * y**2) * (5.0 / 4.0)
    elif (n, m) == (4, -3):
        res = z * (-(y**3) + 3 * y * x**2) * (5.0 / 4.0)
    elif (n, m) == (4, 4):
        res = (x**4 - 6.0 * x**2 * y**2 + y**4) * (15.0 / 32.0)
    elif (n, m) == (4, -4):
        res = (4.0 * x**3 * y - 4.0 * x * y**3) * (15.0 / 32.0)
    else:
        res = []
    return res


def test_multipole_harmonics_unchanged():
    rng = numpy.random.default_rng(7)
    theta = rng.uniform(0., numpy.pi, 30)
    phi = rng.uniform(0., 2. * numpy.pi, 30)
    for l in range(5):
        for m in range(-l, l + 1):
            assert numpy.allclose(realsphharm2df(l, m, theta, phi),
                                  old_realsphharm2df(l, m, theta, phi), rtol=1e-12, atol=1e-14)
    assert realsphharm2df(5, 0, theta, phi) == []
