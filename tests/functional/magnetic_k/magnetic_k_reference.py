"""Independent real-space reference for magnetic structures with k.

It does not use CrysPy: it reads msCIF superspace operations by itself and
gives every atom of N1 x N2 x N3 cells its moment by the msCIF rule

    M_nu(h.x + epsilon x4 + tau) = theta det(R) R M_mu(x4),
    M_mu(x4) = Mcos cos(2 pi x4) + Msin sin(2 pi x4),   x4 = k.r,

checking that all operations giving the same site give the same moment.
"""
import itertools
import re
from fractions import Fraction

import numpy

MAGNETIC_SCATTERING_LENGTH = 0.2695  # 10^-12 cm per μB


def parse_operation(algebraic):
    """(4x4 matrix, 4 translation, theta) of 'x1+1/2,-x2,x3,-x4+1/4,+1'."""
    parts = [part.strip() for part in algebraic.split(",")]
    names = ("x1", "x2", "x3", "x4")
    matrix, translation = numpy.zeros((4, 4)), numpy.zeros(4)
    for i_row, part in enumerate(parts[:4]):
        expression = re.sub(r"(\d)(x)", r"\1*\2", part).replace("/", "*1./")

        def value(point, expression=expression):
            return eval(expression, {}, dict(zip(names, point)))

        translation[i_row] = value((0., 0., 0., 0.))
        for i_col in range(4):
            point = numpy.zeros(4)
            point[i_col] = 1.
            matrix[i_row, i_col] = value(point) - translation[i_row]
    return matrix, translation, int(parts[4])


def cartesian_basis(length, angle):
    """Columns a, b, c in the Cartesian system x || a*, z || c (Å)."""
    a, b, c = length
    alpha, beta, gamma = numpy.radians(angle)
    vector_c = numpy.array([0., 0., c])
    vector_b = numpy.array([0., b * numpy.sin(alpha), b * numpy.cos(alpha)])
    a_z = a * numpy.cos(beta)
    a_y = a * (numpy.cos(gamma) - numpy.cos(beta) * numpy.cos(alpha)) / numpy.sin(alpha)
    vector_a = numpy.array([numpy.sqrt(a * a - a_y * a_y - a_z * a_z), a_y, a_z])
    return numpy.stack([vector_a, vector_b, vector_c], axis=1)


def take_cells(k_vector):
    """Number of cells along a, b, c that make the crystal periodic."""
    return numpy.array([Fraction(value).limit_denominator(1000).denominator
                        for value in k_vector])


def generate_moments(operations, centrings, k_vector, atoms, n_cell,
                     atol: float = 2e-3):
    """Positions and moments of all atoms in N1 x N2 x N3 cells.

    Parameters
    ----------
    operations, centrings : list of str
        Algebraic superspace operations.
    k_vector : sequence of float
        (3,) propagation vector.
    atoms : dict
        {label: (position (3,), Fourier components (3, 2): cos and sin
        along a, b, c in μB)}.
    n_cell : numpy.ndarray
        (3,) number of cells.
    atol : float
        Tolerance of the check of equal moments at the same site (μB).

    Returns
    -------
    dict
        {label: list of (position in the large cell in units of the small
        cell, moment (3,) along a, b, c)}.
    """
    k_float = numpy.array([float(value) for value in k_vector])
    l_op = [parse_operation(text) for text in operations]
    l_cen = [parse_operation(text) for text in centrings] or [
        parse_operation("x1,x2,x3,x4,+1")]
    all_operations = [
        (m_c @ m_o, m_c @ t_o + t_c, theta_c * theta_o)
        for (m_c, t_c, theta_c) in l_cen for (m_o, t_o, theta_o) in l_op]
    lattice = numpy.array(list(itertools.product(*[range(n) for n in n_cell])))
    out = {}
    for label, (position, values) in atoms.items():
        sites = {}
        for matrix, translation, theta in all_operations:
            rotation, h_vector = matrix[:3, :3], matrix[3, :3]
            epsilon, tau = matrix[3, 3], translation[3]
            x = numpy.asarray(position, dtype=float)[numpy.newaxis, :] + lattice
            x_new = x @ rotation.transpose() + translation[:3]
            # x4 of the source atom giving the moment at the actual x4 = k.x_new
            x4 = epsilon * (x_new @ k_float - x @ h_vector - tau)
            angle = 2. * numpy.pi * x4
            values = numpy.asarray(values, dtype=float)
            moment = numpy.cos(angle)[:, numpy.newaxis] * values[:, 0] + \
                numpy.sin(angle)[:, numpy.newaxis] * values[:, 1]
            moment = theta * numpy.linalg.det(rotation) * moment @ rotation.transpose()
            x_new = numpy.mod(x_new, n_cell)
            x_new = numpy.where(numpy.abs(x_new - n_cell) < 1e-6, 0., x_new)
            for site, site_moment in zip(x_new, moment):
                key = tuple(numpy.round(site, 5))
                if key in sites:
                    assert numpy.allclose(sites[key][1], site_moment, atol=atol)
                else:
                    sites[key] = (site, site_moment)
        out[label] = list(sites.values())
    return out


def calc_f_m_perp(sites, n_cell, length, angle, index_hkl):
    """F_M_perp (3, n_hkl) per small cell by the direct sum, form factor 1,
    B = 0, in 10^-12 cm in the Cartesian system x || a*, z || c."""
    basis = cartesian_basis(length, angle)
    unit_vectors = basis / numpy.linalg.norm(basis, axis=0)[numpy.newaxis, :]
    l_site = [site for l_value in sites.values() for site in l_value]
    position = numpy.stack([site[0] for site in l_site], axis=1)
    moment = unit_vectors @ numpy.stack([site[1] for site in l_site], axis=1)
    f_m = MAGNETIC_SCATTERING_LENGTH * (
        moment[:, numpy.newaxis, :] *
        numpy.exp(2j * numpy.pi * (index_hkl.transpose() @ position))[numpy.newaxis]
        ).sum(axis=2) / numpy.prod(n_cell)
    q_vector = numpy.linalg.inv(basis).transpose() @ index_hkl
    e_q = q_vector / numpy.linalg.norm(q_vector, axis=0)
    return f_m - e_q * (e_q * f_m).sum(axis=0)
