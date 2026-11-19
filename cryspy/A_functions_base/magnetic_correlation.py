"""Finite magnetic correlation length in powder diffraction.

The magnetic order of a powder grain can be correlated over a shorter
distance than the crystal structure (magnetic domains, short-range order).
For the exponential decay of the spin correlations exp(-r/xi) the magnetic
Bragg peaks are convolved with a Lorentzian of the half width at half
maximum

    kappa = 1 / xi        (in Q = 4 pi sin(theta) / lambda, inverse angstrom).

In 2theta (degrees) the full width at half maximum of this Lorentzian is

    FWHM = Y_mag / cos(theta),    Y_mag = 180 lambda / (pi^2 xi),

which is the Lorentzian size term Y/cos(theta) of the pseudo-Voigt profile
of FullProf (Thompson, Cox, Hastings, J. Appl. Cryst. 20 (1987) 79). The
apparent Scherrer size given by FullProf for this Y is
D = 360 lambda / (pi^2 Y) = 2 xi (FullProf manual, size broadening; the
integral breadth of the Lorentzian in 1/d is 1 / (2 xi)). So the option is
the "magnetic phase with its own size" of FullProf, but with the nuclear
and the magnetic intensities of the same phase.

Anisotropic correlation length
------------------------------
For an ellipsoid of correlation lengths the length along the scattering
vector H = h a* + k b* + l c* is

    xi(h)^2 = sum_ij h_i h_j a*_i a*_j Xi_ij / |H|^2,

with the symmetric tensor Xi_ij (angstrom^2) written in the convention of
the anisotropic displacement parameters U_ij (Trueblood et al., Acta Cryst.
A52 (1996) 770). Xi_11, Xi_22, Xi_33 are the squares of the correlation
lengths along a*, b*, c*. The isotropic length xi is the tensor
Xi_ij = xi^2 cos(angle between a*_i and a*_j). The tensor obeys the
symmetry of the crystal class as the U_ij of an atom at a site with the
full point symmetry.

Components of the tensors are stored as (11, 22, 33, 12, 13, 23) along
axis 0, as the U_ij and beta_ij in CrysPy.
"""
import numpy

from cryspy.A_functions_base.unit_cell import (
    calc_reciprocal_by_unit_cell_parameters,
    calc_m_b_by_unit_cell_parameters,
    calc_sthovl_by_unit_cell_parameters,
)

na = numpy.newaxis

# (i, j) of the components 11, 22, 33, 12, 13, 23
INDEX_IJ = ((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2))


def calc_y_mag_by_correlation_length(correlation_length, wavelength):
    """Lorentzian size term Y_mag (degrees) of the correlation length.

    The full width at half maximum of the magnetic Lorentzian in 2theta is
    Y_mag / cos(theta), Y_mag = 180 lambda / (pi^2 xi): the half width in
    Q is 1/xi. A length <= 0 means an infinite correlation length (no
    broadening), Y_mag = 0.

    Parameters
    ----------
    correlation_length : float or array
        Correlation length xi, angstrom.
    wavelength : float
        Neutron wavelength, angstrom.

    Returns
    -------
    float or array
        Y_mag, degrees.
    """
    xi = numpy.asarray(correlation_length, dtype=float)
    safe_xi = numpy.where(xi > 0., xi, 1.)
    return numpy.where(xi > 0., 180. * wavelength / (numpy.square(numpy.pi) * safe_xi), 0.)


def calc_aniso_by_isotropic_length(correlation_length, unit_cell_parameters):
    """Tensor Xi_ij (angstrom^2) of an isotropic correlation length.

    Xi_ij = xi^2 cos(angle between a*_i and a*_j), as U_ij = U_iso cos.

    Parameters
    ----------
    correlation_length : float
        Correlation length xi, angstrom.
    unit_cell_parameters : array (6,)
        a, b, c (angstrom), alpha, beta, gamma (radians).

    Returns
    -------
    array (6,)
        Xi_11, Xi_22, Xi_33, Xi_12, Xi_13, Xi_23.
    """
    reciprocal = calc_reciprocal_by_unit_cell_parameters(unit_cell_parameters)[0]
    c_ia, c_ib, c_ig = numpy.cos(reciprocal[3:6])
    xi_sq = float(correlation_length) ** 2
    return xi_sq * numpy.array([1., 1., 1., c_ig, c_ib, c_ia], dtype=float)


def form_b_by_aniso(aniso, unit_cell_parameters):
    """Quadratic form B (3, 3) in hkl: B_ij = a*_i a*_j Xi_ij."""
    reciprocal = calc_reciprocal_by_unit_cell_parameters(unit_cell_parameters)[0]
    a_star = numpy.asarray(reciprocal[:3], dtype=float)
    m_b = numpy.zeros((3, 3), dtype=float)
    for value, (i, j) in zip(aniso, INDEX_IJ):
        m_b[i, j] = value * a_star[i] * a_star[j]
        m_b[j, i] = m_b[i, j]
    return m_b


def form_aniso_by_b(m_b, unit_cell_parameters):
    """Tensor Xi_ij (6,) of the quadratic form B in hkl (inverse of form_b_by_aniso)."""
    reciprocal = calc_reciprocal_by_unit_cell_parameters(unit_cell_parameters)[0]
    a_star = numpy.asarray(reciprocal[:3], dtype=float)
    return numpy.array([m_b[i, j] / (a_star[i] * a_star[j]) for i, j in INDEX_IJ],
                       dtype=float)


def calc_correlation_length_hkl(index_hkl, unit_cell_parameters, aniso):
    """Correlation length along the scattering vector of the reflections.

    xi(h)^2 = sum_ij h_i h_j a*_i a*_j Xi_ij / |H|^2, |H| = 1/d.

    Parameters
    ----------
    index_hkl : array (3, n_hkl)
        Miller indices (fractional for the satellites).
    unit_cell_parameters : array (6,)
        a, b, c (angstrom), alpha, beta, gamma (radians).
    aniso : array (6,)
        Xi_11, Xi_22, Xi_33, Xi_12, Xi_13, Xi_23, angstrom^2.

    Returns
    -------
    array (n_hkl,)
        Correlation length, angstrom; 0 where the quadratic form is not
        positive (not a physical ellipsoid: no broadening there).
    """
    index_hkl = numpy.asarray(index_hkl, dtype=float)
    m_b = form_b_by_aniso(aniso, unit_cell_parameters)
    q_hkl = numpy.einsum("in,ij,jn->n", index_hkl, m_b, index_hkl)
    sthovl = calc_sthovl_by_unit_cell_parameters(index_hkl, unit_cell_parameters)[0]
    inv_d_sq = numpy.square(2. * sthovl)
    xi_sq = numpy.where(inv_d_sq > 0., q_hkl / numpy.where(inv_d_sq > 0., inv_d_sq, 1.), 0.)
    return numpy.sqrt(numpy.where(xi_sq > 0., xi_sq, 0.))


def symmetrize_aniso(aniso, rotations, unit_cell_parameters):
    """Average the tensor Xi_ij over the rotations of the crystal class.

    The quadratic form in hkl transforms as the beta_ij of the atomic
    displacements: B' = 1/N sum_R R B R^T for the rotations acting on hkl
    as h' = R^T h. The result is invariant under the class.

    Parameters
    ----------
    aniso : array (6,)
        Xi_ij, angstrom^2.
    rotations : array (n, 3, 3)
        Rotations of the Laue class acting on hkl as h' = R^T h.
    unit_cell_parameters : array (6,)
        a, b, c (angstrom), alpha, beta, gamma (radians).

    Returns
    -------
    array (6,)
        Symmetrized Xi_ij.
    """
    m_b = form_b_by_aniso(aniso, unit_cell_parameters)
    m_b_sym = numpy.einsum("nij,jk,nlk->il", rotations, m_b, rotations) / rotations.shape[0]
    return form_aniso_by_b(m_b_sym, unit_cell_parameters)


def calc_independent_aniso_components(rotations, unit_cell_parameters):
    """Independent components of Xi_ij under the crystal class.

    The components are taken in the order 11, 22, 33, 12, 13, 23; a
    component is independent when it is not fixed by the previous ones.

    Returns
    -------
    projection : array (6, 6)
        Symmetrization as a matrix acting on the vector Xi_ij.
    independent : array (6,) of bool
        True for the independent components.
    """
    projection = numpy.stack([
        symmetrize_aniso(numpy.eye(6)[i_comp], rotations, unit_cell_parameters)
        for i_comp in range(6)], axis=1)
    independent = numpy.zeros(6, dtype=bool)
    l_column = []
    rank = 0
    for i_comp in range(6):
        trial = numpy.stack(l_column + [projection[:, i_comp]], axis=1)
        if numpy.linalg.matrix_rank(trial, tol=1e-8) > rank:
            l_column.append(projection[:, i_comp])
            rank += 1
            independent[i_comp] = True
    return projection, independent


def calc_aniso_from_independent(aniso, projection, independent):
    """Tensor Xi_ij given by its independent components.

    The dependent components are found so that the symmetrized tensor keeps
    the values of the independent ones: Xi = P Xi, with Xi_ind fixed.
    """
    aniso = numpy.asarray(aniso, dtype=float)
    # the symmetric tensors are the image of P; take the one that has the
    # given independent components
    basis = projection[:, independent]
    coefficient = numpy.linalg.lstsq(basis[independent, :], aniso[independent],
                                     rcond=None)[0]
    return basis @ coefficient


def calc_principal_correlation_lengths(aniso, unit_cell_parameters):
    """Principal correlation lengths and axes of the ellipsoid.

    Parameters
    ----------
    aniso : array (6,)
        Xi_ij, angstrom^2.
    unit_cell_parameters : array (6,)
        a, b, c (angstrom), alpha, beta, gamma (radians).

    Returns
    -------
    lengths : array (3,)
        Principal correlation lengths, angstrom (0 for a non-positive
        eigenvalue), in increasing order.
    axes : array (3, 3)
        Unit vectors of the principal axes as columns, in the Cartesian
        system x || a*, z || c.
    """
    m_b = form_b_by_aniso(aniso, unit_cell_parameters)
    # H_cart = M h, so h^T B h = H^T (M^-T B M^-1) H
    m_m = numpy.asarray(calc_m_b_by_unit_cell_parameters(unit_cell_parameters)[0],
                        dtype=float).reshape(3, 3)
    m_m_inv = numpy.linalg.inv(m_m)
    xi_ccs = m_m_inv.transpose() @ m_b @ m_m_inv
    eigenvalues, axes = numpy.linalg.eigh(0.5 * (xi_ccs + xi_ccs.transpose()))
    lengths = numpy.sqrt(numpy.where(eigenvalues > 0., eigenvalues, 0.))
    return lengths, axes
