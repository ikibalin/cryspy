"""Constraint matrix of a symmetric tensor (susceptibility, exchange).

calc_constr_matrix turns the symmetry projector of a site into the matrix
that sets every component from the independent ones. The old version left
chi_33 = (chi_11 + chi_22)/2 on a cubic site, so the values and the sigmas
changed every time the constraints were applied.
"""
import numpy

from cryspy.C_item_loop_classes.cl_1_atom_site_susceptibility import (
    calc_constr_matrix,
    calc_independent_lines,
)


def projector_cubic():
    # chi_11 = chi_22 = chi_33 and chi_12 = chi_13 = chi_23 (site 3m along [111])
    block = numpy.full((3, 3), 1. / 3.)
    return numpy.block([[block, numpy.zeros((3, 3))], [numpy.zeros((3, 3)), block]])


def projector_tetragonal():
    # chi_11 = chi_22, chi_33 free, off-diagonal zero
    res = numpy.zeros((6, 6))
    res[:2, :2] = 0.5
    res[2, 2] = 1.
    return res


def test_cubic_site_sets_every_component_from_chi_11():
    matrix = calc_constr_matrix(projector_cubic())
    chi = numpy.array([1., 2., 5., 0.1, 0.2, 0.3])
    numpy.testing.assert_allclose(matrix @ chi, [1., 1., 1., 0.1, 0.1, 0.1])
    assert calc_independent_lines(matrix) == [0, 3]


def test_constraint_is_idempotent_for_values_and_sigmas():
    for projector in (projector_cubic(), projector_tetragonal()):
        matrix = calc_constr_matrix(projector)
        numpy.testing.assert_allclose(matrix @ matrix, matrix, atol=1e-12)
        # a tensor allowed by the symmetry is kept
        numpy.testing.assert_allclose(matrix @ projector, projector, atol=1e-12)
        sigma = numpy.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
        sigma_1 = numpy.sqrt(numpy.square(matrix) @ numpy.square(sigma))
        sigma_2 = numpy.sqrt(numpy.square(matrix) @ numpy.square(sigma_1))
        numpy.testing.assert_allclose(sigma_2, sigma_1)


def test_tetragonal_site():
    matrix = calc_constr_matrix(projector_tetragonal())
    chi = numpy.array([1., 2., 5., 0.1, 0.2, 0.3])
    numpy.testing.assert_allclose(matrix @ chi, [1., 1., 5., 0., 0., 0.])
    assert calc_independent_lines(matrix) == [0, 2]
