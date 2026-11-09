"""sin(theta)/lambda and the volume of the unit cell come from one library,
unit_cell.py; the old functions of function_2_crystallography_base are
deprecated wrappers around it."""
import numpy
import pytest

import cryspy
from cryspy.A_functions_base.function_2_crystallography_base import (
    calc_inverse_d_by_hkl_abc_angles,
    calc_sthovl_by_hkl_abc_cosines,
    calc_volume_uc_by_abc_angles,
)

# triclinic cell: a, b, c in Å, angles in degrees
CELL = (5.1, 6.3, 7.7, 81.0, 97.5, 104.0)
INDEX_HKL = numpy.array([[1, 0, 0, 2, -1], [0, 1, 0, 1, 3], [0, 0, 1, -2, 2]])


def reference_sthovl_and_volume():
    """From the metric tensor G: 1/d^2 = h^T G^-1 h, V = sqrt(det G)."""
    a, b, c = CELL[:3]
    alpha, beta, gamma = numpy.radians(CELL[3:])
    m_g = numpy.array([
        [a*a, a*b*numpy.cos(gamma), a*c*numpy.cos(beta)],
        [a*b*numpy.cos(gamma), b*b, b*c*numpy.cos(alpha)],
        [a*c*numpy.cos(beta), b*c*numpy.cos(alpha), c*c]])
    inv_d_sq = numpy.einsum("in,ij,jn->n", INDEX_HKL, numpy.linalg.inv(m_g), INDEX_HKL)
    return 0.5*numpy.sqrt(inv_d_sq), numpy.sqrt(numpy.linalg.det(m_g))


def test_cell_calc_sthovl_against_metric_tensor():
    cell = cryspy.Cell(length_a=CELL[0], length_b=CELL[1], length_c=CELL[2],
                       angle_alpha=CELL[3], angle_beta=CELL[4], angle_gamma=CELL[5])
    sthovl_ref, volume_ref = reference_sthovl_and_volume()
    assert numpy.allclose(cell.calc_sthovl(*INDEX_HKL), sthovl_ref, rtol=1e-12)
    assert numpy.isclose(cell.calc_sthovl(1, 0, 0), sthovl_ref[0], rtol=1e-12)
    assert numpy.isclose(cell.volume, volume_ref, rtol=1e-10)


def test_old_functions_are_deprecated_wrappers():
    sthovl_ref, volume_ref = reference_sthovl_and_volume()
    a, b, c = CELL[:3]
    angles = numpy.radians(CELL[3:])
    with pytest.warns(DeprecationWarning):
        sthovl = calc_sthovl_by_hkl_abc_cosines(*INDEX_HKL, a, b, c, *numpy.cos(angles))
    assert numpy.allclose(sthovl, sthovl_ref, rtol=1e-12)
    with pytest.warns(DeprecationWarning):
        inverse_d = calc_inverse_d_by_hkl_abc_angles(*INDEX_HKL, a, b, c, *angles)
    assert numpy.allclose(inverse_d, 2*sthovl_ref, rtol=1e-12)
    with pytest.warns(DeprecationWarning):
        assert numpy.isclose(calc_volume_uc_by_abc_angles(a, b, c, *angles), volume_ref, rtol=1e-12)
