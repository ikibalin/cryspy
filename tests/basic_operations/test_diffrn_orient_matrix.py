"""UB matrix of DiffrnOrientMatrix from U and the cell in the notations of
the instruments."""
import numpy
import pytest

from cryspy.C_item_loop_classes.cl_1_cell import Cell
from cryspy.C_item_loop_classes.cl_2_diffrn_orient_matrix import DiffrnOrientMatrix

ABC_ANGLES = (4.26, 25.16, 4.28, 90., 97., 90.)


def calc_rotation(angle: float) -> numpy.ndarray:
    cos_angle, sin_angle = numpy.cos(angle), numpy.sin(angle)
    return numpy.array([[cos_angle, -sin_angle, 0.], [sin_angle, cos_angle, 0.], [0., 0., 1.]])


def test_calc_ub_ccsl_is_b():
    """U = 1 in CCSL gives UB = B of the cell (used to fail with
    AttributeError: get_notation)."""
    orient_matrix = DiffrnOrientMatrix.calc_ub(numpy.eye(3), "CCSL", ABC_ANGLES)
    cell = Cell(length_a=4.26, length_b=25.16, length_c=4.28, angle_alpha=90.,
                angle_beta=97., angle_gamma=90.)
    assert numpy.allclose(orient_matrix.ub, cell.m_b)
    assert numpy.allclose(orient_matrix.u, numpy.eye(3), atol=1e-10)


@pytest.mark.parametrize("notation", ["CCSL", "6T2@LLB", "5C1@LLB", "BusingLevy", "-YXZ"])
def test_calc_ub_notations(notation):
    """UB in a notation, taken to CCSL, is U_ccsl B with U_ccsl = T U, and the
    U of the object is U_ccsl; the notation and back give the same matrix."""
    u = calc_rotation(0.3)
    orient_matrix = DiffrnOrientMatrix.calc_ub(u, notation, ABC_ANGLES)
    transform = orient_matrix.take_transform_to_ccsl()
    cell = Cell(length_a=4.26, length_b=25.16, length_c=4.28, angle_alpha=90.,
                angle_beta=97., angle_gamma=90.)
    assert numpy.allclose(orient_matrix.calc_ub_ccsl(), transform @ u @ cell.m_b)
    assert numpy.allclose(orient_matrix.u, transform @ u, atol=1e-10)
    matrix = numpy.arange(9.).reshape(3, 3)
    assert numpy.allclose(
        orient_matrix.take_transform_to_ccsl() @ orient_matrix.calc_matrix_from_ccsl(matrix),
        matrix)


def test_busing_levy_axes():
    """CCSL rows (y, -x, z) of Busing-Levy (calc_ub_ccsl): the CCSL x axis
    is the Busing-Levy y axis."""
    orient_matrix = DiffrnOrientMatrix(type="BusingLevy")
    assert numpy.allclose(orient_matrix.calc_matrix_from_ccsl(numpy.eye(3))[:, 0], [0., 1., 0.])
