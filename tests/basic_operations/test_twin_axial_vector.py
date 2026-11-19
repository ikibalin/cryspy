"""Laboratory axial vectors in the Cartesian system of twin domains."""
import numpy

from cryspy.A_functions_base.unit_cell import calc_axial_vector_by_twin_domains


def test_rotation_about_c():
    # tetragonal metric: 90 degrees about c, h_2 = (-k, h, l)
    unit_cell_parameters = numpy.array([4., 4., 6., 0.5 * numpy.pi, 0.5 * numpy.pi, 0.5 * numpy.pi])
    twin_matrices = numpy.array([0., -1., 0., 1., 0., 0., 0., 0., 1.])[:, numpy.newaxis]
    vector = numpy.array([1., 2., 3.])
    res = calc_axial_vector_by_twin_domains(vector, twin_matrices, unit_cell_parameters, 2)
    expected = numpy.array([[1., 1., -2., -2.], [2., 2., 1., 1.], [3., 3., 3., 3.]])
    assert numpy.allclose(res, expected)


def test_improper_law_and_length():
    # mirror l -> -l: an axial vector keeps its z component and flips x, y;
    # a law a <-> b of an orthorhombic cell (a != b) keeps the length
    unit_cell_parameters = numpy.array([5.36, 5.41, 13.1, 0.5 * numpy.pi, 0.5 * numpy.pi, 0.5 * numpy.pi])
    twin_matrices = numpy.array([
        [1., 0., 0., 0., 1., 0., 0., 0., -1.],
        [0., 1., 0., 1., 0., 0., 0., 0., -1.]]).transpose()
    vector = numpy.array([0.3, -0.2, 0.9])
    res = calc_axial_vector_by_twin_domains(vector, twin_matrices, unit_cell_parameters, 1)
    assert numpy.allclose(res[:, 1], [-0.3, 0.2, 0.9])
    assert numpy.isclose(numpy.linalg.norm(res[:, 2]), numpy.linalg.norm(vector))
    assert numpy.allclose(res[:, 2], [-0.2, 0.3, -0.9], atol=0.01)
