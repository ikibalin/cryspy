import numpy

from cryspy.A_functions_base.magnetic_k_symmetry import calc_supercell_magnetic_space_group


def test_antiferromagnetic_chain_k_half():
    """Moments +m, -m at x = 0, 1/2 of the supercell 2a: the half
    translation of the supercell is an operation with time reversal."""
    rotation = numpy.eye(3).reshape(9, 1)
    translation = numpy.zeros((3, 1))
    fract_xyz = numpy.array([[0., 0.5], [0., 0.], [0., 0.]])
    moment = numpy.array([[0., 0.], [0., 0.], [1.5, -1.5]])
    group = calc_supercell_magnetic_space_group(
        rotation, translation, numpy.array([2, 1, 1]), fract_xyz,
        numpy.array([0, 0]), moment, numpy.array([8., 4., 4.]))
    assert group["theta"].tolist() == [1, -1]
    assert numpy.allclose(group["translation"][:, 1], [0.5, 0., 0.])
    assert group["representative"].tolist() == [True, False]


def test_ferromagnetic_cell_has_no_time_reversal():
    rotation = numpy.eye(3).reshape(9, 1)
    translation = numpy.zeros((3, 1))
    fract_xyz = numpy.array([[0., 0.5], [0., 0.], [0., 0.]])
    moment = numpy.array([[0., 0.], [0., 0.], [1.5, 1.5]])
    group = calc_supercell_magnetic_space_group(
        rotation, translation, numpy.array([2, 1, 1]), fract_xyz,
        numpy.array([0, 0]), moment, numpy.array([8., 4., 4.]))
    assert group["theta"].tolist() == [1, 1]
