import numpy
import pytest

from cryspy.A_functions_base.structure_factor_magnetic_k import (
    calc_f_m_k,
    calc_f_m_perp_k,
    calc_full_symm_elems_by_k_operations,
    calc_generated_multiplicity,
    calc_satellite_index_hkl,
    calc_superspace_operation_arrays,
    is_k_equivalent_to_minus_k,
    parse_operation_xyz,
    parse_ssg_operation,
)

UCP = numpy.array([4.0, 5.0, 6.0, 0.5 * numpy.pi, 0.5 * numpy.pi, 0.5 * numpy.pi])


def test_parse_operation_xyz():
    rotation, translation = parse_operation_xyz("-x+1/2, y, -z+3/4")
    numpy.testing.assert_array_equal(
        rotation, [-1, 0, 0, 0, 1, 0, 0, 0, -1])
    numpy.testing.assert_allclose(translation, [0.5, 0.0, 0.75])
    rotation, _ = parse_operation_xyz("x-y,x,z")
    numpy.testing.assert_array_equal(rotation, [1, -1, 0, 1, 0, 0, 0, 0, 1])


def test_parse_ssg_operation():
    operation = parse_ssg_operation("-x1+1/2,-x2,x3+1/2,-x1-x2+x4+1/4,-1")
    numpy.testing.assert_array_equal(
        operation["rotation"], [-1, 0, 0, 0, -1, 0, 0, 0, 1])
    numpy.testing.assert_allclose(operation["translation"], [0.5, 0., 0.5])
    numpy.testing.assert_allclose(operation["h"], [-1., -1., 0.])
    assert (operation["epsilon"], operation["tau"], operation["theta"]) == (1, 0.25, -1)
    operation = parse_ssg_operation("x2,x1,-x3,3x3-x4+1/2,+1")
    numpy.testing.assert_allclose(operation["h"], [0., 0., 3.])
    assert operation["epsilon"] == -1


def test_superspace_operation_arrays():
    """phi = k.t - tau, moment matrix theta det(R) R, centring x operation."""
    k_vector = numpy.array([0., 0.3, 0.])
    arrays = calc_superspace_operation_arrays(
        ["x1,x2,x3,x4,+1", "-x1,x2+1/2,-x3,x4+1/4,-1"],
        ["x1,x2,x3,x4,+1", "x1,x2,x3,x4+1/2,-1"], k_vector)
    numpy.testing.assert_allclose(arrays["phase"], [0., 0.15 - 0.25, -0.5, 0.15 - 0.75])
    numpy.testing.assert_array_equal(arrays["moment"][:, 1], [1, 0, 0, 0, -1, 0, 0, 0, 1])
    numpy.testing.assert_array_equal(arrays["conjugate"], False)
    with pytest.raises(ValueError):
        # R^T k = -k needs -x4
        calc_superspace_operation_arrays(["-x1,-x2,-x3,x4,+1"], [], k_vector)


def test_is_k_equivalent_to_minus_k():
    assert is_k_equivalent_to_minus_k([0.5, 0.0, 0.5])
    assert is_k_equivalent_to_minus_k([0.0, 0.0, 0.0])
    assert not is_k_equivalent_to_minus_k([1. / 3., 0.0, 0.0])


def test_satellites_and_centring():
    index_hkl, sign_k = calc_satellite_index_hkl(
        0.01, 0.25, UCP, [0.5, 0.0, 0.5])
    assert numpy.all(sign_k == 1)
    numpy.testing.assert_allclose(numpy.mod(index_hkl[0], 1.), 0.5)
    i_centring = numpy.array([[0., 0.5], [0., 0.5], [0., 0.5]])
    index_i, _ = calc_satellite_index_hkl(
        0.01, 0.25, UCP, [0.5, 0.0, 0.5], i_centring)
    h = numpy.round(index_i - numpy.array([0.5, 0.0, 0.5])[:, None])
    assert numpy.all(numpy.mod(h.sum(axis=0), 2) == 0)
    assert index_i.shape[1] < index_hkl.shape[1]
    index_g, sign_g = calc_satellite_index_hkl(
        0.01, 0.25, UCP, [1. / 3., 0.0, 0.0])
    assert set(sign_g) == {1, -1}


def test_generated_multiplicity():
    rotation = numpy.stack([numpy.eye(3).flatten(), -numpy.eye(3).flatten()], axis=1)
    translation = numpy.zeros((3, 2))
    centring = numpy.array([[0., 0.5], [0., 0.5], [0., 0.5]])
    xyz = numpy.array([[0.0, 0.1], [0.0, 0.2], [0.0, 0.3]])
    numpy.testing.assert_array_equal(
        calc_generated_multiplicity(xyz, rotation, translation, centring), [2, 4])


def test_full_symm_elems_by_k_operations():
    rotation = numpy.stack([numpy.eye(3).flatten(), -numpy.eye(3).flatten()], axis=1)
    translation = numpy.array([[0., 0.5], [0., 0.], [0., 0.]])
    centring = numpy.array([[0., 0.5], [0., 0.5], [0., 0.5]])
    elems = calc_full_symm_elems_by_k_operations(rotation, translation, centring)
    assert elems.shape == (14, 4)
    # -x+1/2 with centring 1/2 gives translation 0 (mod 1)
    numpy.testing.assert_array_equal(elems[:4, 3], [0, 1, 1, 2])
    numpy.testing.assert_array_equal(elems[13], 1)


def test_f_m_k_single_atom():
    """One atom at the origin, moment along c, F_M = 0.2695 m f."""
    index_hkl = numpy.array([[0.5], [0.], [0.]])
    identity = numpy.eye(3).flatten()[:, None]
    f_m = calc_f_m_k(
        index_hkl, numpy.array([0.5, 0., 0.]), UCP, identity,
        numpy.zeros((3, 1)), identity, numpy.zeros(1), numpy.zeros((3, 1)),
        numpy.zeros((3, 1)), numpy.ones(1), numpy.zeros(1),
        numpy.full((1, 1), 0.8), numpy.array([[0.], [0.], [2.]], dtype=complex))
    numpy.testing.assert_allclose(f_m[:, 0], [0., 0., 0.2695 * 2. * 0.8], atol=1e-12)
    # Q || a*: moment along c is fully perpendicular
    f_perp = calc_f_m_perp_k(index_hkl, UCP, f_m)
    numpy.testing.assert_allclose(f_perp, f_m, atol=1e-12)
    # moment along a is parallel to Q and does not scatter
    f_m_a = calc_f_m_k(
        index_hkl, numpy.array([0.5, 0., 0.]), UCP, identity,
        numpy.zeros((3, 1)), identity, numpy.zeros(1), numpy.zeros((3, 1)),
        numpy.zeros((3, 1)), numpy.ones(1), numpy.zeros(1),
        numpy.full((1, 1), 0.8), numpy.array([[2.], [0.], [0.]], dtype=complex))
    numpy.testing.assert_allclose(calc_f_m_perp_k(index_hkl, UCP, f_m_a), 0., atol=1e-12)


def test_moments_at_positions_change_from_cell_to_cell():
    """k = (1/2, 0, 0): the moment of the atom at the origin changes sign
    from cell to cell; other positions have no moment."""
    from cryspy.A_functions_base.structure_factor_magnetic_k import calc_moments_at_positions
    fract_xyz = numpy.array([[0., 1., 2., 0.5, -1.], [0., 0., 0., 0., 0.], [0., 0., 0., 0., 1.]])
    moment = calc_moments_at_positions(
        fract_xyz, numpy.zeros((3, 1)), numpy.array([[0.], [0.], [2. + 1j]]),
        numpy.array([0.5, 0., 0.]))
    assert numpy.allclose(moment[2, [0, 1, 2, 4]], [2., -2., 2., -2.])
    assert numpy.all(numpy.isnan(moment[:, 3]))
