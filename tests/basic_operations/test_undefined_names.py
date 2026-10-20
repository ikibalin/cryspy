"""Code paths that used names which were never defined (ruff F821)."""
import numpy

import cryspy
from cryspy.A_functions_base.matrix_operations import calc_vt_q_v, transpose_m
from cryspy.C_item_loop_classes.cl_1_tof_intensity_incident import TOFIntensityIncident


def test_transpose_m_complex_with_derivatives():
    m = numpy.arange(9, dtype=float) + 1j * numpy.arange(9, dtype=float)[::-1]
    res, dder = transpose_m(m, flag_m=True)
    assert numpy.allclose(res, m.reshape(3, 3).transpose().flatten())
    assert dder["m_real"].shape == (9, 9)
    assert dder["m_imag"].dtype == float


def test_calc_vt_q_v_derivative_over_q():
    # q is a quadratic form stored as 11, 22, 33, 12, 13, 23
    q = numpy.array([8.74, 6.7, 1.7, 0.74, 7.7, -6.7])
    v = numpy.array([0.74, 7.7, -6.7])
    res, dder = calc_vt_q_v(q, v, flag_q=True)

    q_full = numpy.array([
        [q[0], q[3], q[4]],
        [q[3], q[1], q[5]],
        [q[4], q[5], q[2]]])
    assert numpy.isclose(res, v @ q_full @ v)

    step = 1e-6
    der_numerical = numpy.zeros(6)
    for i_q in range(6):
        q_plus, q_minus = q.copy(), q.copy()
        q_plus[i_q] += step
        q_minus[i_q] -= step
        der_numerical[i_q] = (calc_vt_q_v(q_plus, v)[0] - calc_vt_q_v(q_minus, v)[0]) / (2. * step)
    assert dder["q"].shape == (6, )
    assert numpy.allclose(dder["q"], der_numerical, rtol=1e-6)


def test_tof_intensity_incident_get_dictionary():
    item = TOFIntensityIncident(
        a0=1., a1=2., a2=3., a3=4., a4=5., a5=6., a6=7., a7=8., a8=9.,
        spectrum="Maxwell")
    res = item.get_dictionary()
    assert res["spectrum_incident_type"] == "Maxwell"
    assert numpy.allclose(res["spectrum_incident_coefficients"], numpy.arange(1., 10.))
    assert not numpy.any(res["flags_spectrum_incident_coefficients"])


def test_magnetization_ellipsoid_without_susceptibility():
    crystal = cryspy.Crystal(data_name="no_susceptibility")
    assert crystal.calc_magnetization_ellipsoid() == []
