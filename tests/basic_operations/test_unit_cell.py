import os
import math
import numpy

from cryspy.A_functions_base.unit_cell import \
    calc_reciprocal_by_unit_cell_parameters, \
    calc_phi_sq_by_unit_cell_parameters, \
    calc_phi_by_unit_cell_parameters, \
    calc_volume_uc_by_unit_cell_parameters, \
    calc_volume_ruc_by_unit_cell_parameters, \
    calc_m_g_by_unit_cell_parameters, \
    calc_m_g_norm_by_unit_cell_parameters, \
    calc_m_reciprocal_g_by_unit_cell_parameters, \
    calc_m_reciprocal_g_norm_by_unit_cell_parameters, \
    calc_m_b_by_unit_cell_parameters, \
    calc_m_inv_b_by_unit_cell_parameters, \
    calc_m_b_norm_by_unit_cell_parameters, \
    calc_m_inv_b_norm_by_unit_cell_parameters, \
    calc_m_m_by_unit_cell_parameters, \
    calc_m_inv_m_by_unit_cell_parameters, \
    calc_m_m_norm_by_unit_cell_parameters, \
    calc_m_inv_m_norm_by_unit_cell_parameters, \
    calc_q_ccs_by_unit_cell_parameters, \
    calc_eq_ccs_by_unit_cell_parameters, \
    calc_inv_d_by_unit_cell_parameters, \
    calc_sthovl_by_unit_cell_parameters, \
    transform_quadratic_form_reciprocal_to_ccs, \
    transform_quadratic_form_ccs_to_reciprocal, \
    transform_quadratic_form_reciprocal_norm_to_ccs, \
    transform_quadratic_form_ccs_to_reciprocal_norm, \
    transform_quadratic_form_direct_to_ccs, \
    transform_quadratic_form_ccs_to_direct, \
    transform_quadratic_form_direct_norm_to_ccs, \
    transform_quadratic_form_ccs_to_direct_norm, \
    calc_m_inv_m_b_norm_by_unit_cell_parameters, \
    calc_m_inv_b_norm_m_by_unit_cell_parameters, \
    calc_matrix_t

abs_tol = 0.00001

unit_cell_parameters = numpy.array(
    [3.4, 7.5, 12,
    87.0*numpy.pi/180., 102.*numpy.pi/180., 124.*numpy.pi/180.], dtype=float)

flag_unit_cell_parameters = True

reciprocal_unit_cell_parameters = numpy.array(
    [0.36332942, 0.16133114, 0.08546101, 1.49188214, 1.35330423, 0.97480368], dtype=float)

volume_uc = 247.36961584064142
volume_ruc = 0.00404
phi_sq = 0.653506416775028
phi = numpy.sqrt(phi_sq)

m_direct_g_ij = numpy.array([11.56, 56.25, 144.0, -14.25942, -8.48280, 4.71024], dtype=float)
m_direct_g_norm_ij = numpy.array([1., 1., 1., -0.5591929, -0.20791169, 0.05233596], dtype=float)
m_reciprocal_g_ij = numpy.array([0.13201, 0.02603, 0.00730, 0.03290, 0.00670, 0.00109], dtype=float)
m_reciprocal_g_norm_ij = numpy.array([1., 1., 1., 0.56133054, 0.21578149, 0.07883231], dtype=float)

m_b_ij = numpy.array([ 0.36332942,  0.0905601,  0.0184409,   0., 0.13351631, -0.00436731, 0., 0., 0.08333333], dtype=float)
m_inv_b_ij = numpy.array([2.75232, -1.86682, -0.70690, 0.0, 7.48972, 0.39252, 0., 0., 12.0], dtype=float)
m_b_norm_ij = numpy.array([1.0, 0.56133, 0.21578, 0., 0.82759, -0.05110, 0., 0., 0.97510], dtype=float)
m_inv_b_norm_ij = numpy.array([1.0, -0.67827, -0.25684, 0., 1.20833, 0.06333, 0., 0., 1.02553], dtype=float)

m_m_ij = numpy.array([2.75232, 0., 0., -1.86682, 7.48972, 0., -0.70690, 0.39252, 12.0], dtype=float)
m_inv_m_ij = numpy.array([0.36333, 0., 0., 0.09056, 0.13352, 0., 0.01844, -0.00437, 0.08333], dtype=float)
m_m_norm_ij = numpy.array([0.80951, 0., 0., -0.54906, 0.99863, 0., -0.20791, 0.05234, 1.0], dtype=float)
m_inv_m_norm_ij = numpy.array([1.23532, 0., 0., 0.67920, 1.00137, 0., 0.22129, -0.05241, 1.0], dtype=float)

m_q_reciprocal_norm_ij = numpy.array([3.1, 2.9, 5.7, 0.4, -0.7, 5.2], dtype=float)
m_q_xyz_ij = numpy.array([3.1, 5.00464377,  7.24204832, -1.61930656, -1.48873842,  7.55124895], dtype=float)
m_q_reciprocal_ij = numpy.array([0.40922562,  0.07548044,  0.04163043,  0.02344654, -0.02173535,  0.07169512], dtype=float)
m_q_direct_ij = numpy.array([86.90688, 326.25521, 1042.85496, -152.48631, -279.76439, 712.79278], dtype=float)
m_q_direct_norm_ij = numpy.array([7.51789858,  5.80008894,  7.24205,    -5.9798562,  -6.8569725,   7.91992089], dtype=float)

index_hkl = numpy.array([
    [0,0,2,0,1],
    [0,2,0,0,2],
    [2,0,0,0,3]], dtype=int)

index_hkl_2 = numpy.array([
    [0,0,2,-8,1, 8, 3],
    [0,2,0, 9,2,-9,-9],
    [2,0,0, 1,3, 1,17]], dtype=int)

q_ccs = numpy.array([
    [0.03688181,  0.18112019,  0.72665883,   0.0,  0.59977232],
    [-0.00873463,  0.26703263,  0.,          0.0,    0.25393068],
    [ 0.16666667,  0.,          0.,          0.0,  0.25]], dtype=float)

eq_ccs = numpy.array([
    [0.21578149, 0.56133054, 1.,  0.0, 0.85971072],
    [-0.05110301, 0.8275917, 0., 0.0, 0.363983],
    [0.9751035,   0., 0., 0.0, 0.35834878]], dtype=float)

inv_d = numpy.array([0.17092203, 0.32266228, 0.72665883, 0.0, 0.68715675], dtype=float)

sthovl = numpy.array([0.08546101, 0.16133114, 0.36332942, 0.0, 0.34357837], dtype=float)

m_inv_m_b_norm_ij = numpy.array([
    0.36332942, 2.03947897e-01, 7.83997613e-02,
    9.05600957e-02, 0.16133114, 1.27181067e-02,
    1.84409045e-02, 6.73708925e-03, 0.08546101], dtype=float)

m_inv_b_norm_m_ij = numpy.array([
     4.20008805, -5.1808664,  -3.08204968,
     -2.30048833, 9.07487662, 0.75990775,
    -0.72494843, 0.40254155, 12.30638596], dtype=float)


def test_calc_reciprocal_by_unit_cell_parameters():
    res, dder = calc_reciprocal_by_unit_cell_parameters(
        unit_cell_parameters,
        flag_unit_cell_parameters=flag_unit_cell_parameters)

    print("res: ", res)
    print(reciprocal_unit_cell_parameters)
    assert math.isclose(res[0], reciprocal_unit_cell_parameters[0], abs_tol=abs_tol)
    assert math.isclose(res[1], reciprocal_unit_cell_parameters[1], abs_tol=abs_tol)
    assert math.isclose(res[2], reciprocal_unit_cell_parameters[2], abs_tol=abs_tol)
    assert math.isclose(res[3], reciprocal_unit_cell_parameters[3], abs_tol=abs_tol)
    assert math.isclose(res[4], reciprocal_unit_cell_parameters[4], abs_tol=abs_tol)
    assert math.isclose(res[5], reciprocal_unit_cell_parameters[5], abs_tol=abs_tol)


def test_calc_phi_sq_by_unit_cell_parameters():
    res, dder = calc_phi_sq_by_unit_cell_parameters(
        unit_cell_parameters,
        flag_unit_cell_parameters=flag_unit_cell_parameters)

    print("res: ", res)
    print(phi_sq)
    assert math.isclose(res, phi_sq, abs_tol=abs_tol)

def test_calc_phi_by_unit_cell_parameters():
    res, dder = calc_phi_by_unit_cell_parameters(
        unit_cell_parameters,
        flag_unit_cell_parameters=flag_unit_cell_parameters)

    print("res: ", res)
    print(phi)
    assert math.isclose(res, phi, abs_tol=abs_tol)

def test_calc_volume_uc_by_unit_cell_parameters():
    res, dder = calc_volume_uc_by_unit_cell_parameters(
        unit_cell_parameters,
        flag_unit_cell_parameters=flag_unit_cell_parameters)

    print("res: ", res)
    print(volume_uc)
    assert math.isclose(res, volume_uc, abs_tol=abs_tol)

def test_calc_volume_ruc_by_unit_cell_parameters():
    res, dder = calc_volume_ruc_by_unit_cell_parameters(
        unit_cell_parameters,
        flag_unit_cell_parameters=flag_unit_cell_parameters)

    print("res: ", res)
    print(volume_ruc)
    assert math.isclose(res, volume_ruc, abs_tol=abs_tol)

def test_calc_m_g_by_unit_cell_parameters():
    res, dder = calc_m_g_by_unit_cell_parameters(
        unit_cell_parameters,
        flag_unit_cell_parameters=flag_unit_cell_parameters)

    print("res: ", res)
    print(m_direct_g_ij)
    assert math.isclose(res[0], m_direct_g_ij[0], abs_tol=abs_tol)
    assert math.isclose(res[1], m_direct_g_ij[1], abs_tol=abs_tol)
    assert math.isclose(res[2], m_direct_g_ij[2], abs_tol=abs_tol)
    assert math.isclose(res[3], m_direct_g_ij[3], abs_tol=abs_tol)
    assert math.isclose(res[4], m_direct_g_ij[4], abs_tol=abs_tol)
    assert math.isclose(res[5], m_direct_g_ij[5], abs_tol=abs_tol)

def test_calc_m_g_norm_by_unit_cell_parameters():
    res, dder = calc_m_g_norm_by_unit_cell_parameters(
        unit_cell_parameters,
        flag_unit_cell_parameters=flag_unit_cell_parameters)

    print("res: ", res)
    print(m_direct_g_norm_ij)
    assert math.isclose(res[0], m_direct_g_norm_ij[0], abs_tol=abs_tol)
    assert math.isclose(res[1], m_direct_g_norm_ij[1], abs_tol=abs_tol)
    assert math.isclose(res[2], m_direct_g_norm_ij[2], abs_tol=abs_tol)
    assert math.isclose(res[3], m_direct_g_norm_ij[3], abs_tol=abs_tol)
    assert math.isclose(res[4], m_direct_g_norm_ij[4], abs_tol=abs_tol)
    assert math.isclose(res[5], m_direct_g_norm_ij[5], abs_tol=abs_tol)

def test_calc_m_reciprocal_g_by_unit_cell_parameters():
    res, dder = calc_m_reciprocal_g_by_unit_cell_parameters(
        unit_cell_parameters,
        flag_unit_cell_parameters=flag_unit_cell_parameters)

    print("res: ", res)
    print(m_reciprocal_g_ij)
    assert math.isclose(res[0], m_reciprocal_g_ij[0], abs_tol=abs_tol)
    assert math.isclose(res[1], m_reciprocal_g_ij[1], abs_tol=abs_tol)
    assert math.isclose(res[2], m_reciprocal_g_ij[2], abs_tol=abs_tol)
    assert math.isclose(res[3], m_reciprocal_g_ij[3], abs_tol=abs_tol)
    assert math.isclose(res[4], m_reciprocal_g_ij[4], abs_tol=abs_tol)
    assert math.isclose(res[5], m_reciprocal_g_ij[5], abs_tol=abs_tol)


def test_calc_m_reciprocal_g_norm_by_unit_cell_parameters():
    res, dder = calc_m_reciprocal_g_norm_by_unit_cell_parameters(
        unit_cell_parameters,
        flag_unit_cell_parameters=flag_unit_cell_parameters)

    print("res: ", res)
    print(m_reciprocal_g_norm_ij)
    assert math.isclose(res[0], m_reciprocal_g_norm_ij[0], abs_tol=abs_tol)
    assert math.isclose(res[1], m_reciprocal_g_norm_ij[1], abs_tol=abs_tol)
    assert math.isclose(res[2], m_reciprocal_g_norm_ij[2], abs_tol=abs_tol)
    assert math.isclose(res[3], m_reciprocal_g_norm_ij[3], abs_tol=abs_tol)
    assert math.isclose(res[4], m_reciprocal_g_norm_ij[4], abs_tol=abs_tol)
    assert math.isclose(res[5], m_reciprocal_g_norm_ij[5], abs_tol=abs_tol)


def test_calc_m_b_by_unit_cell_parameters():
    res, dder = calc_m_b_by_unit_cell_parameters(
        unit_cell_parameters,
        flag_unit_cell_parameters=flag_unit_cell_parameters)

    print("res: ", res)
    print(m_b_ij)
    assert math.isclose(res[0], m_b_ij[0], abs_tol=abs_tol)
    assert math.isclose(res[1], m_b_ij[1], abs_tol=abs_tol)
    assert math.isclose(res[2], m_b_ij[2], abs_tol=abs_tol)
    assert math.isclose(res[3], m_b_ij[3], abs_tol=abs_tol)
    assert math.isclose(res[4], m_b_ij[4], abs_tol=abs_tol)
    assert math.isclose(res[5], m_b_ij[5], abs_tol=abs_tol)
    assert math.isclose(res[6], m_b_ij[6], abs_tol=abs_tol)
    assert math.isclose(res[7], m_b_ij[7], abs_tol=abs_tol)
    assert math.isclose(res[8], m_b_ij[8], abs_tol=abs_tol)


def test_calc_m_inv_b_by_unit_cell_parameters():
    res, dder = calc_m_inv_b_by_unit_cell_parameters(
        unit_cell_parameters,
        flag_unit_cell_parameters=flag_unit_cell_parameters)

    print("res: ", res)
    print(m_inv_b_ij)
    assert math.isclose(res[0], m_inv_b_ij[0], abs_tol=abs_tol)
    assert math.isclose(res[1], m_inv_b_ij[1], abs_tol=abs_tol)
    assert math.isclose(res[2], m_inv_b_ij[2], abs_tol=abs_tol)
    assert math.isclose(res[3], m_inv_b_ij[3], abs_tol=abs_tol)
    assert math.isclose(res[4], m_inv_b_ij[4], abs_tol=abs_tol)
    assert math.isclose(res[5], m_inv_b_ij[5], abs_tol=abs_tol)
    assert math.isclose(res[6], m_inv_b_ij[6], abs_tol=abs_tol)
    assert math.isclose(res[7], m_inv_b_ij[7], abs_tol=abs_tol)
    assert math.isclose(res[8], m_inv_b_ij[8], abs_tol=abs_tol)


def test_calc_m_b_norm_by_unit_cell_parameters():
    res, dder = calc_m_b_norm_by_unit_cell_parameters(
        unit_cell_parameters,
        flag_unit_cell_parameters=flag_unit_cell_parameters)

    print("res: ", res)
    print(m_b_norm_ij)
    assert math.isclose(res[0], m_b_norm_ij[0], abs_tol=abs_tol)
    assert math.isclose(res[1], m_b_norm_ij[1], abs_tol=abs_tol)
    assert math.isclose(res[2], m_b_norm_ij[2], abs_tol=abs_tol)
    assert math.isclose(res[3], m_b_norm_ij[3], abs_tol=abs_tol)
    assert math.isclose(res[4], m_b_norm_ij[4], abs_tol=abs_tol)
    assert math.isclose(res[5], m_b_norm_ij[5], abs_tol=abs_tol)
    assert math.isclose(res[6], m_b_norm_ij[6], abs_tol=abs_tol)
    assert math.isclose(res[7], m_b_norm_ij[7], abs_tol=abs_tol)
    assert math.isclose(res[8], m_b_norm_ij[8], abs_tol=abs_tol)


def test_calc_m_inv_b_norm_by_unit_cell_parameters():
    res, dder = calc_m_inv_b_norm_by_unit_cell_parameters(
        unit_cell_parameters,
        flag_unit_cell_parameters=flag_unit_cell_parameters)

    print("res: ", res)
    print(m_inv_b_norm_ij)
    assert math.isclose(res[0], m_inv_b_norm_ij[0], abs_tol=abs_tol)
    assert math.isclose(res[1], m_inv_b_norm_ij[1], abs_tol=abs_tol)
    assert math.isclose(res[2], m_inv_b_norm_ij[2], abs_tol=abs_tol)
    assert math.isclose(res[3], m_inv_b_norm_ij[3], abs_tol=abs_tol)
    assert math.isclose(res[4], m_inv_b_norm_ij[4], abs_tol=abs_tol)
    assert math.isclose(res[5], m_inv_b_norm_ij[5], abs_tol=abs_tol)
    assert math.isclose(res[6], m_inv_b_norm_ij[6], abs_tol=abs_tol)
    assert math.isclose(res[7], m_inv_b_norm_ij[7], abs_tol=abs_tol)
    assert math.isclose(res[8], m_inv_b_norm_ij[8], abs_tol=abs_tol)


def test_calc_m_m_by_unit_cell_parameters():
    res, dder = calc_m_m_by_unit_cell_parameters(
        unit_cell_parameters,
        flag_unit_cell_parameters=flag_unit_cell_parameters)

    print("res: ", res)
    print(m_m_ij)
    assert math.isclose(res[0], m_m_ij[0], abs_tol=abs_tol)
    assert math.isclose(res[1], m_m_ij[1], abs_tol=abs_tol)
    assert math.isclose(res[2], m_m_ij[2], abs_tol=abs_tol)
    assert math.isclose(res[3], m_m_ij[3], abs_tol=abs_tol)
    assert math.isclose(res[4], m_m_ij[4], abs_tol=abs_tol)
    assert math.isclose(res[5], m_m_ij[5], abs_tol=abs_tol)
    assert math.isclose(res[6], m_m_ij[6], abs_tol=abs_tol)
    assert math.isclose(res[7], m_m_ij[7], abs_tol=abs_tol)
    assert math.isclose(res[8], m_m_ij[8], abs_tol=abs_tol)

def test_calc_m_inv_m_by_unit_cell_parameters():
    res, dder = calc_m_inv_m_by_unit_cell_parameters(
        unit_cell_parameters,
        flag_unit_cell_parameters=flag_unit_cell_parameters)

    print("res: ", res)
    print(m_inv_m_ij)
    assert math.isclose(res[0], m_inv_m_ij[0], abs_tol=abs_tol)
    assert math.isclose(res[1], m_inv_m_ij[1], abs_tol=abs_tol)
    assert math.isclose(res[2], m_inv_m_ij[2], abs_tol=abs_tol)
    assert math.isclose(res[3], m_inv_m_ij[3], abs_tol=abs_tol)
    assert math.isclose(res[4], m_inv_m_ij[4], abs_tol=abs_tol)
    assert math.isclose(res[5], m_inv_m_ij[5], abs_tol=abs_tol)
    assert math.isclose(res[6], m_inv_m_ij[6], abs_tol=abs_tol)
    assert math.isclose(res[7], m_inv_m_ij[7], abs_tol=abs_tol)
    assert math.isclose(res[8], m_inv_m_ij[8], abs_tol=abs_tol)


def test_calc_m_m_norm_by_unit_cell_parameters():
    res, dder = calc_m_m_norm_by_unit_cell_parameters(
        unit_cell_parameters,
        flag_unit_cell_parameters=flag_unit_cell_parameters)

    print("res: ", res)
    print(m_m_norm_ij)
    assert math.isclose(res[0], m_m_norm_ij[0], abs_tol=abs_tol)
    assert math.isclose(res[1], m_m_norm_ij[1], abs_tol=abs_tol)
    assert math.isclose(res[2], m_m_norm_ij[2], abs_tol=abs_tol)
    assert math.isclose(res[3], m_m_norm_ij[3], abs_tol=abs_tol)
    assert math.isclose(res[4], m_m_norm_ij[4], abs_tol=abs_tol)
    assert math.isclose(res[5], m_m_norm_ij[5], abs_tol=abs_tol)
    assert math.isclose(res[6], m_m_norm_ij[6], abs_tol=abs_tol)
    assert math.isclose(res[7], m_m_norm_ij[7], abs_tol=abs_tol)
    assert math.isclose(res[8], m_m_norm_ij[8], abs_tol=abs_tol)


def test_calc_m_inv_m_norm_by_unit_cell_parameters():
    res, dder = calc_m_inv_m_norm_by_unit_cell_parameters(
        unit_cell_parameters,
        flag_unit_cell_parameters=flag_unit_cell_parameters)

    print("res: ", res)
    print(m_inv_m_norm_ij)
    assert math.isclose(res[0], m_inv_m_norm_ij[0], abs_tol=abs_tol)
    assert math.isclose(res[1], m_inv_m_norm_ij[1], abs_tol=abs_tol)
    assert math.isclose(res[2], m_inv_m_norm_ij[2], abs_tol=abs_tol)
    assert math.isclose(res[3], m_inv_m_norm_ij[3], abs_tol=abs_tol)
    assert math.isclose(res[4], m_inv_m_norm_ij[4], abs_tol=abs_tol)
    assert math.isclose(res[5], m_inv_m_norm_ij[5], abs_tol=abs_tol)
    assert math.isclose(res[6], m_inv_m_norm_ij[6], abs_tol=abs_tol)
    assert math.isclose(res[7], m_inv_m_norm_ij[7], abs_tol=abs_tol)
    assert math.isclose(res[8], m_inv_m_norm_ij[8], abs_tol=abs_tol)

def test_calc_q_ccs_by_unit_cell_parameters():
    res, dder = calc_q_ccs_by_unit_cell_parameters(
        index_hkl, unit_cell_parameters,
        flag_unit_cell_parameters=False)

    print("res: ", res)
    print(q_ccs)
    assert math.isclose(res[0][0], q_ccs[0][0], abs_tol=abs_tol)
    assert math.isclose(res[0][1], q_ccs[0][1], abs_tol=abs_tol)
    assert math.isclose(res[0][2], q_ccs[0][2], abs_tol=abs_tol)
    assert math.isclose(res[0][3], q_ccs[0][3], abs_tol=abs_tol)
    assert math.isclose(res[0][4], q_ccs[0][4], abs_tol=abs_tol)

    assert math.isclose(res[1][0], q_ccs[1][0], abs_tol=abs_tol)
    assert math.isclose(res[1][1], q_ccs[1][1], abs_tol=abs_tol)
    assert math.isclose(res[1][2], q_ccs[1][2], abs_tol=abs_tol)
    assert math.isclose(res[1][3], q_ccs[1][3], abs_tol=abs_tol)
    assert math.isclose(res[1][4], q_ccs[1][4], abs_tol=abs_tol)

    assert math.isclose(res[2][0], q_ccs[2][0], abs_tol=abs_tol)
    assert math.isclose(res[2][1], q_ccs[2][1], abs_tol=abs_tol)
    assert math.isclose(res[2][2], q_ccs[2][2], abs_tol=abs_tol)
    assert math.isclose(res[2][3], q_ccs[2][3], abs_tol=abs_tol)
    assert math.isclose(res[2][4], q_ccs[2][4], abs_tol=abs_tol)


def test_calc_eq_ccs_by_unit_cell_parameters():
    res, dder = calc_eq_ccs_by_unit_cell_parameters(
        index_hkl, unit_cell_parameters,
        flag_unit_cell_parameters=False)

    print("res: ", res)
    print(eq_ccs)
    assert math.isclose(res[0][0], eq_ccs[0][0], abs_tol=abs_tol)
    assert math.isclose(res[0][1], eq_ccs[0][1], abs_tol=abs_tol)
    assert math.isclose(res[0][2], eq_ccs[0][2], abs_tol=abs_tol)
    assert math.isclose(res[0][3], eq_ccs[0][3], abs_tol=abs_tol)
    assert math.isclose(res[0][4], eq_ccs[0][4], abs_tol=abs_tol)

    assert math.isclose(res[1][0], eq_ccs[1][0], abs_tol=abs_tol)
    assert math.isclose(res[1][1], eq_ccs[1][1], abs_tol=abs_tol)
    assert math.isclose(res[1][2], eq_ccs[1][2], abs_tol=abs_tol)
    assert math.isclose(res[1][3], eq_ccs[1][3], abs_tol=abs_tol)
    assert math.isclose(res[1][4], eq_ccs[1][4], abs_tol=abs_tol)

    assert math.isclose(res[2][0], eq_ccs[2][0], abs_tol=abs_tol)
    assert math.isclose(res[2][1], eq_ccs[2][1], abs_tol=abs_tol)
    assert math.isclose(res[2][2], eq_ccs[2][2], abs_tol=abs_tol)
    assert math.isclose(res[2][3], eq_ccs[2][3], abs_tol=abs_tol)
    assert math.isclose(res[2][4], eq_ccs[2][4], abs_tol=abs_tol)


def test_calc_inv_d_by_unit_cell_parameters():
    res, dder = calc_inv_d_by_unit_cell_parameters(
        index_hkl, unit_cell_parameters, 
        flag_unit_cell_parameters=flag_unit_cell_parameters)

    print("res: ", res)
    print(inv_d)
    assert math.isclose(res[0], inv_d[0], abs_tol=abs_tol)
    assert math.isclose(res[1], inv_d[1], abs_tol=abs_tol)
    assert math.isclose(res[2], inv_d[2], abs_tol=abs_tol)
    assert math.isclose(res[3], inv_d[3], abs_tol=abs_tol)
    assert math.isclose(res[4], inv_d[4], abs_tol=abs_tol)


def test_calc_sthovl_by_unit_cell_parameters():
    res, dder = calc_sthovl_by_unit_cell_parameters(
        index_hkl, unit_cell_parameters,
        flag_unit_cell_parameters=flag_unit_cell_parameters)

    print("res: ", res)
    print(sthovl)
    assert math.isclose(res[0], sthovl[0], abs_tol=abs_tol)
    assert math.isclose(res[1], sthovl[1], abs_tol=abs_tol)
    assert math.isclose(res[2], sthovl[2], abs_tol=abs_tol)
    assert math.isclose(res[3], sthovl[3], abs_tol=abs_tol)
    assert math.isclose(res[4], sthovl[4], abs_tol=abs_tol)


def test_transform_quadratic_form_reciprocal_to_ccs():
    res, dder = transform_quadratic_form_reciprocal_to_ccs(
        m_q_reciprocal_ij, unit_cell_parameters)

    print("res: ", res)
    print(m_q_xyz_ij)
    assert math.isclose(res[0], m_q_xyz_ij[0], abs_tol=abs_tol)
    assert math.isclose(res[1], m_q_xyz_ij[1], abs_tol=abs_tol)
    assert math.isclose(res[2], m_q_xyz_ij[2], abs_tol=abs_tol)
    assert math.isclose(res[3], m_q_xyz_ij[3], abs_tol=abs_tol)
    assert math.isclose(res[4], m_q_xyz_ij[4], abs_tol=abs_tol)
    assert math.isclose(res[5], m_q_xyz_ij[5], abs_tol=abs_tol)


def test_transform_quadratic_form_ccs_to_reciprocal():
    res, dder = transform_quadratic_form_ccs_to_reciprocal(
        m_q_xyz_ij, unit_cell_parameters)

    print("res: ", res)
    print(m_q_reciprocal_ij)
    assert math.isclose(res[0], m_q_reciprocal_ij[0], abs_tol=abs_tol)
    assert math.isclose(res[1], m_q_reciprocal_ij[1], abs_tol=abs_tol)
    assert math.isclose(res[2], m_q_reciprocal_ij[2], abs_tol=abs_tol)
    assert math.isclose(res[3], m_q_reciprocal_ij[3], abs_tol=abs_tol)
    assert math.isclose(res[4], m_q_reciprocal_ij[4], abs_tol=abs_tol)
    assert math.isclose(res[5], m_q_reciprocal_ij[5], abs_tol=abs_tol)


def test_transform_quadratic_form_reciprocal_norm_to_ccs():
    res, dder = transform_quadratic_form_reciprocal_norm_to_ccs(
        m_q_reciprocal_norm_ij, unit_cell_parameters)

    print("res: ", res)
    print(m_q_xyz_ij)
    assert math.isclose(res[0], m_q_xyz_ij[0], abs_tol=abs_tol)
    assert math.isclose(res[1], m_q_xyz_ij[1], abs_tol=abs_tol)
    assert math.isclose(res[2], m_q_xyz_ij[2], abs_tol=abs_tol)
    assert math.isclose(res[3], m_q_xyz_ij[3], abs_tol=abs_tol)
    assert math.isclose(res[4], m_q_xyz_ij[4], abs_tol=abs_tol)
    assert math.isclose(res[5], m_q_xyz_ij[5], abs_tol=abs_tol)


def test_transform_quadratic_form_ccs_to_reciprocal_norm():
    res, dder = transform_quadratic_form_ccs_to_reciprocal_norm(
        m_q_xyz_ij, unit_cell_parameters)

    print("res: ", res)
    print(m_q_reciprocal_norm_ij)
    assert math.isclose(res[0], m_q_reciprocal_norm_ij[0], abs_tol=abs_tol)
    assert math.isclose(res[1], m_q_reciprocal_norm_ij[1], abs_tol=abs_tol)
    assert math.isclose(res[2], m_q_reciprocal_norm_ij[2], abs_tol=abs_tol)
    assert math.isclose(res[3], m_q_reciprocal_norm_ij[3], abs_tol=abs_tol)
    assert math.isclose(res[4], m_q_reciprocal_norm_ij[4], abs_tol=abs_tol)
    assert math.isclose(res[5], m_q_reciprocal_norm_ij[5], abs_tol=abs_tol)


def test_transform_quadratic_form_direct_to_ccs():
    res, dder = transform_quadratic_form_direct_to_ccs(
        m_q_direct_ij, unit_cell_parameters)

    print("res: ", res)
    print(m_q_xyz_ij)
    assert math.isclose(res[0], m_q_xyz_ij[0], abs_tol=abs_tol)
    assert math.isclose(res[1], m_q_xyz_ij[1], abs_tol=abs_tol)
    assert math.isclose(res[2], m_q_xyz_ij[2], abs_tol=abs_tol)
    assert math.isclose(res[3], m_q_xyz_ij[3], abs_tol=abs_tol)
    assert math.isclose(res[4], m_q_xyz_ij[4], abs_tol=abs_tol)
    assert math.isclose(res[5], m_q_xyz_ij[5], abs_tol=abs_tol)


def test_transform_quadratic_form_ccs_to_direct():
    res, dder = transform_quadratic_form_ccs_to_direct(
        m_q_xyz_ij, unit_cell_parameters)

    print("res: ", res)
    print(m_q_direct_ij)
    assert math.isclose(res[0], m_q_direct_ij[0], abs_tol=abs_tol)
    assert math.isclose(res[1], m_q_direct_ij[1], abs_tol=abs_tol)
    assert math.isclose(res[2], m_q_direct_ij[2], abs_tol=abs_tol)
    assert math.isclose(res[3], m_q_direct_ij[3], abs_tol=abs_tol)
    assert math.isclose(res[4], m_q_direct_ij[4], abs_tol=abs_tol)
    assert math.isclose(res[5], m_q_direct_ij[5], abs_tol=abs_tol)


def test_transform_quadratic_form_direct_norm_to_ccs():
    res, dder = transform_quadratic_form_direct_norm_to_ccs(
        m_q_direct_norm_ij, unit_cell_parameters)

    print("res: ", res)
    print(m_q_xyz_ij)
    assert math.isclose(res[0], m_q_xyz_ij[0], abs_tol=abs_tol)
    assert math.isclose(res[1], m_q_xyz_ij[1], abs_tol=abs_tol)
    assert math.isclose(res[2], m_q_xyz_ij[2], abs_tol=abs_tol)
    assert math.isclose(res[3], m_q_xyz_ij[3], abs_tol=abs_tol)
    assert math.isclose(res[4], m_q_xyz_ij[4], abs_tol=abs_tol)
    assert math.isclose(res[5], m_q_xyz_ij[5], abs_tol=abs_tol)


def test_transform_quadratic_form_ccs_to_direct_norm():
    res, dder = transform_quadratic_form_ccs_to_direct_norm(
        m_q_xyz_ij, unit_cell_parameters)

    print("res: ", res)
    print(m_q_direct_norm_ij)
    assert math.isclose(res[0], m_q_direct_norm_ij[0], abs_tol=abs_tol)
    assert math.isclose(res[1], m_q_direct_norm_ij[1], abs_tol=abs_tol)
    assert math.isclose(res[2], m_q_direct_norm_ij[2], abs_tol=abs_tol)
    assert math.isclose(res[3], m_q_direct_norm_ij[3], abs_tol=abs_tol)
    assert math.isclose(res[4], m_q_direct_norm_ij[4], abs_tol=abs_tol)
    assert math.isclose(res[5], m_q_direct_norm_ij[5], abs_tol=abs_tol)


def test_calc_m_inv_m_b_norm_by_unit_cell_parameters():
    res, dder = calc_m_inv_m_b_norm_by_unit_cell_parameters(
        unit_cell_parameters, flag_unit_cell_parameters=flag_unit_cell_parameters)

    print("res: ", res)
    print(m_inv_m_b_norm_ij)
    assert math.isclose(res[0], m_inv_m_b_norm_ij[0], abs_tol=abs_tol)
    assert math.isclose(res[1], m_inv_m_b_norm_ij[1], abs_tol=abs_tol)
    assert math.isclose(res[2], m_inv_m_b_norm_ij[2], abs_tol=abs_tol)
    assert math.isclose(res[3], m_inv_m_b_norm_ij[3], abs_tol=abs_tol)
    assert math.isclose(res[4], m_inv_m_b_norm_ij[4], abs_tol=abs_tol)
    assert math.isclose(res[5], m_inv_m_b_norm_ij[5], abs_tol=abs_tol)
    assert math.isclose(res[6], m_inv_m_b_norm_ij[6], abs_tol=abs_tol)
    assert math.isclose(res[7], m_inv_m_b_norm_ij[7], abs_tol=abs_tol)
    assert math.isclose(res[8], m_inv_m_b_norm_ij[8], abs_tol=abs_tol)


def test_calc_m_inv_b_norm_m_by_unit_cell_parameters():
    res, dder = calc_m_inv_b_norm_m_by_unit_cell_parameters(
        unit_cell_parameters, flag_unit_cell_parameters=False)

    print("res: ", res)
    print(m_inv_b_norm_m_ij)
    assert math.isclose(res[0], m_inv_b_norm_m_ij[0], abs_tol=abs_tol)
    assert math.isclose(res[1], m_inv_b_norm_m_ij[1], abs_tol=abs_tol)
    assert math.isclose(res[2], m_inv_b_norm_m_ij[2], abs_tol=abs_tol)
    assert math.isclose(res[3], m_inv_b_norm_m_ij[3], abs_tol=abs_tol)
    assert math.isclose(res[4], m_inv_b_norm_m_ij[4], abs_tol=abs_tol)
    assert math.isclose(res[5], m_inv_b_norm_m_ij[5], abs_tol=abs_tol)
    assert math.isclose(res[6], m_inv_b_norm_m_ij[6], abs_tol=abs_tol)
    assert math.isclose(res[7], m_inv_b_norm_m_ij[7], abs_tol=abs_tol)
    assert math.isclose(res[8], m_inv_b_norm_m_ij[8], abs_tol=abs_tol)

def test_calc_matrix_t():
    res, dder_m_t = calc_matrix_t(
        index_hkl_2, unit_cell_parameters, flag_unit_cell_parameters=False)
    
    eq_ccs, dder_eq_ccs = calc_eq_ccs_by_unit_cell_parameters(
        index_hkl_2, unit_cell_parameters, flag_unit_cell_parameters=False)
    eq_x = res[0]*eq_ccs[0]+res[1]*eq_ccs[1]+res[2]*eq_ccs[2]
    eq_y = res[3]*eq_ccs[0]+res[4]*eq_ccs[1]+res[5]*eq_ccs[2]
    eq_z = res[6]*eq_ccs[0]+res[7]*eq_ccs[1]+res[8]*eq_ccs[2]

    assert numpy.all(numpy.isclose(eq_x, 0.))
    assert numpy.all(numpy.isclose(eq_y, 0.))
    assert numpy.all(numpy.isclose(eq_z, 1.))
