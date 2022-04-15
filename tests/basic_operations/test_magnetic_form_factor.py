import os
import math
import numpy

from cryspy.A_functions_base.magnetic_form_factor import \
    get_j0_j2_parameters, \
    calc_j0,\
    calc_j2,\
    calc_form_factor

na = numpy.newaxis

symbols = numpy.array(["Co3+", "O2-", "Co2+", "HH"], dtype=str)

lande_factor = numpy.array([1.8, 1.9, 2.0, 1.4], dtype=float)
kappa = numpy.array([1., 1.1, 0.9, 1.4], dtype=float)
sthovl = numpy.array([0., 0.1, 0.7], dtype=float)



j0_parameters = numpy.array([
    [0.3902,   0.99895,   0.4332,   0.],
    [12.5078,  12.09652, 14.3553,   0.],
    [ 0.6324,   0.28854,  0.5857,   0.],
    [ 4.4574,   0.12914,  4.6077,   0.],
    [-0.15,     0.11425, -0.0382,   0.],
    [ 0.0343,  -0.22968,  0.1338,   0.],
    [ 0.1272,  -0.40685,  0.0179,   0.]], dtype=float)

j2_parameters = numpy.array([
    [ 1.70580e+00,  0.00000e+00,  1.90490e+00,  0.],
    [ 8.85950e+00,  0.00000e+00,  1.16444e+01,  0.],
    [ 1.14090e+00,  0.00000e+00,  1.31590e+00,  0.],
    [ 3.30860e+00,  0.00000e+00,  4.35740e+00,  0.],
    [ 1.47400e-01,  0.00000e+00,  3.14600e-01,  0.],
    [ 1.08990e+00,  0.00000e+00,  1.64530e+00,  0.],
    [-2.50000e-03,  0.00000e+00,  1.70000e-03,  0.]], dtype=float)

j0_coeff = numpy.array([
    [ 9.99800000e-01,  9.94890000e-01,  9.98600000e-01,  0.00000000e+00],
    [ 9.26405285e-01,  8.99762622e-01,  8.95919537e-01,  0.00000000e+00],
    [ 5.17433999e-02, -1.76190696e-04,  1.88124323e-02,  0.00000000e+00]], dtype=float)


j2_coeff = numpy.array([
    [0., 0.02808248, 0.16249887],
    [0., 0.        , 0.        ],
    [0., 0.03958993, 0.12941214],
    [0., 0.        , 0.        ]], dtype=float)




form_factor = numpy.array([
    [ 9.99800000e-01,  9.29525560e-01,  6.97988300e-02],
    [ 9.94890000e-01,  8.99762622e-01, -1.76190696e-04],
    [ 9.98600000e-01,  8.95919537e-01,  1.88124323e-02],
    [ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00]], dtype=float)


def test_get_j0_j2_parameters():
    res_1, res_2 = get_j0_j2_parameters(symbols)
    print("res_1:", res_1)
    print(j0_parameters)
    print("res_2:", res_2)
    print(j2_parameters)

    assert numpy.all(numpy.isclose(j0_parameters, res_1))
    assert numpy.all(numpy.isclose(j2_parameters, res_2))


def test_calc_j0():
    res, dder = calc_j0(sthovl[:, na], kappa[na,:], j0_parameters[:, na, :], flag_sthovl=True, flag_kappa=True)
    print("res:", res)
    print(j0_coeff)

    assert numpy.all(numpy.isclose(j0_coeff, res))


def test_calc_j2():
    res, dder = calc_j2(sthovl[na, :], kappa[:, na], j2_parameters[:, :, na], flag_sthovl=True, flag_kappa=True)
    print("res:", res)
    print(j2_coeff)

    assert numpy.all(numpy.isclose(j2_coeff, res))


def test_calc_form_factor():
    res, dder = calc_form_factor(sthovl[na, :], lande_factor[:, na], kappa[:, na], j0_parameters[:, :, na],  j2_parameters[:, :, na],
        flag_only_orbital=False, flag_sthovl=True, flag_kappa=True, flag_lande_factor=True)
    print("res:", res)
    print(form_factor)

    assert numpy.all(numpy.isclose(form_factor, res))




