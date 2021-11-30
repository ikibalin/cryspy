import os
import math
import numpy

from cryspy.A_functions_base.local_susceptibility import \
    calc_magnetization_ellipsoid_as_u, \
    calc_magnetization_ellipsoid_axes,\
    calc_ellipsoid_factor, \
    calc_chi_direct_norm_with_symmetry

na = numpy.newaxis

abs_tol = 0.00001

unit_cell_parameters = numpy.array(
    [3.4, 7.5, 12,
    87.0*numpy.pi/180., 102.*numpy.pi/180., 124.*numpy.pi/180.], dtype=float)

flag_unit_cell_parameters = True

susceptibility_ij = numpy.array([3.1, 2.9, 5.7, 0.4, -0.7, 5.2], dtype=float)

reciprocal_u_norm = numpy.array([10.67146199, 40.18752035, 58.88475685,  7.31809478,  8.62285557, 47.60932752], dtype=float)

ellipsoid_factor = -0.006104671061403277
# reciprocal_u_norm2 = numpy.array([1.09999359,  11.4156045,  113.944081,  0., 0., 0.], dtype=float)

susceptibility_ij_ccs = numpy.array([
    3.17348518,  1.65365136,  2.92133142,
    0.36680879,  2.33010953,  4.15557402,
    -0.68257245,  6.58982797,  6.19640529], dtype=float)

susceptibility_ij_direct = numpy.array([
    1.30155586,  4.91659845, 12.73686766,
   -0.51378685,  3.77335933,  9.83271568,
   -1.45326383,  4.48188431,  6.62508481], dtype=float)


eigen_val = numpy.array([1.04880579,  3.37869863, 10.67445929], dtype=float)

eigen_fields = numpy.array([
    [-0.32444153,  0.94414243, -0.05773004],
    [-0.69411522, -0.27909669, -0.66355791],
    [ 0.64260544,  0.17521444, -0.74589419]], dtype=float)

eigen_moments = numpy.array([
    [-0.30017151,  3.04655286, -3.459503  ],
    [ 0.93402198,  0.42411048, -4.66695702],
    [-0.3708013 , -1.39794508, -8.95519025]], dtype=float)


index_hkl = numpy.array([
    [0, 0, 2, 2, 1],
    [0, 2, 0, 2, 2],
    [2, 0, 0, 2, 3]], dtype=int)


symm_elems_r = numpy.array([
    [1,-1,-1, 1],
    [0, 0, 0, 0],
    [0, 0, 0, 0],
    [0, 0, 0, 0],
    [1,-1, 1,-1],
    [0, 0, 0, 0],
    [0, 0, 0, 0],
    [0, 0, 0, 0],
    [1,-1, 1,-1]], dtype=int)

susceptibility_2_ij = numpy.array([
    [3.1, 2.9, 5.7, 0.4, -0.7, 5.2],
    [1.5,-0.9, 1.7, 0.0,  0.9, 4.2]], dtype=float).transpose()

chi_s_direct = numpy.array(
[[[  1.30155586,   1.29548633],
  [  1.30155586,   1.29548633],
  [  1.30155586,   1.29548633],
  [  1.30155586,   1.29548633]],
 [[  4.91659845,  -1.33656834],
  [  4.91659845,  -1.33656834],
  [ -4.91659845,   1.33656834],
  [ -4.91659845,   1.33656834]],
 [[ 12.73686766,  14.4193239 ],
  [ 12.73686766,  14.4193239 ],
  [-12.73686766, -14.4193239 ],
  [-12.73686766, -14.4193239 ]],
 [[ -0.51378685,   0.2637808 ],
  [ -0.51378685,   0.2637808 ],
  [  0.51378685,  -0.2637808 ],
  [  0.51378685,  -0.2637808 ]],
 [[  3.77335933,  -1.281711  ],
  [  3.77335933,  -1.281711  ],
  [  3.77335933,  -1.281711  ],
  [  3.77335933,  -1.281711  ]],
 [[  9.83271568,   9.08409822],
  [  9.83271568,   9.08409822],
  [  9.83271568,   9.08409822],
  [  9.83271568,   9.08409822]],
 [[ -1.45326383,  -0.51041841],
  [ -1.45326383,  -0.51041841],
  [  1.45326383,   0.51041841],
  [  1.45326383,   0.51041841]],
 [[  4.48188431,   2.73703601],
  [  4.48188431,   2.73703601],
  [  4.48188431,   2.73703601],
  [  4.48188431,   2.73703601]],
 [[  6.62508481,   2.28622467],
  [  6.62508481,   2.28622467],
  [  6.62508481,   2.28622467],
  [  6.62508481,   2.28622467]]], dtype=float)

def test_calc_magnetization_ellipsoid_as_u():
    res = calc_magnetization_ellipsoid_as_u(
        susceptibility_ij, unit_cell_parameters)[0]

    print("res: ", res)
    print(reciprocal_u_norm)
    assert math.isclose(res[0], reciprocal_u_norm[0], abs_tol=abs_tol)
    assert math.isclose(res[1], reciprocal_u_norm[1], abs_tol=abs_tol)
    assert math.isclose(res[2], reciprocal_u_norm[2], abs_tol=abs_tol)
    assert math.isclose(res[3], reciprocal_u_norm[3], abs_tol=abs_tol)
    assert math.isclose(res[4], reciprocal_u_norm[4], abs_tol=abs_tol)
    assert math.isclose(res[5], reciprocal_u_norm[5], abs_tol=abs_tol)


def test_calc_magnetization_ellipsoid_axes():
    res_1, res_2, res_3 = calc_magnetization_ellipsoid_axes(
        susceptibility_ij, unit_cell_parameters)

    print("res_1: ", res_1)
    print(eigen_val)

    print("res_2: ", res_2)
    print(eigen_fields)

    print("res_3: ", res_3)
    print(eigen_moments)

    assert math.isclose(res_1[0], eigen_val[0], abs_tol=abs_tol)
    assert math.isclose(res_1[1], eigen_val[1], abs_tol=abs_tol)
    assert math.isclose(res_1[2], eigen_val[2], abs_tol=abs_tol)

    assert math.isclose(res_2[0][0], eigen_fields[0][0], abs_tol=abs_tol)
    assert math.isclose(res_2[0][1], eigen_fields[0][1], abs_tol=abs_tol)
    assert math.isclose(res_2[0][2], eigen_fields[0][2], abs_tol=abs_tol)
    assert math.isclose(res_2[1][0], eigen_fields[1][0], abs_tol=abs_tol)
    assert math.isclose(res_2[1][1], eigen_fields[1][1], abs_tol=abs_tol)
    assert math.isclose(res_2[1][2], eigen_fields[1][2], abs_tol=abs_tol)
    assert math.isclose(res_2[2][0], eigen_fields[2][0], abs_tol=abs_tol)
    assert math.isclose(res_2[2][1], eigen_fields[2][1], abs_tol=abs_tol)
    assert math.isclose(res_2[2][2], eigen_fields[2][2], abs_tol=abs_tol)

    assert math.isclose(res_3[0][0], eigen_moments[0][0], abs_tol=abs_tol)
    assert math.isclose(res_3[0][1], eigen_moments[0][1], abs_tol=abs_tol)
    assert math.isclose(res_3[0][2], eigen_moments[0][2], abs_tol=abs_tol)
    assert math.isclose(res_3[1][0], eigen_moments[1][0], abs_tol=abs_tol)
    assert math.isclose(res_3[1][1], eigen_moments[1][1], abs_tol=abs_tol)
    assert math.isclose(res_3[1][2], eigen_moments[1][2], abs_tol=abs_tol)
    assert math.isclose(res_3[2][0], eigen_moments[2][0], abs_tol=abs_tol)
    assert math.isclose(res_3[2][1], eigen_moments[2][1], abs_tol=abs_tol)
    assert math.isclose(res_3[2][2], eigen_moments[2][2], abs_tol=abs_tol)


def test_calc_ellipsoid_factor():
    res = calc_ellipsoid_factor(
        susceptibility_ij, unit_cell_parameters)

    print("res: ", res)
    print(ellipsoid_factor)
    assert math.isclose(res, ellipsoid_factor, abs_tol=abs_tol)


def test_calc_chi_direct_norm_with_symmetry():
    res = calc_chi_direct_norm_with_symmetry(
        susceptibility_2_ij[:, na, :], symm_elems_r[:, :, na], unit_cell_parameters)[0]

    print("res: ", res)
    print(chi_s_direct)
    assert numpy.all(numpy.isclose(res, chi_s_direct))
