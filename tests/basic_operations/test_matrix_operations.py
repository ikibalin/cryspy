import math
import numpy
import numpy.linalg

from cryspy.A_functions_base.matrix_operations import \
    calc_m_q_mt, calc_m1_m2, calc_q1_q2_q1, calc_m_q_inv_m, calc_mt_m, \
    calc_det_m, calc_inv_m, calc_m1_m2_inv_m1, calc_m1_m2_m1t, \
    calc_vt_m_v

m_11, m_12, m_13 =  1., -1., 5.
m_21, m_22, m_23 =  17., 8., 1.2
m_31, m_32, m_33 =  -11.2, -19., 5.7

m2_11, m2_12, m2_13 =  -1.7, -11.5, 5.44
m2_21, m2_22, m2_23 =  7.7, -7.3, 7.0
m2_31, m2_32, m2_33 =  -1.3, -0.7, -7.1

q_11, q_22, q_33 = 8.74, 6.7, 1.7
q_12, q_13, q_23 = 0.74, 7.7, -6.7

q2_11, q2_22, q2_33 = 1.74, 16.7, 23.7
q2_12, q2_13, q2_23 = 8.74, 7.7, -6.7

v_1, v_2, v_3 = 0.74, 7.7, -6.7

m_ij = numpy.stack([m_11, m_12, m_13, m_21, m_22, m_23, m_31, m_32, m_33], axis=0)
m2_ij = numpy.stack([m2_11, m2_12, m2_13, m2_21, m2_22, m2_23, m2_31, m2_32, m2_33], axis=0)
q_ij = numpy.stack([q_11, q_22, q_33, q_12, q_13, q_23], axis=0)
q2_ij = numpy.stack([q2_11, q2_22, q2_33, q2_12, q2_13, q2_23], axis=0)

v_i = numpy.stack([v_1, v_2, v_3], axis=0)

m_np = numpy.array([[m_11, m_12, m_13], [m_21, m_22, m_23], [m_31, m_32, m_33]], dtype=float)
m2_np = numpy.array([[m2_11, m2_12, m2_13], [m2_21, m2_22, m2_23], [m2_31, m2_32, m2_33]], dtype=float)
q_np = numpy.array([[q_11, q_12, q_13], [q_12, q_22, q_23], [q_13, q_23, q_33]], dtype=float)
q2_np = numpy.array([[q2_11, q2_12, q2_13], [q2_12, q2_22, q2_23], [q2_13, q2_23, q2_33]], dtype=float)


def test_calc_m_q_mt():
    o_ij = calc_m_q_mt(m_ij, q_ij)[0]
    o_np = numpy.matmul(numpy.matmul(m_np, q_np), numpy.transpose(m_np))

    print("o_np: ", o_np)
    print("o_ij: ", o_ij)
    assert math.isclose(o_ij[0], o_np[0, 0])
    assert math.isclose(o_ij[3], o_np[0, 1])
    assert math.isclose(o_ij[4], o_np[0, 2])
    assert math.isclose(o_ij[3], o_np[1, 0])
    assert math.isclose(o_ij[1], o_np[1, 1])
    assert math.isclose(o_ij[5], o_np[1, 2])
    assert math.isclose(o_ij[4], o_np[2, 0])
    assert math.isclose(o_ij[5], o_np[2, 1])
    assert math.isclose(o_ij[2], o_np[2, 2])


def test_calc_m1_m2():
    o_ij = calc_m1_m2(m_ij, m2_ij)[0]
    o_np = numpy.matmul(m_np, m2_np)

    print("o_np: ", o_np)
    print("o_ij: ", o_ij)
    assert math.isclose(o_ij[0], o_np[0, 0])
    assert math.isclose(o_ij[1], o_np[0, 1])
    assert math.isclose(o_ij[2], o_np[0, 2])
    assert math.isclose(o_ij[3], o_np[1, 0])
    assert math.isclose(o_ij[4], o_np[1, 1])
    assert math.isclose(o_ij[5], o_np[1, 2])
    assert math.isclose(o_ij[6], o_np[2, 0])
    assert math.isclose(o_ij[7], o_np[2, 1])
    assert math.isclose(o_ij[8], o_np[2, 2])


def test_calc_q1_q2_q1():
    o_ij = calc_q1_q2_q1(q_ij, q2_ij)[0]
    o_np = numpy.matmul(numpy.matmul(q_np, q2_np), q_np)

    print("o_np: ", o_np)
    print("o_ij: ", o_ij)
    assert math.isclose(o_ij[0], o_np[0, 0])
    assert math.isclose(o_ij[3], o_np[0, 1])
    assert math.isclose(o_ij[4], o_np[0, 2])
    assert math.isclose(o_ij[3], o_np[1, 0])
    assert math.isclose(o_ij[1], o_np[1, 1])
    assert math.isclose(o_ij[5], o_np[1, 2])
    assert math.isclose(o_ij[4], o_np[2, 0])
    assert math.isclose(o_ij[5], o_np[2, 1])
    assert math.isclose(o_ij[2], o_np[2, 2])


def test_calc_m_q_inv_m():
    o_ij = calc_m_q_inv_m(m_ij, q_ij)[0]
    inv_m_np = numpy.linalg.inv(m_np)
    o_np = numpy.matmul(numpy.matmul(m_np, q_np), inv_m_np)

    print("o_np: ", o_np)
    print("o_ij: ", o_ij)
    assert math.isclose(o_ij[0], o_np[0, 0])
    assert math.isclose(o_ij[1], o_np[0, 1])
    assert math.isclose(o_ij[2], o_np[0, 2])
    assert math.isclose(o_ij[3], o_np[1, 0])
    assert math.isclose(o_ij[4], o_np[1, 1])
    assert math.isclose(o_ij[5], o_np[1, 2])
    assert math.isclose(o_ij[6], o_np[2, 0])
    assert math.isclose(o_ij[7], o_np[2, 1])
    assert math.isclose(o_ij[8], o_np[2, 2])


def test_calc_mt_m():
    o_ij = calc_mt_m(m_ij)[0]
    o_np = numpy.matmul(numpy.transpose(m_np), m_np)

    print("o_np: ", o_np)
    print("o_ij: ", o_ij)
    assert math.isclose(o_ij[0], o_np[0, 0])
    assert math.isclose(o_ij[3], o_np[0, 1])
    assert math.isclose(o_ij[4], o_np[0, 2])
    assert math.isclose(o_ij[3], o_np[1, 0])
    assert math.isclose(o_ij[1], o_np[1, 1])
    assert math.isclose(o_ij[5], o_np[1, 2])
    assert math.isclose(o_ij[4], o_np[2, 0])
    assert math.isclose(o_ij[5], o_np[2, 1])
    assert math.isclose(o_ij[2], o_np[2, 2])

def test_calc_det_m():
    res = calc_det_m(m_ij)[0]
    det = numpy.linalg.det(m_np)

    print("res: ", res)
    print("det: ", det)
    assert math.isclose(res, det)

def test_calc_inv_m():
    o_ij = calc_inv_m(m_ij)[0]
    o_np = numpy.linalg.inv(m_np)

    print("o_np: ", o_np)
    print("o_ij: ", o_ij)
    assert math.isclose(o_ij[0], o_np[0, 0])
    assert math.isclose(o_ij[1], o_np[0, 1])
    assert math.isclose(o_ij[2], o_np[0, 2])
    assert math.isclose(o_ij[3], o_np[1, 0])
    assert math.isclose(o_ij[4], o_np[1, 1])
    assert math.isclose(o_ij[5], o_np[1, 2])
    assert math.isclose(o_ij[6], o_np[2, 0])
    assert math.isclose(o_ij[7], o_np[2, 1])
    assert math.isclose(o_ij[8], o_np[2, 2])


def test_calc_m1_m2_inv_m1():
    o_ij = calc_m1_m2_inv_m1(m_ij, m2_ij)[0]
    inv_m_np = numpy.linalg.inv(m_np)
    o_np = numpy.matmul(numpy.matmul(m_np, m2_np), inv_m_np)

    print("o_np: ", o_np)
    print("o_ij: ", o_ij)
    assert math.isclose(o_ij[0], o_np[0, 0])
    assert math.isclose(o_ij[1], o_np[0, 1])
    assert math.isclose(o_ij[2], o_np[0, 2])
    assert math.isclose(o_ij[3], o_np[1, 0])
    assert math.isclose(o_ij[4], o_np[1, 1])
    assert math.isclose(o_ij[5], o_np[1, 2])
    assert math.isclose(o_ij[6], o_np[2, 0])
    assert math.isclose(o_ij[7], o_np[2, 1])
    assert math.isclose(o_ij[8], o_np[2, 2])

def test_calc_m1_m2_m1t():
    o_ij = calc_m1_m2_m1t(m_ij, m2_ij)[0]
    o_np = numpy.matmul(numpy.matmul(m_np, m2_np), numpy.transpose(m_np))

    print("o_np: ", o_np)
    print("o_ij: ", o_ij)
    assert math.isclose(o_ij[0], o_np[0, 0])
    assert math.isclose(o_ij[1], o_np[0, 1])
    assert math.isclose(o_ij[2], o_np[0, 2])
    assert math.isclose(o_ij[3], o_np[1, 0])
    assert math.isclose(o_ij[4], o_np[1, 1])
    assert math.isclose(o_ij[5], o_np[1, 2])
    assert math.isclose(o_ij[6], o_np[2, 0])
    assert math.isclose(o_ij[7], o_np[2, 1])
    assert math.isclose(o_ij[8], o_np[2, 2])


def test_calc_vt_m_v():
    o_ij = calc_vt_m_v(m_ij, v_i)[0]
    o_np = (v_i*m_np.dot(v_i)).sum()

    print("o_np: ", o_np)
    print("o_ij: ", o_ij)
    assert math.isclose(o_ij, o_np)

