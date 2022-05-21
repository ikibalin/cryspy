# -*- coding: utf-8 -*-
"""
Expressions for matrix operations expressed as a full matrix (9 elements) or as a quadratic form (6 elements).

Order of 9 element matrix is 11, 12, 13, 21, 22, 23, 31, 32, 33
Order of 6 element quadratic form is 11, 22, 33, 12, 13, 23

For details see documentation.
"""

import numpy 
import numpy.linalg
na = numpy.newaxis
np_ol = numpy.ones_like
np_zl = numpy.zeros_like


def transpose_m(m, flag_m: bool = False):
    """Calculate determinant of the matrix.
    """
    m_11, m_12, m_13 = m[0], m[1], m[2]
    m_21, m_22, m_23 = m[3], m[4], m[5]
    m_31, m_32, m_33 = m[6], m[7], m[8]

    res = numpy.stack([m_11, m_21, m_31, m_12, m_22, m_32, m_13, m_23, m_33], axis=0)
    dder = {}
    if flag_m:
        # FIXME correct derivatives
        if m.dtype == complex:
            dder["m_real"] = numpy.zeros(m.shape[:1] + m.shape, dtype=real) 
            dder["m_imag"] = numpy.zeros(m.shape[:1] + m.shape, dtype=real)
        else:
            dder["m"] = numpy.zeros(m.shape[:1] + m.shape, dtype=m.dtype)
    return res, dder

def calc_norm_v(v, flag_v: bool = False):
    norm = numpy.sqrt(numpy.square(numpy.abs(v)).sum(axis=0))
    dder = {}
    if flag_v:
        flag_0 = numpy.isclose(norm, 0.)
        norm_inv = 1./numpy.where(flag_0, 1., norm)
        if v.dtype == complex:
            dder["v_real"] = numpy.where(flag_0, np_ol(v.real), v.real * norm_inv)
            dder["v_imag"] = numpy.where(flag_0, np_ol(v.real), v.imag * norm_inv)
        else:
            dder["v"] = numpy.where(flag_0, np_ol(v), v * norm_inv)
    return norm, dder

def calc_unity_v(v, flag_v: bool = False):
    norm, dder_norm = calc_norm_v(v, flag_v = flag_v)
    flag_0 = numpy.isclose(norm, 0.)
    norm_inv = 1./numpy.where(flag_0, 1., norm)
    unity = v*norm_inv
    dder = {}
    if flag_v:
        shape_v = v.shape
        shape_dder = v.shape[:1]+shape_v
        delta_ij = numpy.zeros(shape_dder, dtype=float)
        delta_ij[range(shape_v[0]), range(shape_v[0])] = 1. # numpy.ones(shape_v[1:], dtype=float)
        if v.dtype == complex:
            dder["v_real"] = delta_ij/norm - numpy.expand_dims(dder_norm["v_real"], axis=0)/numpy.square(norm) 
            dder["v_imag"] = 1j*delta_ij/norm - numpy.expand_dims(dder_norm["v_imag"], axis=0)/numpy.square(norm)
        else:
            dder["v"] = delta_ij/norm - numpy.expand_dims(dder_norm["v"], axis=0)/numpy.square(norm)
    return unity, dder


def calc_vv_as_v1_v2_v1(v1, flag_v1: bool = False):
    v1_1, v1_2, v1_3 = v1[0], v1[1], v1[2]

    v1_1_sq = numpy.square(v1_1)
    v1_2_sq = numpy.square(v1_2)
    v1_3_sq = numpy.square(v1_3)

    vv = numpy.stack([
        numpy.stack([(v1_3_sq  + v1_2_sq), - v1_1 * v1_2, - v1_1 * v1_3], axis=0),
        numpy.stack([- v1_2 * v1_1, (v1_1_sq  + v1_3_sq), - v1_2 * v1_3], axis=0),
        numpy.stack([- v1_3 * v1_1, - v1_3 * v1_2, (v1_2_sq  + v1_1_sq)], axis=0)
    ], axis=0)
    dder = {}
    if flag_v1:
        zeros = numpy.zeros_like(v1_1)
        dd_11_123 = numpy.stack([zeros, 2.*v1_2, 2*v1_3], axis=0)
        dd_12_123 = numpy.stack([- v1_2, - v1_1, zeros], axis=0)
        dd_13_123 = numpy.stack([- v1_3, zeros, - v1_1], axis=0)

        dd_21_123 = numpy.stack([- v1_2, - v1_1, zeros], axis=0)
        dd_22_123 = numpy.stack([2.*v1_1, zeros, 2*v1_3], axis=0)
        dd_23_123 = numpy.stack([zeros, - v1_3, - v1_2], axis=0)

        dd_31_123 = numpy.stack([- v1_3, zeros, - v1_1], axis=0)
        dd_32_123 = numpy.stack([zeros, - v1_3, - v1_2], axis=0)
        dd_33_123 = numpy.stack([2.*v1_1, 2*v1_2, zeros], axis=0)

        if v1.dtype == complex:
            dder["v1_real"] = numpy.stack([
                numpy.stack([dd_11_123, dd_12_123, dd_13_123], axis=0),
                numpy.stack([dd_21_123, dd_22_123, dd_23_123], axis=0),
                numpy.stack([dd_31_123, dd_32_123, dd_33_123], axis=0)
            ], axis=0)
            dder["v1_imag"] = numpy.stack([
                numpy.stack([dd_11_123*1j, dd_12_123*1j, dd_13_123*1j], axis=0),
                numpy.stack([dd_21_123*1j, dd_22_123*1j, dd_23_123*1j], axis=0),
                numpy.stack([dd_31_123*1j, dd_32_123*1j, dd_33_123*1j], axis=0)
            ], axis=0)
        else:
            dder["v1"] = numpy.stack([
                numpy.stack([dd_11_123, dd_12_123, dd_13_123], axis=0),
                numpy.stack([dd_21_123, dd_22_123, dd_23_123], axis=0),
                numpy.stack([dd_31_123, dd_32_123, dd_33_123], axis=0)
            ], axis=0)
    return vv, dder


def calc_vector_product_v1_v2_v1(v1, v2, flag_v1: bool = False, flag_v2: bool = False):
    """Calculated vector product v1 x v2 x v1.
    """
    vv, dder_vv = calc_vv_as_v1_v2_v1(v1, flag_v1=flag_v1)
    res = (vv*numpy.expand_dims(v2, axis=0)).sum(axis=1)

    # v1_1, v1_2, v1_3 = v1[0], v1[1], v1[2]
    # v2_1, v2_2, v2_3 = v2[0], v2[1], v2[2]

    dder = {}
    if flag_v1:
        if v1.dtype == complex:
            dder["v1_real"] = dder_vv["v1_real"]*numpy.expand_dims(numpy.expand_dims(v2, axis=0), axis=2).sum(axis=1)
            dder["v1_imag"] = dder_vv["v1_imag"]*numpy.expand_dims(numpy.expand_dims(v2, axis=0), axis=2).sum(axis=1)
        else:
            dder["v1"] = dder_vv["v1"]*numpy.expand_dims(numpy.expand_dims(v2, axis=0), axis=2).sum(axis=1)

    if flag_v2:
        d_2 = vv*numpy.ones((1, )+ v2.shape, dtype=float)
        if v2.dtype == complex:
            dder["v2_real"] = d_2
            dder["v2_imag"] = d_2*1j
        else:
            dder["v2"] = d_2
    return res, dder


def calc_vector_product_v1_v2(v1, v2, flag_v1: bool=False, flag_v2: bool=False):
    v1_1, v1_2, v1_3 = v1[0], v1[1], v1[2]
    v2_1, v2_2, v2_3 = v2[0], v2[1], v2[2]

    o_1 = v1_2*v2_3 - v1_3*v2_2
    o_2 = v1_3*v2_1 - v1_1*v2_3 
    o_3 = v1_1*v2_2 - v1_2*v2_1
    res = numpy.stack([o_1, o_2, o_3], axis=0)
    dder = {}
    if flag_v1:
        if v1.dtype == complex:
            dder["v1_real"] = None
            dder["v1_imag"] = None
        else:
            dder["v1"] = None
    if flag_v2:
        if v2.dtype == complex:
            dder["v2_real"] = None
            dder["v2_imag"] = None
        else:
            dder["v2"] = None
    return res, dder


def calc_m_q_mt(m_ij, q_ij, flag_m: bool = False, flag_q: bool = False):
    """
    q is quadratic form q_11, q_22, q_33, q_12, q_13, q_23
    m is matrix m_11, m_12, m_13, m_21, m_22, m_23, m_31, m_32, m_33

    Output is quadratic form o = M * Q * MT:
    o is o_11, o_22, o_33, o_12, o_13, o_23 
    """
    m_11, m_12, m_13 = m_ij[0], m_ij[1], m_ij[2]
    m_21, m_22, m_23 = m_ij[3], m_ij[4], m_ij[5]
    m_31, m_32, m_33 = m_ij[6], m_ij[7], m_ij[8]

    q_11, q_22, q_33 = q_ij[0], q_ij[1], q_ij[2]
    q_12, q_13, q_23 = q_ij[3], q_ij[4], q_ij[5]

    o_11 = m_11*(m_11*q_11 + m_12*q_12 + m_13*q_13) + m_12*(m_11*q_12 + m_12*q_22 + m_13*q_23) + m_13*(m_11*q_13 + m_12*q_23 + m_13*q_33)
    o_22 = m_21*(m_21*q_11 + m_22*q_12 + m_23*q_13) + m_22*(m_21*q_12 + m_22*q_22 + m_23*q_23) + m_23*(m_21*q_13 + m_22*q_23 + m_23*q_33)
    o_33 = m_31*(m_31*q_11 + m_32*q_12 + m_33*q_13) + m_32*(m_31*q_12 + m_32*q_22 + m_33*q_23) + m_33*(m_31*q_13 + m_32*q_23 + m_33*q_33)
    o_12 = m_21*(m_11*q_11 + m_12*q_12 + m_13*q_13) + m_22*(m_11*q_12 + m_12*q_22 + m_13*q_23) + m_23*(m_11*q_13 + m_12*q_23 + m_13*q_33)
    o_13 = m_31*(m_11*q_11 + m_12*q_12 + m_13*q_13) + m_32*(m_11*q_12 + m_12*q_22 + m_13*q_23) + m_33*(m_11*q_13 + m_12*q_23 + m_13*q_33)
    o_23 = m_31*(m_21*q_11 + m_22*q_12 + m_23*q_13) + m_32*(m_21*q_12 + m_22*q_22 + m_23*q_23) + m_33*(m_21*q_13 + m_22*q_23 + m_23*q_33)
    o_ij = numpy.stack([o_11, o_22, o_33, o_12, o_13, o_23], axis=0)
    dder = {}
    if flag_m:
        dder["m"] = numpy.zeros((6, )+ m_ij.shape, dtype=m_ij.dtype)
    if flag_q:
        dder["q"] = numpy.zeros((6, )+ q_ij.shape, dtype=q_ij.dtype)
    return o_ij, dder


def calc_mt_q_m(m, q, flag_m: bool = False, flag_q: bool = False):
    """
    q is quadratic form q_11, q_22, q_33, q_12, q_13, q_23
    m is matrix m_11, m_12, m_13, m_21, m_22, m_23, m_31, m_32, m_33

    Output is quadratic form o = MT * Q * M:
    o is o_11, o_22, o_33, o_12, o_13, o_23 
    """
    m_11, m_12, m_13 = m[0], m[1], m[2]
    m_21, m_22, m_23 = m[3], m[4], m[5]
    m_31, m_32, m_33 = m[6], m[7], m[8]

    q_11, q_22, q_33 = q[0], q[1], q[2]
    q_12, q_13, q_23 = q[3], q[4], q[5]

    qo_11 = m_11*(m_11*q_11 + m_21*q_12 + m_31*q_13) + m_21*(m_11*q_12 + m_21*q_22 + m_31*q_23) + m_31*(m_11*q_13 + m_21*q_23 + m_31*q_33)
    qo_22 = m_12*(m_12*q_11 + m_22*q_12 + m_32*q_13) + m_22*(m_12*q_12 + m_22*q_22 + m_32*q_23) + m_32*(m_12*q_13 + m_22*q_23 + m_32*q_33)
    qo_33 = m_13*(m_13*q_11 + m_23*q_12 + m_33*q_13) + m_23*(m_13*q_12 + m_23*q_22 + m_33*q_23) + m_33*(m_13*q_13 + m_23*q_23 + m_33*q_33)
    qo_12 = m_12*(m_11*q_11 + m_21*q_12 + m_31*q_13) + m_22*(m_11*q_12 + m_21*q_22 + m_31*q_23) + m_32*(m_11*q_13 + m_21*q_23 + m_31*q_33)
    qo_13 = m_13*(m_11*q_11 + m_21*q_12 + m_31*q_13) + m_23*(m_11*q_12 + m_21*q_22 + m_31*q_23) + m_33*(m_11*q_13 + m_21*q_23 + m_31*q_33)
    qo_23 = m_13*(m_12*q_11 + m_22*q_12 + m_32*q_13) + m_23*(m_12*q_12 + m_22*q_22 + m_32*q_23) + m_33*(m_12*q_13 + m_22*q_23 + m_32*q_33)

    qo = numpy.stack([qo_11, qo_22, qo_33, qo_12, qo_13, qo_23], axis=0)
    dder = {}
    if flag_m:
        zero = numpy.zeros((m_11*q_11).shape, dtype=float)
        hh = numpy.stack([
            numpy.stack([2*(m_11*q_11 + m_21*q_12 + m_31*q_13), zero, zero, 2*(m_11*q_12 + m_21*q_22 + m_31*q_23), zero, zero, 2*(m_11*q_13 + m_21*q_23 + m_31*q_33), zero, zero], axis=0),
            numpy.stack([zero, 2*(m_12*q_11 + m_22*q_12 + m_32*q_13), zero, zero, 2*(m_12*q_12 + m_22*q_22 + m_32*q_23), zero, zero, 2*(m_12*q_13 + m_22*q_23 + m_32*q_33), zero], axis=0),
            numpy.stack([zero, zero, 2*(m_13*q_11 + m_23*q_12 + m_33*q_13), zero, zero, 2*(m_13*q_12 + m_23*q_22 + m_33*q_23), zero, zero, 2*(m_13*q_13 + m_23*q_23 + m_33*q_33)], axis=0),
            numpy.stack([m_12*q_11+m_22*q_12+m_32*q_13, (m_11*q_11 + m_21*q_12 + m_31*q_13), zero, m_12*q_12+m_22*q_22+m_32*q_23, (m_11*q_12 + m_21*q_22 + m_31*q_23), zero, m_12*q_13+m_22*q_23+m_32*q_33, (m_11*q_13 + m_21*q_23 + m_31*q_33), zero], axis=0),
            numpy.stack([m_13*q_11+m_23*q_12+m_33*q_13, zero, (m_11*q_11 + m_21*q_12 + m_31*q_13), m_13*q_12+m_23*q_22+m_33*q_23, zero, (m_11*q_12 + m_21*q_22 + m_31*q_23), m_13*q_13+m_23*q_23+m_33*q_33, zero, (m_11*q_13 + m_21*q_23 + m_31*q_33)], axis=0),
            numpy.stack([zero, m_13*q_11+m_23*q_12+m_33*q_13, (m_12*q_11 + m_22*q_12 + m_32*q_13), zero, m_13*q_12+m_23*q_22+m_33*q_23, (m_12*q_12 + m_22*q_22 + m_32*q_23), zero, m_13*q_13+m_23*q_23+m_33*q_33, (m_12*q_13 + m_22*q_23 + m_32*q_33)], axis=0)
            ], axis=0)
        if m.dtype == complex:
            dder["m_real"] = hh
            dder["m_imag"] = 1j*hh
        else:
            dder["m"] = hh
    if flag_q:
        oq = numpy.ones(q_11.shape, dtype=float)
        hh = numpy.stack([
            numpy.stack([m_11*m_11*oq, m_21*m_21*oq, m_31*m_31*oq, (m_11*m_21+m_21*m_11)*oq,  (m_11*m_31+m_31*m_11)*oq, (m_21*m_31+m_31*m_21)*oq], axis=0),
            numpy.stack([m_12*m_12*oq, m_22*m_22*oq, m_32*m_32*oq, (m_12*m_22+m_22*m_12)*oq,  (m_12*m_32+m_32*m_12)*oq, (m_22*m_32+m_32*m_22)*oq], axis=0),
            numpy.stack([m_13*m_13*oq, m_23*m_23*oq, m_33*m_33*oq, (m_13*m_23+m_23*m_13)*oq,  (m_13*m_33+m_33*m_13)*oq, (m_23*m_33+m_33*m_23)*oq], axis=0),
            numpy.stack([m_12*m_11*oq, m_22*m_21*oq, m_32*m_31*oq, (m_12*m_21+m_22*m_11)*oq,  (m_12*m_31+m_32*m_11)*oq, (m_22*m_31+m_32*m_21)*oq], axis=0),
            numpy.stack([m_13*m_11*oq, m_23*m_21*oq, m_33*m_31*oq, (m_13*m_21+m_23*m_11)*oq,  (m_13*m_31+m_33*m_11)*oq, (m_23*m_31+m_33*m_21)*oq], axis=0),
            numpy.stack([m_13*m_12*oq, m_23*m_22*oq, m_33*m_32*oq, (m_13*m_22+m_23*m_12)*oq,  (m_13*m_32+m_33*m_12)*oq, (m_23*m_32+m_33*m_22)*oq], axis=0)],
            axis=0)
        if q.dtype == complex:
            dder["q_real"] = hh
            dder["q_imag"] = 1j*hh
        else:
            dder["q"] = hh
    return qo, dder


def calc_m1_m2(m1_ij, m2_ij, flag_m1: bool = False, flag_m2: bool = False):
    """
    m1, m2 are matrix m_11, m_12, m_13, m_21, m_22, m_23, m_31, m_32, m_33

    Output is matrix o = M1 * M2
    """
    m1_11, m1_12, m1_13 = m1_ij[0], m1_ij[1], m1_ij[2]
    m1_21, m1_22, m1_23 = m1_ij[3], m1_ij[4], m1_ij[5]
    m1_31, m1_32, m1_33 = m1_ij[6], m1_ij[7], m1_ij[8]

    m2_11, m2_12, m2_13 = m2_ij[0], m2_ij[1], m2_ij[2]
    m2_21, m2_22, m2_23 = m2_ij[3], m2_ij[4], m2_ij[5]
    m2_31, m2_32, m2_33 = m2_ij[6], m2_ij[7], m2_ij[8]

    o_11 = m1_11*m2_11 + m1_12*m2_21 + m1_13*m2_31
    o_12 = m1_11*m2_12 + m1_12*m2_22 + m1_13*m2_32
    o_13 = m1_11*m2_13 + m1_12*m2_23 + m1_13*m2_33
    o_21 = m1_21*m2_11 + m1_22*m2_21 + m1_23*m2_31
    o_22 = m1_21*m2_12 + m1_22*m2_22 + m1_23*m2_32
    o_23 = m1_21*m2_13 + m1_22*m2_23 + m1_23*m2_33
    o_31 = m1_31*m2_11 + m1_32*m2_21 + m1_33*m2_31
    o_32 = m1_31*m2_12 + m1_32*m2_22 + m1_33*m2_32
    o_33 = m1_31*m2_13 + m1_32*m2_23 + m1_33*m2_33
    o_ij = numpy.stack([o_11, o_12, o_13, o_21, o_22, o_23, o_31, o_32, o_33], axis=0)
    dder = {}
    if flag_m1:
        dder["m1"] = numpy.zeros((9, )+ m1_ij.shape, dtype=m1_ij.dtype)
    if flag_m2:
        dder["m2"] = numpy.zeros((9, )+ m2_ij.shape, dtype=m2_ij.dtype)
    return o_ij, dder


def calc_q1_q2_q1(q1_ij, q2_ij, flag_q1: bool = False, flag_q2: bool = False):
    """
    q is quadratic form q_11, q_22, q_33, q_12, q_13, q_23
    m is matrix m_11, m_12, m_13, m_21, m_22, m_23, m_31, m_32, m_33

    Output is quadratic form o = M * Q * MT:
    o is o_11, o_22, o_33, o_12, o_13, o_23 
    """
    q1_11, q1_22, q1_33 = q1_ij[0], q1_ij[1], q1_ij[2]
    q1_12, q1_13, q1_23 = q1_ij[3], q1_ij[4], q1_ij[5]

    q2_11, q2_22, q2_33 = q2_ij[0], q2_ij[1], q2_ij[2]
    q2_12, q2_13, q2_23 = q2_ij[3], q2_ij[4], q2_ij[5]

    o_11 = q1_11*(q1_11*q2_11 + q1_12*q2_12 + q1_13*q2_13) + q1_12*(q1_11*q2_12 + q1_12*q2_22 + q1_13*q2_23) + q1_13*(q1_11*q2_13 + q1_12*q2_23 + q1_13*q2_33)
    o_22 = q1_12*(q1_12*q2_11 + q1_22*q2_12 + q1_23*q2_13) + q1_22*(q1_12*q2_12 + q1_22*q2_22 + q1_23*q2_23) + q1_23*(q1_12*q2_13 + q1_22*q2_23 + q1_23*q2_33)
    o_33 = q1_13*(q1_13*q2_11 + q1_23*q2_12 + q1_33*q2_13) + q1_23*(q1_13*q2_12 + q1_23*q2_22 + q1_33*q2_23) + q1_33*(q1_13*q2_13 + q1_23*q2_23 + q1_33*q2_33)
    o_12 = q1_12*(q1_11*q2_11 + q1_12*q2_12 + q1_13*q2_13) + q1_22*(q1_11*q2_12 + q1_12*q2_22 + q1_13*q2_23) + q1_23*(q1_11*q2_13 + q1_12*q2_23 + q1_13*q2_33)
    o_13 = q1_13*(q1_11*q2_11 + q1_12*q2_12 + q1_13*q2_13) + q1_23*(q1_11*q2_12 + q1_12*q2_22 + q1_13*q2_23) + q1_33*(q1_11*q2_13 + q1_12*q2_23 + q1_13*q2_33)
    o_23 = q1_13*(q1_12*q2_11 + q1_22*q2_12 + q1_23*q2_13) + q1_23*(q1_12*q2_12 + q1_22*q2_22 + q1_23*q2_23) + q1_33*(q1_12*q2_13 + q1_22*q2_23 + q1_23*q2_33)
    o_ij = numpy.stack([o_11, o_22, o_33, o_12, o_13, o_23], axis=0)
    dder = {}
    if flag_q1:
        dder["q1"] = numpy.zeros((6, )+ q1_ij.shape, dtype=q1_ij.dtype)
    if flag_q2:
        dder["q2"] = numpy.zeros((6, )+ q2_ij.shape, dtype=q1_ij.dtype)
    return o_ij, dder


def calc_m_q_inv_m(m, q, flag_m: bool = False, flag_q: bool = False):
    """
    q is quadratic form q_11, q_22, q_33, q_12, q_13, q_23
    m is matrix m_11, m_12, m_13, m_21, m_22, m_23, m_31, m_32, m_33

    Output is matrix  o
    """
    
    m_11, m_12, m_13 = m[0], m[1], m[2]
    m_21, m_22, m_23 = m[3], m[4], m[5]
    m_31, m_32, m_33 = m[6], m[7], m[8]
    
    q_11, q_22, q_33 = q[0], q[1], q[2]
    q_12, q_13, q_23 = q[3], q[4], q[5]

    det, dder_det = calc_det_m(m, flag_m=flag_m)
    
    o_11  = ((m_21*m_32 - m_22*m_31)*(q_13*m_11 + q_23*m_12 + q_33*m_13) - (m_21*m_33 - m_23*m_31)*(q_12*m_11 + q_22*m_12 + q_23*m_13) + (m_22*m_33 - m_23*m_32)*(q_11*m_11 + q_12*m_12 + q_13*m_13))/det
    o_12 = (-(m_11*m_32 - m_12*m_31)*(q_13*m_11 + q_23*m_12 + q_33*m_13) + (m_11*m_33 - m_13*m_31)*(q_12*m_11 + q_22*m_12 + q_23*m_13) - (m_12*m_33 - m_13*m_32)*(q_11*m_11 + q_12*m_12 + q_13*m_13))/det
    o_13 = ((m_11*m_22 - m_12*m_21)*(q_13*m_11 + q_23*m_12 + q_33*m_13) - (m_11*m_23 - m_13*m_21)*(q_12*m_11 + q_22*m_12 + q_23*m_13) + (m_12*m_23 - m_13*m_22)*(q_11*m_11 + q_12*m_12 + q_13*m_13))/det
    o_21 = ((m_21*m_32 - m_22*m_31)*(q_13*m_21 + q_23*m_22 + q_33*m_23) - (m_21*m_33 - m_23*m_31)*(q_12*m_21 + q_22*m_22 + q_23*m_23) + (m_22*m_33 - m_23*m_32)*(q_11*m_21 + q_12*m_22 + q_13*m_23))/det
    o_22 = (-(m_11*m_32 - m_12*m_31)*(q_13*m_21 + q_23*m_22 + q_33*m_23) + (m_11*m_33 - m_13*m_31)*(q_12*m_21 + q_22*m_22 + q_23*m_23) - (m_12*m_33 - m_13*m_32)*(q_11*m_21 + q_12*m_22 + q_13*m_23))/det
    o_23 = ((m_11*m_22 - m_12*m_21)*(q_13*m_21 + q_23*m_22 + q_33*m_23) - (m_11*m_23 - m_13*m_21)*(q_12*m_21 + q_22*m_22 + q_23*m_23) + (m_12*m_23 - m_13*m_22)*(q_11*m_21 + q_12*m_22 + q_13*m_23))/det
    o_31 = ((m_21*m_32 - m_22*m_31)*(q_13*m_31 + q_23*m_32 + q_33*m_33) - (m_21*m_33 - m_23*m_31)*(q_12*m_31 + q_22*m_32 + q_23*m_33) + (m_22*m_33 - m_23*m_32)*(q_11*m_31 + q_12*m_32 + q_13*m_33))/det
    o_32 = (-(m_11*m_32 - m_12*m_31)*(q_13*m_31 + q_23*m_32 + q_33*m_33) + (m_11*m_33 - m_13*m_31)*(q_12*m_31 + q_22*m_32 + q_23*m_33) - (m_12*m_33 - m_13*m_32)*(q_11*m_31 + q_12*m_32 + q_13*m_33))/det
    o_33 = ((m_11*m_22 - m_12*m_21)*(q_13*m_31 + q_23*m_32 + q_33*m_33) - (m_11*m_23 - m_13*m_21)*(q_12*m_31 + q_22*m_32 + q_23*m_33) + (m_12*m_23 - m_13*m_22)*(q_11*m_31 + q_12*m_32 + q_13*m_33))/det
    o_ij = numpy.stack([o_11, o_12, o_13, o_21, o_22, o_23, o_31, o_32, o_33], axis=0)
    dder = {}
    if flag_m:
        dder["m"] = numpy.zeros((9, )+ m.shape, dtype=m.dtype)
    if flag_q:
        oq = numpy.ones(q_11.shape, dtype=float)/det
        dder_q = numpy.stack([
            numpy.stack([+(m_22*m_33-m_23*m_32)*m_11*oq, -(m_21*m_33-m_23*m_31)*m_12*oq, +(m_21*m_32-m_22*m_31)*m_13*oq, -(m_21*m_33-m_23*m_31)*m_11*oq+(m_22*m_33-m_23*m_32)*m_12*oq, +(m_21*m_32-m_22*m_31)*m_11*oq+(m_22*m_33-m_23*m_32)*m_13*oq, +(m_21*m_32-m_22*m_31)*m_12-(m_21*m_33-m_23*m_31)*m_13*oq], axis=0),
            numpy.stack([-(m_12*m_33-m_13*m_32)*m_11*oq, +(m_11*m_33-m_13*m_31)*m_12*oq, -(m_11*m_32-m_12*m_31)*m_13*oq, +(m_11*m_33-m_13*m_31)*m_11*oq-(m_12*m_33-m_13*m_32)*m_12*oq, -(m_11*m_32-m_12*m_31)*m_11*oq-(m_12*m_33-m_13*m_32)*m_13*oq, -(m_11*m_32-m_12*m_31)*m_12+(m_11*m_33-m_13*m_31)*m_13*oq], axis=0),
            numpy.stack([+(m_12*m_23-m_13*m_22)*m_11*oq, -(m_11*m_23-m_13*m_21)*m_12*oq, +(m_11*m_22-m_12*m_21)*m_13*oq, -(m_11*m_23-m_13*m_21)*m_11*oq+(m_12*m_23-m_13*m_22)*m_12*oq, +(m_11*m_22-m_12*m_21)*m_11*oq+(m_12*m_23-m_13*m_22)*m_13*oq, +(m_11*m_22-m_12*m_21)*m_12-(m_11*m_23-m_13*m_21)*m_13*oq], axis=0),
            numpy.stack([+(m_22*m_33-m_23*m_32)*m_21*oq, -(m_21*m_33-m_23*m_31)*m_22*oq, +(m_21*m_32-m_22*m_31)*m_23*oq, -(m_21*m_33-m_23*m_31)*m_21*oq+(m_22*m_33-m_23*m_32)*m_22*oq, +(m_21*m_32-m_22*m_31)*m_21*oq+(m_22*m_33-m_23*m_32)*m_23*oq, +(m_21*m_32-m_22*m_31)*m_22-(m_21*m_33-m_23*m_31)*m_23*oq], axis=0),
            numpy.stack([-(m_12*m_33-m_13*m_32)*m_21*oq, +(m_11*m_33-m_13*m_31)*m_22*oq, -(m_11*m_32-m_12*m_31)*m_23*oq, +(m_11*m_33-m_13*m_31)*m_21*oq-(m_12*m_33-m_13*m_32)*m_22*oq, -(m_11*m_32-m_12*m_31)*m_21*oq-(m_12*m_33-m_13*m_32)*m_23*oq, -(m_11*m_32-m_12*m_31)*m_22+(m_11*m_33-m_13*m_31)*m_23*oq], axis=0),
            numpy.stack([+(m_12*m_23-m_13*m_22)*m_21*oq, -(m_11*m_23-m_13*m_21)*m_22*oq, +(m_11*m_22-m_12*m_21)*m_23*oq, -(m_11*m_23-m_13*m_21)*m_21*oq+(m_12*m_23-m_13*m_22)*m_22*oq, +(m_11*m_22-m_12*m_21)*m_21*oq+(m_12*m_23-m_13*m_22)*m_23*oq, +(m_11*m_22-m_12*m_21)*m_22-(m_11*m_23-m_13*m_21)*m_23*oq], axis=0),
            numpy.stack([+(m_22*m_33-m_23*m_32)*m_31*oq, -(m_21*m_33-m_23*m_31)*m_32*oq, +(m_21*m_32-m_22*m_31)*m_33*oq, -(m_21*m_33-m_23*m_31)*m_31*oq+(m_22*m_33-m_23*m_32)*m_32*oq, +(m_21*m_32-m_22*m_31)*m_31*oq+(m_22*m_33-m_23*m_32)*m_33*oq, +(m_21*m_32-m_22*m_31)*m_32-(m_21*m_33-m_23*m_31)*m_33*oq], axis=0),
            numpy.stack([-(m_12*m_33-m_13*m_32)*m_31*oq, +(m_11*m_33-m_13*m_31)*m_32*oq, -(m_11*m_32-m_12*m_31)*m_33*oq, +(m_11*m_33-m_13*m_31)*m_31*oq-(m_12*m_33-m_13*m_32)*m_32*oq, -(m_11*m_32-m_12*m_31)*m_31*oq-(m_12*m_33-m_13*m_32)*m_33*oq, -(m_11*m_32-m_12*m_31)*m_32+(m_11*m_33-m_13*m_31)*m_33*oq], axis=0),
            numpy.stack([+(m_12*m_23-m_13*m_22)*m_31*oq, -(m_11*m_23-m_13*m_21)*m_32*oq, +(m_11*m_22-m_12*m_21)*m_33*oq, -(m_11*m_23-m_13*m_21)*m_31*oq+(m_12*m_23-m_13*m_22)*m_32*oq, +(m_11*m_22-m_12*m_21)*m_31*oq+(m_12*m_23-m_13*m_22)*m_33*oq, +(m_11*m_22-m_12*m_21)*m_32-(m_11*m_23-m_13*m_21)*m_33*oq], axis=0)], axis=0)
        if q.dtype == float:
            dder["q"] = dder_q
        else:
            dder["q_real"] = dder_q
            dder["q_imag"] = 1j*dder_q

    return o_ij, dder


def calc_mt_m(m_ij, flag_m: bool = False):
    """
    m is matrix m_11, m_12, m_13, m_21, m_22, m_23, m_31, m_32, m_33

    Output is quadratic form q
    """
    m_11, m_12, m_13 = m_ij[0], m_ij[1], m_ij[2]
    m_21, m_22, m_23 = m_ij[3], m_ij[4], m_ij[5]
    m_31, m_32, m_33 = m_ij[6], m_ij[7], m_ij[8]

    q_11 = numpy.square(m_11) + numpy.square(m_21) + numpy.square(m_31) 
    q_22 = numpy.square(m_12) + numpy.square(m_22) + numpy.square(m_32) 
    q_33 = numpy.square(m_13) + numpy.square(m_23) + numpy.square(m_33) 
    q_12 = m_11 * m_12 + m_21 * m_22 + m_31 * m_32
    q_13 = m_11 * m_13 + m_21 * m_23 + m_31 * m_33
    q_23 = m_12 * m_13 + m_22 * m_23 + m_32 * m_33

    q_ij = numpy.stack([q_11, q_22, q_33, q_12, q_13, q_23], axis=0)
    dder = {}
    if flag_m:
        dder["m"] = numpy.zeros((6, )+ m_ij.shape, dtype=m_ij.dtype)
    return q_ij, dder



def calc_q_sq(q, flag_q: bool = False):
    """
    q is matrix q_11, q_22, q_33, q_12, q_13, q_23

    Output is quadratic form
    """
    q_11, q_22, q_33 = q[0], q[1], q[2]
    q_12, q_13, q_23 = q[3], q[4], q[5]

    qo_11 = numpy.square(q_11) + numpy.square(q_12) + numpy.square(q_13) 
    qo_22 = numpy.square(q_12) + numpy.square(q_22) + numpy.square(q_23) 
    qo_33 = numpy.square(q_13) + numpy.square(q_23) + numpy.square(q_33) 
    qo_12 = q_11 * q_12 + q_12 * q_22 + q_13 * q_23
    qo_13 = q_11 * q_13 + q_12 * q_23 + q_13 * q_33
    qo_23 = q_12 * q_13 + q_22 * q_23 + q_23 * q_33

    qo = numpy.stack([qo_11, qo_22, qo_33, qo_12, qo_13, qo_23], axis=0)
    dder = {}
    if flag_q:
        zero = numpy.zeros(q_11.shape, dtype=float)
        dder_11 = numpy.stack([2*q_11, zero, zero, 2*q_12, 2*q_13, zero], axis=0)
        dder_22 = numpy.stack([zero, 2*q_22, zero, 2*q_12, zero, 2*q_23], axis=0)
        dder_33 = numpy.stack([zero, zero, 2*q_33, zero, 2*q_13, 2*q_23], axis=0)
        dder_12 = numpy.stack([q_12, q_12, zero, q_11+q_22, q_23, q_13], axis=0)
        dder_13 = numpy.stack([q_13, zero, q_13, q_23, q_11+q_33, q_12], axis=0)
        dder_23 = numpy.stack([zero, q_23, q_23, q_13, q_12, q_22+q_33], axis=0)
        if q.dtype == complex:
            dder["q_real"] = numpy.stack([dder_11, dder_22, dder_33, dder_12, dder_13, dder_23], axis=0)
            dder["q_imag"] = 1j*dder["q_real"]
        else:
            dder["q"] = numpy.stack([dder_11, dder_22, dder_33, dder_12, dder_13, dder_23], axis=0)
    return qo, dder


def calc_m_sq(m, flag_m: bool = False):
    """
    m is matrix m_11, m_12, m_13, m_21, m_22, m_23, m_31, m_32, m_33

    Output is m_sq matrix
    """
    m_11, m_12, m_13 = m[0], m[1], m[2]
    m_21, m_22, m_23 = m[3], m[4], m[5]
    m_31, m_32, m_33 = m[6], m[7], m[8]

    mo_11 = numpy.square(m_11) + m_12*m_21 + m_13*m_31 
    mo_22 = m_12*m_21 + numpy.square(m_22) + m_23*m_32 
    mo_33 = m_13*m_31 + m_23*m_32 + numpy.square(m_33) 
    mo_12 = m_11 * m_12 + m_12 * m_22 + m_13 * m_32
    mo_13 = m_11 * m_13 + m_12 * m_23 + m_13 * m_33
    mo_23 = m_13 * m_21 + m_22 * m_23 + m_23 * m_33
    mo_21 = m_11 * m_21 + m_21 * m_22 + m_23 * m_31
    mo_31 = m_11 * m_31 + m_21 * m_32 + m_33 * m_33
    mo_32 = m_12 * m_31 + m_22 * m_32 + m_32 * m_33

    mo = numpy.stack([
        mo_11, mo_12, mo_13,
        mo_21, mo_22, mo_23,
        mo_31, mo_32, mo_33], axis=0)
    dder = {}
    if flag_m:
        zero = numpy.zeros(m_11.shape, dtype=float)
        dder_11 = numpy.stack([2*m_11 , m_21 , m_31 , m_12 , zero , zero , m_13 , zero , zero], axis=0)
        dder_12 = numpy.stack([m_12 , m_11 + m_22 , m_32 , zero , m_12 , zero , zero , m_13 , zero], axis=0)
        dder_13 = numpy.stack([m_13 , m_23 , m_11 + m_33 , zero , zero , m_12 , zero , zero , m_13], axis=0)
        dder_21 = numpy.stack([m_21 , zero , zero , m_11 + m_22 , m_21 , m_31 , m_23 , zero , zero], axis=0)
        dder_22 = numpy.stack([zero , m_21 , zero , m_12 , 2*m_22 , m_32 , zero , m_23 , zero], axis=0)
        dder_23 = numpy.stack([zero , zero , m_21 , m_13 , m_23 , m_22 + m_33 , zero , zero , m_23], axis=0)
        dder_31 = numpy.stack([m_31 , zero , zero , m_32 , zero , zero , m_11 + m_33 , m_21 , m_31], axis=0)
        dder_32 = numpy.stack([zero , m_31 , zero , zero , m_32 , zero , m_12 , m_22 + m_33 , m_32], axis=0)
        dder_33 = numpy.stack([zero , zero , m_31 , zero , zero , m_32 , m_13 , m_23 , 2*m_33], axis=0)
        if m.dtype == complex:
            dder["m_real"] = numpy.stack([dder_11, dder_12, dder_13, dder_21, dder_22, dder_23, dder_31, dder_32, dder_33], axis=0)
            dder["m_imag"] = 1j*dder["m_real"]
        else:
            dder["m"] = numpy.stack([dder_11, dder_12, dder_13, dder_21, dder_22, dder_23, dder_31, dder_32, dder_33], axis=0)
    return mo, dder


def calc_m1_m2_inv_m1(m1_ij, m2_ij, flag_m1: bool = False, flag_m2: bool = False):
    m1_11, m1_12, m1_13 = m1_ij[0], m1_ij[1], m1_ij[2]
    m1_21, m1_22, m1_23 = m1_ij[3], m1_ij[4], m1_ij[5]
    m1_31, m1_32, m1_33 = m1_ij[6], m1_ij[7], m1_ij[8]

    m2_11, m2_12, m2_13 = m2_ij[0], m2_ij[1], m2_ij[2]
    m2_21, m2_22, m2_23 = m2_ij[3], m2_ij[4], m2_ij[5]
    m2_31, m2_32, m2_33 = m2_ij[6], m2_ij[7], m2_ij[8]
    det = (m1_11*m1_22*m1_33 - m1_11*m1_23*m1_32 - m1_12*m1_21*m1_33 + m1_12*m1_23*m1_31 + m1_13*m1_21*m1_32 - m1_13*m1_22*m1_31)

    m1_m2_inv_m1_11 = ((m1_21*m1_32 - m1_22*m1_31)*(m1_11*m2_13 + m1_12*m2_23 + m1_13*m2_33) -
                       (m1_21*m1_33 - m1_23*m1_31)*(m1_11*m2_12 + m1_12*m2_22 + m1_13*m2_32) +
                       (m1_22*m1_33 - m1_23*m1_32)*(m1_11*m2_11 + m1_12*m2_21 + m1_13*m2_31))/det
    m1_m2_inv_m1_12 = (-(m1_11*m1_32 - m1_12*m1_31)*(m1_11*m2_13 + m1_12*m2_23 + m1_13*m2_33) +
                        (m1_11*m1_33 - m1_13*m1_31)*(m1_11*m2_12 + m1_12*m2_22 + m1_13*m2_32) -
                        (m1_12*m1_33 - m1_13*m1_32)*(m1_11*m2_11 + m1_12*m2_21 + m1_13*m2_31))/det
    m1_m2_inv_m1_13 = ((m1_11*m1_22 - m1_12*m1_21)*(m1_11*m2_13 + m1_12*m2_23 + m1_13*m2_33) -
                       (m1_11*m1_23 - m1_13*m1_21)*(m1_11*m2_12 + m1_12*m2_22 + m1_13*m2_32) +
                       (m1_12*m1_23 - m1_13*m1_22)*(m1_11*m2_11 + m1_12*m2_21 + m1_13*m2_31))/det
    m1_m2_inv_m1_21 = ((m1_21*m1_32 - m1_22*m1_31)*(m1_21*m2_13 + m1_22*m2_23 + m1_23*m2_33) -
                       (m1_21*m1_33 - m1_23*m1_31)*(m1_21*m2_12 + m1_22*m2_22 + m1_23*m2_32) +
                       (m1_22*m1_33 - m1_23*m1_32)*(m1_21*m2_11 + m1_22*m2_21 + m1_23*m2_31))/det
    m1_m2_inv_m1_22 = (-(m1_11*m1_32 - m1_12*m1_31)*(m1_21*m2_13 + m1_22*m2_23 + m1_23*m2_33) +
                       (m1_11*m1_33 - m1_13*m1_31)*(m1_21*m2_12 + m1_22*m2_22 + m1_23*m2_32) -
                       (m1_12*m1_33 - m1_13*m1_32)*(m1_21*m2_11 + m1_22*m2_21 + m1_23*m2_31))/det
    m1_m2_inv_m1_23 = ((m1_11*m1_22 - m1_12*m1_21)*(m1_21*m2_13 + m1_22*m2_23 + m1_23*m2_33) -
                       (m1_11*m1_23 - m1_13*m1_21)*(m1_21*m2_12 + m1_22*m2_22 + m1_23*m2_32) +
                       (m1_12*m1_23 - m1_13*m1_22)*(m1_21*m2_11 + m1_22*m2_21 + m1_23*m2_31))/det
    m1_m2_inv_m1_31 = ((m1_21*m1_32 - m1_22*m1_31)*(m1_31*m2_13 + m1_32*m2_23 + m1_33*m2_33) -
                       (m1_21*m1_33 - m1_23*m1_31)*(m1_31*m2_12 + m1_32*m2_22 + m1_33*m2_32) +
                       (m1_22*m1_33 - m1_23*m1_32)*(m1_31*m2_11 + m1_32*m2_21 + m1_33*m2_31))/det
    m1_m2_inv_m1_32 = (-(m1_11*m1_32 - m1_12*m1_31)*(m1_31*m2_13 + m1_32*m2_23 + m1_33*m2_33) +
                       (m1_11*m1_33 - m1_13*m1_31)*(m1_31*m2_12 + m1_32*m2_22 + m1_33*m2_32) -
                       (m1_12*m1_33 - m1_13*m1_32)*(m1_31*m2_11 + m1_32*m2_21 + m1_33*m2_31))/det
    m1_m2_inv_m1_33 = ((m1_11*m1_22 - m1_12*m1_21)*(m1_31*m2_13 + m1_32*m2_23 + m1_33*m2_33) -
                       (m1_11*m1_23 - m1_13*m1_21)*(m1_31*m2_12 + m1_32*m2_22 + m1_33*m2_32) +
                       (m1_12*m1_23 - m1_13*m1_22)*(m1_31*m2_11 + m1_32*m2_21 + m1_33*m2_31))/det

    m1_m2_inv_m1 = numpy.stack([
        m1_m2_inv_m1_11, m1_m2_inv_m1_12, m1_m2_inv_m1_13,
        m1_m2_inv_m1_21, m1_m2_inv_m1_22, m1_m2_inv_m1_23,
        m1_m2_inv_m1_31, m1_m2_inv_m1_32, m1_m2_inv_m1_33], axis=0)
    dder = {}
    if flag_m1:
        dder["m1"] = numpy.zeros((9, )+ m1_ij.shape, dtype=m1_ij.dtype)
    if flag_m2:
        dder["m2"] = numpy.zeros((9, )+ m2_ij.shape, dtype=m2_ij.dtype)
    return m1_m2_inv_m1, dder

def calc_m1_m2_m1t(m1_ij, m2_ij, flag_m1: bool = False, flag_m2: bool = False):
    m1_11, m1_12, m1_13 = m1_ij[0], m1_ij[1], m1_ij[2]
    m1_21, m1_22, m1_23 = m1_ij[3], m1_ij[4], m1_ij[5]
    m1_31, m1_32, m1_33 = m1_ij[6], m1_ij[7], m1_ij[8]

    m2_11, m2_12, m2_13 = m2_ij[0], m2_ij[1], m2_ij[2]
    m2_21, m2_22, m2_23 = m2_ij[3], m2_ij[4], m2_ij[5]
    m2_31, m2_32, m2_33 = m2_ij[6], m2_ij[7], m2_ij[8]

    m1_m2_m1t_11 = m1_11*(m1_11*m2_11 + m1_12*m2_21 + m1_13*m2_31) + \
        m1_12*(m1_11*m2_12 + m1_12*m2_22 + m1_13*m2_32) + \
        m1_13*(m1_11*m2_13 + m1_12*m2_23 + m1_13*m2_33)
    m1_m2_m1t_12 = m1_21*(m1_11*m2_11 + m1_12*m2_21 + m1_13*m2_31) + \
        m1_22*(m1_11*m2_12 + m1_12*m2_22 + m1_13*m2_32) + \
        m1_23*(m1_11*m2_13 + m1_12*m2_23 + m1_13*m2_33)
    m1_m2_m1t_13 = m1_31*(m1_11*m2_11 + m1_12*m2_21 + m1_13*m2_31) + \
        m1_32*(m1_11*m2_12 + m1_12*m2_22 + m1_13*m2_32) + \
        m1_33*(m1_11*m2_13 + m1_12*m2_23 + m1_13*m2_33)
    m1_m2_m1t_21 = m1_11*(m1_21*m2_11 + m1_22*m2_21 + m1_23*m2_31) + \
        m1_12*(m1_21*m2_12 + m1_22*m2_22 + m1_23*m2_32) + \
        m1_13*(m1_21*m2_13 + m1_22*m2_23 + m1_23*m2_33)
    m1_m2_m1t_22 = m1_21*(m1_21*m2_11 + m1_22*m2_21 + m1_23*m2_31) + \
        m1_22*(m1_21*m2_12 + m1_22*m2_22 + m1_23*m2_32) + \
        m1_23*(m1_21*m2_13 + m1_22*m2_23 + m1_23*m2_33)
    m1_m2_m1t_23 = m1_31*(m1_21*m2_11 + m1_22*m2_21 + m1_23*m2_31) + \
        m1_32*(m1_21*m2_12 + m1_22*m2_22 + m1_23*m2_32) + \
        m1_33*(m1_21*m2_13 + m1_22*m2_23 + m1_23*m2_33)
    m1_m2_m1t_31 = m1_11*(m1_31*m2_11 + m1_32*m2_21 + m1_33*m2_31) + \
        m1_12*(m1_31*m2_12 + m1_32*m2_22 + m1_33*m2_32) + \
        m1_13*(m1_31*m2_13 + m1_32*m2_23 + m1_33*m2_33)
    m1_m2_m1t_32 = m1_21*(m1_31*m2_11 + m1_32*m2_21 + m1_33*m2_31) + \
        m1_22*(m1_31*m2_12 + m1_32*m2_22 + m1_33*m2_32) + \
        m1_23*(m1_31*m2_13 + m1_32*m2_23 + m1_33*m2_33)
    m1_m2_m1t_33 = m1_31*(m1_31*m2_11 + m1_32*m2_21 + m1_33*m2_31) + \
        m1_32*(m1_31*m2_12 + m1_32*m2_22 + m1_33*m2_32) + \
        m1_33*(m1_31*m2_13 + m1_32*m2_23 + m1_33*m2_33)

    m1_m2_m1t = numpy.stack([
        m1_m2_m1t_11, m1_m2_m1t_12, m1_m2_m1t_13,
        m1_m2_m1t_21, m1_m2_m1t_22, m1_m2_m1t_23,
        m1_m2_m1t_31, m1_m2_m1t_32, m1_m2_m1t_33], axis=0)
    dder = {}
    if flag_m1:
        dder["m1"] = numpy.zeros((9, )+ m1_ij.shape, dtype=m1_ij.dtype)
    if flag_m2:
        dder["m2"] = numpy.zeros((9, )+ m2_ij.shape, dtype=m2_ij.dtype)
    return m1_m2_m1t, dder


def calc_vt_m_v(m_ij, v_i, flag_m: bool = False, flag_v: bool = False):
    m_11, m_12, m_13 = m_ij[0], m_ij[1], m_ij[2]
    m_21, m_22, m_23 = m_ij[3], m_ij[4], m_ij[5]
    m_31, m_32, m_33 = m_ij[6], m_ij[7], m_ij[8]
    v_1 = v_i[0]
    v_2 = v_i[1]
    v_3 = v_i[2]
    res = (m_11*v_1*v_1 + m_22*v_2*v_2 + m_33*v_3*v_3 +
           (m_12+m_21)*v_1*v_2 + (m_13+m_31)*v_1*v_3+(m_23+m_32)*v_2*v_3)
    dder = {}
    if flag_m:
        dder["m"] = numpy.stack([
            v_1*v_1*numpy.ones_like(m_11), v_1*v_2*numpy.ones_like(m_12),
            v_1*v_3*numpy.ones_like(m_13), v_2*v_1*numpy.ones_like(m_21),
            v_2*v_2*numpy.ones_like(m_22), v_2*v_3*numpy.ones_like(m_23),
            v_3*v_1*numpy.ones_like(m_31), v_3*v_2*numpy.ones_like(m_32),
            v_3*v_3*numpy.ones_like(m_33)], axis=0)
    if flag_v:
        dder["v"] = numpy.stack([
            (2*m_11*v_1+(m_12+m_21)*v_2+(m_13+m_31)*v_3)*numpy.ones_like(v_1),
            (2*m_22*v_2+(m_12+m_21)*v_1+(m_23+m_32)*v_3)*numpy.ones_like(v_2),
            (2*m_33*v_3+(m_13+m_31)*v_1+(m_23+m_32)*v_2)*numpy.ones_like(v_3)
        ], axis=0)
    return res, dder


def calc_vt_q_v(q, v, flag_q: bool = False, flag_v: bool = False):
    q_11, q_22, q_33 = q[0], q[1], q[2]
    q_12, q_13, q_23 = q[3], q[4], q[5]
    v_1 = v[0]
    v_2 = v[1]
    v_3 = v[2]
    res = (q_11*v_1*v_1 + q_22*v_2*v_2 + q_33*v_3*v_3 +
           2*q_12*v_1*v_2 + 2*q_13*v_1*v_3+2*q_23*v_2*v_3)
    dder = {}
    if flag_q:
        dder["q"] = numpy.stack([
            v_1*v_1*numpy.ones_like(q_11), v_1*v_2*numpy.ones_like(m_12),
            v_1*v_3*numpy.ones_like(q_13), v_2*v_1*numpy.ones_like(m_12),
            v_2*v_2*numpy.ones_like(q_22), v_2*v_3*numpy.ones_like(m_23),
            v_3*v_1*numpy.ones_like(q_13), v_3*v_2*numpy.ones_like(m_23),
            v_3*v_3*numpy.ones_like(q_33)], axis=0)
    if flag_v:
        dder["v"] = numpy.stack([
            (2*q_11*v_1+(q_12+q_12)*v_2+(q_13+q_13)*v_3)*numpy.ones_like(v_1),
            (2*q_22*v_2+(q_12+q_12)*v_1+(q_23+q_23)*v_3)*numpy.ones_like(v_2),
            (2*q_33*v_3+(q_13+q_13)*v_1+(q_23+q_23)*v_2)*numpy.ones_like(v_3)
        ], axis=0)
    return res, dder


def calc_q_v(q, v, flag_q: bool = False, flag_v: bool = False):
    q_11, q_22, q_33 = q[0], q[1], q[2]
    q_12, q_13, q_23 = q[3], q[4], q[5]
    v_1 = v[0]
    v_2 = v[1]
    v_3 = v[2]
    o_1 = q_11*v_1 + q_12*v_2 + q_13*v_3
    o_2 = q_12*v_1 + q_22*v_2 + q_23*v_3
    o_3 = q_13*v_1 + q_23*v_2 + q_33*v_3
    res = numpy.stack([o_1, o_2, o_3], axis=0)
    dder = {}
    if flag_q:
        np_ol_v = np_ol(v_1)
        if q.dtype == complex:
            dder["q_real"] = numpy.stack([
                numpy.stack([v_1*np_ol(q_11.real), np_zl(q_22.real), np_zl(q_33.real), v_2*np_ol(q_12.real), v_3*np_ol(q_13.real), np_zl(q_23.real)], axis=0),
                numpy.stack([np_zl(q_11.real), v_2*np_ol(q_22.real), np_zl(q_33.real), v_1*np_zl(q_12.real), np_zl(q_13.real), v_3*np_zl(q_23.real)], axis=0),
                numpy.stack([np_zl(q_11.real), np_zl(q_22.real), v_3*np_ol(q_33.real), np_zl(q_12.real), v_1*np_zl(q_13.real), v_2*np_zl(q_23.real)], axis=0)
            ], axis=0)
            dder["q_imag"] = numpy.stack([
                numpy.stack([v_1*1j*np_ol(q_11.real), np_zl(q_22.real), np_zl(q_33.real), v_2*1j*np_ol(q_12.real), v_3*1j*np_ol(q_13.real), np_zl(q_23.real)], axis=0),
                numpy.stack([np_zl(q_11.real), v_2*1j*np_ol(q_22.real), np_zl(q_33.real), v_1*1j*np_zl(q_12.real), np_zl(q_13.real), v_3*1j*np_zl(q_23.real)], axis=0),
                numpy.stack([np_zl(q_11.real), np_zl(q_22.real), v_3*1j*np_ol(q_33.real), np_zl(q_12.real), v_1*1j*np_zl(q_13.real), v_2*1j*np_zl(q_23.real)], axis=0)
            ], axis=0)
        else:
            dder["q"] = numpy.stack([
                numpy.stack([v_1*np_ol(q_11.real), np_zl(q_22.real), np_zl(q_33.real), v_2*np_ol(q_12.real), v_3*np_ol(q_13.real), np_zl(q_23.real)], axis=0),
                numpy.stack([np_zl(q_11.real), v_2*np_ol(q_22.real), np_zl(q_33.real), v_1*np_zl(q_12.real), np_zl(q_13.real), v_3*np_zl(q_23.real)], axis=0),
                numpy.stack([np_zl(q_11.real), np_zl(q_22.real), v_3*np_ol(q_33.real), np_zl(q_12.real), v_1*np_zl(q_13.real), v_2*np_zl(q_23.real)], axis=0)
            ], axis=0)
    if flag_v:
        if v.dtype == complex:
            dder["v_real"] = numpy.stack([
                numpy.stack([q_11*np_ol(v_1.real), q_12*np_ol(v_2.real), q_13*np_ol(v_3.real)], axis=0),
                numpy.stack([q_12*np_ol(v_1.real), q_22*np_ol(v_2.real), q_23*np_ol(v_3.real)], axis=0),
                numpy.stack([q_13*np_ol(v_1.real), q_23*np_ol(v_2.real), q_33*np_ol(v_3.real)], axis=0)
            ], axis=0)
            dder["v_imag"] = numpy.stack([
                numpy.stack([q_11*1j*np_ol(v_1.real), q_12*1j*np_ol(v_2.real), q_13*1j*np_ol(v_3.real)], axis=0),
                numpy.stack([q_12*1j*np_ol(v_1.real), q_22*1j*np_ol(v_2.real), q_23*1j*np_ol(v_3.real)], axis=0),
                numpy.stack([q_13*1j*np_ol(v_1.real), q_23*1j*np_ol(v_2.real), q_33*1j*np_ol(v_3.real)], axis=0)
            ], axis=0)
        else:
            dder["v"] = numpy.stack([
                numpy.stack([q_11*np_ol(v_1), q_12*np_ol(v_2), q_13*np_ol(v_3)], axis=0),
                numpy.stack([q_12*np_ol(v_1), q_22*np_ol(v_2), q_23*np_ol(v_3)], axis=0),
                numpy.stack([q_13*np_ol(v_1), q_23*np_ol(v_2), q_33*np_ol(v_3)], axis=0)
            ], axis=0)
    return res, dder


def calc_m_v(m, v, flag_m: bool = False, flag_v: bool = False):
    m_11, m_12, m_13 = m[0], m[1], m[2]
    m_21, m_22, m_23 = m[3], m[4], m[5]
    m_31, m_32, m_33 = m[6], m[7], m[8]
    v_1 = v[0]
    v_2 = v[1]
    v_3 = v[2]
    o_1 = m_11*v_1 + m_12*v_2 + m_13*v_3
    o_2 = m_21*v_1 + m_22*v_2 + m_23*v_3
    o_3 = m_31*v_1 + m_32*v_2 + m_33*v_3
    res = numpy.stack([o_1, o_2, o_3], axis=0)
    dder = {}
    if flag_m:
        np_ol_v = np_ol(v_1)
        if m.dtype == complex:
            dder["m_real"] = numpy.stack([
                numpy.stack([v_1*np_ol(m_11.real), v_2*np_ol(m_12.real), v_3*np_ol(m_13.real), np_zl(m_21.real), np_zl(m_22.real), np_zl(m_23.real), np_zl(m_31.real), np_zl(m_32.real), np_zl(m_33.real)], axis=0),
                numpy.stack([np_zl(m_11.real), np_zl(m_12.real), np_zl(m_13.real), v_1*np_ol(m_21.real), v_2*np_ol(m_22.real), v_3*np_ol(m_23.real), np_zl(m_31.real), np_zl(m_32.real), np_zl(m_33.real)], axis=0),
                numpy.stack([np_zl(m_11.real), np_zl(m_12.real), np_zl(m_13.real), np_zl(m_21.real), np_zl(m_22.real), np_zl(m_23.real), v_1*np_ol(m_31.real), v_2*np_ol(m_32.real), v_3*np_ol(m_33.real)], axis=0)
            ], axis=0)
            dder["m_imag"] = numpy.stack([
                numpy.stack([v_1*1j*np_ol(m_11.real), v_2*1j*np_ol(m_12.real), v_3*1j*np_ol(m_13.real), np_zl(m_21.real), np_zl(m_22.real), np_zl(m_23.real), np_zl(m_31.real), np_zl(m_32.real), np_zl(m_33.real)], axis=0),
                numpy.stack([np_zl(m_11.real), np_zl(m_12.real), np_zl(m_13.real), v_1*1j*np_ol(m_21.real), v_2*1j*np_ol(m_22.real), v_3*1j*np_ol(m_23.real), np_zl(m_31.real), np_zl(m_32.real), np_zl(m_33.real)], axis=0),
                numpy.stack([np_zl(m_11.real), np_zl(m_12.real), np_zl(m_13.real), np_zl(m_21.real), np_zl(m_22.real), np_zl(m_23.real), v_1*1j*np_ol(m_31.real), v_2*1j*np_ol(m_32.real), v_3*1j*np_ol(m_33.real)], axis=0)
            ], axis=0)
        else:
            dder["m"] = numpy.stack([
                numpy.stack([v_1*np_ol(m_11), v_2*np_ol(m_12), v_3*np_ol(m_13), np_ol_v*np_zl(m_21), np_ol_v*np_zl(m_22), np_ol_v*np_zl(m_23), np_ol_v*np_zl(m_31), np_ol_v*np_zl(m_32), np_ol_v*np_zl(m_33)], axis=0),
                numpy.stack([np_ol_v*np_zl(m_11), np_ol_v*np_zl(m_12), np_ol_v*np_zl(m_13), v_1*np_ol(m_21), v_2*np_ol(m_22), v_3*np_ol(m_23), np_ol_v*np_zl(m_31), np_ol_v*np_zl(m_32), np_ol_v*np_zl(m_33)], axis=0),
                numpy.stack([np_ol_v*np_zl(m_11), np_ol_v*np_zl(m_12), np_ol_v*np_zl(m_13), np_ol_v*np_zl(m_21), np_ol_v*np_zl(m_22), np_ol_v*np_zl(m_23), v_1*np_ol(m_31), v_2*np_ol(m_32), v_3*np_ol(m_33)], axis=0)
            ], axis=0)
    if flag_v:
        if v.dtype == complex:
            dder["v_real"] = numpy.stack([
                numpy.stack([m_11*np_ol(v_1.real), m_12*np_ol(v_2.real), m_13*np_ol(v_3.real)], axis=0),
                numpy.stack([m_21*np_ol(v_1.real), m_22*np_ol(v_2.real), m_23*np_ol(v_3.real)], axis=0),
                numpy.stack([m_31*np_ol(v_1.real), m_32*np_ol(v_2.real), m_33*np_ol(v_3.real)], axis=0)
            ], axis=0)
            dder["v_imag"] = numpy.stack([
                numpy.stack([m_11*1j*np_ol(v_1.real), m_12*1j*np_ol(v_2.real), m_13*1j*np_ol(v_3.real)], axis=0),
                numpy.stack([m_21*1j*np_ol(v_1.real), m_22*1j*np_ol(v_2.real), m_23*1j*np_ol(v_3.real)], axis=0),
                numpy.stack([m_31*1j*np_ol(v_1.real), m_32*1j*np_ol(v_2.real), m_33*1j*np_ol(v_3.real)], axis=0)
            ], axis=0)
        else:
            dder["v"] = numpy.stack([
                numpy.stack([m_11*np_ol(v_1), m_12*np_ol(v_2), m_13*np_ol(v_3)], axis=0),
                numpy.stack([m_21*np_ol(v_1), m_22*np_ol(v_2), m_23*np_ol(v_3)], axis=0),
                numpy.stack([m_31*np_ol(v_1), m_32*np_ol(v_2), m_33*np_ol(v_3)], axis=0)
            ], axis=0)
    return res, dder


def calc_det_m(m, flag_m: bool = False):
    """Calculate determinant of the matrix.
    """
    m_11, m_12, m_13 = m[0], m[1], m[2]
    m_21, m_22, m_23 = m[3], m[4], m[5]
    m_31, m_32, m_33 = m[6], m[7], m[8]

    det = m_11*m_22*m_33 - m_11*m_23*m_32 - m_12*m_21*m_33 + m_12*m_23*m_31 + m_13*m_21*m_32 - m_13*m_22*m_31

    dder = {}
    if flag_m:
        hh = numpy.stack([
            (m_22*m_33 - m_23*m_32), (- m_21*m_33 + m_23*m_31), (m_21*m_32 - m_22*m_31),
            (- m_12*m_33 + m_13*m_32), (m_11*m_33 - m_13*m_31), (- m_11*m_32 + m_12*m_31),
            (m_12*m_23 - m_13*m_22), (m_11*m_23 + m_13*m_21), (m_11*m_22 - m_12*m_21)], axis=0)
        if m.type == complex:
            dder["m_real"] = hh
            dder["m_imag"] = 1j*hh
        else:
            dder["m"] = hh
    return det, dder


def calc_inv_m(m_ij, flag_m: bool = False):
    """
    Calculate inversed matrix
    """
    eps = 1.0e-12
    det, dder_det = calc_det_m(m_ij, flag_m=flag_m)

    if numpy.all(det==1):
        inv_det = 1
    else:
        inv_det = numpy.where(numpy.abs(det) < eps, 0., 1./det)

    m_11, m_12, m_13 = m_ij[0], m_ij[1], m_ij[2]
    m_21, m_22, m_23 = m_ij[3], m_ij[4], m_ij[5]
    m_31, m_32, m_33 = m_ij[6], m_ij[7], m_ij[8]

    inv_m_11 = +(m_22*m_33-m_23*m_32) * inv_det
    inv_m_12 = -(m_12*m_33-m_13*m_32) * inv_det
    inv_m_13 = +(m_12*m_23-m_13*m_22) * inv_det
    inv_m_21 = -(m_21*m_33-m_23*m_31) * inv_det
    inv_m_22 = +(m_11*m_33-m_13*m_31) * inv_det
    inv_m_23 = -(m_11*m_23-m_13*m_21) * inv_det
    inv_m_31 = +(m_21*m_32-m_22*m_31) * inv_det
    inv_m_32 = -(m_11*m_32-m_12*m_31) * inv_det
    inv_m_33 = +(m_11*m_22-m_12*m_21) * inv_det

    inv_m_ij = numpy.stack([inv_m_11, inv_m_12, inv_m_13, inv_m_21, inv_m_22, inv_m_23, inv_m_31, inv_m_32, inv_m_33], axis=0)
    dder = {}
    if flag_m:
        inv_det_sq = numpy.square(inv_det)
        d11_11 = (                - (m_22*m_33-m_23*m_32) * inv_det_sq * dder_det["m_ij"][0]) * numpy.ones_like(m_11)
        d11_12 = (                - (m_22*m_33-m_23*m_32) * inv_det_sq * dder_det["m_ij"][1]) * numpy.ones_like(m_12)
        d11_13 = (                - (m_22*m_33-m_23*m_32) * inv_det_sq * dder_det["m_ij"][2]) * numpy.ones_like(m_13)
        d11_21 = (                - (m_22*m_33-m_23*m_32) * inv_det_sq * dder_det["m_ij"][3]) * numpy.ones_like(m_21)
        d11_22 = ( m_33 * inv_det - (m_22*m_33-m_23*m_32) * inv_det_sq * dder_det["m_ij"][4]) * numpy.ones_like(m_22)
        d11_23 = (-m_32 * inv_det - (m_22*m_33-m_23*m_32) * inv_det_sq * dder_det["m_ij"][5]) * numpy.ones_like(m_23)
        d11_31 = (                - (m_22*m_33-m_23*m_32) * inv_det_sq * dder_det["m_ij"][6]) * numpy.ones_like(m_31)
        d11_32 = (-m_23 * inv_det - (m_22*m_33-m_23*m_32) * inv_det_sq * dder_det["m_ij"][7]) * numpy.ones_like(m_32)
        d11_33 = ( m_22 * inv_det - (m_22*m_33-m_23*m_32) * inv_det_sq * dder_det["m_ij"][8]) * numpy.ones_like(m_33)
        d12_11 = (                + (m_12*m_33-m_13*m_32) * inv_det_sq * dder_det["m_ij"][0]) *numpy.ones_like(m_11)
        d12_12 = (-m_33 * inv_det + (m_12*m_33-m_13*m_32) * inv_det_sq * dder_det["m_ij"][1]) *numpy.ones_like(m_12)
        d12_13 = ( m_32 * inv_det + (m_12*m_33-m_13*m_32) * inv_det_sq * dder_det["m_ij"][2]) *numpy.ones_like(m_13)
        d12_21 = (                + (m_12*m_33-m_13*m_32) * inv_det_sq * dder_det["m_ij"][3]) *numpy.ones_like(m_21)
        d12_22 = (                + (m_12*m_33-m_13*m_32) * inv_det_sq * dder_det["m_ij"][4]) *numpy.ones_like(m_22)
        d12_23 = (                + (m_12*m_33-m_13*m_32) * inv_det_sq * dder_det["m_ij"][5]) *numpy.ones_like(m_23)
        d12_31 = (                + (m_12*m_33-m_13*m_32) * inv_det_sq * dder_det["m_ij"][6]) *numpy.ones_like(m_31)
        d12_32 = ( m_13 * inv_det + (m_12*m_33-m_13*m_32) * inv_det_sq * dder_det["m_ij"][7]) *numpy.ones_like(m_32)
        d12_33 = (-m_12 * inv_det + (m_12*m_33-m_13*m_32) * inv_det_sq * dder_det["m_ij"][8]) *numpy.ones_like(m_33)
        d13_11 = (                - (m_12*m_23-m_13*m_22) * inv_det_sq * dder_det["m_ij"][0]) *numpy.ones_like(m_11)
        d13_12 = ( m_23 * inv_det - (m_12*m_23-m_13*m_22) * inv_det_sq * dder_det["m_ij"][1]) *numpy.ones_like(m_12)
        d13_13 = (-m_22 * inv_det - (m_12*m_23-m_13*m_22) * inv_det_sq * dder_det["m_ij"][2]) *numpy.ones_like(m_13)
        d13_21 = (                - (m_12*m_23-m_13*m_22) * inv_det_sq * dder_det["m_ij"][3]) *numpy.ones_like(m_21)
        d13_22 = (-m_13 * inv_det - (m_12*m_23-m_13*m_22) * inv_det_sq * dder_det["m_ij"][4]) *numpy.ones_like(m_22)
        d13_23 = ( m_12 * inv_det - (m_12*m_23-m_13*m_22) * inv_det_sq * dder_det["m_ij"][5]) *numpy.ones_like(m_23)
        d13_31 = (                - (m_12*m_23-m_13*m_22) * inv_det_sq * dder_det["m_ij"][6]) *numpy.ones_like(m_31)
        d13_32 = (                - (m_12*m_23-m_13*m_22) * inv_det_sq * dder_det["m_ij"][7]) *numpy.ones_like(m_32)
        d13_33 = (                - (m_12*m_23-m_13*m_22) * inv_det_sq * dder_det["m_ij"][8]) *numpy.ones_like(m_33)
        d21_11 = (                + (m_21*m_33-m_23*m_31) * inv_det_sq * dder_det["m_ij"][0]) *numpy.ones_like(m_11)
        d21_12 = (                + (m_21*m_33-m_23*m_31) * inv_det_sq * dder_det["m_ij"][1]) *numpy.ones_like(m_12)
        d21_13 = (                + (m_21*m_33-m_23*m_31) * inv_det_sq * dder_det["m_ij"][2]) *numpy.ones_like(m_13)
        d21_21 = (-m_33 * inv_det + (m_21*m_33-m_23*m_31) * inv_det_sq * dder_det["m_ij"][3]) *numpy.ones_like(m_21)
        d21_22 = (                + (m_21*m_33-m_23*m_31) * inv_det_sq * dder_det["m_ij"][4]) *numpy.ones_like(m_22)
        d21_23 = ( m_31 * inv_det + (m_21*m_33-m_23*m_31) * inv_det_sq * dder_det["m_ij"][5]) *numpy.ones_like(m_23)
        d21_31 = ( m_23 * inv_det + (m_21*m_33-m_23*m_31) * inv_det_sq * dder_det["m_ij"][6]) *numpy.ones_like(m_31)
        d21_32 = (                + (m_21*m_33-m_23*m_31) * inv_det_sq * dder_det["m_ij"][7]) *numpy.ones_like(m_32)
        d21_33 = (-m_21 * inv_det + (m_21*m_33-m_23*m_31) * inv_det_sq * dder_det["m_ij"][8]) *numpy.ones_like(m_33)
        d22_11 = ( m_33 * inv_det - (m_11*m_33-m_13*m_31) * inv_det_sq * dder_det["m_ij"][0]) *numpy.ones_like(m_11)
        d22_12 = (                - (m_11*m_33-m_13*m_31) * inv_det_sq * dder_det["m_ij"][1]) *numpy.ones_like(m_12)
        d22_13 = (-m_31 * inv_det - (m_11*m_33-m_13*m_31) * inv_det_sq * dder_det["m_ij"][2]) *numpy.ones_like(m_13)
        d22_21 = (                - (m_11*m_33-m_13*m_31) * inv_det_sq * dder_det["m_ij"][3]) *numpy.ones_like(m_21)
        d22_22 = (                - (m_11*m_33-m_13*m_31) * inv_det_sq * dder_det["m_ij"][4]) *numpy.ones_like(m_22)
        d22_23 = (                - (m_11*m_33-m_13*m_31) * inv_det_sq * dder_det["m_ij"][5]) *numpy.ones_like(m_23)
        d22_31 = (-m_13 * inv_det - (m_11*m_33-m_13*m_31) * inv_det_sq * dder_det["m_ij"][6]) *numpy.ones_like(m_31)
        d22_32 = (                - (m_11*m_33-m_13*m_31) * inv_det_sq * dder_det["m_ij"][7]) *numpy.ones_like(m_32)
        d22_33 = ( m_11 * inv_det - (m_11*m_33-m_13*m_31) * inv_det_sq * dder_det["m_ij"][8]) *numpy.ones_like(m_33)
        d23_11 = (-m_23 * inv_det + (m_11*m_23-m_13*m_21) * inv_det_sq * dder_det["m_ij"][0]) *numpy.ones_like(m_11)
        d23_12 = (                + (m_11*m_23-m_13*m_21) * inv_det_sq * dder_det["m_ij"][1]) *numpy.ones_like(m_12)
        d23_13 = ( m_21 * inv_det + (m_11*m_23-m_13*m_21) * inv_det_sq * dder_det["m_ij"][2]) *numpy.ones_like(m_13)
        d23_21 = ( m_13 * inv_det + (m_11*m_23-m_13*m_21) * inv_det_sq * dder_det["m_ij"][3]) *numpy.ones_like(m_21)
        d23_22 = (                + (m_11*m_23-m_13*m_21) * inv_det_sq * dder_det["m_ij"][4]) *numpy.ones_like(m_22)
        d23_23 = (-m_11 * inv_det + (m_11*m_23-m_13*m_21) * inv_det_sq * dder_det["m_ij"][5]) *numpy.ones_like(m_23)
        d23_31 = (                + (m_11*m_23-m_13*m_21) * inv_det_sq * dder_det["m_ij"][6]) *numpy.ones_like(m_31)
        d23_32 = (                + (m_11*m_23-m_13*m_21) * inv_det_sq * dder_det["m_ij"][7]) *numpy.ones_like(m_32)
        d23_33 = (                + (m_11*m_23-m_13*m_21) * inv_det_sq * dder_det["m_ij"][8]) *numpy.ones_like(m_33)
        d31_11 = (                - (m_21*m_32-m_22*m_31) * inv_det_sq * dder_det["m_ij"][0]) *numpy.ones_like(m_11)
        d31_12 = (                - (m_21*m_32-m_22*m_31) * inv_det_sq * dder_det["m_ij"][1]) *numpy.ones_like(m_12)
        d31_13 = (                - (m_21*m_32-m_22*m_31) * inv_det_sq * dder_det["m_ij"][2]) *numpy.ones_like(m_13)
        d31_21 = ( m_32 * inv_det - (m_21*m_32-m_22*m_31) * inv_det_sq * dder_det["m_ij"][3]) *numpy.ones_like(m_21)
        d31_22 = (-m_31 * inv_det - (m_21*m_32-m_22*m_31) * inv_det_sq * dder_det["m_ij"][4]) *numpy.ones_like(m_22)
        d31_23 = (                - (m_21*m_32-m_22*m_31) * inv_det_sq * dder_det["m_ij"][5]) *numpy.ones_like(m_23)
        d31_31 = (-m_22 * inv_det - (m_21*m_32-m_22*m_31) * inv_det_sq * dder_det["m_ij"][6]) *numpy.ones_like(m_31)
        d31_32 = ( m_21 * inv_det - (m_21*m_32-m_22*m_31) * inv_det_sq * dder_det["m_ij"][7]) *numpy.ones_like(m_32)
        d31_33 = (                - (m_21*m_32-m_22*m_31) * inv_det_sq * dder_det["m_ij"][8]) *numpy.ones_like(m_33)
        d32_11 = (-m_32 * inv_det + (m_11*m_32-m_12*m_31) * inv_det_sq * dder_det["m_ij"][0]) *numpy.ones_like(m_11)
        d32_12 = ( m_31 * inv_det + (m_11*m_32-m_12*m_31) * inv_det_sq * dder_det["m_ij"][1]) *numpy.ones_like(m_12)
        d32_13 = (                + (m_11*m_32-m_12*m_31) * inv_det_sq * dder_det["m_ij"][2]) *numpy.ones_like(m_13)
        d32_21 = (                + (m_11*m_32-m_12*m_31) * inv_det_sq * dder_det["m_ij"][3]) *numpy.ones_like(m_21)
        d32_22 = (                + (m_11*m_32-m_12*m_31) * inv_det_sq * dder_det["m_ij"][4]) *numpy.ones_like(m_22)
        d32_23 = (                + (m_11*m_32-m_12*m_31) * inv_det_sq * dder_det["m_ij"][5]) *numpy.ones_like(m_23)
        d32_31 = ( m_12 * inv_det + (m_11*m_32-m_12*m_31) * inv_det_sq * dder_det["m_ij"][6]) *numpy.ones_like(m_31)
        d32_32 = (-m_11 * inv_det + (m_11*m_32-m_12*m_31) * inv_det_sq * dder_det["m_ij"][7]) *numpy.ones_like(m_32)
        d32_33 = (                + (m_11*m_32-m_12*m_31) * inv_det_sq * dder_det["m_ij"][8]) *numpy.ones_like(m_33)
        d33_11 = ( m_22 * inv_det - (m_11*m_22-m_12*m_21) * inv_det_sq * dder_det["m_ij"][0]) *numpy.ones_like(m_11)
        d33_12 = (-m_21 * inv_det - (m_11*m_22-m_12*m_21) * inv_det_sq * dder_det["m_ij"][1]) *numpy.ones_like(m_12)
        d33_13 = (                - (m_11*m_22-m_12*m_21) * inv_det_sq * dder_det["m_ij"][2]) *numpy.ones_like(m_13)
        d33_21 = (-m_12 * inv_det - (m_11*m_22-m_12*m_21) * inv_det_sq * dder_det["m_ij"][3]) *numpy.ones_like(m_21)
        d33_22 = (+m_11 * inv_det - (m_11*m_22-m_12*m_21) * inv_det_sq * dder_det["m_ij"][4]) *numpy.ones_like(m_22)
        d33_23 = (                - (m_11*m_22-m_12*m_21) * inv_det_sq * dder_det["m_ij"][5]) *numpy.ones_like(m_23)
        d33_31 = (                - (m_11*m_22-m_12*m_21) * inv_det_sq * dder_det["m_ij"][6]) *numpy.ones_like(m_31)
        d33_32 = (                - (m_11*m_22-m_12*m_21) * inv_det_sq * dder_det["m_ij"][7]) *numpy.ones_like(m_32)
        d33_33 = (                - (m_11*m_22-m_12*m_21) * inv_det_sq * dder_det["m_ij"][8]) *numpy.ones_like(m_33)
        dder["m"] = numpy.stack([
            numpy.stack([d11_11, d11_12, d11_13, d11_21, d11_22, d11_23, d11_31, d11_32, d11_33], axis=0),
            numpy.stack([d12_11, d12_12, d12_13, d12_21, d12_22, d12_23, d12_31, d12_32, d12_33], axis=0),
            numpy.stack([d13_11, d13_12, d13_13, d13_21, d13_22, d13_23, d13_31, d13_32, d13_33], axis=0),
            numpy.stack([d21_11, d21_12, d21_13, d21_21, d21_22, d21_23, d21_31, d21_32, d21_33], axis=0),
            numpy.stack([d22_11, d22_12, d22_13, d22_21, d22_22, d22_23, d22_31, d22_32, d22_33], axis=0),
            numpy.stack([d23_11, d23_12, d23_13, d23_21, d23_22, d23_23, d23_31, d23_32, d23_33], axis=0),
            numpy.stack([d31_11, d31_12, d31_13, d31_21, d31_22, d31_23, d31_31, d31_32, d31_33], axis=0),
            numpy.stack([d32_11, d32_12, d32_13, d32_21, d32_22, d32_23, d32_31, d32_32, d32_33], axis=0),
            numpy.stack([d33_11, d33_12, d33_13, d33_21, d33_22, d33_23, d33_31, d33_32, d33_33], axis=0),
        ], axis=0)
    return inv_m_ij, dder


def calc_mm_as_m_q_inv_m(m, flag_m: bool = False):
    r_11, r_12, r_13 = m[0], m[1], m[2]
    r_21, r_22, r_23 = m[3], m[4], m[5]
    r_31, r_32, r_33 = m[6], m[7], m[8]

    det_r = numpy.expand_dims(numpy.expand_dims(
        r_11*r_22*r_33 - r_11*r_23*r_32  - r_12*r_21*r_33 + r_12*r_23*r_31  + r_13*r_21*r_32 - r_13*r_22*r_31,
        axis=0), axis=0)

    mm = numpy.array([
        [r_11*(r_22*r_33 - r_23*r_32), r_12*(-r_21*r_33 + r_23*r_31), r_13*(r_21*r_32 - r_22*r_31), r_11*(-r_21*r_33 + r_23*r_31) + r_12*(r_22*r_33 - r_23*r_32), r_11*(r_21*r_32 - r_22*r_31) + r_13*(r_22*r_33 - r_23*r_32), r_12*(r_21*r_32 - r_22*r_31) + r_13*(-r_21*r_33 + r_23*r_31)],
        [r_11*(-r_12*r_33 + r_13*r_32), r_12*(r_11*r_33 - r_13*r_31), r_13*(-r_11*r_32 + r_12*r_31), r_11*(r_11*r_33 - r_13*r_31) + r_12*(-r_12*r_33 + r_13*r_32), r_11*(-r_11*r_32 + r_12*r_31) + r_13*(-r_12*r_33 + r_13*r_32), r_12*(-r_11*r_32 + r_12*r_31) + r_13*(r_11*r_33 - r_13*r_31)],
        [r_11*(r_12*r_23 - r_13*r_22), r_12*(-r_11*r_23 + r_13*r_21), r_13*(r_11*r_22 - r_12*r_21), r_11*(-r_11*r_23 + r_13*r_21) + r_12*(r_12*r_23 - r_13*r_22), r_11*(r_11*r_22 - r_12*r_21) + r_13*(r_12*r_23 - r_13*r_22), r_12*(r_11*r_22 - r_12*r_21) + r_13*(-r_11*r_23 + r_13*r_21)],
        [r_21*(r_22*r_33 - r_23*r_32), r_22*(-r_21*r_33 + r_23*r_31), r_23*(r_21*r_32 - r_22*r_31), r_21*(-r_21*r_33 + r_23*r_31) + r_22*(r_22*r_33 - r_23*r_32), r_21*(r_21*r_32 - r_22*r_31) + r_23*(r_22*r_33 - r_23*r_32), r_22*(r_21*r_32 - r_22*r_31) + r_23*(-r_21*r_33 + r_23*r_31)],
        [r_21*(-r_12*r_33 + r_13*r_32), r_22*(r_11*r_33 - r_13*r_31), r_23*(-r_11*r_32 + r_12*r_31), r_21*(r_11*r_33 - r_13*r_31) + r_22*(-r_12*r_33 + r_13*r_32), r_21*(-r_11*r_32 + r_12*r_31) + r_23*(-r_12*r_33 + r_13*r_32), r_22*(-r_11*r_32 + r_12*r_31) + r_23*(r_11*r_33 - r_13*r_31)],
        [r_21*(r_12*r_23 - r_13*r_22), r_22*(-r_11*r_23 + r_13*r_21), r_23*(r_11*r_22 - r_12*r_21), r_21*(-r_11*r_23 + r_13*r_21) + r_22*(r_12*r_23 - r_13*r_22), r_21*(r_11*r_22 - r_12*r_21) + r_23*(r_12*r_23 - r_13*r_22), r_22*(r_11*r_22 - r_12*r_21) + r_23*(-r_11*r_23 + r_13*r_21)],
        [r_31*(r_22*r_33 - r_23*r_32), r_32*(-r_21*r_33 + r_23*r_31), r_33*(r_21*r_32 - r_22*r_31), r_31*(-r_21*r_33 + r_23*r_31) + r_32*(r_22*r_33 - r_23*r_32), r_31*(r_21*r_32 - r_22*r_31) + r_33*(r_22*r_33 - r_23*r_32), r_32*(r_21*r_32 - r_22*r_31) + r_33*(-r_21*r_33 + r_23*r_31)],
        [r_31*(-r_12*r_33 + r_13*r_32), r_32*(r_11*r_33 - r_13*r_31), r_33*(-r_11*r_32 + r_12*r_31), r_31*(r_11*r_33 - r_13*r_31) + r_32*(-r_12*r_33 + r_13*r_32), r_31*(-r_11*r_32 + r_12*r_31) + r_33*(-r_12*r_33 + r_13*r_32), r_32*(-r_11*r_32 + r_12*r_31) + r_33*(r_11*r_33 - r_13*r_31)],
        [r_31*(r_12*r_23 - r_13*r_22), r_32*(-r_11*r_23 + r_13*r_21), r_33*(r_11*r_22 - r_12*r_21), r_31*(-r_11*r_23 + r_13*r_21) + r_32*(r_12*r_23 - r_13*r_22), r_31*(r_11*r_22 - r_12*r_21) + r_33*(r_12*r_23 - r_13*r_22), r_32*(r_11*r_22 - r_12*r_21) + r_33*(-r_11*r_23 + r_13*r_21)]],
        dtype=m.dtype)
    mm /= det_r
    dder_mm = {}
    return mm, dder_mm


def calc_mm_as_m1_m2_inv_m1(m1, flag_m1: bool = False):
    r_11, r_12, r_13 = m1[0], m1[1], m1[2]
    r_21, r_22, r_23 = m1[3], m1[4], m1[5]
    r_31, r_32, r_33 = m1[6], m1[7], m1[8]

    det_r = numpy.expand_dims(numpy.expand_dims(
        r_11*r_22*r_33 - r_11*r_23*r_32  - r_12*r_21*r_33 + r_12*r_23*r_31  + r_13*r_21*r_32 - r_13*r_22*r_31,
        axis=0), axis=0)

    mm = numpy.array([
        [r_11*(r_22*r_33 - r_23*r_32), r_11*(-r_21*r_33 + r_23*r_31), r_11*(r_21*r_32 - r_22*r_31), r_12*(r_22*r_33 - r_23*r_32), r_12*(-r_21*r_33 + r_23*r_31), r_12*(r_21*r_32 - r_22*r_31), r_13*(r_22*r_33 - r_23*r_32), r_13*(-r_21*r_33 + r_23*r_31), r_13*(r_21*r_32 - r_22*r_31)],
        [r_11*(-r_12*r_33 + r_13*r_32), r_11*(r_11*r_33 - r_13*r_31), r_11*(-r_11*r_32 + r_12*r_31), r_12*(-r_12*r_33 + r_13*r_32), r_12*(r_11*r_33 - r_13*r_31), r_12*(-r_11*r_32 + r_12*r_31), r_13*(-r_12*r_33 + r_13*r_32), r_13*(r_11*r_33 - r_13*r_31), r_13*(-r_11*r_32 + r_12*r_31)],
        [r_11*(r_12*r_23 - r_13*r_22), r_11*(-r_11*r_23 + r_13*r_21), r_11*(r_11*r_22 - r_12*r_21), r_12*(r_12*r_23 - r_13*r_22), r_12*(-r_11*r_23 + r_13*r_21), r_12*(r_11*r_22 - r_12*r_21), r_13*(r_12*r_23 - r_13*r_22), r_13*(-r_11*r_23 + r_13*r_21), r_13*(r_11*r_22 - r_12*r_21)],
        [r_21*(r_22*r_33 - r_23*r_32), r_21*(-r_21*r_33 + r_23*r_31), r_21*(r_21*r_32 - r_22*r_31), r_22*(r_22*r_33 - r_23*r_32), r_22*(-r_21*r_33 + r_23*r_31), r_22*(r_21*r_32 - r_22*r_31), r_23*(r_22*r_33 - r_23*r_32), r_23*(-r_21*r_33 + r_23*r_31), r_23*(r_21*r_32 - r_22*r_31)],
        [r_21*(-r_12*r_33 + r_13*r_32), r_21*(r_11*r_33 - r_13*r_31), r_21*(-r_11*r_32 + r_12*r_31), r_22*(-r_12*r_33 + r_13*r_32), r_22*(r_11*r_33 - r_13*r_31), r_22*(-r_11*r_32 + r_12*r_31), r_23*(-r_12*r_33 + r_13*r_32), r_23*(r_11*r_33 - r_13*r_31), r_23*(-r_11*r_32 + r_12*r_31)],
        [r_21*(r_12*r_23 - r_13*r_22), r_21*(-r_11*r_23 + r_13*r_21), r_21*(r_11*r_22 - r_12*r_21), r_22*(r_12*r_23 - r_13*r_22), r_22*(-r_11*r_23 + r_13*r_21), r_22*(r_11*r_22 - r_12*r_21), r_23*(r_12*r_23 - r_13*r_22), r_23*(-r_11*r_23 + r_13*r_21), r_23*(r_11*r_22 - r_12*r_21)],
        [r_31*(r_22*r_33 - r_23*r_32), r_31*(-r_21*r_33 + r_23*r_31), r_31*(r_21*r_32 - r_22*r_31), r_32*(r_22*r_33 - r_23*r_32), r_32*(-r_21*r_33 + r_23*r_31), r_32*(r_21*r_32 - r_22*r_31), r_33*(r_22*r_33 - r_23*r_32), r_33*(-r_21*r_33 + r_23*r_31), r_33*(r_21*r_32 - r_22*r_31)],
        [r_31*(-r_12*r_33 + r_13*r_32), r_31*(r_11*r_33 - r_13*r_31), r_31*(-r_11*r_32 + r_12*r_31), r_32*(-r_12*r_33 + r_13*r_32), r_32*(r_11*r_33 - r_13*r_31), r_32*(-r_11*r_32 + r_12*r_31), r_33*(-r_12*r_33 + r_13*r_32), r_33*(r_11*r_33 - r_13*r_31), r_33*(-r_11*r_32 + r_12*r_31)],
        [r_31*(r_12*r_23 - r_13*r_22), r_31*(-r_11*r_23 + r_13*r_21), r_31*(r_11*r_22 - r_12*r_21), r_32*(r_12*r_23 - r_13*r_22), r_32*(-r_11*r_23 + r_13*r_21), r_32*(r_11*r_22 - r_12*r_21), r_33*(r_12*r_23 - r_13*r_22), r_33*(-r_11*r_23 + r_13*r_21), r_33*(r_11*r_22 - r_12*r_21)]],
        dtype=m1.dtype)
    mm /= det_r
    dder_mm = {}
    return mm, dder_mm