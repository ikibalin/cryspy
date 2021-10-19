# -*- coding: utf-8 -*-
"""
Expression for unit cell parameters are presented.
For details see documentation.
"""
import warnings
import numpy 

from .matrix_operations import calc_m_q_mt


def calc_reciprocal_by_unit_cell_parameters(
        unit_cell_parameters, flag_unit_cell_parameters: bool = False):
    a, b = unit_cell_parameters[0], unit_cell_parameters[1]
    c = unit_cell_parameters[2]

    alpha, beta = unit_cell_parameters[3], unit_cell_parameters[4]
    gamma = unit_cell_parameters[5]
    c_a, c_b, c_g = numpy.cos(alpha), numpy.cos(beta), numpy.cos(gamma)
    s_a, s_b, s_g = numpy.sin(alpha), numpy.sin(beta), numpy.sin(gamma)
    vol, dder_vol = \
        calc_volume_uc_by_unit_cell_parameters(
            unit_cell_parameters, flag_unit_cell_parameters=flag_unit_cell_parameters)

    ialpha = numpy.arccos((c_b*c_g-c_a)/(s_b*s_g))
    ibeta = numpy.arccos((c_g*c_a-c_b)/(s_g*s_a))
    igamma = numpy.arccos((c_a*c_b-c_g)/(s_a*s_b))
    ia, ib, ic = b*c*s_a/vol, c*a*s_b/vol, a*b*s_g/vol
    reciprocal_unit_cell_parameters = numpy.stack([
        ia, ib, ic, ialpha, ibeta, igamma], axis=0)
    
    dder = {}
    if flag_unit_cell_parameters:
        d_v = dder_vol["unit_cell_parameters"]
        zeros = numpy.zeros(unit_cell_parameters[0].shape, dtype=float)
        i_vol_sq = 1./numpy.square(vol)
        der_ia = numpy.stack([b*c*s_a*d_v[0]*i_vol_sq,
                              c*s_a/vol+b*c*s_a*d_v[1]*i_vol_sq,
                              b*s_a/vol+b*c*s_a*d_v[2]*i_vol_sq,
                              b*c*c_a/vol+b*c*s_a*d_v[3]*i_vol_sq,
                              b*c*s_a*d_v[4]*i_vol_sq,
                              b*c*s_a*d_v[5]*i_vol_sq], axis=0)
        der_ib = numpy.stack([c*s_b/vol+c*a*s_b*i_vol_sq*d_v[0],
                              c*a*s_b*i_vol_sq*d_v[1],
                              a*s_b/vol+c*a*s_b*i_vol_sq*d_v[2],
                              c*a*s_b*i_vol_sq*d_v[3],
                              c*a*c_b/vol+c*a*s_b*i_vol_sq*d_v[4],
                              c*a*s_b*i_vol_sq*d_v[5]], axis=0)
        der_ic = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_ial = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_ibe = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_iga = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)

        dder["unit_cell_parameters"] = numpy.stack([
            der_ia, der_ib, der_ic, der_ial, der_ibe, der_iga], axis=0)
    return reciprocal_unit_cell_parameters, dder


def calc_phi_sq_by_unit_cell_parameters(
        unit_cell_parameters, flag_unit_cell_parameters: bool = False):
    """
    Calculate phi_sq
    phi_sq = 1 - cos^2 alpha - cos^2 beta - cos^2 gamma - 2 cos alpha cos beta cos gamma
    """
    alpha, beta = unit_cell_parameters[3], unit_cell_parameters[4]
    gamma = unit_cell_parameters[5]
    c_a, c_b, c_g = numpy.cos(alpha), numpy.cos(beta), numpy.cos(gamma)
    phi_sq =  (
        1. - numpy.square(c_a) - numpy.square(c_b) - numpy.square(c_g) +
        2.*c_a*c_b*c_g)
    dder = {}
    if flag_unit_cell_parameters:
        der_a = numpy.zeros_like(unit_cell_parameters[0])
        s_a, s_b, s_g = numpy.sin(alpha), numpy.sin(beta), numpy.sin(gamma)
        der_al = 2*s_a*(c_a-c_b*c_g)
        der_be = 2*s_b*(c_b-c_g*c_a)
        der_ga = 2*s_g*(c_g-c_a*c_b)
        dder["unit_cell_parameters"] = numpy.stack([
            der_a, der_a, der_a, der_al, der_be, der_ga], axis=0)
    return phi_sq, dder


def calc_phi_by_unit_cell_parameters(
        unit_cell_parameters, flag_unit_cell_parameters: bool = False):
    """
    Calculate phi
    phi = (1 - cos^2 alpha - cos^2 beta - cos^2 gamma - 2 cos alpha cos beta cos gamma)^0.5
    """
    phi_sq, dder_phi_sq = calc_phi_sq_by_unit_cell_parameters(
        unit_cell_parameters, flag_unit_cell_parameters)
    phi = numpy.sqrt(phi_sq)
    dder = {}
    if flag_unit_cell_parameters:
        dder["unit_cell_parameters"] = 0.5*dder_phi_sq["unit_cell_parameters"]/phi 
    return phi, dder


def calc_volume_uc_by_unit_cell_parameters(
        unit_cell_parameters, flag_unit_cell_parameters: bool = False):
    """
    Calculate the volume of unit cell defined as 
    a, b, c, alpha, beta, gamma.
    Angles should be given in radians.
    """
    a, b = unit_cell_parameters[0], unit_cell_parameters[1]
    c = unit_cell_parameters[2]
    alpha, beta = unit_cell_parameters[3], unit_cell_parameters[4]
    gamma = unit_cell_parameters[5]
    c_a, c_b, c_g = numpy.cos(alpha), numpy.cos(beta), numpy.cos(gamma)
    c_a_sq, c_b_sq = numpy.square(c_a), numpy.square(c_b)
    c_g_sq = numpy.square(c_g)
    phi, dder_phi = calc_phi_by_unit_cell_parameters(unit_cell_parameters, flag_unit_cell_parameters)
    volume_uc = a*b*c*phi
    dder = {}
    if flag_unit_cell_parameters:
        s_a, s_b, s_g = numpy.sin(alpha), numpy.sin(beta), numpy.sin(gamma)
        der_a = b*c*phi
        der_b = a*c*phi
        der_c = a*b*phi
        der_al = -a*b*c*0.5*(-2*c_a_sq+2.*c_b*c_g)*s_a/phi
        der_be = -a*b*c*0.5*(-2*c_b_sq+2.*c_a*c_g)*s_b/phi
        der_ga = -a*b*c*0.5*(-2*c_g_sq+2.*c_a*c_b)*s_g/phi
        dder["unit_cell_parameters"] = numpy.stack([
            der_a, der_b, der_c, der_al, der_be, der_ga], axis=0)
    return volume_uc, dder


def calc_volume_ruc_by_unit_cell_parameters(
        unit_cell_parameters, flag_unit_cell_parameters: bool = False):
    """
    Calculate the volume of reciprocal unit cell 
    """
    volume_uc , dder_volume_uc = \
        calc_volume_uc_by_unit_cell_parameters(
            unit_cell_parameters,
            flag_unit_cell_parameters=flag_unit_cell_parameters)
    volume_ruc = 1./volume_uc
    dder = {}
    if flag_unit_cell_parameters:
        dder["unit_cell_parameters"] = -numpy.square(volume_ruc)*dder_volume_uc["unit_cell_parameters"]
    return volume_ruc, dder


def calc_m_g_by_unit_cell_parameters(
        unit_cell_parameters, flag_unit_cell_parameters: bool = False):
    a, b = unit_cell_parameters[0], unit_cell_parameters[1]
    c = unit_cell_parameters[2]
    alpha, beta = unit_cell_parameters[3], unit_cell_parameters[4]
    gamma = unit_cell_parameters[5]
    c_a, c_b, c_g = numpy.cos(alpha), numpy.cos(beta), numpy.cos(gamma)

    m_g_ij = numpy.stack([a*a,b*b, c*c, a*b*c_g, a*c*c_b,b*c*c_a])
    dder = {}
    if flag_unit_cell_parameters:
        s_a, s_b, s_g = numpy.sin(alpha), numpy.sin(beta), numpy.sin(gamma)
        zeros = numpy.zeros_like(unit_cell_parameters[0])
        der_11 = numpy.stack([2*a, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_12 = numpy.stack([b*c_g, a*c_g, zeros, zeros, zeros, -a*b*s_g],
                             axis=0)
        der_13 = numpy.stack([c*c_b, zeros, a*c_b, zeros, -a*c*s_b, zeros],
                             axis=0)
        der_22 = numpy.stack([zeros, 2*b, zeros, zeros, zeros, zeros], axis=0)
        der_23 = numpy.stack([zeros, c*c_a, b*c_a, -b*c*s_a, zeros, zeros], axis=0)
        der_33 = numpy.stack([zeros, zeros, 2*c, zeros, zeros, zeros], axis=0)

        dder["unit_cell_parameters"] = numpy.stack([
            der_11, der_22, der_33, der_12, der_13, der_23], axis=0)
    return m_g_ij, dder


def calc_m_g_norm_by_unit_cell_parameters(
        unit_cell_parameters, flag_unit_cell_parameters: bool = False):
    alpha, beta = unit_cell_parameters[3], unit_cell_parameters[4]
    gamma = unit_cell_parameters[5]
    c_a, c_b, c_g = numpy.cos(alpha), numpy.cos(beta), numpy.cos(gamma)

    ones = numpy.ones_like(c_a)

    m_g_norm_ij = numpy.stack([ones,ones,ones, c_g, c_b, c_a])
    dder = {}
    if flag_unit_cell_parameters:
        zeros = numpy.zeros_like(c_a)

        s_a, s_b, s_g = numpy.sin(alpha), numpy.sin(beta), numpy.sin(gamma)
        der_11 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_12 = numpy.stack([c_g, c_g, zeros, zeros, zeros, -s_g], axis=0)
        der_13 = numpy.stack([c_b, zeros, c_b, zeros, -s_b, zeros], axis=0)
        der_22 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_23 = numpy.stack([zeros, c_a, c_a, -s_a, zeros, zeros], axis=0)
        der_33 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)

        dder["unit_cell_parameters"] = numpy.stack([
            der_11, der_22, der_33, der_12, der_13, der_23], axis=0)
    return m_g_norm_ij, dder


def calc_m_reciprocal_g_by_unit_cell_parameters(
        unit_cell_parameters, flag_unit_cell_parameters: bool = False):
    a, b = unit_cell_parameters[0], unit_cell_parameters[1]
    c = unit_cell_parameters[2]
    alpha, beta = unit_cell_parameters[3], unit_cell_parameters[4]
    gamma = unit_cell_parameters[5]
    c_a, c_b, c_g = numpy.cos(alpha), numpy.cos(beta), numpy.cos(gamma)
    s_a, s_b, s_g = numpy.sin(alpha), numpy.sin(beta), numpy.sin(gamma)
    phi_sq, dder_phi_sq = calc_phi_sq_by_unit_cell_parameters(
        unit_cell_parameters, flag_unit_cell_parameters=flag_unit_cell_parameters) 

    inv_g_11 = numpy.square(s_a) / (numpy.square(a)* phi_sq)
    inv_g_22 = numpy.square(s_b) / (numpy.square(b)* phi_sq)
    inv_g_33 = numpy.square(s_g) / (numpy.square(c)* phi_sq)
    inv_g_12 = (c_a*c_b-c_g) / (a*b* phi_sq)
    inv_g_13 = (c_a*c_g-c_b) / (a*c* phi_sq)
    inv_g_23 = (c_b*c_g-c_a) / (b*c* phi_sq)

    # G matrix for reciprocal unit cell
    m_reciprocal_g_ij = numpy.stack(
        [inv_g_11, inv_g_22, inv_g_33,
         inv_g_12, inv_g_13, inv_g_23],
        axis=0)
    dder = {}
    if flag_unit_cell_parameters:
        zeros = numpy.zeros_like(unit_cell_parameters[0])
        der_11 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_12 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_13 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_22 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_23 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_33 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)

        dder["unit_cell_parameters"] = numpy.stack([
            der_11, der_22, der_33, der_12, der_13, der_23], axis=0)
    return m_reciprocal_g_ij, dder


def calc_m_reciprocal_g_norm_by_unit_cell_parameters(
        unit_cell_parameters, flag_unit_cell_parameters: bool = False):
    a, b = unit_cell_parameters[0], unit_cell_parameters[1]
    c = unit_cell_parameters[2]
    alpha, beta = unit_cell_parameters[3], unit_cell_parameters[4]
    gamma = unit_cell_parameters[5]
    c_a, c_b, c_g = numpy.cos(alpha), numpy.cos(beta), numpy.cos(gamma)
    s_a, s_b, s_g = numpy.sin(alpha), numpy.sin(beta), numpy.sin(gamma)
    phi_sq, dder_phi_sq = calc_phi_sq_by_unit_cell_parameters(
        unit_cell_parameters, flag_unit_cell_parameters=flag_unit_cell_parameters) 

    ones = numpy.ones_like(c_a)

    inv_g_norm_11 = ones
    inv_g_norm_22 = ones
    inv_g_norm_33 = ones
    inv_g_norm_12 = (c_a*c_b-c_g) / (s_a*s_b)
    inv_g_norm_13 = (c_a*c_g-c_b) / (s_a*s_g)
    inv_g_norm_23 = (c_b*c_g-c_a) / (s_b*s_g)

    # G matrix for reciprocal unit cell
    m_reciprocal_g_norm_ij = numpy.stack(
        [inv_g_norm_11, inv_g_norm_22, inv_g_norm_33,
         inv_g_norm_12, inv_g_norm_13, inv_g_norm_23],
        axis=0)
    dder = {}
    if flag_unit_cell_parameters:
        zero = numpy.zeros_like(a)
        s_a_sq = numpy.square(s_a)
        s_b_sq = numpy.square(s_b)
        s_g_sq = numpy.square(s_g)

        der_11 = numpy.stack([zero, zero, zero, zero, zero, zero], axis=0)
        der_22 = numpy.stack([zero, zero, zero, zero, zero, zero], axis=0)
        der_33 = numpy.stack([zero, zero, zero, zero, zero, zero], axis=0)
        der_12 = numpy.stack([zero, zero, zero, -(c_a*c_b - c_g)*c_a/(s_a_sq*s_b) - c_b/s_b, -(c_a*c_b - c_g)*c_b/(s_a*s_b_sq) - c_a/s_a, s_g/(s_a*s_b)], axis=0)
        der_13 = numpy.stack([zero, zero, zero, -(c_a*c_g - c_b)*c_a/(s_a_sq*s_g) - c_g/s_g, s_b/(s_a*s_g), -(c_a*c_g - c_b)*c_g/(s_a*s_g_sq) - c_a/s_a], axis=0)
        der_23 = numpy.stack([zero, zero, zero, s_a/(s_b*s_g), -(-c_a + c_b*c_g)*c_b/(s_b_sq*s_g) - c_g/s_g, -(-c_a + c_b*c_g)*c_g/(s_b*s_g_sq) - c_b/s_b], axis=0)

        dder["unit_cell_parameters"] = numpy.stack([
            der_11, der_22, der_33, der_12, der_13, der_23], axis=0)
    return m_reciprocal_g_norm_ij, dder


def calc_m_b_by_unit_cell_parameters(
        unit_cell_parameters, flag_unit_cell_parameters: bool = False):
    a, b = unit_cell_parameters[0], unit_cell_parameters[1]
    c = unit_cell_parameters[2]
    alpha, beta = unit_cell_parameters[3], unit_cell_parameters[4]
    gamma = unit_cell_parameters[5]
    c_a, c_b, c_g = numpy.cos(alpha), numpy.cos(beta), numpy.cos(gamma)
    s_a, s_b, s_g = numpy.sin(alpha), numpy.sin(beta), numpy.sin(gamma)

    phi, dder_phi = calc_phi_by_unit_cell_parameters(
        unit_cell_parameters, flag_unit_cell_parameters=flag_unit_cell_parameters) 

    zeros = numpy.zeros_like(unit_cell_parameters[0])

    m_b_ij = numpy.stack([
        s_a/(a*phi), (c_a*c_b-c_g)/(b*phi*s_a), (c_a*c_g-c_b)/(c*phi*s_a),
              zeros, 1/(b*s_a), -c_a/(c*s_a),
              zeros,     zeros, 1./c], axis=0)
    dder = {}
    if flag_unit_cell_parameters:
        der_11_a = - m_b_ij[0]/a - m_b_ij[0]/phi*dder_phi["unit_cell_parameters"][0]
        der_11_b = - m_b_ij[0]/phi*dder_phi["unit_cell_parameters"][1]
        der_11_c = - m_b_ij[0]/phi*dder_phi["unit_cell_parameters"][2]
        der_11_al = c_a/(a*phi) - m_b_ij[0]/phi*dder_phi["unit_cell_parameters"][3]
        der_11_be = - m_b_ij[0]/phi*dder_phi["unit_cell_parameters"][4]
        der_11_ga = - m_b_ij[0]/phi*dder_phi["unit_cell_parameters"][5]
        der_11 = numpy.stack([der_11_a, der_11_b, der_11_c, der_11_al, der_11_be, der_11_ga], axis=0)

        der_12_a = - m_b_ij[1]/phi*dder_phi["unit_cell_parameters"][0]
        der_12_b = - m_b_ij[1]/b - m_b_ij[1]/phi*dder_phi["unit_cell_parameters"][1]
        der_12_c = - m_b_ij[1]/phi*dder_phi["unit_cell_parameters"][2]
        der_12_al = (c_a*c_g-c_b)/(numpy.square(s_a)*b*phi) - m_b_ij[1]/phi*dder_phi["unit_cell_parameters"][3]
        der_12_be = c_a*s_b/(b*phi*s_a) - m_b_ij[1]/phi*dder_phi["unit_cell_parameters"][4]
        der_12_ga = s_g/(b*phi*s_a) - m_b_ij[1]/phi*dder_phi["unit_cell_parameters"][5]
        der_12 = numpy.stack([der_12_a, der_12_b, der_12_c, der_12_al, der_12_be, der_12_ga], axis=0)

        der_13_a = - m_b_ij[2]/phi*dder_phi["unit_cell_parameters"][0]
        der_13_b = - m_b_ij[2]/phi*dder_phi["unit_cell_parameters"][1]
        der_13_c = -m_b_ij[2]/c - m_b_ij[2]/phi*dder_phi["unit_cell_parameters"][2]
        der_13_al = (c_a*c_b-c_g)/(numpy.square(s_a)*c*phi) - m_b_ij[2]/phi*dder_phi["unit_cell_parameters"][3]
        der_13_be = s_b/(c*phi*s_a) - m_b_ij[2]/phi*dder_phi["unit_cell_parameters"][4]
        der_13_ga = -s_g*c_a/(c*phi*s_a) - m_b_ij[2]/phi*dder_phi["unit_cell_parameters"][5]
        der_13 = numpy.stack([der_13_a, der_13_b, der_13_c, der_13_al, der_13_be, der_13_ga], axis=0)

        der_22_b =  -m_b_ij[4]/b
        der_22_al = -m_b_ij[4]*c_a/numpy.square(s_a)
        der_22 = numpy.stack([zeros, der_22_b, zeros, der_22_al, zeros, zeros], axis=0)

        der_23_c =  -m_b_ij[5]/c
        der_23_al = 1/(numpy.square(s_a)*c)
        der_23 = numpy.stack([zeros, zeros, der_23_c, der_23_al, zeros, zeros], axis=0)

        der_33_c =  -1/numpy.square(c)
        der_33 = numpy.stack([zeros, zeros, der_33_c, zeros, zeros, zeros], axis=0)

        der_21 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_31 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_32 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)

        dder["unit_cell_parameters"] = numpy.stack([
            der_11, der_12, der_13, der_21, der_22, der_23, der_31, der_32,
            der_33], axis=0)
    return m_b_ij, dder


def calc_m_inv_b_by_unit_cell_parameters(
        unit_cell_parameters, flag_unit_cell_parameters: bool = False):
    a, b = unit_cell_parameters[0], unit_cell_parameters[1]
    c = unit_cell_parameters[2]
    alpha, beta = unit_cell_parameters[3], unit_cell_parameters[4]
    gamma = unit_cell_parameters[5]
    c_a, c_b, c_g = numpy.cos(alpha), numpy.cos(beta), numpy.cos(gamma)
    s_a, s_b, s_g = numpy.sin(alpha), numpy.sin(beta), numpy.sin(gamma)

    phi, dder_phi = calc_phi_by_unit_cell_parameters(
        unit_cell_parameters, flag_unit_cell_parameters=flag_unit_cell_parameters) 

    zeros = numpy.zeros_like(unit_cell_parameters[0])

    m_inv_b_ij = numpy.stack([
        (a*phi)/s_a, a*(c_g-c_a*c_b)/s_a, a*c_b,
              zeros,               b*s_a, b*c_a,
              zeros,               zeros, c], axis=0)
    dder = {}
    if flag_unit_cell_parameters:
        der_11 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_12 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_13 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_21 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_22 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_23 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_31 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_32 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_33 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)

        dder["unit_cell_parameters"] = numpy.stack([
            der_11, der_12, der_13, der_21, der_22, der_23, der_31, der_32,
            der_33], axis=0)
    return m_inv_b_ij, dder


def calc_m_b_norm_by_unit_cell_parameters(
        unit_cell_parameters, flag_unit_cell_parameters: bool = False):
    a, b = unit_cell_parameters[0], unit_cell_parameters[1]
    c = unit_cell_parameters[2]
    alpha, beta = unit_cell_parameters[3], unit_cell_parameters[4]
    gamma = unit_cell_parameters[5]
    c_a, c_b, c_g = numpy.cos(alpha), numpy.cos(beta), numpy.cos(gamma)
    s_a, s_b, s_g = numpy.sin(alpha), numpy.sin(beta), numpy.sin(gamma)

    phi, dder_phi = calc_phi_by_unit_cell_parameters(
        unit_cell_parameters, flag_unit_cell_parameters=flag_unit_cell_parameters) 

    zeros = numpy.zeros_like(unit_cell_parameters[0])
    ones = numpy.ones_like(unit_cell_parameters[0])


    m_b_norm_ij = numpy.stack(
        [ones, (c_a*c_b-c_g)/(s_a*s_b), (c_a*c_g-c_b)/(s_a*s_g),
         zeros, phi/(s_a*s_b), -phi*c_a/(s_a*s_g),
         zeros, zeros, phi/s_g], axis=0)
 
    dder = {}
    if flag_unit_cell_parameters:
        der_11 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_12 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_13 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_21 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_22 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_23 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_31 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_32 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_33 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)

        dder["unit_cell_parameters"] = numpy.stack([
            der_11, der_12, der_13, der_21, der_22, der_23, der_31, der_32,
            der_33], axis=0)
    return m_b_norm_ij, dder


def calc_m_inv_b_norm_by_unit_cell_parameters(
        unit_cell_parameters, flag_unit_cell_parameters: bool = False):
    a, b = unit_cell_parameters[0], unit_cell_parameters[1]
    c = unit_cell_parameters[2]
    alpha, beta = unit_cell_parameters[3], unit_cell_parameters[4]
    gamma = unit_cell_parameters[5]
    c_a, c_b, c_g = numpy.cos(alpha), numpy.cos(beta), numpy.cos(gamma)
    s_a, s_b, s_g = numpy.sin(alpha), numpy.sin(beta), numpy.sin(gamma)

    phi, dder_phi = calc_phi_by_unit_cell_parameters(
        unit_cell_parameters, flag_unit_cell_parameters=flag_unit_cell_parameters) 

    zeros = numpy.zeros_like(unit_cell_parameters[0])
    ones = numpy.ones_like(unit_cell_parameters[0])


    m_inv_b_norm = numpy.stack(
        [ones, (c_g-c_a*c_b)/phi, (s_a*c_b)/phi,
         zeros, (s_a*s_b)/phi, (c_a*s_b)/phi,
         zeros, zeros, s_g/phi], axis=0)
 
    dder = {}
    if flag_unit_cell_parameters:

        der_11_1 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_12_1 = numpy.stack([zeros, zeros, zeros, s_a*c_b/phi, s_b*c_a/phi, -s_g/phi], axis=0)
        der_13_1 = numpy.stack([zeros, zeros, zeros, c_a*c_b/phi, -s_a*s_b/phi, zeros], axis=0)
        der_22_1 = numpy.stack([zeros, zeros, zeros, s_b*c_a/phi, s_a*c_b/phi, zeros], axis=0)
        der_23_1 = numpy.stack([zeros, zeros, zeros, -s_a*s_b/phi, c_a*c_b/phi, zeros], axis=0)
        der_33_1 = numpy.stack([zeros, zeros, zeros, zeros, zeros, c_g/phi], axis=0)
        der_21_1 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_31_1 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_32_1 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)

        dder_2 = -numpy.expand_dims(m_inv_b_norm/phi, axis=1)*numpy.expand_dims(dder_phi["uni_cell_parameters"], axis=0)

        dder["unit_cell_parameters"] = numpy.stack([
            der_11_1, der_12_1, der_13_1, der_21_1, der_22_1, der_23_1, der_31_1, der_32_1,
            der_33_1], axis=0) + dder_2
    return m_inv_b_norm, dder


def calc_m_m_by_unit_cell_parameters(
        unit_cell_parameters, flag_unit_cell_parameters: bool = False):
    a, b = unit_cell_parameters[0], unit_cell_parameters[1]
    c = unit_cell_parameters[2]
    alpha, beta = unit_cell_parameters[3], unit_cell_parameters[4]
    gamma = unit_cell_parameters[5]
    c_a, c_b, c_g = numpy.cos(alpha), numpy.cos(beta), numpy.cos(gamma)
    s_a, s_b, s_g = numpy.sin(alpha), numpy.sin(beta), numpy.sin(gamma)

    phi, dder_phi = calc_phi_by_unit_cell_parameters(
        unit_cell_parameters, flag_unit_cell_parameters=flag_unit_cell_parameters) 

    zeros = numpy.zeros_like(unit_cell_parameters[0])

    # M matrix
    m_m_ij = numpy.stack([
        a*phi/s_a, zeros, zeros,
        a*(c_g-c_a*c_b)/s_a, b*s_a, zeros,
        a*c_b, b*c_a, c], axis=0)
    dder = {}
    if flag_unit_cell_parameters:
        der_11 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_12 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_13 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_21 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_22 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_23 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_31 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_32 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_33 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)

        dder["unit_cell_parameters"] = numpy.stack([
            der_11, der_12, der_13, der_21, der_22, der_23, der_31, der_32,
            der_33], axis=0)
    return m_m_ij, dder


def calc_m_inv_m_by_unit_cell_parameters(
        unit_cell_parameters, flag_unit_cell_parameters: bool = False):
    a, b = unit_cell_parameters[0], unit_cell_parameters[1]
    c = unit_cell_parameters[2]
    alpha, beta = unit_cell_parameters[3], unit_cell_parameters[4]
    gamma = unit_cell_parameters[5]
    c_a, c_b, c_g = numpy.cos(alpha), numpy.cos(beta), numpy.cos(gamma)
    s_a, s_b, s_g = numpy.sin(alpha), numpy.sin(beta), numpy.sin(gamma)

    phi, dder_phi = calc_phi_by_unit_cell_parameters(
        unit_cell_parameters, flag_unit_cell_parameters=flag_unit_cell_parameters) 

    zeros = numpy.zeros_like(unit_cell_parameters[0])

    # inversed M matrix
    m_inv_m_ij = numpy.stack([
        s_a/(a*phi), zeros, zeros,
        (c_a*c_b-c_g)/(b*phi*s_a), 1/(b*s_a), zeros,
        (c_a*c_g-c_b)/(c*phi*s_a), -c_a/(c*s_a), 1/c], axis=0)
    dder = {}
    if flag_unit_cell_parameters:
        der_11 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_12 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_13 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_21 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_22 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_23 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_31 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_32 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_33 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)

        dder["unit_cell_parameters"] = numpy.stack([
            der_11, der_12, der_13, der_21, der_22, der_23, der_31, der_32,
            der_33], axis=0)
    return m_inv_m_ij, dder


def calc_m_m_norm_by_unit_cell_parameters(
        unit_cell_parameters, flag_unit_cell_parameters: bool = False):
    """nM matrix."""
    a, b = unit_cell_parameters[0], unit_cell_parameters[1]
    c = unit_cell_parameters[2]
    alpha, beta = unit_cell_parameters[3], unit_cell_parameters[4]
    gamma = unit_cell_parameters[5]
    c_a, c_b, c_g = numpy.cos(alpha), numpy.cos(beta), numpy.cos(gamma)
    s_a, s_b, s_g = numpy.sin(alpha), numpy.sin(beta), numpy.sin(gamma)

    phi, dder_phi = calc_phi_by_unit_cell_parameters(
        unit_cell_parameters, flag_unit_cell_parameters=flag_unit_cell_parameters) 

    zeros = numpy.zeros_like(unit_cell_parameters[0])
    ones = numpy.ones_like(unit_cell_parameters[0])

    m_m_norm_ij = numpy.stack([
        phi/s_a, zeros, zeros,
        (c_g-c_a*c_b)/s_a, s_a, zeros,
        c_b, c_a, ones], axis=0)

    dder = {}
    if flag_unit_cell_parameters:
        s_a_sq = numpy.square(s_a)
        c_a_sq = numpy.square(c_a)
        der_11 = numpy.stack([
            zeros, zeros, zeros,
            -0.5/phi * (2*c_a-2.*c_b*c_g) - phi*c_a/s_a_sq,
            -0.5/phi * (2*c_b-2.*c_a*c_g)*s_b/s_a,
            -0.5/phi * (2*c_g-2.*c_a*c_b)*s_g/s_a], axis=0)
        der_12 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_13 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_21 = numpy.stack([
            zeros, zeros, zeros, c_b-(c_g*c_a-c_a_sq*c_b)/s_a_sq, c_a*s_b/s_a,
            -s_g/s_a], axis=0)
        der_22 = numpy.stack([zeros, zeros, zeros, c_a, zeros, zeros], axis=0)
        der_23 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_31 = numpy.stack([zeros, zeros, zeros, zeros, -1*s_b, zeros], axis=0)
        der_32 = numpy.stack([zeros, zeros, zeros, -1*s_a, zeros, zeros], axis=0)
        der_33 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)

        dder["unit_cell_parameters"] = numpy.stack([
            der_11, der_12, der_13, der_21, der_22, der_23, der_31, der_32,
            der_33], axis=0)
    return m_m_norm_ij, dder


def calc_m_inv_m_norm_by_unit_cell_parameters(
        unit_cell_parameters, flag_unit_cell_parameters: bool = False):
    """nM matrix."""
    a, b = unit_cell_parameters[0], unit_cell_parameters[1]
    c = unit_cell_parameters[2]
    alpha, beta = unit_cell_parameters[3], unit_cell_parameters[4]
    gamma = unit_cell_parameters[5]
    c_a, c_b, c_g = numpy.cos(alpha), numpy.cos(beta), numpy.cos(gamma)
    s_a, s_b, s_g = numpy.sin(alpha), numpy.sin(beta), numpy.sin(gamma)

    phi, dder_phi = calc_phi_by_unit_cell_parameters(
        unit_cell_parameters, flag_unit_cell_parameters=flag_unit_cell_parameters) 

    zeros = numpy.zeros_like(unit_cell_parameters[0])
    ones = numpy.ones_like(unit_cell_parameters[0])

    m_inv_m_norm_ij = numpy.stack([
        s_a/phi, zeros, zeros,
        (c_a*c_b-c_g)/(phi*s_a), 1./s_a, zeros,
        (c_a*c_g-c_b)/(phi*s_a), -c_a/s_a, ones], axis=0)

    dder = {}
    if flag_unit_cell_parameters:
        der_11 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_12 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_13 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_21 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_22 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_23 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_31 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_32 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_33 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)

        dder["unit_cell_parameters"] = numpy.stack([
            der_11, der_12, der_13, der_21, der_22, der_23, der_31, der_32,
            der_33], axis=0)
    return m_inv_m_norm_ij, dder


def calc_q_ccs_by_unit_cell_parameters(
        index_hkl, unit_cell_parameters,
        flag_unit_cell_parameters: bool = False):
    """
    Calculate scattering vector in Cartesian coordinate system.
    Coordinate system is defined as (x||a*, z||c).

    Keyword arguments:

        index_hkl: Miller indices (3, refln)

    Output arguments:

        k_x, k_y, k_z: 1D numpy array of x, y, z components of unity
                       scattering vector
    """
    m_b, dder_m_b = calc_m_b_by_unit_cell_parameters(
        unit_cell_parameters, flag_unit_cell_parameters = flag_unit_cell_parameters)
    q_xyz = numpy.stack(
        [m_b[0]*index_hkl[0]+m_b[1]*index_hkl[1]+m_b[2]*index_hkl[2],
         m_b[3]*index_hkl[0]+m_b[4]*index_hkl[1]+m_b[5]*index_hkl[2],
         m_b[6]*index_hkl[0]+m_b[7]*index_hkl[1]+m_b[8]*index_hkl[2]],
        axis=0)

    dder = {}
    if flag_unit_cell_parameters:
        dder["unit_cell_parameters"] = numpy.stack([
            numpy.expand_dims(index_hkl[0], axis=-1)*numpy.expand_dims(dder_m_b["unit_cell_parameters"][0], axis=0) + 
            numpy.expand_dims(index_hkl[1], axis=-1)*numpy.expand_dims(dder_m_b["unit_cell_parameters"][1], axis=0) +
            numpy.expand_dims(index_hkl[2], axis=-1)*numpy.expand_dims(dder_m_b["unit_cell_parameters"][2], axis=0),
            numpy.expand_dims(index_hkl[0], axis=-1)*numpy.expand_dims(dder_m_b["unit_cell_parameters"][3], axis=0) + 
            numpy.expand_dims(index_hkl[1], axis=-1)*numpy.expand_dims(dder_m_b["unit_cell_parameters"][4], axis=0) +
            numpy.expand_dims(index_hkl[2], axis=-1)*numpy.expand_dims(dder_m_b["unit_cell_parameters"][5], axis=0),
            numpy.expand_dims(index_hkl[0], axis=-1)*numpy.expand_dims(dder_m_b["unit_cell_parameters"][6], axis=0) + 
            numpy.expand_dims(index_hkl[1], axis=-1)*numpy.expand_dims(dder_m_b["unit_cell_parameters"][7], axis=0) +
            numpy.expand_dims(index_hkl[2], axis=-1)*numpy.expand_dims(dder_m_b["unit_cell_parameters"][8], axis=0)
        ], axis=0)
    return q_xyz, dder


def calc_eq_ccs_by_unit_cell_parameters(
        index_hkl, unit_cell_parameters,
        flag_unit_cell_parameters: bool = False):
    """
    Calculate unity scattering vector in Cartesian coordinate system.
    Coordinate system is defined as (x||a*, z||c).

    Keyword arguments:

        index_hkl: Miller indices (3, refln)

    Output arguments:

        k_x, k_y, k_z: 1D numpy array of x, y, z components of unity
                       scattering vector
    """
    q_xyz, dder_qxyz = calc_q_ccs_by_unit_cell_parameters(
        index_hkl, unit_cell_parameters, flag_unit_cell_parameters=flag_unit_cell_parameters)

    q_norm = numpy.sqrt(numpy.sum(numpy.square(q_xyz), axis=0))

    q_norm = numpy.where(q_norm != 0., q_norm, 1.)

    eq_xyz = q_xyz/q_norm
    dder = {}
    if flag_unit_cell_parameters:
        q_norm_exp = numpy.expand_dims(numpy.expand_dims(q_norm, axis=0), axis=-1)
        q_xyz_exp = numpy.expand_dims(q_xyz, axis=-1)
        hh = (q_xyz_exp*dder_qxyz["unit_cell_parameters"]).sum(axis=0)
        hh_exp = numpy.expand_dims(hh, axis=0)
        dder["unit_cell_parameters"] = (
            dder_qxyz["unit_cell_parameters"]/q_norm_exp-
            q_xyz_exp/numpy.square(q_norm_exp)*hh_exp/q_norm_exp)
    return eq_xyz, dder


def calc_inv_d_by_unit_cell_parameters(
        index_hkl, unit_cell_parameters,
        flag_unit_cell_parameters: bool = False):
    """
    Calculate 1/d for given reflections hkl
    and unit cell parameters defined as [a, b, c, alpha, beta, gamma].
    """
    a, b = unit_cell_parameters[0], unit_cell_parameters[1]
    c = unit_cell_parameters[2]
    c_a = numpy.cos(unit_cell_parameters[3])
    c_b = numpy.cos(unit_cell_parameters[4])
    c_g = numpy.cos(unit_cell_parameters[5])
    c_a_sq, c_b_sq = numpy.square(c_a), numpy.square(c_b)
    c_g_sq = numpy.square(c_g)
    s_a_sq, s_b_sq, s_g_sq = (1.-c_a_sq), (1.-c_b_sq), (1.-c_g_sq)

    A = (1.-c_a_sq-c_b_sq-c_g_sq+2.*c_a*c_b*c_g)
    h = index_hkl[0]
    k = index_hkl[1]
    l = index_hkl[2]
    B1 = s_a_sq*numpy.square(h*1./a)+\
         s_b_sq*numpy.square(k*1./b)+\
         s_g_sq*numpy.square(l*1./c)
    B2 = 2.*(k*l*(c_a-c_b*c_g))/(b*c)+\
         2.*(h*l*(c_b-c_a*c_g))/(a*c)+\
         2.*(h*k*(c_g-c_a*c_b))/(a*b)
    B = B1-B2
    inv_d = numpy.sqrt(B*1./A)
    dder = {}
    if flag_unit_cell_parameters:
        coeff = 1./(inv_d * A)
        s_a = numpy.sin(unit_cell_parameters[3])
        s_b = numpy.sin(unit_cell_parameters[4])
        s_g = numpy.sin(unit_cell_parameters[5])
        dder_a = coeff * (h/numpy.square(a)) * ((c_b*l)/c + (c_g*k)/b - (s_a_sq*h)/a)*numpy.ones_like(a)
        dder_b = coeff * (k/numpy.square(b)) * ((c_g*h)/a + (c_a*l)/c - (s_b_sq*k)/b)*numpy.ones_like(b)
        dder_c = coeff * (l/numpy.square(c)) * ((c_a*k)/b + (c_b*h)/a - (s_g_sq*l)/c)*numpy.ones_like(c)
        dder_al = (s_a/A) * ((c_a*numpy.square(h/a) + (k*l)/(b*c))/inv_d - inv_d*(c_a - c_b * c_g)/A)*numpy.ones_like(c_a)
        dder_be = (s_b/A) * ((c_b*numpy.square(k/b) + (l*h)/(c*a))/inv_d - inv_d*(c_b - c_g * c_a)/A)*numpy.ones_like(c_b)
        dder_ga = (s_g/A) * ((c_g*numpy.square(l/c) + (h*k)/(a*b))/inv_d - inv_d*(c_g - c_a * c_b)/A)*numpy.ones_like(c_g)
        dder["unit_cell_parameters"] = numpy.stack([dder_a, dder_b, dder_c, dder_al, dder_be, dder_ga], axis=0)
    return inv_d, dder


def calc_sthovl_by_unit_cell_parameters(index_hkl, unit_cell_parameters,
        flag_unit_cell_parameters: bool = False):
    """
    Calculate sin(theta)/lambda for given reflections h, k, l
    and unit cell parameters defined as a, b, c, cos(alpha), cos(beta), cos(gamma).
    """
    inv_d, dder_inv_d = calc_inv_d_by_unit_cell_parameters(
        index_hkl, unit_cell_parameters,
        flag_unit_cell_parameters=flag_unit_cell_parameters)
    sthovl = 0.5*inv_d
    dder = {}
    if flag_unit_cell_parameters:
        dder["unit_cell_parameters"] = 0.5*dder_inv_d["unit_cell_parameters"]
    return sthovl, dder



def transform_quadratic_form_reciprocal_to_ccs(
    quadratic_form_reciprocal, unit_cell_parameters,
    flag_quadratic_form_reciprocal:bool=False,
    flag_unit_cell_parameters:bool=False):
    """Transform quadratic form from reciprocal system
    ($a^*$, $b^*$, $c^*$) to Cartesian coordinate system (x||a*, z||c)
    """
    m_m_ij, dder_minvb = calc_m_inv_b_by_unit_cell_parameters(
        unit_cell_parameters, flag_unit_cell_parameters=flag_unit_cell_parameters)
    m_mt_ij = numpy.stack([
        m_m_ij[0], m_m_ij[3], m_m_ij[6],
        m_m_ij[1], m_m_ij[4], m_m_ij[7],
        m_m_ij[2], m_m_ij[5], m_m_ij[8]], axis=0)
    quadratic_form_ccs, dder_q = calc_m_q_mt(m_mt_ij, quadratic_form_reciprocal,
        flag_m=flag_unit_cell_parameters, flag_q=flag_quadratic_form_reciprocal)
    dder = {}
    if flag_quadratic_form_reciprocal:
        dder["quadratic_form_reciprocal"] = numpy.zeros((6, )+quadratic_form_reciprocal.shape, dtype=quadratic_form_reciprocal.dtype)
    if flag_unit_cell_parameters:
        dder["unit_cell_parameters"] = numpy.zeros((6, )+unit_cell_parameters.shape, dtype=unit_cell_parameters.dtype)
    return quadratic_form_ccs, dder


def transform_quadratic_form_ccs_to_reciprocal(
    quadratic_form_ccs, unit_cell_parameters,
    flag_quadratic_form_ccs:bool=False,
    flag_unit_cell_parameters:bool=False):
    """Transform quadratic form from Cartesian coordinate system (x||a*, z||c)
   to reciprocal system ($a^*$, $b^*$, $c^*$)
    """
    m_m_ij, dder = calc_m_b_by_unit_cell_parameters(
        unit_cell_parameters, flag_unit_cell_parameters=False)
    m_mt_ij = numpy.stack([
        m_m_ij[0], m_m_ij[3], m_m_ij[6],
        m_m_ij[1], m_m_ij[4], m_m_ij[7],
        m_m_ij[2], m_m_ij[5], m_m_ij[8]], axis=0)
    quadratic_form_reciprocal = calc_m_q_mt(m_mt_ij, quadratic_form_ccs, flag_m=False, flag_q=False)[0]
    dder = {}
    if flag_quadratic_form_ccs:
        dder["quadratic_form_ccs"] = numpy.zeros((6, )+quadratic_form_ccs.shape, dtype=quadratic_form_ccs.dtype)
    if flag_unit_cell_parameters:
        dder["unit_cell_parameters"] = numpy.zeros((6, )+unit_cell_parameters.shape, dtype=unit_cell_parameters.dtype)
    return quadratic_form_reciprocal, dder


def transform_quadratic_form_reciprocal_norm_to_ccs(
    quadratic_form_reciprocal_norm, unit_cell_parameters,
    flag_quadratic_form_reciprocal_norm:bool=False,
    flag_unit_cell_parameters:bool=False):
    """Transform quadratic form from normalized reciprocal system
    ($a^*/a^*$, $b^*/b^*$, $c^*/c^*$) to Cartesian coordinate system (x||a*, z||c)
    """
    m_m_ij, dder = calc_m_inv_b_norm_by_unit_cell_parameters(
        unit_cell_parameters, flag_unit_cell_parameters=False)
    m_mt_ij = numpy.stack([
        m_m_ij[0], m_m_ij[3], m_m_ij[6],
        m_m_ij[1], m_m_ij[4], m_m_ij[7],
        m_m_ij[2], m_m_ij[5], m_m_ij[8]], axis=0)
    quadratic_form_ccs = calc_m_q_mt(m_mt_ij, quadratic_form_reciprocal_norm, flag_m=False, flag_q=False)[0]
    dder = {}
    if flag_quadratic_form_reciprocal_norm:
        dder["quadratic_form_reciprocal_norm"] = numpy.zeros((6, )+quadratic_form_reciprocal_norm.shape, dtype=quadratic_form_reciprocal_norm.dtype)
    if flag_unit_cell_parameters:
        dder["unit_cell_parameters"] = numpy.zeros((6, )+unit_cell_parameters.shape, dtype=unit_cell_parameters.dtype)
    return quadratic_form_ccs, dder


def transform_quadratic_form_ccs_to_reciprocal_norm(
    quadratic_form_ccs, unit_cell_parameters,
    flag_quadratic_form_ccs:bool=False,
    flag_unit_cell_parameters:bool=False):
    """Transform quadratic form from Cartesian coordinate system (x||a*, z||c)
   to normalized reciprocal system ($a^*/a^*$, $b^*/b^*$, $c^*/c^*$)
    """
    m_m_ij, dder = calc_m_b_norm_by_unit_cell_parameters(
        unit_cell_parameters, flag_unit_cell_parameters=False)
    m_mt_ij = numpy.stack([
        m_m_ij[0], m_m_ij[3], m_m_ij[6],
        m_m_ij[1], m_m_ij[4], m_m_ij[7],
        m_m_ij[2], m_m_ij[5], m_m_ij[8]], axis=0)
    quadratic_form_reciprocal_norm = calc_m_q_mt(m_mt_ij, quadratic_form_ccs, flag_m=False, flag_q=False)[0]
    dder = {}
    if flag_quadratic_form_ccs:
        dder["quadratic_form_ccs"] = numpy.zeros((6, )+quadratic_form_ccs.shape, dtype=quadratic_form_ccs.dtype)
    if flag_unit_cell_parameters:
        dder["unit_cell_parameters"] = numpy.zeros((6, )+unit_cell_parameters.shape, dtype=unit_cell_parameters.dtype)
    return quadratic_form_reciprocal_norm, dder


def transform_quadratic_form_direct_to_ccs(
    quadratic_form_direct, unit_cell_parameters,
    flag_quadratic_form_direct:bool=False,
    flag_unit_cell_parameters:bool=False):
    """Transform quadratic form from direct system
    ($a$, $b$, $c$) to Cartesian coordinate system (x||a*, z||c)
    """
    m_m_ij, dder = calc_m_inv_m_by_unit_cell_parameters(
        unit_cell_parameters, flag_unit_cell_parameters=False)
    m_mt_ij = numpy.stack([
        m_m_ij[0], m_m_ij[3], m_m_ij[6],
        m_m_ij[1], m_m_ij[4], m_m_ij[7],
        m_m_ij[2], m_m_ij[5], m_m_ij[8]], axis=0)
    quadratic_form_ccs = calc_m_q_mt(m_mt_ij, quadratic_form_direct, flag_m=False, flag_q=False)[0]
    dder = {}
    if flag_quadratic_form_direct:
        dder["quadratic_form_direct"] = numpy.zeros((6, )+quadratic_form_direct.shape, dtype=quadratic_form_direct.dtype)
    if flag_unit_cell_parameters:
        dder["unit_cell_parameters"] = numpy.zeros((6, )+unit_cell_parameters.shape, dtype=unit_cell_parameters.dtype)
    return quadratic_form_ccs, dder


def transform_quadratic_form_ccs_to_direct(
    quadratic_form_ccs, unit_cell_parameters,
    flag_quadratic_form_ccs:bool=False,
    flag_unit_cell_parameters:bool=False):
    """Transform quadratic form from Cartesian coordinate system (x||a*, z||c)
   to direct system ($a$, $b$, $c$)
    """
    m_m_ij, dder = calc_m_m_by_unit_cell_parameters(
        unit_cell_parameters, flag_unit_cell_parameters=False)
    m_mt_ij = numpy.stack([
        m_m_ij[0], m_m_ij[3], m_m_ij[6],
        m_m_ij[1], m_m_ij[4], m_m_ij[7],
        m_m_ij[2], m_m_ij[5], m_m_ij[8]], axis=0)
    quadratic_form_direct = calc_m_q_mt(m_mt_ij, quadratic_form_ccs, flag_m=False, flag_q=False)[0]
    dder = {}
    if flag_quadratic_form_ccs:
        dder["quadratic_form_ccs"] = numpy.zeros((6, )+quadratic_form_ccs.shape, dtype=quadratic_form_ccs.dtype)
    if flag_unit_cell_parameters:
        dder["unit_cell_parameters"] = numpy.zeros((6, )+unit_cell_parameters.shape, dtype=unit_cell_parameters.dtype)
    return quadratic_form_direct, dder


def transform_quadratic_form_direct_norm_to_ccs(
    quadratic_form_direct_norm, unit_cell_parameters,
    flag_quadratic_form_direct_norm:bool=False,
    flag_unit_cell_parameters:bool=False):
    """Transform quadratic form from normalized direct system
    ($a/a$, $b/b$, $c/c$) to Cartesian coordinate system (x||a*, z||c)
    """
    m_m_ij, dder = calc_m_inv_m_norm_by_unit_cell_parameters(
        unit_cell_parameters, flag_unit_cell_parameters=False)
    m_mt_ij = numpy.stack([
        m_m_ij[0], m_m_ij[3], m_m_ij[6],
        m_m_ij[1], m_m_ij[4], m_m_ij[7],
        m_m_ij[2], m_m_ij[5], m_m_ij[8]], axis=0)
    quadratic_form_ccs = calc_m_q_mt(m_mt_ij, quadratic_form_direct_norm, flag_m=False, flag_q=False)[0]
    dder = {}
    if flag_quadratic_form_direct_norm:
        dder["quadratic_form_direct_norm"] = numpy.zeros((6, )+quadratic_form_direct_norm.shape, dtype=quadratic_form_direct_norm.dtype)
    if flag_unit_cell_parameters:
        dder["unit_cell_parameters"] = numpy.zeros((6, )+unit_cell_parameters.shape, dtype=unit_cell_parameters.dtype)
    return quadratic_form_ccs, dder


def transform_quadratic_form_ccs_to_direct_norm(
    quadratic_form_ccs, unit_cell_parameters,
    flag_quadratic_form_ccs:bool=False,
    flag_unit_cell_parameters:bool=False):
    """Transform quadratic form from Cartesian coordinate system (x||a*, z||c)
   to normalized direct system ($a/a$, $b/b$, $c/c$)
    """
    m_m_ij, dder = calc_m_m_norm_by_unit_cell_parameters(
        unit_cell_parameters, flag_unit_cell_parameters=False)
    m_mt_ij = numpy.stack([
        m_m_ij[0], m_m_ij[3], m_m_ij[6],
        m_m_ij[1], m_m_ij[4], m_m_ij[7],
        m_m_ij[2], m_m_ij[5], m_m_ij[8]], axis=0)
    quadratic_form_direct_norm = calc_m_q_mt(m_mt_ij, quadratic_form_ccs, flag_m=False, flag_q=False)[0]
    dder = {}
    if flag_quadratic_form_ccs:
        dder["quadratic_form_ccs"] = numpy.zeros((6, )+quadratic_form_ccs.shape, dtype=quadratic_form_ccs.dtype)
    if flag_unit_cell_parameters:
        dder["unit_cell_parameters"] = numpy.zeros((6, )+unit_cell_parameters.shape, dtype=unit_cell_parameters.dtype)
    return quadratic_form_direct_norm, dder

def calc_m_inv_m_b_norm_by_unit_cell_parameters(
        unit_cell_parameters, flag_unit_cell_parameters: bool = False):
    """
    Calculation M^-1 B_norm
    """
    a, b, c = unit_cell_parameters[0], unit_cell_parameters[1], unit_cell_parameters[2]
    al, be, ga = unit_cell_parameters[3], unit_cell_parameters[4], unit_cell_parameters[5]
    c_a, c_b, c_g = numpy.cos(al), numpy.cos(be), numpy.cos(ga)
    s_a, s_b, s_g = numpy.sin(al), numpy.sin(be), numpy.sin(ga)
    phi, dder_phi = calc_phi_by_unit_cell_parameters(
        unit_cell_parameters,
        flag_unit_cell_parameters=flag_unit_cell_parameters)

    o_11 = s_a / (a * phi)
    o_12 = (c_a * c_b - c_g) / (a * phi * s_b)
    o_13 = (c_a * c_g - c_b) / (a * phi * s_g)
    o_21 = (c_a * c_b - c_g) / (b * phi * s_a)
    o_22 = s_b / (b * phi)
    o_23 = (c_b * c_g - c_a) / (b * phi * s_g)
    o_31 = (c_a * c_g - c_b) / (c * phi * s_a)
    o_32 = (c_b * c_g - c_a) / (c * phi * s_b)
    o_33 = s_g / (c * phi)

    inv_m_b_norm = numpy.stack([
        o_11, o_12, o_13, o_21, o_22, o_23, o_31, o_32, o_33],
        axis=0)
    dder = {}
    if flag_unit_cell_parameters:
        zeros = numpy.zeros_like(a)
        der_11 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_12 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_13 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_21 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_22 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_23 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_31 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_32 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_33 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        dder["unit_cell_parameters"] = numpy.stack([
            der_11, der_12, der_13, der_21, der_22, der_23, der_31, der_32, der_33], axis=0)
    return inv_m_b_norm, dder

def calc_m_inv_b_norm_m_by_unit_cell_parameters(
        unit_cell_parameters, flag_unit_cell_parameters: bool = False):
    """
    Calculation B^-1_norm M
    """
    a, b, c = unit_cell_parameters[0], unit_cell_parameters[1], unit_cell_parameters[2]
    al, be, ga = unit_cell_parameters[3], unit_cell_parameters[4], unit_cell_parameters[5]
    c_a, c_b, c_g = numpy.cos(al), numpy.cos(be), numpy.cos(ga)
    s_a, s_b, s_g = numpy.sin(al), numpy.sin(be), numpy.sin(ga)
    phi, dder_phi = calc_phi_by_unit_cell_parameters(
        unit_cell_parameters,
        flag_unit_cell_parameters=flag_unit_cell_parameters)
    o_11 = a * s_a / phi
    o_12 = b * s_a * c_g / phi
    o_13 = c * s_a * c_b / phi
    o_21 = a * s_b * c_g / phi
    o_22 = b * s_b / phi
    o_23 = c * s_b * c_a / phi
    o_31 = a * s_g * c_b / phi
    o_32 = b * s_g * c_a / phi
    o_33 = c * s_g / phi

    inv_b_norm_m = numpy.stack([
        o_11, o_12, o_13, o_21, o_22, o_23, o_31, o_32, o_33],
        axis=0)
    dder = {}
    if flag_unit_cell_parameters:
        zeros = numpy.zeros_like(a)
        der_11 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_12 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_13 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_21 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_22 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_23 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_31 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_32 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        der_33 = numpy.stack([zeros, zeros, zeros, zeros, zeros, zeros], axis=0)
        dder["unit_cell_parameters"] = numpy.stack([
            der_11, der_12, der_13, der_21, der_22, der_23, der_31, der_32, der_33], axis=0)
    return inv_b_norm_m, dder


def calc_matrix_t(index_hkl, unit_cell_parameters, flag_unit_cell_parameters: bool = False):
    """Calculate matrix T to have new z axis along scattering vector q = h a* + k b* + l c* from CCS in which (x||a*, z||c).

    m is vector in Cartesian coordinate system (x||a*, z||c)
    M is vector in Cartesian coordinate system (Z||q, axes X and Y are not important)

    M = T * m
    """
    eq_ccs, dder_eq_ccs = calc_eq_ccs_by_unit_cell_parameters(
        index_hkl, unit_cell_parameters, flag_unit_cell_parameters=flag_unit_cell_parameters)

    k_x = eq_ccs[0]
    k_y = eq_ccs[1]
    k_z = eq_ccs[2]
    
    al = numpy.zeros_like(k_x, dtype=float)
    be = numpy.arccos(k_z)
    sb = numpy.sin(be)
    flag = (sb != 0.)
    sa1 = k_x[flag]*1./sb[flag]
    ca2 = -1*k_y[flag]*1./sb[flag]
    sa1[sa1 > 1] = 1.
    sa1[sa1 < -1] = -1.
    ca2[ca2 > 1] = 1.
    ca2[ca2 < -1] = -1.
    al1 = numpy.arcsin(sa1)
    al2 = numpy.arccos(ca2)
    al_sh = numpy.copy(al1)
    al_sh[sa1 > 0.] = al2[sa1 > 0.]
    al_sh[sa1 <= 0.] = 2.*numpy.pi-al2[sa1 <= 0.]
    al_sh[numpy.abs(al2-al1) < 0.00001] = al1[numpy.abs(al2-al1) < 0.00001]
    al[flag] = al_sh
    ga = 0.
    ca, cb, cg = numpy.cos(al), numpy.cos(be), numpy.cos(ga)
    sa, sb, sg = numpy.sin(al), numpy.sin(be), numpy.sin(ga)

    # here is taken transposed one (== inversed one) to get relation M_(Z||hkl) = T * m_(x||a*, z||c)
    t_11, t_21, t_31 = ca*cg-sa*cb*sg, -ca*sg-sa*cb*cg,  sa*sb
    t_12, t_22, t_32 = sa*cg+ca*cb*sg, -sa*sg+ca*cb*cg, -ca*sb
    t_13, t_23, t_33 = sb*sg, sb*cg, cb
    flag = (((sa*sb-k_x)**2+(-ca*sb-k_y)**2+(cb-k_z)**2) > 0.0001)
    if any(flag):
        warnings.warn("Mistake with eq\nProgram is stopped",
                      UserWarning, stacklevel=2)
        quit()
    m_t_ij  = numpy.stack([t_11, t_12, t_13, t_21, t_22, t_23, t_31, t_32, t_33], axis=0)

    dder = {}
    if flag_unit_cell_parameters:
        dder["unit_cell_parameters"] = numpy.zeros((9, ) + unit_cell_parameters.shape, dtype=unit_cell_parameters.dtype)
    return m_t_ij, dder