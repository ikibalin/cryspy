# -*- coding: utf-8 -*-
"""
Expression for local susceptibility in different coordinate systems
and magnetization ellipsoid are presented.

For details see documentation.
"""
import numpy 
import numpy.linalg
na = numpy.newaxis

from .matrix_operations import \
    calc_m_q_inv_m, calc_mt_m, calc_m1_m2, calc_m1_m2_inv_m1, calc_q_sq, calc_mt_q_m

from .unit_cell import \
    calc_m_inv_m_norm_by_unit_cell_parameters, \
    calc_m_m_norm_by_unit_cell_parameters, \
    calc_m_g_norm_by_unit_cell_parameters, \
    calc_m_m_by_unit_cell_parameters,\
    calc_m_reciprocal_g_norm_by_unit_cell_parameters, \
    calc_m_inv_b_norm_by_unit_cell_parameters
    

def calc_magnetization_ellipsoid_as_u(
        susceptibility, unit_cell_parameters, 
        flag_susceptibility: bool = False, flag_unit_cell_parameters: bool = False):
    """
    Calculate magnetization ellipsoid through U_ij parameters.
    The magnetization ellipsoid shows the magnitude of induced magnetic
    moment at given orientation of applied field.
    """
    m_inv_g_norm, dder_inv_g_norm = calc_m_reciprocal_g_norm_by_unit_cell_parameters(unit_cell_parameters, flag_unit_cell_parameters=flag_unit_cell_parameters)
    inv_b_norm, dder_inv_b_norm = calc_m_inv_b_norm_by_unit_cell_parameters(unit_cell_parameters, flag_unit_cell_parameters=flag_unit_cell_parameters)

    flag_hh_1 = flag_unit_cell_parameters or flag_susceptibility
    hh_1, dder_hh_1 = calc_m_q_inv_m(inv_b_norm, susceptibility, flag_m=flag_unit_cell_parameters, flag_q=flag_susceptibility)
    res, dder_res = calc_mt_q_m(hh_1, m_inv_g_norm, flag_m=flag_hh_1, flag_q=flag_susceptibility)
    
    dder = {}
    if flag_susceptibility:
        dder["susceptibility"] = numpy.stack([
            numpy.zeros_like(susceptibility),
            numpy.zeros_like(susceptibility),
            numpy.zeros_like(susceptibility),
            numpy.zeros_like(susceptibility),
            numpy.zeros_like(susceptibility),
            numpy.zeros_like(susceptibility),
            numpy.zeros_like(susceptibility),
            numpy.zeros_like(susceptibility),
            numpy.zeros_like(susceptibility)], axis=0)
    if flag_unit_cell_parameters:
        dder["unit_cell_parameters"] = numpy.stack([
            numpy.zeros_like(unit_cell_parameters),
            numpy.zeros_like(unit_cell_parameters),
            numpy.zeros_like(unit_cell_parameters),
            numpy.zeros_like(unit_cell_parameters),
            numpy.zeros_like(unit_cell_parameters),
            numpy.zeros_like(unit_cell_parameters),
            numpy.zeros_like(unit_cell_parameters),
            numpy.zeros_like(unit_cell_parameters),
            numpy.zeros_like(unit_cell_parameters)], axis=0)
    return res, dder


def calc_magnetization_ellipsoid_axes(susceptibility, unit_cell_parameters):
    """
    Calculate eigenvalues and eigenvectors of magnetization ellipsoid
    The magnetization ellipsoid shows the magnitude of induced magnetic
    moment at given orientation of applied field.

    eigenvalues is modulus of induced magnetic moment at field applied along main axes 

    eigen_fields = (eigenvectors_1i, eigenvectors_2i, eigenvectors_3i)

    eigen_moments = chi * eigen_fields

    in Cartesian coordinate system (X||a*, Z||c).
    """
    q_ij = calc_q_sq(susceptibility, flag_q=False)[0]
    np_q_ccs = numpy.array([
        [q_ij[0], q_ij[3], q_ij[4]],
        [q_ij[3], q_ij[1], q_ij[5]],
        [q_ij[4], q_ij[5], q_ij[2]]], dtype=float)
    
    val, eig_fields = numpy.linalg.eigh(np_q_ccs) # np_q_ccs should have quadratic form for eigh

    moments_pos = numpy.sqrt(val)
    np_chi_ccs = numpy.array([
        [susceptibility[0], susceptibility[3], susceptibility[4]],
        [susceptibility[3], susceptibility[1], susceptibility[5]],
        [susceptibility[4], susceptibility[5], susceptibility[2]]], dtype=float)
    eig_moments = numpy.matmul(np_chi_ccs, eig_fields)

    moments_sign = 2*((eig_fields*eig_moments).sum(axis=0) >= 0)-1
    moments = moments_pos * moments_sign
    return moments, eig_fields, eig_moments


def calc_ellipsoid_factor(susceptibility, unit_cell_parameters):
    """
    Calculate ellipsoid factor: a/b - b/c 
    a <= b <= c

    see doi: 10.3389/fendo.2015.00015
    """
    val = calc_magnetization_ellipsoid_axes(susceptibility, unit_cell_parameters)[0]
    val_sort = numpy.sort(val)

    if val_sort[2] == 0.:
        res = 0.
    elif val_sort[1] == 0.:
        res = 1. - (val_sort[1]/val_sort[2])
    else:
        res = (val_sort[0]/val_sort[1]) - (val_sort[1]/val_sort[2])
    return res 


def calc_chi_direct_norm_with_symmetry(susceptibility, symm_elems_r, flag_susceptibility: bool = False):
    """
    Calculate susceptibility tensor in direct space after symmetry operation:
    Rs chi_direct norm Rs^-1
    """
    flag_m2 = flag_susceptibility
    chi_s_direct, dder = calc_m1_m2_inv_m1(symm_elems_r, susceptibility, flag_m1=False, flag_m2=flag_m2)
    dder = {}
    if flag_susceptibility:
        if "m2" in dder.keys():
            dder["susceptibility"] = dder.pop("m2")
        else:
            dder["susceptibility_real"] = dder.pop("m2_real")
            dder["susceptibility_imag"] = dder.pop("m2_imag")
    return chi_s_direct, dder


# def calc_sc_chi(r_ccs):
#     """Calculate symmetry constraint matrix for susceptibility
#     It is supposed that determinant of r matrix is 1
#     """
#     r_11, r_12, r_13 = r_ccs[0], r_ccs[1], r_ccs[2]
#     r_21, r_22, r_23 = r_ccs[3], r_ccs[4], r_ccs[5]
#     r_31, r_32, r_33 = r_ccs[6], r_ccs[7], r_ccs[8]
#     det_r =  r_11*r_22*r_33 - r_11*r_23*r_32  - r_12*r_21*r_33 + r_12*r_23*r_31  + r_13*r_21*r_32 - r_13*r_22*r_31
#     
#     m_mn = numpy.array([
#         [r_11*(r_22*r_33 - r_23*r_32), r_12*(-r_21*r_33 + r_23*r_31), r_13*(r_21*r_32 - r_22*r_31), r_11*(-r_21*r_33 + r_23*r_31) + r_12*(r_22*r_33 - r_23*r_32), r_11*(r_21*r_32 - r_22*r_31) + r_13*(r_22*r_33 - r_23*r_32), r_12*(r_21*r_32 - r_22*r_31) + r_13*(-r_21*r_33 + r_23*r_31)],
#         [r_11*(-r_12*r_33 + r_13*r_32), r_12*(r_11*r_33 - r_13*r_31), r_13*(-r_11*r_32 + r_12*r_31), r_11*(r_11*r_33 - r_13*r_31) + r_12*(-r_12*r_33 + r_13*r_32), r_11*(-r_11*r_32 + r_12*r_31) + r_13*(-r_12*r_33 + r_13*r_32), r_12*(-r_11*r_32 + r_12*r_31) + r_13*(r_11*r_33 - r_13*r_31)],
#         [r_11*(r_12*r_23 - r_13*r_22), r_12*(-r_11*r_23 + r_13*r_21), r_13*(r_11*r_22 - r_12*r_21), r_11*(-r_11*r_23 + r_13*r_21) + r_12*(r_12*r_23 - r_13*r_22), r_11*(r_11*r_22 - r_12*r_21) + r_13*(r_12*r_23 - r_13*r_22), r_12*(r_11*r_22 - r_12*r_21) + r_13*(-r_11*r_23 + r_13*r_21)],
#         [r_21*(r_22*r_33 - r_23*r_32), r_22*(-r_21*r_33 + r_23*r_31), r_23*(r_21*r_32 - r_22*r_31), r_21*(-r_21*r_33 + r_23*r_31) + r_22*(r_22*r_33 - r_23*r_32), r_21*(r_21*r_32 - r_22*r_31) + r_23*(r_22*r_33 - r_23*r_32), r_22*(r_21*r_32 - r_22*r_31) + r_23*(-r_21*r_33 + r_23*r_31)],
#         [r_21*(-r_12*r_33 + r_13*r_32), r_22*(r_11*r_33 - r_13*r_31), r_23*(-r_11*r_32 + r_12*r_31), r_21*(r_11*r_33 - r_13*r_31) + r_22*(-r_12*r_33 + r_13*r_32), r_21*(-r_11*r_32 + r_12*r_31) + r_23*(-r_12*r_33 + r_13*r_32), r_22*(-r_11*r_32 + r_12*r_31) + r_23*(r_11*r_33 - r_13*r_31)],
#         [r_21*(r_12*r_23 - r_13*r_22), r_22*(-r_11*r_23 + r_13*r_21), r_23*(r_11*r_22 - r_12*r_21), r_21*(-r_11*r_23 + r_13*r_21) + r_22*(r_12*r_23 - r_13*r_22), r_21*(r_11*r_22 - r_12*r_21) + r_23*(r_12*r_23 - r_13*r_22), r_22*(r_11*r_22 - r_12*r_21) + r_23*(-r_11*r_23 + r_13*r_21)],
#         [r_31*(r_22*r_33 - r_23*r_32), r_32*(-r_21*r_33 + r_23*r_31), r_33*(r_21*r_32 - r_22*r_31), r_31*(-r_21*r_33 + r_23*r_31) + r_32*(r_22*r_33 - r_23*r_32), r_31*(r_21*r_32 - r_22*r_31) + r_33*(r_22*r_33 - r_23*r_32), r_32*(r_21*r_32 - r_22*r_31) + r_33*(-r_21*r_33 + r_23*r_31)],
#         [r_31*(-r_12*r_33 + r_13*r_32), r_32*(r_11*r_33 - r_13*r_31), r_33*(-r_11*r_32 + r_12*r_31), r_31*(r_11*r_33 - r_13*r_31) + r_32*(-r_12*r_33 + r_13*r_32), r_31*(-r_11*r_32 + r_12*r_31) + r_33*(-r_12*r_33 + r_13*r_32), r_32*(-r_11*r_32 + r_12*r_31) + r_33*(r_11*r_33 - r_13*r_31)],
#         [r_31*(r_12*r_23 - r_13*r_22), r_32*(-r_11*r_23 + r_13*r_21), r_33*(r_11*r_22 - r_12*r_21), r_31*(-r_11*r_23 + r_13*r_21) + r_32*(r_12*r_23 - r_13*r_22), r_31*(r_11*r_22 - r_12*r_21) + r_33*(r_12*r_23 - r_13*r_22), r_32*(r_11*r_22 - r_12*r_21) + r_33*(-r_11*r_23 + r_13*r_21)]],
#         dtype=r_ccs.dtype)
# 
#     flag = numpy.all(numpy.isclose(numpy.abs(det_r),1))
#     if flag:
#         res = m_mn*det_r[na, na, :]
#     else:
#         res = m_mn/det_r[na, na, :]
#     mm = res.sum(axis=2)/res.shape[2]
# 
#     sc_chi = numpy.stack([mm[0,:], mm[4,:], mm[8,:], mm[1,:], mm[2,:], mm[5,:]], axis=0)
#     return sc_chi
# 
# 
# def calc_scm_chi(r_ccs):
#     """Calculate symmetry constraint matrix for susceptibility
#     It is supposed that determinant of r matrix is 1
#     """
#     
#     r_11, r_12, r_13 = r_ccs[0], r_ccs[1], r_ccs[2]
#     r_21, r_22, r_23 = r_ccs[3], r_ccs[4], r_ccs[5]
#     r_31, r_32, r_33 = r_ccs[6], r_ccs[7], r_ccs[8]
#     det_r =  r_11*r_22*r_33 - r_11*r_23*r_32  - r_12*r_21*r_33 + r_12*r_23*r_31  + r_13*r_21*r_32 - r_13*r_22*r_31
#     
#     m_mm = numpy.array([
#         [r_11*(r_22*r_33 - r_23*r_32), r_11*(-r_21*r_33 + r_23*r_31), r_11*(r_21*r_32 - r_22*r_31), r_12*(r_22*r_33 - r_23*r_32), r_12*(-r_21*r_33 + r_23*r_31), r_12*(r_21*r_32 - r_22*r_31), r_13*(r_22*r_33 - r_23*r_32), r_13*(-r_21*r_33 + r_23*r_31), r_13*(r_21*r_32 - r_22*r_31)],
#         [r_11*(-r_12*r_33 + r_13*r_32), r_11*(r_11*r_33 - r_13*r_31), r_11*(-r_11*r_32 + r_12*r_31), r_12*(-r_12*r_33 + r_13*r_32), r_12*(r_11*r_33 - r_13*r_31), r_12*(-r_11*r_32 + r_12*r_31), r_13*(-r_12*r_33 + r_13*r_32), r_13*(r_11*r_33 - r_13*r_31), r_13*(-r_11*r_32 + r_12*r_31)],
#         [r_11*(r_12*r_23 - r_13*r_22), r_11*(-r_11*r_23 + r_13*r_21), r_11*(r_11*r_22 - r_12*r_21), r_12*(r_12*r_23 - r_13*r_22), r_12*(-r_11*r_23 + r_13*r_21), r_12*(r_11*r_22 - r_12*r_21), r_13*(r_12*r_23 - r_13*r_22), r_13*(-r_11*r_23 + r_13*r_21), r_13*(r_11*r_22 - r_12*r_21)],
#         [r_21*(r_22*r_33 - r_23*r_32), r_21*(-r_21*r_33 + r_23*r_31), r_21*(r_21*r_32 - r_22*r_31), r_22*(r_22*r_33 - r_23*r_32), r_22*(-r_21*r_33 + r_23*r_31), r_22*(r_21*r_32 - r_22*r_31), r_23*(r_22*r_33 - r_23*r_32), r_23*(-r_21*r_33 + r_23*r_31), r_23*(r_21*r_32 - r_22*r_31)],
#         [r_21*(-r_12*r_33 + r_13*r_32), r_21*(r_11*r_33 - r_13*r_31), r_21*(-r_11*r_32 + r_12*r_31), r_22*(-r_12*r_33 + r_13*r_32), r_22*(r_11*r_33 - r_13*r_31), r_22*(-r_11*r_32 + r_12*r_31), r_23*(-r_12*r_33 + r_13*r_32), r_23*(r_11*r_33 - r_13*r_31), r_23*(-r_11*r_32 + r_12*r_31)],
#         [r_21*(r_12*r_23 - r_13*r_22), r_21*(-r_11*r_23 + r_13*r_21), r_21*(r_11*r_22 - r_12*r_21), r_22*(r_12*r_23 - r_13*r_22), r_22*(-r_11*r_23 + r_13*r_21), r_22*(r_11*r_22 - r_12*r_21), r_23*(r_12*r_23 - r_13*r_22), r_23*(-r_11*r_23 + r_13*r_21), r_23*(r_11*r_22 - r_12*r_21)],
#         [r_31*(r_22*r_33 - r_23*r_32), r_31*(-r_21*r_33 + r_23*r_31), r_31*(r_21*r_32 - r_22*r_31), r_32*(r_22*r_33 - r_23*r_32), r_32*(-r_21*r_33 + r_23*r_31), r_32*(r_21*r_32 - r_22*r_31), r_33*(r_22*r_33 - r_23*r_32), r_33*(-r_21*r_33 + r_23*r_31), r_33*(r_21*r_32 - r_22*r_31)],
#         [r_31*(-r_12*r_33 + r_13*r_32), r_31*(r_11*r_33 - r_13*r_31), r_31*(-r_11*r_32 + r_12*r_31), r_32*(-r_12*r_33 + r_13*r_32), r_32*(r_11*r_33 - r_13*r_31), r_32*(-r_11*r_32 + r_12*r_31), r_33*(-r_12*r_33 + r_13*r_32), r_33*(r_11*r_33 - r_13*r_31), r_33*(-r_11*r_32 + r_12*r_31)],
#         [r_31*(r_12*r_23 - r_13*r_22), r_31*(-r_11*r_23 + r_13*r_21), r_31*(r_11*r_22 - r_12*r_21), r_32*(r_12*r_23 - r_13*r_22), r_32*(-r_11*r_23 + r_13*r_21), r_32*(r_11*r_22 - r_12*r_21), r_33*(r_12*r_23 - r_13*r_22), r_33*(-r_11*r_23 + r_13*r_21), r_33*(r_11*r_22 - r_12*r_21)]],
#         dtype=r_ccs.dtype)
# 
#     flag = numpy.all(numpy.isclose(numpy.abs(det_r),1))
#     if flag:
#         res = m_mm*det_r[na, na, :]
#     else:
#         res = m_mm/det_r[na, na, :]
#     scm_chi = res.sum(axis=2)/res.shape[2]
#     return scm_chi
# 

def calc_m_r_inv_m(unit_cell_parameters, symm_elems, flag_unit_cell_parameters: bool=False):
    """Calculate M Rs M^-1
    """
    m_r, dder_hh = calc_m_m_by_unit_cell_parameters(unit_cell_parameters, flag_unit_cell_parameters=flag_unit_cell_parameters)
    res, dder_res = calc_m1_m2_inv_m1(m_r, symm_elems[4:, :], flag_m1 = flag_unit_cell_parameters, flag_m2 = False)
    dder = {}
    return res, dder
