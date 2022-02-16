# -*- coding: utf-8 -*-
"""
Expressions to calculate Debye-Waller factor.

For details see documentation.
"""

import numpy

from .unit_cell import calc_reciprocal_by_unit_cell_parameters

na = numpy.newaxis

def calc_beta_by_u_ij(
        u_ij, reciprocal_unit_cell_parameters,
        flag_u_ij: bool = False):
    r"""
    Calculate beta.

    u_11, u_22, u_33, u_12, u_13, u_23 = u_i

    cell is Cell class

    Output:
        $\beta_{11}$, $\beta_{22}$, $\beta_{33}$
        $\beta_{12}$, $\beta_{13}$, $\beta_{23}$
    """
    u_11, u_22, u_33 = u_ij[0], u_ij[1], u_ij[2]
    u_12, u_13, u_23 = u_ij[3], u_ij[4], u_ij[5]

    ia = reciprocal_unit_cell_parameters[0]
    ib = reciprocal_unit_cell_parameters[1]
    ic = reciprocal_unit_cell_parameters[2]
    beta_11 = 2.*numpy.pi**2*u_11*ia**2
    beta_22 = 2.*numpy.pi**2*u_22*ib**2
    beta_33 = 2.*numpy.pi**2*u_33*ic**2
    beta_12 = 2.*numpy.pi**2*u_12*ia*ib
    beta_13 = 2.*numpy.pi**2*u_13*ia*ic
    beta_23 = 2.*numpy.pi**2*u_23*ib*ic
    beta = numpy.stack([beta_11, beta_22, beta_33, beta_12, beta_13, beta_23], axis=0)
    dder = {}
    if flag_u_ij:
        dder["u_ij"] = numpy.stack([
            2.*numpy.pi**2*numpy.ones_like(u_11)*ia**2,
            2.*numpy.pi**2*numpy.ones_like(u_22)*ib**2,
            2.*numpy.pi**2*numpy.ones_like(u_33)*ic**2,
            2.*numpy.pi**2*numpy.ones_like(u_12)*ia*ib,
            2.*numpy.pi**2*numpy.ones_like(u_13)*ia*ic,
            2.*numpy.pi**2*numpy.ones_like(u_23)*ib*ic], axis=0)
    return beta, dder



def calc_u_ij_by_beta(
        beta, reciprocal_unit_cell_parameters,
        flag_beta: bool = False):
    r"""
    Calculate beta.

    beta_11, beta_22, beta_33, beta_12, beta_13, beta_23 = beta_i

    cell is Cell class

    Output:
        $\u_{11}$, $\u_{22}$, $\u_{33}$
        $\u_{12}$, $\u_{13}$, $\u_{23}$
    """
    beta_11, beta_22, beta_33 = beta[0], beta[1], beta[2]
    beta_12, beta_13, beta_23 = beta[3], beta[4], beta[5]
    pi_sq = numpy.pi**2
    ia = reciprocal_unit_cell_parameters[0]
    ib = reciprocal_unit_cell_parameters[1]
    ic = reciprocal_unit_cell_parameters[2]
    u_11 = beta_11/(2. * pi_sq * ia**2)
    u_22 = beta_22/(2. * pi_sq * ib**2)
    u_33 = beta_33/(2. * pi_sq * ic**2)
    u_12 = beta_12/(2. * pi_sq * ia*ib)
    u_13 = beta_13/(2. * pi_sq * ia*ic)
    u_23 = beta_23/(2. * pi_sq * ib*ic)
    u_ij = numpy.stack([u_11, u_22, u_33, u_12, u_13, u_23], axis=0)
    dder = {}
    if flag_beta:
        dder["beta"] = numpy.stack([
            numpy.ones_like(beta_11)/(2. * pi_sq * ia**2),
            numpy.ones_like(beta_22)/(2. * pi_sq * ib**2),
            numpy.ones_like(beta_33)/(2. * pi_sq * ic**2),
            numpy.ones_like(beta_12)/(2. * pi_sq * ia*ib),
            numpy.ones_like(beta_13)/(2. * pi_sq * ia*ic),
            numpy.ones_like(beta_23)/(2. * pi_sq * ib*ic)], axis=0)
    return u_ij, dder


def calc_b_iso_beta(
        unit_cell_parameters, atom_adp_type, atom_iso_param, atom_aniso_param):
    """
    Calculate b_iso and beta_ij based on atom_site and atom_sites.
    For each atom defined in atom_site.
    """
    coeff = float(8.*numpy.square(numpy.pi))

    reciprocal_unit_cell_parameters, dder = calc_reciprocal_by_unit_cell_parameters(
        unit_cell_parameters, flag_unit_cell_parameters=False)

    np_atom_adp_type, np_atom_iso_param, np_atom_aniso_param = \
        numpy.atleast_1d(atom_adp_type, atom_iso_param, atom_aniso_param)
    
    f_uiso = np_atom_adp_type == "Uiso"
    f_biso = np_atom_adp_type == "Biso"
    f_uovl = np_atom_adp_type == "Uovl"
    f_umpe = np_atom_adp_type == "Umpe"
    f_uani = np_atom_adp_type == "Uani"
    f_bovl = np_atom_adp_type == "Bovl"
    f_bani = np_atom_adp_type == "Bani"
    b_iso = numpy.zeros_like(np_atom_iso_param)
    beta = numpy.zeros_like(np_atom_aniso_param)

    b_iso[f_uiso] = coeff*np_atom_iso_param[f_uiso]
    b_iso[f_biso] = np_atom_iso_param[f_biso]
    b_iso[f_uovl] = coeff*np_atom_iso_param[f_uovl] #FIXME: redo it
    b_iso[f_umpe] = coeff*np_atom_iso_param[f_umpe] #FIXME: redo it
    b_iso[f_bovl] = np_atom_iso_param[f_bovl] #FIXME: redo it

    beta[:, f_uani], dder_beta = calc_beta_by_u_ij(
        np_atom_aniso_param[:, f_uani], reciprocal_unit_cell_parameters, flag_u_ij=False)

    beta[:, f_bani] = np_atom_aniso_param[:, f_bani]
    return b_iso, beta


def calc_param_iso_aniso_by_b_iso_beta(unit_cell_parameters, atom_adp_type, b_iso, beta):
    coeff = 1./float(8.*numpy.square(numpy.pi))

    reciprocal_unit_cell_parameters, dder = calc_reciprocal_by_unit_cell_parameters(
        unit_cell_parameters, flag_unit_cell_parameters=False)

    np_atom_adp_type, np_b_iso, np_beta = \
        numpy.atleast_1d(atom_adp_type, b_iso, beta)

    f_uiso = np_atom_adp_type == "Uiso"
    f_biso = np_atom_adp_type == "Biso"
    f_uovl = np_atom_adp_type == "Uovl"
    f_umpe = np_atom_adp_type == "Umpe"
    f_uani = np_atom_adp_type == "Uani"
    f_bovl = np_atom_adp_type == "Bovl"
    f_bani = np_atom_adp_type == "Bani"

    atom_iso_param = numpy.zeros_like(np_b_iso)
    atom_aniso_param = numpy.zeros_like(np_beta)

    atom_iso_param[f_uiso] = coeff*np_b_iso[f_uiso]
    atom_iso_param[f_biso] = b_iso[f_biso]
    atom_iso_param[f_uovl] = coeff*np_b_iso[f_uovl] #FIXME: redo it
    atom_iso_param[f_umpe] = coeff*np_b_iso[f_umpe] #FIXME: redo it
    atom_iso_param[f_bovl] = np_b_iso[f_bovl] #FIXME: redo it

    atom_aniso_param[:, f_bani] = np_beta[:, f_bani]
    atom_aniso_param[:, f_uani] = calc_u_ij_by_beta(
        np_beta[:, f_uani], reciprocal_unit_cell_parameters, flag_beta=False)[0]

    return atom_iso_param, atom_aniso_param


def calc_power_dwf_iso(sthovl, b_iso, flag_sthovl: bool = False,
                       flag_b_iso: bool = False):
    """Isotropic harmonic Debye-Waller factor.
    """
    sthovl_sq = numpy.square(sthovl)
    power_dwf_iso = b_iso*sthovl_sq
    dder = {}
    if flag_b_iso:
        dder["b_iso"] = sthovl_sq * numpy.ones_like(b_iso)
    if flag_sthovl:
        dder["sthovl"] = 2*b_iso*sthovl*numpy.ones_like(sthovl)
    return power_dwf_iso, dder


def calc_power_dwf_aniso(
        index_hkl, beta, symm_elems_r, flag_beta: bool = False):
    """Anisotropic harmonic Debye-Waller factor.
    """
    b_11, b_22, b_33 = beta[0], beta[1], beta[2]
    b_12, b_13, b_23 = beta[3], beta[4], beta[5]

    h = index_hkl[0]
    k = index_hkl[1]
    l = index_hkl[2]

    r_11, r_12, r_13 = symm_elems_r[0], symm_elems_r[1], symm_elems_r[2]
    r_21, r_22, r_23 = symm_elems_r[3], symm_elems_r[4], symm_elems_r[5]
    r_31, r_32, r_33 = symm_elems_r[6], symm_elems_r[7], symm_elems_r[8]

    h_s = h*r_11 + k*r_21 + l*r_31
    k_s = h*r_12 + k*r_22 + l*r_32
    l_s = h*r_13 + k*r_23 + l*r_33

    power_dwf_aniso = (b_11*numpy.square(h_s) + b_22*numpy.square(k_s) +
                       b_33*numpy.square(l_s) + 2.*b_12*h_s*k_s +
                       2.*b_13*h_s*l_s + 2.*b_23*k_s*l_s)
    dder = {}
    if flag_beta:
        ones_b = numpy.ones_like(b_11)
        der_b_11 = ones_b*numpy.square(h_s)
        der_b_22 = ones_b*numpy.square(k_s)
        der_b_33 = ones_b*numpy.square(l_s)
        der_b_12 = ones_b*2.*h_s*k_s
        der_b_13 = ones_b*2.*h_s*l_s
        der_b_23 = ones_b*2.*k_s*l_s
        dder["beta"] = numpy.stack([der_b_11, der_b_22, der_b_33,
                                    der_b_12, der_b_13, der_b_23], axis=0)
    return power_dwf_aniso, dder


def calc_dwf(index_hkl, sthovl, b_iso, beta, symm_elems_r,
             flag_sthovl: bool = False, flag_b_iso: bool = False,
             flag_beta: bool = False):
    """Calculate Debye-Waller factor."""

    # ["hkl", "atoms"]
    power_iso_2d, dder_power_iso_2d = calc_power_dwf_iso(
        sthovl, b_iso, flag_sthovl=flag_sthovl, flag_b_iso=flag_b_iso)
    #    sthovl[:, na], b_iso[na, :], flag_sthovl=flag_sthovl, flag_b_iso=flag_b_iso)

    # ["hkl", "atoms", "symmetry"]
    power_aniso_3d, dder_power_aniso_3d = calc_power_dwf_aniso(
        index_hkl, beta,
        symm_elems_r, flag_beta=flag_beta)
    #    index_hkl[:, :, na, na], beta[:, na, :, na],
    #    symm_elems_r[:, na, na, :], flag_beta=flag_beta)

    power_3d = power_iso_2d + power_aniso_3d
    #power_3d = power_iso_2d[:, :, na] + power_aniso_3d
    dwf = numpy.exp(-power_3d)
    dder = {}
    if flag_sthovl:
        dder["sthovl"] = None
    if flag_b_iso:
        dder["b_iso"] = None
    if flag_beta:
        dder["beta"] = None
    return dwf, dder

