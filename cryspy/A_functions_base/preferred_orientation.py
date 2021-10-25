import numpy

from .unit_cell import calc_eq_ccs_by_unit_cell_parameters
from .powder_diffraction_const_wavelength import calc_gamma_nu_by_ttheta_phi

def calc_cos_ang(cell, h_1, k_1, l_1, h_2, k_2, l_2):
    """Calculate directed cosines."""
    q_1_x, q_1_y, q_1_z = cell.calc_k_loc(h_1, k_1, l_1)
    q_2_x, q_2_y, q_2_z = cell.calc_k_loc(h_2, k_2, l_2)
    q_1_sq = q_1_x*q_1_x + q_1_y*q_1_y + q_1_z*q_1_z
    q_2_sq = q_2_x*q_2_x + q_2_y*q_2_y + q_2_z*q_2_z
    q_12 = q_1_x*q_2_x + q_1_y*q_2_y + q_1_z*q_2_z
    res = q_12/(q_1_sq*q_2_sq)**0.5
    res[res > 1.] = 1.
    return res

na = numpy.newaxis

def calc_preferred_orientation_pd(
        index_hkl, texture_g1, texture_g2, texture_axis, unit_cell_parameters,
        flag_texture_g1: bool = False, flag_texture_g2: bool = False,
        flag_texture_axis: bool = False, flag_unit_cell_parameters: bool = False):
    """Preferred orintation by Modified March-Dollas model.
    """
    preferred_orientation = None 
    eq_axis_ccs, dder_eq_axis_ccs = calc_eq_ccs_by_unit_cell_parameters(
        texture_axis, unit_cell_parameters, flag_unit_cell_parameters=flag_unit_cell_parameters)

    eq_hkl_ccs, dder_eq_hkl_ccs = calc_eq_ccs_by_unit_cell_parameters(
        index_hkl, unit_cell_parameters, flag_unit_cell_parameters=flag_unit_cell_parameters)

    cos_alpha_ax = (eq_axis_ccs[:, 0][:, na] * eq_hkl_ccs).sum(axis=0)

    sin_alpha_ax_sq = numpy.abs(1. - numpy.square(cos_alpha_ax))
    
    cos_alpha_sq = sin_alpha_ax_sq
    hh = 1./texture_g1 + (texture_g1**2 - 1./texture_g1) * cos_alpha_sq
    
    preferred_orientation = texture_g2 + (1. - texture_g2) * numpy.power(hh, -1.5)

    dder_po = {}
    if flag_texture_g2:
        dder_po["texture_g2"] = 1. - numpy.power(hh, -1.5)
    if flag_texture_g1:
        dder_po["texture_g1"] = -1.5*(1. - texture_g2) * numpy.power(hh, -2.5) * \
            (- (1.-cos_alpha_sq)/numpy.square(texture_g1) + 2*cos_alpha_sq*texture_g1)

    return preferred_orientation, dder_po


def calc_preferred_orientation_pd2d(alpha_det,
        index_hkl, texture_g1, texture_g2, texture_axis, unit_cell_parameters,
        flag_texture_g1: bool = False, flag_texture_g2: bool = False,
        flag_texture_axis: bool = False, flag_unit_cell_parameters: bool = False):
    """Preferred orintation by Modified March-Dollas model.
    """
    preferred_orientation = None 
    eq_axis_ccs, dder_eq_axis_ccs = calc_eq_ccs_by_unit_cell_parameters(
        texture_axis, unit_cell_parameters, flag_unit_cell_parameters=flag_unit_cell_parameters)

    eq_hkl_ccs, dder_eq_hkl_ccs = calc_eq_ccs_by_unit_cell_parameters(
        index_hkl, unit_cell_parameters, flag_unit_cell_parameters=flag_unit_cell_parameters)


    cos_alpha_ax = (eq_axis_ccs[:, 0][:, na] * eq_hkl_ccs).sum(axis=0)
    cos_alpha_ax[cos_alpha_ax > 1.] = 1.

    c_help = 1.-cos_alpha_ax**2
    c_help[c_help < 0.] = 0.
    sin_alpha_ax = numpy.sqrt(c_help)

    cos_alpha_det = numpy.cos(alpha_det)
    sin_alpha_det = numpy.sin(alpha_det)

    cos_alpha = cos_alpha_ax[na, na,:] * cos_alpha_det[:,:,na] + sin_alpha_ax[na, na,:] * sin_alpha_det[:,:,na]
    cos_alpha_sq = numpy.square(cos_alpha)
    hh = 1./texture_g1 + (texture_g1**2 - 1./texture_g1) * cos_alpha_sq
    
    preferred_orientation = texture_g2 + (1. - texture_g2) * numpy.power(hh, -1.5)

    dder_po = {}
    if flag_texture_g2:
        dder_po["texture_g2"] = 1. - numpy.power(hh, -1.5)
    if flag_texture_g1:
        dder_po["texture_g1"] = -1.5*(1. - texture_g2) * numpy.power(hh, -2.5) * \
            (- (1.-cos_alpha_sq)/numpy.square(texture_g1) + 2*cos_alpha_sq*texture_g1)

    return preferred_orientation, dder_po


def calc_gamma_nu_for_textured_peaks(eq_axis, eq_ccs, ttheta_hkl, texture_g1):
    c_csc = numpy.sum(eq_ccs*eq_axis, axis=0)
    inv_c_h = 1./numpy.cos(0.5 * ttheta_hkl)

    if texture_g1 <= 1:
        s_phi = inv_c_h*c_csc
        s_phi[s_phi>1] = 1.
        s_phi[s_phi<-1] = -1.
        phi_max = numpy.arcsin(s_phi)
    else:
        s_csc  = numpy.sqrt(1.-numpy.square(c_csc))
        s_csc[s_csc>1.] = 1.
        s_phi = inv_c_h*s_csc
        s_phi[s_phi>1] = 1.
        s_phi[s_phi<-1] = -1.
        phi_max = numpy.arcsin(s_phi)
    gamma_hkl, nu_hkl = calc_gamma_nu_by_ttheta_phi(
        ttheta_hkl, phi_max, flag_ttheta=False, flag_phi=False)[:2]
    return gamma_hkl, nu_hkl 