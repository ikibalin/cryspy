from typing import Callable
import numpy

na = numpy.newaxis

def calc_m_sq_sin_sq_para(tensor_sigma, flag_tensor_sigma: bool = False):
    """Calculate the term P1 for paramagnetic sublattice.
    For details see documentation "Integrated intensity from powder diffraction".
    """
    sigma_11 = tensor_sigma[0]
    sigma_12 = tensor_sigma[1]
    sigma_22 = tensor_sigma[4]
    
    p_1 = 0.5*(numpy.square(numpy.abs(sigma_11)) + numpy.square(numpy.abs(sigma_22))) + \
        numpy.square(numpy.abs(sigma_12))

    dder = {}
    if flag_tensor_sigma:
        ones = numpy.ones_like(sigma_11.real)
        dder_sigma_11_real = sigma_11.real * ones
        dder_sigma_11_imag = sigma_11.imag * ones
        dder_sigma_22_real = sigma_22.real * ones
        dder_sigma_22_imag = sigma_22.imag * ones
        dder_sigma_12_real = 2*sigma_12.real * ones
        dder_sigma_12_imag = 2*sigma_12.imag * ones
        zeros = numpy.zeros_like(dder_sigma_11_real)
        dder["tensor_sigma_real"] = numpy.stack([
            dder_sigma_11_real, dder_sigma_12_real, zeros,
            zeros, dder_sigma_22_real, zeros, zeros, zeros, zeros], axis=0)
        dder["tensor_sigma_imag"] = numpy.stack([
            dder_sigma_11_imag, dder_sigma_12_imag, zeros,
            zeros, dder_sigma_22_imag, zeros, zeros, zeros, zeros], axis=0)
    return p_1, dder


def calc_m_sq_cos_sq_para(tensor_sigma, flag_tensor_sigma: bool = False):
    """Calculate the term P2 for paramagnetic sublattice.
    For details see documentation "Integrated intensity from powder diffraction".
    """
    sigma_13 = tensor_sigma[2]
    sigma_23 = tensor_sigma[5]
    
    p_2 = numpy.square(numpy.abs(sigma_13)) + numpy.square(numpy.abs(sigma_23))

    dder = {}
    if flag_tensor_sigma:
        ones = numpy.ones_like(sigma_13.real)
        dder_sigma_13_real = 2*sigma_13.real * ones
        dder_sigma_13_imag = 2*sigma_13.imag * ones
        dder_sigma_23_real = 2*sigma_23.real * ones
        dder_sigma_23_imag = 2*sigma_23.imag * ones
        zeros = numpy.zeros_like(dder_sigma_13_real)
        dder["tensor_sigma_real"] = numpy.stack([
            zeros, zeros, dder_sigma_13_real,
            zeros, zeros, dder_sigma_23_real, zeros, zeros, zeros], axis=0)
        dder["tensor_sigma_imag"] = numpy.stack([
            zeros, zeros, dder_sigma_13_imag,
            zeros, zeros, dder_sigma_23_imag, zeros, zeros, zeros], axis=0)
    return p_2, dder


def calc_cross_term_para(f_nucl, tensor_sigma, flag_f_nucl: bool = False, flag_tensor_sigma: bool = False):
    """Calculate the term P3 for paramagnetic sublattice.
    For details see documentation "Integrated intensity from powder diffraction".
    """
    sigma_11 = tensor_sigma[0]
    sigma_22 = tensor_sigma[4]
    
    p_3 = f_nucl.real * (sigma_11.real + sigma_22.real) + f_nucl.imag * (sigma_11.imag + sigma_22.imag)

    dder = {}
    if flag_f_nucl:
        dder["f_nucl_real"] = (sigma_11.real + sigma_22.real)*numpy.ones_like(f_nucl.real)
        dder["f_nucl_imag"] = (sigma_11.imag + sigma_22.imag)*numpy.ones_like(f_nucl.imag)

    if flag_tensor_sigma:
        ones = numpy.ones_like(sigma_11.real)
        dder_sigma_11_real = f_nucl.real * ones
        dder_sigma_11_imag = f_nucl.imag * ones
        dder_sigma_22_real = f_nucl.real * ones
        dder_sigma_22_imag = f_nucl.imag * ones
        zeros = numpy.zeros_like(dder_sigma_11_real)
        dder["tensor_sigma_real"] = numpy.stack([
            dder_sigma_11_real, zeros, zeros, 
            zeros, dder_sigma_22_real, zeros, zeros, zeros, zeros], axis=0)
        dder["tensor_sigma_imag"] = numpy.stack([
            dder_sigma_11_imag, zeros, zeros, 
            zeros, dder_sigma_22_imag, zeros, zeros, zeros, zeros], axis=0)
    return p_3, dder


def calc_chiral_term_cos_sin_sq_para(tensor_sigma, flag_tensor_sigma: bool = False):
    """Calculate the term P4 for paramagnetic sublattice.
    For details see documentation "Integrated intensity from powder diffraction".
    """
    sigma_11 = tensor_sigma[0]
    sigma_12 = tensor_sigma[1]
    sigma_21 = tensor_sigma[3]
    sigma_22 = tensor_sigma[4]
    
    p_4 = sigma_21.real*sigma_11.imag - sigma_21.imag*sigma_11.real + \
          sigma_22.real*sigma_12.imag - sigma_22.imag*sigma_12.real

    dder = {}
    if flag_tensor_sigma:
        ones = numpy.ones_like(sigma_11.real)
        dder_sigma_11_real = - sigma_21.imag * ones
        dder_sigma_11_imag = sigma_21.real * ones
        dder_sigma_12_real = - sigma_22.imag * ones
        dder_sigma_12_imag = sigma_22.real * ones
        dder_sigma_21_real = sigma_11.imag * ones
        dder_sigma_21_imag = - sigma_11.real * ones
        dder_sigma_22_real = sigma_12.imag * ones
        dder_sigma_22_imag = - sigma_12.real * ones
        zeros = numpy.zeros_like(dder_sigma_11_real)
        dder["tensor_sigma_real"] = numpy.stack([
            dder_sigma_11_real, dder_sigma_12_real, zeros, 
            dder_sigma_21_real, dder_sigma_22_real, zeros, zeros, zeros, zeros], axis=0)
        dder["tensor_sigma_imag"] = numpy.stack([
            dder_sigma_11_imag, dder_sigma_12_imag, zeros, 
            dder_sigma_21_imag, dder_sigma_22_imag, zeros, zeros, zeros, zeros], axis=0)
    return p_4, dder


def calc_chiral_term_cos_cube_para(tensor_sigma, flag_tensor_sigma: bool = False):
    """Calculate the term P5 for paramagnetic sublattice.
    For details see documentation "Integrated intensity from powder diffraction".
    """
    sigma_13 = tensor_sigma[2]
    sigma_23 = tensor_sigma[5]
    
    p_5 = 2*(sigma_23.real*sigma_13.imag - sigma_23.imag*sigma_13.real)

    dder = {}
    if flag_tensor_sigma:
        ones = numpy.ones_like(sigma_13.real)
        dder_sigma_13_real = -2 * sigma_23.imag * ones
        dder_sigma_13_imag = 2 * sigma_23.real * ones
        dder_sigma_23_real = 2 * sigma_13.imag * ones
        dder_sigma_23_imag = -2 * sigma_13.real * ones
        zeros = numpy.zeros_like(dder_sigma_13_real)
        dder["tensor_sigma_real"] = numpy.stack([
            zeros, zeros, dder_sigma_13_real, 
            zeros, zeros, dder_sigma_23_real, zeros, zeros, zeros], axis=0)
        dder["tensor_sigma_imag"] = numpy.stack([
            zeros, zeros, dder_sigma_13_imag, 
            zeros, zeros, dder_sigma_23_imag, zeros, zeros, zeros], axis=0)
    return p_5, dder


def calc_cross_term_ordered(f_nucl, f_m_perp, flag_f_nucl: bool = False, flag_f_m_perp: bool = False):
    """Calculate the term O1 for the magnetically ordered sublattice.
    For details see documentation "Integrated intensity from powder diffraction".
    """
    f_m_perp_z = f_m_perp[2]
    
    o_1 = 2*(f_nucl.real*f_m_perp_z.real + f_nucl.imag*f_m_perp_z.imag)

    dder = {}
    if flag_f_nucl:
        dder["f_nucl_real"] = 2 * f_m_perp_z.real * numpy.ones_like(f_nucl.real)
        dder["f_nucl_imag"] = 2 * f_m_perp_z.imag * numpy.ones_like(f_nucl.imag)
    if flag_f_m_perp:
        dder_f_m_perp_z_real = 2 * f_nucl.real * numpy.ones_like(f_m_perp_z.real)
        dder_f_m_perp_z_imag = 2 * f_nucl.imag * numpy.ones_like(f_m_perp_z.imag)
        zeros = numpy.zeros_like(dder_f_m_perp_z_real)
        dder["f_m_perp_real"] = numpy.stack([zeros, zeros, dder_f_m_perp_z_real], axis=0)
        dder["f_m_perp_imag"] = numpy.stack([zeros, zeros, dder_f_m_perp_z_imag], axis=0)
    return o_1, dder


def calc_chiral_term_ordered(f_m_perp, flag_f_m_perp: bool = False):
    """Calculate the term O2 for the magnetically ordered sublattice.
    For details see documentation "Integrated intensity from powder diffraction".
    """
    f_m_perp_x = f_m_perp[0]
    f_m_perp_y = f_m_perp[1]
    
    o_2 = 2*(f_m_perp_y.real*f_m_perp_x.imag - f_m_perp_y.imag*f_m_perp_x.real)

    dder = {}
    if flag_f_m_perp:
        dder_f_m_perp_x_real = -2 * f_m_perp_y.imag * numpy.ones_like(f_m_perp_x.real)
        dder_f_m_perp_x_imag = 2 * f_m_perp_y.real * numpy.ones_like(f_m_perp_x.imag)
        dder_f_m_perp_y_real = 2 * f_m_perp_x.imag * numpy.ones_like(f_m_perp_y.real)
        dder_f_m_perp_y_imag = -2 * f_m_perp_x.real * numpy.ones_like(f_m_perp_y.imag)
        zeros = numpy.zeros_like(dder_f_m_perp_x_real)
        dder["f_m_perp_real"] = numpy.stack([dder_f_m_perp_x_real, dder_f_m_perp_y_real, zeros], axis=0)
        dder["f_m_perp_imag"] = numpy.stack([dder_f_m_perp_x_imag, dder_f_m_perp_y_imag, zeros], axis=0)
    return o_2, dder


def calc_m_sq_mix(tensor_sigma, f_m_perp, flag_tensor_sigma: bool = False, flag_f_m_perp: bool = False):
    """Calculate the term M1 for the case of coexistiong paramatic and magnetically ordered sublattices.
    For details see documentation "Integrated intensity from powder diffraction".

    tensor_sigma describe paramagnetic sublattice
    f_m_perp describe ordered sublattice
    """
    sigma_13 = tensor_sigma[2]
    sigma_23 = tensor_sigma[5]
    f_m_perp_x = f_m_perp[0]
    f_m_perp_y = f_m_perp[1]
    
    m_1 = 2*(sigma_13.real*f_m_perp_x.real + sigma_13.imag*f_m_perp_x.imag + 
            sigma_23.real*f_m_perp_y.real + sigma_23.imag*f_m_perp_y.imag)

    dder = {}
    if flag_tensor_sigma:
        dder_sigma_13_real = 2*f_m_perp_x.real*numpy.ones_like(sigma_13.real)
        dder_sigma_13_imag = 2*f_m_perp_x.imag*numpy.ones_like(sigma_13.imag)
        dder_sigma_23_real = 2*f_m_perp_y.real*numpy.ones_like(sigma_23.real)
        dder_sigma_23_imag = 2*f_m_perp_y.imag*numpy.ones_like(sigma_23.imag)
        zeros = numpy.zeros_like(dder_sigma_13_real)
        dder["tensor_sigma_real"] = numpy.stack([
            zeros, zeros, dder_sigma_13_real,
            zeros, zeros, dder_sigma_23_real,
            zeros, zeros, zeros], axis=0)
        dder["tensor_sigma_imag"] = numpy.stack([
            zeros, zeros, dder_sigma_13_imag,
            zeros, zeros, dder_sigma_23_imag,
            zeros, zeros, zeros], axis=0)
    if flag_f_m_perp:
        dder_f_m_perp_x_real = 2 * sigma_13.real * numpy.ones_like(f_m_perp_x.real)
        dder_f_m_perp_x_imag = 2 * sigma_13.imag * numpy.ones_like(f_m_perp_x.imag)
        dder_f_m_perp_y_real = 2 * sigma_23.real * numpy.ones_like(f_m_perp_y.real)
        dder_f_m_perp_y_imag = 2 * sigma_13.imag * numpy.ones_like(f_m_perp_y.imag)
        zeros = numpy.zeros_like(dder_f_m_perp_x_real)
        dder["f_m_perp_real"] = numpy.stack([dder_f_m_perp_x_real, dder_f_m_perp_y_real, zeros], axis=0)
        dder["f_m_perp_imag"] = numpy.stack([dder_f_m_perp_x_imag, dder_f_m_perp_y_imag, zeros], axis=0)
    return m_1, dder


def calc_chiral_term_sin_sq_mix(
        tensor_sigma, f_m_perp, flag_tensor_sigma: bool = False, flag_f_m_perp: bool = False):
    """Calculate the term M2 for the case of coexistiong paramatic and magnetically ordered sublattices.
    For details see documentation "Integrated intensity from powder diffraction".

    tensor_sigma describe paramagnetic sublattice
    f_m_perp describe ordered sublattice
    """
    sigma_12 = tensor_sigma[1]
    sigma_21 = tensor_sigma[3]
    f_m_perp_z = f_m_perp[2]
    
    m_2 = sigma_12.real*f_m_perp_z.imag - sigma_12.imag*f_m_perp_z.real + \
        sigma_21.imag*f_m_perp_z.real - sigma_21.real*f_m_perp_z.imag

    dder = {}
    if flag_tensor_sigma:
        dder_sigma_12_real =  f_m_perp_z.imag*numpy.ones_like(sigma_12.real)
        dder_sigma_12_imag = -f_m_perp_z.real*numpy.ones_like(sigma_12.imag)
        dder_sigma_21_real = -f_m_perp_z.imag*numpy.ones_like(sigma_21.real)
        dder_sigma_21_imag =  f_m_perp_z.real*numpy.ones_like(sigma_21.imag)
        zeros = numpy.zeros_like(dder_sigma_12_real)
        dder["tensor_sigma_real"] = numpy.stack([
            zeros, dder_sigma_12_real, zeros,
            dder_sigma_21_real, zeros, zeros,
            zeros, zeros, zeros], axis=0)
        dder["tensor_sigma_imag"] = numpy.stack([
            zeros, dder_sigma_12_imag, zeros,
            dder_sigma_21_imag, zeros, zeros,
            zeros, zeros, zeros], axis=0)
    if flag_f_m_perp:
        dder_f_m_perp_z_real = (-sigma_12.imag + sigma_21.imag) * numpy.ones_like(f_m_perp_z.real)
        dder_f_m_perp_z_imag = ( sigma_12.real + sigma_21.real) * numpy.ones_like(f_m_perp_z.imag)
        zeros = numpy.zeros_like(dder_f_m_perp_z_real)
        dder["f_m_perp_real"] = numpy.stack([zeros, zeros, dder_f_m_perp_z_real], axis=0)
        dder["f_m_perp_imag"] = numpy.stack([zeros, zeros, dder_f_m_perp_z_imag], axis=0)
    return m_2, dder


def calc_chiral_term_cos_sq_mix(
        tensor_sigma, f_m_perp, flag_tensor_sigma: bool = False, flag_f_m_perp: bool = False):
    """Calculate the term M3 for the case of coexistiong paramatic and magnetically ordered sublattices.
    For details see documentation "Integrated intensity from powder diffraction".

    tensor_sigma describe paramagnetic sublattice
    f_m_perp describe ordered sublattice
    """
    sigma_13 = tensor_sigma[2]
    sigma_23 = tensor_sigma[5]
    f_m_perp_x = f_m_perp[0]
    f_m_perp_y = f_m_perp[1]
    
    m_3 = 2.*(
        sigma_23.real*f_m_perp_x.imag - sigma_23.imag*f_m_perp_x.real + \
        sigma_13.imag*f_m_perp_y.real - sigma_13.real*f_m_perp_y.imag)

    dder = {}
    if flag_tensor_sigma:
        dder_sigma_13_real = -f_m_perp_y.imag*numpy.ones_like(sigma_13.real)
        dder_sigma_13_imag =  f_m_perp_y.real*numpy.ones_like(sigma_13.imag)
        dder_sigma_23_real =  f_m_perp_x.imag*numpy.ones_like(sigma_23.real)
        dder_sigma_23_imag = -f_m_perp_x.real*numpy.ones_like(sigma_23.imag)
        zeros = numpy.zeros_like(dder_sigma_13_real)
        dder["tensor_sigma_real"] = numpy.stack([
            zeros, zeros, dder_sigma_13_real,
            zeros, zeros, dder_sigma_23_real, 
            zeros, zeros, zeros], axis=0)
        dder["tensor_sigma_imag"] = numpy.stack([
            zeros, zeros, dder_sigma_13_imag,
            zeros, zeros, dder_sigma_23_imag,
            zeros, zeros, zeros], axis=0)
    if flag_f_m_perp:
        dder_f_m_perp_x_real = -2 * sigma_23.imag * numpy.ones_like(f_m_perp_x.real)
        dder_f_m_perp_x_imag =  2 * sigma_23.real * numpy.ones_like(f_m_perp_x.imag)
        dder_f_m_perp_y_real =  2 * sigma_13.imag * numpy.ones_like(f_m_perp_y.real)
        dder_f_m_perp_y_imag = -2 * sigma_13.real * numpy.ones_like(f_m_perp_y.imag)
        zeros = numpy.zeros_like(dder_f_m_perp_x_real)
        dder["f_m_perp_real"] = numpy.stack([dder_f_m_perp_x_real, dder_f_m_perp_y_real, zeros], axis=0)
        dder["f_m_perp_imag"] = numpy.stack([dder_f_m_perp_x_imag, dder_f_m_perp_y_imag, zeros], axis=0)
    return m_3, dder


def calc_powder_iint_1d_para(f_nucl, tensor_sigma, polarization, flipper, magnetic_field,
        flag_f_nucl: bool = False, flag_tensor_sigma: bool = False,
        flag_polarization: bool = False, flag_flipper: bool = False):
    """Calculated powderly averaged integrated intensity for paramagnetic sublattice in equatorial plane (alpha=90 deg.)
    For details see documentation "Integrated intensity from powder diffraction".

    Output is two integrated intensities (plus and minus) with derivatives.
    """
    p_u = polarization
    p_d = polarization*(2.*flipper-1.)    

    f_nucl_sq = numpy.square(numpy.abs(f_nucl))
    p_1, dder_p_1 = calc_m_sq_sin_sq_para(tensor_sigma, flag_tensor_sigma=flag_tensor_sigma)
    p_3, dder_p_3 = calc_cross_term_para(f_nucl, tensor_sigma, flag_f_nucl=flag_f_nucl,
        flag_tensor_sigma=flag_tensor_sigma)
    magnetic_field_sq = numpy.square(magnetic_field)
    iint_plus = f_nucl_sq + magnetic_field_sq * p_1 + p_u * magnetic_field*p_3
    iint_minus = f_nucl_sq + magnetic_field_sq * p_1 - p_d * magnetic_field*p_3

    dder_plus = {}
    dder_minus = {}
    if flag_polarization:
        dder_plus["polarization"] = magnetic_field*p_3*numpy.ones_like(polarization)
        dder_minus["polarization"] = -1*magnetic_field*p_3*(2*flipper-1)*numpy.ones_like(polarization)
    if flag_flipper:
        dder_minus["flipper"] = -2*magnetic_field*p_3*polarization*numpy.ones_like(flipper)
    if flag_f_nucl:
        dder_plus["f_nucl_real"] = 2 * f_nucl.real + p_u * magnetic_field * dder_p_3["f_nucl_real"]
        dder_plus["f_nucl_imag"] = 2 * f_nucl.imag + p_u * magnetic_field * dder_p_3["f_nucl_imag"]
        dder_minus["f_nucl_real"] = 2 * f_nucl.real - p_d * magnetic_field * dder_p_3["f_nucl_real"]
        dder_minus["f_nucl_imag"] = 2 * f_nucl.imag - p_d * magnetic_field * dder_p_3["f_nucl_imag"]
    if flag_tensor_sigma:
        dder_plus["tensor_sigma_real"] = magnetic_field_sq * dder_p_1["tensor_sigma_real"] + p_u * magnetic_field*dder_p_3["tensor_sigma_real"]
        dder_plus["tensor_sigma_imag"] = magnetic_field_sq * dder_p_1["tensor_sigma_imag"] + p_u * magnetic_field*dder_p_3["tensor_sigma_imag"]
        dder_minus["tensor_sigma_real"] = magnetic_field_sq * dder_p_1["tensor_sigma_real"] - p_d * magnetic_field*dder_p_3["tensor_sigma_real"]
        dder_minus["tensor_sigma_imag"] = magnetic_field_sq * dder_p_1["tensor_sigma_imag"] - p_d * magnetic_field*dder_p_3["tensor_sigma_imag"]
    return iint_plus, iint_minus, dder_plus, dder_minus


def calc_powder_iint_1d_ordered(f_nucl, f_m_perp,
        flag_f_nucl: bool = False, flag_f_m_perp: bool = False):
    """Calculated powderly averaged integrated intensity for ordered sublattice in equatorial plane (alpha=90 deg.)
    For details see documentation "Integrated intensity from powder diffraction".

    In equatorial plane the scattering is not a function of polarized neutrons.
    Output is integrated intensities with derivatives.
    """
    f_nucl_sq = numpy.square(numpy.abs(f_nucl))
    f_m_perp_sq = numpy.square(numpy.abs(f_m_perp)).sum(axis=0)
    iint = f_nucl_sq + f_m_perp_sq

    dder = {}
    if flag_f_nucl:
        dder["f_nucl_real"] = 2 * f_nucl.real
        dder["f_nucl_imag"] = 2 * f_nucl.imag
    if flag_f_m_perp:
        dder["f_m_perp_real"] = 2*f_m_perp.real
        dder["f_m_perp_imag"] = 2*f_m_perp.imag
    return iint, dder


def calc_powder_iint_1d_mix(f_nucl, f_m_perp, tensor_sigma, polarization, flipper, magnetic_field,
        flag_f_nucl: bool = False, flag_f_m_perp_ordered: bool = False, flag_tensor_sigma: bool = False,
        flag_polarization: bool = False, flag_flipper: bool = False,
        flag_magnetic_field: bool = False):
    """Calculated powderly averaged integrated intensity in case of coexisting paramagnetic and magnetically ordered sublattice at
    scattering in equatorial plane (alpha=90 deg.).
    For details see documentation "Integrated intensity from powder diffraction".

    Output is two integrated intensities (plus and minus) with derivatives.
    """
    p_u = polarization
    p_d = polarization*(2.*flipper-1.)
    iint_plus_para, iint_minus_para, dder_plus_para, dder_minus_para = calc_powder_iint_1d_para(
        f_nucl, tensor_sigma, polarization, flipper, magnetic_field,
        flag_f_nucl=flag_f_nucl, flag_tensor_sigma=flag_tensor_sigma,
        flag_polarization=flag_polarization, flag_flipper=flag_flipper)

    f_m_perp_sq = numpy.square(numpy.abs(f_m_perp)).sum(axis=0)

    m_2, dder_m_2 = calc_chiral_term_sin_sq_mix(
        tensor_sigma, f_m_perp, flag_tensor_sigma=flag_tensor_sigma, flag_f_m_perp=flag_f_m_perp_ordered)

    iint_plus = iint_plus_para + f_m_perp_sq + p_u * magnetic_field * m_2
    iint_minus = iint_minus_para + f_m_perp_sq + p_d * magnetic_field * m_2

    dder_plus = {}
    dder_minus = {}
    if flag_f_nucl:
        dder_plus["f_nucl_real"] = dder_plus_para["f_nucl_real"] 
        dder_minus["f_nucl_real"] = dder_minus_para["f_nucl_real"] 
    if flag_f_m_perp_ordered:
        dder_plus["f_m_perp_ordered_real"] = 2*f_m_perp.real + p_u * magnetic_field * dder_m_2["f_m_perp_real"]
        dder_plus["f_m_perp_ordered_imag"] = 2*f_m_perp.imag + p_u * magnetic_field * dder_m_2["f_m_perp_imag"]
        dder_minus["f_m_perp_ordered_real"] = 2*f_m_perp.real + p_d * magnetic_field * dder_m_2["f_m_perp_real"]
        dder_minus["f_m_perp_ordered_imag"] = 2*f_m_perp.imag + p_d * magnetic_field * dder_m_2["f_m_perp_imag"]
    if flag_tensor_sigma:
        dder_plus["tensor_sigma_real"] = dder_plus_para["tensor_sigma_real"] + p_u * magnetic_field * dder_m_2["tensor_sigma_real"]
        dder_plus["tensor_sigma_imag"] = dder_plus_para["tensor_sigma_imag"] + p_u * magnetic_field * dder_m_2["tensor_sigma_imag"]
        dder_minus["tensor_sigma_real"] = dder_minus_para["tensor_sigma_real"] + p_d * magnetic_field * dder_m_2["tensor_sigma_real"]
        dder_minus["tensor_sigma_imag"] = dder_minus_para["tensor_sigma_imag"] + p_d * magnetic_field * dder_m_2["tensor_sigma_imag"]
    if flag_polarization:
        dder_plus["polarization"] = dder_plus_para["polarization"] + magnetic_field * m_2 * numpy.ones_like(polarization)
        dder_minus["polarization"] = dder_minus_para["polarization"] + magnetic_field * m_2 * (2.*flipper-1.)* numpy.ones_like(polarization)
    if flag_flipper:
        dder_minus["flipper"] = dder_minus_para["flippe"] + magnetic_field * m_2 * 2*polarization* numpy.ones_like(flipper)
    if flag_magnetic_field:
        dder_plus["magnetic_field"] = dder_plus_para["magnetic_field"] + p_u* m_2 * numpy.ones_like(magnetic_field)
        dder_minus["magnetic_field"] = dder_minus_para["magnetic_field"] + p_d* m_2 * numpy.ones_like(magnetic_field)
    return iint_plus, iint_minus, dder_plus, dder_minus


def calc_powder_iint_2d_para(f_nucl, tensor_sigma, beam_polarization, flipper_efficiency, magnetic_field,
        alpha_det, dict_in_out,
        flag_f_nucl: bool = False, flag_tensor_sigma: bool = False,
        flag_polarization: bool = False, flag_flipper: bool = False):
    """Calculated powderly averaged integrated intensity for paramagnetic sublattice in equatorial plane (alpha=90 deg.)
    For details see documentation "Integrated intensity from powder diffraction".

    Output is two integrated intensities (plus and minus) with derivatives.
    """
    p_u = beam_polarization
    p_d = beam_polarization*(2.*flipper_efficiency-1.)    

    s_a_sq = numpy.square(numpy.sin(alpha_det))[:, :, na]
    c_a = numpy.cos(alpha_det)[:, :, na]

    f_nucl_sq = numpy.square(numpy.abs(f_nucl))
    p_1, dder_p_1 = calc_m_sq_sin_sq_para(tensor_sigma, flag_tensor_sigma=flag_tensor_sigma)
    p_2, dder_p_2 = calc_m_sq_cos_sq_para(tensor_sigma, flag_tensor_sigma=flag_tensor_sigma)
    p_3, dder_p_3 = calc_cross_term_para(f_nucl, tensor_sigma, flag_f_nucl=flag_f_nucl,
        flag_tensor_sigma=flag_tensor_sigma)
    p_4, dder_p_4 = calc_chiral_term_cos_sin_sq_para(tensor_sigma, flag_tensor_sigma=flag_tensor_sigma)
    p_5, dder_p_5 = calc_chiral_term_cos_cube_para(tensor_sigma, flag_tensor_sigma=flag_tensor_sigma)
    magnetic_field_sq = numpy.square(magnetic_field)
    dict_in_out["p_1"] = p_1
    dict_in_out["p_2"] = p_2
    dict_in_out["p_3"] = p_3
    dict_in_out["p_4"] = p_4
    dict_in_out["p_5"] = p_5

    iint_plus = (f_nucl_sq + magnetic_field_sq * p_2)[na, na, :] + \
        (magnetic_field_sq * p_1 - magnetic_field_sq * p_2 + p_u * magnetic_field*p_3)[na, na, :] * s_a_sq + \
         (p_u * magnetic_field_sq)[na, na, :] * ((p_4-p_5)*s_a_sq + p_5[na, na, :])*c_a
    iint_minus = (f_nucl_sq + magnetic_field_sq * p_2)[na, na, :] + \
        (magnetic_field_sq * p_1 - magnetic_field_sq * p_2 - p_d * magnetic_field*p_3)[na, na, :] * s_a_sq - \
        (p_d * magnetic_field_sq)[na, na, :] * ((p_4-p_5)[na, na, :]*s_a_sq + p_5[na, na, :])*c_a

    dder_plus = {}
    dder_minus = {}
    if flag_polarization:
        dder_plus["beam_polarization"] = magnetic_field*p_3*s_a_sq + magnetic_field_sq * ((p_4-p_5)*s_a_sq + p_5)*c_a 
        dder_minus["beam_polarization"] = -(2*flipper_efficiency-1)*dder_plus["beam_polarization"]
    if flag_flipper:
        dder_minus["flipper_efficiency"] = -2*beam_polarization*(magnetic_field*p_3*s_a_sq + magnetic_field_sq * ((p_4-p_5)*s_a_sq + p_5)*c_a)
    if flag_f_nucl:
        dder_plus["f_nucl_real"] = 2 * f_nucl.real + p_u * magnetic_field * s_a_sq * dder_p_3["f_nucl_real"]
        dder_plus["f_nucl_imag"] = 2 * f_nucl.imag + p_u * magnetic_field * s_a_sq * dder_p_3["f_nucl_imag"]
        dder_minus["f_nucl_real"] = 2 * f_nucl.real - p_d * magnetic_field * s_a_sq * dder_p_3["f_nucl_real"]
        dder_minus["f_nucl_imag"] = 2 * f_nucl.imag - p_d * magnetic_field * s_a_sq * dder_p_3["f_nucl_imag"]
    if flag_tensor_sigma:
        dder_plus["tensor_sigma_real"] = None
        dder_plus["tensor_sigma_imag"] = None
        dder_minus["tensor_sigma_real"] = None
        dder_minus["tensor_sigma_imag"] = None
    return iint_plus, iint_minus, dder_plus, dder_minus


def calc_powder_iint_2d_ordered(f_nucl, f_m_perp, beam_polarization, flipper_efficiency,
        alpha_det, dict_in_out, 
        flag_f_nucl: bool = False, flag_f_m_perp: bool = False,
        flag_polarization: bool = False, flag_flipper: bool = False):
    """Calculated powderly averaged integrated intensity for paramagnetic sublattice in equatorial plane (alpha=90 deg.)
    For details see documentation "Integrated intensity from powder diffraction".

    Output is two integrated intensities (plus and minus) with derivatives.
    """
    p_u = beam_polarization
    p_d = beam_polarization*(2.*flipper_efficiency-1.)    

    s_a_sq = numpy.square(numpy.sin(alpha_det))[:, :, na]
    c_a = numpy.cos(alpha_det)[:, :, na]

    f_nucl_sq = numpy.square(numpy.abs(f_nucl))
    f_m_perp_sq = ((f_m_perp *numpy.conjugate(f_m_perp)).real).sum(axis=0)
    o_1, dder_o_1 = calc_cross_term_ordered(f_nucl, f_m_perp, flag_f_nucl=flag_f_nucl, flag_f_m_perp=flag_f_m_perp)
    o_2, dder_o_2 = calc_chiral_term_ordered(f_m_perp, flag_f_m_perp=flag_f_m_perp)
    dict_in_out["o_1"] = o_1
    dict_in_out["o_2"] = o_2

    iint_plus = (f_nucl_sq + f_m_perp_sq)[na, na, :] + p_u * (o_1+o_2)[na, na, :] * c_a
    iint_minus = (f_nucl_sq + f_m_perp_sq)[na, na, :] - p_d * (o_1+o_2)[na, na, :] * c_a

    dder_plus = {}
    dder_minus = {}
    if flag_polarization:
        dder_plus["beam_polarization"] = (o_1 + o_2) * c_a
        dder_minus["beam_polarization"] = -(2*flipper_efficiency - 1) * dder_plus["beam_polarization"]
    if flag_flipper:
        dder_minus["flipper_efficiency"] = -beam_polarization * (o_1 + o_2) * c_a
    if flag_f_nucl:
        dder_plus["f_nucl_real"] = 2 * f_nucl.real 
        dder_plus["f_nucl_imag"] = 2 * f_nucl.imag 
        dder_minus["f_nucl_real"] = 2 * f_nucl.real 
        dder_minus["f_nucl_imag"] = 2 * f_nucl.imag 
    if flag_f_m_perp:
        dder_plus["tensor_sigma_real"] = None
        dder_plus["tensor_sigma_imag"] = None
        dder_minus["tensor_sigma_real"] = None
        dder_minus["tensor_sigma_imag"] = None
    return iint_plus, iint_minus, dder_plus, dder_minus



def calc_powder_iint_2d_mix(f_nucl, tensor_sigma, f_m_perp_ordered, beam_polarization, flipper_efficiency, magnetic_field,
        alpha_det, dict_in_out_phase,
        flag_f_nucl: bool = False, flag_tensor_sigma: bool = False, flag_f_m_perp_ordered: bool = False,
        flag_polarization: bool = False, flag_flipper: bool = False):
    """Calculated powderly averaged integrated intensity for paramagnetic sublattice in equatorial plane (alpha=90 deg.)
    For details see documentation "Integrated intensity from powder diffraction".

    Output is two integrated intensities (plus and minus) with derivatives.
    """
    p_u = beam_polarization
    p_d = beam_polarization*(2.*flipper_efficiency-1.)    

    s_a_sq = numpy.square(numpy.sin(alpha_det))[:, :, na]
    c_a = numpy.cos(alpha_det)[:, :, na]

    f_nucl_sq = numpy.square(numpy.abs(f_nucl))
    p_1, dder_p_1 = calc_m_sq_sin_sq_para(tensor_sigma, flag_tensor_sigma=flag_tensor_sigma)
    p_2, dder_p_2 = calc_m_sq_cos_sq_para(tensor_sigma, flag_tensor_sigma=flag_tensor_sigma)
    p_3, dder_p_3 = calc_cross_term_para(f_nucl, tensor_sigma, flag_f_nucl=flag_f_nucl,
        flag_tensor_sigma=flag_tensor_sigma)
    p_4, dder_p_4 = calc_chiral_term_cos_sin_sq_para(tensor_sigma, flag_tensor_sigma=flag_tensor_sigma)
    p_5, dder_p_5 = calc_chiral_term_cos_cube_para(tensor_sigma, flag_tensor_sigma=flag_tensor_sigma)
    magnetic_field_sq = numpy.square(magnetic_field)
    dict_in_out_phase["p_1"] = p_1
    dict_in_out_phase["p_2"] = p_2
    dict_in_out_phase["p_3"] = p_3
    dict_in_out_phase["p_4"] = p_4
    dict_in_out_phase["p_5"] = p_5

    f_m_perp_sq = ((f_m_perp_ordered *numpy.conjugate(f_m_perp_ordered)).real).sum(axis=0)
    o_1, dder_o_1 = calc_cross_term_ordered(f_nucl, f_m_perp_ordered, flag_f_nucl=flag_f_nucl, flag_f_m_perp=flag_f_m_perp_ordered)
    o_2, dder_o_2 = calc_chiral_term_ordered(f_m_perp_ordered, flag_f_m_perp=flag_f_m_perp_ordered)
    dict_in_out_phase["o_1"] = o_1
    dict_in_out_phase["o_2"] = o_2

    m_1, dder_m_1 = calc_m_sq_mix(
        tensor_sigma, f_m_perp_ordered, flag_tensor_sigma=flag_tensor_sigma, flag_f_m_perp=flag_f_m_perp_ordered)
    m_2, dder_m_2 = calc_chiral_term_sin_sq_mix(
        tensor_sigma, f_m_perp_ordered, flag_tensor_sigma=flag_tensor_sigma, flag_f_m_perp=flag_f_m_perp_ordered)
    m_3, dder_m_3 = calc_chiral_term_cos_sq_mix(
        tensor_sigma, f_m_perp_ordered, flag_tensor_sigma=flag_tensor_sigma, flag_f_m_perp=flag_f_m_perp_ordered)

    dict_in_out_phase["m_1"] = m_1
    dict_in_out_phase["m_2"] = m_2
    dict_in_out_phase["m_3"] = m_3

    iint_plus = (f_nucl_sq + magnetic_field_sq * p_2 + f_m_perp_sq)[na, na, :] + p_u*magnetic_field*m_3+\
        (magnetic_field_sq * p_1 - magnetic_field_sq * p_2 + p_u * magnetic_field*p_3 + p_u*magnetic_field*m_2-p_u*magnetic_field*m_3)[na, na, :] * s_a_sq + \
         (p_u * magnetic_field_sq)[na, na, :] * ((p_4-p_5)*s_a_sq + p_5[na, na, :])*c_a + \
         + (p_u*o_1+p_u*o_2+magnetic_field*m_1)[na, na, :] * c_a 
    iint_minus = (f_nucl_sq + magnetic_field_sq * p_2 + f_m_perp_sq)[na, na, :] - p_d*magnetic_field*m_3+\
        (magnetic_field_sq * p_1 - magnetic_field_sq * p_2 - p_d * magnetic_field*p_3 - p_d*magnetic_field*m_2-p_d*magnetic_field*m_3)[na, na, :] * s_a_sq - \
        (p_d * magnetic_field_sq)[na, na, :] * ((p_4-p_5)[na, na, :]*s_a_sq + p_5[na, na, :])*c_a +\
        (-p_d*o_1-p_d*o_2+magnetic_field*m_1)[na, na, :] * c_a

    dder_plus = {}
    dder_minus = {}
    return iint_plus, iint_minus, dder_plus, dder_minus
