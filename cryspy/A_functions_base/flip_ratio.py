from typing import Callable
import numpy


def calc_iint(
        beam_polarization: float, flipper_efficiency: float, f_nucl, f_m_perp, matrix_u, func_extinction: Callable = None,
        flag_beam_polarization: bool = False, flag_flipper_efficiency: bool = False,
        flag_f_nucl: bool = False, flag_f_m_perp: bool = False,
        dict_in_out: dict = None, flag_use_precalculated_data: bool = False):
    """Calculate integrated intensities.

    It is supposed that crystal is not rotate during the calculations 
    (orientation matrix is the same for all reflections)
    """
    if dict_in_out is None:
        flag_dict = False
        dict_in_out_keys = []
    else:
        flag_dict = True
        dict_in_out_keys = dict_in_out.keys()

    p_u = beam_polarization
    p_d = beam_polarization*(2*flipper_efficiency-1)
    flag_f_plus_sq = flag_f_nucl or flag_f_m_perp
    flag_f_minus_sq = flag_f_nucl or flag_f_m_perp
    flag_f_m_perp_xy_sq = flag_f_m_perp

    f_n = numpy.atleast_1d(f_nucl)
    f_m_perp = numpy.atleast_1d(f_m_perp)
    f_n_sq = numpy.square(numpy.abs(f_n))

    f_m_perp_x = f_m_perp[0]*matrix_u[0] + f_m_perp[1]*matrix_u[1] + f_m_perp[2]*matrix_u[2]
    f_m_perp_y = f_m_perp[0]*matrix_u[3] + f_m_perp[1]*matrix_u[4] + f_m_perp[2]*matrix_u[5]
    f_m_perp_z = f_m_perp[0]*matrix_u[6] + f_m_perp[1]*matrix_u[7] + f_m_perp[2]*matrix_u[8]
    
    f_m_perp_x_sq = numpy.square(numpy.abs(f_m_perp_x))
    f_m_perp_y_sq = numpy.square(numpy.abs(f_m_perp_y))
    f_m_perp_z_sq = numpy.square(numpy.abs(f_m_perp_z))

    f_n_f_m_perp_z = 2.*(f_n.real * f_m_perp_z.real + f_n.imag * f_m_perp_z.imag)

    f_plus_sq = f_n_sq + f_m_perp_z_sq + f_n_f_m_perp_z
    f_minus_sq = f_n_sq + f_m_perp_z_sq - f_n_f_m_perp_z
    f_m_perp_xy_sq = f_m_perp_x_sq + f_m_perp_y_sq

    if func_extinction is None:
        dder_y_plus, dder_y_minus, dder_y_m_perp_xy = {}, {}, {}
        y_plus = numpy.ones_like(f_plus_sq)
        y_minus = numpy.ones_like(f_minus_sq)
        y_m_perp_xy = numpy.ones_like(f_m_perp_xy_sq)
    else:
        y_plus, dder_y_plus = func_extinction(f_plus_sq, flag_f_sq=flag_f_plus_sq)
        y_minus, dder_y_minus = func_extinction(f_minus_sq, flag_f_sq=flag_f_minus_sq)
        y_m_perp_xy, dder_y_m_perp_xy = func_extinction(f_m_perp_xy_sq, flag_f_sq=flag_f_m_perp_xy_sq)

    chiral_term = 2.*(f_m_perp_x.imag * f_m_perp_y.real - f_m_perp_x.real * f_m_perp_y.imag)
    if flag_dict:
        dict_in_out["chiral_term"] = chiral_term

    iint_plus = 0.5*((1.+p_u)*y_plus*f_plus_sq +
                     (1.-p_u)*y_minus*f_minus_sq) + \
                    y_m_perp_xy * f_m_perp_xy_sq + \
                    p_u * chiral_term

    iint_minus = 0.5*((1.-p_d)*y_plus*f_plus_sq +
                     (1.+p_d)*y_minus*f_minus_sq) + \
                    y_m_perp_xy * f_m_perp_xy_sq - \
                    p_d * chiral_term

    if flag_dict:
        dict_in_out["iint_plus"] = iint_plus
        dict_in_out["iint_minus"] = iint_minus
        dict_in_out["y_plus"] = y_plus
        dict_in_out["f_plus_sq"] = f_plus_sq
        dict_in_out["y_minus"] = y_minus
        dict_in_out["f_minus_sq"] = f_minus_sq
        dict_in_out["y_m_perp_xy"] = y_m_perp_xy
        dict_in_out["f_m_perp_xy_sq"] = f_m_perp_xy_sq

    dder_plus = {}
    dder_minus = {}
    if flag_beam_polarization:
        dder_plus["beam_polarization"] = 0.5*(y_plus*f_plus_sq - y_minus*f_minus_sq + 2.*chiral_term) * \
            numpy.ones_like(beam_polarization)
        dder_minus["beam_polarization"] = 0.5*(-y_plus*f_plus_sq + y_minus*f_minus_sq - 2.*chiral_term) * \
            numpy.ones_like(beam_polarization)*(2.*flipper_efficiency-1.)

    if flag_flipper_efficiency:
        dder_minus["flipper_efficiency"] = beam_polarization*(-y_plus*f_plus_sq + y_minus*f_minus_sq - 2.*chiral_term) * \
            numpy.ones_like(flipper_efficiency)
    
    if flag_f_nucl:
        f_plus_sq_f_n_real = 2. * (f_n.real + f_m_perp_z.real) * numpy.ones_like(f_n.real)
        f_plus_sq_f_n_imag = 2. * (f_n.imag + f_m_perp_z.imag) * numpy.ones_like(f_n.imag)
        f_minus_sq_f_n_real = 2. * (f_n.real - f_m_perp_z.real)* numpy.ones_like(f_n.real)
        f_minus_sq_f_n_imag = 2. * (f_n.imag - f_m_perp_z.imag) * numpy.ones_like(f_n.imag)

        y_plus_f_n_real, y_minus_f_n_real = 0, 0 
        y_plus_f_n_imag, y_minus_f_n_imag = 0, 0 
        if "f_sq" in dder_y_plus.keys(): 
            y_plus_f_n_real = dder_y_plus["f_sq"]*f_plus_sq_f_n_real
            y_plus_f_n_imag = dder_y_plus["f_sq"]*f_plus_sq_f_n_imag
        if "f_sq" in dder_y_minus.keys(): 
            y_minus_f_n_real = dder_y_minus["f_sq"]*f_minus_sq_f_n_real
            y_minus_f_n_imag = dder_y_minus["f_sq"]*f_minus_sq_f_n_imag

        dder_plus["f_nucl_real"] =  0.5*(
            (1.+p_u)*(y_plus*f_plus_sq_f_n_real+y_plus_f_n_real*f_plus_sq) +
            (1.-p_u)*(y_minus*f_minus_sq_f_n_real+y_minus_f_n_real*f_minus_sq))
        
        dder_plus["f_nucl_imag"] =  0.5*(
            (1.+p_u)*(y_plus*f_plus_sq_f_n_imag+y_plus_f_n_imag*f_plus_sq) +
            (1.-p_u)*(y_minus*f_minus_sq_f_n_imag+y_minus_f_n_imag*f_minus_sq))

        dder_minus["f_nucl_real"] =  0.5*(
            (1.-p_d)*(y_plus*f_plus_sq_f_n_real+y_plus_f_n_real*f_plus_sq) +
            (1.+p_d)*(y_minus*f_minus_sq_f_n_real+y_minus_f_n_real*f_minus_sq))
        
        dder_minus["f_nucl_imag"] =  0.5*(
            (1.-p_d)*(y_plus*f_plus_sq_f_n_imag+y_plus_f_n_imag*f_plus_sq) +
            (1.+p_d)*(y_minus*f_minus_sq_f_n_imag+y_minus_f_n_imag*f_minus_sq))


    if flag_f_m_perp:
        f_plus_sq_f_m_perp_z_real = 2. * (f_n.real + f_m_perp_z.real) * numpy.ones_like(f_m_perp_z.real)
        f_plus_sq_f_m_perp_z_imag = 2. * (f_n.imag + f_m_perp_z.imag) * numpy.ones_like(f_m_perp_z.imag)
        f_minus_sq_f_m_perp_z_real = -2. * (f_n.real - f_m_perp_z.real) * numpy.ones_like(f_m_perp_z.real)
        f_minus_sq_f_m_perp_z_imag = -2. * (f_n.imag - f_m_perp_z.imag) * numpy.ones_like(f_m_perp_z.imag)
        f_m_perp_xy_sq_f_m_perp_x_real = 2 * f_m_perp_x.real * numpy.ones_like(f_m_perp_x.real)
        f_m_perp_xy_sq_f_m_perp_x_imag = 2 * f_m_perp_x.imag * numpy.ones_like(f_m_perp_x.imag)
        f_m_perp_xy_sq_f_m_perp_y_real = 2 * f_m_perp_y.real * numpy.ones_like(f_m_perp_y.real)
        f_m_perp_xy_sq_f_m_perp_y_imag = 2 * f_m_perp_y.imag * numpy.ones_like(f_m_perp_y.imag)
        chiral_term_f_m_perp_x_real = -2 * f_m_perp_y.imag * numpy.ones_like(f_m_perp_x.real)
        chiral_term_f_m_perp_x_imag = 2 * f_m_perp_y.real * numpy.ones_like(f_m_perp_x.imag)
        chiral_term_f_m_perp_y_real = 2 * f_m_perp_x.imag * numpy.ones_like(f_m_perp_y.real)
        chiral_term_f_m_perp_y_imag = -2 * f_m_perp_x.real * numpy.ones_like(f_m_perp_y.imag)

        y_plus_f_m_perp_z_real, y_plus_f_m_perp_z_imag = 0, 0 
        y_minus_f_m_perp_z_real, y_minus_f_m_perp_z_imag = 0, 0 
        y_m_perp_xy_f_m_perp_x_real, y_m_perp_xy_f_m_perp_x_imag = 0, 0 
        y_m_perp_xy_f_m_perp_y_real, y_m_perp_xy_f_m_perp_y_imag = 0, 0 
        if "f_sq" in dder_y_plus.keys(): 
            y_plus_f_m_perp_z_real = dder_y_plus["f_sq"]*f_plus_sq_f_m_perp_z_real
            y_plus_f_m_perp_z_imag = dder_y_plus["f_sq"]*f_plus_sq_f_m_perp_z_imag
        if "f_sq" in dder_y_minus.keys(): 
            y_minus_f_m_perp_z_real = dder_y_minus["f_sq"]*f_minus_sq_f_m_perp_z_real
            y_minus_f_m_perp_z_imag = dder_y_minus["f_sq"]*f_minus_sq_f_m_perp_z_imag
        if "f_sq" in dder_y_m_perp_xy.keys(): 
            y_m_perp_xy_f_m_perp_x_real = dder_y_m_perp_xy["f_sq"]*f_m_perp_xy_sq_f_m_perp_x_real
            y_m_perp_xy_f_m_perp_x_imag = dder_y_m_perp_xy["f_sq"]*f_m_perp_xy_sq_f_m_perp_x_imag
            y_m_perp_xy_f_m_perp_y_real = dder_y_m_perp_xy["f_sq"]*f_m_perp_xy_sq_f_m_perp_y_real
            y_m_perp_xy_f_m_perp_y_imag = dder_y_m_perp_xy["f_sq"]*f_m_perp_xy_sq_f_m_perp_y_imag

        dder_plus_f_m_perp_x_real = \
            y_m_perp_xy_f_m_perp_x_real * f_m_perp_xy_sq + \
            y_m_perp_xy * f_m_perp_xy_sq_f_m_perp_x_real +\
            p_u * chiral_term_f_m_perp_x_real
        dder_plus_f_m_perp_x_imag = \
            y_m_perp_xy_f_m_perp_x_imag * f_m_perp_xy_sq + \
            y_m_perp_xy * f_m_perp_xy_sq_f_m_perp_x_imag +\
            p_u * chiral_term_f_m_perp_x_imag
        dder_plus_f_m_perp_y_real = \
            y_m_perp_xy_f_m_perp_y_real * f_m_perp_xy_sq + \
            y_m_perp_xy * f_m_perp_xy_sq_f_m_perp_y_real +\
            p_u * chiral_term_f_m_perp_y_real
        dder_plus_f_m_perp_y_imag = \
            y_m_perp_xy_f_m_perp_y_imag * f_m_perp_xy_sq + \
            y_m_perp_xy * f_m_perp_xy_sq_f_m_perp_y_imag +\
            p_u * chiral_term_f_m_perp_y_imag
        dder_plus_f_m_perp_z_real =  0.5*(
            (1.+p_u)*(y_plus*f_plus_sq_f_m_perp_z_real+y_plus_f_m_perp_z_real*f_plus_sq) +
            (1.-p_u)*(y_minus*f_minus_sq_f_m_perp_z_real+y_minus_f_m_perp_z_real*f_minus_sq))
        dder_plus_f_m_perp_z_imag =  0.5*(
            (1.+p_u)*(y_plus*f_plus_sq_f_m_perp_z_imag+y_plus_f_m_perp_z_imag*f_plus_sq) +
            (1.-p_u)*(y_minus*f_minus_sq_f_m_perp_z_imag+y_minus_f_m_perp_z_imag*f_minus_sq))

        dder_plus_f_m_perp_1_real = \
            dder_plus_f_m_perp_x_real*matrix_u[0] + \
            dder_plus_f_m_perp_y_real*matrix_u[3] + \
            dder_plus_f_m_perp_z_real*matrix_u[6]
        dder_plus_f_m_perp_1_imag = \
            dder_plus_f_m_perp_x_imag*matrix_u[0] + \
            dder_plus_f_m_perp_y_imag*matrix_u[3] + \
            dder_plus_f_m_perp_z_imag*matrix_u[6]
        dder_plus_f_m_perp_2_real = \
            dder_plus_f_m_perp_x_real*matrix_u[1] + \
            dder_plus_f_m_perp_y_real*matrix_u[4] + \
            dder_plus_f_m_perp_z_real*matrix_u[7]
        dder_plus_f_m_perp_2_imag = \
            dder_plus_f_m_perp_x_imag*matrix_u[1] + \
            dder_plus_f_m_perp_y_imag*matrix_u[4] + \
            dder_plus_f_m_perp_z_imag*matrix_u[7]
        dder_plus_f_m_perp_3_real = \
            dder_plus_f_m_perp_x_real*matrix_u[2] + \
            dder_plus_f_m_perp_y_real*matrix_u[5] + \
            dder_plus_f_m_perp_z_real*matrix_u[8]
        dder_plus_f_m_perp_3_imag = \
            dder_plus_f_m_perp_x_imag*matrix_u[2] + \
            dder_plus_f_m_perp_y_imag*matrix_u[5] + \
            dder_plus_f_m_perp_z_imag*matrix_u[8]

        dder_plus["f_m_perp_real"] = numpy.stack([
            dder_plus_f_m_perp_1_real, dder_plus_f_m_perp_2_real, dder_plus_f_m_perp_3_real],
            axis=0) 
        dder_plus["f_m_perp_imag"] = numpy.stack([
            dder_plus_f_m_perp_1_imag, dder_plus_f_m_perp_2_imag, dder_plus_f_m_perp_3_imag],
            axis=0)


        dder_minus_f_m_perp_x_real = \
            y_m_perp_xy_f_m_perp_x_real * f_m_perp_xy_sq + \
            y_m_perp_xy * f_m_perp_xy_sq_f_m_perp_x_real -\
            p_d * chiral_term_f_m_perp_x_real
        dder_minus_f_m_perp_x_imag = \
            y_m_perp_xy_f_m_perp_x_imag * f_m_perp_xy_sq + \
            y_m_perp_xy * f_m_perp_xy_sq_f_m_perp_x_imag -\
            p_d * chiral_term_f_m_perp_x_imag
        dder_minus_f_m_perp_y_real = \
            y_m_perp_xy_f_m_perp_y_real * f_m_perp_xy_sq + \
            y_m_perp_xy * f_m_perp_xy_sq_f_m_perp_y_real -\
            p_d * chiral_term_f_m_perp_y_real
        dder_minus_f_m_perp_y_imag = \
            y_m_perp_xy_f_m_perp_y_imag * f_m_perp_xy_sq + \
            y_m_perp_xy * f_m_perp_xy_sq_f_m_perp_y_imag -\
            p_d * chiral_term_f_m_perp_y_imag
        dder_minus_f_m_perp_z_real =  0.5*(
            (1.-p_d)*(y_plus*f_plus_sq_f_m_perp_z_real+y_plus_f_m_perp_z_real*f_plus_sq) +
            (1.+p_d)*(y_minus*f_minus_sq_f_m_perp_z_real+y_minus_f_m_perp_z_real*f_minus_sq))
        dder_minus_f_m_perp_z_imag =  0.5*(
            (1.-p_d)*(y_plus*f_plus_sq_f_m_perp_z_imag+y_plus_f_m_perp_z_imag*f_plus_sq) +
            (1.+p_d)*(y_minus*f_minus_sq_f_m_perp_z_imag+y_minus_f_m_perp_z_imag*f_minus_sq))

        dder_minus_f_m_perp_1_real = \
            dder_minus_f_m_perp_x_real*matrix_u[0] + \
            dder_minus_f_m_perp_y_real*matrix_u[3] + \
            dder_minus_f_m_perp_z_real*matrix_u[6]
        dder_minus_f_m_perp_1_imag = \
            dder_minus_f_m_perp_x_imag*matrix_u[0] + \
            dder_minus_f_m_perp_y_imag*matrix_u[3] + \
            dder_minus_f_m_perp_z_imag*matrix_u[6]
        dder_minus_f_m_perp_2_real = \
            dder_minus_f_m_perp_x_real*matrix_u[1] + \
            dder_minus_f_m_perp_y_real*matrix_u[4] + \
            dder_minus_f_m_perp_z_real*matrix_u[7]
        dder_minus_f_m_perp_2_imag = \
            dder_minus_f_m_perp_x_imag*matrix_u[1] + \
            dder_minus_f_m_perp_y_imag*matrix_u[4] + \
            dder_minus_f_m_perp_z_imag*matrix_u[7]
        dder_minus_f_m_perp_3_real = \
            dder_minus_f_m_perp_x_real*matrix_u[2] + \
            dder_minus_f_m_perp_y_real*matrix_u[5] + \
            dder_minus_f_m_perp_z_real*matrix_u[8]
        dder_minus_f_m_perp_3_imag = \
            dder_minus_f_m_perp_x_imag*matrix_u[2] + \
            dder_minus_f_m_perp_y_imag*matrix_u[5] + \
            dder_minus_f_m_perp_z_imag*matrix_u[8]

        dder_minus["f_m_perp_real"] = numpy.stack([
            dder_minus_f_m_perp_1_real, dder_minus_f_m_perp_2_real, dder_minus_f_m_perp_3_real], 
            axis=0) 
        dder_minus["f_m_perp_imag"] = numpy.stack([
            dder_minus_f_m_perp_1_imag, dder_minus_f_m_perp_2_imag, dder_minus_f_m_perp_3_imag], 
            axis=0)

    extinction_keys = dder_y_plus.keys()
    if len(extinction_keys) != 0:
        for key in extinction_keys:
            if key == "f_sq":
                pass
            else:
                dder_plus[key] = \
                    0.5*((1.+p_u)*f_plus_sq*dder_y_plus[key] +
                         (1.-p_u)*f_minus_sq*dder_y_minus[key]) + \
                    f_m_perp_xy_sq * dder_y_m_perp_xy[key]
                dder_minus[key] = \
                    0.5*((1.-p_d)*f_plus_sq*dder_y_plus[key] +
                         (1.+p_d)*f_minus_sq*dder_y_minus[key]) + \
                    f_m_perp_xy_sq * dder_y_m_perp_xy[key]
    return iint_plus, iint_minus, dder_plus, dder_minus


def calc_flip_ratio_by_iint(
        iint_plus, iint_minus, c_lambda2: float = None, iint_2hkl = None,
        flag_iint_plus: bool = False, flag_iint_minus: bool = False, 
        flag_c_lambda2: bool = False, flag_iint_2hkl: bool = False):
    """Calculate flip ratio by given integrated intensities.
    """
    dder = {}
    cond = (c_lambda2 is None) or (iint_2hkl is None)
    if cond:
        signal_plus = iint_plus 
        signal_minus = iint_minus 
    else:
        signal_plus = iint_plus + c_lambda2*iint_2hkl
        signal_minus = iint_minus + c_lambda2*iint_2hkl

    flip_ratio = signal_plus/signal_minus

    if flag_iint_plus:
        dder["iint_plus"] = numpy.ones_like(iint_plus) / signal_minus
    if flag_iint_minus:
        dder["iint_minus"] = -1. * flip_ratio * numpy.ones_like(iint_minus)/ signal_minus
    if flag_c_lambda2:
        dder["c_lambda2"] = numpy.ones_like(c_lambda2)*iint_2hkl / signal_minus - \
            flip_ratio * numpy.ones_like(c_lambda2)*iint_2hkl/ signal_minus
    if flag_iint_2hkl:
        dder["iint_2hkl"] = numpy.ones_like(iint_2hkl)*c_lambda2 / signal_minus - \
            flip_ratio * numpy.ones_like(iint_2hkl)*c_lambda2/ signal_minus
    return flip_ratio, dder


def calc_asymmetry_by_iint(
        iint_plus, iint_minus, c_lambda2: float = None, iint_2hkl = None,
        flag_iint_plus: bool = False, flag_iint_minus: bool = False, 
        flag_c_lambda2: bool = False, flag_iint_2hkl: bool = False):
    """Calculate asymmetry (I^+ - I^-)/(I^+ + I^-).
    """
    dder = {}
    cond = (c_lambda2 is None) or (iint_2hkl is None)
    if cond:
        contamination = 0.
    else:
        contamination = 2.*c_lambda2*iint_2hkl

    denom = (iint_plus + iint_minus + contamination)
    asymmetry = (iint_plus - iint_minus) / denom

    if flag_iint_plus:
        dder["iint_plus"] = (1. / denom - (iint_plus - iint_minus) / numpy.square(denom)) * \
            numpy.ones_like(iint_plus)
    if flag_iint_minus:
        dder["iint_minus"] = (-1. / denom - (iint_plus - iint_minus) / numpy.square(denom)) * \
            numpy.ones_like(iint_minus)
    if flag_c_lambda2:
        dder["c_lambda2"] = 2.*iint_2hkl*( - (iint_plus - iint_minus) / numpy.square(denom)) * \
            numpy.ones_like(c_lambda2)
    if flag_iint_2hkl:
        dder["iint_2hkl"] = 2.*c_lambda2*( - (iint_plus - iint_minus) / numpy.square(denom)) * \
            numpy.ones_like(iint_2hkl)
    return asymmetry, dder


def calc_intensities_by_structure_factors(
        beam_polarization: float, flipper_efficiency: float, f_nucl, f_m_perp, matrix_u, func_extinction: Callable = None,
        c_lambda2: float=None, f_nucl_2hkl=None, f_m_perp_2hkl=None,
        flag_beam_polarization: bool = False, flag_flipper_efficiency: bool = False,
        flag_f_nucl: bool = False, flag_f_m_perp: bool = False,
        flag_c_lambda2: bool = False,
        flag_f_nucl_2hkl: bool = False, flag_f_m_perp_2hkl: bool = False,
        dict_in_out: dict = None, flag_use_precalculated_data: bool = False):
    """Calculate flip ratio by given structure factors
    """
    if dict_in_out is None:
        flag_dict = False
        dict_in_out_keys = []
    else:
        flag_dict = True
        dict_in_out_keys = dict_in_out.keys()
    
    iint_plus, iint_minus, dder_plus_hkl, dder_minus_hkl = calc_iint(
        beam_polarization, flipper_efficiency, f_nucl, f_m_perp, matrix_u, func_extinction=func_extinction,
        flag_beam_polarization=flag_beam_polarization, flag_flipper_efficiency=flag_flipper_efficiency,
        flag_f_nucl=flag_f_nucl, flag_f_m_perp=flag_f_m_perp,
        dict_in_out=dict_in_out, flag_use_precalculated_data=flag_use_precalculated_data)

    cond = (c_lambda2 is None) or (f_nucl_2hkl is None) or (f_m_perp_2hkl is None)
    if isinstance(c_lambda2, numpy.ndarray) and not(cond):
        cond = numpy.any(numpy.isnan(c_lambda2))
    if cond:
        c_lambda2 = None
        iint_2hkl = None
        flag_iint_2hkl = False
        dder_2hkl = {}
    else:
        iint_plus_2hkl, iint_minus_2hkl, dder_plus_2hkl, dder_minus_2hkl = calc_iint(
            0, 0, f_nucl_2hkl, f_m_perp_2hkl, matrix_u, func_extinction=func_extinction,
            flag_beam_polarization=False, flag_flipper_efficiency=False,
            flag_f_nucl=flag_f_nucl_2hkl, flag_f_m_perp=flag_f_m_perp_2hkl,
            dict_in_out=None, flag_use_precalculated_data=False)
        iint_2hkl = iint_plus_2hkl
        dder_2hkl = dder_plus_2hkl
        flag_iint_2hkl = len(dder_2hkl.keys()) > 0
        print("here")
        iint_plus += c_lambda2*iint_2hkl
        iint_minus += c_lambda2*iint_2hkl

        if flag_f_nucl_2hkl:
            dder_2hkl["f_nucl_2hkl_real"] = dder_2hkl.pop("f_nucl_real")
            dder_2hkl["f_nucl_2hkl_imag"] = dder_2hkl.pop("f_nucl_imag")
        if flag_f_m_perp_2hkl:
            dder_2hkl["f_m_perp_2hkl_real"] = dder_2hkl.pop("f_m_perp_real")
            dder_2hkl["f_m_perp_2hkl_imag"] = dder_2hkl.pop("f_m_perp_imag")
    
    dder_plus = {}
    dder_minus = {}
    keys_plus = dder_plus_hkl.keys()
    keys_minus = dder_minus_hkl.keys()
    keys_2h = dder_2hkl.keys()
    set_key = set(keys_plus) | set(keys_minus) | set(keys_2h)
    for key in set_key:
        flag_1 = key in keys_plus
        flag_2 = key in keys_minus
        flag_3 = key in keys_2h
        if flag_1:
            term_1 = dder_plus_hkl[key]
        if flag_2:
            term_2 = dder_minus_hkl[key]
        if flag_3:
            term_3 = dder_2hkl[key]

        if flag_1 and flag_3:
            dder_plus[key] = term_1 + c_lambda2 * term_3
        elif flag_1:
            dder_plus[key] = term_1 
        elif flag_3:
            dder_plus[key] = c_lambda2 * term_3
        
        if flag_2 and flag_3:
            dder_minus[key] = term_2 + c_lambda2 * term_3
        elif flag_2:
            dder_minus[key] = term_2
        elif flag_3:
            dder_minus[key] = c_lambda2 * term_3

    if flag_c_lambda2:
        dder_plus["c_lambda2"] = iint_2hkl
        dder_minus["c_lambda2"] = iint_2hkl
    return iint_plus, iint_minus, dder_plus, dder_minus


def calc_flip_ratio_by_structure_factors(
        beam_polarization: float, flipper_efficiency: float, f_nucl, f_m_perp, matrix_u, func_extinction: Callable = None,
        c_lambda2: float=None, f_nucl_2hkl=None, f_m_perp_2hkl=None,
        flag_beam_polarization: bool = False, flag_flipper_efficiency: bool = False,
        flag_f_nucl: bool = False, flag_f_m_perp: bool = False,
        flag_c_lambda2: bool = False,
        flag_f_nucl_2hkl: bool = False, flag_f_m_perp_2hkl: bool = False,
        dict_in_out: dict = None, flag_use_precalculated_data: bool = False):
    """Calculate flip ratio by given structure factors
    """
    if dict_in_out is None:
        flag_dict = False
        dict_in_out_keys = []
    else:
        flag_dict = True
        dict_in_out_keys = dict_in_out.keys()

    iint_plus, iint_minus, dder_plus, dder_minus = calc_intensities_by_structure_factors(
        beam_polarization, flipper_efficiency, f_nucl, f_m_perp, matrix_u, func_extinction=func_extinction,
        c_lambda2=c_lambda2, f_nucl_2hkl=f_nucl_2hkl, f_m_perp_2hkl=f_m_perp_2hkl,
        flag_beam_polarization=flag_beam_polarization, flag_flipper_efficiency=flag_flipper_efficiency,
        flag_f_nucl=flag_f_nucl, flag_f_m_perp=flag_f_m_perp,
        flag_c_lambda2=flag_c_lambda2,
        flag_f_nucl_2hkl=flag_f_nucl_2hkl, flag_f_m_perp_2hkl=flag_f_m_perp_2hkl,
        dict_in_out=dict_in_out, flag_use_precalculated_data=flag_use_precalculated_data)

    flip_ration = iint_plus/iint_minus

    dder = {}
    keys_plus = dder_plus.keys()
    keys_minus = dder_minus.keys()
    set_key = set(keys_plus) | set(keys_minus)
    for key in set_key:
        flag_1 = key in keys_plus
        flag_2 = key in keys_minus
        if flag_1 and flag_2:
            dder[key] = dder_plus[key]/iint_minus - flip_ration*dder_minus[key]/iint_minus
        elif flag_1:
            dder[key] = dder_plus[key]/iint_minus
        elif flag_2:
            dder[key] = - flip_ration*dder_minus[key]/iint_minus
    return flip_ration, dder


