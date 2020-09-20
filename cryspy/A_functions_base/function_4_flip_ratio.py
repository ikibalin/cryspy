import numpy

from .function_1_algebra import calc_modulus_sq_by_complex_vector,\
    calc_scalar_product_by_complex_vectors

from .function_2_crystallography_base import calc_volume_uc_by_abc_cosines, \
    calc_sthovl_by_hkl_abc_cosines


from .function_3_extinction import calc_extinction_2


def calc_f_plus_sq(f_nucl:complex, f_mag_perp_vert:complex, 
                   flag_f_nucl:bool=False, flag_f_mag_perp_vert:bool=False):
    """
$$
F_{+}^{2} = (F_N + F_{M, \perp, z})^{2}
$$
    """
    flag_derivatives = any([flag_f_nucl, flag_f_mag_perp_vert])
    f_plus = f_nucl + f_mag_perp_vert
    f_plus_sq = (f_plus*f_plus.conjugate()).real
    dder = {}
    if flag_derivatives:
        if flag_f_nucl:
            der_f_plus_sq__f_nucl = 2.*f_plus*(1.+0.*f_nucl)
            dder["der__f_nucl_re"] = der_f_plus_sq__f_nucl.real
            dder["der__f_nucl_im"] = der_f_plus_sq__f_nucl.imag
        if flag_f_mag_perp_vert:
            der_f_plus_sq__f_mag_perp_vert = 2.*f_plus*(1.+0.*f_mag_perp_vert)
            dder["der__f_mag_perp_vert_re"] = der_f_plus_sq__f_mag_perp_vert.real
            dder["der__f_mag_perp_vert_im"] = der_f_plus_sq__f_mag_perp_vert.imag
    return f_plus_sq, dder  


def calc_f_minus_sq(f_nucl:complex, f_mag_perp_vert:complex, 
                   flag_f_nucl:bool=False, flag_f_mag_perp_vert:bool=False):
    """
$$
F_{-}^{2} = (F_N - F_{M, \perp, z})^{2}
$$
    """
    flag_derivatives = any([flag_f_nucl, flag_f_mag_perp_vert])
    f_minus = f_nucl - f_mag_perp_vert
    f_minus_sq = (f_minus*f_minus.conjugate()).real
    dder = {}
    if flag_derivatives:
        if flag_f_nucl:
            der_f_minus_sq__f_nucl = 2.*f_minus*(1.+0.*f_nucl)
            dder["der__f_nucl_re"] = der_f_minus_sq__f_nucl.real
            dder["der__f_nucl_im"] = der_f_minus_sq__f_nucl.imag
        if flag_f_mag_perp_vert:
            der_f_plus_sq__f_mag_perp_vert = 2.*f_minus*(1.+0.*f_mag_perp_vert)
            dder["der__f_mag_perp_vert_re"] = der_f_minus_sq__f_mag_perp_vert.real
            dder["der__f_mag_perp_vert_im"] = der_f_minus_sq__f_mag_perp_vert.imag
    return f_minus_sq, dder


def calc_flip_ratio(f_nucl:float, f_mag_perp:tuple, 
                    e_up:tuple, beam_polarization:float, flipper_efficiency:float, 
                    radius:float, mosaicity:float, model_extinction:str,
                    volume_unit_cell:float, sthovl:float,
                    wavelength:float,
                    flag_f_nucl=False, flag_f_mag_perp=False, flag_e_up=False, 
                    flag_beam_polarization=False, flag_flipper_efficiency=False,
                    flag_radius=False, flag_mosaicity=False,
                    flag_volume_unit_cell=False, flag_sthovl=False,
                    flag_volume_wavelength=False):
    """
    Calculate flip ratio by given nuclear structure factor and 
    perpendicular component of magnetic structure factor 
    (in 10**-12 cm)
    

    cell (only for extinction correction)
    f_nucl [hkl]
    f_m_perp [hkl]
    delta_f_nucl [hkl, parameters]
    delta_f_m_perp [hkl, parameters]
    """
    flag_derivatives = any([flag_f_nucl, flag_f_mag_perp, flag_e_up, 
                    flag_beam_polarization, flag_flipper_efficiency,
                    flag_radius, flag_mosaicity,
                    flag_volume_unit_cell, flag_sthovl,
                    flag_volume_wavelength])
              
    np_b_p = numpy.array(beam_polarization, dtype=float)              
    np_f_e = numpy.array(flipper_efficiency, dtype=float)              
    p_u = np_b_p
    p_d = (2.* np_f_e - 1.) * p_u

    part_u_u, part_u_d = 0.5*(1. + p_u), 0.5*(1. - p_u)
    part_d_u, part_d_d = 0.5*(1. - p_d), 0.5*(1. + p_d)

    f_mag_perp_vert, dder_f_mag_perp_vert = calc_scalar_product_by_complex_vectors(
                    f_mag_perp, e_up, flag_vector_1=flag_f_mag_perp, flag_vector_2=flag_e_up)
    if flag_f_mag_perp:
        der_f_mag_perp_vert__f_mag_perp_re = dder_f_mag_perp_vert["der__vector_1_re"] 
        der_f_mag_perp_vert__f_mag_perp_im = dder_f_mag_perp_vert["der__vector_1_im"] 

    if flag_e_up:
        der_f_mag_perp_vert__e_up = dder_f_mag_perp_vert["der__vector_2_re"] 
            
    flag_f_mag_perp_vert = any([flag_f_mag_perp, flag_e_up])

    f_plus_sq, dder_f_plus_sq = calc_f_plus_sq(f_nucl, f_mag_perp_vert, flag_f_nucl, flag_f_mag_perp_vert)
    f_minus_sq, dder_f_minus_sq = calc_f_minus_sq(f_nucl, f_mag_perp_vert, flag_f_nucl, flag_f_mag_perp_vert)

    if flag_f_nucl:
        der_f_plus_sq__f_nucl_re = dder_f_plus_sq["der__f_nucl_re"]
        der_f_plus_sq__f_nucl_im = dder_f_plus_sq["der__f_nucl_im"]
        der_f_minus_sq__f_nucl_re = dder_f_minus_sq["der__f_nucl_re"]
        der_f_minus_sq__f_nucl_im = dder_f_minus_sq["der__f_nucl_im"]
    if flag_f_mag_perp_vert:
        der_f_plus_sq__f_mag_perp_vert_re = dder_f_plus_sq["der__f_mag_perp_vert_re"]
        der_f_plus_sq__f_mag_perp_vert_im = dder_f_plus_sq["der__f_mag_perp_vert_im"]
        der_f_minus_sq__f_mag_perp_vert_re = dder_f_minus_sq["der__f_mag_perp_vert_re"]
        der_f_minus_sq__f_mag_perp_vert_im = dder_f_minus_sq["der__f_mag_perp_vert_im"]
    
    f_mag_perp_sq, dder_f_mag_perp_sq = calc_modulus_sq_by_complex_vector(f_mag_perp, flag_f_mag_perp)

    if flag_f_mag_perp:
        der_f_mag_perp_sq__f_mag_perp_re = dder_f_mag_perp_sq["der__vector_1_re"]
        der_f_mag_perp_sq__f_mag_perp_im = dder_f_mag_perp_sq["der__vector_1_im"]

    f_mag_perp_vert_sq = (f_mag_perp_vert*f_mag_perp_vert.conjugate()).real
    if flag_f_mag_perp_vert:
        der_f_mag_perp_vert_sq__f_mag_perp_vert_re = 2.*f_mag_perp_vert.real
        der_f_mag_perp_vert_sq__f_mag_perp_vert_im = 2.*f_mag_perp_vert.imag

    f_mag_sq = f_mag_perp_sq - f_mag_perp_vert_sq
    
    if any([flag_f_mag_perp, flag_f_mag_perp_vert]):
        der_f_mag_sq__f_mag_perp_re = der_f_mag_perp_sq__f_mag_perp_re - \
            der_f_mag_perp_vert_sq__f_mag_perp_vert_re * der_f_mag_perp_vert__f_mag_perp_re.real + \
            der_f_mag_perp_vert_sq__f_mag_perp_vert_im * der_f_mag_perp_vert__f_mag_perp_re.imag 
    
        der_f_mag_sq__f_mag_perp_im = der_f_mag_perp_sq__f_mag_perp_im - \
            der_f_mag_perp_vert_sq__f_mag_perp_vert_re * der_f_mag_perp_vert__f_mag_perp_im.real + \
            der_f_mag_perp_vert_sq__f_mag_perp_vert_im * der_f_mag_perp_vert__f_mag_perp_im.imag 
    
    
    #extinction correction
    y_plus, dder_y_plus = calc_extinction_2(radius, mosaicity, model_extinction,
                            f_plus_sq, 
                            volume_unit_cell, sthovl, wavelength)
    y_minus, dder_y_minus = calc_extinction_2(radius, mosaicity, model_extinction,
                            f_minus_sq, 
                            volume_unit_cell, sthovl, wavelength)
    y_mag, dder_y_mag = calc_extinction_2(radius, mosaicity, model_extinction,
                            f_mag_sq, 
                            volume_unit_cell, sthovl, wavelength)
    if flag_derivatives:
        der_y_plus__f_plus_sq = dder_y_plus["der_yext__f_sq"]
        der_y_minus__f_minus_sq = dder_y_minus["der_yext__f_sq"]
        der_y_mag__f_mag_sq = dder_y_mag["der_yext__f_sq"]
                            
    int_up = part_u_u*y_plus*f_plus_sq + part_u_d*y_minus*f_minus_sq + y_mag*f_mag_sq
    int_down = part_d_u*y_plus*f_plus_sq + part_d_d*y_minus*f_minus_sq + y_mag*f_mag_sq

    if flag_derivatives:
        der_int_up__f_plus_sq = part_u_u*(1.+0.*f_plus_sq)*y_plus
        der_int_up__f_minus_sq = part_u_d*(1.+0.*f_minus_sq)*y_minus
        der_int_up__f_mag_sq = (1.+0.*f_mag_sq)*y_mag
        der_int_up__y_plus = part_u_u*f_plus_sq*(1.+0.*y_plus)
        der_int_up__y_minus = part_u_d*f_minus_sq*(1.+0.*y_minus)
        der_int_up__y_mag = f_mag_sq*(1.+0.*y_mag)
        
        der_int_down__f_plus_sq = part_d_u*(1.+0.*f_plus_sq)*y_plus
        der_int_down__f_minus_sq = part_d_d*(1.+0.*f_minus_sq)*y_minus
        der_int_down__f_mag_sq = der_int_up__f_mag_sq
        der_int_down__y_plus = part_d_u*f_plus_sq*(1.+0.*y_plus)
        der_int_down__y_minus = part_d_d*f_minus_sq*(1.+0.*y_minus)
        der_int_down__y_mag = der_int_up__y_mag
        
    flip_ratio = int_up / int_down

    dder = {}
    #derivatives
    if flag_derivatives:
        der_flip_ratio__int_up = (1.+0* int_up) / int_down
        der_flip_ratio__int_down = - flip_ratio / int_down

    if flag_beam_polarization:
        tt = y_plus*f_plus_sq - y_minus*f_minus_sq
        der_int_up__beam_polarization = 0.5*tt*(1.+0.*np_b_p)
        der_int_down__beam_polarization = -0.5*(2.* np_f_e - 1.)*tt*(1.+0.*np_b_p)
        der_flip_ratio__beam_polarization = der_flip_ratio__int_up*der_int_up__beam_polarization + \
                                            der_flip_ratio__int_down*der_int_down__beam_polarization
        dder["der_flip_ratio__beam_polarization"] = der_flip_ratio__beam_polarization
        
    if flag_flipper_efficiency:
        tt = y_plus*f_plus_sq - y_minus*f_minus_sq
        der_int_down__flipper_efficiency = -1.0*(0.* np_f_e + 1.)*tt*np_b_p
        der_flip_ratio__flipper_efficiency = der_flip_ratio__int_down*der_int_down__flipper_efficiency
        dder["der_flip_ratio__flipper_efficiency"] = der_flip_ratio__flipper_efficiency

    if flag_f_mag_perp:
        der_f_plus_sq__f_mag_perp_re =  der_f_plus_sq__f_mag_perp_vert_re * der_f_mag_perp_vert__f_mag_perp_re.real + \
                                        der_f_plus_sq__f_mag_perp_vert_im * der_f_mag_perp_vert__f_mag_perp_re.imag
        
        der_f_minus_sq__f_mag_perp_re =  der_f_minus_sq__f_mag_perp_vert_re * der_f_mag_perp_vert__f_mag_perp_re.real + \
                                         der_f_minus_sq__f_mag_perp_vert_im * der_f_mag_perp_vert__f_mag_perp_re.imag

        der_int_up__f_mag_perp_re = der_int_up__f_plus_sq * der_f_plus_sq__f_mag_perp_re + \
                                    der_int_up__f_minus_sq * der_f_minus_sq__f_mag_perp_re + \
                                    der_int_up__y_plus * der_y_plus__f_plus_sq * der_f_plus_sq__f_mag_perp_re + \
                                    der_int_up__y_minus * der_y_minus__f_minus_sq * der_f_minus_sq__f_mag_perp_re + \
                                    der_int_up__y_mag * der_y_mag__f_mag_sq * der_f_mag_sq__f_mag_perp_re

        der_f_plus_sq__f_mag_perp_im =  der_f_plus_sq__f_mag_perp_vert_re * der_f_mag_perp_vert__f_mag_perp_im.real + \
                                        der_f_plus_sq__f_mag_perp_vert_im * der_f_mag_perp_vert__f_mag_perp_im.imag
        
        der_f_minus_sq__f_mag_perp_im =  der_f_minus_sq__f_mag_perp_vert_re * der_f_mag_perp_vert__f_mag_perp_im.real + \
                                         der_f_minus_sq__f_mag_perp_vert_im * der_f_mag_perp_vert__f_mag_perp_im.imag

        der_int_up__f_mag_perp_im = der_int_up__f_plus_sq * der_f_plus_sq__f_mag_perp_im + \
                                    der_int_up__f_minus_sq * der_f_minus_sq__f_mag_perp_im + \
                                    der_int_up__y_plus * der_y_plus__f_plus_sq * der_f_plus_sq__f_mag_perp_im + \
                                    der_int_up__y_minus * der_y_minus__f_minus_sq * der_f_minus_sq__f_mag_perp_im + \
                                    der_int_up__y_mag * der_y_mag__f_mag_sq * der_f_mag_sq__f_mag_perp_im

        der_int_down__f_mag_perp_re = der_int_down__f_plus_sq * der_f_plus_sq__f_mag_perp_re + \
                                      der_int_down__f_minus_sq * der_f_minus_sq__f_mag_perp_re + \
                                      der_int_down__y_plus * der_y_plus__f_plus_sq * der_f_plus_sq__f_mag_perp_re + \
                                      der_int_down__y_minus * der_y_minus__f_minus_sq * der_f_minus_sq__f_mag_perp_re + \
                                      der_int_down__y_mag * der_y_mag__f_mag_sq * der_f_mag_sq__f_mag_perp_re

        der_int_down__f_mag_perp_im = der_int_down__f_plus_sq * der_f_plus_sq__f_mag_perp_im + \
                                      der_int_down__f_minus_sq * der_f_minus_sq__f_mag_perp_im + \
                                      der_int_down__y_plus * der_y_plus__f_plus_sq * der_f_plus_sq__f_mag_perp_im + \
                                      der_int_down__y_minus * der_y_minus__f_minus_sq * der_f_minus_sq__f_mag_perp_im + \
                                      der_int_down__y_mag * der_y_mag__f_mag_sq * der_f_mag_sq__f_mag_perp_im
        
        der_flip_ratio__f_mag_perp_re = der_flip_ratio__int_up*der_int_up__f_mag_perp_re + \
                                        der_flip_ratio__int_down*der_int_down__f_mag_perp_re
        der_flip_ratio__f_mag_perp_im = der_flip_ratio__int_up*der_int_up__f_mag_perp_im + \
                                        der_flip_ratio__int_down*der_int_down__f_mag_perp_im
        dder["der_flip_ratio__f_mag_perp_re"] = der_flip_ratio__f_mag_perp_re
        dder["der_flip_ratio__f_mag_perp_im"] = der_flip_ratio__f_mag_perp_im
        
    return flip_ratio, dder


"""    
    f_n_sq = (f_nucl * f_nucl.conjugate()).real
    f_m_p_x, f_m_p_y, f_m_p_z = f_m_perp
    
        
    f_m_p_sq = (f_m_p_x*f_m_p_x.conjugate() + 
                f_m_p_y*f_m_p_y.conjugate() + 
                f_m_p_z*f_m_p_z.conjugate()).real


        
    if delta_f_m_perp is not None:
        delta_f_m_p_x, delta_f_m_p_y, delta_f_m_p_z = delta_f_m_perp
        delta_f_m_p_vert = scalar_product(delta_f_m_perp, e_up)
            
        
    two_f_nucl_f_m_p_vert = (f_nucl * f_m_p_vert.conjugate() + 
                             f_m_p_vert * f_nucl.conjugate()).real 

    if delta_f_m_perp is not None:
        delta_two_f_nucl_f_m_p_vert = (
            f_nucl[:, numpy.newaxis]*delta_f_m_p_vert.conjugate() +
            delta_f_m_p_vert*f_nucl.conjugate()[:, numpy.newaxis]).real


        delta_f_m_p_sq = 2.*(
            numpy.real(f_m_p_x)[:,numpy.newaxis]*numpy.real(delta_f_m_p_x)+
            numpy.imag(f_m_p_x)[:,numpy.newaxis]*numpy.imag(delta_f_m_p_x)+
            numpy.real(f_m_p_y)[:,numpy.newaxis]*numpy.real(delta_f_m_p_y)+
            numpy.imag(f_m_p_y)[:,numpy.newaxis]*numpy.imag(delta_f_m_p_y)+
            numpy.real(f_m_p_z)[:,numpy.newaxis]*numpy.real(delta_f_m_p_z)+
            numpy.imag(f_m_p_z)[:,numpy.newaxis]*numpy.imag(delta_f_m_p_z))

    h, k, l = hkl

    # extinction correction
    f_m_p_vert_sq = (f_m_p_vert * f_m_p_vert.conjugate()).real

    if delta_f_m_perp is not None:
        
        delta_f_m_p_vert_sq = 2.*(
         numpy.real(f_m_p_vert)[:,numpy.newaxis]*numpy.real(delta_f_m_p_vert)+
         numpy.imag(f_m_p_vert)[:,numpy.newaxis]*numpy.imag(delta_f_m_p_vert)) 

    fnp = two_f_nucl_f_m_p_vert
    fp_sq = f_n_sq + f_m_p_sq + two_f_nucl_f_m_p_vert
    fm_sq = f_n_sq + f_m_p_sq - two_f_nucl_f_m_p_vert
    fpm_sq = f_m_p_sq - f_m_p_vert_sq

    if delta_f_m_perp is not None:
        delta_fnp = delta_two_f_nucl_f_m_p_vert
        delta_fp_sq = delta_f_m_p_sq + delta_two_f_nucl_f_m_p_vert 
        delta_fm_sq = delta_f_m_p_sq - delta_two_f_nucl_f_m_p_vert 
        delta_fpm_sq = delta_f_m_p_sq - delta_f_m_p_vert_sq 
        
        if extinction is not None:
            yp, delta_yp = extinction.calc_extinction(cell, h, k, l, 
                        fp_sq, wavelength, flag_derivative_f_sq=True)
            ym, delta_ym = extinction.calc_extinction(cell, h, k, l, 
                        fm_sq, wavelength, flag_derivative_f_sq=True)
            ypm, delta_ypm = extinction.calc_extinction(cell, h, k, l, 
                        fpm_sq, wavelength, flag_derivative_f_sq=True)
        else:
            yp = numpy.ones(shape=fp_sq.shape, dtype=float), 
            delta_yp = numpy.zeros(shape=fp_sq.shape, dtype=float)
            ym = numpy.ones(shape=fp_sq.shape, dtype=float)
            delta_ym = numpy.zeros(shape=fp_sq.shape, dtype=float)
            ypm = numpy.ones(shape=fp_sq.shape, dtype=float)
            delta_ypm = numpy.zeros(shape=fp_sq.shape, dtype=float)

        pppl = 0.5 * ((1 + p_u) * yp + (1 - p_u) * ym)
        ppmin = 0.5 * ((1 - p_d) * yp + (1 + p_d) * ym)
        pmpl = 0.5 * ((1 + p_u) * yp - (1 - p_u) * ym)
        pmmin = 0.5 * ((1 - p_d) * yp - (1 + p_d) * ym)

        delta_pppl = 0.5*((1+p_u)*delta_yp[:,numpy.newaxis]*delta_fp_sq+\
                     (1-p_u)*delta_ym[:,numpy.newaxis] * delta_fm_sq)
        delta_ppmin = 0.5*((1-p_d)*delta_yp[:,numpy.newaxis]*delta_fp_sq+\
                     (1+p_d)*delta_ym[:,numpy.newaxis]*delta_fm_sq)
        delta_pmpl = 0.5*((1+p_u)*delta_yp[:,numpy.newaxis]*delta_fp_sq-\
                     (1-p_u)*delta_ym[:,numpy.newaxis]*delta_fm_sq)
        delta_pmmin = 0.5*((1-p_d)*delta_yp[:,numpy.newaxis]*delta_fp_sq-\
                     (1+p_d)*delta_ym[:,numpy.newaxis]*delta_fm_sq)

    iint_u = (f_n_sq+f_m_p_vert_sq)*pppl+pmpl*fnp+ypm*fpm_sq
    iint_d = (f_n_sq+f_m_p_vert_sq)*ppmin+pmmin*fnp+ypm*fpm_sq

    delta_iint_u = (((f_n_sq + f_m_p_vert_sq)[:, numpy.newaxis] * delta_pppl + delta_pmpl * fnp[:,
                           numpy.newaxis] + (delta_ypm * fpm_sq)[:,numpy.newaxis] * delta_fpm_sq) +
                           (delta_f_m_p_vert_sq) * pppl[:, numpy.newaxis] +
                           pmpl[:, numpy.newaxis] * delta_fnp +
                           ypm[:, numpy.newaxis] * delta_fpm_sq)
                       
    delta_iint_d = (((f_n_sq + f_m_p_vert_sq)[:, numpy.newaxis] * delta_ppmin + delta_pmmin * fnp[:,
                                numpy.newaxis] + (delta_ypm * fpm_sq)[:,numpy.newaxis] * delta_fpm_sq) +
                           (delta_f_m_p_vert_sq) * ppmin[:, numpy.newaxis] +
                           pmmin[:, numpy.newaxis] * delta_fnp +
                           ypm[:, numpy.newaxis] * delta_fpm_sq)

    f_r_m = iint_u/iint_d

    delta_f_r_m = (iint_d[:,numpy.newaxis]*delta_iint_u - 
                       iint_u[:,numpy.newaxis]*delta_iint_d) / \
                       (iint_d**2)[:,numpy.newaxis]
                                                                                                  
    return f_r_m, delta_f_r_m
"""
