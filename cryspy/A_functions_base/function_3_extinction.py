"""
Module gives correction on extinction in case of isotropical case.
"""

import numpy

from .function_2_crystallography_base import \
    calc_volume_uc_by_abc_cosines,\
    calc_sthovl_by_hkl_abc_cosines

def calc_extinction(radius:float, mosaicity:float, model:str,
                    a:float, b:float, c:float, alpha:float, beta:float, gamma:float,
                    h:float, k:float, l:float, 
                    f_sq:float, wavelength:float, flag_derivative_f_sq=False):
    """
    Isotropical extinction coorection y:

    $$
    |F|^2_{\text{corrected}} = y \\cdot |F|^2
    $$

    radius primary extinction ???
    mosaisity secondary extinction
    model= "gauss" or "lorentz"
    a,b,c,alpha,beta,gamma are unit cell parameters (in angstrem and radians)
    h, k, l are Miller indices
    f_sq is square of structure factor (in 10-12cm)
    wavelength is neutron wavelength in angstrems
    flag_derivative_radius
    flag_derivative_mosaicity
    flag_derivative_a
    flag_derivative_b
    flag_derivative_c
    flag_derivative_alpha
    flag_derivative_beta
    flag_derivative_gamma
    flag_derivative_f_sq
    flag_derivative_wavelength
    """
    r = float(radius)
    g = float(mosaicity)
    g_sq = numpy.square(g)

    kk = 1.
    c_a, c_b, c_g = numpy.cos(alpha), numpy.cos(beta), numpy.cos(gamma)
    
    volume_unit_cell = calc_volume_uc_by_abc_cosines(a, b, c, c_a, c_b, c_g)
    
    sthovl = calc_sthovl_by_hkl_abc_cosines(h, k, l, a, b, c, c_a, c_b, c_g)

    yext, dder  = calc_extinction_2(radius, mosaicity, model,
                                f_sq, volume_unit_cell, sthovl, wavelength) 
    return yext, dder 
    

def calc_extinction_2(radius:float, mosaicity:float, model:str,
                    f_sq:float, volume_unit_cell:float, sthovl:float, 
                    wavelength:float):
    """
    radius in micro meters
    g in minutes
    f_sq in 10**-12cm
    volume_unit_cell in angstrem**3
    wavelength in angstrem
    sthovl in angstrem**-1
    """
    r = radius
    g = mosaicity*numpy.pi/(180*60) #transfer minutes to radians
    g_sq = numpy.square(g)
    kk = 1.
    volume_sq = numpy.square(volume_unit_cell)
    stheta = sthovl * wavelength

    s2theta = 2.*stheta*numpy.sqrt(1.-numpy.square(stheta))
    c2theta = 1.-2.*numpy.square(stheta)
    c2theta_sq = numpy.square(c2theta)

    wavelength_cube = wavelength**3
    
    q = (f_sq*kk/volume_sq)*(wavelength_cube)*1./s2theta 
    delta_f_sq_q = (kk/volume_sq)*(wavelength_cube)*1./s2theta 

    t = 1.5*r
    alpha = 1.5*r*s2theta*1./wavelength
    x = 2./3*q*alpha*t

    x_sq = numpy.square(x)
    delta_f_sq_x = 2./3*delta_f_sq_q*alpha*t

    A = 0.20 + 0.45 * c2theta
    B = 0.22 - 0.12 * numpy.square(0.5-c2theta)
    yp = 1./(numpy.sqrt(1.+2.*x+(A*x_sq)*1./(1.+B*x)))
    yp_cube = yp**3
    delta_f_sq_yp = -0.5*(yp_cube)*(2.*delta_f_sq_x+
                      2.*A*delta_f_sq_x*x*1./(1.+B*x)-
                      B*delta_f_sq_x*(A*x_sq)*1./(numpy.square(1.+B*x))) 

    ag = 0.*sthovl 
    al = 0.*sthovl

    flag = alpha != 0.
    ag[flag] = alpha[flag]*g*(g_sq+0.5*alpha[flag]**2)**(-0.5)
    al[flag] = alpha[flag]*g*1./(g+alpha[flag]*2./3.)

    if model == "gauss":
        xs = 2./3.*q*ag*t
        xs_sq = numpy.square(xs)
        delta_f_sq_xs = 2./3.*delta_f_sq_q*ag*t
        A = 0.58 + 0.48 * c2theta + 0.24 * c2theta_sq
        B = 0.02 - 0.025 * c2theta
        #print("A, B", A, B)
        ys = 1./numpy.sqrt(1+2.12*xs+(A*xs_sq)*1./(1+B*xs))
        ys_cube = ys**3
        delta_f_sq_ys = -0.5*ys_cube*\
                        (2.12*delta_f_sq_xs+\
                        2.*(A*xs*delta_f_sq_xs)*1./(1+B*xs)-\
                        B*delta_f_sq_xs*(A*xs_sq)*1./(numpy.square(1+B*xs)))
    elif model == "lorentz":
        xs = 2./3.*q*al*t
        xs_sq = numpy.square(xs)
        delta_f_sq_xs = 2./3.*delta_f_sq_q*al*t

        A = 0.025 + 0.285 * c2theta
        B = -0.45 * c2theta
        flag = c2theta>0
        B[flag] = 0.15 - 0.2 * numpy.square(0.75-c2theta[flag])
        ys = 1./numpy.sqrt(1.+2.*xs+(A*xs_sq)*1./(1+B*xs))
        ys_cube = ys**3
        delta_f_sq_ys = -0.5*ys_cube*\
                        (2.*delta_f_sq_xs+\
                        2.*(A*xs*delta_f_sq_xs)*1./(1+B*xs)-\
                        B*delta_f_sq_xs*(A*xs_sq)*1./(numpy.square(1.+B*xs)))
    else:
        ys = 1.
        delta_f_sq_ys = 0.
    #print("ys", ys)
    yext = yp * ys
    delta_f_sq_yext = delta_f_sq_yp * ys + yp * delta_f_sq_ys

    dder = {"der_yext__f_sq": delta_f_sq_yext}

    return yext, dder 
