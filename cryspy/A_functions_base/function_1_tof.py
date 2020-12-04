""" TOF

Expressions are taken from paper

J. Appl. Cryst. (1982). 15, 581-589


Functions
---------
    - 
    
"""
import numpy
import scipy
import scipy.special

inv_8ln2 = 0.18033688011112042591999058512524

na = numpy.newaxis
square = numpy.square
exp = numpy.exp
imag = numpy.imag
erfc = scipy.special.erfc
exp1 = scipy.special.exp1  # E1 = \int_{z}^{\inf} exp(-t)/t dt


def calc_hpv_eta(h_g, h_l):
    """pseudo-Voight function
    
    calculate h_pV and eta based on Gauss and Lorentz Size
    """
    h_g_2, h_l_2 = h_g**2, h_l**2
    h_g_3, h_l_3 = h_g_2*h_g, h_l_2*h_l
    h_g_4, h_l_4 = h_g_3*h_g, h_l_3*h_l
    h_g_5, h_l_5 = h_g_4*h_g, h_l_4*h_l
    c_2, c_3, c_4, c_5 = 2.69269, 2.42843, 4.47163, 0.07842
    h_pv = (h_g_5 + c_2*h_g_4*h_l + c_3*h_g_3*h_l_2 + 
            c_4*h_g_2*h_l_3 + c_5*h_g*h_l_4 + h_l_5)**0.2
    hh = h_l*1./h_pv
    eta = 1.36603*hh - 0.47719*hh**2 + 0.11116*hh**3
    return h_pv, eta


def calc_alpha(alpha0, alpha1, d):
    """alpha0+alpha1/d
    """
    return alpha0+alpha1/d


def calc_beta(beta0, beta1, d):
    """beta0+beta1/d**4
    """
    return beta0+beta1/d**4


def calc_sigma(sigma0, sigma1, d):
    """sigma0+sigma1*d
    """
    return sigma0+sigma1*d


def calc_sigma_gamma(
        d, sigma0, sigma1, sigma2, gamma0, gamma1, gamma2,
        size_g: float = 0., size_l: float = 0., strain_g: float = 0.,
        strain_l: float = 0.):
    """Calculate H_G (sigma) and H_L (gamma)
    
        H_G**2 = (sigma2+size_g)*d**4 + (sigma1+strain_g)*d**2 + sigma0
        H_L = (gamma2+size_l)*d**2 + (sigma1+strain_l)*d + sigma0
    """
    d_sq = numpy.square(d)
    d_sq_sq = numpy.square(d_sq)

    h_g_sq = (sigma2+size_g) * d_sq_sq + (sigma1+strain_l) * d_sq + sigma0
    h_l = (gamma2+size_l) * d_sq + (gamma1+strain_l) * d + gamma0
    h_g = numpy.sqrt(h_g_sq)
    return h_g, h_l


def calc_y_z_u_v(alpha, beta, sigma, delta_2d):
    """Calculate y, z, u, v

    y = (alpha * sigma**2 + delta)/(sigma * 2**0.5)
    z = (beta * sigma**2 - delta)/(sigma * 2**0.5)
    u = 0.5 * alpha * (alpha*sigma**2 + 2 delta)
    v = 0.5 * beta * (beta*sigma**2 - 2 delta)

    """
    sigma_sq = square(sigma)
    y = (alpha * sigma/(2.**0.5))[:, na] + delta_2d/(sigma*2.**0.5)[:, na]
    z = (beta * sigma/(2.**0.5))[:, na] - delta_2d/(sigma*2.**0.5)[:, na]
    u = (0.5*square(alpha)*sigma_sq)[:, na] + delta_2d*alpha[:, na]
    v = (0.5*square(beta)*sigma_sq)[:, na] - delta_2d*beta[:, na]
    return y, z, u, v

def tof_Carpenter():
    pass


def tof_Jorgensen(alpha, beta, sigma, time, time_hkl):
    norm = 0.5*alpha*beta/(alpha+beta)
    time_2d, time_hkl_2d = numpy.meshgrid(time, time_hkl, indexing="ij")
    delta_2d = time_2d-time_hkl_2d
    y, z, u, v = calc_y_z_u_v(alpha, beta, sigma, delta_2d)

    exp_u = exp(u)
    exp_v = exp(v)
    exp_u[numpy.isinf(exp_u)] = 1e200
    exp_v[numpy.isinf(exp_v)] = 1e200

    res_2d = norm[:, na] * (exp_u * erfc(y) + exp_v * erfc(z))
    return res_2d


def tof_Jorgensen_VonDreele(alpha, beta, sigma, gamma, time, time_hkl):
    two_over_pi = 2.*numpy.pi
    norm = 0.5*alpha*beta/(alpha+beta)
    time_2d, time_hkl_2d = numpy.meshgrid(time, time_hkl, indexing="ij")
    delta_2d = time_2d-time_hkl_2d

    # FIXME: it has to be checked
    # sigma = gamma*(inv_8ln2)**0.5
    h_pv, eta = calc_hpv_eta(sigma, gamma)

    y, z, u, v = calc_y_z_u_v(alpha, beta, sigma, delta_2d) 

    exp_u = exp(u)
    exp_v = exp(v)
    exp_u[numpy.isinf(exp_u)] = 1e200
    exp_v[numpy.isinf(exp_v)] = 1e200

    profile_g_2d = norm[:, na] * (exp_u * erfc(y) + exp_v * erfc(z))

    z1_2d = alpha[:, na]*delta_2d + (1j*0.5*alpha*gamma)[:, na]
    z2_2d = -beta[:, na]*delta_2d + (1j*0.5*beta*gamma)[:, na]

    fz1_2d = exp1(z1_2d)
    fz2_2d = exp1(z2_2d)

    # FIXME: check it
    fz1_2d[numpy.isnan(fz1_2d)] = 0.
    fz1_2d[numpy.isinf(fz1_2d)] = 0.
    fz2_2d[numpy.isnan(fz2_2d)] = 0.
    fz2_2d[numpy.isinf(fz2_2d)] = 0.

    oml_a_2d = -imag(fz1_2d) * two_over_pi
    oml_b_2d = -imag(fz2_2d) * two_over_pi
    profile_l_2d = norm[:, na] * (oml_a_2d + oml_b_2d)
    one_e = 1. - eta
    tof_peak_2d = one_e[:, na] * profile_g_2d + eta[:, na] * profile_l_2d
    return tof_peak_2d


