import numpy
import scipy
import scipy.special

na = numpy.newaxis



def calc_lorentz_factor(ttheta, flag_ttheta: bool = False):
    """Angular (bank-angle) part of the TOF Lorentz factor.

    Returns ``sin(ttheta)`` for the fixed detector bank angle ``ttheta``
    (the full scattering angle ``2theta_bank``). The reflection-dependent
    ``d**4`` part is applied separately, so the complete TOF Lorentz
    factor is ``d**4 * sin(ttheta)`` (cf. CrysFML ``lorentz = dsp4 * sina``).
    """
    res = numpy.sin(ttheta)
    dder = {}
    if flag_ttheta:
        dder["ttheta"] = numpy.cos(ttheta)
    return res, dder


def calc_spectrum(time, spectrum_type, spectrum_parameters, flag_spectrum_parameters=False):
    exp = numpy.exp
    time_sq = numpy.square(time)
    time_4 = numpy.square(time_sq)
    a0 = spectrum_parameters[0]
    a1 = spectrum_parameters[1]
    a2 = spectrum_parameters[2]
    a3 = spectrum_parameters[3]
    a4 = spectrum_parameters[4]
    a5 = spectrum_parameters[5]
    a6 = spectrum_parameters[6]
    a7 = spectrum_parameters[7]
    a8 = spectrum_parameters[8]
    if spectrum_type == "Empirical-Exponents":
        res = a0 + a1 * exp(-a2 * time) + \
            a3 * exp(-a4 * time_sq) + \
            a5 * exp(-a6 * time*time_sq) + \
            a7 * exp(-a8 * time_4)
    else:  # spectrum_type == "Maxwell"
        res = a0 + a1 * exp(-a2 * time_sq)/(time*time_4) + \
            a3 * exp(-a4 * time_sq) + \
            a5 * exp(-a6 * time*time_sq) + \
            a7 * exp(-a8 * numpy.square(time_sq))
    dder = {}
    if flag_spectrum_parameters:
        pass
    return res, dder


def calc_time_for_epithermal_neutrons_by_d(d, zero, dtt1, zerot, dtt1t, dtt2t, width, x_cross):
    time_e = zero + dtt1 * d
    time_t = zerot + dtt1t * d - dtt2t / d
    n_cross = 0.5*scipy.special.erfc(width * (x_cross - 1./d))
    time = n_cross * time_e + (1-n_cross) * time_t
    return time


def calc_time_for_thermal_neutrons_by_d(d, zero, dtt1, dtt2):
    time = zero + dtt1 * d + dtt2 * numpy.square(d)
    return time


def calc_d_by_time_for_thermal_neutrons(time, zero, dtt1, dtt2):
    det = numpy.sqrt(numpy.square(dtt1) - 4.*(zero-time)*dtt2)
    if dtt2 < 0.:
        d = (det-dtt1)/(2.*dtt2)
    elif dtt2 > 0.:
        d = (-det-dtt1)/(2.*dtt2)
    else:
        d = (time - zero)/dtt1
    return d


def calc_d_min_max_by_time_thermal_neutrons(time, zero, dtt1, dtt2):
    zero, dtt1, dtt2 = zero.squeeze(), dtt1.squeeze(), dtt2.squeeze()
    time_min = numpy.min(time)
    time_max = numpy.max(time)
    d_min = calc_d_by_time_for_thermal_neutrons(time_min, zero, dtt1, dtt2)
    d_max = calc_d_by_time_for_thermal_neutrons(time_max, zero, dtt1, dtt2)
    return numpy.stack([d_min, d_max], axis=0)

def calc_d_min_max_by_time_epithermal_neutrons(time, zero, dtt1, zerot, dtt1t, dtt2t):
    zero, dtt1, zerot, dtt1t, dtt2t = zero.squeeze(), dtt1.squeeze(), zerot.squeeze(), dtt1t.squeeze(), dtt2t.squeeze()
    d_min_max_t = calc_d_min_max_by_time_thermal_neutrons(time, zerot, dtt1t, dtt2t)
    time_min = numpy.min(time)
    time_max = numpy.max(time)
    d_min_max_e = numpy.stack([(time_min-zero)/dtt1, (time_max-zero)/dtt1], axis=0)
    d_min_max = numpy.array([
        min([d_min_max_t[0], d_min_max_e[0]]),
        max([d_min_max_t[1], d_min_max_e[1]])], dtype=float)
    return d_min_max


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


def tof_Jorgensen(alpha, beta, sigma, time, time_hkl, cutoff_fwhm=0.):
    norm = 0.5*alpha*beta/(alpha+beta)
    time_2d, time_hkl_2d = numpy.meshgrid(time, time_hkl, indexing="ij")
    delta_2d = time_2d-time_hkl_2d
    y, z, u, v = calc_y_z_u_v(alpha, beta, sigma, delta_2d)

    with numpy.errstate(over='ignore'):
        exp_u = exp(u)
        exp_v = exp(v)
    exp_u[numpy.isinf(exp_u)] = 1e200
    exp_v[numpy.isinf(exp_v)] = 1e200

    profile_g_2d = norm[:, na] * (exp_u * erfc(y) + exp_v * erfc(z))
    if cutoff_fwhm > 0.:
        half_width = cutoff_fwhm * (sigma*numpy.sqrt(8.*numpy.log(2.)) + 1./alpha + 1./beta)
        profile_g_2d = profile_g_2d * (numpy.abs(delta_2d) <= half_width[:, na])
    return profile_g_2d


def tof_Jorgensen_VonDreele(alpha, beta, sigma, gamma, time, time_hkl, cutoff_fwhm=0.):
    two_over_pi = 2./numpy.pi
    norm = 0.5*alpha*beta/(alpha+beta)
    time_2d, time_hkl_2d = numpy.meshgrid(time, time_hkl, indexing="ij")
    delta_2d = time_2d-time_hkl_2d

    # Match FullProf/CrysFML: build one Thompson-Cox-Hastings pseudo-Voigt
    # FWHM from the Gaussian FWHM and Lorentzian FWHM, then use that common
    # width for both the Gaussian and Lorentzian components.
    h_g_fwhm = sigma*numpy.sqrt(8.*numpy.log(2.))
    h_pv, eta = calc_hpv_eta(h_g_fwhm, gamma)
    sigma_c = h_pv*numpy.sqrt(inv_8ln2)
    gamma_c = h_pv

    y, z, u, v = calc_y_z_u_v(alpha, beta, sigma_c, delta_2d)

    with numpy.errstate(over='ignore'):
        exp_u = exp(u)
        exp_v = exp(v)
    exp_u[numpy.isinf(exp_u)] = 1e200
    exp_v[numpy.isinf(exp_v)] = 1e200

    profile_g_2d = norm[:, na] * (exp_u * erfc(y) + exp_v * erfc(z))

    z1_2d = alpha[:, na]*delta_2d + (1j*0.5*alpha*gamma_c)[:, na]
    z2_2d = -beta[:, na]*delta_2d + (1j*0.5*beta*gamma_c)[:, na]

    # The Lorentzian term is exp(z) * E1(z); omitting exp(z) breaks
    # the pseudo-Voigt tails for non-zero gamma.
    with numpy.errstate(over='ignore', invalid='ignore'):
        fz1_2d = exp(z1_2d) * exp1(z1_2d)
        fz2_2d = exp(z2_2d) * exp1(z2_2d)

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
    if cutoff_fwhm > 0.:
        half_width = cutoff_fwhm * (h_pv + 1./alpha + 1./beta)
        tof_peak_2d = tof_peak_2d * (numpy.abs(delta_2d) <= half_width[:, na])
    return tof_peak_2d


def tof_non_convoluted_pseudo_voigt(sigma, gamma, time, time_hkl, cutoff_fwhm=0.):
    """Direct symmetric pseudo-Voigt profile in TOF coordinates."""
    h_g = numpy.sqrt(8.*numpy.log(2.))*sigma
    h_pv, eta = calc_hpv_eta(h_g, gamma)
    time_2d, time_hkl_2d = numpy.meshgrid(time, time_hkl, indexing="ij")
    delta_2d = time_2d-time_hkl_2d
    delta_over_h_pv = delta_2d/h_pv[:, na]

    profile_g_2d = (
        2./h_pv*numpy.sqrt(numpy.log(2.)/numpy.pi))[:, na] * numpy.exp(
            -4.*numpy.log(2.)*numpy.square(delta_over_h_pv))
    profile_l_2d = (
        2./(numpy.pi*h_pv))[:, na] / (
            1. + 4.*numpy.square(delta_over_h_pv))

    res_2d = eta[:, na] * profile_l_2d + (1.-eta)[:, na] * profile_g_2d
    if cutoff_fwhm > 0.:
        res_2d = res_2d * (numpy.abs(delta_over_h_pv) <= cutoff_fwhm)
    return res_2d



inv_8ln2 = 0.18033688011112042591999058512524

na = numpy.newaxis
square = numpy.square
exp = numpy.exp
imag = numpy.imag
erfc = scipy.special.erfc
exp1 = scipy.special.exp1  # E1 = \int_{z}^{\inf} exp(-t)/t dt



def calc_alpha(alpha0, alpha1, d):
    """alpha0+alpha1/d
    """
    return alpha0+alpha1/d


def calc_beta(beta0, beta1, d):
    """beta0+beta1/d**4
    """
    return beta0+beta1*numpy.power(d, -4)


def calc_sigma(d, sigma0, sigma1, sigma2, size_g:float=0., strain_g:float=0.):
    """
    H_G**2 = (sigma2+size_g)*d**4 + (sigma1+strain_g)*d**2 + sigma0
    """
    d_sq = numpy.square(d)
    d_sq_sq = numpy.square(d_sq)
    h_g_sq = numpy.abs((sigma2+size_g) * d_sq_sq + (sigma1+strain_g) * d_sq + sigma0)
    h_g = numpy.sqrt(h_g_sq)
    return h_g


def calc_sigma_gamma(
        d, sigma0, sigma1, sigma2, gamma0, gamma1, gamma2,
        size_g: float = 0., size_l: float = 0., strain_g: float = 0.,
        strain_l: float = 0.):
    """Calculate H_G (sigma) and H_L (gamma)

        H_G**2 = (sigma2+size_g)*d**4 + (sigma1+strain_g)*d**2 + sigma0
        H_L = (gamma2+size_l)*d**2 + (sigma1+strain_l)*d + sigma0
    """
    d_sq = numpy.square(d)
    h_g = calc_sigma(d, sigma0, sigma1, sigma2, size_g, strain_g)
    h_l = (gamma2+size_l) * d_sq + (gamma1+strain_l) * d + gamma0
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



def calc_peak_shape_function(alphas, betas, sigmas,
        d, time, time_hkl, gammas = None, size_g: float = 0., strain_g: float = 0.,
        size_l: float = 0., strain_l: float = 0., peak_shape: str = "pseudo-Voigt", cutoff_fwhm: float = 0.):
    """Calculate peak-shape function F(DELTA)
    For Gauss peak-shape:
        F(DELTA) = alpha * beta / (alpha+beta) * [exp(u) erfc(y) + exp(v) erfc(z)]
    For pseudo-Voigt peak-shape:
        F(DELTA) = ...
    alpha = alpha0 + alpha1 / d
    beta = beta0 + beta1 / d**4
    sigma = sigma0 + sigma1 * d
    u = alpha/2 * (alpha * sigma**2 + 2 DELTA)
    v = beta/2 * (beta * sigma**2 - 2 DELTA)
    y = (alpha * sigma**2 + DELTA)/(sigma * 2**0.5)
    z = (beta * sigma**2 - DELTA)/(sigma * 2**0.5)
    DELTA = time - time_hkl
    """
    sigma0, sigma1, sigma2 = sigmas[0], sigmas[1], sigmas[2]

    if peak_shape == "Gauss" or peak_shape == "pseudo-Voigt":
        alpha0, alpha1 = alphas[0], alphas[1]
        beta0, beta1 = betas[0], betas[1]
        alpha = alpha0 + alpha1 / d
        beta = beta0 + beta1 / d**4

        if peak_shape == "Gauss":
            sigma = calc_sigma(
                d, sigma0, sigma1, sigma2, size_g=size_g, strain_g=strain_g)

            res_2d = tof_Jorgensen(alpha, beta, sigma, time, time_hkl, cutoff_fwhm=cutoff_fwhm)

        elif peak_shape == "pseudo-Voigt":
            gamma0, gamma1, gamma2 = gammas[0], gammas[1], gammas[2]

            sigma, gamma = calc_sigma_gamma(
                d, sigma0, sigma1, sigma2, gamma0,
                gamma1, gamma2, size_g=size_g, size_l=size_l,
                strain_g=strain_g, strain_l=strain_l)

            res_2d = tof_Jorgensen_VonDreele(
                alpha, beta, sigma, gamma, time, time_hkl, cutoff_fwhm=cutoff_fwhm)

    elif peak_shape == "non-conv-pseudo-Voigt":
        gamma0, gamma1, gamma2 = gammas[0], gammas[1], gammas[2]

        sigma, gamma = calc_sigma_gamma(
            d, sigma0, sigma1, sigma2, gamma0, gamma1, gamma2,
            size_g=0., size_l=0., strain_g=0., strain_l=0.)

        res_2d = tof_non_convoluted_pseudo_voigt(
            sigma, gamma, time, time_hkl, cutoff_fwhm=cutoff_fwhm)

    else:
        raise ValueError(f"Unknown peak shape: {peak_shape}")

    return res_2d
