import numpy

na = numpy.newaxis


def calc_lorentz_factor(ttheta, k:float=0.0, cthm:float = 0.91, flag_ttheta: bool=False):
    """Lorentz factor for 1D powder diffraction profile.

    k is 0 for neutrons
    k is 0.5 for characteristic x-ray radiation
    k is 0.1 for synchrotron radiation
    cthm = cos**2 (2 theta_M )
    """
    hh = 1. - k  + k * cthm * numpy.square(numpy.cos(ttheta))
    res = hh/(numpy.sin(0.5*ttheta)*numpy.sin(ttheta))
    dder = {}
    if flag_ttheta:
        dder["ttheta"] = - numpy.square(res) * \
            (0.5*numpy.cos(0.5*ttheta)*numpy.sin(ttheta) +
             numpy.sin(0.5*ttheta)*numpy.cos(ttheta)) 
    return res, dder


def func_asymmetry_f_a(z, flag_z: bool = False):
    """Function F_a(z) for asymmetry factor.
    """
    res = 2. * z * numpy.exp(-numpy.square(z))
    dder = {}
    if flag_z:
        dder["z"] = 2.*(1.-2.*numpy.square(z)) * numpy.exp(-numpy.square(z))
    return res, dder


def func_asymmetry_f_b(z, flag_z: bool = False):
    """Function F_b(z) for asymmetry factor.
    """
    f_a , dder_f_a = func_asymmetry_f_a(z, flag_z=flag_z)
    res = 2*(2*numpy.square(z)-3)*f_a
    dder = {}
    if flag_z:
        dder["z"] = 8 * z * f_a + 2*(2*numpy.square(z)-3)*dder_f_a["z"]
    return res, dder


def calc_asymmetry_factor(z, ttheta, p_1, p_2, p_3, p_4, 
        flag_z: bool = False, flag_p_1: bool = False, flag_p_2: bool = False,
        flag_p_3: bool = False, flag_p_4: bool = False):
    """Asymmetry factor for 1D powder diffraction profile.
    """
    f_a , dder_f_a = func_asymmetry_f_a(z, flag_z=flag_z)
    f_b , dder_f_b = func_asymmetry_f_b(z, flag_z=flag_z)
    
    inv_tan_th = numpy.expand_dims(1./numpy.tan(0.5*ttheta), axis=-1)
    inv_tan_tth = numpy.expand_dims(1./numpy.tan(ttheta), axis=-1)
    
    res = 1. + (p_1 * f_a + p_2 * f_b)*inv_tan_th \
         + (p_3 * f_a + p_4 * f_b)*inv_tan_tth
    dder = {}
    if flag_z:
        dder["z"] = (p_1 * dder_f_a["z"] + p_2 * dder_f_b["z"])*inv_tan_th + \
            (p_3 * dder_f_a["z"] + p_4 * dder_f_b["z"])*inv_tan_tth
    if flag_p_1:
        dder["p_1"] = f_a * inv_tan_th
    if flag_p_2:
        dder["p_2"] = f_b * inv_tan_th
    if flag_p_3:
        dder["p_3"] = f_a * inv_tan_tth
    if flag_p_4:
        dder["p_4"] = f_b * inv_tan_tth

    return res, dder


def func_gauss_by_h_pv(z, h_pv, flag_z: bool = False, flag_h_pv: bool = False):
    """Gauss function as function of h_pv
    """
    inv_h_pv = 1./h_pv 
    inv_h_pv_sq = numpy.square(inv_h_pv)
    z_deg = z * 180./numpy.pi
    z_deg_sq = numpy.square(z_deg)
    a_g = 2.*inv_h_pv*numpy.sqrt(numpy.log(2)/numpy.pi)
    b_g = 4.*numpy.log(2)*inv_h_pv_sq
    res = numpy.expand_dims(a_g, axis=-1) * numpy.exp(-numpy.expand_dims(b_g, axis=-1) * z_deg_sq) # numpy.where(val_1 < 5., numpy.exp(-val_1), 0.)
    dder = {}
    if  flag_z:
        dder["z"] = -2.*z_deg*res * 180./numpy.pi
    if flag_h_pv:
        dder["h_pv"] = (-1. + 2.* numpy.expand_dims(b_g, axis=-1)*z_deg_sq)*res*numpy.expand_dims(inv_h_pv, axis=-1)
    return res, dder


def func_lorentz_by_h_pv(z, h_pv, flag_z: bool = False, flag_h_pv: bool = False):
    """Gauss function as function of h_pv
    """
    inv_h_pv = 1./h_pv
    inv_h_pv_sq = numpy.square(inv_h_pv)
    z_deg = z * 180./numpy.pi
    c_a = 2./numpy.pi
    a_l = c_a * inv_h_pv
    b_l = 4.*inv_h_pv_sq
    z_deg_sq = numpy.square(z_deg)
    res = numpy.expand_dims(a_l, axis=-1) /(1+ numpy.expand_dims(b_l, axis=-1) * z_deg_sq)
    dder = {}
    if  flag_z:
        dder["z"] = -2.*z_deg*numpy.expand_dims(b_l,axis=-1)*res/(1.+numpy.expand_dims(b_l, axis=-1)*z_deg_sq) * 180./numpy.pi
    if flag_h_pv:
        dder["h_pv"] = (c_a * (numpy.expand_dims(h_pv, axis=-1) + 4*z_deg_sq) - \
            c_a * numpy.expand_dims(h_pv, axis=-1))/numpy.square(numpy.expand_dims(h_pv, axis=-1) + 4*z_deg_sq)
    return res, dder


def calc_h_g(u, v, w, i_g, ttheta,
        flag_u: bool = False, flag_v: bool = False, flag_w: bool = False, flag_i_g: bool = False):
    """Calculate H_G. 
    For more documentation see documentation module "Powder diffraction at constant wavelenght".
    """
    tan_th  = numpy.tan(0.5*ttheta)
    tan_th_sq = numpy.square(tan_th)
    inv_c_th_sq = 1./numpy.square(numpy.cos(0.5*ttheta))
    hh = u*tan_th_sq + v*tan_th + w + i_g * inv_c_th_sq
    hh[ hh<numpy.finfo(float).eps ] = numpy.finfo(float).eps # protection against negative variables
    res = numpy.sqrt(hh)
    dder = {}
    if flag_u:
        dder["u"] = 0.5*tan_th_sq/res
    if flag_v:
        dder["v"] = 0.5*tan_th/res
    if flag_w:
        dder["w"] = 0.5/res
    if flag_i_g:
        dder["i_g"] = 0.5*inv_c_th_sq/res
    return res, dder


def calc_h_l(x, y, ttheta, flag_x: bool = False, flag_y: bool = False):
    """Calculate H_L. 
    For more documentation see documentation module "Powder diffraction at constant wavelenght".
    """
    tan_th  = numpy.tan(0.5*ttheta)
    inv_c_th = 1./numpy.cos(0.5*ttheta)
    res = x * tan_th + y*inv_c_th
    dder = {}
    if flag_x:
        dder["x"] = tan_th
    if flag_y:
        dder["y"] = inv_c_th
    return res, dder


def calc_h_pv(h_g, h_l, flag_h_g: bool = False, flag_h_l: bool = False):
    """Calculate H_pV. 
    For more documentation see documentation module "Powder diffraction at constant wavelenght".
    """
    h_g_2, h_l_2 = numpy.square(h_g), numpy.square(h_l)
    h_g_3, h_l_3 = h_g_2*h_g, h_l_2*h_l
    h_g_4, h_l_4 = numpy.square(h_g_2), numpy.square(h_l_2)
    h_g_5, h_l_5 = h_g_4*h_g, h_l_4*h_l

    c_2, c_3, c_4, c_5 = 2.69269, 2.42843, 4.47163, 0.07842
    hh = h_g_5 + c_2*h_g_4*h_l + c_3*h_g_3*h_l_2 + c_4*h_g_2*h_l_3 + c_5*h_g*h_l_4 + h_l_5
    res = numpy.power(hh, 0.2)

    dder = {}
    if flag_h_g or flag_h_l:
        c_help = -0.2*numpy.power(hh, -0.8)
    if flag_h_g:
        dder["h_g"] = c_help * (5*h_g_4 + 4*c_2*h_g_3*h_l + 3*c_3*h_g_2*h_l_2 + 2*c_4*h_g*h_l_3 + c_5*h_l_4)
    if flag_h_l:
        dder["h_l"] = c_help * (c_2*h_g_4 + 2*c_3*h_g_3*h_l + 3*c_4*h_g_2*h_l_2 + 4*c_5*h_g*h_l_3 + 5*h_l_4)
    return res, dder


def calc_eta(h_l, h_pv, flag_h_l: bool = False, flag_h_pv: bool = False):
    """Calculate H_pV. 
    For more documentation see documentation module "Powder diffraction at constant wavelenght".
    """
    c_1 = 1.36603
    c_2 = -0.47719
    c_3 = 0.11116

    rel = h_l/h_pv
    rel_sq = numpy.square(rel)
    res = c_1*rel + c_2*rel_sq  + c_3*rel_sq*rel

    dder = {}
    if flag_h_pv or flag_h_l:
        c_help = c_1 + 2*c_2*rel  + 3*c_3*rel_sq
    if flag_h_pv:
        dder["h_pv"] = -1.*c_help * rel / h_pv 
    if flag_h_l:
        dder["h_l"] = c_help/h_pv
    return res, dder


def calc_profile_pseudo_voight(ttheta, ttheta_hkl, u, v, w, i_g, x, y,
        p_1, p_2, p_3, p_4, 
        flag_ttheta: bool=False, 
        flag_ttheta_hkl: bool=False, flag_u: bool=False,
        flag_v: bool=False, flag_w: bool=False, flag_i_g: bool=False, flag_x: bool=False,
        flag_y: bool=False, flag_p_1: bool = False, flag_p_2: bool = False,
        flag_p_3: bool = False, flag_p_4: bool = False):
    """Calculate profile as psevdo-Voight function defined by parameters.
    For more documentation see documentation module "Powder diffraction at constant wavelenght".
    """
    delta_angle = (ttheta[:, na] - ttheta_hkl[na, :])
    flag_delta_angle = flag_ttheta_hkl or flag_ttheta

    flag_h_g = flag_u or flag_v or flag_w or flag_i_g
    h_g, dder_h_g = calc_h_g(
        u,v, w, i_g, ttheta,
        flag_u=flag_u, flag_v=flag_v, flag_w=flag_w, flag_i_g=flag_i_g)

    flag_h_l = flag_x or flag_y
    h_l, dder_h_l = calc_h_l(
        x, y, ttheta, flag_x=flag_x, flag_y=flag_y)

    flag_h_pv = flag_h_g or flag_h_l
    h_pv, dder_h_pv = calc_h_pv(h_g, h_l, flag_h_g=flag_h_g, flag_h_l=flag_h_l)

    eta, dder_eta = calc_eta(h_l, h_pv, flag_h_l=flag_h_l, flag_h_pv=flag_h_pv)

    z = (delta_angle*180./numpy.pi)/numpy.expand_dims(h_pv, axis=1)
    flag_z = flag_h_pv or flag_ttheta_hkl or flag_ttheta
    af, dder_af  = calc_asymmetry_factor(z, ttheta, p_1, p_2, p_3, p_4, 
        flag_z = flag_z, flag_p_1 =  flag_p_1, flag_p_2 =  flag_p_2,
        flag_p_3 = flag_p_3, flag_p_4 = flag_p_4)
    
    profile_l, dder_profile_l = func_lorentz_by_h_pv(
        delta_angle, h_pv, flag_z=flag_delta_angle, flag_h_pv=flag_h_pv)
    profile_g, dder_profile_l = func_gauss_by_h_pv(
        delta_angle, h_pv, flag_z=flag_delta_angle, flag_h_pv=flag_h_pv)

    res = (numpy.expand_dims(eta, axis=1) * profile_l + numpy.expand_dims((1.-eta), axis=1)*profile_g)*af
    dder = {}
    return res, dder


def calc_gamma_nu_by_ttheta_phi(ttheta, phi, flag_ttheta: bool = False, flag_phi: bool = False):
    """See the documentation module "Powder diffraction at constant wavelength".
    """
    s_tth, c_tth = numpy.sin(ttheta), numpy.cos(ttheta)
    s_phi, c_phi = numpy.sin(phi), numpy.cos(phi)

    gamma = numpy.arctan2(s_tth*c_phi, c_tth)
    nu = numpy.arcsin(s_tth*s_phi)

    dder_gamma = {}
    dder_nu = {}
    return gamma, nu, dder_gamma, dder_nu


def calc_ttheta_phi_by_gamma_nu(gamma, nu, flag_gamma: bool = False, flag_nu: bool = False):
    """See the documentation module "Powder diffraction at constant wavelength".
    """
    s_ga, c_ga = numpy.sin(gamma), numpy.cos(gamma)
    s_nu, c_nu = numpy.sin(nu), numpy.cos(nu)

    ttheta = numpy.arccos(c_nu*c_ga)
    phi = numpy.arctan2(s_nu, s_ga*c_nu)

    dder_ttheta = {}
    dder_phi = {}
    return ttheta, phi, dder_ttheta, dder_phi


def calc_profile_pseudo_voight_2d(ttheta, phi,
        ttheta_hkl, u, v, w, i_g, x, y,
        p_1, p_2, p_3, p_4, p_phi,
        flag_ttheta: bool=False, 
        flag_ttheta_hkl: bool=False, flag_u: bool=False,
        flag_v: bool=False, flag_w: bool=False, flag_i_g: bool=False, flag_x: bool=False,
        flag_y: bool=False, flag_p_1: bool = False, flag_p_2: bool = False,
        flag_p_3: bool = False, flag_p_4: bool = False, flag_p_phi: bool = False):
    """Calculate profile as psevdo-Voight function defined by parameters.
    For more documentation see documentation module "Powder diffraction at constant wavelenght".
    """

    delta_ttheta = ttheta[:, :, na] - ttheta_hkl[na, na, :]
    flag_delta_angle = flag_ttheta_hkl or flag_ttheta

    flag_h_g = flag_u or flag_v or flag_w or flag_i_g
    h_g, dder_h_g = calc_h_g(u,v, w, i_g, ttheta, flag_u=flag_u, flag_v=flag_v, flag_w=flag_w, flag_i_g=flag_i_g)

    flag_h_l = flag_x or flag_y
    h_l, dder_h_l = calc_h_l(x, y, ttheta, flag_x=flag_x, flag_y=flag_y)

    flag_h_pv = flag_h_g or flag_h_l
    h_pv, dder_h_pv = calc_h_pv(h_g, h_l, flag_h_g=flag_h_g, flag_h_l=flag_h_l)

    eta, dder_eta = calc_eta(h_l, h_pv, flag_h_l=flag_h_l, flag_h_pv=flag_h_pv)


    coeff_phi = calc_coeff_phi(phi, p_phi, flag_p_phi=flag_p_phi)[0]
    h_pv_phi = h_pv*coeff_phi
    flag_h_pv_phi = flag_h_pv or flag_p_phi

    z = (delta_ttheta*180./numpy.pi)/numpy.expand_dims(h_pv_phi, axis=-1)
    flag_z = flag_h_pv_phi or flag_ttheta_hkl or flag_ttheta
    af, dder_af  = calc_asymmetry_factor(z, ttheta, p_1, p_2, p_3, p_4, 
        flag_z = flag_z, flag_p_1 =  flag_p_1, flag_p_2 =  flag_p_2,
        flag_p_3 = flag_p_3, flag_p_4 = flag_p_4)
    
    profile_l, dder_profile_l = func_lorentz_by_h_pv(delta_ttheta, h_pv_phi, flag_z=flag_delta_angle, flag_h_pv=flag_h_pv_phi)
    profile_g, dder_profile_l = func_gauss_by_h_pv(delta_ttheta, h_pv_phi, flag_z=flag_delta_angle, flag_h_pv=flag_h_pv_phi)

    res = (numpy.expand_dims(eta, axis=-1) * profile_l + numpy.expand_dims((1.-eta), axis=-1)*profile_g)*af
    dder = {}
    return res, dder


def calc_coeff_phi(phi, p_phi, flag_p_phi: bool=False):
    res= 0.5*(p_phi+1.)-numpy.cos(phi)*(p_phi-1.)*0.5
    dder = {}
    if flag_p_phi:
        dder["p_phi"] = 0.5-0.5*numpy.cos(phi)
    return res, dder