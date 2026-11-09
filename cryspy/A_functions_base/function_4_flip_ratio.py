"""Calculation of flip ratio."""
import numpy

from .function_1_algebra import calc_modulus_sq_by_complex_vector,\
    calc_scalar_product_by_complex_vectors

from .function_3_extinction import calc_extinction_2


def calc_f_plus_sq(
        f_nucl: complex, f_mag_perp_vert: complex, flag_f_nucl: bool = False,
        flag_f_mag_perp_vert: bool = False):
    r"""F_{+}^{2} = (F_N + F_{M, \perp, z})^{2}."""
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
            dder["der__f_mag_perp_vert_re"] = \
                der_f_plus_sq__f_mag_perp_vert.real
            dder["der__f_mag_perp_vert_im"] = \
                der_f_plus_sq__f_mag_perp_vert.imag
    return f_plus_sq, dder


def calc_f_minus_sq(
        f_nucl: complex, f_mag_perp_vert: complex, flag_f_nucl: bool = False,
        flag_f_mag_perp_vert: bool = False):
    r"""F_{-}^{2} = (F_N - F_{M, \perp, z})^{2}."""
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
            der_f_minus_sq__f_mag_perp_vert = 2.*f_minus*(
                1.+0.*f_mag_perp_vert)
            dder["der__f_mag_perp_vert_re"] = \
                der_f_minus_sq__f_mag_perp_vert.real
            dder["der__f_mag_perp_vert_im"] = \
                der_f_minus_sq__f_mag_perp_vert.imag
    return f_minus_sq, dder


def calc_flip_ratio(
        f_nucl: float, f_mag_perp: tuple, e_up: tuple,
        beam_polarization: float, flipper_efficiency: float, radius: float,
        mosaicity: float, model_extinction: str, volume_unit_cell: float,
        sthovl: float, wavelength: float,
        flag_f_nucl: bool = False,
        flag_f_mag_perp: bool = False, flag_e_up: bool = False,
        flag_beam_polarization: bool = False,
        flag_flipper_efficiency: bool = False, flag_radius: bool = False,
        flag_mosaicity: bool = False, flag_volume_unit_cell: bool = False,
        flag_sthovl: bool = False, flag_volume_wavelength: bool = False):
    """Calculate flip ratio.

    Calculate flip ratio by given nuclear structure factor and
    perpendicular component of magnetic structure factor
    (in 10**-12 cm)

    cell (only for extinction correction)
    f_nucl [hkl]
    f_m_perp [hkl]
    delta_f_nucl [hkl, parameters]
    delta_f_m_perp [hkl, parameters]
    """
    flag_derivatives = any([
        flag_f_nucl, flag_f_mag_perp, flag_e_up, flag_beam_polarization,
        flag_flipper_efficiency, flag_radius, flag_mosaicity,
        flag_volume_unit_cell, flag_sthovl, flag_volume_wavelength])

    np_b_p = numpy.array(beam_polarization, dtype=float)
    np_f_e = numpy.array(flipper_efficiency, dtype=float)
    p_u = np_b_p
    p_d = (2. * np_f_e - 1.) * p_u

    part_u_u, part_u_d = 0.5*(1. + p_u), 0.5*(1. - p_u)
    part_d_u, part_d_d = 0.5*(1. - p_d), 0.5*(1. + p_d)

    # f_mag_perp_vert = (f_mag_perp, e_up), e_up is a real unit vector
    f_mag_perp_vert = calc_scalar_product_by_complex_vectors(
        f_mag_perp, e_up)[0]
    flag_f_mag_perp_vert = any([flag_f_mag_perp, flag_e_up])

    f_plus_sq, dder_f_plus_sq = calc_f_plus_sq(
        f_nucl, f_mag_perp_vert, flag_f_nucl, flag_f_mag_perp_vert)
    f_minus_sq, dder_f_minus_sq = calc_f_minus_sq(
        f_nucl, f_mag_perp_vert, flag_f_nucl, flag_f_mag_perp_vert)
    f_mag_perp_sq, dder_f_mag_perp_sq = calc_modulus_sq_by_complex_vector(
        f_mag_perp, flag_f_mag_perp)
    f_mag_perp_vert_sq = (f_mag_perp_vert*f_mag_perp_vert.conjugate()).real
    f_mag_sq = f_mag_perp_sq - f_mag_perp_vert_sq

    # extinction correction
    y_plus, dder_y_plus = calc_extinction_2(
        radius, mosaicity, model_extinction, f_plus_sq, volume_unit_cell,
        sthovl, wavelength)
    y_minus, dder_y_minus = calc_extinction_2(
        radius, mosaicity, model_extinction, f_minus_sq, volume_unit_cell,
        sthovl, wavelength)
    y_mag, dder_y_mag = calc_extinction_2(
        radius, mosaicity, model_extinction, f_mag_sq, volume_unit_cell,
        sthovl, wavelength)

    int_up = part_u_u*y_plus*f_plus_sq + part_u_d*y_minus*f_minus_sq + \
        y_mag*f_mag_sq
    int_down = part_d_u*y_plus*f_plus_sq + part_d_d*y_minus*f_minus_sq + \
        y_mag*f_mag_sq

    flip_ratio = int_up / int_down

    dder = {}
    if not flag_derivatives:
        return flip_ratio, dder

    # d(flip_ratio) = (d(int_up) - flip_ratio d(int_down)) / int_down
    def der_flip_ratio(der_int_up, der_int_down):
        return (der_int_up - flip_ratio * der_int_down) / int_down

    # derivatives of y*F^2 over F^2 (extinction depends on F^2)
    der_yf_plus = y_plus + f_plus_sq*dder_y_plus["der_yext__f_sq"]
    der_yf_minus = y_minus + f_minus_sq*dder_y_minus["der_yext__f_sq"]
    der_yf_mag = y_mag + f_mag_sq*dder_y_mag["der_yext__f_sq"]

    def der_flip_ratio_by_f_sq(der_f_plus_sq, der_f_minus_sq, der_f_mag_sq):
        der_int_up = part_u_u*der_yf_plus*der_f_plus_sq + \
            part_u_d*der_yf_minus*der_f_minus_sq + der_yf_mag*der_f_mag_sq
        der_int_down = part_d_u*der_yf_plus*der_f_plus_sq + \
            part_d_d*der_yf_minus*der_f_minus_sq + der_yf_mag*der_f_mag_sq
        return der_flip_ratio(der_int_up, der_int_down)

    if flag_beam_polarization:
        tt = y_plus*f_plus_sq - y_minus*f_minus_sq
        dder["der_flip_ratio__beam_polarization"] = der_flip_ratio(
            0.5*tt, -0.5*(2.*np_f_e - 1.)*tt)

    if flag_flipper_efficiency:
        tt = y_plus*f_plus_sq - y_minus*f_minus_sq
        dder["der_flip_ratio__flipper_efficiency"] = der_flip_ratio(
            0.*tt, -np_b_p*tt)

    if flag_f_nucl:
        # F_M,perp,vert does not depend on F_N
        zeros = numpy.zeros_like(f_plus_sq)
        dder["der_flip_ratio__f_nucl_re"] = der_flip_ratio_by_f_sq(
            dder_f_plus_sq["der__f_nucl_re"], dder_f_minus_sq["der__f_nucl_re"],
            zeros)
        dder["der_flip_ratio__f_nucl_im"] = der_flip_ratio_by_f_sq(
            dder_f_plus_sq["der__f_nucl_im"], dder_f_minus_sq["der__f_nucl_im"],
            zeros)

    f_plus = f_nucl + f_mag_perp_vert
    f_minus = f_nucl - f_mag_perp_vert

    if flag_f_mag_perp:
        # d(F_M,perp,vert)/d(Re F_M,perp,j) = e_j, d/d(Im ...) = i e_j
        l_der_re, l_der_im = [], []
        for e_j, der_re_j, der_im_j in zip(
                e_up, dder_f_mag_perp_sq["der__vector_1_re"],
                dder_f_mag_perp_sq["der__vector_1_im"]):
            l_der_re.append(der_flip_ratio_by_f_sq(
                2.*e_j*f_plus.real, -2.*e_j*f_minus.real,
                der_re_j - 2.*e_j*f_mag_perp_vert.real))
            l_der_im.append(der_flip_ratio_by_f_sq(
                2.*e_j*f_plus.imag, -2.*e_j*f_minus.imag,
                der_im_j - 2.*e_j*f_mag_perp_vert.imag))
        dder["der_flip_ratio__f_mag_perp_re"] = tuple(l_der_re)
        dder["der_flip_ratio__f_mag_perp_im"] = tuple(l_der_im)

    if flag_e_up:
        # d(F_M,perp,vert)/d(e_j) = F_M,perp,j
        l_der = []
        for f_mag_perp_j in f_mag_perp:
            l_der.append(der_flip_ratio_by_f_sq(
                2.*(f_plus.conjugate()*f_mag_perp_j).real,
                -2.*(f_minus.conjugate()*f_mag_perp_j).real,
                -2.*(f_mag_perp_vert.conjugate()*f_mag_perp_j).real))
        dder["der_flip_ratio__e_up"] = tuple(l_der)

    return flip_ratio, dder
