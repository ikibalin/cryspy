"""Calculation of flip ratio.

Functions
---------
    - calc_flip_ratio
    - calc_intensity_plus_down
    - calc_e_up_loc
    - calc_rotated_ub

"""
import numpy

from cryspy.A_functions_base.function_3_extinction import calc_extinction_2



def calc_flip_ratio(
        sthovl: numpy.ndarray, wavelength: float, field_norm: float,
        u_matrix: numpy.ndarray, polarization: float,
        flipper_efficiency: float, volume_unit_cell: float,
        model_extinction: str, radius: float, mosaicity: float,
        f_nucl: numpy.array, fm_perp_loc: tuple,
        ratio_lambdaover2: float = None, f_nucl_2hkl: numpy.array = None,
        fm_perp_loc_2hkl: tuple = None):
    """Calculate flip ratio."""
    dder = {}
    iint_u, iint_d, dder_iint = calc_intensity_plus_down(
        sthovl, wavelength, field_norm, u_matrix, polarization,
        flipper_efficiency, f_nucl, fm_perp_loc, volume_unit_cell,
        model_extinction, radius, mosaicity)

    if ratio_lambdaover2 is not None:
        sthovl_2hkl = 2.0 * sthovl
        wavelength_2hkl = 0.5*wavelength
        polarization_2hkl = 0.

        iint_u_2hkl, iint_d_2hkl, dder_iint_2hkl = calc_intensity_plus_down(
            sthovl_2hkl, wavelength_2hkl, field_norm, u_matrix,
            polarization_2hkl, flipper_efficiency, f_nucl_2hkl,
            fm_perp_loc_2hkl, volume_unit_cell, model_extinction, radius,
            mosaicity)

        flip_ratio = (iint_u+ratio_lambdaover2*iint_u_2hkl) / \
            (iint_d+ratio_lambdaover2*iint_d_2hkl)
    else:
        flip_ratio = iint_u / iint_d

    return flip_ratio, dder


def calc_intensity_plus_down(
        sthovl: numpy.ndarray, wavelength: float, field_norm: float,
        u_matrix: numpy.ndarray, polarization: float,
        flipper_efficiency: float, f_nucl: numpy.array,
        fm_perp_loc: tuple, volume_unit_cell: float,
        model_extinction: str, radius: float, mosaicity: float):
    """Calculate integrated intensity up and down."""
    dder = {}

    p_u = polarization
    p_d = (2.*flipper_efficiency-1.)*polarization

    phi_d, chi_d, omega_d = 0., 0., 0.
    e_up_loc = calc_e_up_loc(phi_d, chi_d, omega_d, u_matrix)

    mag_p_1, mag_p_2, mag_p_3 = fm_perp_loc

    mag_p_sq = abs(mag_p_1*mag_p_1.conjugate() +
                   mag_p_2*mag_p_2.conjugate() +
                   mag_p_3*mag_p_3.conjugate())

    mag_p_e_u = mag_p_1*e_up_loc[0]+mag_p_2*e_up_loc[1]+mag_p_3*e_up_loc[2]

    f_nucl_sq = abs(f_nucl)**2
    mag_p_e_u_sq = abs(mag_p_e_u*mag_p_e_u.conjugate())
    fnp = (mag_p_e_u*f_nucl.conjugate()+mag_p_e_u.conjugate()*f_nucl).real
    # fp_sq = f_nucl_sq + mag_p_sq + fnp
    # fm_sq = f_nucl_sq + mag_p_sq - fnp
    fp_sq = f_nucl_sq + mag_p_e_u_sq + fnp
    fm_sq = f_nucl_sq + mag_p_e_u_sq - fnp
    fpm_sq = mag_p_sq - mag_p_e_u_sq

    l_model_extinction = ["gauss", "lorentz"]
    if model_extinction.lower() in l_model_extinction:
        yp, dder_yp = calc_extinction_2(
            radius, mosaicity, model_extinction, fp_sq, volume_unit_cell,
            sthovl, wavelength)

        ym, dder_ym = calc_extinction_2(
            radius, mosaicity, model_extinction, fm_sq, volume_unit_cell,
            sthovl, wavelength)

        ypm, dder_ypm = calc_extinction_2(
            radius, mosaicity, model_extinction, fpm_sq, volume_unit_cell,
            sthovl, wavelength)

    else:
        yp = 1. + 0.*f_nucl_sq
        ym = yp
        ypm = yp

    iint_u = 0.5*(1.+p_u)*yp*fp_sq + 0.5*(1.-p_u)*ym*fm_sq + ypm*fpm_sq
    iint_d = 0.5*(1.-p_d)*yp*fp_sq + 0.5*(1.+p_d)*ym*fm_sq + ypm*fpm_sq
    # print("iint_u_1\n", iint_u, "\n")
    # print("iint_d_1\n", iint_d, "\n")

    # pppl = 0.5*((1+p_u)*yp+(1-p_u)*ym)
    # ppmin = 0.5*((1-p_d)*yp+(1+p_d)*ym)
    # pmpl = 0.5*((1+p_u)*yp-(1-p_u)*ym)
    # pmmin = 0.5*((1-p_d)*yp-(1+p_d)*ym)

    # # integral intensities and flipping ratios
    # iint_u = (f_nucl_sq+mag_p_e_u_sq)*pppl + pmpl*fnp + ypm*fpm_sq
    # iint_d = (f_nucl_sq+mag_p_e_u_sq)*ppmin + pmmin*fnp + ypm*fpm_sq
    # print("iint_u_2\n", iint_u, "\n")
    # print("iint_d_2\n", iint_d, "\n")

    return iint_u, iint_d, dder


def calc_e_up_loc(phi_d: float, chi_d: float, omega_d: float,
                  u_matrix: numpy.ndarray):
    """Calculate vertical direction in local Cartezian coordinate system.

    Local Cartezian coordinate system is defined such way that
    (x||a*), (z||c).
    """
    m_u_rot = calc_rotated_ub(phi_d, chi_d, omega_d, u_matrix)
    e_u_vec = numpy.array([0., 0., 1.], dtype=float)
    e_u_loc = numpy.matmul(m_u_rot.transpose(), e_u_vec)
    return e_u_loc


def calc_rotated_ub(phi_d: float, chi_d: float, omega_d: float,
                    u_matrix: numpy.ndarray):
    """Rotated U matrix by three angles."""
    m_phi_d = numpy.array([[numpy.cos(phi_d), numpy.sin(phi_d), 0.],
                           [-numpy.sin(phi_d), numpy.cos(phi_d), 0.],
                           [0., 0., 1.]], dtype=float)

    m_omega_d = numpy.array([[numpy.cos(omega_d), numpy.sin(omega_d), 0.],
                             [-numpy.sin(omega_d), numpy.cos(omega_d), 0.],
                             [0., 0., 1.]], dtype=float)

    m_chi_d = numpy.array([[numpy.cos(chi_d), 0., numpy.sin(chi_d)],
                           [0., 1., 0.],
                           [-numpy.sin(chi_d), 0., numpy.cos(chi_d)]],
                          dtype=float)

    m_u_rot = numpy.matmul(m_omega_d,
                           numpy.matmul(m_chi_d,
                                        numpy.matmul(m_phi_d, u_matrix)))
    return m_u_rot
