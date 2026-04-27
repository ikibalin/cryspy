"""Functions for magnetic section.

Functions
---------
    - get_j0_j2_by_symbol
"""
import numpy

from cryspy.A_functions_base.database import DATABASE
from cryspy.A_functions_base.charge_form_factor import get_atomic_symbol_ion_charge_isotope_number_by_ion_symbol

from cryspy.A_functions_base.matrix_operations import calc_m_v, calc_inv_m, calc_m1_m2_inv_m1
from cryspy.A_functions_base.orbital_functions import calc_jl


from cryspy.A_functions_base.unit_cell import calc_m_b_by_unit_cell_parameters, calc_m_m_by_unit_cell_parameters, calc_q_ccs_by_unit_cell_parameters


from cryspy.A_functions_base.multipoles import cartesian_to_spherical, realsphharm2df

def get_j0_j2_parameters(symbols):
    """Get <j0>, <j2> parameters by symbols."""
    j0_parameters = numpy.zeros((7, symbols.shape[0]), dtype=float)
    j2_parameters = numpy.zeros((7, symbols.shape[0]), dtype=float)
    d_mff = DATABASE["Magnetic Form Factor, tabulated"]
    for i_symbol, symbol in enumerate(symbols):
        atom_name, ion_charge, isotope_number = get_atomic_symbol_ion_charge_isotope_number_by_ion_symbol(symbol)
        try:
            j0_parameters[:, i_symbol] = d_mff[("j0_AaBbCcD", atom_name, ion_charge)]
        except KeyError:
            pass
        try:
            j2_parameters[:, i_symbol] = d_mff[("j2_AaBbCcD", atom_name, ion_charge)]
        except KeyError:
            pass
    return j0_parameters, j2_parameters

def calc_j0(sthovl, kappa, j0_parameters,
        flag_sthovl: bool=False, flag_kappa: bool=False):
    """Calculate j0."""
    j0_A = j0_parameters[0]
    j0_a = j0_parameters[1]
    j0_B = j0_parameters[2]
    j0_b = j0_parameters[3]
    j0_C = j0_parameters[4]
    j0_c = j0_parameters[5]
    j0_D = j0_parameters[6]
    _h = numpy.square(sthovl/kappa)
    j0_av = (j0_A*numpy.exp(-j0_a*_h) +
             j0_B*numpy.exp(-j0_b*_h) +
             j0_C*numpy.exp(-j0_c*_h)+j0_D)
    dder = {}
    if flag_sthovl:
        s_k = sthovl/kappa
        dder["sthovl"] = -2.*s_k*(j0_A*numpy.exp(-j0_a*_h)*j0_a +
             j0_B*numpy.exp(-j0_b*_h)*j0_b +
             j0_C*numpy.exp(-j0_c*_h)*j0_c)/kappa
    if flag_kappa:
        s_k = sthovl/kappa
        dder["kappa"] = 2.*_h*(j0_A*numpy.exp(-j0_a*_h)*j0_a +
             j0_B*numpy.exp(-j0_b*_h)*j0_b +
             j0_C*numpy.exp(-j0_c*_h)*j0_c)/kappa
    return j0_av, dder

def calc_j2(sthovl, kappa, j2_parameters,
        flag_sthovl: bool=False, flag_kappa: bool=False):
    """Calculate j2."""
    j2_A = j2_parameters[0]
    j2_a = j2_parameters[1]
    j2_B = j2_parameters[2]
    j2_b = j2_parameters[3]
    j2_C = j2_parameters[4]
    j2_c = j2_parameters[5]
    j2_D = j2_parameters[6]
    _h = numpy.square(sthovl/kappa)
    j2_av = (j2_A*numpy.exp(-j2_a*_h) +
             j2_B*numpy.exp(-j2_b*_h) +
             j2_C*numpy.exp(-j2_c*_h) + j2_D)*_h
    dder = {}
    if flag_sthovl:
        s_k = sthovl/kappa
        dder["sthovl"] = -2.*s_k*_h*(j2_A*numpy.exp(-j2_a*_h)*j2_a +
             j2_B*numpy.exp(-j2_b*_h)*j2_b +
             j2_C*numpy.exp(-j2_c*_h)*j2_c)/kappa + (j2_A*numpy.exp(-j2_a*_h) +
             j2_B*numpy.exp(-j2_b*_h) +
             j2_C*numpy.exp(-j2_c*_h) + j2_D)/kappa
    if flag_kappa:
        s_k = sthovl/kappa
        dder["kappa"] = 2.*_h*_h*(j2_A*numpy.exp(-j2_a*_h)*j2_a +
             j2_B*numpy.exp(-j2_b*_h)*j2_b +
             j2_C*numpy.exp(-j2_c*_h)*j2_c)/kappa - s_k*(j2_A*numpy.exp(-j2_a*_h) +
             j2_B*numpy.exp(-j2_b*_h) +
             j2_C*numpy.exp(-j2_c*_h) + j2_D)/kappa
    return j2_av, dder


def calc_form_factor(
        sthovl, lande_factor, kappa, j0_parameters, j2_parameters,
        flag_only_orbital=False,
        flag_sthovl: bool=False, flag_kappa: bool=False,
        flag_lande_factor: bool=False):
    r"""Calculate magnetic form factor in frame of Spherical model.
    (Int.Tabl.C.p.592)

    Calculation of magnetic form factor
    `<https://journals.aps.org/prb/pdf/10.1103/PhysRevB.79.140405>`
    mismatch with international tables where (1.0-2.0/np_factor_lande)
    """
    # not sure about kappa, it is here just for test, by default it is 1.0
    j2_av, dder_j2 = calc_j2(
        sthovl, kappa, j2_parameters, flag_sthovl=flag_sthovl,
        flag_kappa=flag_kappa)
    form_factor_orbital = (2.0/lande_factor-1.0)*j2_av

    if flag_only_orbital:
        form_factor = form_factor_orbital
    else:
        j0_av, dder_j0 = calc_j0(
            sthovl, kappa, j0_parameters, flag_sthovl=flag_sthovl,
            flag_kappa=flag_kappa)
        form_factor = j0_av + form_factor_orbital
    dder = {}
    if flag_sthovl:
        if flag_only_orbital:
            dder["sthovl"] = dder_j2["sthovl"]*(2.0/lande_factor-1.0)
        else:
            dder["sthovl"] = dder_j0["sthovl"] + dder_j2["sthovl"]*(2.0/lande_factor-1.0)
    if flag_kappa:
        if flag_only_orbital:
            dder["kappa"] = dder_j2["kappa"]*(2.0/lande_factor-1.0)
        else:
            dder["kappa"] = dder_j0["kappa"] + dder_j2["kappa"]*(2.0/lande_factor-1.0)
    if flag_lande_factor:
        dder["lande_factor"]=-2.0*j2_av/numpy.square(lande_factor)
    return form_factor, dder



def calc_magnetic_form_factor_by_multipole_model(index_hkl, r_d, 
        atom_rho_multipole_transformation_matrix, atom_rho_multipole_plm, atom_rho_multipole_lm, 
        l_kappa, l_n, l_zeta, l_coeff, unit_cell_parameters):   
    
    q_ccs = calc_q_ccs_by_unit_cell_parameters(unit_cell_parameters=unit_cell_parameters, index_hkl=index_hkl)[0]
    sthovl = 0.5*numpy.sqrt(numpy.sum(numpy.square(q_ccs), axis=0))
    eq_ccs = q_ccs/numpy.expand_dims(2.*sthovl, axis=0)
    m_m = calc_m_m_by_unit_cell_parameters(unit_cell_parameters=unit_cell_parameters, flag_unit_cell_parameters=False)[0]
    r_xyz = calc_m1_m2_inv_m1(m_m, r_d)[0]

    l_ff = []
    l_max = 4
    for kappa, n, zeta, coeff, transformation_matrix, np_plm in zip(l_kappa, l_n, l_zeta, l_coeff, atom_rho_multipole_transformation_matrix.T, atom_rho_multipole_plm.T):
        inv_r_g_to_l = calc_inv_m(transformation_matrix)[0]
        # for coordinate the inversed one is used
        eq_ccs_local = calc_m_v(inv_r_g_to_l, eq_ccs, flag_m=False, flag_v=False)[0]

        # FIXME: it should be checked
        trt_xyz = calc_m1_m2_inv_m1(transformation_matrix, r_xyz)[0]

        # [3, N_hkl, N_symm]
        eq_ccs_local_s = calc_m_v(
            numpy.expand_dims(trt_xyz, axis=1), 
            numpy.expand_dims(eq_ccs_local, axis=2), 
            flag_m=False, flag_v=False)[0]
        
        r_sph, theta_sph, phi_sph = cartesian_to_spherical(eq_ccs_local_s[0], eq_ccs_local_s[1], eq_ccs_local_s[2])
        np_ff = numpy.zeros(shape=r_sph.shape, dtype=float)
        jl = calc_jl(sthovl, coeff, n, zeta, kappa=kappa, l_max = l_max).transpose()
        for lm, plm in zip(atom_rho_multipole_lm, np_plm):
            if numpy.any(numpy.logical_not(numpy.isclose(plm,0)),axis=0):
                np_orient = plm * realsphharm2df(int(lm[0]), int(lm[1]), theta_sph, phi_sph)
                # FIXME: .real is not correct but for 0,2,4,6 should be ok
                np_ff +=  numpy.power(1j,lm[0]).real * np_orient * numpy.expand_dims(jl[int(lm[0])], axis=1) *  4 * numpy.pi
        l_ff.append(np_ff)
    # [N_hkl, N_symm, N_atom_rho_multipole]
    magnetic_form_factor = numpy.stack(l_ff, axis=2)
    dder = {}
    return magnetic_form_factor, dder
