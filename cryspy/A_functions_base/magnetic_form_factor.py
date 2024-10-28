"""Functions for magnetic section.

Functions
---------
    - get_j0_j2_by_symbol
"""
import numpy

from cryspy.A_functions_base.database import DATABASE
from cryspy.A_functions_base.charge_form_factor import get_atomic_symbol_ion_charge_isotope_number_by_ion_symbol


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