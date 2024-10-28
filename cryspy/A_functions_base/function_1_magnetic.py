"""Functions for magnetic section.

Functions
---------
    - get_j0_j2_by_symbol
"""
from cryspy.A_functions_base.database import DATABASE
from cryspy.A_functions_base.charge_form_factor import get_atomic_symbol_ion_charge_isotope_number_by_ion_symbol



def get_j0_j2_by_symbol(symbol: str):
    """Get <j0>, <j2> parameters by symbol."""

    j0_a2, j0_A1, j0_B1, j0_b2 = 0., 0., 0., 0.
    j0_C1, j0_c2, j0_D = 0., 0., 0.
    j2_a2, j2_A1, j2_B1, j2_b2 = 0., 0., 0., 0.
    j2_C1, j2_c2, j2_D = 0., 0., 0.
    
    d_mff = DATABASE["Magnetic Form Factor, tabulated"]
    atom_name, ion_charge, isotope_number = get_atomic_symbol_ion_charge_isotope_number_by_ion_symbol(symbol)
    try:
        j0 = d_mff[("j0_AaBbCcD", atom_name, ion_charge)]
        j0_A1, j0_a2, j0_B1, j0_b2, j0_C1, j0_c2, j0_D = j0[0], j0[1], j0[2], j0[3], j0[4], j0[5], j0[6] 
    except KeyError:
        pass
    try:
        j2 = d_mff[("j2_AaBbCcD", atom_name, ion_charge)]
        j2_A1, j2_a2, j2_B1, j2_b2, j2_C1, j2_c2, j2_D = j2[0], j2[1], j2[2], j2[3], j2[4], j2[5], j2[6]
    except KeyError:
        pass
    return j0_A1, j0_a2, j0_B1, j0_b2, j0_C1, j0_c2, j0_D, j2_A1, j2_a2, \
        j2_B1, j2_b2, j2_C1, j2_c2, j2_D
