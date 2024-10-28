"""
:Constants:


:Functions:



"""

from cryspy.A_functions_base.database import DATABASE
from cryspy.A_functions_base.charge_form_factor import get_atomic_symbol_ion_charge_isotope_number_by_ion_symbol


def apply_constraint_on_cell_by_type_cell(cell, type_cell:str,
                                          it_coordinate_system_code:str):
    length_a, length_b, length_c = cell.length_a, cell.length_b, cell.length_c
    angle_alpha, angle_beta = cell.angle_alpha, cell.angle_beta
    angle_gamma = cell.angle_gamma
    if length_a is None: 
        cell.length_a = 1.
    if length_b is None: 
        cell.length_b = 1.
    if length_c is None: 
        cell.length_c = 1.
    if angle_alpha is None: 
        cell.angle_alpha = 90.
    if angle_beta is None: 
        cell.angle_beta = 90.
    if angle_gamma is None: 
        cell.angle_gamma = 90.

    if type_cell == "aP":
        pass
    elif type_cell.startswith("m"):
        cell.angle_alpha, cell.angle_alpha_refinement, cell.angle_alpha_constraint = 90., False, True
        cell.angle_gamma, cell.angle_gamma_refinement, cell.angle_gamma_constraint = 90., False, True
    elif type_cell.startswith("o"):
        cell.angle_alpha, cell.angle_alpha_refinement, cell.angle_alpha_constraint = 90., False, True
        cell.angle_beta, cell.angle_beta_refinement, cell.angle_beta_constraint = 90., False, True
        cell.angle_gamma, cell.angle_gamma_refinement, cell.angle_gamma_constraint = 90., False, True
    elif ((type_cell.startswith("t")) | (type_cell == "hP")):
        cell.length_b, cell.length_b_sigma = cell.length_a, cell.length_a_sigma
        cell.length_b_refinement, cell.length_b_constraint = False, True
        cell.angle_alpha, cell.angle_alpha_refinement, cell.angle_alpha_constraint = 90., False, True
        cell.angle_beta, cell.angle_beta_refinement, cell.angle_beta_constraint = 90., False, True
        cell.angle_gamma, cell.angle_gamma_refinement, cell.angle_gamma_constraint = 90., False, True
    elif (type_cell == "hR"):
        if it_coordinate_system_code.lower() == "h":
            cell.length_b, cell.length_b_sigma = cell.length_a, cell.length_a_sigma
            cell.length_b_refinement, cell.length_b_constraint = False, True
            cell.angle_alpha, cell.angle_alpha_refinement, cell.angle_alpha_constraint = 90., False, True
            cell.angle_beta, cell.angle_beta_refinement, cell.angle_beta_constraint = 90., False, True
            cell.angle_gamma, cell.angle_gamma_refinement, cell.angle_gamma_constraint = 120., False, True
        else:
            cell.length_b, cell.length_b_sigma = cell.length_a, cell.length_a_sigma
            cell.length_b.refinement, cell.length_b_constraint = False, True
            cell.length_c, cell.length_c_sigma = cell.length_a, cell.length_a_sigma
            cell.length_c_refinement, cell.length_c_constraint = False, True
            cell.angle_beta, cell.angle_beta_sigma = cell.angle_alpha, cell.angle_alpha_sigma
            cell.angle_beta_refinement, cell.angle_beta_constraint = False, True
            cell.angle_gamma, cell.angle_gamma_sigma = cell.angle_alph, cell.angle_alpha_sigma
            cell.angle_gamma_refinement, cell.angle_gamma_constraint = False, True
    elif type_cell.startswith("c"):
        cell.length_b, cell.length_b_sigma = cell.length_a, cell.length_a_sigma
        cell.length_b.refinement, cell.length_b_constraint = False, True
        cell.length_c, cell.length_c_sigma = cell.length_a, cell.length_a_sigma
        cell.length_c_refinement, cell.length_c_constraint = False, True
        cell.angle_alpha, cell.angle_alpha_refinement, cell.angle_alpha_constraint = 90., False, True
        cell.angle_beta, cell.angle_beta_refinement, cell.angle_beta_constraint = 90., False, True
        cell.angle_gamma, cell.angle_gamma_refinement, cell.angle_gamma_constraint = 90., False, True


def get_scat_length_neutron(type_n: str):
    """
    Take scat_length_neutron.
    """
    atomic_symbol, ion_charge, isotope_number = get_atomic_symbol_ion_charge_isotope_number_by_ion_symbol(type_n)
    if isotope_number is None:
        s_isotope = atomic_symbol
    else:
        s_isotope = f"{isotope_number:}{atomic_symbol:}"
    res = DATABASE["Isotopes"][("b_scat", s_isotope)]
    return res
