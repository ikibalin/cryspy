"""
:Constants:


:Functions:



"""
import os
import numpy
import warnings

import pycifstar
from typing import List, Tuple

F_BSCAT = os.path.join(os.path.dirname(__file__), "bscat.tab")
BSCAT = pycifstar.to_loop(F_BSCAT)


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


def get_scat_length_neutron(type_n):
    """
    Take scat_length_neutron.
    """

    str_1 = type_n.strip().lower()
    str_1 = "".join([hh if hh.isalpha() else ' ' for hh in str_1 ]).split(" ")[0]

    flag = False
    for _1, _2 in zip(BSCAT["_atom_type_symbol"], BSCAT["_atom_type_cohb"]):
        if (_1.lower() == str_1):
            res = 0.1 * complex(_2)  # in 10**-12cm
            flag = True
        elif flag:
            break
    if not(flag):
        res = 0.
        warnings.warn(
            f"Can not find b_scat for '{type_n:}'.\n It is putted as 0.",
            UserWarning, stacklevel=2)
    return res
