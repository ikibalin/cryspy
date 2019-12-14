"""
:Constants:


:Functions:



"""
import os
import numpy
import warnings

import pycifstar
from typing import List, Tuple

from cryspy.common.cl_fitable import Fitable

F_BSCAT = os.path.join(os.path.dirname(__file__), "bscat.tab")
BSCAT = pycifstar.to_loop(F_BSCAT)


def apply_constraint_on_cell_by_type_cell(cell, type_cell:str):
    length_a, length_b, length_c = cell.length_a, cell.length_b, cell.length_c
    angle_alpha, angle_beta, angle_gamma = cell.angle_alpha, cell.angle_beta, cell.angle_gamma
    if length_a is None: 
        length_a = Fitable(1.)
        cell.length_a = length_a
    if length_b is None: 
        length_b = Fitable(1.)
        cell.length_b = length_b
    if length_c is None: 
        length_c = Fitable(1.)
        cell.length_c = length_c
    if angle_alpha is None: 
        angle_alpha = Fitable(90.)
        cell.angle_alpha = angle_alpha
    if angle_beta is None: 
        angle_beta = Fitable(90.)
        cell.angle_beta = angle_beta
    if angle_gamma is None: 
        angle_gamma = Fitable(90.)
        cell.angle_gamma = angle_gamma

    if type_cell == "aP":
        pass
    elif type_cell.startswith("m"):
        angle_alpha.value, angle_alpha.refinement, angle_alpha.constraint_flag = 90., False, True
        angle_gamma.value, angle_gamma.refinement, angle_gamma.constraint_flag = 90., False, True
    elif type_cell.startswith("o"):
        angle_alpha.value, angle_alpha.refinement, angle_alpha.constraint_flag = 90., False, True
        angle_beta.value, angle_beta.refinement, angle_beta.constraint_flag = 90., False, True
        angle_gamma.value, angle_gamma.refinement, angle_gamma.constraint_flag = 90., False, True
    elif ((type_cell.startswith("t")) | (type_cell == "hP")):
        length_b.value, length_b.sigma = length_a.value, length_a.sigma
        length_b.refinement, length_b.constraint_flag = False, True
        angle_alpha.value, angle_alpha.refinement, angle_alpha.constraint_flag = 90., False, True
        angle_beta.value, angle_beta.refinement, angle_beta.constraint_flag = 90., False, True
        angle_gamma.value, angle_gamma.refinement, angle_gamma.constraint_flag = 90., False, True
    elif (type_cell == "hR"):
        length_b.value, length_b.sigma = length_a.value, length_a.sigma
        length_b.refinement, length_b.constraint_flag = False, True
        length_c.value, length_c.sigma = length_a.value, length_a.sigma
        length_c.refinement, length_c.constraint_flag = False, True
        angle_beta.value, angle_beta.sigma = angle_alpha.value, angle_alpha.sigma
        angle_beta.refinement, angle_beta.constraint_flag = False, True
        angle_gamma.value, angle_gamma.sigma = angle_alpha.value, angle_alpha.sigma
        angle_gamma.refinement, angle_gamma.constraint_flag = False, True
    elif type_cell.startswith("c"):
        length_b.value, length_b.sigma = length_a.value, length_a.sigma
        length_b.refinement, length_b.constraint_flag = False, True
        length_c.value, length_c.sigma = length_a.value, length_a.sigma
        length_c.refinement, length_c.constraint_flag = False, True
        angle_alpha.value, angle_alpha.refinement, angle_alpha.constraint_flag = 90., False, True
        angle_beta.value, angle_beta.refinement, angle_beta.constraint_flag = 90., False, True
        angle_gamma.value, angle_gamma.refinement, angle_gamma.constraint_flag = 90., False, True
    return 
        

def get_scat_length_neutron(type_n):
    """
Take scat_length_neutron
    """
    
    str_1 = type_n.strip().lower()
    str_1 = "".join([hh if hh.isalpha() else ' ' for hh in str_1 ]).split(" ")[0]
    
    flag = False
    for _1, _2 in zip(BSCAT["_atom_type_symbol"], BSCAT["_atom_type_cohb"]):
        if (_1.lower() == str_1):
            res = complex(_2)
            flag = True
        elif flag:
            break
    if not(flag):
        res = 0.
        warnings.warn(f"Can not find b_scat for '{type_n:}'.\n It is putted as 0.", UserWarning, stacklevel=2)
    return res 





