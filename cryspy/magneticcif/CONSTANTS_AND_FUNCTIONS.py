"""
:Constants:

:Functions:

"""
import os
import numpy
import warnings
from fractions import Fraction
from typing import List, Tuple
import pycifstar

F_FORMMAG = os.path.join(os.path.dirname(__file__), "formmag.tab")
FORMMAG_data = pycifstar.to_data(F_FORMMAG)
    

def get_j0_j2_by_symbol(symbol:str):
    s_1 = "_atom_type_scat_symbol"
    j0_a2, j0_A1, j0_B1, j0_b2 = 0., 0., 0., 0.
    j0_C1, j0_c2, j0_D = 0., 0., 0. 
    j2_a2, j2_A1, j2_B1, j2_b2 = 0., 0., 0., 0.
    j2_C1, j2_c2, j2_D = 0., 0., 0.
    flag_0, flag_2 = False, False
    for _loop in FORMMAG_data.loops:
        if "_atom_type_scat_neutron_magnetic_j0_a1" in _loop.names:
            hh = [_i1 for _i1, _1 in enumerate(_loop[s_1]) if (_1 == symbol)]
            if len(hh) > 0:
                _ind = hh[0]
                j0_A1 = float(_loop["_atom_type_scat_neutron_magnetic_j0_A1"][_ind])
                j0_a2 = float(_loop["_atom_type_scat_neutron_magnetic_j0_a2"][_ind])
                j0_B1 = float(_loop["_atom_type_scat_neutron_magnetic_j0_B1"][_ind])
                j0_b2 = float(_loop["_atom_type_scat_neutron_magnetic_j0_b2"][_ind])
                j0_C1 = float(_loop["_atom_type_scat_neutron_magnetic_j0_C1"][_ind])
                j0_c2 = float(_loop["_atom_type_scat_neutron_magnetic_j0_c2"][_ind])
                j0_D = float(_loop["_atom_type_scat_neutron_magnetic_j0_D"][_ind])
                flag_0 = True
        if "_atom_type_scat_neutron_magnetic_j2_a1" in _loop.names:
            hh = [_i1 for _i1, _1 in enumerate(_loop[s_1]) if (_1 == symbol)]
            if len(hh) > 0:
                _ind = hh[0]
                j2_A1 = float(_loop["_atom_type_scat_neutron_magnetic_j2_A1"][_ind])
                j2_a2 = float(_loop["_atom_type_scat_neutron_magnetic_j2_a2"][_ind])
                j2_B1 = float(_loop["_atom_type_scat_neutron_magnetic_j2_B1"][_ind])
                j2_b2 = float(_loop["_atom_type_scat_neutron_magnetic_j2_b2"][_ind])
                j2_C1 = float(_loop["_atom_type_scat_neutron_magnetic_j2_C1"][_ind])
                j2_c2 = float(_loop["_atom_type_scat_neutron_magnetic_j2_c2"][_ind])
                j2_D = float(_loop["_atom_type_scat_neutron_magnetic_j2_D"][_ind])
                flag_2 = True
        if all([flag_0, flag_2]):
            break
    return j0_A1, j0_a2, j0_B1, j0_b2, j0_C1, j0_c2, j0_D, j2_A1, j2_a2, j2_B1, j2_b2, j2_C1, j2_c2, j2_D


