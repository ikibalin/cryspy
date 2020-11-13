"""Functions for magnetic section.

Functions
---------
    - get_j0_j2_by_symbol
"""
import os
import pycifstar

F_FORMMAG = os.path.join(os.path.dirname(__file__), "formmag.tab")
FORMMAG_data = pycifstar.to_data(F_FORMMAG)


def get_j0_j2_by_symbol(symbol: str):
    """Get <j0>, <j2> parameters by symbol."""
    s_1 = "_atom_type_scat_symbol"
    j0_a2, j0_A1, j0_B1, j0_b2 = 0., 0., 0., 0.
    j0_C1, j0_c2, j0_D = 0., 0., 0.
    j2_a2, j2_A1, j2_B1, j2_b2 = 0., 0., 0., 0.
    j2_C1, j2_c2, j2_D = 0., 0., 0.
    flag_0, flag_2 = False, False
    for loop in FORMMAG_data.loops:
        if "_atom_type_scat_neutron_magnetic_j0_a1" in loop.names:
            hh = [_i1 for _i1, _1 in enumerate(loop[s_1]) if (_1 == symbol)]
            if len(hh) > 0:
                ind = hh[0]
                j0_A1 = float(loop["_atom_type_scat_neutron_magnetic_j0_A1"][
                    ind])
                j0_a2 = float(loop["_atom_type_scat_neutron_magnetic_j0_a2"][
                    ind])
                j0_B1 = float(loop["_atom_type_scat_neutron_magnetic_j0_B1"][
                    ind])
                j0_b2 = float(loop["_atom_type_scat_neutron_magnetic_j0_b2"][
                    ind])
                j0_C1 = float(loop["_atom_type_scat_neutron_magnetic_j0_C1"][
                    ind])
                j0_c2 = float(loop["_atom_type_scat_neutron_magnetic_j0_c2"][
                    ind])
                j0_D = float(loop["_atom_type_scat_neutron_magnetic_j0_D"][
                    ind])
                flag_0 = True
        if "_atom_type_scat_neutron_magnetic_j2_a1" in loop.names:
            hh = [_i1 for _i1, _1 in enumerate(loop[s_1]) if (_1 == symbol)]
            if len(hh) > 0:
                ind = hh[0]
                j2_A1 = float(loop["_atom_type_scat_neutron_magnetic_j2_A1"][
                    ind])
                j2_a2 = float(loop["_atom_type_scat_neutron_magnetic_j2_a2"][
                    ind])
                j2_B1 = float(loop["_atom_type_scat_neutron_magnetic_j2_B1"][
                    ind])
                j2_b2 = float(loop["_atom_type_scat_neutron_magnetic_j2_b2"][
                    ind])
                j2_C1 = float(loop["_atom_type_scat_neutron_magnetic_j2_C1"][
                    ind])
                j2_c2 = float(loop["_atom_type_scat_neutron_magnetic_j2_c2"][
                    ind])
                j2_D = float(loop["_atom_type_scat_neutron_magnetic_j2_D"][
                    ind])
                flag_2 = True
        if all([flag_0, flag_2]):
            break
    return j0_A1, j0_a2, j0_B1, j0_b2, j0_C1, j0_c2, j0_D, j2_A1, j2_a2, \
        j2_B1, j2_b2, j2_C1, j2_c2, j2_D
