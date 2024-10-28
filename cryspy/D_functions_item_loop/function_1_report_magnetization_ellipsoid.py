from typing import Union
import math
import numpy

from cryspy.A_functions_base.function_1_strings import \
    value_error_mark_to_string

import cryspy.A_functions_base.local_susceptibility as local_susceptibility

from cryspy.C_item_loop_classes.cl_1_cell import Cell
from cryspy.C_item_loop_classes.cl_1_atom_site_susceptibility import \
    AtomSiteSusceptibility, AtomSiteSusceptibilityL

def magnetization_ellipsoid_by_u_ij(
        cell: Cell,
        atom_site_susceptibility: AtomSiteSusceptibility):
    """Calc magnetization ellispoid through thermal parameters U_ij.
    """
    ucp = cell.get_unit_cell_parameters()
    
    q_rn = numpy.array([
        atom_site_susceptibility.chi_11, atom_site_susceptibility.chi_22, atom_site_susceptibility.chi_33,
        atom_site_susceptibility.chi_12, atom_site_susceptibility.chi_13, atom_site_susceptibility.chi_23], dtype=float)
    

    u_ij = local_susceptibility.calc_magnetization_ellipsoid_as_u(
        q_rn, ucp)[0]

    return u_ij


def report_main_axes_of_magnetization_ellipsoids(
        cell: Cell,
        atom_site_susceptibility: Union[AtomSiteSusceptibility, AtomSiteSusceptibilityL]):
    """Give a report about magnetizatin ellipsoid"""

    ls_out = []
    ucp = cell.get_unit_cell_parameters()
    a_s_m_a = atom_site_susceptibility

    if isinstance(atom_site_susceptibility, AtomSiteSusceptibilityL):
        a_s_s_items = atom_site_susceptibility.items
    elif isinstance(atom_site_susceptibility, AtomSiteSusceptibility):
        a_s_s_items = [atom_site_susceptibility, ]
    else:
        return ""

    for a_s_s in a_s_s_items:
        chi_xyz = numpy.array([a_s_s.chi_11, a_s_s.chi_22, a_s_s.chi_33,
        a_s_s.chi_12, a_s_s.chi_13, a_s_s.chi_23], dtype=float)

        moments, moments_sigma, eig_fields, eig_moments = \
            a_s_s.calc_main_axes_of_magnetization_ellipsoid(cell)
        
        u_ij = magnetization_ellipsoid_by_u_ij(cell, a_s_s)

        label = a_s_s.label 
        ls_out.append(f"For **`{label:}`** the principal axes of the atomic susceptibility tensor are:\n")
        # cycle over three main axes
        directions = eig_fields.transpose()

        ef = local_susceptibility.calc_ellipsoid_factor(chi_xyz, ucp)

        ls_out.append("|Axes  (mu_B/T)|Orientation: |X along inv.a| Y is [inv.a, c]|Z along c|")
        ls_out.append("|--------------|-------|----------|----------|----------|")
        for _val1, val_sigma, _direction in zip(
                moments, moments_sigma, directions):
            if math.isclose(val_sigma, 0.):
                s_val = f"{_val1: 9.5f}".rjust(9)
            else:
                s_param = value_error_mark_to_string(_val1, val_sigma, "")
                s_val = f"{s_param:}".rjust(9)
            ls_out.append(
                f"|  {s_val:} | along:| {_direction[0]: 9.5f}| \
{_direction[1]: 9.5f}| {_direction[2]: 9.5f}|")
        ls_out.append(f"\nEllipsoid factor {ef:.5f}.\n")

        ls_out.append("\nUse thermal parameters U_ij to plot ellispoid.\n")
        ls_out.append(f"|      U_11|      U_22|      U_33|      U_12|      \
U_13|      U_23|")
        ls_out.append(f"|----------|----------|----------|----------|------\
----|----------|")
        ls_out.append(
                f"|{u_ij[0]:9.3f}|{u_ij[1]:9.3f}|\
{u_ij[2]:9.3f}|{u_ij[3]:9.3f}|{u_ij[4]:9.3f}|\
{u_ij[5]:9.3f}|\n")
    # ls_out.append("Cartezian coordinate system is x||a*, z||c.")
    if len(ls_out) != 0:
        ls_out.insert(0, "# Magnetization ellipsoids\n")
    return "\n".join(ls_out)
