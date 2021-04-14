from typing import Union
import math
import numpy

from cryspy.A_functions_base.function_1_strings import \
    value_error_to_string

from cryspy.C_item_loop_classes.cl_1_cell import Cell
from cryspy.C_item_loop_classes.cl_1_atom_site_susceptibility import \
    AtomSiteSusceptibility, AtomSiteSusceptibilityL

def magnetization_ellipsoid_by_u_ij(
        cell: Cell,
        atom_site_susceptibility: AtomSiteSusceptibility):
    """Calc magnetization ellispoid through thermal parameters U_ij.
    """
    m_m_norm = cell.m_m_norm
    m_im_norm = numpy.linalg.inv(m_m_norm)
    chi_11 = atom_site_susceptibility.chi_11
    chi_22 = atom_site_susceptibility.chi_22
    chi_33 = atom_site_susceptibility.chi_33
    chi_12 = atom_site_susceptibility.chi_12
    chi_13 = atom_site_susceptibility.chi_13
    chi_23 = atom_site_susceptibility.chi_23

    m_chi_loc = numpy.array(
        [[chi_11, chi_12, chi_13],
         [chi_12, chi_22, chi_23],
         [chi_13, chi_23, chi_33]], dtype=float)

    m_chi_orto = numpy.matmul(numpy.matmul(m_m_norm, m_chi_loc),
                              m_m_norm.transpose())
    val, vec = numpy.linalg.eigh(m_chi_orto)
    
    # D = numpy.diag(numpy.abs(val))  # only positive eigenvalues
    D = numpy.diag(numpy.square(val))  # only positive eigenvalues
    R = numpy.array(vec)
    
    chi_as_u_orto = numpy.matmul(R, numpy.matmul(D, R.transpose()))
    chi_as_u_loc = numpy.matmul(numpy.matmul(m_im_norm, chi_as_u_orto),
                                m_im_norm.transpose())
    return chi_as_u_loc


def report_main_axes_of_magnetization_ellipsoids(
        cell: Cell,
        atom_site_susceptibility: Union[AtomSiteSusceptibility, AtomSiteSusceptibilityL]):
    """Give a report about magnetizatin ellipsoid"""

    ls_out = []

    a_s_m_a = atom_site_susceptibility

    if isinstance(atom_site_susceptibility, AtomSiteSusceptibilityL):
        a_s_s_items = atom_site_susceptibility.items
    elif isinstance(atom_site_susceptibility, AtomSiteSusceptibility):
        a_s_s_items = [atom_site_susceptibility, ]
    else:
        return ""

    for a_s_s in a_s_s_items:
        moments, moments_sigma, rot_matrix = \
            a_s_s.calc_main_axes_of_magnetization_ellipsoid(cell)
        
        chi_as_u = magnetization_ellipsoid_by_u_ij(cell, a_s_s)

        label = a_s_s.label 
        ls_out.append(f"For **`{label:}`** the susceptibility is:\n")
        # cycle over three main axes
        directions = rot_matrix.transpose()

        ls_out.append("|Susceptibility (mu_B/T)|Orientation: |X along inv.a| Y is [inv.a, c]|Z along c|")
        ls_out.append("|--------------|-------|----------|----------|----------|")
        for _val1, val_sigma, _direction in zip(
                moments, moments_sigma, directions):
            if math.isclose(val_sigma, 0.):
                s_val = f"{_val1: 9.5f}".rjust(9)
            else:
                s_param = value_error_to_string(_val1, val_sigma)
                s_val = f"{s_param:}".rjust(9)
            ls_out.append(
                f"|  {s_val:} | along:| {_direction[0]: 9.5f}| \
{_direction[1]: 9.5f}| {_direction[2]: 9.5f}|")
        ls_out.append("\nUse thermal parameters U_ij to plot magn. ellispoid.\n")
        ls_out.append(f"|     U_11|     U_22|     U_33|     U_12|     \
U_13|     U_23|")
        ls_out.append(f"|---------|---------|---------|---------|-----\
----|---------|")
        ls_out.append(
                f"|{chi_as_u[0, 0]:9.2f}|{chi_as_u[1, 1]:9.2f}|\
{chi_as_u[2, 2]:9.2f}|{chi_as_u[0, 1]:9.2f}|{chi_as_u[0, 2]:9.2f}|\
{chi_as_u[1, 2]:9.2f}|\n")
    # ls_out.append("Cartezian coordinate system is x||a*, z||c.")
    if len(ls_out) != 0:
        ls_out.insert(0, "# Magnetization ellipsoids\n")
    return "\n".join(ls_out)
