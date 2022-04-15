"""
Base crystallographic functions.

Functions
---------
    - calc_beta_by_u
    - vibration_constraints,
    - apply_constraint_on_cell_by_type_cell

They are stored in the list FUNCTIONS
"""
from typing import Tuple
import numpy


def calc_beta_by_u(u_i, cell):
    r"""
    Calculate beta.

    u_11, u_22, u_33, u_12, u_13, u_23 = u_i

    cell is Cell class

    Output:
        $\beta_{11}$, $\beta_{22}$, $\beta_{33}$
        $\beta_{12}$, $\beta_{13}$, $\beta_{23}$
    """
    (u_11, u_22, u_33, u_12, u_13, u_23) = u_i

    ia = cell.reciprocal_length_a
    ib = cell.reciprocal_length_b
    ic = cell.reciprocal_length_c
    beta_11 = 2.*numpy.pi**2*u_11*ia**2
    beta_22 = 2.*numpy.pi**2*u_22*ib**2
    beta_33 = 2.*numpy.pi**2*u_33*ic**2
    beta_12 = 2.*numpy.pi**2*u_12*ia*ib
    beta_13 = 2.*numpy.pi**2*u_13*ia*ic
    beta_23 = 2.*numpy.pi**2*u_23*ib*ic
    return beta_11, beta_22, beta_33, beta_12, beta_13, beta_23


def vibration_constraints(numb, param_i, sigma_i, ref_i):
    """
    Constraints on atomic vibrations.

    Parameters
    ----------
    numb : TYPE
        DESCRIPTION.
    param_i : TYPE
        DESCRIPTION.
    sigma_i : TYPE
        DESCRIPTION.
    ref_i : TYPE
        DESCRIPTION.

    Returns
    -------
    param_i : TYPE
        DESCRIPTION.
    sigma_i : TYPE
        DESCRIPTION.
    ref_i : TYPE
        DESCRIPTION.
    constr_i : TYPE
        DESCRIPTION.

    """
    p_11, p_22, p_33, p_12, p_13, p_23 = param_i
    s_11, s_22, s_33, s_12, s_13, s_23 = sigma_i
    p_11_ref, p_22_ref, p_33_ref, p_12_ref, p_13_ref, p_23_ref = ref_i
    p_11_constr, p_22_constr, p_33_constr = False, False, False
    p_12_constr, p_13_constr, p_23_constr = False, False, False

    if numb == 1:
        p_12 = 0.
        s_12 = 0.
        p_12_ref = False
        p_12_constr = True
        p_23 = 0.
        s_23 = 0.
        p_23_ref = False
        p_23_constr = True
    elif numb == 2:
        p_13 = 0.
        s_13 = 0.
        p_13_ref = False
        p_13_constr = True
        p_23 = 0.
        s_23 = 0.
        p_23_ref = False
        p_23_constr = True
    elif numb == 3:
        p_12 = 0.
        s_12 = 0.
        p_12_ref = False
        p_12_constr = True
        p_13 = 0.
        s_13 = 0.
        p_13_ref = False
        p_13_constr = True
    elif numb == 4:
        p_12 = 0.
        s_12 = 0.
        p_12_ref = False
        p_12_constr = True
        p_13 = 0.
        s_13 = 0.
        p_13_ref = False
        p_13_constr = True
        p_23 = 0.
        s_23 = 0.
        p_23_ref = False
        p_23_constr = True
    elif numb == 5:
        p_22 = p_11
        s_22 = s_11
        p_22_ref = False
        p_22_constr = True
        p_13 = 0.
        s_13 = 0.
        p_13_ref = False
        p_13_constr = True
        p_23 = 0.
        s_23 = 0.
        p_23_ref = False
        p_23_constr = True
    elif numb == 6:
        p_22 = p_11
        s_22 = s_11
        p_22_ref = False
        p_22_constr = True
        p_23 = p_13
        s_23 = s_13
        p_23_ref = False
        p_23_constr = True
    elif numb == 7:
        p_22 = p_11
        s_22 = s_11
        p_22_ref = False
        p_22_constr = True
        p_23 = -1.*p_13
        s_23 = s_13
        p_23_ref = False
        p_23_constr = True
    elif numb == 8:
        p_22 = p_11
        s_22 = s_11
        p_22_ref = False
        p_22_constr = True
        p_12 = 0.
        s_12 = 0.
        p_12_ref = False
        p_12_constr = True
        p_13 = 0.
        s_13 = 0.
        p_13_ref = False
        p_13_constr = True
        p_23 = 0.
        s_23 = 0.
        p_23_ref = False
        p_23_constr = True
    elif numb == 9:
        p_33 = p_22
        s_33 = s_22
        p_33_ref = False
        p_33_constr = True
        p_12 = 0.
        s_12 = 0.
        p_12_ref = False
        p_12_constr = True
        p_13 = 0.
        s_13 = 0.
        p_13_ref = False
        p_13_constr = True
    elif numb == 10:
        p_33 = p_22
        s_33 = s_22
        p_33_ref = False
        p_33_constr = True
        p_13 = p_12
        s_13 = s_12
        p_13_ref = False
        p_13_constr = True
    elif numb == 11:
        p_33 = p_22
        s_33 = s_22
        p_33_ref = False
        p_33_constr = True
        p_13 = -1.*p_12
        s_13 = s_12
        p_13_ref = False
        p_13_constr = True
    elif numb == 12:
        p_33 = p_22
        s_33 = s_22
        p_33_ref = False
        p_33_constr = True
        p_12 = 0.
        s_12 = 0.
        p_12_ref = False
        p_12_constr = True
        p_13 = 0.
        s_13 = 0.
        p_13_ref = False
        p_13_constr = True
        p_23 = 0.
        s_23 = 0.
        p_23_ref = False
        p_23_constr = True
    elif numb == 13:
        p_12 = 0.5*p_22         
        s_12 = 0.5*s_22         
        p_12_ref = False
        p_12_constr = True
        p_23 = 0.
        s_23 = 0.
        p_23_ref = False
        p_23_constr = True
    elif numb == 14:
        p_12 = 0.5*p_22
        s_12 = 0.5*s_22
        p_12_ref = False
        p_12_constr = True
        p_13 = 0.
        s_13 = 0.
        p_13_ref = False
        p_13_constr = True
        p_23 = 0.
        s_23 = 0.
        p_23_ref = False
        p_23_constr = True
    elif numb == 15:
        p_12 = 0.5*p_22
        s_12 = 0.5*s_22
        p_12_ref = False
        p_12_constr = True
        p_23 = 2.0*p_13
        s_23 = 2.0*s_13
        p_23_ref = False
        p_23_constr = True
    elif numb == 16:
        p_22 = p_11
        s_22 = s_11
        p_22_ref = False
        p_22_constr = True
        p_12 = 0.5*p_11      
        s_12 = 0.5*s_11      
        p_12_ref = False
        p_12_constr = True
        p_13 = 0.
        s_13 = 0.
        p_13_ref = False
        p_13_constr = True
        p_23 = 0.
        s_23 = 0.
        p_23_ref = False
        p_23_constr = True
    elif numb == 17:
        p_22 = p_11
        s_22 = s_11
        p_22_ref = False
        p_22_constr = True
        p_33 = p_11
        s_33 = s_11
        p_33_ref = False
        p_33_constr = True
        p_12 = 0.
        s_12 = 0.
        p_12_ref = False
        p_12_constr = True
        p_13 = 0.
        s_13 = 0.
        p_13_ref = False
        p_13_constr = True
        p_23 = 0.
        s_23 = 0.
        p_23_ref = False
        p_23_constr = True
    elif numb == 18:
        p_22 = p_11
        s_22 = s_11
        p_22_ref = False
        p_22_constr = True
        p_33 = p_11
        s_33 = s_11
        p_33_ref = False
        p_33_constr = True
        p_13 = p_12
        s_13 = s_12
        p_13_ref = False
        p_13_constr = True
        p_23 = p_12
        s_23 = s_12
        p_23_ref = False
        p_23_constr = True

    param_i = (p_11, p_22, p_33, p_12, p_13, p_23)
    sigma_i = (s_11, s_22, s_33, s_12, s_13, s_23)
    ref_i = (p_11_ref, p_22_ref, p_33_ref, p_12_ref, p_13_ref, p_23_ref)
    constr_i = (p_11_constr, p_22_constr, p_33_constr, p_12_constr,
                p_13_constr, p_23_constr)
    return param_i, sigma_i, ref_i, constr_i


def apply_constraint_on_cell_by_type_cell(
        cell_param: Tuple[float], cell_sigma: Tuple[float],
        cell_ref: Tuple[bool], type_cell: str,
        it_coordinate_system_code: str) -> (Tuple[float], Tuple[float],
                                            Tuple[bool], Tuple[bool]):
    """
       Apply constraints.

    Parameters
    ----------
    cell_param : Tuple[float]
        (a, b, c, alpha, beta, gamma).
    cell_sigma : Tuple[float]
        (a_sigma, b_sigma, c_sigma, alpha_sigma, beta_sigma, gamma_sigma).
    cell_ref : Tuple[bool]
        (a_sigma, b_sigma, c_sigma, alpha_sigma, beta_sigma, gamma_sigma).
    type_cell : str
        DESCRIPTION.
    it_coordinate_system_code : str
        DESCRIPTION.

    Returns
    -------
    (Tuple[float], Tuple[float], Tuple[bool], Tuple[bool])
        cell_param, cell_sigma, cell_ref, cell_constr.

    """
    (a, b, c, alpha, beta, gamma) = cell_param
    (a_sig, b_sig, c_sig, alpha_sig, beta_sig, gamma_sig) = cell_sigma
    (a_ref, b_ref, c_ref, alpha_ref, beta_ref, gamma_ref) = cell_ref
    a_constr, b_constr, c_constr, alpha_constr, beta_constr, gamma_constr = \
        False, False, False, False, False, False

    if type_cell == "aP":
        pass
    elif type_cell.startswith("m"):
        alpha, alpha_sig, alpha_ref, alpha_constr = 90., 0., False, True
        gamma, gamma_sig, gamma_ref, gamma_constr = 90., 0., False, True
    elif type_cell.startswith("o"):
        alpha, alpha_sig, alpha_ref, alpha_constr = 90., 0., False, True
        beta, beta_sig, beta_ref, beta_constr = 90., 0., False, True
        gamma, gamma_sig, gamma_ref, gamma_constr = 90., 0., False, True
    elif ((type_cell.startswith("t"))):  # FIXME: check  | (type_cell == "hP")
        b, b_sig, b_ref, b_constr = a, a_sig, False, True
        alpha, alpha_sig, alpha_ref, alpha_constr = 90., 0., False, True
        beta, beta_sig, beta_ref, beta_constr = 90., 0., False, True
        gamma, gamma_sig, gamma_ref, gamma_constr = 90., 0., False, True
    elif ((type_cell.startswith("hP"))):
        if it_coordinate_system_code.lower() == "h":
            b, b_sig, b_ref, b_constr = a, a_sig, False, True
            alpha, alpha_sig, alpha_ref, alpha_constr = 90., 0., False, True
            beta, beta_sig, beta_ref, beta_constr = 90., 0., False, True
            gamma, gamma_sig, gamma_ref, gamma_constr = 120., 0., False, True
        else:
            b, b_sig, b_ref, b_constr = a, a_sig, False, True
            c, c_sig, c_ref, c_constr = a, a_sig, False, True
            beta, beta_sig, beta_ref, beta_constr = alpha, alpha_sig, \
                False, True
            gamma, gamma_sig, gamma_ref, gamma_constr = alpha, alpha_sig,\
                False, True
    elif (type_cell == "hR"):
        if it_coordinate_system_code.lower() == "h":
            b, b_sig, b_ref, b_constr = a, a_sig, False, True
            alpha, alpha_sig, alpha_ref, alpha_constr = 90., 0., False, True
            beta, beta_sig, beta_ref, beta_constr = 90., 0., False, True
            gamma, gamma_sig, gamma_ref, gamma_constr = 120., 0., False, True
        else:
            b, b_sig, b_ref, b_constr = a, a_sig, False, True
            c, c_sig, c_ref, c_constr = a, a_sig, False, True
            beta, beta_sig, beta_ref, beta_constr = alpha, alpha_sig, \
                False, True
            gamma, gamma_sig, gamma_ref, gamma_constr = alpha, alpha_sig,\
                False, True
    elif type_cell.startswith("c"):
        b, b_sig, b_ref, b_constr = a, a_sig, False, True
        c, c_sig, c_ref, c_constr = a, a_sig, False, True
        alpha, alpha_sig, alpha_ref, alpha_constr = 90., 0., False, True
        beta, beta_sig, beta_ref, beta_constr = 90., 0., False, True
        gamma, gamma_sig, gamma_ref, gamma_constr = 90., 0., False, True
    cell_param = (a, b, c, alpha, beta, gamma)
    cell_sigma = (a_sig, b_sig, c_sig, alpha_sig, beta_sig, gamma_sig)
    cell_ref = (a_ref, b_ref, c_ref, alpha_ref, beta_ref, gamma_ref)
    cell_constr = (a_constr, b_constr, c_constr, alpha_constr, beta_constr,
                   gamma_constr)
    return cell_param, cell_sigma, cell_ref, cell_constr
