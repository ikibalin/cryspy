import numpy

from cryspy.A_functions_base.unit_cell import calc_m_m_by_unit_cell_parameters
from cryspy.A_functions_base.matrix_operations import calc_m1_m2_inv_m1, calc_m_v, calc_mm_as_m_q_inv_m, calc_mm_as_m1_m2_inv_m1, calc_det_m

def calc_sc_fract_sc_b(symm_elems, atom_fract_xyz):
    sc_fract = (symm_elems[4:13]).sum(axis=1)/symm_elems.shape[1]
    b_s = symm_elems[:3]/(symm_elems.shape[1]*numpy.expand_dims(symm_elems[3], axis=0))

    atom_fract_xyz = numpy.mod(atom_fract_xyz, 1)
    r_s_x = calc_m_v(symm_elems[4:13], atom_fract_xyz, flag_m=False, flag_v=False)[0]
    n_s = numpy.expand_dims(atom_fract_xyz, axis=1)-r_s_x-b_s
    sc_b = (b_s + n_s).sum(axis=1)/n_s.shape[1]

    # sc_b = ().sum(axis=1)
    # x_new = calc_m_v(symm_elems[4:13], atom_fract_xyz, flag_m=False, flag_v=False)[0]
    # n_s, x0 = numpy.divmod(x_new,1)
    # n_s = -n_s.sum(axis=1)/n_s.shape[1]
    # sc_b = sc_b + n_s
    return sc_fract, sc_b


def calc_sc_beta(symm_elems):
    """Calculate \beta_av. = 1/N_s R_s \beta R_s^T
    """
    r_11, r_12, r_13 = symm_elems[4], symm_elems[5], symm_elems[6]
    r_21, r_22, r_23 = symm_elems[7], symm_elems[8], symm_elems[9]
    r_31, r_32, r_33 = symm_elems[10], symm_elems[11], symm_elems[12]

    r_r_t = numpy.stack([
        numpy.stack([r_11**2, r_12**2, r_13**2, 2*r_11*r_12, 2*r_11*r_13, 2*r_12*r_13], axis=0), # 11
        numpy.stack([r_21**2, r_22**2, r_23**2, 2*r_21*r_22, 2*r_21*r_23, 2*r_22*r_23], axis=0), # 22
        numpy.stack([r_31**2, r_32**2, r_33**2, 2*r_31*r_32, 2*r_31*r_33, 2*r_32*r_33], axis=0), # 33
        numpy.stack([r_11*r_21, r_12*r_22, r_13*r_23, r_11*r_22 + r_12*r_21, r_11*r_23 + r_13*r_21, r_12*r_23 + r_13*r_22], axis=0), # 12
        numpy.stack([r_11*r_31, r_12*r_32, r_13*r_33, r_11*r_32 + r_12*r_31, r_11*r_33 + r_13*r_31, r_12*r_33 + r_13*r_32], axis=0), # 13
        numpy.stack([r_21*r_31, r_22*r_32, r_23*r_33, r_21*r_32 + r_22*r_31, r_21*r_33 + r_23*r_31, r_22*r_33 + r_23*r_32], axis=0)], axis=0 # 23
    )
    sc_beta = r_r_t.sum(axis=2)/r_11.shape[0]
    return sc_beta


def calc_sc_chi(symm_elems, unit_cell_parameters, flag_unit_cell_parameters: bool = False):
    """Calculate symmetry constraint matrix for susceptibility
    It is supposed that determinant of r matrix is 1
    """
    m_m, dder_m_m = calc_m_m_by_unit_cell_parameters(
        unit_cell_parameters, flag_unit_cell_parameters=False)

    r_direct = symm_elems[4:13]

    r_ccs, dder_r_ccs = calc_m1_m2_inv_m1(m_m, r_direct, flag_m1=False, flag_m2=False)

    if symm_elems.shape[0] == 14:
        theta = numpy.expand_dims(symm_elems[13], axis=(0,1))
        # det_r_ccs = calc_det_m(r_ccs, flag_m=False)[0]
        # res = theta_s*det_r_ccs
    else:
        theta = numpy.expand_dims(numpy.ones_like(symm_elems[0]), axis=(0,1))

    res, dder_mm = calc_mm_as_m_q_inv_m(r_ccs, flag_m= flag_unit_cell_parameters)

    # if r_direct.shape[1] == 2:
    #     print("r_ccs: ", r_ccs)

    mm = (res*theta).sum(axis=2)/res.shape[2]

    sc_chi = numpy.stack([mm[0,:], mm[4,:], mm[8,:], mm[1,:], mm[2,:], mm[5,:]], axis=0)
    dder_sc_chi = {}
    return sc_chi, dder_sc_chi


def calc_sc_chi_full(symm_elems, unit_cell_parameters, flag_unit_cell_parameters:bool=False):
    """Calculate symmetry constraint matrix for susceptibility
    It is supposed that determinant of r matrix is 1
    """
    m_m, dder_m_m = calc_m_m_by_unit_cell_parameters(
        unit_cell_parameters, flag_unit_cell_parameters=False)

    r_direct = symm_elems[4:13]
    if symm_elems.shape[0] == 14:
        theta = numpy.expand_dims(symm_elems[13], axis=(0,1))
    else:
        theta = numpy.expand_dims(numpy.ones_like(symm_elems[0]), axis=(0,1))
    r_ccs, dder_r_ccs = calc_m1_m2_inv_m1(m_m, r_direct, flag_m1=False, flag_m2=False)

    res, dder_mm = calc_mm_as_m1_m2_inv_m1(r_ccs, flag_m1= flag_unit_cell_parameters)

    sc_chi_full = (res*theta).sum(axis=2)/res.shape[2]
    dder_sc_chi_full = {}
    return sc_chi_full, dder_sc_chi_full

