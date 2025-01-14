"""
Base crystallographic functions.
"""
import numpy
from cryspy.A_functions_base.function_1_matrices import calc_mRmCmRT, \
    calc_product_matrix_vector, scalar_product

from .unit_cell import calc_sthovl_by_unit_cell_parameters


def calc_cos_ang(cell, h_1, k_1, l_1, h_2, k_2, l_2):
    """Calculate directed cosines."""
    q_1_x, q_1_y, q_1_z = cell.calc_k_loc(h_1, k_1, l_1)
    q_2_x, q_2_y, q_2_z = cell.calc_k_loc(h_2, k_2, l_2)
    q_1_sq = q_1_x*q_1_x + q_1_y*q_1_y + q_1_z*q_1_z
    q_2_sq = q_2_x*q_2_x + q_2_y*q_2_y + q_2_z*q_2_z
    q_12 = q_1_x*q_2_x + q_1_y*q_2_y + q_1_z*q_2_z
    res = q_12/(q_1_sq*q_2_sq)**0.5
    res[res > 1.] = 1.
    return res


def calc_volume_uc_by_abc_cosines(a:float, b:float, c:float, 
        cos_alpha:float, cos_beta:float, cos_gamma:float):
    """
    Calculate the volume of unit cell defined as 
    a, b, c, cos(alpha), cos(beta), cos(gamma).
    """
    c_a_sq, c_b_sq = numpy.square(cos_alpha), numpy.square(cos_beta)
    c_g_sq = numpy.square(cos_gamma)
    v_sqrt = numpy.sqrt(1.-c_a_sq-c_b_sq-c_g_sq+2.*cos_alpha*cos_beta*cos_gamma)
    volume_uc = a*b*c*v_sqrt
    return volume_uc


def calc_volume_uc_by_abc_angles(a:float, b:float, c:float, 
                          alpha:float, beta:float, gamma:float):
    """
    Calculate the volume of unit cell defined as 
    a, b, c, alpha, beta, gamma.
    Angles should be given in radians.
    """
    c_a, c_b, c_g = numpy.cos(alpha), numpy.cos(beta), numpy.cos(gamma)
    volume_uc = calc_volume_uc_by_abc_cosines(a, b, c, c_a, c_b, c_g)
    return volume_uc


def calc_inverse_d_by_hkl_abc_cosines(h:float, k:float, l:float, 
        a:float, b:float, c:float, 
        cos_alpha:float, cos_beta:float, cos_gamma:float):
    """
    Calculate 1/d for given reflections h, k, l
    and unit cell parameters defined as a, b, c, cos(alpha), cos(beta), cos(gamma).
    """
    c_a_sq, c_b_sq = numpy.square(cos_alpha), numpy.square(cos_beta)
    c_g_sq = numpy.square(cos_gamma)
    s_a_sq, s_b_sq, s_g_sq = (1.-c_a_sq), (1.-c_b_sq), (1.-c_g_sq)

    A = (1.-c_a_sq-c_b_sq-c_g_sq+2.*cos_alpha*cos_beta*cos_gamma)

    B1 = s_a_sq*numpy.square(h*1./a)+\
         s_b_sq*numpy.square(k*1./b)+\
         s_g_sq*numpy.square(l*1./c)
    B2 = 2.*(k*l*(cos_alpha-cos_beta*cos_gamma))/(b*c)+\
         2.*(h*l*(cos_beta-cos_alpha*cos_gamma))/(a*c)+\
         2.*(h*k*(cos_gamma-cos_alpha*cos_beta))/(a*b)

    # B2 = 2.*(k*l*cos_alpha)/(b*c)+\
    #      2.*(h*l*cos_beta)/(a*c)+\
    #      2.*(h*k*cos_gamma)/(a*b)
    #it should be checked, I am not sure
    B = B1-B2
    inverse_d = numpy.sqrt(B*1./A)
    return inverse_d


def calc_inverse_d_by_hkl_abc_angles(
        h:float, k:float, l:float, a:float, b:float, c:float, alpha:float,
        beta:float, gamma:float):
    """
    Calculate 1/d for given reflections h, k, l
    and unit cell parameters defined as a, b, c, alpha, beta, gamma.
    Angles should be given in radians.
    """
    c_a, c_b, c_g = numpy.cos(alpha), numpy.cos(beta), numpy.cos(gamma)
    inverse_d = calc_inverse_d_by_hkl_abc_cosines(h, k, l, a, b, c, c_a, c_b, c_g)
    return inverse_d


def calc_sthovl_by_hkl_abc_cosines(
        h: float, k: float, l: float, a: float, b: float, c: float,
        cos_alpha: float, cos_beta: float, cos_gamma: float):
    """
    Calculate sin(theta)/lambda for given reflections h, k, l
    and unit cell parameters defined as a, b, c, cos(alpha), cos(beta), cos(gamma).
    """
    inv_d = calc_inverse_d_by_hkl_abc_cosines(h, k, l, a, b, c, 
                                              cos_alpha, cos_beta, cos_gamma)
    sthovl = 0.5*inv_d
    return sthovl


def calc_sthovl_by_hkl_abc_angles(
        h: float, k: float, l: float, a: float, b: float, c: float,
        alpha: float, beta: float, gamma: float):
    """
    Calculate sin(theta)/lambda for given reflections h, k, l
    and unit cell parameters defined as a, b, c, alpha, beta, gamma.
    Angles should be given in radians.
    """
    inv_d = calc_inverse_d_by_hkl_abc_angles(h, k, l, a, b, c, alpha, beta, gamma)
    sthovl = 0.5*inv_d
    return sthovl


def calc_phase_3d(hkl, r_ij, b_i, fract_xyz):
    """
Calculation of phases over 3 dimensions (hkl, points, symmetry):

fract_xyz.shape = (pointxs, 3)

r_ij, b_i are elements of symmetry

phase = exp(2\\pi i (h*Rs*x + h*bs))

phase_3d: [hkl, points, symmetry]
    """
    h, k, l = hkl[0], hkl[1], hkl[2]
    x, y, z = fract_xyz[:, 0], fract_xyz[:, 1], fract_xyz[:, 2]
    r_11, r_12, r_13, r_21, r_22, r_23, r_31, r_32, r_33 = r_ij
    b_1, b_2, b_3 = b_i

    # [ind, symm]
    mr_2d_1, mr_2d_2, mr_2d_3 = calc_product_matrix_vector(
        (r_11[numpy.newaxis, :], r_12[numpy.newaxis, :], r_13[numpy.newaxis, :],
         r_21[numpy.newaxis, :], r_22[numpy.newaxis, :], r_23[numpy.newaxis, :],
         r_31[numpy.newaxis, :], r_32[numpy.newaxis, :], r_33[numpy.newaxis, :]), 
        (x[:, numpy.newaxis], y[:, numpy.newaxis], z[:, numpy.newaxis]))

    # [hkl, ind, symm]
    hrr_3d = scalar_product((h[:, numpy.newaxis, numpy.newaxis], 
                                     k[:, numpy.newaxis, numpy.newaxis], 
                                     l[:, numpy.newaxis, numpy.newaxis]),
        (mr_2d_1[numpy.newaxis, :, :], mr_2d_2[numpy.newaxis, :, :], mr_2d_3[numpy.newaxis, :, :]))

    # [hkl, symm]
    hb_2d = scalar_product((h[:, numpy.newaxis], k[:, numpy.newaxis], l[:, numpy.newaxis]), 
                           (b_1[numpy.newaxis, :], b_2[numpy.newaxis, :], b_3[numpy.newaxis, :]))
    phase_3d = numpy.exp(2.*numpy.pi*1j*(hrr_3d+hb_2d[:, numpy.newaxis,:]))
    return phase_3d


def calc_phase_by_hkl_xyz_rb(index_hkl, x, y, z, r_11, r_12, r_13, r_21, r_22,
                             r_23, r_31, r_32, r_33, b_1, b_2, b_3):
    h, k, l = index_hkl[0], index_hkl[1], index_hkl[2]
    np_h, np_x, np_r_11 = numpy.meshgrid(h, x, r_11, indexing="ij")
    np_k, np_y, np_r_22 = numpy.meshgrid(k, y, r_22, indexing="ij")
    np_l, np_z, np_r_33 = numpy.meshgrid(l, z, r_33, indexing="ij")

    np_r_12 = numpy.meshgrid(h, x, r_12, indexing="ij")[2]
    np_r_13 = numpy.meshgrid(k, y, r_13, indexing="ij")[2]
    np_r_23 = numpy.meshgrid(l, z, r_23, indexing="ij")[2]
    np_r_21 = numpy.meshgrid(h, x, r_21, indexing="ij")[2]
    np_r_31 = numpy.meshgrid(k, y, r_31, indexing="ij")[2]
    np_r_32 = numpy.meshgrid(l, z, r_32, indexing="ij")[2]
    np_b_1 = numpy.meshgrid(l, z, b_1, indexing="ij")[2]
    np_b_2 = numpy.meshgrid(l, z, b_2, indexing="ij")[2]
    np_b_3 = numpy.meshgrid(l, z, b_3, indexing="ij")[2]

    np_x_s = np_x*np_r_11 + np_y*np_r_12 + np_z*np_r_13 + np_b_1
    np_y_s = np_x*np_r_21 + np_y*np_r_22 + np_z*np_r_23 + np_b_2
    np_z_s = np_x*np_r_31 + np_y*np_r_32 + np_z*np_r_33 + np_b_3
    hh = (2*numpy.pi*1j*(np_h*np_x_s + np_k*np_y_s+ np_l*np_z_s)).astype(
        complex)
    phase_3d = numpy.exp(hh)
    return phase_3d


def calc_power_dwf_iso(b_iso, sthovl):
    """
isotropic harmonic Debye-Waller factor
    """
    sthovl_sq = sthovl**2
    b_iso_2d, sthovl_sq_2d = numpy.meshgrid(sthovl_sq, b_iso, indexing="ij")
    power_dwf_iso_2d = b_iso_2d*sthovl_sq_2d
    return power_dwf_iso_2d


def calc_power_dwf_aniso(index_hkl, beta, 
    r_11, r_12, r_13, r_21, r_22, r_23, r_31, r_32, r_33):
    """
anisotropic harmonic Debye-Waller factor

h,k,l is 1D (temporary solution)
    """
    h, k, l = index_hkl[0], index_hkl[1], index_hkl[2]
    b_11, b_22, b_33 = beta[:, 0], beta[:, 1], beta[:, 2]
    b_12, b_13, b_23 = beta[:, 3], beta[:, 4], beta[:, 5]

    np_h, np_b_11, np_r_11 = numpy.meshgrid(h, b_11, r_11, indexing="ij")
    np_k, np_b_22, np_r_22 = numpy.meshgrid(k, b_22, r_22, indexing="ij")
    np_l, np_b_33, np_r_33 = numpy.meshgrid(l, b_33, r_33, indexing="ij")
    np_h, np_b_12, np_r_12 = numpy.meshgrid(h, b_12, r_12, indexing="ij")
    np_h, np_b_13, np_r_13 = numpy.meshgrid(h, b_13, r_13, indexing="ij")
    np_h, np_b_23, np_r_23 = numpy.meshgrid(h, b_23, r_23, indexing="ij")
    np_r_21 = numpy.meshgrid(h, b_23, r_21, indexing="ij")[2]
    np_r_31 = numpy.meshgrid(h, b_23, r_31, indexing="ij")[2]
    np_r_32 = numpy.meshgrid(h, b_23, r_32, indexing="ij")[2]

    np_h_s = np_h*np_r_11 + np_k*np_r_21 + np_l*np_r_31
    np_k_s = np_h*np_r_12 + np_k*np_r_22 + np_l*np_r_32
    np_l_s = np_h*np_r_13 + np_k*np_r_23 + np_l*np_r_33

    power_dwf_aniso = (np_b_11*np_h_s**2 + np_b_22*np_k_s**2 +
                       np_b_33*np_l_s**2 + 2.*np_b_12*np_h_s*np_k_s +
                       2.*np_b_13*np_h_s*np_l_s + 2.*np_b_23*np_k_s*np_l_s)
    return power_dwf_aniso


def calc_dwf(cell, index_hkl, b_iso, beta, r_11, r_12, r_13, r_21, r_22, r_23,
             r_31, r_32, r_33):
    """Calculate Debye-Waller factor."""
    unit_cell_parameters = cell.get_unit_cell_parameters()
    sthovl = calc_sthovl_by_unit_cell_parameters(index_hkl, unit_cell_parameters, flag_unit_cell_parameters=False)[0]
    # dimensions (hkl, atoms in assymmetric unit cell)
    power_iso_2d = calc_power_dwf_iso(b_iso, sthovl)

    # dimensions (hkl, atoms in assymmetric unit cell, el.symmetry)
    power_aniso_3d = calc_power_dwf_aniso(
        index_hkl, beta, r_11, r_12, r_13, r_21, r_22, r_23, r_31, r_32, r_33)
    power_3d = power_iso_2d[:, :, numpy.newaxis] + power_aniso_3d
    dwf_3d = numpy.exp(-power_3d)
    return dwf_3d


def calc_form_factor_tensor_susceptibility(
        chi_11, chi_22, chi_33, chi_12, chi_13, chi_23, space_group_symop,
        form_factor, cell, h, k, l):
    """Give components of form factor tensor for susceptibility.

    fft_11, fft_12, fft_13
    fft_21, fft_22, fft_23
    fft_31, fft_32, fft_33

    in 3 dimension (hkl, atoms, symmetry elements).
    """
    ff = numpy.array(form_factor, dtype=float)
    sthovl = cell.calc_sthovl(h, k, l)
    # dimension (hkl, atoms)
    
    r_11 = numpy.array(space_group_symop.r_11, dtype=float)
    r_12 = numpy.array(space_group_symop.r_12, dtype=float)
    r_13 = numpy.array(space_group_symop.r_13, dtype=float)
    r_21 = numpy.array(space_group_symop.r_21, dtype=float)
    r_22 = numpy.array(space_group_symop.r_22, dtype=float)
    r_23 = numpy.array(space_group_symop.r_23, dtype=float)
    r_31 = numpy.array(space_group_symop.r_31, dtype=float)
    r_32 = numpy.array(space_group_symop.r_32, dtype=float)
    r_33 = numpy.array(space_group_symop.r_33, dtype=float)

    chi_21, chi_31, chi_32 = chi_12, chi_13, chi_23

    c11, r11 = numpy.meshgrid(chi_11, r_11, indexing="ij")
    c22, r22 = numpy.meshgrid(chi_22, r_22, indexing="ij")
    c33, r33 = numpy.meshgrid(chi_33, r_33, indexing="ij")
    c12, r12 = numpy.meshgrid(chi_12, r_12, indexing="ij")
    c13, r13 = numpy.meshgrid(chi_13, r_13, indexing="ij")
    c23, r23 = numpy.meshgrid(chi_23, r_23, indexing="ij")
    c21, r21 = numpy.meshgrid(chi_21, r_21, indexing="ij")
    c31, r31 = numpy.meshgrid(chi_31, r_31, indexing="ij")
    c32, r32 = numpy.meshgrid(chi_32, r_32, indexing="ij")

    rcrt_11, rcrt_12, rcrt_13, rcrt_21, rcrt_22, rcrt_23, rcrt_31, rcrt_32, \
        rcrt_33 = calc_mRmCmRT(
            (r11, r12, r13, r21, r22, r23, r31, r32, r33),
            (c11, c12, c13, c21, c22, c23, c31, c32, c33))

    # dimension (hkl, atoms, symmetry)
    n_a = numpy.newaxis
    fft_11 = ff[:, :, n_a] * rcrt_11[n_a, :, :]
    fft_12 = ff[:, :, n_a] * rcrt_12[n_a, :, :]
    fft_13 = ff[:, :, n_a] * rcrt_13[n_a, :, :]
    fft_21 = ff[:, :, n_a] * rcrt_21[n_a, :, :]
    fft_22 = ff[:, :, n_a] * rcrt_22[n_a, :, :]
    fft_23 = ff[:, :, n_a] * rcrt_23[n_a, :, :]
    fft_31 = ff[:, :, n_a] * rcrt_31[n_a, :, :]
    fft_32 = ff[:, :, n_a] * rcrt_32[n_a, :, :]
    fft_33 = ff[:, :, n_a] * rcrt_33[n_a, :, :]
    # ortogonalization should be done
    return fft_11, fft_12, fft_13, fft_21, fft_22, fft_23, fft_31, fft_32, \
        fft_33

def calc_moment_2d_by_susceptibility(r_ij, susc_i, m_norm_ij, h_loc):
    """Recalculate chi_i given according to symmetry elements for each point.

    chi_i  is given in reciprocal unit cell.

    After susceptibility is multiplied in magnetic field defined 

    r_ij:= r_11, r_12, r_13, r_21, r_22, r_23, r_31, r_32, r_33 
    susc_i:= chi_11, chi_22, chi_33, chi_12, chi_13, chi_23

    Matrix m_norm  used to recalculate coordinates  from direct space
    (a1/|a1|, a2/|a2|, a3/|a3|) to Cartesian one (x||a*, z||c).

    x_cart = m_norm * x_direct

    m_norm  = [[(1 - cos**2 alpha1 - cos**2 alpha2 - cos**2 alpha3 + \
                 2 cos alpha1 cos alpha2 cos alpha3)**0.5/sin(alpha1),  0,  0],
           [(cos alpha3 - cos alpha1 cos alpha2) / sin alpha1, sin alpha1,  0],
           [cos alpha2, cos alpha1,  1]]

    Matrix m_norm should be given as

        m_norm_ij = (_11, _12, _13, _21, _22, _23, _31, _32, _33) 

    Output
    ------
        moment_2d: [points, symmetry]
    """
    n_a = numpy.newaxis
    r_11, r_12, r_13, r_21, r_22, r_23, r_31, r_32, r_33 = r_ij
    chi_11, chi_22, chi_33, chi_12, chi_13, chi_23 = susc_i
    # [ind, symm]
    chi_2d_ij = calc_mRmCmRT(
        (r_11[n_a, :], r_12[n_a, :], r_13[n_a, :], r_21[n_a, :], r_22[n_a, :],
         r_23[n_a, :], r_31[n_a, :], r_32[n_a, :], r_33[n_a, :]),
        (chi_11[:, n_a], chi_12[:, n_a], chi_13[:, n_a], chi_12[:, n_a],
         chi_22[:, n_a], chi_23[:, n_a], chi_13[:, n_a], chi_23[:, n_a],
         chi_33[:, n_a]))
    chi_orto_ij = ortogonalize_matrix(chi_2d_ij, m_norm_ij)
    moment_2d = calc_product_matrix_vector(chi_orto_ij, h_loc)
    return moment_2d


def ortogonalize_matrix(m_ij, m_norm_ij):
    """Ortogonalize matrix.

    matrix m_ij is defined in coordinate system (a, b, c).
    It is given as tuple

    Matrix m_norm  used to recalculate coordinates  from direct space
    (a1/|a1|, a2/|a2|, a3/|a3|) to Cartesian one (x||a*, z||c).

    x_cart = m_norm * x_direct

    m_norm  = [[(1 - cos**2 alpha1 - cos**2 alpha2 - cos**2 alpha3 +  \
                 2 cos alpha1 cos alpha2 cos alpha3)**0.5/sin(alpha1),  0,  0],
           [(cos alpha3 - cos alpha1 cos alpha2) / sin alpha1, sin alpha1,  0],
           [cos alpha2, cos alpha1,  1]]

    Should be given as

        m_norm_ij = (_11, _12, _13, _21, _22, _23, _31, _32, _33)

    output matrix s_ij is defined in Cartezian coordinate system defined as
    x||a*, z||c, y= [z x] (right handed)
    """
    s_11, s_12, s_13, s_21, s_22, s_23, s_31, s_32, s_33 = calc_mRmCmRT(
        m_norm_ij, m_ij)
    return s_11, s_12, s_13, s_21, s_22, s_23, s_31, s_32, s_33


def calc_atoms_in_unit_cell(r_ij, b_i, fract_xyz, atom_label):
    """Calculate atoms in unit cell.

    Arguments
    ---------
        - r_ij = (r_11, r_12, r_13, r_21, r_22, r_23, r_31, r_32, r_33)
        - b_i = (b_1, b_2, b_3)
        - fract_xyz = (fract_x, fract_y, fract_z)
        - atom_label

    Output
    ------
        - fract_uc_x
        - fract_uc_y
        - fract_uc_z
        - label_uc
    """
    (r_11, r_12, r_13, r_21, r_22, r_23, r_31, r_32, r_33) = r_ij
    (b_1, b_2, b_3) = b_i
    atom_number = numpy.linspace(0, atom_label.size-1, atom_label.size)
    (fract_x, fract_y, fract_z) = fract_xyz
    na = numpy.newaxis
    f_x = numpy.mod(r_11[:, na]*fract_x[na, :] + r_12[:, na]*fract_y[na, :] +
                    r_13[:, na]*fract_z[na, :] + b_1[:, na], 1.)
    f_y = numpy.mod(r_21[:, na]*fract_x[na, :] + r_22[:, na]*fract_y[na, :] +
                    r_23[:, na]*fract_z[na, :] + b_2[:, na], 1.)
    f_z = numpy.mod(r_31[:, na]*fract_x[na, :] + r_32[:, na]*fract_y[na, :] +
                    r_33[:, na]*fract_z[na, :] + b_3[:, na], 1.)
    r_f_x = numpy.round(f_x, decimals=5)
    r_f_y = numpy.round(f_y, decimals=5)
    r_f_z = numpy.round(f_z, decimals=5)

    a_n = numpy.ones(r_f_x.shape, dtype=float) * atom_number[na, :]
    r_f_xyz = numpy.array([r_f_x.flatten(), r_f_y.flatten(), r_f_z.flatten(),
                           a_n.flatten()],
                          dtype=float)
    u_f_xyz = numpy.unique(r_f_xyz, axis=1)

    fract_uc_x, fract_uc_y, fract_uc_z = u_f_xyz[0], u_f_xyz[1], u_f_xyz[2]
    number_uc = numpy.round(u_f_xyz[3], decimals=0).astype(int)
    label_uc = atom_label[number_uc]
    ind_sort = numpy.argsort(number_uc)

    fract_uc_x_s = numpy.take_along_axis(fract_uc_x, ind_sort, 0)
    fract_uc_y_s = numpy.take_along_axis(fract_uc_y, ind_sort, 0)
    fract_uc_z_s = numpy.take_along_axis(fract_uc_z, ind_sort, 0)
    label_uc_s = numpy.take_along_axis(label_uc, ind_sort, 0)

    return fract_uc_x_s, fract_uc_y_s, fract_uc_z_s, label_uc_s
