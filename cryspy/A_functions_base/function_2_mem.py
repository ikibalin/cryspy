"""Common functions from MEM calculations.
"""

from typing import NoReturn
import numpy

from cryspy.A_functions_base.function_1_matrices import \
    calc_product_matrix_vector, calc_vector_product, calc_mRmCmRT

#FIXME delete as it is doubled in symmetry_elements
def calc_asymmetric_unit_cell_indexes(n_x: int, n_y: int, n_z: int, r_ij, b_i)\
       -> NoReturn:
    """
    Calculate indexes of asymmetric unit cell.

    Input parameters:
        - symmetry elements;
        - points number.

    Multiplication of points number on corresponding symmetry element should
    give integer number.

    """
    r_11, r_12, r_13, r_21, r_22, r_23, r_31, r_32, r_33 = r_ij
    b_1, b_2, b_3 = b_i

    ind_x = numpy.array(range(n_x), dtype=int)
    ind_y = numpy.array(range(n_y), dtype=int)
    ind_z = numpy.array(range(n_z), dtype=int)

    ind_2d_x, ind_2d_y, ind_2d_z = numpy.meshgrid(ind_x, ind_y, ind_z,
                                                  indexing="ij")

    ind_x, ind_y = ind_2d_x.flatten(), ind_2d_y.flatten()
    ind_z = ind_2d_z.flatten()

    r_ij_2d = (
        r_11[numpy.newaxis, :], r_12[numpy.newaxis, :], r_13[numpy.newaxis, :],
        r_21[numpy.newaxis, :], r_22[numpy.newaxis, :], r_23[numpy.newaxis, :],
        r_31[numpy.newaxis, :], r_32[numpy.newaxis, :], r_33[numpy.newaxis, :])
    nb_i_2d = (numpy.around(n_x*b_1, 0).astype(int)[numpy.newaxis, :],
               numpy.around(n_y*b_2, 0).astype(int)[numpy.newaxis, :],
               numpy.around(n_z*b_3, 0).astype(int)[numpy.newaxis, :])

    ind_xyz = (ind_x[:, numpy.newaxis], ind_y[:, numpy.newaxis],
               ind_z[:, numpy.newaxis])

    ind_2d_a, ind_2d_b, ind_2d_c = calc_product_matrix_vector(r_ij_2d, ind_xyz)

    ind_2d_a += nb_i_2d[0]
    ind_2d_b += nb_i_2d[1]
    ind_2d_c += nb_i_2d[2]
    ind_2d_a = numpy.mod(ind_2d_a.astype(int), n_x)
    ind_2d_b = numpy.mod(ind_2d_b.astype(int), n_y)
    ind_2d_c = numpy.mod(ind_2d_c.astype(int), n_z)

    ind_2d_abc = n_z*n_y*ind_2d_a + n_z*ind_2d_b + ind_2d_c

    ind_2d_abc_sorted = numpy.sort(ind_2d_abc, axis=1)
    a, ind_a_u_c, counts_a_u_c = numpy.unique(
        ind_2d_abc_sorted[:, 0], return_index=True, return_counts=True)

    ind_x_a_u_c = ind_x[ind_a_u_c]
    ind_y_a_u_c = ind_y[ind_a_u_c]
    ind_z_a_u_c = ind_z[ind_a_u_c]

    return ind_x_a_u_c, ind_y_a_u_c, ind_z_a_u_c, counts_a_u_c


def calc_index_atom_symmetry_closest_to_fract_xyz(
        fract_xyz, fract_atom_xyz, r_ij, b_i, cell):
    """
    Calculate index of atoms and applied symmetry to have closest atoms.

    Basins are defined as closest points to atom

    fract_xyz = [n_points, 3]
    fract_atom_xyz = [n_atoms, 3]
    r_ij = (r_11, r_12, r_13, r_21, r_22, r_23, r_31, r_32, r_33)
    b_i = (b_1, b_2, b_3)
    cell = Cell(length_a, length_b, length_c, angle_alpha, angle_beta,
                angle_gamma)


    Output:
        - n_atom_index = [n_points]: integers from 0 until n_atoms
        - n_symmetry = [n_points]: integers from 0 until n_symmetry
    """
    r_11, r_12, r_13, r_21, r_22, r_23, r_31, r_32, r_33 = r_ij
    b_1, b_2, b_3 = b_i
    atom_x, atom_y, atom_z = fract_atom_xyz[:, 0], fract_atom_xyz[:, 1], \
        fract_atom_xyz[:, 2]
    fract_x, fract_y, fract_z = numpy.mod(fract_xyz[:, 0], 1.), \
        numpy.mod(fract_xyz[:, 1], 1.), numpy.mod(fract_xyz[:, 2], 1.)

    # 2d dimension [atoms, symmetry]
    f_at_x = numpy.mod(
        r_11[numpy.newaxis, :]*atom_x[:, numpy.newaxis] +
        r_12[numpy.newaxis, :]*atom_y[:, numpy.newaxis] +
        r_13[numpy.newaxis, :]*atom_z[:, numpy.newaxis] +
        b_1[numpy.newaxis, :], 1.)
    f_at_y = numpy.mod(
        r_21[numpy.newaxis, :]*atom_x[:, numpy.newaxis] +
        r_22[numpy.newaxis, :]*atom_y[:, numpy.newaxis] +
        r_23[numpy.newaxis, :]*atom_z[:, numpy.newaxis] +
        b_2[numpy.newaxis, :], 1.)
    f_at_z = numpy.mod(
        r_31[numpy.newaxis, :]*atom_x[:, numpy.newaxis] +
        r_32[numpy.newaxis, :]*atom_y[:, numpy.newaxis] +
        r_33[numpy.newaxis, :]*atom_z[:, numpy.newaxis] +
        b_3[numpy.newaxis, :], 1.)

    val_x = fract_x[:, numpy.newaxis, numpy.newaxis]-f_at_x[numpy.newaxis, :,
                                                            :]
    flag_x_05 = numpy.abs(val_x) > 0.5
    flag_x_sign = numpy.where(val_x >= 0.0, -1.0, 1.0)
    val_x_1 = numpy.where(flag_x_05, val_x+flag_x_sign, val_x)

    val_y = fract_y[:, numpy.newaxis, numpy.newaxis]-f_at_y[numpy.newaxis, :,
                                                            :]
    flag_y_05 = numpy.abs(val_y) > 0.5
    flag_y_sign = numpy.where(val_y >= 0.0, -1.0, 1.0)
    val_y_1 = numpy.where(flag_y_05, val_y+flag_y_sign, val_y)

    val_z = fract_z[:, numpy.newaxis, numpy.newaxis]-f_at_z[numpy.newaxis, :,
                                                            :]
    flag_z_05 = numpy.abs(val_z) > 0.5
    flag_z_sign = numpy.where(val_z >= 0.0, -1.0, 1.0)
    val_z_1 = numpy.where(flag_z_05, val_z+flag_z_sign, val_z)

    # 3d dimension [fractions, atoms, symmetry]
    p_x, p_y, p_z = cell.calc_position_by_coordinate(val_x_1, val_y_1, val_z_1)

    dist_sq_3d = numpy.square(p_x) + numpy.square(p_y) + numpy.square(p_z)
    dist_sq_2d_over_sym = numpy.min(dist_sq_3d, axis=2)
    dist_sq_2d_over_at = numpy.min(dist_sq_3d, axis=1)
    
    ind_at = numpy.argmin(dist_sq_2d_over_sym, axis=1)
    ind_sym = numpy.argmin(dist_sq_2d_over_at, axis=1)
    dist_sq = numpy.min(dist_sq_2d_over_sym, axis=1)
    distance = numpy.sqrt(dist_sq)

    return ind_at, ind_sym, distance


def calc_factor_in_front_of_density_for_fm(mult_i, np, volume,
                                           moment_2d, phase_3d, axis=2):
    r"""
    Calculate factor in front of density for magnetic structure factor.

    [hkl, points, symmetry]

    V_i = V_uc/(Ns * Np) * mult_i * moment_2d[i, s] * phase_3d[hkl, i, s]

    F_M = V_i * den_i

    V_uc is volume of unit cell
    Ns is the number of symmetry elements
    Np is the number of points in unit cell
    i is the points in assymetric unit cell
    s is the elemnt of symmetries
    mult_i is the multiplicity of i point
    den_i is the density of i point
    moment_2d[i, s] is the moment in assymetric point i at element symmetry s
        in local coordinate system
    phase_3d[hkl, i, s] is the phase for reflection hkl in the point i for
        symmetry element s


    Output data:
        - v_i_1, v_i_2, v_i_3: [hkl, points]

    """
    # number of symmetry elements
    ns = phase_3d.shape[2]

    # [ind]
    m_rho = (volume*1./float(ns*np))*mult_i

    # [ind, symm]
    m_2d_1, m_2d_2, m_2d_3 = moment_2d
    t_2d_1, t_2d_2, t_2d_3 = m_rho[:, numpy.newaxis] * m_2d_1, \
        m_rho[:, numpy.newaxis] * m_2d_2, m_rho[:, numpy.newaxis] * m_2d_3

    # [hkl, ind, symm]
    v_hkl_3d_1 = t_2d_1[numpy.newaxis, :, :] * phase_3d[:, :, :]
    v_hkl_3d_2 = t_2d_2[numpy.newaxis, :, :] * phase_3d[:, :, :]
    v_hkl_3d_3 = t_2d_3[numpy.newaxis, :, :] * phase_3d[:, :, :]
    v_hkl_2d_1 = v_hkl_3d_1.sum(axis=axis)
    v_hkl_2d_2 = v_hkl_3d_2.sum(axis=axis)
    v_hkl_2d_3 = v_hkl_3d_3.sum(axis=axis)
    return v_hkl_2d_1, v_hkl_2d_2, v_hkl_2d_3


def calc_moment_perp(k_hkl, moment_2d):
    """
    Calculate perpendicular component to k_hkl.

    k_hkl = (k_hkl_x, k_hkl_y, k_hkl_z)                   [hkl]
    moment_2d = (moment_2d_x, moment_2d_y, moment_2d_z)   [hkl, points]
    moment_perp_2d = k_hkl x moment_2d x k_hkl
    """
    k_1, k_2, k_3 = k_hkl
    k_norm = numpy.sqrt(numpy.square(k_1) + numpy.square(k_2) +
                        numpy.square(k_3))
    k_norm[k_norm == 0.] = 1.
    k_1, k_2, k_3 = k_1/k_norm, k_2/k_norm, k_3/k_norm
    # k_hkl_2d = (k_1[:, newaxis], k_2[:, newaxis], k_3[:, newaxis])
    h_1 = calc_vector_product(moment_2d, k_hkl)
    moment_perp = calc_vector_product(k_hkl, h_1)
    return moment_perp


def calc_fm_3d_by_density(mult_i, den_i, np, volume, moment_2d, phase_3d):
    """
    Calculate magnetic structure factor.

    [hkl, points, symmetry]

    F_M = V_uc / (Ns * Np) mult_i den_i moment_2d[i, s] * phase_3d[hkl, i, s]

    V_uc is volume of unit cell
    Ns is the number of symmetry elements
    Np is the number of points in unit cell
    i is the points in assymetric unit cell
    s is the elemnt of symmetries
    mult_i is the multiplicity of i point
    den_i is the density of i point
    moment_2d[i, s] is the moment in assymetric point i at element symmetry s
        in local coordinate system
    phase_3d[hkl, i, s] is the phase for reflection hkl in the point i for
        symmetry element s


    Output data:
        - f_hkl_1d_1, f_hkl_1d_2, f_hkl_1d_3: [hkl]

    """
    # number of symmetry elements
    ns = phase_3d.shape[2]

    # [ind]
    m_rho = (volume*1./float(ns*np))*den_i*mult_i

    m_2d_1, m_2d_2, m_2d_3 = moment_2d

    # [ind, symm]
    t_2d_1, t_2d_2, t_2d_3 = m_rho[:, numpy.newaxis] * \
        m_2d_1, m_rho[:, numpy.newaxis] * m_2d_2, m_rho[:, numpy.newaxis] * \
        m_2d_3

    # [hkl, ind, symm]
    f_hkl_3d_1 = t_2d_1[numpy.newaxis, :, :] * phase_3d[:, :, :]
    f_hkl_3d_2 = t_2d_2[numpy.newaxis, :, :] * phase_3d[:, :, :]
    f_hkl_3d_3 = t_2d_3[numpy.newaxis, :, :] * phase_3d[:, :, :]
    return f_hkl_3d_1, f_hkl_3d_2, f_hkl_3d_3


def calc_fm_2d_by_density(mult_i, den_i, np, volume, moment_2d, phase_3d,
                          axis: int = 2):
    r"""
    Calculate magnetic structure factor.

    [hkl, points] at axis = 2

    F_M = V_uc / (Ns * Np) mult_i den_i \sum_{s}^{Ns} moment_2d[i, s] *
        phase_3d[hkl, i, s]

    [hkl, symmetry] at axis = 1

    F_M = V_uc / (Ns * Np) \sum_{i}^{a.u.c} mult_i den_i  moment_2d[i, s] *
        phase_3d[hkl, i, s]

    """
    f_hkl_3d_1, f_hkl_3d_2, f_hkl_3d_3 = \
        calc_fm_3d_by_density(mult_i, den_i, np, volume, moment_2d, phase_3d)
    f_hkl_2d_1 = f_hkl_3d_1.sum(axis=axis)
    f_hkl_2d_2 = f_hkl_3d_2.sum(axis=axis)
    f_hkl_2d_3 = f_hkl_3d_3.sum(axis=axis)
    return f_hkl_2d_1, f_hkl_2d_2, f_hkl_2d_3


def calc_fm_by_density(mult_i, den_i, np, volume, moment_2d, phase_3d):
    r"""
    Calculate magnetic structure factor.

    [hkl]

    F_M = V_uc / (Ns * Np) \sum_{i}^{a.u.c} mult_i den_i \sum_{s}^{Ns}
        moment_2d[i, s] * phase_3d[hkl, i, s]

    """
    f_hkl_2d_1, f_hkl_2d_2, f_hkl_2d_3 = \
        calc_fm_2d_by_density(mult_i, den_i, np, volume, moment_2d, phase_3d,
                              axis=2)
    f_hkl_1d_1, f_hkl_1d_2, f_hkl_1d_3 = \
        f_hkl_2d_1.sum(axis=1), f_hkl_2d_2.sum(axis=1), f_hkl_2d_3.sum(axis=1)
    return f_hkl_1d_1, f_hkl_1d_2, f_hkl_1d_3


def transfer_to_density_3d(np_indexes, density, n_xyz, r_ij, b_i):
    """Give 3D array of density from asymetric cell.

    Arguments
    ---------
        - np_xyz is numpy array of integer numbers
        - val_1, ... are numpy array of float numbers.
    """
    (n_x, n_y, n_z) = n_xyz
    (i_x, i_y, i_z) = np_indexes
    den_3d = numpy.zeros(n_xyz)
    for r_11, r_12, r_13, r_21, r_22, r_23, r_31, r_32, r_33, b_1, b_2, b_3 \
            in zip(*r_ij, *b_i):
        np_ind_x = numpy.mod((numpy.around((
            i_x*r_11 + i_y*r_12 + i_z*r_13 + n_x*b_1).astype(float), 0)
            ).astype(int), n_x)
        np_ind_y = numpy.mod((numpy.around((
            i_x*r_21 + i_y*r_22 + i_z*r_23 + n_y*b_2).astype(float), 0)
            ).astype(int), n_y)
        np_ind_z = numpy.mod((numpy.around((
            i_x*r_31 + i_y*r_32 + i_z*r_33 + n_z*b_3).astype(float), 0)
            ).astype(int), n_z)
        den_3d[np_ind_x, np_ind_y, np_ind_z] = density
    return den_3d


def transfer_to_chi_3d(np_indexes, chi_11, chi_22, chi_33, chi_12, chi_13,
                       chi_23, n_xyz, r_ij, b_i):
    """
    Give six 3D arrays of susceptibility.

    Input arguments:
        - np_xyz is numpy array of integer numbers;
        - val_1, ... are numpy array of float numbers.
    """
    (n_x, n_y, n_z) = n_xyz
    (i_x, i_y, i_z) = np_indexes
    chi_3d_11, chi_3d_22 = numpy.zeros(n_xyz), numpy.zeros(n_xyz)
    chi_3d_33, chi_3d_12 = numpy.zeros(n_xyz), numpy.zeros(n_xyz)
    chi_3d_13, chi_3d_23 = numpy.zeros(n_xyz), numpy.zeros(n_xyz)

    for r_11, r_12, r_13, r_21, r_22, r_23, r_31, r_32, r_33, b_1, b_2, b_3 \
            in zip(*r_ij, *b_i):
        np_ind_x = numpy.mod((numpy.around((
            i_x*r_11 + i_y*r_12 + i_z*r_13 + n_x*b_1).astype(float), 0)
            ).astype(int), n_x)
        np_ind_y = numpy.mod((numpy.around((
            i_x*r_21 + i_y*r_22 + i_z*r_23 + n_y*b_2).astype(float), 0)
            ).astype(int), n_y)
        np_ind_z = numpy.mod((numpy.around((
            i_x*r_31 + i_y*r_32 + i_z*r_33 + n_z*b_3).astype(float), 0)
            ).astype(int), n_z)
        chi_out = calc_mRmCmRT(
            (r_11, r_12, r_13, r_21, r_22, r_23, r_31, r_32, r_33),
            (chi_11, chi_12, chi_13, chi_12, chi_22, chi_23, chi_13, chi_23,
             chi_33))
        chi_11_rot, chi_12_rot, chi_13_rot, chi_21_rot, chi_22_rot, \
            chi_23_rot, chi_31_rot, chi_32_rot, chi_33_rot = chi_out
        chi_3d_11[np_ind_x, np_ind_y, np_ind_z] = chi_11_rot
        chi_3d_22[np_ind_x, np_ind_y, np_ind_z] = chi_22_rot
        chi_3d_33[np_ind_x, np_ind_y, np_ind_z] = chi_33_rot
        chi_3d_12[np_ind_x, np_ind_y, np_ind_z] = chi_12_rot
        chi_3d_13[np_ind_x, np_ind_y, np_ind_z] = chi_13_rot
        chi_3d_23[np_ind_x, np_ind_y, np_ind_z] = chi_23_rot
    return chi_3d_11, chi_3d_22, chi_3d_33, chi_3d_12, chi_3d_13, chi_3d_23
