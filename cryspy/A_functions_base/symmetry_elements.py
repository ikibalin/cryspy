"""
Expressions to work with symmetry elements

symm_elems: numpy.shape[13, n_elems]
[numerator_x, numerator_y, numerator_z, denominator_xyz,
 r_11, r_12, r_13,
 r_21, r_22, r_23,
 r_31, r_32, r_33]

mag_symm_elems: [22, n_elems]
[numerator_x, numerator_y, numerator_z, denominator_xyz,
 r_11, r_12, r_13,
 r_21, r_22, r_23,
 r_31, r_32, r_33,
 m_11, m_12, m_13,
 m_21, m_22, m_23,
 m_31, m_32, m_33]

For details see documentation.
"""
import numpy
from .function_1_strings import transform_string_to_digits
from cryspy.A_functions_base.function_2_sym_elems import calc_numerators_and_common_denominator

from cryspy.A_functions_base.matrix_operations import calc_m1_m2, calc_m_v

na = numpy.newaxis


def calc_symm_flags(symm_elems, atom_symm_elems):
    # cond_1: (R_s - E) R_a = 0
    bn_s_1, bn_s_2, bn_s_3, bd_s = symm_elems[0], symm_elems[1], symm_elems[2], symm_elems[3]
    rse_11 = symm_elems[4] - 1
    rse_12 = symm_elems[5]
    rse_13 = symm_elems[6]
    rse_21 = symm_elems[7]
    rse_22 = symm_elems[8] - 1
    rse_23 = symm_elems[9]
    rse_31 = symm_elems[10]
    rse_32 = symm_elems[11]
    rse_33 = symm_elems[12] - 1

    bn_a_1, bn_a_2, bn_a_3, bd_a = atom_symm_elems[0], atom_symm_elems[1], \
        atom_symm_elems[2], atom_symm_elems[3]
    if (atom_symm_elems.shape[0]==4):
        cond_1_11 = numpy.zeros_like(rse_33*bn_a_1) == 0
        cond_1_12 = cond_1_11
        cond_1_13 = cond_1_11
        cond_1_21 = cond_1_11
        cond_1_22 = cond_1_11
        cond_1_23 = cond_1_11
        cond_1_31 = cond_1_11
        cond_1_32 = cond_1_11
        cond_1_33 = cond_1_11
    else:
        ra_11, ra_12, ra_13, ra_21, ra_22, ra_23, ra_31, ra_32, ra_33 = \
            atom_symm_elems[4], atom_symm_elems[5], atom_symm_elems[6], \
            atom_symm_elems[7], atom_symm_elems[8], atom_symm_elems[9], \
            atom_symm_elems[10], atom_symm_elems[11], atom_symm_elems[12]

        cond_1_11 = rse_11*ra_11 + rse_12*ra_21 + rse_13*ra_31 == 0
        cond_1_12 = rse_11*ra_12 + rse_12*ra_22 + rse_13*ra_32 == 0
        cond_1_13 = rse_11*ra_13 + rse_12*ra_23 + rse_13*ra_33 == 0
        cond_1_21 = rse_21*ra_11 + rse_22*ra_21 + rse_23*ra_31 == 0
        cond_1_22 = rse_21*ra_12 + rse_22*ra_22 + rse_23*ra_32 == 0
        cond_1_23 = rse_21*ra_13 + rse_22*ra_23 + rse_23*ra_33 == 0
        cond_1_31 = rse_31*ra_11 + rse_32*ra_21 + rse_33*ra_31 == 0
        cond_1_32 = rse_31*ra_12 + rse_32*ra_22 + rse_33*ra_32 == 0
        cond_1_33 = rse_31*ra_13 + rse_32*ra_23 + rse_33*ra_33 == 0

    # cond_2: (R_s - E) b_a + b_s = int
    
    bd_sa = bd_s*bd_a

    p1_1 = bd_s * (rse_11*bn_a_1 + rse_12*bn_a_2 + rse_13*bn_a_3) + bd_a * bn_s_1
    p1_2 = bd_s * (rse_21*bn_a_1 + rse_22*bn_a_2 + rse_23*bn_a_3) + bd_a * bn_s_2
    p1_3 = bd_s * (rse_31*bn_a_1 + rse_32*bn_a_2 + rse_33*bn_a_3) + bd_a * bn_s_3
    
    cond_2_1 = numpy.mod(p1_1, bd_sa) == 0
    cond_2_2 = numpy.mod(p1_2, bd_sa) == 0
    cond_2_3 = numpy.mod(p1_3, bd_sa) == 0
    
    la = numpy.logical_and
    symm_flag_point_to_point = \
        la(cond_1_11, la(cond_1_12, la(cond_1_13, la(
           cond_1_21, la(cond_1_22, la(cond_1_23, la(
           cond_1_31, la(cond_1_32, la(cond_1_33, la(
           cond_2_1, la(cond_2_2, cond_2_3)))))))))))
    return symm_flag_point_to_point


def calc_multiplicity_by_atom_symm_elems(full_symm_elems, atom_symm_elems):
    """Calculate the multiplicity for atoms defined by symm_elems and atom_symm_atoms.
    
    atom_symm_atoms can be composed as only from atom position and also from atom position and symmetry elements
    """
    symm_flag_point_to_point = calc_symm_flags(
        full_symm_elems[:, :, na], atom_symm_elems[:, na, :])
    point_multiplicity = numpy.floor_divide(
        full_symm_elems.shape[1], numpy.sum(symm_flag_point_to_point, axis=0))
    return point_multiplicity


def calc_asymmetric_unit_cell_indexes(n_abc, full_symm_elems):
    """
    Calculate indexes of asymmetric unit cell.

    Input parameters:
        - symmetry elements;
        - points number.

    Multiplication of points number on corresponding symmetry element should
    give integer number.

    """
    n_a, n_b, n_c = n_abc[0], n_abc[1], n_abc[2]

    point_index = numpy.stack(numpy.meshgrid(
        numpy.arange(n_a), numpy.arange(n_b), numpy.arange(n_c),
        indexing="ij"), axis=0)
    point_index = point_index.reshape(point_index.shape[0], numpy.prod(point_index.shape[1:]))
    
    elem_r = full_symm_elems[4:13]
    elem_b = full_symm_elems[:4]

    r_ind = calc_m_v(
        numpy.expand_dims(elem_r, axis=1),
        numpy.expand_dims(point_index, axis=2), flag_m=False, flag_v=False)[0]

    div, mod = numpy.divmod(numpy.expand_dims(n_abc, axis=1), numpy.expand_dims(elem_b[3], axis=0))
    if not(numpy.all(mod == 0)):
        raise KeyError("Symmetry elements do not match with number of points")
    point_index_s = numpy.mod(r_ind + numpy.expand_dims(div * elem_b[:3], axis=1),
        numpy.expand_dims(numpy.expand_dims(n_abc, axis=1), axis=2))
    value_index_s = n_c*n_b*point_index_s[0] + n_c*point_index_s[1] + point_index_s[2]
    value_index_s_sorted = numpy.sort(value_index_s, axis=1)

    a, ind_a_u_c, counts_a_u_c = numpy.unique(
        value_index_s_sorted[:, 0], return_index=True, return_counts=True)

    point_index_s_a_u_c = point_index[:, ind_a_u_c]

    return point_index_s_a_u_c, counts_a_u_c


def sum_elem_symm_b(elem_symm_b1, elem_symm_b2):
    num_1, denom_1 = elem_symm_b1[:3], elem_symm_b1[3:4]
    num_2, denom_2 = elem_symm_b2[:3], elem_symm_b2[3:4]
    num = num_1 * denom_2 + num_2 * denom_1
    denom = denom_1*denom_2
    gcd = (numpy.gcd(num, denom)).min(axis=0)
    denom_out = numpy.floor_divide(denom, gcd)
    
    num_out = numpy.floor_divide(num, gcd)
    res = numpy.concatenate([num_out, denom_out], axis=0)
    return res


def calc_full_mag_elems(mag_elems_o, mag_elems_c):
    r_o = mag_elems_o[4:13,:]
    r_c = mag_elems_c[4:13,:]
    b_n_o = mag_elems_o[:3,:]
    b_d_o = mag_elems_o[3:4,:]
    b_c = mag_elems_c[:4,:]
    theta_o = mag_elems_o[13:14,:]
    theta_c = mag_elems_c[13:14,:]
    
    
    theta_fs = theta_o[:, :, na] * theta_c[:, na, :]
    r_fs, dder_r_fs = calc_m1_m2(
        r_o[:, :, na], r_c[:, na, :],
        flag_m1=False, flag_m2=False)
    
    rb, dder_rb = calc_m_v(r_c[:, na, :], b_n_o[:, :, na], flag_m=False, flag_v=False)
    b_d_o = numpy.broadcast_arrays(b_c[0, na, :], b_d_o[:, :, na])[1]
    
    se_rb = numpy.concatenate([rb, b_d_o], axis=0)

    b_fs = sum_elem_symm_b(se_rb, b_c[:, na, :])
    res = numpy.concatenate([b_fs, r_fs, theta_fs], axis=0)
    res_2d = res.reshape((res.shape[0], numpy.prod(res.shape[1:])), order="C")
    mag_elems_fs = res_2d # May be it is not enough and unique elements should be given
    # mag_elems_fs = numpy.unique(res_2d, axis=1)
    return mag_elems_fs


def take_symm_elem_by_string(string):
    """
    Define symm_elems from string

    string: x,z+1/2,-y+3/4
    symm_elems: 0 2 3 4 1 0 0 0 0 1 0 -1 0
    """
    labels=("x", "y", "z")
    l_name = "".join(string.strip().split()).lstrip("(").rstrip(")").split(",")
    rij, bi_nd = [], []
    for _name in l_name:
        coefficients, offset = transform_string_to_digits(_name, labels)
        rij.extend([coefficients[0].numerator, coefficients[1].numerator, coefficients[2].numerator])
        bi_nd.extend([offset.numerator, offset.denominator])
    b_nd = calc_numerators_and_common_denominator(*bi_nd)
    symm_elem = numpy.array(b_nd+tuple(rij), dtype=int)
    return symm_elem


def take_symm_elems_by_string(string):
    l_symm_elems = []
    for s_val in string:
        symm_elem = take_symm_elem_by_string(s_val)
        l_symm_elems.append(symm_elem)
    symm_elems = numpy.stack(l_symm_elems).transpose()
    return symm_elems


def apply_symm_elems_to_index_xyz(symm_elems, index_xyz, points_abc):
    """Apply symmetry elements to indexes."""
    b_n_1, b_n_2, b_n_3 = symm_elems[0], symm_elems[1], symm_elems[2]
    b_d = symm_elems[3]
    r_11, r_12, r_13 = symm_elems[4], symm_elems[5], symm_elems[6]
    r_21, r_22, r_23 = symm_elems[7], symm_elems[8], symm_elems[9]
    r_31, r_32, r_33 = symm_elems[10], symm_elems[11], symm_elems[12]
    i_1, i_2, i_3 = index_xyz[0], index_xyz[1], index_xyz[2]
    
    n1, n2, n3 = points_abc[0], points_abc[1], points_abc[2]
    p_1, p_2, p_3 = n1//b_d, n2//b_d, n3//b_d
    
    ni_1 = numpy.mod(r_11*i_1 + r_12*i_2 + r_13*i_3 + b_n_1*p_1, n1)
    ni_2 = numpy.mod(r_21*i_1 + r_22*i_2 + r_23*i_3 + b_n_2*p_2, n2)
    ni_3 = numpy.mod(r_31*i_1 + r_32*i_2 + r_33*i_3 + b_n_3*p_3, n3)
    ni = numpy.stack([ni_1, ni_2, ni_3], axis=0)
    return ni


def sum_b_elems(b_1, b_2):
    b_1_x, b_1_y, b_1_z, b_1_d = b_1[0], b_1[1], b_1[2], b_1[3]
    b_2_x, b_2_y, b_2_z, b_2_d = b_2[0], b_2[1], b_2[2], b_2[3]
    denom_common = numpy.lcm.reduce([b_1_d, b_2_d])
    c_1 = denom_common // b_1_d
    c_2 = denom_common // b_2_d 
    b_s_x = (b_1_x*c_1 + b_2_x*c_2)%denom_common
    b_s_y = (b_1_y*c_1 + b_2_y*c_2)%denom_common
    b_s_z = (b_1_z*c_1 + b_2_z*c_2)%denom_common
    b_s = numpy.stack([b_s_x, b_s_y, b_s_z, denom_common])
    b_sum = b_s // numpy.expand_dims(numpy.gcd.reduce(b_s), axis=0)
    return b_sum


def calc_full_symm_elems_by_reduced(reduced_symm_elems, centrosymmetry, centrosymmetry_position, translation_elems):
    l_fse = []    
    for te in translation_elems.transpose():
        fse_1 = numpy.copy(reduced_symm_elems)
        te_exp = numpy.ones_like(fse_1[:4,:])*numpy.expand_dims(te, axis=1)
        fse_1[:4,:] = sum_b_elems(fse_1[:4,:], te_exp)
        l_fse.append(fse_1)
    fse_a = numpy.concatenate(l_fse, axis=1)

    if centrosymmetry:
        dcp = numpy.copy(centrosymmetry_position)
        dcp[:3] = 2 * dcp[:3]
        fse_c = numpy.copy(fse_a)
        fse_c[:3] = -1*fse_c[:3]
        fse_c[4:] = -1*fse_c[4:]
        dcp_exp = numpy.ones_like(fse_c[:4,:]) * numpy.expand_dims(dcp, axis=1)
        fse_c[:4,:] = sum_b_elems(dcp_exp, fse_c[:4,:])
        full_symm_elems = numpy.concatenate([fse_a, fse_c], axis=1)
    else:
        full_symm_elems = fse_a
    return full_symm_elems


def calc_equivalent_reflections(index_hkl, reduced_symm_elems):
    r_11, r_12, r_13 = numpy.expand_dims(reduced_symm_elems[4], axis=0), numpy.expand_dims(reduced_symm_elems[5], axis=0), numpy.expand_dims(reduced_symm_elems[6], axis=0)
    r_21, r_22, r_23 = numpy.expand_dims(reduced_symm_elems[7], axis=0), numpy.expand_dims(reduced_symm_elems[8], axis=0), numpy.expand_dims(reduced_symm_elems[9], axis=0)
    r_31, r_32, r_33 = numpy.expand_dims(reduced_symm_elems[10], axis=0), numpy.expand_dims(reduced_symm_elems[11], axis=0), numpy.expand_dims(reduced_symm_elems[12], axis=0)

    h, k, l = numpy.expand_dims(index_hkl[0], axis=1), numpy.expand_dims(index_hkl[1], axis=1), numpy.expand_dims(index_hkl[2], axis=1)

    index_hkl_equivalent = numpy.stack([
        r_11*h + r_21*k + r_31*l,
        r_12*h + r_22*k + r_32*l,
        r_13*h + r_23*k + r_33*l], axis=0)

    index_hkl_equivalent = numpy.concatenate([index_hkl_equivalent, -1*index_hkl_equivalent], axis=2)
    return index_hkl_equivalent
