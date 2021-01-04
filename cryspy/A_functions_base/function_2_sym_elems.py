"""Module define operation with functions 
sym_elems  and magn_sym_elems 

sym_elems: numpy.shape[13, n_elems]
[numerator_x, numerator_y, numerator_z, denominator_xyz,
 r_11, r_12, r_13,
 r_21, r_22, r_23,
 r_31, r_32, r_33]

magn_sym_elems: [22, n_elems]
[numerator_x, numerator_y, numerator_z, denominator_xyz,
 r_11, r_12, r_13,
 r_21, r_22, r_23,
 r_31, r_32, r_33,
 m_11, m_12, m_13,
 m_21, m_22, m_23,
 m_31, m_32, m_33]

Functions
---------
    - form_sym_elems_by_b_i_r_ij
    - calc_numerators_denominator_for_b_i
    - calc_common_denominator
    - calc_rational_sum
    - transform_to_p1

"""
import numpy

# from cryspy.A_functions_base.function_1_strings import 


def form_sym_elems_by_b_i_r_ij(b_i, r_ij):
    b_num_1, b_num_2, b_num_3, b_den = calc_numerators_denominator_for_b_i(*b_i)
    (r_11, r_12, r_13, r_21, r_22, r_23, r_31, r_32, r_33) = r_ij

    n_1st = 13
    n_2nd = len(r_11)

    sym_elems = numpy.zeros(shape=(n_1st, n_2nd), dtype=int)

    sym_elems[0, :] = b_num_1
    sym_elems[1, :] = b_num_2
    sym_elems[2, :] = b_num_3
    sym_elems[3, :] = b_den

    sym_elems[4, :] = numpy.array(r_11, dtype=int)
    sym_elems[5, :] = numpy.array(r_12, dtype=int)
    sym_elems[6, :] = numpy.array(r_13, dtype=int)
    sym_elems[7, :] = numpy.array(r_21, dtype=int)
    sym_elems[8, :] = numpy.array(r_22, dtype=int)
    sym_elems[9, :] = numpy.array(r_23, dtype=int)
    sym_elems[10, :] = numpy.array(r_31, dtype=int)
    sym_elems[11, :] = numpy.array(r_32, dtype=int)
    sym_elems[12, :] = numpy.array(r_33, dtype=int)

    return sym_elems


def calc_sym_elems_b_with_common_denominator(*argv):
    """Brings elements of several sym_elems_b to a common denominator

    Parameters
    ----------
    *argv : sym_elems_b [4, n_sym_elems]
        

    Returns
    -------
    normalized sym_elems_b [3, n_sym_elems]
        DESCRIPTION.
    denom_common : int
        common denominator.

    """
    l_norm_sym_elems_b, l_norm_denom = [], []
    for sym_elems_b in argv:
        norm_sym_elems_b, norm_denom = calc_normalized_sym_elems_b(sym_elems_b)
        l_norm_sym_elems_b.append(norm_sym_elems_b)
        l_norm_denom.append(norm_denom)
    denom_common = numpy.lcm.reduce(l_norm_denom)
    t_res = tuple([norm_sym_elems_b * (denom_common // norm_denom) 
                   for norm_sym_elems_b, norm_denom
                   in zip(l_norm_sym_elems_b, l_norm_denom)])
    return *t_res, denom_common
        

def calc_normalized_sym_elems_b(sym_elems_b):
    """Brings elements of sym_elems_b to a common denominator

    Parameters
    ----------
    sym_elems_b : [{"num_x", "num_y", "num_z", "denominator"}, n_symm_elems]
        DESCRIPTION.

    Returns
    -------
    normalized_sym_elems_b : [{"num_x", "num_y", "num_z"}, n_symm_elems],
    common_denominator

    """
    sym_elems_b_denom = sym_elems_b[3, :]
    denom_common = numpy.lcm.reduce(sym_elems_b_denom)
    coeff = denom_common // sym_elems_b_denom
    norm_sym_elems_b = numpy.stack([sym_elems_b[0, :]*coeff,
                                    sym_elems_b[1, :]*coeff,
                                    sym_elems_b[2, :]*coeff], axis=0)
    return norm_sym_elems_b, denom_common


def calc_numerators_denominator_for_b_i(b_1, b_2, b_3):
    # 6*8 = 48
    den_max =  48

    np_b_1 =  numpy.round(den_max * numpy.array(b_1, dtype=float)).astype(int)
    np_b_2 =  numpy.round(den_max * numpy.array(b_2, dtype=float)).astype(int)
    np_b_3 =  numpy.round(den_max * numpy.array(b_3, dtype=float)).astype(int)
    
    np_den_max = den_max * numpy.ones(shape=np_b_1.shape, dtype=int)
    gcd = numpy.gcd.reduce([np_b_1, np_b_2, np_b_3, np_den_max])

    num_x = np_b_1 // gcd
    num_y = np_b_2 // gcd
    num_z = np_b_3 // gcd
    den = den_max // gcd
    return num_x, num_y, num_z, den


def calc_numerators_and_common_denominator(*argv):
    """Calculate common denominator and corresponding numerators
    
    given parameters:
        numerator_1, denominator_1, numerator_2, denominator_2 and so on

    """
    if len(argv)%2 != 0:
        raise ValueError("Wrong number of input arguments.")

    n_numb = len(argv) // 2
    denoms = [argv[2*ind+1] for ind in range(n_numb)]
    den_common = numpy.lcm.reduce(denoms)

    l_res = []
    for ind in range(n_numb):
        num = argv[2*ind]
        denom = argv[2*ind+1]
        coeff = den_common//denom
        l_res.append(num*coeff)
    return *tuple(l_res), den_common


def calc_rational_sum(*argv):
    """Calculate sum of rational numbers given by numerator and denominator
    
    given parameters:
        numerator_1, denominator_1, numerator_2, denominator_2 and so on

    """
    denom_nums = calc_numerators_and_common_denominator(*argv)
    den_common = denom_nums[-1]
    sum_num = numpy.sum(denom_nums[:-1], axis=0)
    gcd = numpy.gcd(sum_num, den_common)
    return sum_num // gcd, den_common // gcd


def transform_to_p1(points_abc, sym_elems, shift, centrosymmetry, indexes_xyz,
                    densities):
    """Transoform indexes_xyz and densities to P1 space group.
    

    Parameters
    ----------
    points_abc : numpy.int [6]
        DESCRIPTION.
    sym_elems : numpy.int [13, n_sym_elems]
        DESCRIPTION.
    shift : numpy.int [4, n_shift]
        DESCRIPTION.
    indexes_xyz : numpy.int [3, n_indexes_xyz]
        DESCRIPTION.
    densities : numpy.int [3, n_densities]
        DESCRIPTION.

    Returns
    -------
    indexes_xyz_p1 : numpy.int [3, n_indexes_xyz_p1]
        DESCRIPTION.
    densities_p1 : numpy.int [3, n_densities_p1]
        DESCRIPTION.

    """
    na = numpy.newaxis

    denom_x = points_abc[0]
    denom_y = points_abc[1]
    denom_z = points_abc[2]

    sym_elems_b_num_1 = sym_elems[0, :]
    sym_elems_b_num_2 = sym_elems[1, :]
    sym_elems_b_num_3 = sym_elems[2, :]
    sym_elems_b_denom = sym_elems[3, :]

    sym_elems_r_11 = sym_elems[4, :]
    sym_elems_r_12 = sym_elems[5, :]
    sym_elems_r_13 = sym_elems[6, :]

    sym_elems_r_21 = sym_elems[7, :]
    sym_elems_r_22 = sym_elems[8, :]
    sym_elems_r_23 = sym_elems[9, :]

    sym_elems_r_31 = sym_elems[10, :]
    sym_elems_r_32 = sym_elems[11, :]
    sym_elems_r_33 = sym_elems[12, :]

    shift_num_1 = shift[0, :]
    shift_num_2 = shift[1, :]
    shift_num_3 = shift[2, :]
    shift_denom = shift[3, :]

    # checking:
    flag_b_x = numpy.alltrue(numpy.mod(denom_x, sym_elems_b_denom) == 0)
    flag_b_y = numpy.alltrue(numpy.mod(denom_y, sym_elems_b_denom) == 0)
    flag_b_z = numpy.alltrue(numpy.mod(denom_z, sym_elems_b_denom) == 0)
    flag_shift_x = numpy.alltrue(numpy.mod(denom_x, shift_denom) == 0)
    flag_shift_y = numpy.alltrue(numpy.mod(denom_y, shift_denom) == 0)
    flag_shift_z = numpy.alltrue(numpy.mod(denom_z, shift_denom) == 0)

    if not(flag_b_x):
        raise ValueError("The number of points along a is not compatible with \
symmetry operations")
    elif not(flag_b_y):
        raise ValueError("The number of points along b is not compatible with \
symmetry operations")
    elif not(flag_b_z):
        raise ValueError("The number of points along c is not compatible with \
symmetry operations")
    elif not(flag_shift_x):
        raise ValueError("The number of points along a is not compatible with \
shifting parameters")
    elif not(flag_shift_y):
        raise ValueError("The number of points along b is not compatible with \
shifting parameters")
    elif not(flag_shift_z):
        raise ValueError("The number of points along c is not compatible with \
shifting parameters")
    
    
    ind_num_x = indexes_xyz[0, :]
    ind_num_y = indexes_xyz[1, :]
    ind_num_z = indexes_xyz[2, :]
    
    denom_xyz = numpy.lcm.reduce([denom_x, denom_y, denom_z])
    num_x = ind_num_x * denom_xyz // denom_x
    num_y = ind_num_y * denom_xyz // denom_y
    num_z = ind_num_z * denom_xyz // denom_z

    denom_new = sym_elems_b_denom * denom_xyz
    
    num_new_x = ((sym_elems_r_11[na, :]*num_x[:, na] +
                  sym_elems_r_12[na, :]*num_y[:, na] +
                  sym_elems_r_13[na, :]*num_z[:, na]) *
                 sym_elems_b_denom + sym_elems_b_num_1[na, :] *
                 denom_xyz)[:, :, na] + \
        (shift_num_1[na, :]*denom_new[:, na]//shift_denom[na, :])[na, :, :]
        
    num_new_y = ((sym_elems_r_21[na, :]*num_x[:, na] +
                  sym_elems_r_22[na, :]*num_y[:, na] +
                  sym_elems_r_23[na, :]*num_z[:, na]) *
                 sym_elems_b_denom + sym_elems_b_num_2[na, :] *
                 denom_xyz)[:, :, na] + \
        (shift_num_2[na, :]*denom_new[:, na]//shift_denom[na, :])[na, :, :]

    num_new_z = ((sym_elems_r_31[na, :]*num_x[:, na] +
                  sym_elems_r_32[na, :]*num_y[:, na] +
                  sym_elems_r_33[na, :]*num_z[:, na]) *
                 sym_elems_b_denom + sym_elems_b_num_3[na, :] *
                 denom_xyz)[:, :, na] + \
        (shift_num_3[na, :]*denom_new[:, na]//shift_denom[na, :])[na, :, :]

    
    ind_new_x = ((num_new_x * denom_x) // denom_new[na, :, na]) % denom_x
    ind_new_y = ((num_new_y * denom_y) // denom_new[na, :, na]) % denom_y
    ind_new_z = ((num_new_z * denom_z) // denom_new[na, :, na]) % denom_z

    if centrosymmetry:
        ind_new_x = numpy.concatenate((ind_new_x, (-1*ind_new_x) % denom_x),
                                      axis=1)
        ind_new_y = numpy.concatenate((ind_new_y, (-1*ind_new_y) % denom_y),
                                      axis=1)
        ind_new_z = numpy.concatenate((ind_new_z, (-1*ind_new_z) % denom_z),
                                      axis=1)

    
    densities_4d_p1 = densities[:, :, na, na]*numpy.ones(shape=(3, )+
                                                     ind_new_x.shape)
    indexes_xyz_4d_p1 = numpy.stack([ind_new_x, ind_new_y, ind_new_z], axis=0)
    
        
    sh_h = densities_4d_p1.shape
    densities_p1 = densities_4d_p1.reshape((sh_h[0], sh_h[1]*sh_h[2]*sh_h[3]))
    
    sh_h = indexes_xyz_4d_p1.shape
    indexes_xyz_p1 = indexes_xyz_4d_p1.reshape((sh_h[0], sh_h[1]*sh_h[2]*
                                                sh_h[3]))
    
    indexes_xyz_p1_uniq, uniq_ind = numpy.unique(
        indexes_xyz_p1, return_index=True, axis=1)
    densities_p1_uniq = densities_p1[:, uniq_ind]

    return indexes_xyz_p1_uniq, densities_p1_uniq
