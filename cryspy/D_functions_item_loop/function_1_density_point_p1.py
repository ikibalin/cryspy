"""Transformation of density points to P1 space_group

Functions
---------
    -
"""
import numpy

from cryspy.A_functions_base.function_2_mem import transfer_to_density_3d

from cryspy.C_item_loop_classes.cl_1_space_group_symop import \
    SpaceGroupSymopL

from cryspy.C_item_loop_classes.cl_3_density_point import DensityPoint, \
    DensityPointL


def transform_density_point_to_p1(
        density_point : DensityPointL, space_group_symop: SpaceGroupSymopL(),
        points_a: int, points_b: int, points_c: int):
    """Transform density_point to p1.
    
    FIXME: should be deleted
    """
    r_11 = space_group_symop.numpy_r_11.astype(int)
    r_12 = space_group_symop.numpy_r_12.astype(int)
    r_13 = space_group_symop.numpy_r_13.astype(int)

    r_21 = space_group_symop.numpy_r_21.astype(int)
    r_22 = space_group_symop.numpy_r_22.astype(int)
    r_23 = space_group_symop.numpy_r_23.astype(int)

    r_31 = space_group_symop.numpy_r_31.astype(int)
    r_32 = space_group_symop.numpy_r_32.astype(int)
    r_33 = space_group_symop.numpy_r_33.astype(int)

    b_1 = space_group_symop.numpy_b_1.astype(int)
    b_2 = space_group_symop.numpy_b_2.astype(int)
    b_3 = space_group_symop.numpy_b_3.astype(int)

    index_xyz = (numpy.array(density_point.index_x, dtype=int),
                 numpy.array(density_point.index_y, dtype=int),
                 numpy.array(density_point.index_z, dtype=int))

    r_ij = (r_11, r_12, r_13, r_21, r_22, r_23, r_31, r_32, r_33)
    b_i = (b_1, b_2, b_3)

    points_abc = (points_a, points_b, points_c)

    den_i = numpy.array(density_point.density, dtype=float)
    # den_ferro_i = numpy.array(density_point.density_ferro, dtype=float)
    # den_aferro_i = numpy.array(density_point.density_antiferro, dtype=float)

    den_3d = transfer_to_density_3d(index_xyz, den_i, points_abc, r_ij, b_i)

    den_1d = den_3d.flatten()
    p_bc = points_b * points_c
    l_item = []
    for i_den, den in enumerate(den_1d):
        ind_x = i_den // p_bc
        h_y = i_den - p_bc*ind_x
        ind_y = h_y // points_c
        ind_z = h_y % points_c
        item = DensityPoint(
            index_x=ind_x, index_y=ind_y, index_z=ind_z,
            fract_x=float(ind_x)/float(points_a),
            fract_y=float(ind_y)/float(points_b),
            fract_z=float(ind_z)/float(points_c),
            density=den)
        l_item.append(item)

    res_obj = DensityPointL()
    res_obj.items = l_item
    return res_obj
