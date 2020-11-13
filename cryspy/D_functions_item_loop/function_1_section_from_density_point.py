"""Calculation section from density_point.

Functions
---------
    - calc_section_from_density_point
"""
import numpy

from cryspy.A_functions_base.function_1_matrices import \
    tri_linear_interpolation

from cryspy.A_functions_base.function_2_mem import transfer_to_density_3d

from cryspy.C_item_loop_classes.cl_1_atom_site import AtomSiteL
from cryspy.C_item_loop_classes.cl_1_atom_site_susceptibility import \
    AtomSiteSusceptibilityL
from cryspy.C_item_loop_classes.cl_1_cell import Cell
from cryspy.C_item_loop_classes.cl_1_mem_parameters import MEMParameters
from cryspy.C_item_loop_classes.cl_1_space_group_symop import SpaceGroupSymopL

from cryspy.C_item_loop_classes.cl_2_section import Section

from cryspy.C_item_loop_classes.cl_3_density_point import DensityPointL


def calc_section_from_density_point(
        section: Section, density_point: DensityPointL,
        mem_parameters: MEMParameters, cell: Cell,
        space_group_symop: SpaceGroupSymopL, atom_site: AtomSiteL,
        atom_site_susceptibility: AtomSiteSusceptibilityL):
    """Calculate section in density point.

    Arguments
    ---------
        - section
        - density_point
        - cell
        - space_group_symop

    Output:
        - matrix_chi, matrix_b
    """
    r_11 = numpy.array(space_group_symop.r_11, dtype=int)
    r_12 = numpy.array(space_group_symop.r_12, dtype=int)
    r_13 = numpy.array(space_group_symop.r_13, dtype=int)
    r_21 = numpy.array(space_group_symop.r_21, dtype=int)
    r_22 = numpy.array(space_group_symop.r_22, dtype=int)
    r_23 = numpy.array(space_group_symop.r_23, dtype=int)
    r_31 = numpy.array(space_group_symop.r_31, dtype=int)
    r_32 = numpy.array(space_group_symop.r_32, dtype=int)
    r_33 = numpy.array(space_group_symop.r_33, dtype=int)

    b_1 = numpy.array(space_group_symop.b_1, dtype=float)
    b_2 = numpy.array(space_group_symop.b_2, dtype=float)
    b_3 = numpy.array(space_group_symop.b_3, dtype=float)

    r_ij = (r_11, r_12, r_13, r_21, r_22, r_23, r_31, r_32, r_33)
    b_i = (b_1, b_2, b_3)

    np_ind_x = density_point.numpy_index_x
    np_ind_y = density_point.numpy_index_y
    np_ind_z = density_point.numpy_index_z
    np_indexes = (np_ind_x, np_ind_y, np_ind_z)

    den_chi = density_point.calc_density_chi(cell, atom_site_susceptibility)
    den_f = numpy.array(density_point.density_ferro, dtype=float)
    den_a = numpy.array(density_point.density_antiferro, dtype=float)

    points_a = mem_parameters.points_a
    points_b = mem_parameters.points_b
    points_c = mem_parameters.points_c
    chi_f = mem_parameters.chi_ferro
    chi_a = mem_parameters.chi_antiferro

    n_xyz = (points_a, points_b, points_c)
    den_chi_3d = transfer_to_density_3d(np_indexes, den_chi, n_xyz, r_ij, b_i)

    den_b = chi_f*den_f + chi_a*den_a
    den_b_3d = transfer_to_density_3d(np_indexes, den_b, n_xyz, r_ij, b_i)

    fract_x, fract_y, fract_z = section.calc_fractions(cell, atom_site)
    fract_x_uc, fract_y_uc = fract_x % 1., fract_y % 1.
    fract_z_uc = fract_z % 1.

    ind_xyz = numpy.array([fract_x_uc*points_a, fract_y_uc*points_b,
                           fract_z_uc*points_c], dtype=float).transpose()

    den_chi_section = tri_linear_interpolation(ind_xyz, den_chi_3d)
    den_b_section = tri_linear_interpolation(ind_xyz, den_b_3d)

    points_x = section.points_x
    points_y = section.points_y
    den_chi_section = den_chi_section.reshape((points_x, points_y))
    den_b_section = den_b_section.reshape((points_x, points_y))
    return den_chi_section, den_b_section
