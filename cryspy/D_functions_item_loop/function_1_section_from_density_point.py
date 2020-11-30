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
from cryspy.C_item_loop_classes.cl_1_space_group_symop import \
    SpaceGroupSymop, SpaceGroupSymopL

from cryspy.C_item_loop_classes.cl_2_section import Section

from cryspy.C_item_loop_classes.cl_3_density_point import DensityPoint, \
    DensityPointL

def read_file_den(file_name: str):
    """Read density file.
    
    Arguments
    ---------
        -   file_name is a file name
    
    Output
    ------
        - cell
        - space_group_symop
        - density_point
    """
    with open(file_name, "r") as fid:
        l_content = fid.readlines()

    number_lines = int(l_content[1])

    hh = l_content[number_lines+2].strip().split()
    cell = Cell(length_a=float(hh[0]), length_b=float(hh[1]),
                length_c=float(hh[2]), angle_alpha=float(hh[3]),
                angle_beta=float(hh[4]), angle_gamma=float(hh[5]))

    number_points = [int(hh) for hh in 
                     l_content[number_lines+3][:-1].split()[:3]]
    [n_el_symm, centr, n_orig] = [int(hh) for hh in 
                                  l_content[number_lines+3][:-1].split()[3:6]]
    
    l_item_dp = []
    for line in l_content[2:number_lines+2]:
        hh = line.strip().split()
        ind_x, ind_y, ind_z = int(hh[0]), int(hh[1]), int(hh[2])
        den, den_f, den_a = 0., 0., 0.
        if len(hh) == 6:
            den_f = float(hh[4])
            den_a = float(hh[5])
        else:
            den_u = float(hh[3])
            if den_u >= 0.:
                den_f = den_u
            else:
                den_a = den_u
        item = DensityPoint(
            index_x=ind_x, index_y=ind_y, index_z=ind_z,
            fract_x=float(ind_x)/float(number_points[0]),
            fract_y=float(ind_y)/float(number_points[1]),
            fract_z=float(ind_z)/float(number_points[2]),
            density=den, density_ferro=den_f, density_antiferro=den_a)
        l_item_dp.append(item)
    density_point = DensityPointL()
    density_point.items = l_item_dp


    l_el_symm = []
    for line in l_content[(number_lines+4):(number_lines+4+n_el_symm)]:
        hh = line[:-1].strip().split()
        l_el_symm.append([float(hh[9]), int(hh[0]), int(hh[1]), int(hh[2]),
                          float(hh[10]), int(hh[3]), int(hh[4]), int(hh[5]),
                          float(hh[11]), int(hh[6]), int(hh[7]), int(hh[8])])
    l_item_sg = []
    for line in l_content[(number_lines+4+n_el_symm):
                          (number_lines+4+n_el_symm+n_orig)]:
        orig = [float(hh) for hh in line[:-1].split()[:3]]
        for el_symm in l_el_symm:
            el_symm_full = [
                (el_symm[0]+orig[0])%1., el_symm[1], el_symm[2], el_symm[3],
                (el_symm[4]+orig[1])%1., el_symm[5], el_symm[6], el_symm[7],
                (el_symm[8]+orig[2])%1., el_symm[9], el_symm[10], el_symm[11]]
            item = SpaceGroupSymop.define_by_el_symm(el_symm_full)
            l_item_sg.append(item)
    space_group_symop = SpaceGroupSymopL()
    space_group_symop.items = l_item_sg
    return cell, space_group_symop, density_point, number_points

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
