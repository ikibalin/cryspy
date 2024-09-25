#!/usr/bin/env python
"""
CrysPy is a crystallographic library for neutron data analysis. It allows to
refine polarized neutron diffraction experiments performed with single crystals
as well as with powder magnetic compounds. The program has been mainly
developed for Rietveld analysis of polarized neutron powder diffraction data.
The program can be also used as a MEM tool to reconstruct magnetization
(or spin) density. Time-of-flight (TOF) neutron data analysis is also
available.

List of functions is accessible by command:
    >>> cryspy.functions()

List of classes is accessible by command:
    >>> cryspy.classes()

To search the specific function or the class key words should be used:
    >>> cryspy.functions("flip ratio")
    >>> cryspy.classes("tof")

To load data from ".cif" file use the function:
    >>> data = cryspy.load_file(file_name)

where `file_name` is the address of the file.
    
    
Features
  - Analysis of the polarized neutron scattering on crystals;
  - Diffraction data refinement for single crystals or powder by RhoChi class;
  - MEM for reconstruction of magnetization density by MEM class;
  - Powder refinement of magnetic space groups.
  - TOF neutron data analysis

For more information please see the site: https://www.cryspy.fr/
"""
__author__ = 'Iurii KIBALIN'
__copyright__   = "Copyright 2024, "
__credits__ = ["Iurii KIBALIN", "Andrew SAZONOV", "Arsen GOUKASSOV"]
__license__ = "GPL"
__version__ = "0.7.7"
__maintainer__ = "Iurii KIBALIN"
__email__ = "iurii.kibalin@cea.fr"
__status__ = "Development"
__date__ = "12.09.2024"
name = "cryspy"


from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN
from cryspy.B_parent_classes.cl_3_data import DataN
from cryspy.B_parent_classes.cl_4_global import GlobalN


from cryspy.A_functions_base.function_1_algebra import \
    calc_scalar_product_by_vectors,\
    calc_scalar_product_by_complex_vectors,\
    calc_modulus_sq_by_complex_vector

from cryspy.A_functions_base.function_1_atomic_vibrations import calc_beta_by_u, vibration_constraints, apply_constraint_on_cell_by_type_cell
from cryspy.A_functions_base.function_1_error_simplex import error_estimation_simplex
from cryspy.A_functions_base.function_1_gamma_nu import gammanu_to_tthphi, tthphi_to_gammanu, recal_int_to_tthphi_grid, recal_int_to_gammanu_grid

from cryspy.A_functions_base.function_1_inversed_hessian import estimate_inversed_hessian_matrix
from cryspy.A_functions_base.function_1_magnetic import get_j0_j2_by_symbol
from cryspy.A_functions_base.function_1_markdown import md_to_html

from cryspy.A_functions_base.function_1_matrices import \
    calc_chi_sq, tri_linear_interpolation, transform_string_to_r_b, transform_string_to_digits,\
    transform_fraction_with_label_to_string, transform_digits_to_string, transform_r_b_to_string, \
    calc_mRmCmRT, calc_rotation_matrix_ij_by_euler_angles,\
    calc_euler_angles_by_rotation_matrix_ij, calc_determinant_matrix_ij, \
    calc_inverse_matrix_ij, calc_rotation_matrix_ij_around_axis, \
    calc_product_matrices, calc_product_matrix_vector, calc_vector_angle,\
    calc_vector_product, scalar_product, calc_rotation_matrix_by_two_vectors, \
    ortogonalize_matrix, calc_moment_2d_by_susceptibility, calc_phase_3d
from cryspy.A_functions_base.function_1_objects import get_functions_of_objet, \
    variable_name_to_string, change_variable_name, get_table_html_for_variables
from cryspy.A_functions_base.function_1_roots import calc_roots
from cryspy.A_functions_base.function_1_scat_length_neutron import \
    apply_constraint_on_cell_by_type_cell, get_scat_length_neutron

from cryspy.A_functions_base.function_1_strings import value_error_mark_to_string, \
    ciftext_to_html, find_prefix, common_string, string_to_value_error_mark, \
    transform_string_to_r_b, transform_string_to_digits,\
    transform_fraction_with_label_to_string, transform_digits_to_string,\
    transform_r_b_to_string

from cryspy.A_functions_base.powder_diffraction_tof import tof_Jorgensen, \
    tof_Jorgensen_VonDreele


from cryspy.A_functions_base.function_2_crystallography_base import \
    calc_volume_uc_by_abc_cosines,\
    calc_volume_uc_by_abc_angles,\
    calc_inverse_d_by_hkl_abc_cosines,\
    calc_inverse_d_by_hkl_abc_angles,\
    calc_sthovl_by_hkl_abc_cosines,\
    calc_sthovl_by_hkl_abc_angles,\
    calc_phase_3d,\
    calc_moment_2d_by_susceptibility,\
    ortogonalize_matrix

from cryspy.A_functions_base.function_2_mem import \
    calc_asymmetric_unit_cell_indexes, calc_index_atom_symmetry_closest_to_fract_xyz, \
    calc_factor_in_front_of_density_for_fm, calc_moment_perp,\
    transfer_to_density_3d, transfer_to_chi_3d

from cryspy.A_functions_base.function_2_space_group import \
    get_crystal_system_by_it_number, get_name_hm_short_by_it_number, \
    get_name_schoenflies_by_it_number, get_name_hall_by_it_number, \
    get_name_hm_extended_by_it_number_it_coordinate_system_code, \
    get_name_hm_full_by_it_number

from cryspy.A_functions_base.function_2_sym_elems import \
    form_symm_elems_by_b_i_r_ij, transform_to_p1, get_string_by_symm_elem

from cryspy.A_functions_base.function_3_den_file import \
    read_den_file, save_to_den_file

from cryspy.A_functions_base.function_3_extinction import \
    calc_extinction,\
    calc_extinction_2

from cryspy.A_functions_base.function_3_mcif import \
    calc_fract_by_sym_elem, calc_moment_by_sym_elem

from cryspy.A_functions_base.function_4_flip_ratio import \
    calc_f_plus_sq, calc_f_minus_sq, calc_flip_ratio

from cryspy.B_parent_classes.cl_2_loop import \
 get_prefix_of_loop

from cryspy.D_functions_item_loop.function_1_section_from_density_point \
    import calc_section_from_density_point


from cryspy.D_functions_item_loop.function_1_report_magnetization_ellipsoid import \
    magnetization_ellipsoid_by_u_ij, report_main_axes_of_magnetization_ellipsoids

# from cryspy.D_functions_item_loop.structure_factor import \
#     calculate_nuclear_structure_factor, calculate_structure_factor_tensor

from cryspy.H_functions_global.function_1_cryspy_objects import \
    str_to_globaln, file_to_globaln, str_to_items

from cryspy.H_functions_global.function_1_cryspy_objects import \
    file_to_globaln as load_file

from cryspy.H_functions_global.function_1_cryspy_objects import \
    L_GLOBAL_CLASS, L_DATA_CLASS, L_LOOP_CLASS, L_ITEM_CLASS, load_packages, \
    add_package, packages, delete_package, L_FUNCTION_ADD

from cryspy.H_functions_global.powder_experiments import \
    report_powder_experiments

from cryspy.procedure_rhochi.rhochi import rhochi_rietveld_refinement, \
    rhochi_rietveld_refinement_with_parameters, \
    rhochi_no_refinement, rhochi_inversed_hessian


from cryspy.procedure_mempy.mempy import  mempy_magnetization_density_reconstruction, \
    mempy_spin_density_reconstruction,\
    mempy_reconstruction_with_parameters,\
    mempy_cycle_with_parameters




def repr_function(function, flag_long: bool = False):
    ls_out = []
    name = function.__name__
    code = function.__code__
    defaults_val = function.__defaults__
    if defaults_val is None:
        len_defaults_val = 0
    else:
        len_defaults_val = len(defaults_val)
    argcount = code.co_argcount
    varnames_in = code.co_varnames[:argcount]
    if len_defaults_val == 0:
        s_in_1 = ", ".join([name for name in varnames_in])
    else:
        s_in_1 = ", ".join([name for name in varnames_in[:-len_defaults_val]])
        
    if len_defaults_val == 0:
        s_in_2 = ""
    else:
        s_in_2 = ", ".join([f"{name:}=..." for name in varnames_in[-len_defaults_val:]])

    if ((s_in_1 == "") & (s_in_2 == "")):
        ls_out.append(f"Function: {name:}()")
    elif ((s_in_1 == "") & (s_in_2 != "")):
        ls_out.append(f"Function: {name:}({s_in_2:})")
    elif ((s_in_1 != "") & (s_in_2 == "")):
        ls_out.append(f"Function: {name:}({s_in_1:})")
    else:
        ls_out.append(f"Function: {name:}({s_in_1:}, {s_in_2:})")
    doc = function.__doc__
    if doc is not None:
        if flag_long:
            ls_out.exted([f"    {line:}" for line in doc.split("\n")])
        else:
            
            ls_out.append("    "+doc.strip().split("\n")[0])
    return "\n".join(ls_out)


L_FUNCTION = [
    calc_scalar_product_by_vectors,
    calc_scalar_product_by_complex_vectors,
    calc_modulus_sq_by_complex_vector,
    calc_beta_by_u, vibration_constraints, apply_constraint_on_cell_by_type_cell,
    error_estimation_simplex,
    gammanu_to_tthphi, tthphi_to_gammanu, recal_int_to_tthphi_grid, recal_int_to_gammanu_grid,
    estimate_inversed_hessian_matrix,
    get_j0_j2_by_symbol,
    md_to_html,
    calc_chi_sq, tri_linear_interpolation, transform_string_to_r_b, transform_string_to_digits,
    transform_fraction_with_label_to_string, transform_digits_to_string, transform_r_b_to_string,
    calc_mRmCmRT, calc_rotation_matrix_ij_by_euler_angles,
    calc_euler_angles_by_rotation_matrix_ij, calc_determinant_matrix_ij,
    calc_inverse_matrix_ij, calc_rotation_matrix_ij_around_axis,
    calc_product_matrices, calc_product_matrix_vector, calc_vector_angle,
    calc_vector_product, scalar_product, calc_rotation_matrix_by_two_vectors,
    ortogonalize_matrix, calc_moment_2d_by_susceptibility, calc_phase_3d,
    get_functions_of_objet,
    variable_name_to_string, change_variable_name, get_table_html_for_variables,
    calc_roots,
    apply_constraint_on_cell_by_type_cell, get_scat_length_neutron,
    value_error_mark_to_string, 
    ciftext_to_html, find_prefix, common_string, string_to_value_error_mark,
    transform_string_to_r_b, transform_string_to_digits,
    transform_fraction_with_label_to_string, transform_digits_to_string,
    transform_r_b_to_string,
    tof_Jorgensen,
    tof_Jorgensen_VonDreele,
    calc_volume_uc_by_abc_cosines,
    calc_volume_uc_by_abc_angles,
    calc_inverse_d_by_hkl_abc_cosines,
    calc_inverse_d_by_hkl_abc_angles,
    calc_sthovl_by_hkl_abc_cosines,
    calc_sthovl_by_hkl_abc_angles,
    calc_phase_3d,
    calc_moment_2d_by_susceptibility,
    ortogonalize_matrix,
    calc_asymmetric_unit_cell_indexes, calc_index_atom_symmetry_closest_to_fract_xyz, 
    calc_factor_in_front_of_density_for_fm, calc_moment_perp,
    transfer_to_density_3d, transfer_to_chi_3d,
    get_crystal_system_by_it_number, get_name_hm_short_by_it_number, 
    get_name_schoenflies_by_it_number, get_name_hall_by_it_number, 
    get_name_hm_extended_by_it_number_it_coordinate_system_code, 
    get_name_hm_full_by_it_number,
    form_symm_elems_by_b_i_r_ij, transform_to_p1, get_string_by_symm_elem,
    read_den_file, save_to_den_file,
    calc_extinction,
    calc_extinction_2,
    calc_fract_by_sym_elem, calc_moment_by_sym_elem,
    calc_f_plus_sq, calc_f_minus_sq, calc_flip_ratio,
    calc_section_from_density_point,
    str_to_globaln, file_to_globaln, str_to_items,
    load_file,
    load_packages,
    add_package, packages, delete_package,
    repr_function,
    magnetization_ellipsoid_by_u_ij,
    report_main_axes_of_magnetization_ellipsoids,
    report_powder_experiments,
    rhochi_rietveld_refinement, 
    rhochi_rietveld_refinement_with_parameters, 
    rhochi_no_refinement,
    rhochi_inversed_hessian,
    mempy_magnetization_density_reconstruction, 
    mempy_spin_density_reconstruction,
    mempy_reconstruction_with_parameters,
    mempy_cycle_with_parameters,
    get_prefix_of_loop]

def functions(s_name: str = "", flag_long: bool = False):
    l_function = L_FUNCTION + L_FUNCTION_ADD
    ls_out = []
    for function in l_function:
        cond_1 = function.__name__.lower().find(s_name.lower()) != -1
        if function.__doc__ is None:
            cond_2 = False
        else:
            cond_2 = function.__doc__.lower().find(s_name.lower()) != -1
        if cond_1:
            ls_out.append(repr_function(function, flag_long=flag_long))
        elif cond_2:
            ls_out.append(repr_function(function, flag_long=flag_long))
    s_out = "\n\n".join(ls_out)
    print(s_out)

def repr_class(obj_class):
    s_name = obj_class.__name__
    doc = obj_class.__doc__
    if doc is None:
        s_doc = ""
    else:
        s_doc = doc.strip().split("\n")[0]
    return f"|{s_name:30}|{s_doc:70}|"
    
def classes(s_name: str = ""):
    l_class = L_GLOBAL_CLASS+L_DATA_CLASS+L_LOOP_CLASS+L_ITEM_CLASS
    ls_out = ["|class                         |description                                                           |"]
    ls_out.append("|------------------------------|----------------------------------------------------------------------|")
    for class_ in l_class:
        if class_.__name__.lower().find(s_name.lower()) != -1:
            ls_out.append(repr_class(class_))
    s_out = "\n".join(ls_out)
    print(s_out)


def __getattr__(name):
    for obj_class in L_ITEM_CLASS + L_LOOP_CLASS + L_DATA_CLASS + \
            L_GLOBAL_CLASS + L_FUNCTION_ADD:
        if name == obj_class.__name__:
            return obj_class
    raise AttributeError(f"module 'cryspy' has no attribute '{name:}'")

    