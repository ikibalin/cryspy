## Constitution of the library

The library could be separadet on eight categories. 
The functions and objects from each category can 
import objects predefiened earlier.

1. Function, used the standard modules.
2. Common object used as parent classes
3 Item/Loop objects based on the parent classes with defined methods.
4. Functions, which use the Item/Loop objects.
5. Data objects, which includes Item/Loop objects
6. Functions, which use the Data objects
7. Global objects, which includes Item/Loop objects and Data objects
8. Functions, which use the Global objects

Import of the objects is performed consecutively from 1. until 8.


## Naming

The functions are defined in files named as `functions_#level_#tag.py` where `#tag` 
is short description of defined functions, `#level` marks dependences from object defined 
in the same category (it is 1 if there are no dependencies, it is 2 if there is dependence
from functions defined with level 1, it is 3 if there is dependence from functions 
defined with level 2 and so on)

The class for each Item/Loop/Data/Global object is defined in files named as 
`cl_#level_#class_name.py`, where `#class_name.py` is the name of the class.

## Directory

A_functions_base
B_parent_classes
C_item_loop_classes
D_functions_item_loop
E_data_classes
F_functions_data
G_global_classes
H_functions_global

## Function names

The function should be named as `#verb_#output_by_#input`

## Contnt

A_functions_base
    function_1_algebra
    function_1_atomic_vibrations
    function_1_error_simplex
    function_1_gamma_nu
    function_1_magnetic
    function_1_matrices
    function_1_rhocif
    function_1_roots
    function_1_scat_length_neutron
    function_1_strings
    function_2_crystallography_base
    function_2_mem
    function_2_space_group
    function_3_extinction
    function_3_mcif
    function_4_flip_ratio (test)
B_parent_classes
    cl_1_item
    cl_2_loop
    cl_3_data
    cl_4_global
C_item_loop_classes
    cl_1_atom_electron_configuration
    cl_1_atom_local_axes
    cl_1_atom_site_aniso
    cl_1_atom_site_moment
    cl_1_atom_site_susceptibility
    cl_1_atom_site
    cl_1_atom_type_scat
    cl_1_atom_type
    cl_1_cell
    cl_1_chi2
    cl_1_diffrn_radiation
    cl_1_diffrn_refln
    cl_1_exclude
    cl_1_extinction
    cl_1_mem_parameters
    cl_1_pd_background
    cl_1_pd_instr_reflex_asymmetry
    cl_1_pd_instr_resolution
    cl_1_pd_meas
    cl_1_pd_peak
    cl_1_pd_proc
    cl_1_pd2d_background
    cl_1_pd2d_instr_reflex_asymmetry
    cl_1_pd2d_instr_resolution
    cl_1_pd2d_meas
    cl_1_pd2d_peak
    cl_1_pd2d_proc
    cl_1_phase
    cl_1_range
    cl_1_refine_ls
    cl_1_refln_susceptibility
    cl_1_refln
    cl_1_setup
    cl_1_space_group_symop_magn_centering
    cl_1_space_group_symop
    cl_1_space_group_wyckoff
    cl_1_texture
    cl_2_atom_rho_orbital_radial_slater
    cl_2_atom_site_scat
    cl_2_diffrn_orient_matrix
    cl_2_section
    cl_2_space_group_symop_magn_operation
    cl_2_space_group
    cl_3_density_point
D_functions_item_loop
    function_1_calc_for_magcrystal
    function_1_flip_ratio
E_data_classes
    cl_1_crystal
    cl_1_mag_crystal
    cl_1_diffrn
    cl_1_pd
    cl_1_pd2d
F_functions_data
    script_1_mem
H_functions_global
    function_1_cryspy_objects
