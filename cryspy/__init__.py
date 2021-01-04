#!/usr/bin/env python
"""
Crystallographic library for polarized neutron diffraction.

Content
-------
    Section A: Base functions
    =========================
        1. Algebra
        ~~~~~~~~~~
            - calc_scalar_product_by_vectors
            - calc_scalar_product_by_complex_vectors
            - calc_modulus_sq_by_complex_vector
        2. Crystallography base
        ~~~~~~~~~~~~~~~~~~~~~~~
            - calc_volume_uc_by_abc_cosines
            - calc_volume_uc_by_abc_angles
            - calc_inverse_d_by_hkl_abc_cosines
            - calc_inverse_d_by_hkl_abc_angles
            - calc_sthovl_by_hkl_abc_cosines
            - calc_sthovl_by_hkl_abc_angles
            - calc_phase_3d
            - calc_moment_2d_by_susceptibility
            - ortogonalize_matrix
        3. Extinction
        ~~~~~~~~~~~~~
            - calc_extinction
            - calc_extinction_2
        4. Flip ratio
        ~~~~~~~~~~~~~
            - calc_f_plus_sq
            - calc_f_minus_sq
            - calc_flip_ratio
       5. Matrices
            - tri_linear_interpolation

    Section B: Parent classes
    =========================
        1. Internal classes (not for user)
            - ItemN, LoopN, DataN, GlobalN
    Section C: Item/Loop classes
    ============================
        1. Classes to describe crystal strucutres
            - unit_cell: Cell, CellL
        2. Classes to describe experiments
    Section D: Functions operating with item, loop
    ==============================================
        1.

    Section E: Data classes
    =======================
        1.

    Section F: Functions operating with Data classes
    ================================================
        1.

    Section G: Global classes
    =========================
        1.

    Section H: Functions operating with Global classes
    ==================================================
        1. str_to_globaln, file_to_globaln, str_to_items
"""
__author__ = 'Iurii KIBALIN'
__copyright__   = "Copyright 2021, "
__credits__ = ["Iurii KIBALIN", "Andrew SAZONOV", "Arsen GOUKASSOV"]
__license__ = "GPL"
__version__ = "0.5.1"
__maintainer__ = "Iurii KIBALIN"
__email__ = "iurii.kibalin@cea.fr"
__status__ = "Development"
__date__ = "04.01.2021"
name = "cryspy"

from .A_functions_base.function_1_algebra import \
    calc_scalar_product_by_vectors,\
    calc_scalar_product_by_complex_vectors,\
    calc_modulus_sq_by_complex_vector

from .A_functions_base.function_1_strings import value_error_to_string

from .A_functions_base.function_1_matrices import tri_linear_interpolation
from .A_functions_base.function_1_markdown import md_to_html

from .A_functions_base.function_2_crystallography_base import \
    calc_volume_uc_by_abc_cosines,\
    calc_volume_uc_by_abc_angles,\
    calc_inverse_d_by_hkl_abc_cosines,\
    calc_inverse_d_by_hkl_abc_angles,\
    calc_sthovl_by_hkl_abc_cosines,\
    calc_sthovl_by_hkl_abc_angles,\
    calc_phase_3d,\
    calc_moment_2d_by_susceptibility,\
    ortogonalize_matrix

from .A_functions_base.function_3_extinction import \
    calc_extinction,\
    calc_extinction_2

from .A_functions_base.function_4_flip_ratio import \
    calc_f_plus_sq, calc_f_minus_sq, calc_flip_ratio

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN
from cryspy.B_parent_classes.cl_3_data import DataN
from cryspy.B_parent_classes.cl_4_global import GlobalN

from cryspy.D_functions_item_loop.function_1_section_from_density_point \
    import calc_section_from_density_point

from cryspy.H_functions_global.function_1_cryspy_objects import \
    str_to_globaln, file_to_globaln, str_to_items

from cryspy.H_functions_global.function_1_cryspy_objects import \
    file_to_globaln as load_file

from cryspy.H_functions_global.function_1_cryspy_objects import \
    L_GLOBAL_CLASS, L_DATA_CLASS, L_LOOP_CLASS, L_ITEM_CLASS, load_packages, \
    add_package, packages, delete_package

def __getattr__(name):
    for obj_class in L_ITEM_CLASS+L_LOOP_CLASS+L_DATA_CLASS+L_GLOBAL_CLASS:
        if name == obj_class.__name__:
            return obj_class
    raise AttributeError(f"module 'cryspy' has no attribute '{name:}'")
    