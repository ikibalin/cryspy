#!/usr/bin/env python
"""
CrysPy is a crystallographic library for neutron data analysis. It allows to
refine polarized neutron diffraction experiments performed with single crystals
as well as with powder magnetic compounds. The program has been mainly
developed for Rietveld analysis of polarized neutron powder diffraction data.
The program can be also used as a MEM tool to reconstruct magnetization
(or spin) density. Time-of-flight (TOF) neutron data analysis is also
available.

Features
  - Analysis of the polarized neutron scattering on crystals by the library
    CrysPy;
  - Diffraction data refinement for single crystals or powder by RhoChi;
  - MEM for reconstruction of magnetization density;
  - Powder refinement of magnetic space groups (under development).
  - TOF neutron data analysis (under development)

For more information please see the site: https://ikibalin.github.io/cryspy/
"""
__author__ = 'Iurii KIBALIN'
__copyright__   = "Copyright 2021, "
__credits__ = ["Iurii KIBALIN", "Andrew SAZONOV", "Arsen GOUKASSOV"]
__license__ = "GPL"
__version__ = "0.5.5"
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
    