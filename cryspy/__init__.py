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
__copyright__   = "Copyright 2020, "
__credits__ = ["Iurii KIBALIN", "Andrew SAZONOV", "Arsen GOUKASSOV"]
__license__ = "GPL"
__version__ = "0.4.20"
__maintainer__ = "Iurii KIBALIN"
__email__ = "iurii.kibalin@cea.fr"
__status__ = "Development"
__date__ = "30.11.2020"
name = "cryspy"

from .A_functions_base.function_1_algebra import \
    calc_scalar_product_by_vectors,\
    calc_scalar_product_by_complex_vectors,\
    calc_modulus_sq_by_complex_vector

from .A_functions_base.function_1_strings import value_error_to_string

from .A_functions_base.function_1_matrices import tri_linear_interpolation

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

from cryspy.C_item_loop_classes.cl_1_cell import Cell, CellL
from cryspy.C_item_loop_classes.cl_1_atom_site import AtomSite, \
    AtomSiteL
from cryspy.C_item_loop_classes.cl_1_atom_type import AtomType, \
    AtomTypeL
from cryspy.C_item_loop_classes.cl_1_atom_site_aniso import \
    AtomSiteAniso, AtomSiteAnisoL
from cryspy.C_item_loop_classes.cl_1_inversed_hessian import InversedHessian
from cryspy.C_item_loop_classes.cl_1_refln import Refln, ReflnL
from cryspy.C_item_loop_classes.cl_1_atom_site_susceptibility import \
    AtomSiteSusceptibility, AtomSiteSusceptibilityL
from cryspy.C_item_loop_classes.cl_1_atom_type_scat import \
    AtomTypeScat, AtomTypeScatL
from cryspy.C_item_loop_classes.cl_1_refln_susceptibility import \
    ReflnSusceptibility, ReflnSusceptibilityL
from cryspy.C_item_loop_classes.cl_1_atom_local_axes import \
    AtomLocalAxes, AtomLocalAxesL
from cryspy.C_item_loop_classes.cl_1_atom_electron_configuration \
    import AtomElectronConfiguration, AtomElectronConfigurationL
from cryspy.C_item_loop_classes.cl_1_chi2 import Chi2, Chi2L
from cryspy.C_item_loop_classes.cl_1_exclude import Exclude, ExcludeL
from cryspy.C_item_loop_classes.cl_1_diffrn_radiation import \
    DiffrnRadiation, DiffrnRadiationL
from cryspy.C_item_loop_classes.cl_1_diffrn_refln import DiffrnRefln, \
    DiffrnReflnL
from cryspy.C_item_loop_classes.cl_1_extinction import Extinction, \
    ExtinctionL
from cryspy.C_item_loop_classes.cl_1_mem_parameters import \
    MEMParameters, MEMParametersL
from cryspy.C_item_loop_classes.cl_1_phase import Phase, PhaseL
from cryspy.C_item_loop_classes.cl_1_pd_background import \
    PdBackground, PdBackgroundL
from cryspy.C_item_loop_classes.cl_1_pd_instr_reflex_asymmetry import \
    PdInstrReflexAsymmetry, PdInstrReflexAsymmetryL
from cryspy.C_item_loop_classes.cl_1_pd_instr_resolution import \
    PdInstrResolution, PdInstrResolutionL
from cryspy.C_item_loop_classes.cl_1_pd_meas import PdMeas, PdMeasL
from cryspy.C_item_loop_classes.cl_1_pd_proc import PdProc, PdProcL
from cryspy.C_item_loop_classes.cl_1_pd_peak import PdPeak, PdPeakL
from cryspy.C_item_loop_classes.cl_1_pd2d_background import \
    Pd2dBackground
from cryspy.C_item_loop_classes.cl_1_pd2d_instr_reflex_asymmetry import\
    Pd2dInstrReflexAsymmetry, Pd2dInstrReflexAsymmetryL
from cryspy.C_item_loop_classes.cl_1_pd2d_instr_resolution import \
    Pd2dInstrResolution, Pd2dInstrResolutionL
from cryspy.C_item_loop_classes.cl_1_pd2d_meas import Pd2dMeas
from cryspy.C_item_loop_classes.cl_1_pd2d_proc import Pd2dProc
from cryspy.C_item_loop_classes.cl_1_pd2d_peak import Pd2dPeak, \
    Pd2dPeakL
from cryspy.C_item_loop_classes.cl_1_range import Range, RangeL
from cryspy.C_item_loop_classes.cl_1_refln import Refln, ReflnL
from cryspy.C_item_loop_classes.cl_1_refln_susceptibility import \
    ReflnSusceptibility, ReflnSusceptibilityL
from cryspy.C_item_loop_classes.cl_1_refine_ls import RefineLs, \
    RefineLsL
from cryspy.C_item_loop_classes.cl_1_setup import Setup, SetupL
from cryspy.C_item_loop_classes.cl_1_texture import Texture, TextureL
from cryspy.C_item_loop_classes.cl_1_tof_background import TOFBackground, \
    TOFBackgroundL
from cryspy.C_item_loop_classes.cl_1_tof_meas import TOFMeas, TOFMeasL
from cryspy.C_item_loop_classes.cl_1_tof_parameters import TOFParameters, \
    TOFParametersL
from cryspy.C_item_loop_classes.cl_1_tof_peak import TOFPeak, TOFPeakL
from cryspy.C_item_loop_classes.cl_1_tof_proc import TOFProc, TOFProcL
from cryspy.C_item_loop_classes.cl_1_tof_profile import TOFProfile, TOFProfileL

from cryspy.C_item_loop_classes.cl_2_atom_site_scat import \
    AtomSiteScat, AtomSiteScatL
from cryspy.C_item_loop_classes.cl_2_diffrn_orient_matrix import \
    DiffrnOrientMatrix, DiffrnOrientMatrixL
from cryspy.C_item_loop_classes.cl_2_section import Section, SectionL
from cryspy.C_item_loop_classes.cl_2_space_group import SpaceGroup
from cryspy.C_item_loop_classes.cl_2_space_group_symop_magn_operation import \
    SpaceGroupSymopMagnOperation, SpaceGroupSymopMagnOperationL

from cryspy.C_item_loop_classes.cl_3_density_point import DensityPointL

from cryspy.D_functions_item_loop.function_1_section_from_density_point \
    import calc_section_from_density_point

from cryspy.E_data_classes.cl_1_crystal import Crystal
from cryspy.E_data_classes.cl_2_diffrn import Diffrn
from cryspy.E_data_classes.cl_2_pd import Pd
from cryspy.E_data_classes.cl_2_pd2d import Pd2d
from cryspy.E_data_classes.cl_2_tof import TOF

from cryspy.G_global_classes.cl_1_rhochi import RhoChi
from cryspy.G_global_classes.cl_1_mem import MEM

from cryspy.H_functions_global.function_1_cryspy_objects import \
    str_to_globaln, file_to_globaln, str_to_items

from cryspy.H_functions_global.function_1_cryspy_objects import \
    L_GLOBAL_CLASS, L_DATA_CLASS, L_LOOP_CLASS, L_ITEM_CLASS

