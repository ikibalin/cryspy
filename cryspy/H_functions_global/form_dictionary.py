import numpy
from cryspy.A_functions_base.magnetic_form_factor import get_j0_j2_parameters

from cryspy.B_parent_classes.cl_3_data import DataN
from cryspy.B_parent_classes.preocedures import take_items_by_class
from cryspy.C_item_loop_classes.cl_1_cell import Cell
from cryspy.C_item_loop_classes.cl_1_atom_site import AtomSiteL
from cryspy.C_item_loop_classes.cl_1_atom_site_aniso import AtomSiteAnisoL
from cryspy.C_item_loop_classes.cl_2_atom_site_scat import AtomSiteScatL
from cryspy.C_item_loop_classes.cl_1_atom_type import AtomTypeL
from cryspy.C_item_loop_classes.cl_1_atom_type_scat import AtomTypeScatL
from cryspy.C_item_loop_classes.cl_1_atom_site_susceptibility import AtomSiteSusceptibilityL
from cryspy.C_item_loop_classes.cl_1_atom_electron_configuration import AtomElectronConfigurationL
from cryspy.C_item_loop_classes.cl_2_space_group import SpaceGroup

from cryspy.C_item_loop_classes.cl_1_setup import Setup
from cryspy.C_item_loop_classes.cl_1_diffrn_refln import DiffrnReflnL
from cryspy.C_item_loop_classes.cl_1_extinction import Extinction
from cryspy.C_item_loop_classes.cl_1_diffrn_radiation import DiffrnRadiation
from cryspy.C_item_loop_classes.cl_2_diffrn_orient_matrix import DiffrnOrientMatrix


from cryspy.C_item_loop_classes.cl_1_refln import ReflnL
from cryspy.C_item_loop_classes.cl_1_refln_susceptibility import ReflnSusceptibilityL


from cryspy.E_data_classes.cl_1_crystal import Crystal
from cryspy.E_data_classes.cl_2_diffrn import Diffrn
from cryspy.E_data_classes.cl_2_pd import Pd



def form_dictionary_by_crystal(data_obj) -> dict:
    """
    Based on cryspy objects form dictionary which contents whole
    information about crystal structure
    """
    if isinstance(data_obj, Crystal):
        ddict = data_obj.get_dictionary()
    else:
        ddict = {}

    return ddict



def form_dictionary_by_diffrn(data_obj) -> dict:
    """
    Based on cryspy objects form dictionary from numpy arrays needed for
    calculations

    Parameters
    ----------
    global_obj : GlobalN of DataN containers
        data container with needed information.

    Returns
    -------
    ddict : dictionary
        output keys: .

    """
    if isinstance(data_obj, Diffrn):
        ddict = data_obj.get_dictionary()
    else:
        ddict = {}

    return ddict



def form_dictionary_by_pd(global_obj) -> dict:
    """
    Based on cryspy objects form dictionary from numpy arrays needed for
    calculations

    Parameters
    ----------
    global_obj : GlobalN of DataN containers
        data container with needed information.

    Returns
    -------
    ddict : dictionary
        output keys: .

    """
    if isinstance(global_obj, Pd):
        global_obj.form_object()

    ddict = {}
    chi2, diffrn_radiation = None, None
    exclude, pd_background = None, None
    pd_instr_reflex_asymmetry, pd_instr_resolution = None, None,
    pd_meas, pd_peak = None, None
    pd_proc, phase = None, None
    range_, refine_ls = None, None
    refln, refln_susceptibility = None, None
    setup = None

    l_obj = take_items_by_class(global_obj, (Setup, ))
    if len(l_obj) > 0:
        setup = l_obj[0]

    l_obj = take_items_by_class(global_obj, (DiffrnRadiation, ))
    if len(l_obj) > 0:
        diffrn_radiation = l_obj[0]

    l_obj = take_items_by_class(global_obj, (ReflnL, ))
    if len(l_obj) > 0:
        refln = l_obj

    l_obj = take_items_by_class(global_obj, (ReflnSusceptibilityL, ))
    if len(l_obj) > 0:
        refln_susceptibility = l_obj


    ddict["name"] = global_obj.get_name()
    if setup is not None:
        ddict["magnetic_field"] = numpy.atleast_1d(setup.field)
        ddict["offset_ttheta"] = numpy.atleast_1d(setup.offset_ttheta)
        ddict["wavelength"] = numpy.atleast_1d(setup.wavelength)

    if diffrn_radiation is not None:
        ddict["beam_polarization"] = numpy.atleast_1d(diffrn_radiation.polarization)
        ddict["flipper_efficiency"] = numpy.atleast_1d(diffrn_radiation.efficiency)

    if refln is not None:
        for refln_phase in refln:
            phase_name = refln_phase.loop_name
            if (refln_phase.is_attribute("index_h") and refln_phase.is_attribute("index_k") and refln_phase.is_attribute("index_l")):
                index_hkl = numpy.array([refln_phase.index_h, refln_phase.index_k, refln_phase.index_l], dtype=int)
                ddict[f"index_hkl_{phase_name:}"] = index_hkl
            if refln_phase.is_attribute("f_calc"):
                f_calc = numpy.array(refln_phase.f_calc, dtype=complex)
                ddict[f"f_nucl_{phase_name:}"] = f_calc
            if (refln_phase.is_attribute("a_calc") and refln_phase.is_attribute("b_calc")):
                a_calc = numpy.array(refln_phase.a_calc, dtype=complex)
                b_calc = numpy.array(refln_phase.b_calc, dtype=complex)
                ddict[f"f_nucl_{phase_name:}"] = a_calc + 1j*b_calc

    if refln_susceptibility is not None:
        for refln_phase in refln_susceptibility:
            phase_name = refln_phase.loop_name
            if (refln_phase.is_attribute("index_h") and refln_phase.is_attribute("index_k") and refln_phase.is_attribute("index_l")):
                index_hkl = numpy.array([refln_phase.index_h, refln_phase.index_k, refln_phase.index_l], dtype=int)
                ddict[f"index_hkl_{phase_name:}"] = index_hkl
            if refln_phase.is_attribute("chi_11_calc"):
                chi_11 = numpy.array(refln_phase.chi_11_calc, dtype=complex)
                chi_12 = numpy.array(refln_phase.chi_12_calc, dtype=complex)
                chi_13 = numpy.array(refln_phase.chi_13_calc, dtype=complex)
                chi_21 = numpy.array(refln_phase.chi_21_calc, dtype=complex)
                chi_22 = numpy.array(refln_phase.chi_22_calc, dtype=complex)
                chi_23 = numpy.array(refln_phase.chi_23_calc, dtype=complex)
                chi_31 = numpy.array(refln_phase.chi_31_calc, dtype=complex)
                chi_32 = numpy.array(refln_phase.chi_32_calc, dtype=complex)
                chi_33 = numpy.array(refln_phase.chi_33_calc, dtype=complex)

                ddict[f"sft_ccs_{phase_name:}"] = numpy.stack([
                    chi_11, chi_12, chi_13, chi_21, chi_22, chi_23, chi_31, chi_32, chi_33], axis=0)

    return ddict

