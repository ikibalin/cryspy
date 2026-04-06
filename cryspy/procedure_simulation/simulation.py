import os
from typing import List, Union
import numpy
import scipy
import scipy.optimize

from warnings import warn


from cryspy.A_functions_base.function_1_inversed_hessian import \
    estimate_inversed_hessian_matrix
from cryspy.A_functions_base.function_1_error_simplex import \
    error_estimation_simplex

from cryspy.B_parent_classes.cl_3_data import DataN
from cryspy.B_parent_classes.cl_4_global import GlobalN

from cryspy.C_item_loop_classes.cl_1_phase import Phase, PhaseL
from cryspy.C_item_loop_classes.cl_1_pd_meas import PdMeas, PdMeasL
from cryspy.C_item_loop_classes.cl_1_diffrn_refln import DiffrnRefln, DiffrnReflnL
from cryspy.E_data_classes.cl_1_crystal import Crystal
from cryspy.A_functions_base.unit_cell import calc_m_b_by_unit_cell_parameters
# from cryspy.E_data_classes.cl_1_mag_crystal import MagCrystal
from cryspy.E_data_classes.cl_2_diffrn import Diffrn
from cryspy.E_data_classes.cl_2_pd import Pd
from cryspy.E_data_classes.cl_2_pd2d import Pd2d
from cryspy.E_data_classes.cl_2_tof import TOF

from cryspy.procedure_rhochi.rhochi import rhochi_no_refinement

from cryspy.A_functions_base.structure_factor import \
    calc_index_hkl_multiplicity_in_range
        
import cryspy

na = numpy.newaxis

DIR = os.path.dirname(__file__)


def simulation_paramagnetic_crystal(cryspy_object: cryspy.GlobalN):
    """Add example of paramagnetic crystal to the cryspy object
    """
    f_crystal_para = os.path.join(DIR, "crystal_Ho2Ti2O7.rcif")
    rcif_crystal = cryspy.load_file(f_crystal_para).crystal_ho2Ti2o7
    cryspy_object.add_items([rcif_crystal, ])


def simulation_magnetic_crystal(cryspy_object: cryspy.GlobalN):
    """Add example of magnetic crystal to the cryspy object
    """
    f_crystal_mag = os.path.join(DIR, "crystal_URu0.96Rh0.04Si2.rcif")
    rcif_crystal = cryspy.load_file(f_crystal_mag).crystal_ururhsi
    cryspy_object.add_items([rcif_crystal, ])


def simulation_polarized_neutron_powder_diffraction(cryspy_object: cryspy.GlobalN):
    """Simulate PNPD experiment with default parameters
    """
    f_pd_pnpd = os.path.join(DIR, "pd_pnpd.rcif")
    rcif_pd = cryspy.load_file(f_pd_pnpd).pd_pnpd

    l_crystal = []
    l_name = []
    for item in cryspy_object.items:
        if isinstance(item, cryspy.Crystal):
            l_crystal.append(item)
            l_name.append(item.data_name)

    l_item_phase = []
    for name in l_name:
        item_phase = Phase(label=name, scale=1.)
        l_item_phase.append(item_phase)
    loop_phase = PhaseL()
    loop_phase.items=l_item_phase
    rcif_pd.phase = loop_phase
    cryspy_object.add_items([rcif_pd, ])
    rhochi_no_refinement(cryspy_object)

    np_tth = rcif_pd.pd_proc.numpy_ttheta
    np_plus_net = rcif_pd.pd_proc.numpy_intensity_plus_net
    np_minus_net = rcif_pd.pd_proc.numpy_intensity_minus_net
    np_bkg_calc = rcif_pd.pd_proc.numpy_intensity_bkg_calc
    np_sigma = numpy.ones_like(np_plus_net) * 0.1

    l_item_meas = []
    for tth, plus, minus, bkg_calc, sigma in zip(np_tth, np_plus_net, np_minus_net, np_bkg_calc, np_sigma):
        item_meas = PdMeas(ttheta=tth, intensity_plus=plus + 0.5*bkg_calc, intensity_plus_sigma=sigma,intensity_minus=minus + 0.5*bkg_calc, intensity_minus_sigma=sigma,)
        l_item_meas.append(item_meas)
    loop_meas = PdMeasL()
    loop_meas.items = l_item_meas
    rcif_pd.pd_meas = loop_meas
    rhochi_no_refinement(cryspy_object)
    return     

def simulation_unpolarized_neutron_powder_diffraction(cryspy_object: cryspy.GlobalN):
    """Simulate UNPD experiment with default parameters
    """
    f_pd_unpd = os.path.join(DIR, "pd_unpd.rcif")
    rcif_pd = cryspy.load_file(f_pd_unpd).pd_unpd

    l_crystal = []
    l_name = []
    for item in cryspy_object.items:
        if isinstance(item, cryspy.Crystal):
            l_crystal.append(item)
            l_name.append(item.data_name)

    l_item_phase = []
    for name in l_name:
        item_phase = Phase(label=name, scale=1.)
        l_item_phase.append(item_phase)
    loop_phase = PhaseL()
    loop_phase.items=l_item_phase
    rcif_pd.phase = loop_phase
    cryspy_object.add_items([rcif_pd, ])
    rhochi_no_refinement(cryspy_object)

    np_tth = rcif_pd.pd_proc.numpy_ttheta
    np_plus_net = rcif_pd.pd_proc.numpy_intensity_plus_net
    np_minus_net = rcif_pd.pd_proc.numpy_intensity_minus_net
    np_bkg_calc = rcif_pd.pd_proc.numpy_intensity_bkg_calc
    np_sigma = numpy.ones_like(np_plus_net) * 0.1

    l_item_meas = []
    for tth, plus, minus, bkg_calc, sigma in zip(np_tth, np_plus_net, np_minus_net, np_bkg_calc, np_sigma):
        item_meas = PdMeas(ttheta=tth, intensity=plus + minus + bkg_calc, intensity_sigma=sigma,)
        l_item_meas.append(item_meas)
    loop_meas = PdMeasL()
    loop_meas.items = l_item_meas
    rcif_pd.pd_meas = loop_meas
    rhochi_no_refinement(cryspy_object)
    return  

def simulation_single_crystal(cryspy_object: cryspy.GlobalN):
    """Simulate single crystal experiment with default parameters
    """
    f_diffrn_sc = os.path.join(DIR, "diffrn_sc.rcif")

    l_diffrn = []
    l_name = []
    for item in cryspy_object.items:
        if isinstance(item, cryspy.Crystal):
            d_crystal = item.get_dictionary()

            rcif_diffrn = cryspy.load_file(f_diffrn_sc).diffrn_sc
            rcif_diffrn.phase.label = item.data_name
            wavelength = rcif_diffrn.setup.wavelength

            unit_cell_parameters = d_crystal["unit_cell_parameters"]
            
            m_b = calc_m_b_by_unit_cell_parameters(unit_cell_parameters)[0]
            reduced_symm_elems = d_crystal.get("reduced_symm_elems", numpy.array([0,0,0,1,1,0,0,0,1,0,0,0,1],dtype=int))
            translation_elems = d_crystal.get("translation_elems", numpy.array([0,0,0,1],dtype=int))
            centrosymmetry = d_crystal.get("centrosymmetry", True)
            sthovl_min = 0.01/wavelength
            sthovl_max = 0.85/wavelength
            np_hkl = calc_index_hkl_multiplicity_in_range(
                sthovl_min=sthovl_min, 
                sthovl_max=sthovl_max, unit_cell_parameters=unit_cell_parameters, 
                reduced_symm_elems=reduced_symm_elems, translation_elems=translation_elems, 
                centrosymmetry=centrosymmetry)[0]
            
            rcif_diffrn.diffrn_orient_matrix.ub_11 = float(numpy.round(m_b[0],6))
            rcif_diffrn.diffrn_orient_matrix.ub_12 = float(numpy.round(m_b[1],6))
            rcif_diffrn.diffrn_orient_matrix.ub_13 = float(numpy.round(m_b[2],6))
            rcif_diffrn.diffrn_orient_matrix.ub_21 = float(numpy.round(m_b[3],6))
            rcif_diffrn.diffrn_orient_matrix.ub_22 = float(numpy.round(m_b[4],6))
            rcif_diffrn.diffrn_orient_matrix.ub_23 = float(numpy.round(m_b[5],6))
            rcif_diffrn.diffrn_orient_matrix.ub_31 = float(numpy.round(m_b[6],6))
            rcif_diffrn.diffrn_orient_matrix.ub_32 = float(numpy.round(m_b[7],6))
            rcif_diffrn.diffrn_orient_matrix.ub_33 = float(numpy.round(m_b[8],6))

            l_item_diffrn_refln = []
            for hkl in np_hkl.T:
                item_diffrn_refln = DiffrnRefln(
                    index_h=hkl[0], index_k=hkl[1], index_l=hkl[2], 
                    intensity=1., intensity_sigma=0.01, 
                    )
                l_item_diffrn_refln.append(item_diffrn_refln)
            loop_diffrn_refln = DiffrnReflnL()
            loop_diffrn_refln.items = l_item_diffrn_refln
            rcif_diffrn.diffrn_refln = loop_diffrn_refln
            l_diffrn.append(rcif_diffrn)

    cryspy_object.add_items(l_diffrn)
    rhochi_no_refinement(cryspy_object)

    for rcif_diffrn in l_diffrn:
        for item in rcif_diffrn.diffrn_refln.items:
            item.intensity = float(numpy.round(item.intensity_calc, 5))
    rhochi_no_refinement(cryspy_object)


def simulation_fliping_ratio(cryspy_object: cryspy.GlobalN):
    """Simulate FR experiment with default parameters
    """
    f_diffrn_fr = os.path.join(DIR, "diffrn_fr.rcif")

    l_diffrn = []
    l_name = []
    for item in cryspy_object.items:
        if isinstance(item, cryspy.Crystal):
            d_crystal = item.get_dictionary()

            rcif_diffrn = cryspy.load_file(f_diffrn_fr).diffrn_fr
            rcif_diffrn.phase.label = item.data_name
            wavelength = rcif_diffrn.setup.wavelength

            unit_cell_parameters = d_crystal["unit_cell_parameters"]
            
            m_b = calc_m_b_by_unit_cell_parameters(unit_cell_parameters)[0]

            reduced_symm_elems = d_crystal.get("reduced_symm_elems", numpy.array([0,0,0,1,1,0,0,0,1,0,0,0,1],dtype=int))
            translation_elems = d_crystal.get("translation_elems", numpy.array([0,0,0,1],dtype=int))
            centrosymmetry = d_crystal.get("centrosymmetry", True)
            sthovl_min = 0.01/wavelength
            sthovl_max = 0.85/wavelength
            np_hkl = calc_index_hkl_multiplicity_in_range(
                sthovl_min=sthovl_min, 
                sthovl_max=sthovl_max, unit_cell_parameters=unit_cell_parameters, 
                reduced_symm_elems=reduced_symm_elems, translation_elems=translation_elems, 
                centrosymmetry=centrosymmetry)[0]
            
            rcif_diffrn.diffrn_orient_matrix.ub_11 = float(numpy.round(m_b[0],6))
            rcif_diffrn.diffrn_orient_matrix.ub_12 = float(numpy.round(m_b[1],6))
            rcif_diffrn.diffrn_orient_matrix.ub_13 = float(numpy.round(m_b[2],6))
            rcif_diffrn.diffrn_orient_matrix.ub_21 = float(numpy.round(m_b[3],6))
            rcif_diffrn.diffrn_orient_matrix.ub_22 = float(numpy.round(m_b[4],6))
            rcif_diffrn.diffrn_orient_matrix.ub_23 = float(numpy.round(m_b[5],6))
            rcif_diffrn.diffrn_orient_matrix.ub_31 = float(numpy.round(m_b[6],6))
            rcif_diffrn.diffrn_orient_matrix.ub_32 = float(numpy.round(m_b[7],6))
            rcif_diffrn.diffrn_orient_matrix.ub_33 = float(numpy.round(m_b[8],6))

            l_item_diffrn_refln = []
            for hkl in np_hkl.T:
                item_diffrn_refln = DiffrnRefln(
                    index_h=hkl[0], index_k=hkl[1], index_l=hkl[2], 
                    fr=1., fr_sigma=0.01, 
                    )
                l_item_diffrn_refln.append(item_diffrn_refln)
            loop_diffrn_refln = DiffrnReflnL()
            loop_diffrn_refln.items = l_item_diffrn_refln
            rcif_diffrn.diffrn_refln = loop_diffrn_refln
            l_diffrn.append(rcif_diffrn)

    cryspy_object.add_items(l_diffrn)
    rhochi_no_refinement(cryspy_object)

    for rcif_diffrn in l_diffrn:
        for item in rcif_diffrn.diffrn_refln.items:
            item.fr = float(numpy.round(item.fr_calc, 5))
    rhochi_no_refinement(cryspy_object)

def simulation_current_experiments(cryspy_object: cryspy.GlobalN):
    """Simulate existing experiments
    """
    rhochi_no_refinement(cryspy_object)
    l_diffrn = []
    l_pd = []
    for item in cryspy_object.items:
        if isinstance(item, cryspy.Diffrn):
            l_diffrn.append(item)
        elif isinstance(item, cryspy.Pd):
            l_pd.append(item)

    for rcif_diffrn in l_diffrn:
        for item in rcif_diffrn.diffrn_refln.items:
            if item.is_attribute('fr_calc'):
                item.fr = float(numpy.round(item.fr_calc, 5))
            elif item.is_attribute('intensity_calc'):
                item.intensity = float(numpy.round(item.intensity_calc, 5))


    for rcif_pd in l_pd:
        np_tth = rcif_pd.pd_proc.numpy_ttheta
        np_plus_net = rcif_pd.pd_proc.numpy_intensity_plus_net
        np_minus_net = rcif_pd.pd_proc.numpy_intensity_minus_net
        np_bkg_calc = rcif_pd.pd_proc.numpy_intensity_bkg_calc
        np_sigma = numpy.ones_like(np_plus_net) * 0.1

        flag_polarized = rcif_pd.pd_meas.items[0].is_attribute("intensity_plus") and rcif_pd.pd_meas.items[0].is_attribute("intensity_minus")
        flag_unpolarized = rcif_pd.pd_meas.items[0].is_attribute("intensity")

        if flag_polarized:  
            l_item_meas = []
            for tth, plus, minus, bkg_calc, sigma in zip(np_tth, np_plus_net, np_minus_net, np_bkg_calc, np_sigma):
                item_meas = PdMeas(ttheta=tth, intensity_plus=plus + 0.5*bkg_calc, intensity_plus_sigma=sigma,intensity_minus=minus + 0.5*bkg_calc, intensity_minus_sigma=sigma,)
                l_item_meas.append(item_meas)
        elif flag_unpolarized:              
            l_item_meas = []
            for tth, plus, minus, bkg_calc, sigma in zip(np_tth, np_plus_net, np_minus_net, np_bkg_calc, np_sigma):
                item_meas = PdMeas(ttheta=tth, intensity=plus + minus + bkg_calc, intensity_sigma=sigma,)
                l_item_meas.append(item_meas)
        loop_meas = PdMeasL()
        loop_meas.items = l_item_meas
        rcif_pd.pd_meas = loop_meas
    rhochi_no_refinement(cryspy_object)
    pass
