import numpy

from cryspy.B_parent_classes.cl_4_global import GlobalN

from cryspy.C_item_loop_classes.cl_1_chi2 import Chi2
from cryspy.C_item_loop_classes.cl_1_pd_background import PdBackgroundL
from cryspy.C_item_loop_classes.cl_1_pd_instr_reflex_asymmetry import PdInstrReflexAsymmetry
from cryspy.C_item_loop_classes.cl_1_pd_instr_resolution import PdInstrResolution
from cryspy.C_item_loop_classes.cl_1_pd_meas import PdMeasL
from cryspy.C_item_loop_classes.cl_1_diffrn_radiation import DiffrnRadiation
from cryspy.C_item_loop_classes.cl_1_range import Range
from cryspy.C_item_loop_classes.cl_1_phase import PhaseL
from cryspy.C_item_loop_classes.cl_1_setup import Setup

from cryspy.E_data_classes.cl_1_crystal import Crystal
from cryspy.E_data_classes.cl_2_pd import Pd

from cryspy.procedure_rhochi.rhochi import rhochi_no_refinement

def simulation_pd(l_crystal:list,
                  ttheta_start:float=4.,
                  ttheta_end:float=120.,
                  ttheta_step:float=0.2,
                  flag_polarized:bool=False):
    l_crystal = [item for item in l_crystal if isinstance(item, Crystal)]
    n_step = int(numpy.round((ttheta_end-ttheta_start)/ttheta_step+1.))
    tth_deg = numpy.round(numpy.linspace(ttheta_start, ttheta_end, n_step), 3)
    bckg_pd = PdBackgroundL(
        intensity=[0.,0.],
        ttheta=[ttheta_start, ttheta_end])
    instr_res_pd = PdInstrResolution(w=0.1)
    phase_pd = PhaseL(scale=[1. for hh in l_crystal], label=[hh.data_name for hh in l_crystal])
    setup_pd = Setup(wavelength=1.4, offset_ttheta=0.)
    instr_asym_pd = PdInstrReflexAsymmetry(p1=0, p2=0, p3=0, p4=0)
    range_pd = Range(ttheta_min=ttheta_start, ttheta_max=ttheta_end)
    chi_sq = Chi2(sum=True, diff=False, up=False, down=False)
    l_tem_pd = [bckg_pd, instr_res_pd, phase_pd, setup_pd, instr_asym_pd, range_pd, chi_sq]
    if flag_polarized:
        meas_pd = PdMeasL(ttheta=tth_deg,
                                 intensity_plus=0*tth_deg,
                                 intensity_plus_sigma=0*tth_deg+0.1,
                                 intensity_minus=0*tth_deg,
                                 intensity_minus_sigma=0*tth_deg+0.1,)
        diffrn_rad = DiffrnRadiation(polarization=1., efficiency=1.)
        setup_pd.field = 1.
        l_tem_pd.extend([meas_pd, diffrn_rad])
    else:
        meas_pd = PdMeasL(ttheta=tth_deg, intensity=0*tth_deg, intensity_sigma=0*tth_deg+0.1)
        l_tem_pd.append(meas_pd)

    pd = Pd(data_name="sim_1d")
    pd.add_items(l_tem_pd)

    l_item_global = [pd, ] + l_crystal
    classes = set([type(item) for item in l_item_global])
    global_obj_2 = GlobalN.make_container((), tuple(classes), "global")
    global_obj_2.global_name = ""
    global_obj_2.add_items(l_item_global)

    rhochi_no_refinement(global_obj_2)

    pd_proc = pd.pd_proc
    pd_meas = pd.pd_meas
    if flag_polarized:
        pd_proc.numpy_intensity_plus = numpy.round(pd_proc.numpy_intensity_plus_net + 0.5 * pd_proc.numpy_intensity_bkg_calc, 3)
        pd_proc.numpy_intensity_minus = numpy.round(pd_proc.numpy_intensity_plus_net + 0.5 * pd_proc.numpy_intensity_bkg_calc, 3)
        pd_meas.numpy_intensity_plus = numpy.round(pd_proc.numpy_intensity_plus, 3)
        pd_meas.numpy_intensity_minus = numpy.round(pd_proc.numpy_intensity_minus, 3)
    else:
        pd_proc.numpy_intensity = numpy.round(pd_proc.numpy_intensity_plus_net + pd_proc.numpy_intensity_minus_net + pd_proc.numpy_intensity_bkg_calc, 3)
        pd_meas.numpy_intensity = numpy.round(pd_proc.numpy_intensity, 3)

    pd_proc.numpy_to_items()
    pd_meas.numpy_to_items()
    return pd