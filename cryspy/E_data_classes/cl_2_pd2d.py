"""Description of Pd2d class."""
__author__ = 'ikibalin'
__version__ = "2020_08_19"
from warnings import warn
import numpy
from typing import List, NoReturn

from cryspy.A_functions_base.function_2_crystallography_base import \
    calc_cos_ang

from cryspy.A_functions_base.matrix_operations import calc_m1_m2_m1t
from cryspy.A_functions_base.unit_cell import calc_matrix_t
from cryspy.A_functions_base.powder_diffraction_const_wavelength import \
    calc_ttheta_phi_by_gamma_nu

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN
from cryspy.B_parent_classes.cl_3_data import DataN
from cryspy.B_parent_classes.preocedures import take_items_by_class

from cryspy.C_item_loop_classes.cl_1_setup import Setup
from cryspy.C_item_loop_classes.cl_1_diffrn_radiation import \
    DiffrnRadiation
from cryspy.C_item_loop_classes.cl_1_extinction import Extinction
from cryspy.C_item_loop_classes.cl_1_phase import PhaseL
from cryspy.C_item_loop_classes.cl_1_pd_peak import PdPeakL
from cryspy.C_item_loop_classes.cl_1_refln import ReflnL
from cryspy.C_item_loop_classes.cl_1_refln_susceptibility import \
    ReflnSusceptibilityL
from cryspy.C_item_loop_classes.cl_1_refine_ls import RefineLs
from cryspy.C_item_loop_classes.cl_1_pd2d_background import \
    Pd2dBackground
from cryspy.C_item_loop_classes.cl_1_pd2d_instr_reflex_asymmetry import\
    Pd2dInstrReflexAsymmetry
from cryspy.C_item_loop_classes.cl_1_pd2d_instr_resolution import \
    Pd2dInstrResolution
from cryspy.C_item_loop_classes.cl_1_pd2d_meas import Pd2dMeas
from cryspy.C_item_loop_classes.cl_1_pd2d_proc import Pd2dProc
from cryspy.C_item_loop_classes.cl_1_pd2d_peak import Pd2dPeakL
from cryspy.C_item_loop_classes.cl_1_chi2 import Chi2
from cryspy.C_item_loop_classes.cl_1_range import Range
from cryspy.C_item_loop_classes.cl_1_texture import TextureL
from cryspy.C_item_loop_classes.cl_1_exclude import ExcludeL

from cryspy.E_data_classes.cl_1_crystal import Crystal

na = numpy.newaxis

class Pd2d(DataN):
    """
    Powder diffraction experiment with polarized or unpolarized neutrons (2d).

    Data items in the DIFFRN category record details about
    2d powder diffraction measurements.

    Methods
    -------
        - calc_iint_u_d_flip_ratio
        - calc_fr
        - calc_fm_perp_loc
        - calc_chi_sq
        - params_to_cif
        - data_to_cif
        - calc_to_cif
        - estimate_FM


    Attributes
    ----------
        - setup (mandatory)
        - pd2d_instr_resolution (mandatory)
        - phase (mandatory)
        - pd2d_background (mandatory)
        - pd_meas (mandatory)
        - diffrn_radiation
        - chi2
        - range
        - extinction
        - pd2d_instr_reflex_asymmetry
        - texture
        - exclude
        - pd2d_proc
        - pd2d_peak
        - refine_ls
        - refln_#phase_name
        - refln_susceptibility_#phase_name
    """

    CLASSES_MANDATORY = (Pd2dInstrResolution, PhaseL, DiffrnRadiation,
                         Setup, Range, Pd2dBackground, Pd2dMeas)
    CLASSES_OPTIONAL = (Extinction, ExcludeL, Chi2, Pd2dInstrReflexAsymmetry,
                        TextureL, RefineLs, ReflnL, ReflnSusceptibilityL,
                        Pd2dPeakL, Pd2dProc, PdPeakL)
    # CLASSES_INTERNAL = ()

    CLASSES = CLASSES_MANDATORY + CLASSES_OPTIONAL

    PREFIX = "pd2d"

    # default values for the parameters
    D_DEFAULT = {}

    def __init__(self, data_name=None, **kwargs) -> NoReturn:
        super(Pd2d, self).__init__()

        self.__dict__["items"] = []
        self.__dict__["data_name"] = data_name

        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)
    #
    # def calc_profile(self, tth, phi, l_crystal: List[Crystal], l_peak_in=None,
    #                  l_refln_in=None, l_refln_susceptibility_in=None,
    #                  l_dd_in=None, flag_internal: bool = True):
    #     """Calculate intensity for the given diffraction angle."""
    #     proc = Pd2dProc()  # it is output
    #     proc.ttheta = tth
    #     proc.phi = phi
    #
    #     background = self.pd2d_background
    #     int_bkgd = background.interpolate_by_points(tth, phi)
    #     proc.intensity_bkg_calc = int_bkgd
    #
    #     tth_rad = tth*numpy.pi/180.
    #     phi_rad = phi*numpy.pi/180.
    #     cos_theta_1d = numpy.cos(0.5*tth_rad)
    #     sin_phi_1d = numpy.sin(phi_rad)
    #
    #     setup = self.setup
    #     wavelength = setup.wavelength
    #     diffrn_radiation = self.diffrn_radiation
    #     if setup.offset_phi is None:
    #         setup.offset_phi = 0.
    #
    #     phi_0 = setup.offset_phi
    #     phi_rad = (phi-phi_0)*numpy.pi/180.
    #     sin_phi_1d = numpy.sin(phi_rad)
    #
    #     p_u = float(diffrn_radiation.polarization)
    #     p_d = (2.*float(diffrn_radiation.efficiency)-1.)*p_u
    #
    #     tth_min = tth.min()
    #     tth_max = tth.max()+3.
    #     if tth_max > 180.:
    #         tth_max = 180.
    #
    #     sthovl_min = numpy.sin(0.5*tth_min*numpy.pi/180.)/wavelength
    #     sthovl_max = numpy.sin(0.5*tth_max*numpy.pi/180.)/wavelength
    #
    #     res_u_2d = numpy.zeros((tth.shape[0], phi.shape[0]), dtype=float)
    #     res_d_2d = numpy.zeros((tth.shape[0], phi.shape[0]), dtype=float)
    #
    #     try:
    #         texture = self.texture
    #         h_ax, k_ax = float(texture.h_ax[0]), float(texture.k_ax[0]) # FIXME
    #         l_ax = float(texture.l_ax[0])
    #         g_1, g_2 = float(texture.g_1[0]), float(texture.g_2[0])
    #     except AttributeError:
    #         texture = None
    #
    #     phase = self.phase
    #     _h = len(phase.label)
    #     if l_peak_in is None:
    #         l_peak_in = len(phase.label) * [None]
    #     else:
    #         if _h != len(l_peak_in):
    #             l_peak_in = len(phase.label) * [None]
    #     if l_refln_in is None:
    #         l_refln_in = len(phase.label) * [None]
    #     else:
    #         if _h != len(l_refln_in):
    #             l_refln_in = len(phase.label) * [None]
    #     if l_refln_susceptibility_in is None:
    #         l_refln_susceptibility_in = len(phase.label) * [None]
    #     else:
    #         if _h != len(l_refln_susceptibility_in):
    #             l_refln_susceptibility_in = len(phase.label) * [None]
    #     if l_dd_in is None:
    #         l_dd_in = len(phase.label) * [None]
    #     else:
    #         if _h != len(l_dd_in):
    #             l_dd_in = len(phase.label) * [None]
    #
    #     l_peak, l_dd_out = [], []
    #     l_refln, l_refln_s = [], []
    #     for item_phase, peak_in, refln_in, refln_susceptibility_in, dd_in in \
    #             zip(phase.items, l_peak_in, l_refln_in,
    #                 l_refln_susceptibility_in, l_dd_in):
    #         phase_label = item_phase.label
    #         phase_scale = item_phase.scale
    #         try:
    #             phase_igsize = item_phase.igsize
    #             if phase_igsize is None: phase_igsize = 0. # temporary solution
    #         except AttributeError:
    #             phase_igsize = 0.
    #         try:
    #             phase_u = item_phase.u
    #             if phase_u is None: phase_u = 0. # temporary solution
    #         except AttributeError:
    #             phase_u = 0.
    #         try:
    #             phase_v = item_phase.v
    #             if phase_v is None: phase_v = 0. # temporary solution
    #         except AttributeError:
    #             phase_v = 0.
    #         try:
    #             phase_w = item_phase.w
    #             if phase_w is None: phase_w = 0. # temporary solution
    #         except AttributeError:
    #             phase_w = 0.
    #         try:
    #             phase_x = item_phase.x
    #             if phase_x is None: phase_x = 0. # temporary solution
    #         except AttributeError:
    #             phase_x = 0.
    #         try:
    #             phase_y = item_phase.y
    #             if phase_y is None: phase_y = 0. # temporary solution
    #         except AttributeError:
    #             phase_y = 0.
    #
    #         dd_out = {}
    #         for i_crystal, crystal in enumerate(l_crystal):
    #             if crystal.data_name.lower() == phase_label.lower():
    #                 ind_cry = i_crystal
    #                 break
    #         if ind_cry is None:
    #             warn(f"Crystal with name '{phase_label:}' is not found.",
    #                  UserWarning)
    #             return
    #
    #         crystal = l_crystal[ind_cry]
    #
    #         cell = crystal.cell
    #         # space_group = crystal.space_group
    #
    #         if peak_in is not None:
    #             index_h = peak_in.numpy_index_h
    #             index_k = peak_in.numpy_index_k
    #             index_l = peak_in.numpy_index_l
    #             mult = peak_in.numpy_index_mult
    #         else:
    #             if texture is None:
    #                 ind_hkl_mult = crystal.calc_hkl(sthovl_min, sthovl_max)
    #             else:
    #                 ind_hkl_mult = crystal.calc_hkl_in_range(sthovl_min, sthovl_max)
    #             index_h, index_k, index_l, mult = ind_hkl_mult[0], ind_hkl_mult[1], ind_hkl_mult[2], ind_hkl_mult[3]
    #
    #         peak = Pd2dPeakL(loop_name=phase_label)
    #         peak.numpy_index_h = index_h
    #         peak.numpy_index_k = index_k
    #         peak.numpy_index_l = index_l
    #         peak.numpy_index_mult = mult
    #         peak.numpy_to_items()
    #
    #         cond_1 = len(crystal.get_variable_names()) == 0
    #         cond_2 = (peak_in is not None) & (refln_in is not None)
    #
    #         if (cond_1 & cond_2):
    #             f_nucl_sq = peak_in.numpy_f_nucl_sq
    #             f_m_p_sin_sq = peak_in.numpy_f_m_p_sin_sq
    #             f_m_p_cos_sq = peak_in.numpy_f_m_p_cos_sq
    #             cross_sin = peak_in.numpy_cross_sin
    #             refln = refln_in
    #             refln_s = refln_susceptibility_in
    #         else:
    #             f_nucl_sq, f_m_p_sin_sq, f_m_p_cos_sq, cross_sin, refln, \
    #                 refln_s = self.calc_for_iint(index_h, index_k, index_l,
    #                                              crystal)
    #             refln.loop_name = phase_label
    #             refln_s.loop_name = phase_label
    #         l_refln.append(refln)
    #         l_refln_s.append(refln_s)
    #
    #         peak.numpy_f_nucl_sq = f_nucl_sq
    #         peak.numpy_f_m_p_sin_sq = f_m_p_sin_sq
    #         peak.numpy_f_m_p_cos_sq = f_m_p_cos_sq
    #         peak.numpy_cross_sin = cross_sin
    #
    #         cond_1 = dd_in is not None
    #         cond_2 = ((len(crystal.get_variable_names()) == 0) &
    #                   (not(diffrn_radiation.polarization_refinement)) &
    #                   (not(diffrn_radiation.efficiency_refinement)))
    #         if cond_1 & cond_2:
    #             iint_u_3d, iint_d_3d = dd_in["iint_u_3d"], dd_in["iint_d_3d"]
    #             cos_theta_3d = dd_in["cos_theta_3d"]
    #             sin_phi_3d = dd_in["sin_phi_3d"]
    #         else:
    #             cos_theta_3d, sin_phi_3d, mult_f_n_3d = numpy.meshgrid(
    #                 cos_theta_1d, sin_phi_1d, mult*f_nucl_sq, indexing="ij")
    #
    #             mult_f_m_c_3d = numpy.meshgrid(
    #                 tth_rad, phi_rad, mult*f_m_p_cos_sq, indexing="ij")[2]
    #
    #             hh_u_s_3d = numpy.meshgrid(
    #                 tth_rad, phi_rad, mult*(f_m_p_sin_sq+p_u*cross_sin),
    #                 indexing="ij")[2]
    #             hh_d_s_3d = numpy.meshgrid(
    #                 tth_rad, phi_rad, mult*(f_m_p_sin_sq-p_d*cross_sin),
    #                 indexing="ij")[2]
    #
    #             c_a_sq_3d = (cos_theta_3d * sin_phi_3d)**2
    #             s_a_sq_3d = 1.-c_a_sq_3d
    #             iint_u_3d = (mult_f_n_3d + hh_u_s_3d*s_a_sq_3d +
    #                          mult_f_m_c_3d*c_a_sq_3d)
    #             iint_d_3d = (mult_f_n_3d + hh_d_s_3d*s_a_sq_3d +
    #                          mult_f_m_c_3d*c_a_sq_3d)
    #         dd_out["cos_theta_3d"] = cos_theta_3d
    #         dd_out["sin_phi_3d"] = sin_phi_3d
    #         dd_out["iint_u_3d"] = iint_u_3d
    #         dd_out["iint_d_3d"] = iint_d_3d
    #
    #         sthovl_hkl = cell.calc_sthovl(index_h, index_k, index_l)
    #         tth_hkl_rad = numpy.where(sthovl_hkl*wavelength < 1.,
    #                                   2.*numpy.arcsin(sthovl_hkl*wavelength),
    #                                   numpy.pi)
    #         tth_hkl = tth_hkl_rad*180./numpy.pi
    #
    #         cond_1 = dd_in is not None
    #         flag_phase_uvwxy_ref = (
    #             item_phase.igsize_refinement | item_phase.u_refinement |
    #             item_phase.v_refinement | item_phase.w_refinement |
    #             item_phase.x_refinement | item_phase.y_refinement)
    #         cond_2 = ((not(flag_phase_uvwxy_ref)) &
    #                   (not(self.pd2d_instr_resolution.is_variables())) &
    #                   (not(setup.is_variables())))
    #         if cond_1 & cond_2:
    #             profile_3d, tth_zs = dd_in["profile_3d"], dd_in["tth_zs"]
    #             h_pv = dd_in["h_pv"]
    #         else:
    #             profile_3d, tth_zs, h_pv = self.calc_shape_profile(
    #                 tth, phi, tth_hkl, phase_igsize=phase_igsize,
    #                 phase_u=phase_u, phase_v=phase_v, phase_w=phase_w,
    #                 phase_x=phase_x, phase_y=phase_y)
    #         dd_out["profile_3d"] = profile_3d
    #         dd_out["tth_zs"] = tth_zs
    #         dd_out["h_pv"] = h_pv
    #
    #         peak.numpy_ttheta = tth_hkl + setup.offset_ttheta
    #         peak.width_ttheta = h_pv
    #         peak.numpy_to_items()
    #
    #         # texture
    #         if texture is not None:
    #             cond_1 = dd_in is not None
    #             cond_2 = ((not(setup.offset_phi_refinement)))
    #             if cond_1 & cond_2:
    #                 cos_alpha_ang_3d = dd_in["cos_alpha_ang_3d"]
    #                 sin_alpha_ang_3d = dd_in["sin_alpha_ang_3d"]
    #             else:
    #                 cos_alpha_ang_3d = cos_theta_3d * sin_phi_3d
    #                 sin_alpha_ang_3d = numpy.sqrt(1.-cos_alpha_ang_3d**2)
    #             dd_out["cos_alpha_ang_3d"] = cos_alpha_ang_3d
    #             dd_out["sin_alpha_ang_3d"] = sin_alpha_ang_3d
    #             cond_2 = (cond_2 & (not(texture.is_variables())) &
    #                       (not(cell.is_variables())))
    #             if cond_1 & cond_2:
    #                 texture_3d = dd_in["texture_3d"]
    #             else:
    #                 cos_alpha_ax = calc_cos_ang(cell, h_ax, k_ax, l_ax,
    #                                             index_h, index_k, index_l)
    #
    #                 c_help = 1.-cos_alpha_ax**2
    #                 c_help[c_help < 0.] = 0.
    #                 sin_alpha_ax = numpy.sqrt(c_help)
    #                 cos_alpha_ax_3d = numpy.meshgrid(tth, phi, cos_alpha_ax,
    #                                                  indexing="ij")[2]
    #                 sin_alpha_ax_3d = numpy.meshgrid(tth, phi, sin_alpha_ax,
    #                                                  indexing="ij")[2]
    #                 cos_alpha_3d = cos_alpha_ax_3d*cos_alpha_ang_3d + \
    #                     sin_alpha_ax_3d*sin_alpha_ang_3d
    #                 texture_3d = g_2 + (1.-g_2) * \
    #                     (1./g_1 + (g_1**2-1./g_1)*cos_alpha_3d**2)**(-1.5)
    #             dd_out["texture_3d"] = texture_3d
    #
    #             profile_3d = profile_3d*texture_3d
    #
    #         res_u_3d = profile_3d*iint_u_3d
    #         res_d_3d = profile_3d*iint_d_3d
    #
    #         # 0.5 to have the same meaning for scale factor as in FullProf
    #         res_u_2d += 0.5*phase_scale*res_u_3d.sum(axis=2)
    #         res_d_2d += 0.5*phase_scale*res_d_3d.sum(axis=2)
    #         l_peak.append(peak)
    #         l_dd_out.append(dd_out)
    #
    #     proc.ttheta_corrected = tth_zs
    #     proc.intensity_plus_net = res_u_2d
    #     proc.intensity_minus_net = res_d_2d
    #     proc.intensity_plus_total = res_u_2d+int_bkgd
    #     proc.intensity_minus_total = res_d_2d+int_bkgd
    #     if flag_internal:
    #         if background.is_variables():
    #             background.form_ttheta_phi_intensity()
    #         proc.form_ttheta_phi_intensity_plus_net()
    #         proc.form_ttheta_phi_intensity_minus_net()
    #         proc.form_ttheta_phi_intensity_plus_total()
    #         proc.form_ttheta_phi_intensity_minus_total()
    #     l_calc_objs = l_refln + l_refln_s + l_peak
    #     l_calc_objs.append(proc)
    #     self.add_items(l_calc_objs)
    #
    #     return proc, l_peak, l_refln, l_dd_out
    #
    # def calc_chi_sq(self, l_crystal, flag_internal: bool = True):
    #     """
    #     Calculate chi square.
    #
    #     Arguments
    #     ---------
    #         - l_crystal: a list of Crystal objects of cryspy library
    #         - flag_internal: a flag to calculate internal objects
    #           (default is True)
    #
    #     Output arguments
    #     ----------------
    #         - chi_sq_val: chi square of flip ratio
    #           (Sum_i ((y_e_i - y_m_i) / sigma_i)**2)
    #         - n: number of measured reflections
    #     """
    #     meas = self.pd2d_meas
    #
    #     tth = meas.ttheta
    #     phi = meas.phi
    #     int_u_exp = meas.intensity_plus
    #     sint_u_exp = meas.intensity_plus_sigma
    #     int_d_exp = meas.intensity_minus
    #     sint_d_exp = meas.intensity_minus_sigma
    #
    #     l_peak_in, l_refln_in = [], []
    #     l_refln_susceptibility_in = []
    #     l_dd_in = []
    #     flag_1 = not(flag_internal)
    #
    #     try:
    #         dd = self.dd
    #         flag_2 = True
    #     except AttributeError:
    #         flag_2 = False
    #     if (flag_1 & flag_2):
    #         for phase_item in self.phase.items:
    #             crystal = None
    #             for cryst in l_crystal:
    #                 if cryst.data_name.lower() == phase_item.label.lower():
    #                     crystal = cryst
    #                     break
    #             attr_peak = f"pd2d_peak_{crystal.data_name:}"
    #             attr_refln = f"refln_{crystal.data_name:}"
    #             attr_refln_s = f"refln_susceptibility_{crystal.data_name:}"
    #             l_peak_in.append(getattr(self, attr_peak))
    #             l_refln_in.append(getattr(self, attr_refln))
    #             l_refln_susceptibility_in.append(getattr(self, attr_refln_s))
    #
    #     cond_tth_in = numpy.ones(tth.size, dtype=bool)
    #     cond_phi_in = numpy.ones(phi.size, dtype=bool)
    #     try:
    #         range_ = self.range
    #         cond_tth_in = numpy.logical_and(cond_tth_in, tth >=
    #                                         range_.ttheta_min)
    #         cond_tth_in = numpy.logical_and(cond_tth_in, tth <=
    #                                         range_.ttheta_max)
    #
    #         cond_phi_in = numpy.logical_and(cond_phi_in, phi >= range_.phi_min)
    #         cond_phi_in = numpy.logical_and(cond_phi_in, phi <= range_.phi_max)
    #
    #     except AttributeError:
    #         pass
    #     # cond_1_in, cond_2_in = numpy.meshgrid(cond_tth_in, cond_phi_in,
    #     #                                       indexing="ij")
    #     # cond_in = numpy.logical_and(cond_1_in, cond_2_in)
    #     tth_in = tth[cond_tth_in]
    #     phi_in = phi[cond_phi_in]
    #     int_u_exp_in = int_u_exp[cond_tth_in, :][:, cond_phi_in]
    #     sint_u_exp_in = sint_u_exp[cond_tth_in, :][:, cond_phi_in]
    #     int_d_exp_in = int_d_exp[cond_tth_in, :][:, cond_phi_in]
    #     sint_d_exp_in = sint_d_exp[cond_tth_in, :][:, cond_phi_in]
    #
    #     proc, l_peak, l_refln, l_dd_out = self.calc_profile(
    #         tth_in, phi_in, l_crystal, l_peak_in=l_peak_in,
    #         l_refln_in=l_refln_in,
    #         l_refln_susceptibility_in=l_refln_susceptibility_in,
    #         l_dd_in=l_dd_in, flag_internal=flag_internal)
    #     proc.intensity_plus = int_u_exp_in
    #     proc.intensity_plus_sigma = sint_u_exp_in
    #     proc.intensity_minus = int_d_exp_in
    #     proc.intensity_minus_sigma = sint_d_exp_in
    #
    #     # self.proc = proc
    #     # self.peaks = l_peak
    #     # self.reflns = l_refln
    #     self.dd = l_dd_out
    #
    #     int_u_mod = proc.intensity_plus_total
    #     int_d_mod = proc.intensity_minus_total
    #
    #     sint_sum_exp_in = (sint_u_exp_in**2 + sint_d_exp_in**2)**0.5
    #
    #     chi_sq_u = ((int_u_mod-int_u_exp_in)/sint_u_exp_in)**2
    #     chi_sq_d = ((int_d_mod-int_d_exp_in)/sint_d_exp_in)**2
    #     chi_sq_sum = ((int_u_mod+int_d_mod-int_u_exp_in-int_d_exp_in) /
    #                   sint_sum_exp_in)**2
    #     chi_sq_dif = ((int_u_mod-int_d_mod-int_u_exp_in+int_d_exp_in) /
    #                   sint_sum_exp_in)**2
    #
    #     cond_u = numpy.logical_not(numpy.isnan(chi_sq_u))
    #     cond_d = numpy.logical_not(numpy.isnan(chi_sq_d))
    #     cond_sum = numpy.logical_not(numpy.isnan(chi_sq_sum))
    #     cond_dif = numpy.logical_not(numpy.isnan(chi_sq_dif))
    #
    #     # exclude region
    #     try:
    #         exclude = self.exclude
    #         l_excl_tth_min = exclude.numpy_ttheta_min
    #         l_excl_tth_max = exclude.numpy_ttheta_max
    #         l_excl_phi_min = exclude.numpy_phi_min
    #         l_excl_phi_max = exclude.numpy_phi_max
    #         for excl_tth_min, excl_tth_max, excl_phi_min, excl_phi_max in \
    #             zip(l_excl_tth_min, l_excl_tth_max, l_excl_phi_min,
    #                 l_excl_phi_max):
    #             cond_1 = numpy.logical_or(tth_in < 1.*excl_tth_min,
    #                                       tth_in > 1.*excl_tth_max)
    #             cond_2 = numpy.logical_or(phi_in < 1.*excl_phi_min,
    #                                       phi_in > 1.*excl_phi_max)
    #             cond_11, cond_22 = numpy.meshgrid(cond_1, cond_2,
    #                                               indexing="ij")
    #             cond_12 = numpy.logical_or(cond_11, cond_22)
    #             cond_u = numpy.logical_and(cond_u, cond_12)
    #             cond_d = numpy.logical_and(cond_d, cond_12)
    #             cond_sum = numpy.logical_and(cond_sum, cond_12)
    #     except AttributeError:
    #         pass
    #
    #     chi_sq_u_val = (chi_sq_u[cond_u]).sum()
    #     n_u = cond_u.sum()
    #
    #     chi_sq_d_val = (chi_sq_d[cond_d]).sum()
    #     n_d = cond_d.sum()
    #
    #     chi_sq_sum_val = (chi_sq_sum[cond_sum]).sum()
    #     n_sum = cond_sum.sum()
    #
    #     chi_sq_dif_val = (chi_sq_dif[cond_dif]).sum()
    #     n_dif = cond_dif.sum()
    #     chi2 = self.chi2
    #     flag_u = chi2.up
    #     flag_d = chi2.down
    #     flag_sum = chi2.sum
    #     flag_dif = chi2.diff
    #
    #     chi_sq_val = (int(flag_u)*chi_sq_u_val + int(flag_d)*chi_sq_d_val +
    #                   int(flag_sum)*chi_sq_sum_val +
    #                   int(flag_dif)*chi_sq_dif_val)
    #     n = (int(flag_u)*n_u + int(flag_d)*n_d + int(flag_sum)*n_sum +
    #          int(flag_dif)*n_dif)
    #     # print(f"chi_sq_val/n: {chi_sq_val/n:.2f}    \
    #     #     chi_sq_val: {chi_sq_val: .2f}")
    #     # d_exp_out = {"chi_sq_val": chi_sq_val, "n": n}
    #     # d_exp_out.update(d_exp_prof_out)
    #
    #     if flag_internal:
    #         refine_ls = RefineLs(number_reflns=n,
    #                              goodness_of_fit_all=chi_sq_val/float(n),
    #                              weighting_scheme="sigma")
    #         self.refine_ls = refine_ls
    #         proc.form_ttheta_phi_intensity_bkg_calc()
    #         proc.form_ttheta_phi_intensity_plus()
    #         proc.form_ttheta_phi_intensity_plus_sigma()
    #         proc.form_ttheta_phi_intensity_minus()
    #         proc.form_ttheta_phi_intensity_minus_sigma()
    #     return chi_sq_val, n
    #
    # def calc_for_iint(self, index_h, index_k, index_l, crystal,
    #                   flag_internal: bool = True):
    #     """Calculate the integral intensity for h, k, l reflections."""
    #     index_hkl = numpy.stack([index_h, index_k, index_l], axis=0)
    #     setup = self.setup
    #     field = float(setup.field)
    #
    #     refln = crystal.calc_refln(index_hkl)
    #     f_nucl = refln.numpy_f_calc
    #
    #     f_nucl_sq = abs(f_nucl*f_nucl.conjugate())
    #
    #     refln_s = crystal.calc_refln_susceptibility(index_hkl)
    #     sft_11 = refln_s.numpy_chi_11_calc
    #     sft_12 = refln_s.numpy_chi_12_calc
    #     sft_13 = refln_s.numpy_chi_13_calc
    #     sft_21 = refln_s.numpy_chi_21_calc
    #     sft_22 = refln_s.numpy_chi_22_calc
    #     sft_23 = refln_s.numpy_chi_23_calc
    #     sft_31 = refln_s.numpy_chi_31_calc
    #     sft_32 = refln_s.numpy_chi_32_calc
    #     sft_33 = refln_s.numpy_chi_33_calc
    #     _ij = numpy.stack([sft_11, sft_12, sft_13, sft_21, sft_22, sft_23, sft_31, sft_32, sft_33], axis=0)
    #
    #     cell = crystal.cell
    #     unit_cell_parameters = cell.get_unit_cell_parameters()
    #     t_ij = calc_matrix_t(index_hkl, unit_cell_parameters, flag_unit_cell_parameters=False)[0]
    #     #SIGMA = T CHI T^-1 = T chi T^T (because T is rotation matrix, therefore T^-1 = T^T)
    #     th_ij = calc_m1_m2_m1t(t_ij, _ij)[0]
    #     th_11, th_12, th_13, th_22, th_23 = th_ij[0], th_ij[1], th_ij[2], th_ij[4], th_ij[5]
    #
    #     # f_m_p_sin_sq = (field**2)*abs(0.5*(th_11*th_11.conjugate()+\
    #     #                th_22*th_22.conjugate())+th_12*th_12.conjugate())
    #     # f_m_p_cos_sq = (field**2)*abs(th_13*th_13.conjugate()+\
    #     #                 th_23*th_23.conjugate())
    #     # f_m_p_field = 0.5*field*(th_11+th_22)
    #
    #     f_m_p_sin_sq = abs(0.5 * (th_11 * th_11.conjugate()+th_22 *
    #                               th_22.conjugate())+th_12 * th_12.conjugate())
    #     f_m_p_cos_sq = abs(th_13 * th_13.conjugate() +
    #                        th_23 * th_23.conjugate())
    #     f_m_p_field = 0.5 * (th_11+th_22)
    #     cross_sin = 2. * (f_nucl.real * f_m_p_field.real +
    #                       f_nucl.imag * f_m_p_field.imag)
    #
    #     return f_nucl_sq, f_m_p_sin_sq, f_m_p_cos_sq, cross_sin, refln, refln_s
    #
    # def _gauss_pd(self, tth_2d):
    #     """One dimensional gauss powder diffraction."""
    #     ag, bg = self.ag, self.bg
    #     val_1 = bg*tth_2d**2
    #     val_2 = numpy.where(val_1 < 5., numpy.exp(-val_1), 0.)
    #     self.gauss_pd = ag*val_2
    #
    # def _lor_pd(self, tth_2d):
    #     """One dimensional lorentz powder diffraction."""
    #     al, bl = self.al, self.bl
    #     self.lor_pd = al*1./(1.+bl*tth_2d**2)
    #
    # def calc_shape_profile(
    #         self, tth, phi, tth_hkl, phase_igsize: float = 0.,
    #         phase_u: float = 0., phase_v: float = 0., phase_w: float = 0.,
    #         phase_x: float = 0., phase_y: float = 0.):
    #     """
    #     Calculate profile in the range ttheta.
    #
    #     For reflections placed on
    #     ttheta_hkl with i_g parameter by default equal to zero
    #
    #     tth, phi, tth_hkl in degrees
    #     """
    #     setup = self.setup
    #     zero_shift = float(setup.offset_ttheta)
    #     tth_zs = tth-zero_shift
    #
    #     resolution = self.pd2d_instr_resolution
    #
    #     h_pv, eta, h_g, h_l, a_g, b_g, a_l, b_l = resolution.calc_resolution(
    #         tth_hkl, phase_igsize=phase_igsize, phase_u=phase_u,
    #         phase_v=phase_v, phase_w=phase_w, phase_x=phase_x, phase_y=phase_y)
    #
    #     tth_3d, phi_3d, tth_hkl_3d = numpy.meshgrid(tth_zs, phi, tth_hkl,
    #                                                 indexing="ij")
    #
    #     self.ag = numpy.meshgrid(tth_zs, phi, a_g, indexing="ij")[2]
    #     self.bg = numpy.meshgrid(tth_zs, phi, b_g, indexing="ij")[2]
    #     self.al = numpy.meshgrid(tth_zs, phi, a_l, indexing="ij")[2]
    #     self.bl = numpy.meshgrid(tth_zs, phi, b_l, indexing="ij")[2]
    #     eta_3d = numpy.meshgrid(tth_zs, phi, eta, indexing="ij")[2]
    #     self.eta = eta_3d
    #
    #     self._gauss_pd(tth_3d-tth_hkl_3d)
    #     self._lor_pd(tth_3d-tth_hkl_3d)
    #     g_pd2d_3d = self.gauss_pd
    #     l_pd2d_3d = self.lor_pd
    #
    #     np_shape_3d = eta_3d * l_pd2d_3d + (1.-eta_3d) * g_pd2d_3d
    #
    #     asymmetry = self.pd2d_instr_reflex_asymmetry
    #     np_ass_2d = asymmetry.calc_asymmetry(tth_zs, tth_hkl, h_pv)
    #     np_ass_3d = np_ass_2d[:, numpy.newaxis, :] * numpy.ones(
    #         phi.size, dtype=float)[numpy.newaxis, :, numpy.newaxis]
    #
    #     # Lorentz factor
    #     tth_rad = tth_zs*numpy.pi/180.
    #     np_lor_1d = 1./(numpy.sin(tth_rad)*numpy.sin(0.5*tth_rad))
    #     np_lor_3d = numpy.meshgrid(np_lor_1d, phi, tth_hkl, indexing="ij")[0]
    #
    #     profile_3d = np_shape_3d*np_ass_3d*np_lor_3d
    #
    #     return profile_3d, tth_zs, h_pv
    #
    # def params_to_cif(self, separator="_", flag: bool = False,
    #                   flag_minimal: bool = True) -> str:
    #     """Save parameters to cif format."""
    #     ls_out = []
    #     l_cls = (Pd2dBackground, Pd2dInstrResolution, PhaseL, DiffrnRadiation,
    #              Setup, Range, Chi2, Extinction, ExcludeL,
    #              Pd2dInstrReflexAsymmetry, TextureL)
    #     l_obj = [item for item in self.items if type(item) in l_cls]
    #     l_itemn = [item for item in l_obj if isinstance(item, ItemN)]
    #     l_loopn = [item for item in l_obj if isinstance(item, LoopN)]
    #     ls_out.extend([_.to_cif(separator=separator)+"\n" for _ in l_itemn])
    #     ls_out.extend([_.to_cif(separator=separator)+"\n" for _ in l_loopn])
    #     return "\n".join(ls_out)
    #
    # def data_to_cif(self, separator="_", flag: bool = False,
    #                 flag_minimal: bool = True) -> str:
    #     """Save data to cif format."""
    #     ls_out = []
    #     l_cls = (Pd2dMeas, )
    #     l_obj = [item for item in self.items if type(item) in l_cls]
    #     l_itemn = [item for item in l_obj if isinstance(item, ItemN)]
    #     l_loopn = [item for item in l_obj if isinstance(item, LoopN)]
    #     ls_out.extend([_.to_cif(separator=separator)+"\n" for _ in l_itemn])
    #     ls_out.extend([_.to_cif(separator=separator)+"\n" for _ in l_loopn])
    #     return "\n".join(ls_out)
    #
    # def calc_to_cif(self, separator="_", flag=False, flag_minimal=True) -> str:
    #     """Save calculations to cif format."""
    #     ls_out = []
    #     l_cls = (RefineLs, ReflnL, ReflnSusceptibilityL, Pd2dPeakL, Pd2dProc)
    #     l_obj = [item for item in self.items if type(item) in l_cls]
    #     l_itemn = [item for item in l_obj if isinstance(item, ItemN)]
    #     l_loopn = [item for item in l_obj if isinstance(item, LoopN)]
    #     ls_out.extend([_.to_cif(separator=separator)+"\n" for _ in l_itemn])
    #     ls_out.extend([_.to_cif(separator=separator)+"\n" for _ in l_loopn])
    #     return "\n".join(ls_out)

    def apply_constraints(self):
        """Apply constraints."""
        pass

    def plots(self):
        if self.is_attribute("pd2d_proc"):
            pd2d_proc = self.pd2d_proc
            fig = pd2d_proc.plot_gamma_nu()[0]
            l_pd_peak = [item for item in self.items if isinstance(item, PdPeakL)]
            [ax1, ax2] = fig.axes[:2]
            gamma_min, gamma_max = ax1.get_xlim()
            nu_min, nu_max = ax1.get_ylim()
            for pd_peak in l_pd_peak:
                if pd_peak.is_attribute("gamma"):
                    index_h = numpy.array(pd_peak.index_h, dtype=int)
                    index_k = numpy.array(pd_peak.index_k, dtype=int)
                    index_l = numpy.array(pd_peak.index_l, dtype=int)
                    gamma = numpy.array(pd_peak.gamma, dtype=float)
                    nu = numpy.array(pd_peak.nu, dtype=float)

                    flag_gn = numpy.logical_and(
                        numpy.logical_and(gamma>=gamma_min, gamma<=gamma_max),
                        numpy.logical_and(nu>=nu_min, nu<=nu_max))

                    if (pd_peak.is_attribute("intensity_plus") and pd_peak.is_attribute("intensity_minus")):
                        i_plus = numpy.array(pd_peak.intensity_plus, dtype=float)
                        i_minus = numpy.array(pd_peak.intensity_minus, dtype=float)
                        inv_lf = numpy.sin(gamma*numpy.pi/180)*numpy.sin(0.5*gamma*numpy.pi/180)
                        i_sum = numpy.abs(i_plus*flag_gn+i_minus*flag_gn)
                        i_difference = numpy.abs(i_plus*flag_gn-i_minus*flag_gn)
                        flag_i_plus = (i_sum-i_sum.min()) > 0.01 * (i_sum.max() - i_sum.min())
                        flag_i_difference = i_difference  > 0.01 * i_difference.max()
                        flag_i = numpy.logical_or(flag_i_plus, flag_i_difference)
                        flag_gn = numpy.logical_and(flag_gn, flag_i)

                    ax1.plot(gamma[flag_gn], 0.5*(nu[flag_gn]+nu_max), "k.", alpha=0.3)
                    ax2.plot(gamma[flag_gn], 0.5*(nu[flag_gn]+nu_max), "k.", alpha=0.3)
                    ax1.plot(gamma[flag_gn], nu_min+0.5*(-nu[flag_gn]+nu_max), "k.", alpha=0.3)
                    ax2.plot(gamma[flag_gn], nu_min+0.5*(-nu[flag_gn]+nu_max), "k.", alpha=0.3)

                    if (pd_peak.is_attribute("intensity_plus") and pd_peak.is_attribute("intensity_minus")):
                        i_plus = numpy.array(pd_peak.intensity_plus, dtype=float)
                        i_minus = numpy.array(pd_peak.intensity_minus, dtype=float)
                        inv_lf = numpy.sin(gamma*numpy.pi/180)*numpy.sin(0.5*gamma*numpy.pi/180)
                        i_sum = numpy.abs(i_plus*flag_gn+i_minus*flag_gn)
                        i_difference = numpy.abs(i_plus*flag_gn-i_minus*flag_gn)
                        flag_i_plus = (i_sum-i_sum.min()) > 0.2 * (i_sum.max() - i_sum.min())
                        flag_i_difference = i_difference  > 0.2 * i_difference.max()
                        flag_i = numpy.logical_or(flag_i_plus, flag_i_difference)
                        flag_gn = numpy.logical_and(flag_gn, flag_i)

                    l_gn = []
                    for i_h, i_k, i_l, i_g, i_n in zip(
                            index_h[flag_gn], index_k[flag_gn], index_l[flag_gn], gamma[flag_gn], nu[flag_gn]):
                        if not((i_g, i_n) in l_gn):
                            s_text = f"({i_h:}{i_k:}{i_l:})"
                            ax1.text(i_g, 0.5*(+i_n+nu_max), s_text, alpha=0.5, color="k")
                            ax2.text(i_g, 0.5*(+i_n+nu_max), s_text, alpha=0.5, color="k")
                            ax1.text(i_g, 0.5*(-i_n+nu_max)+nu_min, s_text, alpha=0.5, color="k")
                            ax2.text(i_g, 0.5*(-i_n+nu_max)+nu_min, s_text, alpha=0.5, color="k")
                            l_gn.append((i_g, i_n))
            return [(fig, ax1), ]
            if self.is_attribute("chi2"):
                flag_up = self.chi2.up
                flag_down = self.chi2.down
                flag_sum = self.chi2.sum
                flag_diff = self.chi2.diff
            else:
                flag_up, flag_down, flag_sum = False, False, True
                flag_diff = False

            if flag_sum:
                fig_s, ax_s = pd2d_proc.plot_projection_sum()
                ax_s.set_title(self.data_name + " - "+ax_s.title.get_text())
                y_min_s, y_max_s = ax_s.get_ylim()
                y_dist_s = y_max_s-y_min_s
                y_step_s = 0.

            if flag_diff:
                fig_d, ax_d = pd2d_proc.plot_projection_diff()
                ax_d.set_title(self.data_name + " - "+ax_d.title.get_text())
                y_min_d, y_max_d = ax_d.get_ylim()
                y_dist_d = y_max_d-y_min_d
                y_step_d = 0.

            for item in self.items:
                if isinstance(item, Pd2dPeakL):
                    np_tth = item.numpy_ttheta
                    if flag_sum:
                        ax_s.plot(np_tth, 0.*np_tth+y_min_s-y_step_s, "|",
                                  label=item.loop_name)
                        y_step_s += 0.05*y_dist_s
                    if flag_diff:
                        ax_d.plot(np_tth, 0.*np_tth+y_min_d-y_step_d, "|",
                                  label=item.loop_name)
                        y_step_d += 0.05*y_dist_d
            res = []
            if flag_sum:
                ax_s.legend(loc='upper right')
                res.append((fig_s, ax_s))
            if flag_diff:
                ax_d.legend(loc='upper right')
                res.append((fig_d, ax_d))
            return res
        elif self.is_attribute("pd2d_meas"):
            return self.pd2d_meas.plots()
        return []

    def get_dictionary(self):
        """Form dictionary. See documentation moduel CrysPy using Jupyter notebook.
        """
        ddict = super(Pd2d, self).get_dictionary()
        pd2d_meas, range_, exclude = None, None, None
        # self.form_object()
        # ddict = {}
        # setup, pd2d_meas, resolution = None, None, None
        # pd2d_background, range_, exclude = None, None, None
        # asymmetry, diffrn_radiation = None, None
        # phase, texture, chi2 = None, None, None
        #
        # l_obj = take_items_by_class(self, (Setup, ))
        # if len(l_obj) > 0:
        #     setup = l_obj[0]

        l_obj = take_items_by_class(self, (Pd2dMeas, ))
        if len(l_obj) > 0:
            pd2d_meas = l_obj[0]

        # l_obj = take_items_by_class(self, (Pd2dInstrResolution, ))
        # if len(l_obj) > 0:
        #     resolution = l_obj[0]
        #
        # l_obj = take_items_by_class(self, (Pd2dBackground, ))
        # if len(l_obj) > 0:
        #     pd2d_background = l_obj[0]

        l_obj = take_items_by_class(self, (Range, ))
        if len(l_obj) > 0:
            range_ = l_obj[0]

        l_obj = take_items_by_class(self, (ExcludeL, ))
        if len(l_obj) > 0:
            exclude = l_obj[0]

        # l_obj = take_items_by_class(self, (Pd2dInstrReflexAsymmetry, ))
        # if len(l_obj) > 0:
        #     asymmetry = l_obj[0]
        #
        # l_obj = take_items_by_class(self, (DiffrnRadiation, ))
        # if len(l_obj) > 0:
        #     diffrn_radiation = l_obj[0]
        #
        # l_obj = take_items_by_class(self, (PhaseL, ))
        # if len(l_obj) > 0:
        #     phase = l_obj[0]
        #
        # l_obj = take_items_by_class(self, (TextureL, ))
        # if len(l_obj) > 0:
        #     texture = l_obj[0]
        #
        # l_obj = take_items_by_class(self, (Chi2, ))
        # if len(l_obj) > 0:
        #     chi2 = l_obj[0]
        #
        # ddict["name"] = self.data_name
        # ddict["type_name"] = self.get_name()
        # if setup is not None:
        #     ddict["magnetic_field"] = numpy.array([setup.field], dtype=float)
        #     ddict["wavelength"] = numpy.array([setup.wavelength], dtype=float)
        #     ddict["flags_wavelength"] = numpy.array([setup.wavelength_refinement], dtype=bool)
        #     ddict["offset_gamma"] = numpy.array([setup.offset_gamma * numpy.pi/180.], dtype=float) # FIXME
        #     ddict["flags_offset_gamma"] = numpy.array([setup.offset_gamma_refinement], dtype=bool) # FIXME
        #     ddict["offset_nu"] = numpy.array([setup.offset_nu * numpy.pi/180.], dtype=float) # FIXME
        #     ddict["flags_offset_nu"] = numpy.array([setup.offset_nu_refinement], dtype=bool) # FIXME

        # if chi2 is not None:
        #     ddict["flag_chi_sq_sum"] = chi2.sum
        #     ddict["flag_chi_sq_difference"] = chi2.diff

        if pd2d_meas is not None:
            gamma_deg = numpy.array(pd2d_meas.gamma, dtype=float) #
            gamma_min_deg = range_.gamma_min
            gamma_max_deg = range_.gamma_max
            nu_deg = numpy.array(pd2d_meas.nu, dtype=float) #
            nu_min_deg = range_.nu_min
            nu_max_deg = range_.nu_max

            flag_gamma_in = numpy.logical_and(
                gamma_deg >= gamma_min_deg,
                gamma_deg <= gamma_max_deg)

            flag_nu_in = numpy.logical_and(
                nu_deg >= nu_min_deg,
                nu_deg <= nu_max_deg)

            flag_in = flag_gamma_in[:, na] * flag_nu_in[na, :]

            # ddict["ttheta"] = ttheta_deg[flag_in] * numpy.pi/180.
            ddict["gamma"] = gamma_deg[flag_gamma_in] * numpy.pi/180.
            ddict["nu"] = nu_deg[flag_nu_in] * numpy.pi/180.

            if pd2d_meas.is_attribute("gn_intensity_plus"):
                int_plus = numpy.array(pd2d_meas.gn_intensity_plus, dtype=float)[flag_gamma_in][:, flag_nu_in]
                s_int_plus = numpy.array(pd2d_meas.gn_intensity_plus_sigma, dtype=float)[flag_gamma_in][:, flag_nu_in]
                int_minus = numpy.array(pd2d_meas.gn_intensity_minus, dtype=float)[flag_gamma_in][:, flag_nu_in]
                s_int_minus = numpy.array(pd2d_meas.gn_intensity_minus_sigma, dtype=float)[flag_gamma_in][:, flag_nu_in]
                ddict["signal_exp_plus"] = numpy.stack([int_plus, s_int_plus], axis=0)
                ddict["signal_exp_minus"] = numpy.stack([int_minus, s_int_minus], axis=0)
            else:
                int_sum = numpy.array(pd2d_meas.gn_intensity, dtype=float)[flag_gamma_in][:, flag_nu_in]
                s_int_sum = numpy.array(pd2d_meas.gn_intensity_sigma, dtype=float)[flag_gamma_in][:, flag_nu_in]
                ddict["signal_exp"] = numpy.stack([int_sum, s_int_sum], axis=0)

            gamma_in_range = gamma_deg[flag_gamma_in]
            nu_in_range = nu_deg[flag_nu_in]
            ttheta_in_range, phi_in_range = calc_ttheta_phi_by_gamma_nu(
                gamma_in_range[:, na]*numpy.pi/180. - ddict["offset_gamma"],
                nu_in_range[na, :]*numpy.pi/180. - ddict["offset_nu"],
                flag_gamma=False, flag_nu=False)[:2]
            ttheta_in_range = ttheta_in_range*180./numpy.pi
            phi_in_range = phi_in_range*180./numpy.pi
            flag_exclude = numpy.zeros(gamma_in_range.shape + nu_in_range.shape, dtype=bool)
            if exclude is not None:
                for item_e in exclude.items:
                    if item_e.is_attribute("gamma_min"):
                        if item_e.gamma_min is not None:
                            flag_gamma_in_1 = numpy.logical_and(
                                gamma_in_range >= item_e.gamma_min,
                                gamma_in_range <= item_e.gamma_max)
                            flag_nu_in_1 = numpy.logical_and(
                                nu_in_range >= item_e.nu_min,
                                nu_in_range <= item_e.nu_max)
                            flag_in_1 = flag_gamma_in_1[:, na] * flag_nu_in_1[na, :]
                            flag_exclude = numpy.logical_or(flag_exclude, flag_in_1)
                    if item_e.is_attribute("ttheta_min"):
                        if item_e.ttheta_min is not None:
                            flag_ttheta_in_1 = numpy.logical_and(
                                ttheta_in_range >= item_e.ttheta_min,
                                ttheta_in_range <= item_e.ttheta_max)
                            flag_phi_in_1 = numpy.logical_and(
                                phi_in_range >= item_e.phi_min,
                                phi_in_range <= item_e.phi_max)
                            flag_in_1 = numpy.logical_and(
                                flag_ttheta_in_1,
                                flag_phi_in_1)
                            flag_exclude = numpy.logical_or(flag_exclude, flag_in_1)

            ddict["excluded_points"] = flag_exclude

        # if resolution is not None:
        #     ddict["resolution_parameters"] = numpy.array([
        #         resolution.u, resolution.v, resolution.w,
        #         resolution.x, resolution.y], dtype=float)
        #
        #     ddict["flags_resolution_parameters"] = numpy.array([
        #         resolution.u_refinement, resolution.v_refinement, resolution.w_refinement,
        #         resolution.x_refinement, resolution.y_refinement], dtype=bool)
        #
        #     ddict["resolution_phi_parameter"] = numpy.array([
        #         resolution.phi,], dtype=float)
        #     ddict["flags_resolution_phi_parameter"] = numpy.array([
        #         resolution.phi_refinement,], dtype=float)
        #
        # if asymmetry is not None:
        #     ddict["asymmetry_parameters"] = numpy.array([
        #         asymmetry.p1, asymmetry.p2, asymmetry.p3,
        #         asymmetry.p4], dtype=float)
        #
        #     ddict["flags_asymmetry_parameters"] = numpy.array([
        #         asymmetry.p1_refinement, asymmetry.p2_refinement, asymmetry.p3_refinement,
        #         asymmetry.p4_refinement], dtype=bool)
        # else:
        #     ddict["asymmetry_parameters"] = numpy.array([
        #         0, 0, 0, 0], dtype=float)
        #
        #     ddict["flags_asymmetry_parameters"] = numpy.array([
        #         False, False, False, False], dtype=bool)
        #
        # if diffrn_radiation is not None:
        #     beam_polarization = diffrn_radiation.polarization
        #     flipper_efficiency = diffrn_radiation.efficiency
        #     ddict["beam_polarization"] = numpy.array([beam_polarization], dtype=float)
        #     ddict["flipper_efficiency"] = numpy.array([flipper_efficiency], dtype=float)
        #     ddict["flags_beam_polarization"] = numpy.array([diffrn_radiation.polarization_refinement], dtype=bool)
        #     ddict["flags_flipper_efficiency"] = numpy.array([diffrn_radiation.efficiency_refinement], dtype=bool)
        #
        # if phase is not None:
        #     ddict["phase_name"] = numpy.array(phase.label, dtype=str)
        #     if phase.is_attribute("u"):
        #         p_u = numpy.array(phase.u, dtype=float)
        #         r_u = numpy.array(phase.u_refinement, dtype=bool)
        #     else:
        #         p_u = numpy.zeros((len(phase.items),), dtype=float)
        #         r_u = numpy.zeros((len(phase.items),), dtype=bool)
        #
        #     if phase.is_attribute("v"):
        #         p_v = numpy.array(phase.v, dtype=float)
        #         r_v = numpy.array(phase.v_refinement, dtype=bool)
        #     else:
        #         p_v = numpy.zeros((len(phase.items),), dtype=float)
        #         r_v = numpy.zeros((len(phase.items),), dtype=bool)
        #
        #     if phase.is_attribute("w"):
        #         p_w = numpy.array(phase.w, dtype=float)
        #         r_w = numpy.array(phase.w_refinement, dtype=bool)
        #     else:
        #         p_w = numpy.zeros((len(phase.items),), dtype=float)
        #         r_w = numpy.zeros((len(phase.items),), dtype=bool)
        #
        #     if phase.is_attribute("x"):
        #         p_x = numpy.array(phase.x, dtype=float)
        #         r_x = numpy.array(phase.x_refinement, dtype=bool)
        #     else:
        #         p_x = numpy.zeros((len(phase.items),), dtype=float)
        #         r_x = numpy.zeros((len(phase.items),), dtype=bool)
        #
        #     if phase.is_attribute("y"):
        #         p_y = numpy.array(phase.y, dtype=float)
        #         r_y = numpy.array(phase.y_refinement, dtype=bool)
        #     else:
        #         p_y = numpy.zeros((len(phase.items),), dtype=float)
        #         r_y = numpy.zeros((len(phase.items),), dtype=bool)
        #
        #     ddict["phase_resolution_parameters"] = numpy.stack([p_u, p_v, p_w, p_x, p_y], axis=0)
        #     ddict["flags_phase_resolution_parameters"] = numpy.stack([r_u, r_v, r_w, r_x, r_y], axis=0)
        #
        #     if phase.is_attribute("igsize"):
        #         ddict["phase_ig"] = numpy.array(phase.igsize, dtype=float)
        #         ddict["flags_phase_ig"] = numpy.array(phase.igsize_refinement, dtype=bool)
        #     else:
        #         ddict["phase_ig"] = numpy.zeros((len(phase.items),), dtype=float)
        #         ddict["flags_phase_ig"] = numpy.zeros((len(phase.items),), dtype=bool)
        #
        #     if phase.is_attribute("scale"):
        #         ddict["phase_scale"] = numpy.array(phase.scale, dtype=float)
        #         ddict["flags_phase_scale"] = numpy.array(phase.scale_refinement, dtype=bool)
        #     else:
        #         ddict["phase_scale"] = numpy.zeros((len(phase.items),), dtype=float)
        #         ddict["flags_phase_scale"] = numpy.zeros((len(phase.items),), dtype=bool)
        #
        # if texture is not None:
        #     ddict["texture_name"] = numpy.array(texture.label, dtype=str)
        #     ddict["texture_g1"] = numpy.array(texture.g_1, dtype=float)
        #     ddict["texture_g2"] = numpy.array(texture.g_2, dtype=float)
        #     ddict["texture_axis"] = numpy.array(
        #         [texture.h_ax, texture.k_ax, texture.l_ax], dtype=float)
        #     ddict["flags_texture_g1"] = numpy.array(texture.g_1_refinement, dtype=bool)
        #     ddict["flags_texture_g2"] = numpy.array(texture.g_2_refinement, dtype=bool)
        #     ddict["flags_texture_axis"] = numpy.array(
        #         [texture.h_ax_refinement,
        #          texture.k_ax_refinement,
        #          texture.l_ax_refinement], dtype=bool)
        #
        # if pd2d_background is not None:
        #     ddict["background_gamma"] = numpy.array(pd2d_background.gamma, dtype=float)*numpy.pi/180.
        #     ddict["background_nu"] = numpy.array(pd2d_background.nu, dtype=float)*numpy.pi/180.
        #     ddict["background_intensity"] = numpy.array(pd2d_background.intensity, dtype=float)
        #     ddict["flags_background_intensity"] = numpy.array(pd2d_background.intensity_refinement, dtype=bool)

        return ddict

    def take_parameters_from_dictionary(self, ddict_diffrn, l_parameter_name: list=None, l_sigma: list=None):
        keys = ddict_diffrn.keys()
        if l_parameter_name is not None:
            parameter_label = [hh[0] for hh in l_parameter_name]
        else:
            parameter_label = []

        if "beam_polarization" in keys:
            self.diffrn_radiation.polarization = ddict_diffrn["beam_polarization"]

        if "flipper_efficiency" in keys:
            self.diffrn_radiation.efficiency = ddict_diffrn["flipper_efficiency"]

        if "resolution_parameters" in keys:
            np_resolution = ddict_diffrn["resolution_parameters"]
            resolution = self.pd2d_instr_resolution
            resolution.u = float(np_resolution[0])
            resolution.v = float(np_resolution[1])
            resolution.w = float(np_resolution[2])
            resolution.x = float(np_resolution[3])
            resolution.y = float(np_resolution[4])

        if "resolution_phi_parameter" in keys:
            np_resolution = ddict_diffrn["resolution_phi_parameter"]
            resolution = self.pd2d_instr_resolution
            resolution.phi = float(np_resolution[0])

        if "background_intensity" in keys:
            background_intensity = ddict_diffrn["background_intensity"]
            self.pd2d_background.intensity = numpy.round(background_intensity, decimals=5)
            self.pd2d_background.form_gamma_nu_intensity()

        if "texture_g1" in keys:
            hh_g1 = ddict_diffrn["texture_g1"]
            hh_g2 = ddict_diffrn["texture_g2"]
            hh_axis = ddict_diffrn["texture_axis"]
            for i_item, item in enumerate(self.texture):
                item.g_1  = float(hh_g1[i_item])
                item.g_2  = float(hh_g2[i_item])
                item.h_ax = float(hh_axis[0, i_item])
                item.k_ax = float(hh_axis[1, i_item])
                item.l_ax = float(hh_axis[2, i_item])

        if "asymmetry_parameters" in keys:
            hh = ddict_diffrn["asymmetry_parameters"]
            if self.is_attribute("pd2d_instr_reflex_asymmetry"):
                asymmetry = self.pd2d_instr_reflex_asymmetry
            else:
                asymmetry = Pd2dInstrReflexAsymmetry()
                self.items.append(asymmetry)
            asymmetry.p1 = float(hh[0])
            asymmetry.p2 = float(hh[1])
            asymmetry.p3 = float(hh[2])
            asymmetry.p4 = float(hh[3])


        if "phase_scale" in keys:
            hh = ddict_diffrn["phase_scale"]
            for i_item, item in enumerate(self.phase.items):
                item.scale = float(hh[i_item])

        if "phase_ig" in keys:
            hh = ddict_diffrn["phase_ig"]
            for i_item, item in enumerate(self.phase.items):
                item.igsize = float(hh[i_item])

        if "offset_gamma" in keys:
            self.setup.offset_gamma = float(ddict_diffrn["offset_gamma"]) * 180./numpy.pi

        if "offset_nu" in keys:
            self.setup.offset_nu = float(ddict_diffrn["offset_nu"]) * 180./numpy.pi

        if "wavelength" in keys:
            self.setup.wavelength = float(ddict_diffrn["wavelength"])

        if len(parameter_label) > 0:
            for name, sigma in zip(l_parameter_name, l_sigma):
                if name[0] == "background_intensity":
                    self.pd2d_background.intensity_sigma[name[1]] = float(sigma)
                if name[0] == "phase_scale":
                    self.phase.items[name[1][0]].scale_sigma = float(sigma)
                if name[0] == "phase_ig":
                    self.phase.items[name[1][0]].igsize_sigma = float(sigma)
                if name[0] == "wavelength":
                    self.setup.wavelength_sigma = float(sigma)
                if name[0] == "offset_gamma":
                    self.setup.offset_gamma_sigma = float(sigma) * 180./numpy.pi
                if name[0] == "offset_nu":
                    self.setup.offset_nu_sigma = float(sigma) * 180./numpy.pi
                if name[0] == "beam_polarization":
                    self.diffrn_radiation.polarization_sigma = float(sigma)
                if name[0] == "flipper_efficiency":
                    self.diffrn_radiation.efficiency_sigma = float(sigma)
                if name[0] == "asymmetry_parameters":
                    asymmetry = self.pd2d_instr_reflex_asymmetry
                    if name[1][0] == 0:
                        asymmetry.p1_sigma = float(sigma)
                    elif name[1][0] == 1:
                        asymmetry.p2_sigma = float(sigma)
                    elif name[1][0] == 2:
                        asymmetry.p3_sigma = float(sigma)
                    elif name[1][0] == 3:
                        asymmetry.p4_sigma = float(sigma)
                if name[0] == "resolution_parameters":
                    resolution = self.pd2d_instr_resolution
                    if name[1][0] == 0:
                        resolution.u_sigma = float(sigma)
                    elif name[1][0] == 1:
                        resolution.v_sigma = float(sigma)
                    elif name[1][0] == 2:
                        resolution.w_sigma = float(sigma)
                    elif name[1][0] == 3:
                        resolution.x_sigma = float(sigma)
                    elif name[1][0] == 4:
                        resolution.y_sigma = float(sigma)
                if name[0] == "resolution_phi_parameter":
                    resolution = self.pd2d_instr_resolution
                    resolution.phi_sigma = float(sigma)
                if name[0] == "texture_g1":
                    self.texture.items[name[1][0]].g1_sigma = float(sigma)
                if name[0] == "texture_g2":
                    self.texture.items[name[1][0]].g2_sigma = float(sigma)
                if name[0] == "texture_axis":
                    ind_param, ind_a = name[1]
                    if ind_param == 0:
                        self.texture.items[ind_a].h_ax_sigma = float(sigma)
                    if ind_param == 1:
                        self.texture.items[ind_a].k_ax_sigma = float(sigma)
                    if ind_param == 2:
                        self.texture.items[ind_a].l_ax_sigma = float(sigma)

            if "background_intensity" in parameter_label:
                self.pd2d_background.form_gamma_nu_intensity()

        if (("signal_plus" in keys) and ("signal_minus" in keys) and ("signal_background" in keys)):
            pd2d_proc = Pd2dProc()
            pd2d_proc.gamma = numpy.round(ddict_diffrn["gamma"] * 180./numpy.pi, decimals=5)
            pd2d_proc.nu = numpy.round(ddict_diffrn["nu"] * 180./numpy.pi, decimals=5)
            pd2d_proc.intensity_plus_net = numpy.round(ddict_diffrn["signal_plus"], decimals=5)
            pd2d_proc.intensity_minus_net = numpy.round(ddict_diffrn["signal_minus"], decimals=5)
            pd2d_proc.intensity_bkg_calc = numpy.round(ddict_diffrn["signal_background"], decimals=5)
            pd2d_proc.excluded_points = ddict_diffrn["excluded_points"][:, :]
            pd2d_proc.form_gamma_nu_intensity_plus_net()
            pd2d_proc.form_gamma_nu_intensity_minus_net()
            pd2d_proc.form_gamma_nu_intensity_bkg_calc()
            pd2d_proc.form_gamma_nu_excluded_points()

            if "signal_exp_plus" in keys:
                pd2d_proc.intensity_plus = numpy.round(ddict_diffrn["signal_exp_plus"][0, :, :], decimals=5)
                pd2d_proc.intensity_plus_sigma = numpy.round(ddict_diffrn["signal_exp_plus"][1, :, :], decimals=5)
                pd2d_proc.intensity_minus = numpy.round(ddict_diffrn["signal_exp_minus"][0, :, :], decimals=5)
                pd2d_proc.intensity_minus_sigma = numpy.round(ddict_diffrn["signal_exp_minus"][1, :, :], decimals=5)
                pd2d_proc.form_gamma_nu_intensity_plus()
                pd2d_proc.form_gamma_nu_intensity_plus_sigma()
                pd2d_proc.form_gamma_nu_intensity_minus()
                pd2d_proc.form_gamma_nu_intensity_minus_sigma()
            elif "signal_exp" in keys:
                pd2d_proc.intensity = numpy.round(ddict_diffrn["signal_exp"][0, :, :], decimals=5)
                pd2d_proc.intensity_sigma = numpy.round(ddict_diffrn["signal_exp"][1, :, :], decimals=5)
                pd2d_proc.form_gamma_nu_intensity()
                pd2d_proc.form_gamma_nu_intensity_sigma()

            self.pd2d_proc = pd2d_proc

        l_pd_peak = []
        for item_phase in self.phase.items:
            s_label = f"dict_in_out_{item_phase.label:}"
            if s_label in keys:
                dict_crystal = ddict_diffrn[s_label]
                dict_crystal_keys = dict_crystal.keys()
                if (("index_hkl" in dict_crystal_keys) and  ("ttheta_hkl" in dict_crystal_keys) and  ("f_nucl" in dict_crystal_keys)):
                    pd_peak = PdPeakL(loop_name = item_phase.label)
                    int_max = numpy.square(numpy.abs(dict_crystal["f_nucl"]))

                    if "f_m_perp_o" in dict_crystal_keys:
                        int_max += numpy.sum(numpy.square(numpy.abs(dict_crystal["f_m_perp_o"])), axis=0)
                    if "p_1" in dict_crystal_keys:
                        int_max += dict_crystal["p_1"]

                    if "p_3" in dict_crystal_keys:
                        int_plus_max = int_max+dict_crystal["p_3"]
                        int_minus_max = int_max-dict_crystal["p_3"]
                    else:
                        int_plus_max = int_max
                        int_minus_max = int_max

                    index_hkl = dict_crystal["index_hkl"]
                    ttheta_hkl = dict_crystal["ttheta_hkl"]
                    lf = numpy.sin(ttheta_hkl)*numpy.sin(0.5*ttheta_hkl)
                    lf[lf<0] = 1.
                    pd_peak.numpy_index_h = index_hkl[0]
                    pd_peak.numpy_index_k = index_hkl[1]
                    pd_peak.numpy_index_l = index_hkl[2]
                    pd_peak.numpy_ttheta = numpy.round(ttheta_hkl * 180./numpy.pi, decimals=3)
                    pd_peak.numpy_intensity_plus = int_plus_max/lf
                    pd_peak.numpy_intensity_minus = int_minus_max/lf
                    if "multiplicity_hkl" in dict_crystal_keys:
                        pd_peak.numpy_index_mult = dict_crystal["multiplicity_hkl"]
                    if "ttheta_hkl" in dict_crystal_keys:
                        pd_peak.numpy_ttheta = numpy.round(dict_crystal["ttheta_hkl"] * 180./numpy.pi, decimals=3)
                    if (("gamma_hkl" in dict_crystal_keys) and ("nu_hkl" in dict_crystal_keys)):
                        gamma_hkl = dict_crystal["gamma_hkl"]
                        nu_hkl = dict_crystal["nu_hkl"]
                        pd_peak.numpy_gamma = numpy.round(gamma_hkl * 180./numpy.pi, decimals=3)
                        pd_peak.numpy_nu = numpy.round(nu_hkl * 180./numpy.pi, decimals=3)
                    pd_peak.numpy_to_items()
                    l_pd_peak.append(pd_peak)
        if len(l_pd_peak) > 0:
            self.add_items(l_pd_peak)

    def estimate_background(self):
        pd2d_background = self.pd2d_background
        intensity_refinement = numpy.copy(pd2d_background.intensity_refinement)
        pd2d_background.intensity_refinement = numpy.ones_like(intensity_refinement, dtype=bool)
        setup=self.setup
        pd2d_proc = self.pd2d_proc
        res = pd2d_proc.estimate_background(
            pd2d_background,
            offset_gamma=setup.offset_gamma,
            offset_nu = setup.offset_nu)
        pd2d_background.intensity_refinement = intensity_refinement
        pd2d_background.form_gamma_nu_intensity()
        return res
