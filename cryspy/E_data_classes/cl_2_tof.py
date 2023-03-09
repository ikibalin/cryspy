"""Description of Pd class."""
__author__ = 'ikibalin'
__version__ = "2020_08_19"

import numpy
from typing import NoReturn

from cryspy.A_functions_base.function_2_crystallography_base import \
    calc_cos_ang

from cryspy.A_functions_base.matrix_operations import calc_m1_m2_m1t
from cryspy.A_functions_base.unit_cell import calc_matrix_t

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN
from cryspy.B_parent_classes.cl_3_data import DataN
from cryspy.B_parent_classes.preocedures import take_items_by_class

from cryspy.C_item_loop_classes.cl_1_diffrn_radiation import \
    DiffrnRadiation
from cryspy.C_item_loop_classes.cl_1_extinction import Extinction
from cryspy.C_item_loop_classes.cl_1_phase import PhaseL
from cryspy.C_item_loop_classes.cl_1_refln import ReflnL
from cryspy.C_item_loop_classes.cl_1_refln_susceptibility import \
    ReflnSusceptibilityL
from cryspy.C_item_loop_classes.cl_1_refine_ls import RefineLs
from cryspy.C_item_loop_classes.cl_1_chi2 import Chi2
from cryspy.C_item_loop_classes.cl_1_range import Range
from cryspy.C_item_loop_classes.cl_1_setup import Setup
from cryspy.C_item_loop_classes.cl_1_texture import Texture, TextureL
from cryspy.C_item_loop_classes.cl_1_exclude import ExcludeL

from cryspy.C_item_loop_classes.cl_1_tof_background import TOFBackground
from cryspy.C_item_loop_classes.cl_1_tof_intensity_incident import \
    TOFIntensityIncident
from cryspy.C_item_loop_classes.cl_1_tof_meas import TOFMeasL
from cryspy.C_item_loop_classes.cl_1_tof_parameters import TOFParameters
from cryspy.C_item_loop_classes.cl_1_tof_peak import TOFPeakL
from cryspy.C_item_loop_classes.cl_1_tof_proc import TOFProcL
from cryspy.C_item_loop_classes.cl_1_tof_profile import TOFProfile

from cryspy.E_data_classes.cl_1_crystal import Crystal
from cryspy.E_data_classes.cl_1_mag_crystal import MagCrystal


class TOF(DataN):
    """
    Time-of-flight powder diffraction (polarized or unpolarized neutrons, 1d).

    Data items in the TOF category record details about
    powder diffraction measurements by Time of FLight.

    Methods
    -------
        - calc_profile

    Attributes
    ----------
        - setup, tof_profile, phase, tof_background, tof_meas (mandatory)
        - diffrn_radiation, chi2, range, extinction,
          texture, exclude, tof_proc, tof_peak, refine_ls, refln_#phase_name
          refln_susceptibility_#phase_name (optional)
    """

    CLASSES_MANDATORY = (TOFParameters, TOFProfile, PhaseL, TOFBackground,
                         TOFMeasL)
    CLASSES_OPTIONAL = (DiffrnRadiation, Chi2, Range,
                        Texture, TOFIntensityIncident, ExcludeL, TOFProcL,
                        TOFPeakL, RefineLs, ReflnL, ReflnSusceptibilityL, Setup)
    # CLASSES_INTERNAL = ()

    CLASSES = CLASSES_MANDATORY + CLASSES_OPTIONAL

    PREFIX = "tof"

    # default values for the parameters
    D_DEFAULT = {}

    def __init__(self, data_name=None, **kwargs) -> NoReturn:
        super(TOF, self).__init__()

        self.__dict__["items"] = []
        self.__dict__["data_name"] = data_name

        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)

    # def calc_profile(self, time, l_crystal, flag_internal: bool = True,
    #                  flag_polarized: bool = True):
    #     """Calculate intensity for the given diffraction angles.
    #
    #     Arguments
    #     ---------
    #         - time: 1D numpy array of time in micro seconds (???)
    #         - l_crystal: a list of Crystal objects of cryspy library
    #         - l_peak_in: precalculated data about integrated intensities
    #           (used to speed up the calculations)
    #         - l_refln_in: precalculated data about nuclear structure factors
    #           (used to speed up the calculations)
    #         - l_refln_susceptibility_in: precalculated data about
    #           (used to speed up the calculations)
    #         - flag_internal: a flag to calculate or to use internal objects.
    #                It should be True if user call the function.
    #                It's True by default.
    #
    #     Output
    #     ------
    #         - proc: output profile
    #         - l_peak: data about peaks
    #         - l_refln: data about nuclear structure factor
    #     """
    #     if flag_internal:
    #         d_internal_val = {}
    #         self.d_internal_val = d_internal_val
    #     else:
    #         try:
    #             d_internal_val = self.d_internal_val
    #         except AttributeError:
    #             d_internal_val = {}
    #             self.d_internal_val = d_internal_val
    #
    #     tof_proc = TOFProcL()
    #     tof_proc.numpy_time = time
    #
    #     try:
    #         tof_background = self.tof_background
    #         int_bkgd = tof_background.calc_background(time)
    #     except AttributeError:
    #         int_bkgd = 0.*time
    #
    #     tof_proc.numpy_intensity_bkg_calc = int_bkgd
    #
    #     tof_parameters = self.tof_parameters
    #     # tof_profile = self.tof_profile
    #
    #     d = tof_parameters.calc_d_by_time(time)
    #     d_min, d_max = tof_parameters.calc_d_min_max(time)
    #     sthovl_min = 0.5/d_max
    #     sthovl_max = 0.5/d_min
    #
    #     try:
    #         texture = self.texture
    #         h_ax, k_ax, l_ax = texture.h_ax, texture.k_ax, texture.l_ax
    #         g_1, g_2 = texture.g_1, texture.g_2
    #     except AttributeError:
    #         texture = None
    #
    #     res_u_1d = numpy.zeros(time.shape[0], dtype=float)
    #     res_d_1d = numpy.zeros(time.shape[0], dtype=float)
    #
    #     phase = self.phase
    #
    #     for phase_item in phase.items:
    #         phase_label = phase_item.label
    #         phase_scale = phase_item.scale
    #         try:
    #             phase_igsize = phase_item.igsize
    #             if phase_igsize is None: phase_igsize = 0. # temporary solution
    #         except AttributeError:
    #             phase_igsize = 0.
    #
    #         for i_crystal, crystal in enumerate(l_crystal):
    #             if crystal.data_name.lower() == phase_label.lower():
    #                 ind_cry = i_crystal
    #                 break
    #         if ind_cry is None:
    #             raise AttributeError(
    #                 f"Crystal with name '{crystal.data_name:}' is not found.")
    #             return
    #         crystal = l_crystal[ind_cry]
    #
    #         scale = phase_scale
    #
    #         try:
    #             peak = d_internal_val[f"peak_{crystal.data_name:}"]
    #             index_h = peak.numpy_index_h
    #             index_k = peak.numpy_index_k
    #             index_l = peak.numpy_index_l
    #             mult = peak.numpy_index_mult
    #             cond_2 = True
    #         except KeyError:
    #             if texture is None:
    #                 index_h, index_k, index_l, mult = crystal.calc_hkl(
    #                     sthovl_min, sthovl_max)
    #             else:
    #                 index_h, index_k, index_l, mult = \
    #                     crystal.calc_hkl_in_range(sthovl_min, sthovl_max)
    #             peak = TOFPeakL(loop_name=phase_label)
    #             peak.numpy_index_h = numpy.array(index_h, dtype=int)
    #             peak.numpy_index_k = numpy.array(index_k, dtype=int)
    #             peak.numpy_index_l = numpy.array(index_l, dtype=int)
    #             peak.numpy_index_mult = numpy.array(mult, dtype=int)
    #             d_internal_val[f"peak_{crystal.data_name:}"] = peak
    #             cond_2 = False
    #
    #         cond_1 = not(crystal.is_variables())
    #         if cond_1 & cond_2:
    #             np_iint_u = peak.numpy_intensity_plus
    #             np_iint_d = peak.numpy_intensity_minus
    #         else:
    #             np_iint_u, np_iint_d = self.calc_iint(
    #                 index_h, index_k, index_l, crystal,
    #                 flag_internal=flag_internal)
    #             peak.numpy_intensity_plus = np_iint_u
    #             peak.numpy_intensity_minus = np_iint_d
    #
    #         cell = crystal.cell
    #         sthovl_hkl = cell.calc_sthovl(index_h, index_k, index_l)
    #         d_hkl = 0.5/sthovl_hkl
    #         time_hkl = tof_parameters.calc_time_by_d(d_hkl)
    #
    #
    #         profile_2d = self.calc_shape_profile(
    #             d, d_hkl, phase_igsize=phase_igsize)
    #
    #         peak.numpy_time = time_hkl
    #
    #         tth_rad = tof_parameters.ttheta_bank * numpy.pi/180.
    #         wavelength = 2.*d_hkl*numpy.sin(0.5*tth_rad)
    #         wavelength_4 = wavelength**4
    #
    #         iint_u_0 = np_iint_u * mult * wavelength_4
    #         iint_d_0 = np_iint_d * mult * wavelength_4
    #
    #         if tof_parameters.is_attribute("extinction"):
    #             tof_extinction = tof_parameters.extinction
    #             iint_u_0 = (1. - tof_extinction*iint_u_0)*mult*wavelength_4
    #             iint_d_0 = (1. - tof_extinction*iint_d_0)*mult*wavelength_4
    #
    #         np_iint_u_2d = numpy.meshgrid(time, iint_u_0, indexing="ij")[1]
    #         np_iint_d_2d = numpy.meshgrid(time, iint_d_0, indexing="ij")[1]
    #
    #         # texture
    #         if texture is not None:
    #             cond_2 = ((not(texture.is_variables())) &
    #                       (not(cell.is_variables())))
    #             if cond_2:
    #                 try:
    #                     texture_2d = d_internal_val[
    #                         f"texture_2d_{crystal.data_name:}"]
    #                     cond_1 = False
    #                 except KeyError:
    #                     cond_1 = True
    #             if cond_1:
    #                 cos_alpha_ax = calc_cos_ang(cell, h_ax, k_ax, l_ax,
    #                                             index_h, index_k, index_l)
    #                 c_help = 1.-cos_alpha_ax**2
    #                 c_help[c_help < 0.] = 0.
    #                 sin_alpha_ax = numpy.sqrt(c_help)
    #                 cos_alpha_ax_2d = numpy.meshgrid(time, cos_alpha_ax,
    #                                                  indexing="ij")[1]
    #                 sin_alpha_ax_2d = numpy.meshgrid(time, sin_alpha_ax,
    #                                                  indexing="ij")[1]
    #                 # cos_alpha_2d = cos_alpha_ax_2d*cos_alpha_ang_2d + \
    #                 #                sin_alpha_ax_2d*sin_alpha_ang_2d
    #                 # cos_alpha_ang_2d, sin_alpha_ang_2d
    #                 cos_alpha_2d = cos_alpha_ax_2d*1.+sin_alpha_ax_2d*0.
    #                 texture_2d = g_2 + (1. - g_2) * (1./g_1 + (g_1**2 - 1./g_1)
    #                                                  * cos_alpha_2d**2)**(-1.5)
    #                 d_internal_val[f"texture_2d_{crystal.data_name:}"] = \
    #                     texture_2d
    #
    #             profile_2d = profile_2d*texture_2d
    #
    #         res_u_2d = profile_2d*np_iint_u_2d
    #         res_d_2d = profile_2d*np_iint_d_2d
    #
    #         # 0.5 to have the same meaning for scale factor as in FullProf
    #         res_u_1d += 0.5*scale*res_u_2d.sum(axis=1)
    #         res_d_1d += 0.5*scale*res_d_2d.sum(axis=1)
    #
    #         if flag_internal:
    #             peak.numpy_to_items()
    #
    #     tof_proc.numpy_d_spacing = d
    #     tof_proc.numpy_intensity_plus_net = res_u_1d
    #     tof_proc.numpy_intensity_minus_net = res_d_1d
    #     tof_proc.numpy_intensity_net = res_u_1d+res_d_1d
    #     tof_proc.numpy_intensity_diff_total = res_u_1d-res_d_1d
    #     if flag_polarized:
    #         tof_proc.numpy_intensity_plus_total = res_u_1d+int_bkgd
    #         tof_proc.numpy_intensity_minus_total = res_d_1d+int_bkgd
    #         tof_proc.numpy_intensity_total = res_u_1d+res_d_1d+int_bkgd+int_bkgd
    #     else:
    #         tof_proc.numpy_intensity_plus_total = res_u_1d+0.5*int_bkgd
    #         tof_proc.numpy_intensity_minus_total = res_d_1d+0.5*int_bkgd
    #         tof_proc.numpy_intensity_total = res_u_1d+res_d_1d+int_bkgd
    #
    #     if flag_internal:
    #         tof_proc.numpy_to_items()
    #         l_calc_objs = []
    #         for crystal in l_crystal:
    #             try:
    #                 obj = d_internal_val[f"refln_{crystal.data_name}"]
    #                 l_calc_objs.append(obj)
    #             except KeyError:
    #                 pass
    #             try:
    #                 obj = d_internal_val[
    #                     f"refln_susceptibility_{crystal.data_name}"]
    #                 l_calc_objs.append(obj)
    #             except KeyError:
    #                 pass
    #             try:
    #                 obj = d_internal_val[f"peak_{crystal.data_name}"]
    #                 l_calc_objs.append(obj)
    #             except KeyError:
    #                 pass
    #         l_calc_objs.append(tof_proc)
    #         self.add_items(l_calc_objs)
    #     return tof_proc
    #
    # def calc_chi_sq(self, l_crystal, flag_internal=True):
    #     """
    #     Calculate chi square.
    #
    #     Arguments
    #     ---------
    #         - l_crystal: a list of Crystal objects of cryspy library
    #         - flag_internal: a flag to calculate or to use internal objects.
    #                It should be True if user call the function.
    #                It's True by default.
    #
    #     Output
    #     ------
    #         - chi_sq_val: chi square of flip ratio
    #           (Sum_i ((y_e_i - y_m_i) / sigma_i)**2)
    #         - n: number of measured reflections
    #     """
    #     tof_meas = self.tof_meas
    #     flag_polarized = tof_meas.is_polarized()
    #
    #     np_time = tof_meas.numpy_time
    #     if flag_polarized:
    #         int_u_exp = tof_meas.numpy_intensity_plus
    #         sint_u_exp = tof_meas.numpy_intensity_plus_sigma
    #         int_d_exp = tof_meas.numpy_intensity_minus
    #         sint_d_exp = tof_meas.numpy_intensity_minus_sigma
    #     else:
    #         int_exp = tof_meas.numpy_intensity
    #         sint_exp = tof_meas.numpy_intensity_sigma
    #
    #     cond_in = numpy.ones(np_time.shape, dtype=bool)
    #     try:
    #         range_ = self.range
    #         time_min = numpy.array(range_.time_min, dtype=float)
    #         time_max = numpy.array(range_.time_max, dtype=float)
    #
    #         cond_in = numpy.logical_and(cond_in, np_time >= time_min)
    #         cond_in = numpy.logical_and(cond_in, np_time <= time_max)
    #     except AttributeError:
    #         pass
    #
    #     np_time_in = np_time[cond_in]
    #     if flag_polarized:
    #         int_u_exp_in = int_u_exp[cond_in]
    #         sint_u_exp_in = sint_u_exp[cond_in]
    #         int_d_exp_in = int_d_exp[cond_in]
    #         sint_d_exp_in = sint_d_exp[cond_in]
    #     else:
    #         int_exp_in = int_exp[cond_in]
    #         sint_exp_in = sint_exp[cond_in]
    #
    #     tof_proc = self.calc_profile(
    #         np_time_in, l_crystal, flag_internal=flag_internal,
    #         flag_polarized=flag_polarized)
    #
    #     if flag_polarized:
    #         tof_proc.numpy_intensity_plus = int_u_exp_in
    #         tof_proc.numpy_intensity_plus_sigma = sint_u_exp_in
    #         tof_proc.numpy_intensity_minus = int_d_exp_in
    #         tof_proc.numpy_intensity_minus_sigma = sint_d_exp_in
    #         tof_proc.numpy_intensity = int_u_exp_in+int_d_exp_in
    #         tof_proc.numpy_intensity_sigma = numpy.sqrt(
    #             numpy.square(sint_u_exp_in) + numpy.square(sint_d_exp_in))
    #     else:
    #         tof_proc.numpy_intensity = int_exp_in
    #         tof_proc.numpy_intensity_sigma = sint_exp_in
    #
    #     int_u_mod = tof_proc.numpy_intensity_plus_total
    #     int_d_mod = tof_proc.numpy_intensity_minus_total
    #
    #     if flag_polarized:
    #         sint_sum_exp_in = (sint_u_exp_in**2 + sint_d_exp_in**2)**0.5
    #         chi_sq_u = ((int_u_mod-int_u_exp_in)/sint_u_exp_in)**2
    #         chi_sq_d = ((int_d_mod-int_d_exp_in)/sint_d_exp_in)**2
    #         chi_sq_sum = ((int_u_mod+int_d_mod-int_u_exp_in-int_d_exp_in) /
    #                       sint_sum_exp_in)**2
    #         chi_sq_dif = ((int_u_mod-int_d_mod-int_u_exp_in+int_d_exp_in) /
    #                       sint_sum_exp_in)**2
    #
    #         cond_u = numpy.logical_not(numpy.isnan(chi_sq_u))
    #         cond_d = numpy.logical_not(numpy.isnan(chi_sq_d))
    #         cond_sum = numpy.logical_not(numpy.isnan(chi_sq_sum))
    #         cond_dif = numpy.logical_not(numpy.isnan(chi_sq_dif))
    #     else:
    #         chi_sq_sum = ((int_u_mod+int_d_mod-int_exp_in)/sint_exp_in)**2
    #         cond_sum = numpy.logical_not(numpy.isnan(chi_sq_sum))
    #
    #     # exclude region
    #     try:
    #         exclude = self.exclude
    #         l_excl_time_min = exclude.time_low
    #         l_excl_time_max = exclude.time_high
    #         if flag_polarized:
    #             for excl_time_min, excl_time_max in zip(l_excl_time_min,
    #                                                     l_excl_time_max):
    #                 cond_1 = numpy.logical_or(np_time_in < 1.*excl_time_min,
    #                                           np_time_in > 1.*excl_time_max)
    #                 cond_u = numpy.logical_and(cond_u, cond_1)
    #                 cond_d = numpy.logical_and(cond_d, cond_1)
    #                 cond_sum = numpy.logical_and(cond_sum, cond_1)
    #         else:
    #             for excl_time_min, excl_time_max in zip(l_excl_time_min,
    #                                                     l_excl_time_max):
    #                 cond_1 = numpy.logical_or(np_time_in < 1.*excl_time_min,
    #                                           np_time_in > 1.*excl_time_max)
    #                 cond_sum = numpy.logical_and(cond_sum, cond_1)
    #     except AttributeError:
    #         pass
    #
    #     tof_proc.numpy_excluded = numpy.logical_not(cond_sum)
    #     chi_sq_sum_val = (chi_sq_sum[cond_sum]).sum()
    #     n_sum = cond_sum.sum()
    #
    #     if flag_polarized:
    #
    #         chi_sq_u_val = (chi_sq_u[cond_u]).sum()
    #         n_u = cond_u.sum()
    #
    #         chi_sq_d_val = (chi_sq_d[cond_d]).sum()
    #         n_d = cond_d.sum()
    #
    #         chi_sq_dif_val = (chi_sq_dif[cond_dif]).sum()
    #         n_dif = cond_dif.sum()
    #
    #         chi2 = self.chi2
    #
    #         flag_u = chi2.up
    #         flag_d = chi2.down
    #         flag_sum = chi2.sum
    #         flag_dif = chi2.diff
    #
    #         chi_sq_val = (int(flag_u)*chi_sq_u_val + int(flag_d)*chi_sq_d_val +
    #                       int(flag_sum)*chi_sq_sum_val +
    #                       int(flag_dif)*chi_sq_dif_val)
    #         n = (int(flag_u)*n_u + int(flag_d)*n_d + int(flag_sum)*n_sum +
    #              int(flag_dif)*n_dif)
    #     else:
    #         chi_sq_val = chi_sq_sum_val
    #         n = n_sum
    #
    #     if flag_internal:
    #         refine_ls = RefineLs(number_reflns=n,
    #                              goodness_of_fit_all=chi_sq_val/float(n),
    #                              weighting_scheme="sigma")
    #         self.refine_ls = refine_ls
    #         tof_proc.numpy_to_items()
    #     return chi_sq_val, n

    # def simulation(self, l_crystal, ttheta_start: float = 4.,
    #                ttheta_end: float = 120., ttheta_step: float = 0.1,
    #                flag_polarized: bool = True) -> NoReturn:
    #     """Simulate."""
    #     n_points = int(round((ttheta_end-ttheta_start)/ttheta_step))
    #     tth_in = numpy.linspace(ttheta_start, ttheta_end, n_points)
    #     if isinstance(l_crystal, Crystal):
    #         l_crystal = [l_crystal]
    #     self.calc_profile(tth_in, l_crystal, l_peak_in=None, l_refln_in=None,
    #                       l_refln_susceptibility_in=None, l_dd_in=None,
    #                       flag_internal=True, flag_polarized=flag_polarized)
    #     return
    #
    # def calc_iint(self, index_h, index_k, index_l, crystal: Crystal,
    #               flag_internal: bool = True):
    #     """Calculate the integrated intensity for h, k, l reflections.
    #
    #     Arguments
    #     ---------
    #         - h, k, l: 1D numpy array of Miller indexes, dtype = int32
    #         - l_crystal: a list of Crystal objects of cryspy library
    #         - flag_internal: a flag to calculate or to use internal objects.
    #                It should be True if user call the function.
    #                It's True by default.
    #
    #     Output
    #     ------
    #         - iint_u: 1D numpy array of integrated intensity up, dtype = float
    #         - iint_d: 1D numpy array of integrated intensity up, dtype = float
    #         - refln: ReflnL object of cryspy library (nuclear structure factor)
    #         - refln_s: ReflnSusceptibilityL object of cryspy library
    #           (nuclear structure factor)
    #     """
    #     index_hkl = numpy.stack([index_h, index_k, index_l], axis=0)
    #
    #     if flag_internal:
    #         try:
    #             d_internal_val = self.d_internal_val
    #         except AttributeError:
    #             d_internal_val = {}
    #             self.d_internal_val = d_internal_val
    #     else:
    #         d_internal_val = {}
    #         self.d_internal_val = d_internal_val
    #
    #     tof_parameters = self.tof_parameters
    #     try:
    #         field = tof_parameters.field
    #     except AttributeError:
    #         field = 0.
    #
    #     try:
    #         diffrn_radiation = self.diffrn_radiation
    #         p_u = float(diffrn_radiation.polarization)
    #         p_d = (2.*float(diffrn_radiation.efficiency)-1)*p_u
    #     except AttributeError:
    #         p_u = 0.0
    #         p_d = 0.0
    #
    #     try:
    #         if (not(flag_internal) | crystal.is_variables()):
    #             raise KeyError
    #         refln = d_internal_val[f"refln_{crystal.data_name:}"]
    #     except KeyError:
    #         refln = crystal.calc_refln(index_h, index_k, index_l,
    #                                    flag_internal=flag_internal)
    #         d_internal_val[f"refln_{crystal.data_name:}"] = refln
    #     f_nucl = refln.numpy_f_calc
    #     f_nucl_sq = abs(f_nucl*f_nucl.conjugate())
    #
    #     if isinstance(crystal, Crystal):
    #         try:
    #             if (not(flag_internal) | crystal.is_variables()):
    #                 raise KeyError
    #             refln_s = d_internal_val[
    #                 f"refln_susceptibility_{crystal.data_name:}"]
    #         except KeyError:
    #             refln_s = crystal.calc_refln_susceptibility(
    #                 index_h, index_k, index_l, flag_internal=flag_internal)
    #             d_internal_val[f"refln_susceptibility_{crystal.data_name:}"] \
    #                 = refln_s
    #
    #         sft_11 = refln_s.numpy_chi_11_calc
    #         sft_12 = refln_s.numpy_chi_12_calc
    #         sft_13 = refln_s.numpy_chi_13_calc
    #         sft_21 = refln_s.numpy_chi_21_calc
    #         sft_22 = refln_s.numpy_chi_22_calc
    #         sft_23 = refln_s.numpy_chi_23_calc
    #         sft_31 = refln_s.numpy_chi_31_calc
    #         sft_32 = refln_s.numpy_chi_32_calc
    #         sft_33 = refln_s.numpy_chi_33_calc
    #
    #
    #         _11, _12 = field*sft_11, field*sft_12
    #         _21, _13 = field*sft_21, field*sft_13
    #         _22, _23 = field*sft_22, field*sft_23
    #         _31, _32 = field*sft_31, field*sft_32
    #         _33 = field*sft_33
    #         _ij = numpy.stack([_11, _12, _13, _21, _22, _23, _31, _32, _33], axis=0)
    #         cell = crystal.cell
    #
    #         unit_cell_parameters = numpy.array([
    #             cell.length_a, cell.length_b, cell.length_c,
    #             cell.angle_alpha*numpy.pi/180.,
    #             cell.angle_beta*numpy.pi/180.,
    #             cell.angle_gamma*numpy.pi/180.], dtype=float)
    #
    #         t_ij = calc_matrix_t(index_hkl, unit_cell_parameters, flag_unit_cell_parameters=False)[0]
    #         #SIGMA = T CHI T^-1 = T chi T^T (because T is rotation matrix, therefore T^-1 = T^T)
    #         th_ij = calc_m1_m2_m1t(t_ij, _ij)[0]
    #         th_11, th_12, th_22 = th_ij[0], th_ij[1], th_ij[4]
    #
    #         # fm_p_sq = (field**2)*abs(0.5*(th_11*th_11.conjugate()+
    #         #            th_22*th_22.conjugate())+th_12*th_12.conjugate())
    #         # fm_p_field = field*0.5*(th_11+th_22)
    #         fm_p_sq = abs(0.5 * (th_11 * th_11.conjugate() +
    #                              th_22 * th_22.conjugate()) +
    #                       th_12 * th_12.conjugate())
    #         fm_p_field = 0.5*(th_11 + th_22)
    #         cross = 2.*(f_nucl.real*fm_p_field.real +
    #                     f_nucl.imag*fm_p_field.imag)
    #
    #         iint_u = f_nucl_sq + fm_p_sq + p_u*cross
    #         iint_d = f_nucl_sq + fm_p_sq - p_d*cross
    #
    #     elif isinstance(crystal, MagCrystal):
    #         try:
    #             if (not(flag_internal) | crystal.is_variables()):
    #                 raise KeyError
    #             f_mag_perp = d_internal_val[f"f_mag_perp_{crystal.data_name:}"]
    #
    #         except KeyError:
    #             f_mag_perp = crystal.calc_f_mag_perp(index_h, index_k, index_l)
    #             d_internal_val[f"f_mag_perp_{crystal.data_name:}"] = f_mag_perp
    #         f_mag_perp_sq = abs(f_mag_perp*f_mag_perp.conjugate()).sum(axis=0)
    #         iint_u = f_nucl_sq + f_mag_perp_sq
    #         iint_d = f_nucl_sq + f_mag_perp_sq
    #
    #     return iint_u, iint_d
    #
    # def _gauss_pd(self, time_2d):
    #     """One dimensional gauss powder diffraction."""
    #     ag, bg = self.ag, self.bg
    #     val_1 = bg*time_2d**2
    #     val_2 = numpy.where(val_1 < 5., numpy.exp(-val_1), 0.)
    #     self.gauss_pd = ag*val_2
    #
    # def _lor_pd(self, time_2d):
    #     """One dimensional lorentz powder diffraction."""
    #     al, bl = self.al, self.bl
    #     self.lor_pd = al*1./(1.+bl*time_2d**2)
    #
    # def calc_shape_profile(self, d, d_hkl, phase_igsize: float = 0.):
    #     """
    #     Calculate shape profile.
    #
    #     Calculate profile in the range of time for reflections placed on
    #     time_hkl with i_g parameter by default equal to zero
    #
    #     d, d_hkl in angstrems
    #     """
    #     tof_parameters = self.tof_parameters
    #     tof_profile = self.tof_profile
    #
    #     time = tof_parameters.calc_time_by_d(d)
    #     time_hkl = tof_parameters.calc_time_by_d(d_hkl)
    #
    #     # FIXME: strain_g, size_g, size_l
    #     np_shape_2d = tof_profile.calc_peak_shape_function(
    #         d, time, time_hkl, size_g=phase_igsize, strain_g=0.,
    #         size_l=0., strain_l=0.)
    #
    #     # Lorentz factor
    #     tth_rad = tof_parameters.ttheta_bank*numpy.pi/180.
    #     lorentz_factor = 1./(numpy.sin(tth_rad)*numpy.sin(0.5*tth_rad))
    #
    #     profile_2d = np_shape_2d*lorentz_factor
    #
    #     return profile_2d
    #
    # def params_to_cif(self, separator="_", flag: bool = False,
    #                   flag_minimal: bool = True) -> str:
    #     """Save parameters to cif format."""
    #     ls_out = []
    #     l_cls = (DiffrnRadiation, Chi2, Range, Extinction, Texture,
    #              TOFIntensityIncident, TOFParameters, TOFProfile, PhaseL,
    #              ExcludeL, TOFBackground)
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
    #     l_cls = (TOFMeasL, )
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
    #     l_cls = (RefineLs, ReflnL, ReflnSusceptibilityL, TOFPeakL, TOFProcL)
    #     l_obj = [item for item in self.items if type(item) in l_cls]
    #     l_itemn = [item for item in l_obj if isinstance(item, ItemN)]
    #     l_loopn = [item for item in l_obj if isinstance(item, LoopN)]
    #     ls_out.extend([_.to_cif(separator=separator)+"\n" for _ in l_itemn])
    #     ls_out.extend([_.to_cif(separator=separator)+"\n" for _ in l_loopn])
    #     return "\n".join(ls_out)

    def apply_constraints(self):
        """Apply constraints."""
        pass
        # if self.pd_meas is not None:
        #     self.pd_meas.apply_constraints()

    def plots(self):
        if self.is_attribute("tof_proc"):
            tof_proc = self.tof_proc
            fig_s, ax_s = tof_proc.plot_sum_diff()
            ax_s = fig_s.axes[0]
            ax_hkl = fig_s.axes[1]
            ax_s.set_title(self.data_name + " - "+ax_s.title.get_text())
            y_min_s, y_max_s = ax_hkl.get_ylim()
            y_dist_s = y_max_s-y_min_s
            y_step_s = 0.
            for item in self.items:
                if isinstance(item, TOFPeakL):
                    np_time = item.numpy_time
                    ax_hkl.plot(np_time, 0.*np_time+y_min_s-y_step_s, "|", label=item.loop_name)
                    y_step_s += 0.05*y_dist_s
            res = []
            ax_s.legend(loc='upper right')
            res.append((fig_s, ax_s))
            return res
        elif self.is_attribute("tof_meas"):
            return self.tof_meas.plots()

    def get_dictionary(self):
        self.form_object()
        dict_tof = {}

        phase, range_, tof_background = None, None, None
        exclude = None
        tof_meas, tof_parameters, tof_profile = None, None, None
        diffrn_radiation, texture = None, None
        setup = None
        tof_intensity_incident = None

        l_obj = take_items_by_class(self, (PhaseL, ))
        if len(l_obj) > 0:
            phase = l_obj[0]

        l_obj = take_items_by_class(self, (Range, ))
        if len(l_obj) > 0:
            range_ = l_obj[0]

        l_obj = take_items_by_class(self, (TOFBackground, ))
        if len(l_obj) > 0:
            tof_background = l_obj[0]

        l_obj = take_items_by_class(self, (ExcludeL, ))
        if len(l_obj) > 0:
            exclude = l_obj[0]

        l_obj = take_items_by_class(self, (TOFMeasL, ))
        if len(l_obj) > 0:
            tof_meas = l_obj[0]

        l_obj = take_items_by_class(self, (TOFParameters, ))
        if len(l_obj) > 0:
            tof_parameters = l_obj[0]

        l_obj = take_items_by_class(self, (TOFProfile, ))
        if len(l_obj) > 0:
            tof_profile = l_obj[0]

        l_obj = take_items_by_class(self, (DiffrnRadiation, ))
        if len(l_obj) > 0:
            diffrn_radiation = l_obj[0]

        l_obj = take_items_by_class(self, (TextureL, ))
        if len(l_obj) > 0:
            texture = l_obj[0]

        l_obj = take_items_by_class(self, (Setup, ))
        if len(l_obj) > 0:
            setup = l_obj[0]

        l_obj = take_items_by_class(self, (TOFIntensityIncident, ))
        if len(l_obj) > 0:
            tof_intensity_incident = l_obj[0]


        dict_tof["name"] = self.data_name
        dict_tof["type_name"] = self.get_name()

        if phase is not None:
            dict_tof["phase_name"] = numpy.array(phase.label, dtype=str)

            if phase.is_attribute("igsize"):
                dict_tof["phase_ig"] = numpy.array(phase.igsize, dtype=float)
                dict_tof["flags_phase_ig"] = numpy.array(phase.igsize_refinement, dtype=bool)
            else:
                dict_tof["phase_ig"] = numpy.zeros((len(phase.items),), dtype=float)
                dict_tof["flags_phase_ig"] = numpy.zeros((len(phase.items),), dtype=bool)

            if phase.is_attribute("scale"):
                dict_tof["phase_scale"] = numpy.array(phase.scale, dtype=float)
                dict_tof["flags_phase_scale"] = numpy.array(phase.scale_refinement, dtype=bool)
            else:
                dict_tof["phase_scale"] = numpy.zeros((len(phase.items),), dtype=float)
                dict_tof["flags_phase_scale"] = numpy.zeros((len(phase.items),), dtype=bool)

        if tof_background is not None:
            dict_tof["background_time_max"] = tof_background.time_max
            dict_tof["background_coefficients"] = tof_background.get_coefficients()
            dict_tof["flags_background_coefficients"] = tof_background.get_flags_coefficients()


        if tof_meas is not None:
            time = numpy.array(tof_meas.time, dtype=float)

            time_min = range_.time_min
            time_max = range_.time_max
            dict_tof["time_min"] = time_min
            dict_tof["time_max"] = time_max


            flag_in = numpy.logical_and(
                time >= time_min,
                time <= time_max)

            dict_tof["time"] = time[flag_in]
            if tof_meas.is_attribute("intensity_plus"):
                int_plus = numpy.array(tof_meas.intensity_plus, dtype=float)[flag_in]
                s_int_plus = numpy.array(tof_meas.intensity_plus_sigma, dtype=float)[flag_in]
                int_minus = numpy.array(tof_meas.intensity_minus, dtype=float)[flag_in]
                s_int_minus = numpy.array(tof_meas.intensity_minus_sigma, dtype=float)[flag_in]
                dict_tof["signal_exp_plus"] = numpy.stack([int_plus, s_int_plus], axis=0)
                dict_tof["signal_exp_minus"] = numpy.stack([int_minus, s_int_minus], axis=0)
            else:
                int_sum = numpy.array(tof_meas.intensity, dtype=float)[flag_in]
                s_int_sum = numpy.array(tof_meas.intensity_sigma, dtype=float)[flag_in]
                dict_tof["signal_exp"] = numpy.stack([int_sum, s_int_sum], axis=0)

            time_in_range = time[flag_in]
            flag_exclude = numpy.zeros(time_in_range.shape, dtype=bool)
            if exclude is not None:
                for item_e in exclude.items:
                    flag_in_1 = numpy.logical_and(
                        time_in_range >= item_e.time_min,
                        time_in_range <= item_e.time_max)
                    flag_exclude = numpy.logical_or(flag_exclude, flag_in_1)
            dict_tof["excluded_points"] = flag_exclude


        if tof_parameters is not None:
            dict_tof["ttheta_bank"] = tof_parameters.ttheta_bank * numpy.pi/180
            dict_tof["neutron_type"] = tof_parameters.neutrons
            dict_tof["zero"] = numpy.array([getattr(tof_parameters, "zero"), ], dtype=float)
            dict_tof["flags_zero"] = numpy.array([getattr(tof_parameters, "zero_refinement"), ], dtype=bool)
            dict_tof["dtt1"] = numpy.array([getattr(tof_parameters, "dtt1"), ], dtype=float)
            dict_tof["flags_dtt1"] = numpy.array([getattr(tof_parameters, "dtt1_refinement"), ], dtype=bool)
            if tof_parameters.is_attribute("dtt2"):
                dict_tof["dtt2"] = numpy.array([getattr(tof_parameters, "dtt2"), ], dtype=float)
                dict_tof["flags_dtt2"] = numpy.array([getattr(tof_parameters, "dtt2_refinement"), ], dtype=bool)
            if tof_parameters.is_attribute("zerot"):
                dict_tof["zerot"] = numpy.array([getattr(tof_parameters, "zerot"), ], dtype=float)
                dict_tof["flags_zerot"] = numpy.array([getattr(tof_parameters, "zerot_refinement"), ], dtype=bool)
            if tof_parameters.is_attribute("dtt1t"):
                dict_tof["dtt1t"] = numpy.array([getattr(tof_parameters, "dtt1t"), ], dtype=float)
                dict_tof["flags_dtt1t"] = numpy.array([getattr(tof_parameters, "dtt1t_refinement"), ], dtype=bool)
            if tof_parameters.is_attribute("dtt2t"):
                dict_tof["dtt2t"] = numpy.array([getattr(tof_parameters, "dtt2t"), ], dtype=float)
                dict_tof["flags_dtt2t"] = numpy.array([getattr(tof_parameters, "dtt2t_refinement"), ], dtype=bool)


        if tof_profile is not None:
            if tof_profile.is_attribute("peak_shape"):
                peak_shape = tof_profile.peak_shape
            else:
                peak_shape = "pseudo-Voigt"
            dict_tof["profile_peak_shape"] = peak_shape


            l_sigma = []
            l_sigma_refinement = []
            for numb in range(3):
                if tof_profile.is_attribute(f"sigma{numb:}"):
                    l_sigma.append(getattr(tof_profile, f"sigma{numb:}"))
                    l_sigma_refinement.append(getattr(tof_profile, f"sigma{numb:}_refinement"))
                else:
                    break
            dict_tof["profile_sigmas"] = numpy.array(l_sigma, dtype=float)
            dict_tof["flags_profile_sigmas"] = numpy.array(l_sigma_refinement, dtype=float)

            if peak_shape in ["pseudo-Voigt", "type0m"]:
                l_gamma = []
                l_gamma_refinement = []
                for numb in range(3):
                    if tof_profile.is_attribute(f"gamma{numb:}"):
                        l_gamma.append(getattr(tof_profile, f"gamma{numb:}"))
                        l_gamma_refinement.append(getattr(tof_profile, f"gamma{numb:}_refinement"))
                    else:
                        break
                dict_tof["profile_gammas"] = numpy.array(l_gamma, dtype=float)
                dict_tof["flags_profile_gammas"] = numpy.array(l_gamma_refinement, dtype=float)


            if peak_shape in ["pseudo-Voigt", "Gauss", ]:
                l_alpha = []
                l_alpha_refinement = []
                for numb in range(2):
                    if tof_profile.is_attribute(f"alpha{numb:}"):
                        l_alpha.append(getattr(tof_profile, f"alpha{numb:}"))
                        l_alpha_refinement.append(getattr(tof_profile, f"alpha{numb:}_refinement"))
                    else:
                        break

                l_beta = []
                l_beta_refinement = []
                for numb in range(2):
                    if tof_profile.is_attribute(f"beta{numb:}"):
                        l_beta.append(getattr(tof_profile, f"beta{numb:}"))
                        l_beta_refinement.append(getattr(tof_profile, f"beta{numb:}_refinement"))
                    else:
                        break

                dict_tof["profile_alphas"] = numpy.array(l_alpha, dtype=float)
                dict_tof["flags_profile_alphas"] = numpy.array(l_alpha_refinement, dtype=float)
                dict_tof["profile_betas"] = numpy.array(l_beta, dtype=float)
                dict_tof["flags_profile_betas"] = numpy.array(l_beta_refinement, dtype=float)
            elif peak_shape in ["type0m", ]:
                profile_alphas = numpy.array([
                    getattr(tof_profile, "alpha1"), getattr(tof_profile, "alpha2")
                    ], dtype=float)
                flags_profile_alphas = numpy.array([
                    getattr(tof_profile, "alpha1_refinement"), getattr(tof_profile, "alpha2_refinement")
                    ], dtype=bool)
                profile_betas = numpy.array([
                    getattr(tof_profile, "beta00"), getattr(tof_profile, "beta01"), getattr(tof_profile, "beta10")
                    ], dtype=float)
                flags_profile_betas = numpy.array([
                    getattr(tof_profile, "beta00_refinement"),
                    getattr(tof_profile, "beta01_refinement"),
                    getattr(tof_profile, "beta10_refinement")
                    ], dtype=bool)
                profile_rs = numpy.array([
                    getattr(tof_profile, "r01"), getattr(tof_profile, "r02"), getattr(tof_profile, "r03")
                    ], dtype=float)
                flags_profile_rs = numpy.array([
                    getattr(tof_profile, "r01_refinement"),
                    getattr(tof_profile, "r02_refinement"),
                    getattr(tof_profile, "r03_refinement")
                    ], dtype=bool)

                dict_tof["profile_alphas"] = profile_alphas
                dict_tof["flags_profile_alphas"] = flags_profile_alphas
                dict_tof["profile_betas"] = profile_betas
                dict_tof["flags_profile_betas"] = flags_profile_betas
                dict_tof["profile_rs"] = profile_rs
                dict_tof["flags_profile_rs"] = flags_profile_rs


        if diffrn_radiation is not None:
            beam_polarization = diffrn_radiation.polarization
            flipper_efficiency = diffrn_radiation.efficiency
            dict_tof["beam_polarization"] = numpy.array([beam_polarization], dtype=float)
            dict_tof["flipper_efficiency"] = numpy.array([flipper_efficiency], dtype=float)
            dict_tof["flags_beam_polarization"] = numpy.array([diffrn_radiation.polarization_refinement], dtype=bool)
            dict_tof["flags_flipper_efficiency"] = numpy.array([diffrn_radiation.efficiency_refinement], dtype=bool)


        if texture is not None:
            dict_tof["texture_name"] = numpy.array(texture.label, dtype=str)
            dict_tof["texture_g1"] = numpy.array(texture.g_1, dtype=float)
            dict_tof["texture_g2"] = numpy.array(texture.g_2, dtype=float)
            dict_tof["texture_axis"] = numpy.array(
                [texture.h_ax, texture.k_ax, texture.l_ax], dtype=float)
            dict_tof["flags_texture_g1"] = numpy.array(texture.g_1_refinement, dtype=bool)
            dict_tof["flags_texture_g2"] = numpy.array(texture.g_2_refinement, dtype=bool)
            dict_tof["flags_texture_axis"] = numpy.array(
                [texture.h_ax_refinement,
                 texture.k_ax_refinement,
                 texture.l_ax_refinement], dtype=bool)


        if setup is not None:
            dict_tof["magnetic_field"] = numpy.array([setup.field], dtype=float)

        if tof_intensity_incident is not None:
            dict_tof["spectrum_incident_type"] = tof_intensity_incident.spectrum
            dict_tof["spectrum_incident_coefficients"] = tof_intensity_incident.get_coefficients()
            dict_tof["flags_spectrum_incident_coefficients"] = tof_intensity_incident.get_flags_coefficients()

        return dict_tof


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

        if "phase_scale" in keys:
            hh = ddict_diffrn["phase_scale"]
            for i_item, item in enumerate(self.phase.items):
                item.scale = float(hh[i_item])

        if "phase_ig" in keys:
            hh = ddict_diffrn["phase_ig"]
            for i_item, item in enumerate(self.phase.items):
                item.igsize = float(hh[i_item])

        if "spectrum_incident_coefficients" in keys:
            tof_intensity_incident = self.tof_intensity_incident
            hh = ddict_diffrn["spectrum_incident_coefficients"]
            tof_intensity_incident.a0 = float(hh[0])
            tof_intensity_incident.a1 = float(hh[1])
            tof_intensity_incident.a2 = float(hh[2])
            tof_intensity_incident.a3 = float(hh[3])
            tof_intensity_incident.a4 = float(hh[4])
            tof_intensity_incident.a5 = float(hh[5])
            tof_intensity_incident.a6 = float(hh[6])
            tof_intensity_incident.a7 = float(hh[7])
            tof_intensity_incident.a8 = float(hh[8])

        if "zero" in keys:
            tof_parameters = self.tof_parameters
            setattr(tof_parameters, "zero", float(ddict_diffrn["zero"]))
            setattr(tof_parameters, "dtt1", float(ddict_diffrn["dtt1"]))
            if "dtt2" in keys:
                setattr(tof_parameters, "dtt2", float(ddict_diffrn["dtt2"]))
            if "zerot" in keys:
                setattr(tof_parameters, "zerot", float(ddict_diffrn["zerot"]))
            if "dtt1t" in keys:
                setattr(tof_parameters, "dtt1t", float(ddict_diffrn["dtt1t"]))
            if "dtt2t" in keys:
                setattr(tof_parameters, "dtt2t", float(ddict_diffrn["dtt2t"]))

        if "background_coefficients" in keys:
            tof_background = self.tof_background
            hh = ddict_diffrn["background_coefficients"]
            for i_hh in range(hh.size):
                setattr(tof_background, f"coeff{i_hh+1:}", float(hh[i_hh]))


        if len(parameter_label) > 0:
            for name, sigma in zip(l_parameter_name, l_sigma):
                if name[0] == "phase_scale":
                    self.phase.items[name[1][0]].scale_sigma = float(sigma)
                if name[0] == "phase_ig":
                    self.phase.items[name[1][0]].igsize_sigma = float(sigma)
                if name[0] == "beam_polarization":
                    self.diffrn_radiation.polarization_sigma = float(sigma)
                if name[0] == "flipper_efficiency":
                    self.diffrn_radiation.efficiency_sigma = float(sigma)
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

        if (("signal_plus" in keys) and ("signal_minus" in keys)):
            tof_proc = TOFProcL()
            tof_proc.numpy_time = numpy.round(ddict_diffrn["time"], decimals=5)
            tof_proc.numpy_intensity_plus_net = numpy.round(ddict_diffrn["signal_plus"], decimals=5)
            tof_proc.numpy_intensity_minus_net = numpy.round(ddict_diffrn["signal_minus"], decimals=5)
            if "signal_exp_plus" in keys:
                tof_proc.numpy_intensity_plus = numpy.round(ddict_diffrn["signal_exp_plus"][0], decimals=5)
                tof_proc.numpy_intensity_plus_sigma = numpy.round(ddict_diffrn["signal_exp_plus"][1], decimals=5)
                tof_proc.numpy_intensity_minus = numpy.round(ddict_diffrn["signal_exp_minus"][0], decimals=5)
                tof_proc.numpy_intensity_minus_sigma = numpy.round(ddict_diffrn["signal_exp_minus"][1], decimals=5)
            else:
                tof_proc.numpy_intensity = numpy.round(ddict_diffrn["signal_exp"][0], decimals=5)
                tof_proc.numpy_intensity_sigma = numpy.round(ddict_diffrn["signal_exp"][1], decimals=5)
            tof_proc.numpy_intensity_bkg_calc = numpy.round(ddict_diffrn["signal_background"], decimals=5)
            tof_proc.numpy_excluded = ddict_diffrn["excluded_points"]
            tof_proc.numpy_to_items()
            self.tof_proc = tof_proc

        l_tof_peak = []
        l_refln = []
        for item_phase in self.phase.items:
            s_label = f"dict_in_out_{item_phase.label:}"
            if s_label in keys:
                dict_crystal = ddict_diffrn[s_label]
                dict_crystal_keys = dict_crystal.keys()
                if (("index_hkl" in dict_crystal_keys) and  ("time_hkl" in dict_crystal_keys)):
                    index_hkl = dict_crystal["index_hkl"]
                    time_hkl = dict_crystal["time_hkl"]

                    tof_peak = TOFPeakL(loop_name = item_phase.label)
                    int_plus_max, int_minus_max = None, None
                    if "iint_plus_with_factors" in dict_crystal_keys:
                        # int_plus_max = dict_crystal["iint_plus_with_factors"]
                        int_plus_max = dict_crystal["iint_plus"]

                    if "iint_minus_with_factors" in dict_crystal_keys:
                        # int_minus_max = dict_crystal["iint_minus_with_factors"]
                        int_minus_max = dict_crystal["iint_minus"]

                    if "f_nucl":
                        refln = ReflnL(loop_name = item_phase.label)
                        refln.numpy_index_h = index_hkl[0]
                        refln.numpy_index_k = index_hkl[1]
                        refln.numpy_index_l = index_hkl[2]
                        refln.numpy_a_calc = dict_crystal["f_nucl"].real
                        refln.numpy_b_calc = dict_crystal["f_nucl"].imag
                        refln.numpy_to_items()
                        l_refln.append(refln)

                    if not((int_plus_max is None) and (int_minus_max is None)):
                        tof_peak.numpy_index_h = index_hkl[0]
                        tof_peak.numpy_index_k = index_hkl[1]
                        tof_peak.numpy_index_l = index_hkl[2]
                        tof_peak.numpy_index_mult = dict_crystal["multiplicity_hkl"]
                        tof_peak.numpy_time = numpy.round(time_hkl, decimals=3)
                        tof_peak.numpy_intensity_plus = int_plus_max
                        tof_peak.numpy_intensity_minus = int_minus_max
                        tof_peak.numpy_to_items()
                        l_tof_peak.append(tof_peak)

        if len(l_tof_peak) > 0:
            self.add_items(l_tof_peak)
        if len(l_refln) > 0:
            self.add_items(l_refln)

    # def estimate_background(self):
    #     tof_background = self.tof_background
    #     tof_proc = self.tof_proc
    #     tof_proc.estimate_background(tof_background)

