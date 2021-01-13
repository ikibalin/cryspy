"""Description of Pd class."""
__author__ = 'ikibalin'
__version__ = "2020_08_19"

import numpy
from typing import NoReturn

from cryspy.A_functions_base.function_2_crystallography_base import \
    calc_cos_ang
from cryspy.A_functions_base.function_1_matrices import calc_mRmCmRT

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN
from cryspy.B_parent_classes.cl_3_data import DataN

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
from cryspy.C_item_loop_classes.cl_1_texture import Texture
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
                        TOFPeakL, RefineLs, ReflnL, ReflnSusceptibilityL)
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

    def calc_profile(self, time, l_crystal, flag_internal: bool = True,
                     flag_polarized: bool = True):
        """Calculate intensity for the given diffraction angles.

        Arguments
        ---------
            - time: 1D numpy array of time in micro seconds (???)
            - l_crystal: a list of Crystal objects of cryspy library
            - l_peak_in: precalculated data about integrated intensities
              (used to speed up the calculations)
            - l_refln_in: precalculated data about nuclear structure factors
              (used to speed up the calculations)
            - l_refln_susceptibility_in: precalculated data about
              (used to speed up the calculations)
            - flag_internal: a flag to calculate or to use internal objects.
                   It should be True if user call the function.
                   It's True by default.

        Output
        ------
            - proc: output profile
            - l_peak: data about peaks
            - l_refln: data about nuclear structure factor
        """
        if flag_internal:
            d_internal_val = {}
            self.d_internal_val = d_internal_val
        else:
            try:
                d_internal_val = self.d_internal_val
            except AttributeError:
                d_internal_val = {}
                self.d_internal_val = d_internal_val

        tof_proc = TOFProcL()
        tof_proc.numpy_time = time

        try:
            tof_background = self.tof_background
            int_bkgd = tof_background.calc_background(time)
        except AttributeError:
            int_bkgd = 0.*time

        tof_proc.numpy_intensity_bkg_calc = int_bkgd

        tof_parameters = self.tof_parameters
        # tof_profile = self.tof_profile
        
        d = tof_parameters.calc_d_by_time(time)
        d_min, d_max = tof_parameters.calc_d_min_max(time)
        sthovl_min = 0.5/d_max
        sthovl_max = 0.5/d_min

        try:
            texture = self.texture
            h_ax, k_ax, l_ax = texture.h_ax, texture.k_ax, texture.l_ax
            g_1, g_2 = texture.g_1, texture.g_2
        except AttributeError:
            texture = None

        res_u_1d = numpy.zeros(time.shape[0], dtype=float)
        res_d_1d = numpy.zeros(time.shape[0], dtype=float)

        phase = self.phase

        for phase_item in phase.items:
            phase_label = phase_item.label
            phase_scale = phase_item.scale
            try:
                phase_igsize = phase_item.igsize
                if phase_igsize is None: phase_igsize = 0. # temporary solution
            except AttributeError:
                phase_igsize = 0.

            for i_crystal, crystal in enumerate(l_crystal):
                if crystal.data_name.lower() == phase_label.lower():
                    ind_cry = i_crystal
                    break
            if ind_cry is None:
                raise AttributeError(
                    f"Crystal with name '{crystal.data_name:}' is not found.")
                return
            crystal = l_crystal[ind_cry]

            scale = phase_scale

            try:
                peak = d_internal_val[f"peak_{crystal.data_name:}"]
                index_h = peak.numpy_index_h
                index_k = peak.numpy_index_k
                index_l = peak.numpy_index_l
                mult = peak.numpy_index_multiplicity
                cond_2 = True
            except KeyError:
                if texture is None:
                    index_h, index_k, index_l, mult = crystal.calc_hkl(
                        sthovl_min, sthovl_max)
                else:
                    index_h, index_k, index_l, mult = \
                        crystal.calc_hkl_in_range(sthovl_min, sthovl_max)
                peak = TOFPeakL(loop_name=phase_label)
                peak.numpy_index_h = numpy.array(index_h, dtype=int)
                peak.numpy_index_k = numpy.array(index_k, dtype=int)
                peak.numpy_index_l = numpy.array(index_l, dtype=int)
                peak.numpy_index_multiplicity = numpy.array(mult, dtype=int)
                d_internal_val[f"peak_{crystal.data_name:}"] = peak
                cond_2 = False

            cond_1 = not(crystal.is_variables())
            if cond_1 & cond_2:
                np_iint_u = peak.numpy_intensity_up
                np_iint_d = peak.numpy_intensity_down
            else:
                np_iint_u, np_iint_d = self.calc_iint(
                    index_h, index_k, index_l, crystal,
                    flag_internal=flag_internal)
                peak.numpy_intensity_up = np_iint_u
                peak.numpy_intensity_down = np_iint_d

            cell = crystal.cell
            sthovl_hkl = cell.calc_sthovl(index_h, index_k, index_l)
            d_hkl = 0.5/sthovl_hkl
            time_hkl = tof_parameters.calc_time_by_d(d_hkl)
            

            profile_2d = self.calc_shape_profile(
                d, d_hkl, phase_igsize=phase_igsize)

            peak.numpy_time = time_hkl

            tth_rad = tof_parameters.ttheta_bank * numpy.pi/180.
            wavelength = 2.*d_hkl*numpy.sin(0.5*tth_rad)
            wavelength_4 = wavelength**4
            
            iint_u_0 = np_iint_u * mult * wavelength_4
            iint_d_0 = np_iint_d * mult * wavelength_4

            if tof_parameters.is_attribute("extinction"):
                tof_extinction = tof_parameters.extinction
                iint_u_0 = (1. - tof_extinction*iint_u_0)*mult*wavelength_4
                iint_d_0 = (1. - tof_extinction*iint_d_0)*mult*wavelength_4

            np_iint_u_2d = numpy.meshgrid(time, iint_u_0, indexing="ij")[1]
            np_iint_d_2d = numpy.meshgrid(time, iint_d_0, indexing="ij")[1]

            # texture
            if texture is not None:
                cond_2 = ((not(texture.is_variables())) &
                          (not(cell.is_variables())))
                if cond_2:
                    try:
                        texture_2d = d_internal_val[
                            f"texture_2d_{crystal.data_name:}"]
                        cond_1 = False
                    except KeyError:
                        cond_1 = True
                if cond_1:
                    cos_alpha_ax = calc_cos_ang(cell, h_ax, k_ax, l_ax,
                                                index_h, index_k, index_l)
                    c_help = 1.-cos_alpha_ax**2
                    c_help[c_help < 0.] = 0.
                    sin_alpha_ax = numpy.sqrt(c_help)
                    cos_alpha_ax_2d = numpy.meshgrid(time, cos_alpha_ax,
                                                     indexing="ij")[1]
                    sin_alpha_ax_2d = numpy.meshgrid(time, sin_alpha_ax,
                                                     indexing="ij")[1]
                    # cos_alpha_2d = cos_alpha_ax_2d*cos_alpha_ang_2d + \
                    #                sin_alpha_ax_2d*sin_alpha_ang_2d
                    # cos_alpha_ang_2d, sin_alpha_ang_2d
                    cos_alpha_2d = cos_alpha_ax_2d*1.+sin_alpha_ax_2d*0.
                    texture_2d = g_2 + (1. - g_2) * (1./g_1 + (g_1**2 - 1./g_1)
                                                     * cos_alpha_2d**2)**(-1.5)
                    d_internal_val[f"texture_2d_{crystal.data_name:}"] = \
                        texture_2d

                profile_2d = profile_2d*texture_2d

            res_u_2d = profile_2d*np_iint_u_2d
            res_d_2d = profile_2d*np_iint_d_2d

            # 0.5 to have the same meaning for scale factor as in FullProf
            res_u_1d += 0.5*scale*res_u_2d.sum(axis=1)
            res_d_1d += 0.5*scale*res_d_2d.sum(axis=1)

            if flag_internal:
                peak.numpy_to_items()

        tof_proc.numpy_d_spacing = d
        tof_proc.numpy_intensity_up_net = res_u_1d
        tof_proc.numpy_intensity_down_net = res_d_1d
        tof_proc.numpy_intensity_net = res_u_1d+res_d_1d
        tof_proc.numpy_intensity_diff_total = res_u_1d-res_d_1d
        if flag_polarized:
            tof_proc.numpy_intensity_up_total = res_u_1d+int_bkgd
            tof_proc.numpy_intensity_down_total = res_d_1d+int_bkgd
            tof_proc.numpy_intensity_total = res_u_1d+res_d_1d+int_bkgd+int_bkgd
        else:
            tof_proc.numpy_intensity_up_total = res_u_1d+0.5*int_bkgd
            tof_proc.numpy_intensity_down_total = res_d_1d+0.5*int_bkgd
            tof_proc.numpy_intensity_total = res_u_1d+res_d_1d+int_bkgd

        if flag_internal:
            tof_proc.numpy_to_items()
            l_calc_objs = []
            for crystal in l_crystal:
                try:
                    obj = d_internal_val[f"refln_{crystal.data_name}"]
                    l_calc_objs.append(obj)
                except KeyError:
                    pass
                try:
                    obj = d_internal_val[
                        f"refln_susceptibility_{crystal.data_name}"]
                    l_calc_objs.append(obj)
                except KeyError:
                    pass
                try:
                    obj = d_internal_val[f"peak_{crystal.data_name}"]
                    l_calc_objs.append(obj)
                except KeyError:
                    pass
            l_calc_objs.append(tof_proc)
            self.add_items(l_calc_objs)
        return tof_proc

    def calc_chi_sq(self, l_crystal, flag_internal=True):
        """
        Calculate chi square.

        Arguments
        ---------
            - l_crystal: a list of Crystal objects of cryspy library
            - flag_internal: a flag to calculate or to use internal objects.
                   It should be True if user call the function.
                   It's True by default.

        Output
        ------
            - chi_sq_val: chi square of flip ratio
              (Sum_i ((y_e_i - y_m_i) / sigma_i)**2)
            - n: number of measured reflections
        """
        tof_meas = self.tof_meas
        flag_polarized = tof_meas.is_polarized()

        np_time = tof_meas.numpy_time
        if flag_polarized:
            int_u_exp = tof_meas.numpy_intensity_up
            sint_u_exp = tof_meas.numpy_intensity_up_sigma
            int_d_exp = tof_meas.numpy_intensity_down
            sint_d_exp = tof_meas.numpy_intensity_down_sigma
        else:
            int_exp = tof_meas.numpy_intensity
            sint_exp = tof_meas.numpy_intensity_sigma

        cond_in = numpy.ones(np_time.shape, dtype=bool)
        try:
            range_ = self.range
            time_min = numpy.array(range_.time_min, dtype=float)
            time_max = numpy.array(range_.time_max, dtype=float)

            cond_in = numpy.logical_and(cond_in, np_time >= time_min)
            cond_in = numpy.logical_and(cond_in, np_time <= time_max)
        except AttributeError:
            pass

        np_time_in = np_time[cond_in]
        if flag_polarized:
            int_u_exp_in = int_u_exp[cond_in]
            sint_u_exp_in = sint_u_exp[cond_in]
            int_d_exp_in = int_d_exp[cond_in]
            sint_d_exp_in = sint_d_exp[cond_in]
        else:
            int_exp_in = int_exp[cond_in]
            sint_exp_in = sint_exp[cond_in]

        tof_proc = self.calc_profile(
            np_time_in, l_crystal, flag_internal=flag_internal,
            flag_polarized=flag_polarized)

        if flag_polarized:
            tof_proc.numpy_intensity_up = int_u_exp_in
            tof_proc.numpy_intensity_up_sigma = sint_u_exp_in
            tof_proc.numpy_intensity_down = int_d_exp_in
            tof_proc.numpy_intensity_down_sigma = sint_d_exp_in
            tof_proc.numpy_intensity = int_u_exp_in+int_d_exp_in
            tof_proc.numpy_intensity_sigma = numpy.sqrt(
                numpy.square(sint_u_exp_in) + numpy.square(sint_d_exp_in))
        else:
            tof_proc.numpy_intensity = int_exp_in
            tof_proc.numpy_intensity_sigma = sint_exp_in

        int_u_mod = tof_proc.numpy_intensity_up_total
        int_d_mod = tof_proc.numpy_intensity_down_total

        if flag_polarized:
            sint_sum_exp_in = (sint_u_exp_in**2 + sint_d_exp_in**2)**0.5
            chi_sq_u = ((int_u_mod-int_u_exp_in)/sint_u_exp_in)**2
            chi_sq_d = ((int_d_mod-int_d_exp_in)/sint_d_exp_in)**2
            chi_sq_sum = ((int_u_mod+int_d_mod-int_u_exp_in-int_d_exp_in) /
                          sint_sum_exp_in)**2
            chi_sq_dif = ((int_u_mod-int_d_mod-int_u_exp_in+int_d_exp_in) /
                          sint_sum_exp_in)**2

            cond_u = numpy.logical_not(numpy.isnan(chi_sq_u))
            cond_d = numpy.logical_not(numpy.isnan(chi_sq_d))
            cond_sum = numpy.logical_not(numpy.isnan(chi_sq_sum))
            cond_dif = numpy.logical_not(numpy.isnan(chi_sq_dif))
        else:
            chi_sq_sum = ((int_u_mod+int_d_mod-int_exp_in)/sint_exp_in)**2
            cond_sum = numpy.logical_not(numpy.isnan(chi_sq_sum))

        # exclude region
        try:
            exclude = self.exclude
            l_excl_time_min = exclude.time_low
            l_excl_time_max = exclude.time_high
            if flag_polarized:
                for excl_time_min, excl_time_max in zip(l_excl_time_min,
                                                        l_excl_time_max):
                    cond_1 = numpy.logical_or(np_time_in < 1.*excl_time_min,
                                              np_time_in > 1.*excl_time_max)
                    cond_u = numpy.logical_and(cond_u, cond_1)
                    cond_d = numpy.logical_and(cond_d, cond_1)
                    cond_sum = numpy.logical_and(cond_sum, cond_1)
            else:
                for excl_time_min, excl_time_max in zip(l_excl_time_min,
                                                        l_excl_time_max):
                    cond_1 = numpy.logical_or(np_time_in < 1.*excl_time_min,
                                              np_time_in > 1.*excl_time_max)
                    cond_sum = numpy.logical_and(cond_sum, cond_1)
        except AttributeError:
            pass

        tof_proc.numpy_excluded = numpy.logical_not(cond_sum)
        chi_sq_sum_val = (chi_sq_sum[cond_sum]).sum()
        n_sum = cond_sum.sum()

        if flag_polarized:

            chi_sq_u_val = (chi_sq_u[cond_u]).sum()
            n_u = cond_u.sum()

            chi_sq_d_val = (chi_sq_d[cond_d]).sum()
            n_d = cond_d.sum()

            chi_sq_dif_val = (chi_sq_dif[cond_dif]).sum()
            n_dif = cond_dif.sum()

            chi2 = self.chi2

            flag_u = chi2.up
            flag_d = chi2.down
            flag_sum = chi2.sum
            flag_dif = chi2.diff

            chi_sq_val = (int(flag_u)*chi_sq_u_val + int(flag_d)*chi_sq_d_val +
                          int(flag_sum)*chi_sq_sum_val +
                          int(flag_dif)*chi_sq_dif_val)
            n = (int(flag_u)*n_u + int(flag_d)*n_d + int(flag_sum)*n_sum +
                 int(flag_dif)*n_dif)
        else:
            chi_sq_val = chi_sq_sum_val
            n = n_sum

        if flag_internal:
            refine_ls = RefineLs(number_reflns=n,
                                 goodness_of_fit_all=chi_sq_val/float(n),
                                 weighting_scheme="sigma")
            self.refine_ls = refine_ls
            tof_proc.numpy_to_items()
        return chi_sq_val, n

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

    def calc_iint(self, index_h, index_k, index_l, crystal: Crystal,
                  flag_internal: bool = True):
        """Calculate the integrated intensity for h, k, l reflections.

        Arguments
        ---------
            - h, k, l: 1D numpy array of Miller indexes, dtype = int32
            - l_crystal: a list of Crystal objects of cryspy library
            - flag_internal: a flag to calculate or to use internal objects.
                   It should be True if user call the function.
                   It's True by default.

        Output
        ------
            - iint_u: 1D numpy array of integrated intensity up, dtype = float
            - iint_d: 1D numpy array of integrated intensity up, dtype = float
            - refln: ReflnL object of cryspy library (nuclear structure factor)
            - refln_s: ReflnSusceptibilityL object of cryspy library
              (nuclear structure factor)
        """
        if flag_internal:
            try:
                d_internal_val = self.d_internal_val
            except AttributeError:
                d_internal_val = {}
                self.d_internal_val = d_internal_val
        else:
            d_internal_val = {}
            self.d_internal_val = d_internal_val

        tof_parameters = self.tof_parameters
        try:
            field = tof_parameters.field
        except AttributeError:
            field = 0.

        try:
            diffrn_radiation = self.diffrn_radiation
            p_u = float(diffrn_radiation.polarization)
            p_d = (2.*float(diffrn_radiation.efficiency)-1)*p_u
        except AttributeError:
            p_u = 0.0
            p_d = 0.0

        try:
            if (not(flag_internal) | crystal.is_variables()):
                raise KeyError
            refln = d_internal_val[f"refln_{crystal.data_name:}"]
        except KeyError:
            refln = crystal.calc_refln(index_h, index_k, index_l,
                                       flag_internal=flag_internal)
            d_internal_val[f"refln_{crystal.data_name:}"] = refln
        f_nucl = refln.numpy_f_calc
        f_nucl_sq = abs(f_nucl*f_nucl.conjugate())

        if isinstance(crystal, Crystal):
            try:
                if (not(flag_internal) | crystal.is_variables()):
                    raise KeyError
                refln_s = d_internal_val[
                    f"refln_susceptibility_{crystal.data_name:}"]
            except KeyError:
                refln_s = crystal.calc_refln_susceptibility(
                    index_h, index_k, index_l, flag_internal=flag_internal)
                d_internal_val[f"refln_susceptibility_{crystal.data_name:}"] \
                    = refln_s

            sft_11 = refln_s.numpy_chi_11_calc
            sft_12 = refln_s.numpy_chi_12_calc
            sft_13 = refln_s.numpy_chi_13_calc
            sft_21 = refln_s.numpy_chi_21_calc
            sft_22 = refln_s.numpy_chi_22_calc
            sft_23 = refln_s.numpy_chi_23_calc
            sft_31 = refln_s.numpy_chi_31_calc
            sft_32 = refln_s.numpy_chi_32_calc
            sft_33 = refln_s.numpy_chi_33_calc

            sftm_11 = refln_s.numpy_moment_11_calc
            sftm_12 = refln_s.numpy_moment_12_calc
            sftm_13 = refln_s.numpy_moment_13_calc
            sftm_21 = refln_s.numpy_moment_21_calc
            sftm_22 = refln_s.numpy_moment_22_calc
            sftm_23 = refln_s.numpy_moment_23_calc
            sftm_31 = refln_s.numpy_moment_31_calc
            sftm_32 = refln_s.numpy_moment_32_calc
            sftm_33 = refln_s.numpy_moment_33_calc

            _11, _12 = sftm_11+field*sft_11, sftm_12+field*sft_12
            _21, _13 = sftm_21+field*sft_21, sftm_13+field*sft_13
            _22, _23 = sftm_22+field*sft_22, sftm_23+field*sft_23
            _31, _32 = sftm_31+field*sft_31, sftm_32+field*sft_32
            _33 = sftm_33+field*sft_33
            _ij = (_11, _12, _13, _21, _22, _23, _31, _32, _33)
            cell = crystal.cell

            # k_loc = cell.calc_k_loc(h, k, l)
            t_ij = cell.calc_m_t(index_h, index_k, index_l)
            # FIXME: I would like to recheck the expression for T
            #        and expression SIGMA = T^T CHI T
            t_tr_ij = (t_ij[0], t_ij[3], t_ij[6],
                       t_ij[1], t_ij[4], t_ij[7],
                       t_ij[2], t_ij[5], t_ij[8])
            th_11, th_12, th_13, th_21, th_22, th_23, th_31, th_32, th_33 = \
                calc_mRmCmRT(t_tr_ij, _ij)

            # fm_p_sq = (field**2)*abs(0.5*(th_11*th_11.conjugate()+
            #            th_22*th_22.conjugate())+th_12*th_12.conjugate())
            # fm_p_field = field*0.5*(th_11+th_22)
            fm_p_sq = abs(0.5 * (th_11 * th_11.conjugate() +
                                 th_22 * th_22.conjugate()) +
                          th_12 * th_12.conjugate())
            fm_p_field = 0.5*(th_11 + th_22)
            cross = 2.*(f_nucl.real*fm_p_field.real +
                        f_nucl.imag*fm_p_field.imag)

            iint_u = f_nucl_sq + fm_p_sq + p_u*cross
            iint_d = f_nucl_sq + fm_p_sq - p_d*cross

        elif isinstance(crystal, MagCrystal):
            try:
                if (not(flag_internal) | crystal.is_variables()):
                    raise KeyError
                f_mag_perp = d_internal_val[f"f_mag_perp_{crystal.data_name:}"]

            except KeyError:
                f_mag_perp = crystal.calc_f_mag_perp(index_h, index_k, index_l)
                d_internal_val[f"f_mag_perp_{crystal.data_name:}"] = f_mag_perp
            f_mag_perp_sq = abs(f_mag_perp*f_mag_perp.conjugate()).sum(axis=0)
            iint_u = f_nucl_sq + f_mag_perp_sq
            iint_d = f_nucl_sq + f_mag_perp_sq

        return iint_u, iint_d

    def _gauss_pd(self, time_2d):
        """One dimensional gauss powder diffraction."""
        ag, bg = self.ag, self.bg
        val_1 = bg*time_2d**2
        val_2 = numpy.where(val_1 < 5., numpy.exp(-val_1), 0.)
        self.gauss_pd = ag*val_2

    def _lor_pd(self, time_2d):
        """One dimensional lorentz powder diffraction."""
        al, bl = self.al, self.bl
        self.lor_pd = al*1./(1.+bl*time_2d**2)

    def calc_shape_profile(self, d, d_hkl, phase_igsize: float = 0.):
        """
        Calculate shape profile.

        Calculate profile in the range of time for reflections placed on
        time_hkl with i_g parameter by default equal to zero

        d, d_hkl in angstrems
        """
        tof_parameters = self.tof_parameters
        tof_profile = self.tof_profile
        
        time = tof_parameters.calc_time_by_d(d)
        time_hkl = tof_parameters.calc_time_by_d(d_hkl)

        # FIXME: strain_g, size_g, size_l
        np_shape_2d = tof_profile.calc_peak_shape_function(
            d, time, time_hkl, size_g=phase_igsize, strain_g=0.,
            size_l=0., strain_l=0.)

        # Lorentz factor
        tth_rad = tof_parameters.ttheta_bank*numpy.pi/180.
        lorentz_factor = 1./(numpy.sin(tth_rad)*numpy.sin(0.5*tth_rad))

        profile_2d = np_shape_2d*lorentz_factor

        return profile_2d

    def params_to_cif(self, separator="_", flag: bool = False,
                      flag_minimal: bool = True) -> str:
        """Save parameters to cif format."""
        ls_out = []
        l_cls = (DiffrnRadiation, Chi2, Range, Extinction, Texture,
                 TOFIntensityIncident, TOFParameters, TOFProfile, PhaseL,
                 ExcludeL, TOFBackground)
        l_obj = [item for item in self.items if type(item) in l_cls]
        l_itemn = [item for item in l_obj if isinstance(item, ItemN)]
        l_loopn = [item for item in l_obj if isinstance(item, LoopN)]
        ls_out.extend([_.to_cif(separator=separator)+"\n" for _ in l_itemn])
        ls_out.extend([_.to_cif(separator=separator)+"\n" for _ in l_loopn])
        return "\n".join(ls_out)

    def data_to_cif(self, separator="_", flag: bool = False,
                    flag_minimal: bool = True) -> str:
        """Save data to cif format."""
        ls_out = []
        l_cls = (TOFMeasL, )
        l_obj = [item for item in self.items if type(item) in l_cls]
        l_itemn = [item for item in l_obj if isinstance(item, ItemN)]
        l_loopn = [item for item in l_obj if isinstance(item, LoopN)]
        ls_out.extend([_.to_cif(separator=separator)+"\n" for _ in l_itemn])
        ls_out.extend([_.to_cif(separator=separator)+"\n" for _ in l_loopn])
        return "\n".join(ls_out)

    def calc_to_cif(self, separator="_", flag=False, flag_minimal=True) -> str:
        """Save calculations to cif format."""
        ls_out = []
        l_cls = (RefineLs, ReflnL, ReflnSusceptibilityL, TOFPeakL, TOFProcL)
        l_obj = [item for item in self.items if type(item) in l_cls]
        l_itemn = [item for item in l_obj if isinstance(item, ItemN)]
        l_loopn = [item for item in l_obj if isinstance(item, LoopN)]
        ls_out.extend([_.to_cif(separator=separator)+"\n" for _ in l_itemn])
        ls_out.extend([_.to_cif(separator=separator)+"\n" for _ in l_loopn])
        return "\n".join(ls_out)

    def apply_constraints(self):
        """Apply constraints."""
        if self.pd_meas is not None:
            self.pd_meas.apply_constraints()

    def plots(self):
        if self.is_attribute("tof_proc"):
            tof_proc = self.tof_proc
            fig_s, ax_s = tof_proc.plot_sum()
            fig_d_ax_d = tof_proc.plot_diff()
            flag_d = fig_d_ax_d is not None
            if flag_d:
                fig_d, ax_d = fig_d_ax_d
                y_min_d, y_max_d = ax_d.get_ylim()
                y_dist_d = y_max_d-y_min_d
                y_step_d = 0.

            y_min_s, y_max_s = ax_s.get_ylim()
            y_dist_s = y_max_s-y_min_s
            y_step_s = 0.
            for item in self.items:
                if isinstance(item, TOFPeakL):
                    np_tth = item.numpy_time
                    ax_s.plot(np_tth, 0.*np_tth+y_min_s-y_step_s, "|", label=item.loop_name)
                    y_step_s += 0.05*y_dist_s
                    if flag_d:
                        ax_d.plot(np_tth, 0.*np_tth+y_min_d-y_step_d, "|", label=item.loop_name)
                        y_step_d += 0.05*y_dist_d
            ax_s.set_title(self.data_name)
            ax_s.legend(loc='upper right')
            if flag_d:
                ax_d.set_title(self.data_name)
                ax_d.legend(loc='upper right')
            if flag_d:
                res = [(fig_s, ax_s), (fig_d, ax_d)]
            else:
                res = [(fig_s, ax_s)]
            return res
            
            return self.tof_proc.plots()
        elif self.is_attribute("tof_meas"):
            return self.tof_meas.plots()

# s_cont = """
#   data_tof

# _tof_profile_sigma2 0.00000
# _tof_profile_sigma1 61.67600
# _tof_profile_sigma0 6.19600
# _tof_profile_gamma2 0.00000
# _tof_profile_gamma1 0.60400
# _tof_profile_gamma0 0.00000

# _range_time_min 3001.589
# _range_time_max 19000.0

# _tof_background_time_max 19000.0
# _tof_background_coeff1 24832.850
# _tof_background_coeff2 6139.244
# _tof_background_coeff3 8063.472
# _tof_background_coeff4 3125.050
# _tof_background_coeff5 2566.956
# _tof_background_coeff6 311.077
# _tof_background_coeff7 837.348
# _tof_background_coeff8 -103.742
# _tof_background_coeff9 -11.806

# loop_
# _phase_label
# _phase_scale
# _phase_igsize
# cecual   0.024678(78)   0.0  

# loop_
# _tof_meas_time
# _tof_meas_intensity
# _tof_meas_intensity_sigma
# 3001.589    40409    00462  
# 3003.090    40171    00460  

# _tof_parameters_zero 2.921
# _tof_parameters_dtt1 6167.247
# _tof_parameters_dtt1 -2.280
# _tof_parameters_2theta_bank 145.


# loop_
# _exclude_time_low
# _exclude_time_high
# 2000.1    3000.0   
# 19000.0   20000.0  

# """
# obj = TOF.from_cif(s_cont)
# print(obj)
# for var_name in obj.get_variable_names():
#     print(var_name)
