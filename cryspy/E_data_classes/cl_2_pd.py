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

from cryspy.C_item_loop_classes.cl_1_setup import Setup
from cryspy.C_item_loop_classes.cl_1_diffrn_radiation import \
    DiffrnRadiation
from cryspy.C_item_loop_classes.cl_1_extinction import Extinction
from cryspy.C_item_loop_classes.cl_1_phase import PhaseL
from cryspy.C_item_loop_classes.cl_1_refln import ReflnL
from cryspy.C_item_loop_classes.cl_1_refln_susceptibility import \
    ReflnSusceptibilityL
from cryspy.C_item_loop_classes.cl_1_refine_ls import RefineLs
from cryspy.C_item_loop_classes.cl_1_pd_background import PdBackgroundL
from cryspy.C_item_loop_classes.cl_1_pd_instr_reflex_asymmetry import \
    PdInstrReflexAsymmetry
from cryspy.C_item_loop_classes.cl_1_pd_instr_resolution import \
    PdInstrResolution
from cryspy.C_item_loop_classes.cl_1_pd_meas import PdMeasL
from cryspy.C_item_loop_classes.cl_1_pd_proc import PdProcL
from cryspy.C_item_loop_classes.cl_1_pd_peak import PdPeakL
from cryspy.C_item_loop_classes.cl_1_chi2 import Chi2
from cryspy.C_item_loop_classes.cl_1_range import Range
from cryspy.C_item_loop_classes.cl_1_texture import Texture
from cryspy.C_item_loop_classes.cl_1_exclude import ExcludeL

from cryspy.E_data_classes.cl_1_crystal import Crystal
from cryspy.E_data_classes.cl_1_mag_crystal import MagCrystal


class Pd(DataN):
    """
    Powder diffraction experiment with polarized or unpolarized neutrons (1d).

    Data items in the DIFFRN category record details about
    powder diffraction measurements.

    Methods
    -------
        - calc_profile
        - calc_chi_sq
        - simmulation
        - calc_iint
        - _gauss_pd
        - _lor_pd
        - calc_shape_profile
        - params_to_cif
        - data_to_cif
        - calc_to_cif
        - apply_constraints

    Attributes
    ----------
        - setup, pd_instr_resolution, phase, pd_background, pd_meas (mandatory)
        - diffrn_radiation, chi2, range, extinction, pd_instr_reflex_asymmetry
          texture, exclude, pd_proc, pd_peak, refine_ls, refln_#phase_name
          refln_susceptibility_#phase_name (optional)
    """

    CLASSES_MANDATORY = (Setup, PdInstrResolution, PhaseL, PdBackgroundL,
                         PdMeasL)
    CLASSES_OPTIONAL = (DiffrnRadiation, Chi2, Range, Extinction,
                        PdInstrReflexAsymmetry, Texture, ExcludeL, PdProcL,
                        PdPeakL, RefineLs, ReflnL, ReflnSusceptibilityL)
    # CLASSES_INTERNAL = ()

    CLASSES = CLASSES_MANDATORY + CLASSES_OPTIONAL

    PREFIX = "pd"

    # default values for the parameters
    D_DEFAULT = {}

    def __init__(self, data_name=None, **kwargs) -> NoReturn:
        super(Pd, self).__init__()

        self.__dict__["items"] = []
        self.__dict__["data_name"] = data_name

        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)

    def calc_profile(self, tth, l_crystal, flag_internal: bool = True,
                     flag_polarized: bool = True):
        """Calculate intensity for the given diffraction angles.

        Arguments
        ---------
            - tth: 1D numpy array of 2theta in degrees
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

        proc = PdProcL()
        proc.numpy_ttheta = tth

        try:
            background = self.pd_background
            int_bkgd = background.interpolate_by_points(tth)
        except AttributeError:
            int_bkgd = 0.*tth

        proc.numpy_intensity_bkg_calc = int_bkgd

        setup = self.setup
        wavelength = float(setup.wavelength)

        tth_min = tth.min()
        tth_max = tth.max() + 3.
        if tth_max > 180.:
            tth_max = 180.
        sthovl_min = numpy.sin(0.5*tth_min*numpy.pi/180.)/wavelength
        sthovl_max = numpy.sin(0.5*tth_max*numpy.pi/180.)/wavelength

        try:
            texture = self.texture
            h_ax, k_ax, l_ax = texture.h_ax, texture.k_ax, texture.l_ax
            g_1, g_2 = texture.g_1, texture.g_2
        except AttributeError:
            texture = None

        res_u_1d = numpy.zeros(tth.shape[0], dtype=float)
        res_d_1d = numpy.zeros(tth.shape[0], dtype=float)

        phase = self.phase

        for phase_item in phase.items:
            phase_label = phase_item.label
            phase_scale = phase_item.scale
            try:
                phase_igsize = phase_item.igsize
                if phase_igsize is None: phase_igsize = 0. # temporary solution
            except AttributeError:
                phase_igsize = 0.
            try:
                phase_u = phase_item.u
                if phase_u is None: phase_u = 0. # temporary solution
            except AttributeError:
                phase_u = 0.
            try:
                phase_v = phase_item.v
                if phase_v is None: phase_v = 0. # temporary solution
            except AttributeError:
                phase_v = 0.
            try:
                phase_w= phase_item.w
                if phase_w is None: phase_w = 0. # temporary solution
            except AttributeError:
                phase_w= 0.
            try:
                phase_x = phase_item.x
                if phase_x is None: phase_x = 0. # temporary solution
            except AttributeError:
                phase_x = 0.
            try:
                phase_y = phase_item.y
                if phase_y is None: phase_y = 0. # temporary solution
            except AttributeError:
                phase_y = 0.

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
                peak = PdPeakL(loop_name=phase_label)
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

            tth_hkl_rad = numpy.where(sthovl_hkl*wavelength < 1.,
                                      2.*numpy.arcsin(sthovl_hkl*wavelength),
                                      numpy.pi)
            tth_hkl = tth_hkl_rad*180./numpy.pi

            profile_2d, tth_zs, h_pv = self.calc_shape_profile(
                tth, tth_hkl, phase_igsize=phase_igsize, phase_u=phase_u,
                phase_v=phase_v, phase_w=phase_w, phase_x=phase_x,
                phase_y=phase_y)

            peak.numpy_ttheta = tth_hkl+self.setup.offset_ttheta
            peak.numpy_width_ttheta = h_pv

            np_iint_u_2d = numpy.meshgrid(tth, np_iint_u*mult,
                                          indexing="ij")[1]
            np_iint_d_2d = numpy.meshgrid(tth, np_iint_d*mult,
                                          indexing="ij")[1]

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
                    cos_alpha_ax_2d = numpy.meshgrid(tth, cos_alpha_ax,
                                                     indexing="ij")[1]
                    sin_alpha_ax_2d = numpy.meshgrid(tth, sin_alpha_ax,
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

        proc.numpy_ttheta_corrected = tth_zs
        proc.numpy_intensity_up_net = res_u_1d
        proc.numpy_intensity_down_net = res_d_1d
        proc.numpy_intensity_net = res_u_1d+res_d_1d
        proc.numpy_intensity_diff_total = res_u_1d-res_d_1d
        if flag_polarized:
            proc.numpy_intensity_up_total = res_u_1d+int_bkgd
            proc.numpy_intensity_down_total = res_d_1d+int_bkgd
            proc.numpy_intensity_total = res_u_1d+res_d_1d+int_bkgd+int_bkgd
        else:
            proc.numpy_intensity_up_total = res_u_1d+0.5*int_bkgd
            proc.numpy_intensity_down_total = res_d_1d+0.5*int_bkgd
            proc.numpy_intensity_total = res_u_1d+res_d_1d+int_bkgd

        if flag_internal:
            proc.numpy_to_items()
            l_calc_objs = []
            for crystal in l_crystal:
                try:
                    obj = d_internal_val[f"refln_{crystal.data_name:}"]
                    l_calc_objs.append(obj)
                except KeyError:
                    pass
                try:
                    obj = d_internal_val[
                        f"refln_susceptibility_{crystal.data_name:}"]
                    l_calc_objs.append(obj)
                except KeyError:
                    pass
                try:
                    obj = d_internal_val[f"peak_{crystal.data_name:}"]
                    l_calc_objs.append(obj)
                except KeyError:
                    pass
            l_calc_objs.append(proc)
            self.add_items(l_calc_objs)
        return proc

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
        meas = self.pd_meas
        flag_polarized = meas.is_polarized()

        tth = meas.numpy_ttheta
        if flag_polarized:
            int_u_exp = meas.numpy_intensity_up
            sint_u_exp = meas.numpy_intensity_up_sigma
            int_d_exp = meas.numpy_intensity_down
            sint_d_exp = meas.numpy_intensity_down_sigma
        else:
            int_exp = meas.numpy_intensity
            sint_exp = meas.numpy_intensity_sigma

        cond_in = numpy.ones(tth.shape, dtype=bool)
        try:
            range_ = self.range
            range_min = numpy.array(range_.ttheta_min, dtype=float)
            range_max = numpy.array(range_.ttheta_max, dtype=float)

            cond_in = numpy.logical_and(cond_in, tth >= range_min)
            cond_in = numpy.logical_and(cond_in, tth <= range_max)
        except AttributeError:
            pass

        tth_in = tth[cond_in]
        if flag_polarized:
            int_u_exp_in = int_u_exp[cond_in]
            sint_u_exp_in = sint_u_exp[cond_in]
            int_d_exp_in = int_d_exp[cond_in]
            sint_d_exp_in = sint_d_exp[cond_in]
        else:
            int_exp_in = int_exp[cond_in]
            sint_exp_in = sint_exp[cond_in]

        proc = self.calc_profile(
            tth_in, l_crystal, flag_internal=flag_internal,
            flag_polarized=flag_polarized)

        if flag_polarized:
            proc.numpy_intensity_up = int_u_exp_in
            proc.numpy_intensity_up_sigma = sint_u_exp_in
            proc.numpy_intensity_down = int_d_exp_in
            proc.numpy_intensity_down_sigma = sint_d_exp_in
            proc.numpy_intensity = int_u_exp_in+int_d_exp_in
            proc.numpy_intensity_sigma = numpy.sqrt(
                numpy.square(sint_u_exp_in) + numpy.square(sint_d_exp_in))
        else:
            proc.numpy_intensity = int_exp_in
            proc.numpy_intensity_sigma = sint_exp_in

        int_u_mod = proc.numpy_intensity_up_total
        int_d_mod = proc.numpy_intensity_down_total

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
            l_excl_tth_min = exclude.ttheta_min
            l_excl_tth_max = exclude.ttheta_max
            if flag_polarized:
                for excl_tth_min, excl_tth_max in zip(l_excl_tth_min,
                                                      l_excl_tth_max):
                    cond_1 = numpy.logical_or(tth_in < 1.*excl_tth_min,
                                              tth_in > 1.*excl_tth_max)
                    cond_u = numpy.logical_and(cond_u, cond_1)
                    cond_d = numpy.logical_and(cond_d, cond_1)
                    cond_sum = numpy.logical_and(cond_sum, cond_1)
            else:
                for excl_tth_min, excl_tth_max in zip(l_excl_tth_min,
                                                      l_excl_tth_max):
                    cond_1 = numpy.logical_or(tth_in < 1.*excl_tth_min,
                                              tth_in > 1.*excl_tth_max)
                    cond_sum = numpy.logical_and(cond_sum, cond_1)
        except AttributeError:
            pass

        proc.numpy_excluded = numpy.logical_not(cond_sum)
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
            proc.numpy_to_items()
        return chi_sq_val, n

    def simulation(self, l_crystal, ttheta_start: float = 4.,
                   ttheta_end: float = 120., ttheta_step: float = 0.1,
                   flag_polarized: bool = True) -> NoReturn:
        """Simulate."""
        n_points = int(round((ttheta_end-ttheta_start)/ttheta_step))
        tth_in = numpy.linspace(ttheta_start, ttheta_end, n_points)
        if isinstance(l_crystal, Crystal):
            l_crystal = [l_crystal]
        self.calc_profile(tth_in, l_crystal, l_peak_in=None, l_refln_in=None,
                          l_refln_susceptibility_in=None, l_dd_in=None,
                          flag_internal=True, flag_polarized=flag_polarized)
        return

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

        setup = self.setup
        try:
            field = setup.field
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

    def _gauss_pd(self, tth_2d):
        """One dimensional gauss powder diffraction."""
        ag, bg = self.ag, self.bg
        val_1 = bg*tth_2d**2
        val_2 = numpy.where(val_1 < 5., numpy.exp(-val_1), 0.)
        self.gauss_pd = ag*val_2

    def _lor_pd(self, tth_2d):
        """One dimensional lorentz powder diffraction."""
        al, bl = self.al, self.bl
        self.lor_pd = al*1./(1.+bl*tth_2d**2)

    def calc_shape_profile(
            self, tth, tth_hkl, phase_igsize: float = 0., phase_u: float = 0.,
            phase_v: float = 0., phase_w: float = 0., phase_x: float = 0.,
            phase_y: float = 0.):
        """
        Calculate shape profile.

        Calculate profile in the range ttheta for reflections placed on
        ttheta_hkl with i_g parameter by default equal to zero

        tth, tth_hkl in degrees
        """
        zero_shift = float(self.setup.offset_ttheta)
        tth_zs = tth-zero_shift

        resolution = self.pd_instr_resolution

        h_pv, eta, h_g, h_l, a_g, b_g, a_l, b_l = resolution.calc_resolution(
            tth_hkl, phase_igsize=phase_igsize, phase_u=phase_u,
            phase_v=phase_v, phase_w=phase_w, phase_x=phase_x, phase_y=phase_y)

        tth_2d, tth_hkl_2d = numpy.meshgrid(tth_zs, tth_hkl, indexing="ij")

        self.ag = numpy.meshgrid(tth_zs, a_g, indexing="ij")[1]
        self.bg = numpy.meshgrid(tth_zs, b_g, indexing="ij")[1]
        self.al = numpy.meshgrid(tth_zs, a_l, indexing="ij")[1]
        self.bl = numpy.meshgrid(tth_zs, b_l, indexing="ij")[1]
        eta_2d = numpy.meshgrid(tth_zs, eta, indexing="ij")[1]
        self.eta = eta_2d

        self._gauss_pd(tth_2d-tth_hkl_2d)
        self._lor_pd(tth_2d-tth_hkl_2d)
        g_pd_2d = self.gauss_pd
        l_pd_2d = self.lor_pd

        np_shape_2d = eta_2d * l_pd_2d + (1.-eta_2d) * g_pd_2d

        try:
            asymmetry = self.asymmetry
            np_ass_2d = asymmetry.calc_asymmetry(tth_zs, tth_hkl, h_pv)
        except AttributeError:
            np_ass_2d = numpy.ones(shape=np_shape_2d.shape)

        # Lorentz factor
        tth_rad = tth_zs*numpy.pi/180.
        np_lor_1d = 1./(numpy.sin(tth_rad)*numpy.sin(0.5*tth_rad))

        np_lor_2d = numpy.meshgrid(np_lor_1d, tth_hkl, indexing="ij")[0]

        profile_2d = np_shape_2d*np_ass_2d*np_lor_2d

        return profile_2d, tth_zs, h_pv

    def params_to_cif(self, separator="_", flag: bool = False,
                      flag_minimal: bool = True) -> str:
        """Save parameters to cif format."""
        ls_out = []
        l_cls = (Setup, PdInstrResolution, DiffrnRadiation, Chi2, Range,
                 Extinction, PdInstrReflexAsymmetry, Texture, PhaseL, ExcludeL,
                 PdBackgroundL)
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
        l_cls = (PdMeasL, )
        l_obj = [item for item in self.items if type(item) in l_cls]
        l_itemn = [item for item in l_obj if isinstance(item, ItemN)]
        l_loopn = [item for item in l_obj if isinstance(item, LoopN)]
        ls_out.extend([_.to_cif(separator=separator)+"\n" for _ in l_itemn])
        ls_out.extend([_.to_cif(separator=separator)+"\n" for _ in l_loopn])
        return "\n".join(ls_out)

    def calc_to_cif(self, separator="_", flag=False, flag_minimal=True) -> str:
        """Save calculations to cif format."""
        ls_out = []
        l_cls = (RefineLs, ReflnL, ReflnSusceptibilityL, PdPeakL, PdProcL)
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
        if self.is_attribute("pd_proc"):
            pd_proc = self.pd_proc
            if self.is_attribute("chi2"):
                flag_up = self.chi2.up
                flag_down = self.chi2.down
                flag_sum = self.chi2.sum
                flag_diff = self.chi2.diff
            else:
                flag_up, flag_down, flag_sum = False, False, True
                flag_diff = False
            if flag_sum:
                fig_s, ax_s = pd_proc.plot_sum()
                ax_s.set_title(self.data_name + " - "+ax_s.title.get_text())
                y_min_s, y_max_s = ax_s.get_ylim()
                y_dist_s = y_max_s-y_min_s
                y_step_s = 0.
            
            flag_d = False
            if flag_diff:
                fig_d_ax_d = pd_proc.plot_diff()
                flag_d = fig_d_ax_d is not None
                if flag_d:
                    fig_d, ax_d = fig_d_ax_d
                    ax_d.set_title(self.data_name + " - "+ax_d.title.get_text())
                    y_min_d, y_max_d = ax_d.get_ylim()
                    y_dist_d = y_max_d-y_min_d
                    y_step_d = 0.

            for item in self.items:
                if isinstance(item, PdPeakL):
                    np_tth = item.numpy_ttheta
                    if flag_sum:
                        ax_s.plot(np_tth, 0.*np_tth+y_min_s-y_step_s, "|",
                                  label=item.loop_name)
                        y_step_s += 0.05*y_dist_s
                    if (flag_d & flag_diff):
                        ax_d.plot(np_tth, 0.*np_tth+y_min_d-y_step_d, "|",
                                  label=item.loop_name)
                        y_step_d += 0.05*y_dist_d
            res = []
            if flag_sum:
                ax_s.legend(loc='upper right')
                res.append((fig_s, ax_s))
            if (flag_d & flag_diff):
                ax_d.legend(loc='upper right')
                res.append((fig_d, ax_d))
            return res
        elif self.is_attribute("pd_meas"):
            return self.pd_meas.plots()

# s_cont = """
#  data_pnd
#  _setup_wavelength     0.840
#  _setup_field          1.000
#  _setup_offset_2theta -0.385(12)
#  _chi2_sum  True
#  _chi2_diff False
#  _chi2_up   False
#  _chi2_down False
#  _range_2theta_min     4.000
#  _range_2theta_max    80.000
#  loop_
#  _pd_background_2theta
#  _pd_background_intensity
#   4.5 256.0
#  40.0 158.0
#  80.0  65.0
#  loop_
#  _exclude_2theta_min
#  _exclude_2theta_max
#  0.0 1.0
#  _pd_instr_reflex_asymmetry_p1 0.0
#  _pd_instr_reflex_asymmetry_p2 0.0
#  _pd_instr_reflex_asymmetry_p3 0.0
#  _pd_instr_reflex_asymmetry_p4 0.0
#  _diffrn_radiation_polarization 0.0
#  _diffrn_radiation_efficiency   1.0
#  _pd_instr_resolution_u 16.9776
#  _pd_instr_resolution_v -2.8357
#  _pd_instr_resolution_w  0.5763
#  _pd_instr_resolution_x  0.0
#  _pd_instr_resolution_y  0.0
#  loop_
#  _phase_label
#  _phase_scale
#  _phase_igsize
#  Fe3O4 0.02381 0.0
#  loop_
#  _pd_meas_2theta
#  _pd_meas_intensity_up
#  _pd_meas_intensity_up_sigma
#  _pd_meas_intensity_down
#  _pd_meas_intensity_down_sigma
#  4.0 465.80 128.97 301.88 129.30
#  4.2 323.78 118.22 206.06 120.00
#  4.4 307.14 115.90 230.47 116.53
#  loop_
#  _pd_peak_index_h
#  _pd_peak_index_k
#  _pd_peak_index_l
#  _pd_peak_mult
#  _pd_peak_ttheta
#  _pd_peak_intensity_up
#  _pd_peak_intensity_down
#  _pd_peak_width_2theta
#  1 1 1  8  9.748 128.15576 128.15576 0.677
#  2 0 0  6 11.260   0.00000   0.00000 0.680
#  2 2 0 12 15.950  94.21107  94.21107 0.716
#  loop_
#  _pd_proc_2theta
#  _pd_proc_2theta_corrected
#  _pd_proc_intensity_up_net
#  _pd_proc_intensity_down_net
#  _pd_proc_intensity_up_total
#  _pd_proc_intensity_down_total
#  _pd_proc_intensity_bkg_calc
#  _pd_proc_intensity_up
#  _pd_proc_intensity_up_sigma
#  _pd_proc_intensity_down
#  _pd_proc_intensity_down_sigma
#  4.000 4.385 0.0 0.0 256.0 256.0 256.0 465.8 128.97000 301.88000 129.30000
#  4.200 4.585 0.0 0.0 256.0 256.0 256.0 323.78 118.22000 206.06000 120.00000
#  4.400 4.785 0.0 0.0 256.0 256.0 256.0 307.14 115.90000 230.47000 116.53000
# """
# obj = Pd.from_cif(s_cont)
# print(obj)
# for var_name in obj.get_variable_names():
#     print(var_name)
