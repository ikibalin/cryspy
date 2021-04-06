"""Description of Pd2d class."""
__author__ = 'ikibalin'
__version__ = "2020_08_19"
from warnings import warn
import numpy
from typing import List, NoReturn

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
from cryspy.C_item_loop_classes.cl_1_texture import Texture
from cryspy.C_item_loop_classes.cl_1_exclude import ExcludeL

from cryspy.E_data_classes.cl_1_crystal import Crystal


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
                         Setup, Range, Chi2, Pd2dBackground, Pd2dMeas)
    CLASSES_OPTIONAL = (Extinction, ExcludeL, Pd2dInstrReflexAsymmetry,
                        Texture, RefineLs, ReflnL, ReflnSusceptibilityL,
                        Pd2dPeakL, Pd2dProc)
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

    def calc_profile(self, tth, phi, l_crystal: List[Crystal], l_peak_in=None,
                     l_refln_in=None, l_refln_susceptibility_in=None,
                     l_dd_in=None, flag_internal: bool = True):
        """Calculate intensity for the given diffraction angle."""
        proc = Pd2dProc()  # it is output
        proc.ttheta = tth
        proc.phi = phi

        background = self.pd2d_background
        int_bkgd = background.interpolate_by_points(tth, phi)
        proc.intensity_bkg_calc = int_bkgd

        tth_rad = tth*numpy.pi/180.
        phi_rad = phi*numpy.pi/180.
        cos_theta_1d = numpy.cos(0.5*tth_rad)
        sin_phi_1d = numpy.sin(phi_rad)

        setup = self.setup
        wavelength = setup.wavelength
        diffrn_radiation = self.diffrn_radiation
        if setup.offset_phi is None:
            setup.offset_phi = 0.

        phi_0 = setup.offset_phi
        phi_rad = (phi-phi_0)*numpy.pi/180.
        sin_phi_1d = numpy.sin(phi_rad)

        p_u = float(diffrn_radiation.polarization)
        p_d = (2.*float(diffrn_radiation.efficiency)-1.)*p_u

        tth_min = tth.min()
        tth_max = tth.max()+3.
        if tth_max > 180.:
            tth_max = 180.

        sthovl_min = numpy.sin(0.5*tth_min*numpy.pi/180.)/wavelength
        sthovl_max = numpy.sin(0.5*tth_max*numpy.pi/180.)/wavelength

        res_u_2d = numpy.zeros((tth.shape[0], phi.shape[0]), dtype=float)
        res_d_2d = numpy.zeros((tth.shape[0], phi.shape[0]), dtype=float)

        try:
            texture = self.texture
            h_ax, k_ax = float(texture.h_ax), float(texture.k_ax)
            l_ax = float(texture.l_ax)
            g_1, g_2 = float(texture.g_1), float(texture.g_2)
        except AttributeError:
            texture = None

        phase = self.phase
        _h = len(phase.label)
        if l_peak_in is None:
            l_peak_in = len(phase.label) * [None]
        else:
            if _h != len(l_peak_in):
                l_peak_in = len(phase.label) * [None]
        if l_refln_in is None:
            l_refln_in = len(phase.label) * [None]
        else:
            if _h != len(l_refln_in):
                l_refln_in = len(phase.label) * [None]
        if l_refln_susceptibility_in is None:
            l_refln_susceptibility_in = len(phase.label) * [None]
        else:
            if _h != len(l_refln_susceptibility_in):
                l_refln_susceptibility_in = len(phase.label) * [None]
        if l_dd_in is None:
            l_dd_in = len(phase.label) * [None]
        else:
            if _h != len(l_dd_in):
                l_dd_in = len(phase.label) * [None]

        l_peak, l_dd_out = [], []
        l_refln, l_refln_s = [], []
        for item_phase, peak_in, refln_in, refln_susceptibility_in, dd_in in \
                zip(phase.items, l_peak_in, l_refln_in,
                    l_refln_susceptibility_in, l_dd_in):
            phase_label = item_phase.label
            phase_scale = item_phase.scale
            try:
                phase_igsize = item_phase.igsize
                if phase_igsize is None: phase_igsize = 0. # temporary solution
            except AttributeError:
                phase_igsize = 0.
            try:
                phase_u = item_phase.u
                if phase_u is None: phase_u = 0. # temporary solution
            except AttributeError:
                phase_u = 0.
            try:
                phase_v = item_phase.v
                if phase_v is None: phase_v = 0. # temporary solution
            except AttributeError:
                phase_v = 0.
            try:
                phase_w = item_phase.w
                if phase_w is None: phase_w = 0. # temporary solution
            except AttributeError:
                phase_w = 0.
            try:
                phase_x = item_phase.x
                if phase_x is None: phase_x = 0. # temporary solution
            except AttributeError:
                phase_x = 0.
            try:
                phase_y = item_phase.y
                if phase_y is None: phase_y = 0. # temporary solution
            except AttributeError:
                phase_y = 0.

            dd_out = {}
            for i_crystal, crystal in enumerate(l_crystal):
                if crystal.data_name.lower() == phase_label.lower():
                    ind_cry = i_crystal
                    break
            if ind_cry is None:
                warn(f"Crystal with name '{phase_label:}' is not found.",
                     UserWarning)
                return

            crystal = l_crystal[ind_cry]

            cell = crystal.cell
            space_group = crystal.space_group

            if peak_in is not None:
                index_h = peak_in.numpy_index_h
                index_k = peak_in.numpy_index_k
                index_l = peak_in.numpy_index_l
                mult = peak_in.numpy_index_mult
            else:
                if texture is None:
                    index_h, index_k, index_l, mult = cell.calc_hkl(
                        space_group, sthovl_min, sthovl_max)
                else:
                    index_h, index_k, index_l, mult = cell.calc_hkl_in_range(
                        sthovl_min, sthovl_max)

            peak = Pd2dPeakL(loop_name=phase_label)
            peak.numpy_index_h = index_h
            peak.numpy_index_k = index_k
            peak.numpy_index_l = index_l
            peak.numpy_index_mult = mult
            peak.numpy_to_items()

            cond_1 = len(crystal.get_variable_names()) == 0
            cond_2 = (peak_in is not None) & (refln_in is not None)

            if (cond_1 & cond_2):
                f_nucl_sq = peak_in.numpy_f_nucl_sq
                f_m_p_sin_sq = peak_in.numpy_f_m_p_sin_sq
                f_m_p_cos_sq = peak_in.numpy_f_m_p_cos_sq
                cross_sin = peak_in.numpy_cross_sin
                refln = refln_in
                refln_s = refln_susceptibility_in
            else:
                f_nucl_sq, f_m_p_sin_sq, f_m_p_cos_sq, cross_sin, refln, \
                    refln_s = self.calc_for_iint(index_h, index_k, index_l,
                                                 crystal)
                refln.loop_name = phase_label
                refln_s.loop_name = phase_label
            l_refln.append(refln)
            l_refln_s.append(refln_s)

            peak.numpy_f_nucl_sq = f_nucl_sq
            peak.numpy_f_m_p_sin_sq = f_m_p_sin_sq
            peak.numpy_f_m_p_cos_sq = f_m_p_cos_sq
            peak.numpy_cross_sin = cross_sin

            cond_1 = dd_in is not None
            cond_2 = ((len(crystal.get_variable_names()) == 0) &
                      (not(diffrn_radiation.polarization_refinement)) &
                      (not(diffrn_radiation.efficiency_refinement)))
            if cond_1 & cond_2:
                iint_u_3d, iint_d_3d = dd_in["iint_u_3d"], dd_in["iint_d_3d"]
                cos_theta_3d = dd_in["cos_theta_3d"]
                sin_phi_3d = dd_in["sin_phi_3d"]
            else:
                cos_theta_3d, sin_phi_3d, mult_f_n_3d = numpy.meshgrid(
                    cos_theta_1d, sin_phi_1d, mult*f_nucl_sq, indexing="ij")

                mult_f_m_c_3d = numpy.meshgrid(
                    tth_rad, phi_rad, mult*f_m_p_cos_sq, indexing="ij")[2]

                hh_u_s_3d = numpy.meshgrid(
                    tth_rad, phi_rad, mult*(f_m_p_sin_sq+p_u*cross_sin),
                    indexing="ij")[2]
                hh_d_s_3d = numpy.meshgrid(
                    tth_rad, phi_rad, mult*(f_m_p_sin_sq-p_d*cross_sin),
                    indexing="ij")[2]

                c_a_sq_3d = (cos_theta_3d * sin_phi_3d)**2
                s_a_sq_3d = 1.-c_a_sq_3d
                iint_u_3d = (mult_f_n_3d + hh_u_s_3d*s_a_sq_3d +
                             mult_f_m_c_3d*c_a_sq_3d)
                iint_d_3d = (mult_f_n_3d + hh_d_s_3d*s_a_sq_3d +
                             mult_f_m_c_3d*c_a_sq_3d)
            dd_out["cos_theta_3d"] = cos_theta_3d
            dd_out["sin_phi_3d"] = sin_phi_3d
            dd_out["iint_u_3d"] = iint_u_3d
            dd_out["iint_d_3d"] = iint_d_3d

            sthovl_hkl = cell.calc_sthovl(index_h, index_k, index_l)
            tth_hkl_rad = numpy.where(sthovl_hkl*wavelength < 1.,
                                      2.*numpy.arcsin(sthovl_hkl*wavelength),
                                      numpy.pi)
            tth_hkl = tth_hkl_rad*180./numpy.pi

            cond_1 = dd_in is not None
            flag_phase_uvwxy_ref = (
                item_phase.igsize_refinement | item_phase.u_refinement |
                item_phase.v_refinement | item_phase.w_refinement |
                item_phase.x_refinement | item_phase.y_refinement)
            cond_2 = ((not(flag_phase_uvwxy_ref)) &
                      (not(self.pd2d_instr_resolution.is_variables())) &
                      (not(setup.is_variables())))
            if cond_1 & cond_2:
                profile_3d, tth_zs = dd_in["profile_3d"], dd_in["tth_zs"]
                h_pv = dd_in["h_pv"]
            else:
                profile_3d, tth_zs, h_pv = self.calc_shape_profile(
                    tth, phi, tth_hkl, phase_igsize=phase_igsize,
                    phase_u=phase_u, phase_v=phase_v, phase_w=phase_w,
                    phase_x=phase_x, phase_y=phase_y)
            dd_out["profile_3d"] = profile_3d
            dd_out["tth_zs"] = tth_zs
            dd_out["h_pv"] = h_pv

            peak.numpy_ttheta = tth_hkl + setup.offset_ttheta
            peak.width_ttheta = h_pv
            peak.numpy_to_items()

            # texture
            if texture is not None:
                cond_1 = dd_in is not None
                cond_2 = ((not(setup.offset_phi_refinement)))
                if cond_1 & cond_2:
                    cos_alpha_ang_3d = dd_in["cos_alpha_ang_3d"]
                    sin_alpha_ang_3d = dd_in["sin_alpha_ang_3d"]
                else:
                    cos_alpha_ang_3d = cos_theta_3d * sin_phi_3d
                    sin_alpha_ang_3d = numpy.sqrt(1.-cos_alpha_ang_3d**2)
                dd_out["cos_alpha_ang_3d"] = cos_alpha_ang_3d
                dd_out["sin_alpha_ang_3d"] = sin_alpha_ang_3d
                cond_2 = (cond_2 & (not(texture.is_variables())) &
                          (not(cell.is_variables())))
                if cond_1 & cond_2:
                    texture_3d = dd_in["texture_3d"]
                else:
                    cos_alpha_ax = calc_cos_ang(cell, h_ax, k_ax, l_ax,
                                                index_h, index_k, index_l)

                    c_help = 1.-cos_alpha_ax**2
                    c_help[c_help < 0.] = 0.
                    sin_alpha_ax = numpy.sqrt(c_help)
                    cos_alpha_ax_3d = numpy.meshgrid(tth, phi, cos_alpha_ax,
                                                     indexing="ij")[2]
                    sin_alpha_ax_3d = numpy.meshgrid(tth, phi, sin_alpha_ax,
                                                     indexing="ij")[2]
                    cos_alpha_3d = cos_alpha_ax_3d*cos_alpha_ang_3d + \
                        sin_alpha_ax_3d*sin_alpha_ang_3d
                    texture_3d = g_2 + (1.-g_2) * \
                        (1./g_1 + (g_1**2-1./g_1)*cos_alpha_3d**2)**(-1.5)
                dd_out["texture_3d"] = texture_3d

                profile_3d = profile_3d*texture_3d

            res_u_3d = profile_3d*iint_u_3d
            res_d_3d = profile_3d*iint_d_3d

            # 0.5 to have the same meaning for scale factor as in FullProf
            res_u_2d += 0.5*phase_scale*res_u_3d.sum(axis=2)
            res_d_2d += 0.5*phase_scale*res_d_3d.sum(axis=2)
            l_peak.append(peak)
            l_dd_out.append(dd_out)

        proc.ttheta_corrected = tth_zs
        proc.intensity_up_net = res_u_2d
        proc.intensity_down_net = res_d_2d
        proc.intensity_up_total = res_u_2d+int_bkgd
        proc.intensity_down_total = res_d_2d+int_bkgd
        if flag_internal:
            if background.is_variables():
                background.form_ttheta_phi_intensity()
            proc.form_ttheta_phi_intensity_up_net()
            proc.form_ttheta_phi_intensity_down_net()
            proc.form_ttheta_phi_intensity_up_total()
            proc.form_ttheta_phi_intensity_down_total()
        l_calc_objs = l_refln + l_refln_s + l_peak
        l_calc_objs.append(proc)
        self.add_items(l_calc_objs)

        return proc, l_peak, l_refln, l_dd_out

    def calc_chi_sq(self, l_crystal, flag_internal: bool = True):
        """
        Calculate chi square.

        Arguments
        ---------
            - l_crystal: a list of Crystal objects of cryspy library
            - flag_internal: a flag to calculate internal objects
              (default is True)

        Output arguments
        ----------------
            - chi_sq_val: chi square of flip ratio
              (Sum_i ((y_e_i - y_m_i) / sigma_i)**2)
            - n: number of measured reflections
        """
        meas = self.pd2d_meas

        tth = meas.ttheta
        phi = meas.phi
        int_u_exp = meas.intensity_up
        sint_u_exp = meas.intensity_up_sigma
        int_d_exp = meas.intensity_down
        sint_d_exp = meas.intensity_down_sigma

        l_peak_in, l_refln_in = [], []
        l_refln_susceptibility_in = []
        l_dd_in = []
        flag_1 = not(flag_internal)

        try:
            dd = self.dd
            flag_2 = True
        except AttributeError:
            flag_2 = False
        if (flag_1 & flag_2):
            for phase_item in self.phase.items:
                crystal = None
                for cryst in l_crystal:
                    if cryst.data_name.lower() == phase_item.label.lower():
                        crystal = cryst
                        break
                attr_peak = f"pd2d_peak_{crystal.data_name:}"
                attr_refln = f"refln_{crystal.data_name:}"
                attr_refln_s = f"refln_susceptibility_{crystal.data_name:}"
                l_peak_in.append(getattr(self, attr_peak))
                l_refln_in.append(getattr(self, attr_refln))
                l_refln_susceptibility_in.append(getattr(self, attr_refln_s))

        cond_tth_in = numpy.ones(tth.size, dtype=bool)
        cond_phi_in = numpy.ones(phi.size, dtype=bool)
        try:
            range_ = self.range
            cond_tth_in = numpy.logical_and(cond_tth_in, tth >=
                                            range_.ttheta_min)
            cond_tth_in = numpy.logical_and(cond_tth_in, tth <=
                                            range_.ttheta_max)

            cond_phi_in = numpy.logical_and(cond_phi_in, phi >= range_.phi_min)
            cond_phi_in = numpy.logical_and(cond_phi_in, phi <= range_.phi_max)

        except AttributeError:
            pass
        # cond_1_in, cond_2_in = numpy.meshgrid(cond_tth_in, cond_phi_in,
        #                                       indexing="ij")
        # cond_in = numpy.logical_and(cond_1_in, cond_2_in)
        tth_in = tth[cond_tth_in]
        phi_in = phi[cond_phi_in]
        int_u_exp_in = int_u_exp[cond_tth_in, :][:, cond_phi_in]
        sint_u_exp_in = sint_u_exp[cond_tth_in, :][:, cond_phi_in]
        int_d_exp_in = int_d_exp[cond_tth_in, :][:, cond_phi_in]
        sint_d_exp_in = sint_d_exp[cond_tth_in, :][:, cond_phi_in]

        proc, l_peak, l_refln, l_dd_out = self.calc_profile(
            tth_in, phi_in, l_crystal, l_peak_in=l_peak_in,
            l_refln_in=l_refln_in,
            l_refln_susceptibility_in=l_refln_susceptibility_in,
            l_dd_in=l_dd_in, flag_internal=flag_internal)
        proc.intensity_up = int_u_exp_in
        proc.intensity_up_sigma = sint_u_exp_in
        proc.intensity_down = int_d_exp_in
        proc.intensity_down_sigma = sint_d_exp_in

        # self.proc = proc
        # self.peaks = l_peak
        # self.reflns = l_refln
        self.dd = l_dd_out

        int_u_mod = proc.intensity_up_total
        int_d_mod = proc.intensity_down_total

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

        # exclude region
        try:
            exclude = self.exclude
            l_excl_tth_min = exclude.numpy_ttheta_min
            l_excl_tth_max = exclude.numpy_ttheta_max
            l_excl_phi_min = exclude.numpy_phi_min
            l_excl_phi_max = exclude.numpy_phi_max
            for excl_tth_min, excl_tth_max, excl_phi_min, excl_phi_max in \
                zip(l_excl_tth_min, l_excl_tth_max, l_excl_phi_min,
                    l_excl_phi_max):
                cond_1 = numpy.logical_or(tth_in < 1.*excl_tth_min,
                                          tth_in > 1.*excl_tth_max)
                cond_2 = numpy.logical_or(phi_in < 1.*excl_phi_min,
                                          phi_in > 1.*excl_phi_max)
                cond_11, cond_22 = numpy.meshgrid(cond_1, cond_2,
                                                  indexing="ij")
                cond_12 = numpy.logical_or(cond_11, cond_22)
                cond_u = numpy.logical_and(cond_u, cond_12)
                cond_d = numpy.logical_and(cond_d, cond_12)
                cond_sum = numpy.logical_and(cond_sum, cond_12)
        except AttributeError:
            pass

        chi_sq_u_val = (chi_sq_u[cond_u]).sum()
        n_u = cond_u.sum()

        chi_sq_d_val = (chi_sq_d[cond_d]).sum()
        n_d = cond_d.sum()

        chi_sq_sum_val = (chi_sq_sum[cond_sum]).sum()
        n_sum = cond_sum.sum()

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
        # print(f"chi_sq_val/n: {chi_sq_val/n:.2f}    \
        #     chi_sq_val: {chi_sq_val: .2f}")
        # d_exp_out = {"chi_sq_val": chi_sq_val, "n": n}
        # d_exp_out.update(d_exp_prof_out)

        if flag_internal:
            refine_ls = RefineLs(number_reflns=n,
                                 goodness_of_fit_all=chi_sq_val/float(n),
                                 weighting_scheme="sigma")
            self.refine_ls = refine_ls
            proc.form_ttheta_phi_intensity_bkg_calc()
            proc.form_ttheta_phi_intensity_up()
            proc.form_ttheta_phi_intensity_up_sigma()
            proc.form_ttheta_phi_intensity_down()
            proc.form_ttheta_phi_intensity_down_sigma()
        return chi_sq_val, n

    def calc_for_iint(self, index_h, index_k, index_l, crystal,
                      flag_internal: bool = True):
        """Calculate the integral intensity for h, k, l reflections."""
        setup = self.setup
        field = float(setup.field)

        refln = crystal.calc_refln(index_h, index_k, index_l)

        refln_s = crystal.calc_refln_susceptibility(index_h, index_k, index_l)

        f_nucl = refln.numpy_f_calc

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
        # k_loc = cell.calc_k_loc(index_h, index_k, index_l)
        t_ij = cell.calc_m_t(index_h, index_k, index_l)
        # FIXME: I would like to recheck the expression for T
        #        and expression SIGMA = T^T CHI T

        t_tr_ij = (t_ij[0], t_ij[3], t_ij[6],
                   t_ij[1], t_ij[4], t_ij[7],
                   t_ij[2], t_ij[5], t_ij[8])
        th_11, th_12, th_13, th_21, th_22, th_23, th_31, th_32, th_33 = \
            calc_mRmCmRT(t_tr_ij, _ij)

        # f_m_p_sin_sq = (field**2)*abs(0.5*(th_11*th_11.conjugate()+\
        #                th_22*th_22.conjugate())+th_12*th_12.conjugate())
        # f_m_p_cos_sq = (field**2)*abs(th_13*th_13.conjugate()+\
        #                 th_23*th_23.conjugate())
        # f_m_p_field = 0.5*field*(th_11+th_22)

        f_nucl_sq = abs(f_nucl*f_nucl.conjugate())
        f_m_p_sin_sq = abs(0.5 * (th_11 * th_11.conjugate()+th_22 *
                                  th_22.conjugate())+th_12 * th_12.conjugate())
        f_m_p_cos_sq = abs(th_13 * th_13.conjugate() +
                           th_23 * th_23.conjugate())
        f_m_p_field = 0.5 * (th_11+th_22)
        cross_sin = 2. * (f_nucl.real * f_m_p_field.real +
                          f_nucl.imag * f_m_p_field.imag)

        return f_nucl_sq, f_m_p_sin_sq, f_m_p_cos_sq, cross_sin, refln, refln_s

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
            self, tth, phi, tth_hkl, phase_igsize: float = 0.,
            phase_u: float = 0., phase_v: float = 0., phase_w: float = 0.,
            phase_x: float = 0., phase_y: float = 0.):
        """
        Calculate profile in the range ttheta.

        For reflections placed on
        ttheta_hkl with i_g parameter by default equal to zero

        tth, phi, tth_hkl in degrees
        """
        setup = self.setup
        zero_shift = float(setup.offset_ttheta)
        tth_zs = tth-zero_shift

        resolution = self.pd2d_instr_resolution

        h_pv, eta, h_g, h_l, a_g, b_g, a_l, b_l = resolution.calc_resolution(
            tth_hkl, phase_igsize=phase_igsize, phase_u=phase_u,
            phase_v=phase_v, phase_w=phase_w, phase_x=phase_x, phase_y=phase_y)

        tth_3d, phi_3d, tth_hkl_3d = numpy.meshgrid(tth_zs, phi, tth_hkl,
                                                    indexing="ij")

        self.ag = numpy.meshgrid(tth_zs, phi, a_g, indexing="ij")[2]
        self.bg = numpy.meshgrid(tth_zs, phi, b_g, indexing="ij")[2]
        self.al = numpy.meshgrid(tth_zs, phi, a_l, indexing="ij")[2]
        self.bl = numpy.meshgrid(tth_zs, phi, b_l, indexing="ij")[2]
        eta_3d = numpy.meshgrid(tth_zs, phi, eta, indexing="ij")[2]
        self.eta = eta_3d

        self._gauss_pd(tth_3d-tth_hkl_3d)
        self._lor_pd(tth_3d-tth_hkl_3d)
        g_pd2d_3d = self.gauss_pd
        l_pd2d_3d = self.lor_pd

        np_shape_3d = eta_3d * l_pd2d_3d + (1.-eta_3d) * g_pd2d_3d

        asymmetry = self.pd2d_instr_reflex_asymmetry
        np_ass_2d = asymmetry.calc_asymmetry(tth_zs, tth_hkl, h_pv)
        np_ass_3d = np_ass_2d[:, numpy.newaxis, :] * numpy.ones(
            phi.size, dtype=float)[numpy.newaxis, :, numpy.newaxis]

        # Lorentz factor
        tth_rad = tth_zs*numpy.pi/180.
        np_lor_1d = 1./(numpy.sin(tth_rad)*numpy.sin(0.5*tth_rad))
        np_lor_3d = numpy.meshgrid(np_lor_1d, phi, tth_hkl, indexing="ij")[0]

        profile_3d = np_shape_3d*np_ass_3d*np_lor_3d

        return profile_3d, tth_zs, h_pv

    def params_to_cif(self, separator="_", flag: bool = False,
                      flag_minimal: bool = True) -> str:
        """Save parameters to cif format."""
        ls_out = []
        l_cls = (Pd2dBackground, Pd2dInstrResolution, PhaseL, DiffrnRadiation,
                 Setup, Range, Chi2, Extinction, ExcludeL,
                 Pd2dInstrReflexAsymmetry, Texture)
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
        l_cls = (Pd2dMeas, )
        l_obj = [item for item in self.items if type(item) in l_cls]
        l_itemn = [item for item in l_obj if isinstance(item, ItemN)]
        l_loopn = [item for item in l_obj if isinstance(item, LoopN)]
        ls_out.extend([_.to_cif(separator=separator)+"\n" for _ in l_itemn])
        ls_out.extend([_.to_cif(separator=separator)+"\n" for _ in l_loopn])
        return "\n".join(ls_out)

    def calc_to_cif(self, separator="_", flag=False, flag_minimal=True) -> str:
        """Save calculations to cif format."""
        ls_out = []
        l_cls = (RefineLs, ReflnL, ReflnSusceptibilityL, Pd2dPeakL, Pd2dProc)
        l_obj = [item for item in self.items if type(item) in l_cls]
        l_itemn = [item for item in l_obj if isinstance(item, ItemN)]
        l_loopn = [item for item in l_obj if isinstance(item, LoopN)]
        ls_out.extend([_.to_cif(separator=separator)+"\n" for _ in l_itemn])
        ls_out.extend([_.to_cif(separator=separator)+"\n" for _ in l_loopn])
        return "\n".join(ls_out)

    def apply_constraints(self):
        """Apply constraints."""
        pass

    def plots(self):
        if self.is_attribute("pd2d_proc"):
            pd2d_proc = self.pd2d_proc
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
        
# s_cont = """
#  data_powder2d
#  _setup_wavelength     0.840
#  _setup_field          1.000
#  _setup_offset_2theta -0.385
#  _chi2_sum  True
#  _chi2_diff False
#  _chi2_up   False
#  _chi2_down False
#  _range_2theta_min     4.000
#  _range_2theta_max    80.000
#  _range_phi_min       -2.000
#  _range_phi_max       40.000
#  _pd2d_background_2theta_phi_intensity
#  ;
#  2 5 80
#  5 2. 2.(1)
#  40 2. 2.
#  ;
#  loop_
#  _exclude_2theta_min
#  _exclude_2theta_max
#  _exclude_phi_min
#  _exclude_phi_max
#  0.0 1.0 0. 1.0
#  _pd2d_instr_reflex_asymmetry_p1 0.0
#  _pd2d_instr_reflex_asymmetry_p2 0.0
#  _pd2d_instr_reflex_asymmetry_p3 0.0
#  _pd2d_instr_reflex_asymmetry_p4 0.0
#  _diffrn_radiation_polarization 0.0
#  _diffrn_radiation_efficiency   1.0

#  _pd2d_instr_resolution_u 16.9776
#  _pd2d_instr_resolution_v -2.8357(1)
#  _pd2d_instr_resolution_w  0.5763
#  _pd2d_instr_resolution_x  0.0
#  _pd2d_instr_resolution_y  0.0
#  loop_
#  _phase_label
#  _phase_scale
#  _phase_igsize
#  Fe3O4 0.02381() 0.0
#  _pd2d_meas_2theta_phi_intensity_up
#  ;
#  2 5 80
#  5 2. 2.
#  40 2. 2.
#  ;
#  _pd2d_meas_2theta_phi_intensity_up_sigma
#  ;
#  2 5 80
#  5 2. 2.
#  40 2. 2.
#  ;
#  _pd2d_meas_2theta_phi_intensity_down
#  ;
#  2 5 80
#  5 2. 2.
#  40 2. 2.
#  ;
#  _pd2d_meas_2theta_phi_intensity_down_sigma
#  ;
#  2 5 80
#  5 2. 2.
#  40 2. 2.
#  ;
# """
# obj = Pd2d.from_cif(s_cont)
# for item in obj.items:
#     print(item, end="\n\n")
# print(50*"*", "\n", obj.is_defined())
# for var_name in obj.get_variable_names():
#     print(var_name)
