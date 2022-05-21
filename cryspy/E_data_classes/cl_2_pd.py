"""Description of Pd class."""
__author__ = 'ikibalin'
__version__ = "2020_08_19"

import numpy
from typing import NoReturn

from cryspy.A_functions_base.function_2_crystallography_base import \
    calc_cos_ang

from cryspy.A_functions_base.matrix_operations import calc_m1_m2_m1t
from cryspy.A_functions_base.unit_cell import calc_matrix_t

from cryspy.A_functions_base.integrated_intensity_powder_diffraction import calc_powder_iint_1d_para, calc_powder_iint_1d_ordered

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN
from cryspy.B_parent_classes.cl_3_data import DataN
from cryspy.B_parent_classes.preocedures import take_items_by_class

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
from cryspy.C_item_loop_classes.cl_1_texture import Texture, TextureL
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
                        PdInstrReflexAsymmetry, TextureL, ExcludeL, PdProcL,
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
                
                mult = peak.numpy_index_mult
                cond_2 = True
                if self.is_attribute("diffrn_radiation"):
                    diffrn_radiation = self.diffrn_radiation
                    cond_21 = not(
                        diffrn_radiation.polarization_refinement | 
                        diffrn_radiation.efficiency_refinement)
                    cond_2 = cond_2 & cond_21
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
                peak.numpy_index_mult = numpy.array(mult, dtype=int)
                d_internal_val[f"peak_{crystal.data_name:}"] = peak
                cond_2 = False
            
            index_hkl = numpy.stack([index_h, index_k, index_l], axis=0)
            cond_1 = not(crystal.is_variables())
            if cond_1 & cond_2:
                np_iint_u = peak.numpy_intensity_plus
                np_iint_d = peak.numpy_intensity_minus
            else:
                np_iint_u, np_iint_d = self.calc_iint(index_hkl, crystal, flag_internal=flag_internal)
                peak.numpy_intensity_plus = np_iint_u
                peak.numpy_intensity_minus = np_iint_d

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
                cond_1 = True
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
                    cos_alpha_2d = cos_alpha_ax_2d*0.+sin_alpha_ax_2d*1.
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
        proc.numpy_intensity_plus_net = res_u_1d
        proc.numpy_intensity_minus_net = res_d_1d
        proc.numpy_intensity_net = res_u_1d+res_d_1d
        proc.numpy_intensity_diff_total = res_u_1d-res_d_1d
        if flag_polarized:
            proc.numpy_intensity_plus_total = res_u_1d+int_bkgd
            proc.numpy_intensity_minus_total = res_d_1d+int_bkgd
            proc.numpy_intensity_total = res_u_1d+res_d_1d+int_bkgd+int_bkgd
        else:
            proc.numpy_intensity_plus_total = res_u_1d+0.5*int_bkgd
            proc.numpy_intensity_minus_total = res_d_1d+0.5*int_bkgd
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
            int_u_exp = meas.numpy_intensity_plus
            sint_u_exp = meas.numpy_intensity_plus_sigma
            int_d_exp = meas.numpy_intensity_minus
            sint_d_exp = meas.numpy_intensity_minus_sigma
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
            proc.numpy_intensity_plus = int_u_exp_in
            proc.numpy_intensity_plus_sigma = sint_u_exp_in
            proc.numpy_intensity_minus = int_d_exp_in
            proc.numpy_intensity_minus_sigma = sint_d_exp_in
            proc.numpy_intensity = int_u_exp_in+int_d_exp_in
            proc.numpy_intensity_sigma = numpy.sqrt(
                numpy.square(sint_u_exp_in) + numpy.square(sint_d_exp_in))
        else:
            proc.numpy_intensity = int_exp_in
            proc.numpy_intensity_sigma = sint_exp_in

        int_u_mod = proc.numpy_intensity_plus_total
        int_d_mod = proc.numpy_intensity_minus_total

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

    def calc_iint(self, index_hkl, crystal: Crystal,
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
            magnetic_field = setup.field
            if magnetic_field is None:
                magnetic_field = 0.
        except AttributeError:
            magnetic_field = 0.

        try:
            diffrn_radiation = self.diffrn_radiation
            polarization = diffrn_radiation.polarization
            flipper = diffrn_radiation.efficiency
        except AttributeError:
            polarization = 0.0
            flipper = 0.0

        try:
            if (not(flag_internal) | crystal.is_variables()):
                raise KeyError
            refln = d_internal_val[f"refln_{crystal.data_name:}"]
        except KeyError:
            refln = crystal.calc_refln(index_hkl, flag_internal=flag_internal)
            d_internal_val[f"refln_{crystal.data_name:}"] = refln
        f_nucl = refln.numpy_f_calc

        if isinstance(crystal, Crystal):
            try:
                if (not(flag_internal) | crystal.is_variables()):
                    raise KeyError
                refln_s = d_internal_val[
                    f"refln_susceptibility_{crystal.data_name:}"]
            except KeyError:
                refln_s = crystal.calc_refln_susceptibility(
                    index_hkl, flag_internal=flag_internal)
                d_internal_val[f"refln_susceptibility_{crystal.data_name:}"] \
                    = refln_s

            sft_ccs_11 = refln_s.numpy_chi_11_calc
            sft_ccs_12 = refln_s.numpy_chi_12_calc
            sft_ccs_13 = refln_s.numpy_chi_13_calc
            sft_ccs_21 = refln_s.numpy_chi_21_calc
            sft_ccs_22 = refln_s.numpy_chi_22_calc
            sft_ccs_23 = refln_s.numpy_chi_23_calc
            sft_ccs_31 = refln_s.numpy_chi_31_calc
            sft_ccs_32 = refln_s.numpy_chi_32_calc
            sft_ccs_33 = refln_s.numpy_chi_33_calc
            sft_ccs = numpy.stack([
                sft_ccs_11, sft_ccs_12, sft_ccs_13,
                sft_ccs_21, sft_ccs_22, sft_ccs_23,
                sft_ccs_31, sft_ccs_32, sft_ccs_33], axis=0)

            cell = crystal.cell
            unit_cell_parameters = cell.get_unit_cell_parameters()
            
            matrix_t = calc_matrix_t(index_hkl, unit_cell_parameters, flag_unit_cell_parameters=False)[0]
            #SIGMA = T CHI T^-1 = T chi T^T (because T is rotation matrix, therefore T^-1 = T^T)
            tensor_sigma, dder_tensor_sigma = calc_m1_m2_m1t(matrix_t, sft_ccs, flag_m1=False, flag_m2=False)
            iint_plus, iint_minus, dder_iint_plus, dder_iint_minus = calc_powder_iint_1d_para(f_nucl, tensor_sigma, polarization, flipper, magnetic_field,
                flag_f_nucl=False, flag_tensor_sigma=False,
                flag_polarization=False, flag_flipper=False)

        elif isinstance(crystal, MagCrystal):
            try:
                if (not(flag_internal) | crystal.is_variables()):
                    raise KeyError
                f_mag_perp = d_internal_val[f"f_mag_perp_{crystal.data_name:}"]

            except KeyError:
                f_mag_perp = crystal.calc_f_mag_perp(index_hkl)
                d_internal_val[f"f_mag_perp_{crystal.data_name:}"] = f_mag_perp

            iint_plus, dder_iint_plus, = calc_powder_iint_1d_ordered(
                f_nucl, f_mag_perp)
            iint_minus = numpy.copy(iint_plus) # For 1D case I_plus = I_minus
        return iint_plus, iint_minus

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
            
            fig_s, ax_s = pd_proc.plot_sum_diff()
            ax_s = fig_s.axes[0]
            ax_hkl = fig_s.axes[1]
            ax_s.set_title(self.data_name + " - "+ax_s.title.get_text())
            y_min_s, y_max_s = ax_hkl.get_ylim()
            y_dist_s = y_max_s-y_min_s
            y_step_s = 0.
            
            flag_d = False
            # if flag_diff:
            #     fig_d_ax_d = pd_proc.plot_diff()
            #     flag_d = fig_d_ax_d is not None
            #     if flag_d:
            #         fig_d, ax_d = fig_d_ax_d
            #         ax_d.set_title(self.data_name + " - "+ax_d.title.get_text())
            #         y_min_d, y_max_d = ax_d.get_ylim()
            #         y_dist_d = y_max_d-y_min_d
            #         y_step_d = 0.

            for item in self.items:
                if isinstance(item, PdPeakL):
                    np_tth = item.numpy_ttheta
                    ax_hkl.plot(np_tth, 0.*np_tth+y_min_s-y_step_s, "|", label=item.loop_name)
                    y_step_s += 0.05*y_dist_s
            res = []
            ax_s.legend(loc='upper right')
            res.append((fig_s, ax_s))
            return res
        elif self.is_attribute("pd_meas"):
            return self.pd_meas.plots()


    def get_dictionary(self):
        """Form dictionary. See documentation module CrysPy using Jupyter notebook.
        """
        self.form_object()
        ddict = {}
        setup, pd_meas, resolution = None, None, None
        pd_background, range_, exclude = None, None, None
        asymmetry, diffrn_radiation = None, None
        phase, texture, chi2 = None, None, None

        l_obj = take_items_by_class(self, (Setup, ))
        if len(l_obj) > 0:
            setup = l_obj[0]

        l_obj = take_items_by_class(self, (PdMeasL, ))
        if len(l_obj) > 0:
            pd_meas = l_obj[0]

        l_obj = take_items_by_class(self, (PdInstrResolution, ))
        if len(l_obj) > 0:
            resolution = l_obj[0]

        l_obj = take_items_by_class(self, (PdBackgroundL, ))
        if len(l_obj) > 0:
            pd_background = l_obj[0]

        l_obj = take_items_by_class(self, (Range, ))
        if len(l_obj) > 0:
            range_ = l_obj[0]

        l_obj = take_items_by_class(self, (ExcludeL, ))
        if len(l_obj) > 0:
            exclude = l_obj[0]

        l_obj = take_items_by_class(self, (PdInstrReflexAsymmetry, ))
        if len(l_obj) > 0:
            asymmetry = l_obj[0]

        l_obj = take_items_by_class(self, (DiffrnRadiation, ))
        if len(l_obj) > 0:
            diffrn_radiation = l_obj[0]

        l_obj = take_items_by_class(self, (PhaseL, ))
        if len(l_obj) > 0:
            phase = l_obj[0]

        l_obj = take_items_by_class(self, (TextureL, ))
        if len(l_obj) > 0:
            texture = l_obj[0]

        l_obj = take_items_by_class(self, (Chi2, ))
        if len(l_obj) > 0:
            chi2 = l_obj[0]

        ddict["name"] = self.data_name
        ddict["type_name"] = self.get_name()
        if setup is not None:
            if setup.is_attribute("field"):
                ddict["magnetic_field"] = numpy.array([setup.field], dtype=float)
            if setup.is_attribute("temperature"):
                ddict["temperature"] = numpy.array([setup.temperature], dtype=float)
            ddict["wavelength"] = numpy.array([setup.wavelength], dtype=float)
            ddict["flags_wavelength"] = numpy.array([setup.wavelength_refinement], dtype=bool)
            ddict["offset_ttheta"] = numpy.array([setup.offset_ttheta * numpy.pi/180.], dtype=float)
            ddict["flags_offset_ttheta"] = numpy.array([setup.offset_ttheta_refinement], dtype=bool)
            ddict["radiation"] = numpy.array([setup.radiation], dtype=str)
            ddict["k"] = numpy.array([setup.k], dtype=float)
            ddict["cthm"] = numpy.array([setup.cthm], dtype=float)

        if chi2 is not None:
            ddict["flag_chi_sq_sum"] = chi2.sum
            ddict["flag_chi_sq_difference"] = chi2.diff

        if pd_meas is not None:
            ttheta_deg = numpy.array(pd_meas.ttheta, dtype=float)

            ttheta_min_deg = range_.ttheta_min
            ttheta_max_deg = range_.ttheta_max
            flag_in = numpy.logical_and(
                ttheta_deg >= ttheta_min_deg,
                ttheta_deg <= ttheta_max_deg)
            
            ddict["ttheta"] = ttheta_deg[flag_in] * numpy.pi/180.
            if pd_meas.is_attribute("intensity_plus"):
                int_plus = numpy.array(pd_meas.intensity_plus, dtype=float)[flag_in]
                s_int_plus = numpy.array(pd_meas.intensity_plus_sigma, dtype=float)[flag_in]
                int_minus = numpy.array(pd_meas.intensity_minus, dtype=float)[flag_in]
                s_int_minus = numpy.array(pd_meas.intensity_minus_sigma, dtype=float)[flag_in]
                ddict["signal_exp_plus"] = numpy.stack([int_plus, s_int_plus], axis=0)
                ddict["signal_exp_minus"] = numpy.stack([int_minus, s_int_minus], axis=0)
            else:
                int_sum = numpy.array(pd_meas.intensity, dtype=float)[flag_in]
                s_int_sum = numpy.array(pd_meas.intensity_sigma, dtype=float)[flag_in]
                ddict["signal_exp"] = numpy.stack([int_sum, s_int_sum], axis=0)

            ttheta_in_range = ttheta_deg[flag_in]
            flag_exclude = numpy.zeros(ttheta_in_range.shape, dtype=bool)
            if exclude is not None:
                for item_e in exclude.items:
                    flag_in_1 = numpy.logical_and(
                        ttheta_in_range >= item_e.ttheta_min,
                        ttheta_in_range <= item_e.ttheta_max)
                    flag_exclude = numpy.logical_or(flag_exclude, flag_in_1)
            ddict["excluded_points"] = flag_exclude

        if resolution is not None:
            ddict["resolution_parameters"] = numpy.array([
                resolution.u, resolution.v, resolution.w,
                resolution.x, resolution.y], dtype=float)

            ddict["flags_resolution_parameters"] = numpy.array([
                resolution.u_refinement, resolution.v_refinement, resolution.w_refinement,
                resolution.x_refinement, resolution.y_refinement], dtype=bool)

        if asymmetry is not None:
            ddict["asymmetry_parameters"] = numpy.array([
                asymmetry.p1, asymmetry.p2, asymmetry.p3,
                asymmetry.p4], dtype=float)

            ddict["flags_asymmetry_parameters"] = numpy.array([
                asymmetry.p1_refinement, asymmetry.p2_refinement, asymmetry.p3_refinement,
                asymmetry.p4_refinement], dtype=bool)
        else:
            ddict["asymmetry_parameters"] = numpy.array([
                0, 0, 0, 0], dtype=float)

            ddict["flags_asymmetry_parameters"] = numpy.array([
                False, False, False, False], dtype=bool)

        if diffrn_radiation is not None:
            beam_polarization = diffrn_radiation.polarization
            flipper_efficiency = diffrn_radiation.efficiency
            ddict["beam_polarization"] = numpy.array([beam_polarization], dtype=float)
            ddict["flipper_efficiency"] = numpy.array([flipper_efficiency], dtype=float)
            ddict["flags_beam_polarization"] = numpy.array([diffrn_radiation.polarization_refinement], dtype=bool)
            ddict["flags_flipper_efficiency"] = numpy.array([diffrn_radiation.efficiency_refinement], dtype=bool)

        if phase is not None:
            ddict["phase_name"] = numpy.array(phase.label, dtype=str)
            if phase.is_attribute("u"):
                p_u = numpy.array(phase.u, dtype=float)
                r_u = numpy.array(phase.u_refinement, dtype=bool)
            else:
                p_u = numpy.zeros((len(phase.items),), dtype=float)
                r_u = numpy.zeros((len(phase.items),), dtype=bool)

            if phase.is_attribute("v"):
                p_v = numpy.array(phase.v, dtype=float)
                r_v = numpy.array(phase.v_refinement, dtype=bool)
            else:
                p_v = numpy.zeros((len(phase.items),), dtype=float)
                r_v = numpy.zeros((len(phase.items),), dtype=bool)

            if phase.is_attribute("w"):
                p_w = numpy.array(phase.w, dtype=float)
                r_w = numpy.array(phase.w_refinement, dtype=bool)
            else:
                p_w = numpy.zeros((len(phase.items),), dtype=float)
                r_w = numpy.zeros((len(phase.items),), dtype=bool)

            if phase.is_attribute("x"):
                p_x = numpy.array(phase.x, dtype=float)
                r_x = numpy.array(phase.x_refinement, dtype=bool)
            else:
                p_x = numpy.zeros((len(phase.items),), dtype=float)
                r_x = numpy.zeros((len(phase.items),), dtype=bool)

            if phase.is_attribute("y"):
                p_y = numpy.array(phase.y, dtype=float)
                r_y = numpy.array(phase.y_refinement, dtype=bool)
            else:
                p_y = numpy.zeros((len(phase.items),), dtype=float)
                r_y = numpy.zeros((len(phase.items),), dtype=bool)

            ddict["phase_resolution_parameters"] = numpy.stack([p_u, p_v, p_w, p_x, p_y], axis=0)
            ddict["flags_phase_resolution_parameters"] = numpy.stack([r_u, r_v, r_w, r_x, r_y], axis=0)

            if phase.is_attribute("igsize"):
                ddict["phase_ig"] = numpy.array(phase.igsize, dtype=float)
                ddict["flags_phase_ig"] = numpy.array(phase.igsize_refinement, dtype=bool)
            else:
                ddict["phase_ig"] = numpy.zeros((len(phase.items),), dtype=float)
                ddict["flags_phase_ig"] = numpy.zeros((len(phase.items),), dtype=bool)

            if phase.is_attribute("scale"):
                ddict["phase_scale"] = numpy.array(phase.scale, dtype=float)
                ddict["flags_phase_scale"] = numpy.array(phase.scale_refinement, dtype=bool)
            else:
                ddict["phase_scale"] = numpy.zeros((len(phase.items),), dtype=float)
                ddict["flags_phase_scale"] = numpy.zeros((len(phase.items),), dtype=bool)

        if texture is not None:
            ddict["texture_name"] = numpy.array(texture.label, dtype=str)
            ddict["texture_g1"] = numpy.array(texture.g_1, dtype=float)
            ddict["texture_g2"] = numpy.array(texture.g_2, dtype=float)
            ddict["texture_axis"] = numpy.array(
                [texture.h_ax, texture.k_ax, texture.l_ax], dtype=float)
            ddict["flags_texture_g1"] = numpy.array(texture.g_1_refinement, dtype=bool)
            ddict["flags_texture_g2"] = numpy.array(texture.g_2_refinement, dtype=bool)
            ddict["flags_texture_axis"] = numpy.array(
                [texture.h_ax_refinement, 
                 texture.k_ax_refinement, 
                 texture.l_ax_refinement], dtype=bool)

        if pd_background is not None:
            ddict["background_ttheta"] = numpy.array(pd_background.ttheta, dtype=float)*numpy.pi/180.
            ddict["background_intensity"] = numpy.array(pd_background.intensity, dtype=float)
            ddict["flags_background_intensity"] = numpy.array(pd_background.intensity_refinement, dtype=bool)

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

        if "phase_scale" in keys:
            hh = ddict_diffrn["phase_scale"]
            for i_item, item in enumerate(self.phase.items):
                item.scale = float(hh[i_item])

        if "phase_ig" in keys:
            hh = ddict_diffrn["phase_ig"]
            for i_item, item in enumerate(self.phase.items):
                item.igsize = float(hh[i_item])

        if "resolution_parameters" in keys:
            hh = ddict_diffrn["resolution_parameters"]
            resolution = self.pd_instr_resolution 
            resolution.u = float(hh[0]) 
            resolution.v = float(hh[1]) 
            resolution.w = float(hh[2]) 
            resolution.x = float(hh[3]) 
            resolution.y = float(hh[4]) 

        if "asymmetry_parameters" in keys:
            hh = ddict_diffrn["asymmetry_parameters"]
            if self.is_attribute("pd_instr_reflex_asymmetry"):
                asymmetry = self.pd_instr_reflex_asymmetry
            else:
                asymmetry = PdInstrReflexAsymmetry()
                self.items.append(asymmetry)
            asymmetry.p1 = float(hh[0]) 
            asymmetry.p2 = float(hh[1]) 
            asymmetry.p3 = float(hh[2]) 
            asymmetry.p4 = float(hh[3])

        if "background_intensity" in keys:
            hh = ddict_diffrn["background_intensity"]
            for i_item, item in enumerate(self.pd_background.items):
                item.intensity = float(hh[i_item])

        if "offset_ttheta" in keys:
            self.setup.offset_ttheta = float(ddict_diffrn["offset_ttheta"]) * 180./numpy.pi

        if "wavelength" in keys:
            self.setup.wavelength = float(ddict_diffrn["wavelength"])

        if len(parameter_label) > 0:
            for name, sigma in zip(l_parameter_name, l_sigma):
                if name[0] == "background_intensity":
                    self.pd_background.items[name[1][0]].intensity_sigma = float(sigma)
                if name[0] == "phase_scale":
                    self.phase.items[name[1][0]].scale_sigma = float(sigma)
                if name[0] == "phase_ig":
                    self.phase.items[name[1][0]].igsize_sigma = float(sigma)
                if name[0] == "wavelength":
                    self.setup.wavelength_sigma = float(sigma)
                if name[0] == "offset_ttheta":
                    self.setup.offset_ttheta_sigma = float(sigma) * 180./numpy.pi
                if name[0] == "beam_polarization":
                    self.diffrn_radiation.polarization_sigma = float(sigma)
                if name[0] == "flipper_efficiency":
                    self.diffrn_radiation.efficiency_sigma = float(sigma)
                if name[0] == "asymmetry_parameters":
                    asymmetry = self.pd_instr_reflex_asymmetry
                    if name[1][0] == 0:
                        asymmetry.p1_sigma = float(sigma)
                    elif name[1][0] == 1:
                        asymmetry.p2_sigma = float(sigma)
                    elif name[1][0] == 2:
                        asymmetry.p3_sigma = float(sigma)
                    elif name[1][0] == 3:
                        asymmetry.p4_sigma = float(sigma)
                if name[0] == "resolution_parameters":
                    resolution = self.pd_instr_resolution 
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
            pd_proc = PdProcL()
            pd_proc.numpy_ttheta = numpy.round(ddict_diffrn["ttheta"] * 180./numpy.pi, decimals=5)
            pd_proc.numpy_ttheta_corrected = numpy.round(ddict_diffrn["ttheta_corrected"] * 180./numpy.pi, decimals=5)
            pd_proc.numpy_intensity_plus_net = numpy.round(ddict_diffrn["signal_plus"], decimals=5)
            pd_proc.numpy_intensity_minus_net = numpy.round(ddict_diffrn["signal_minus"], decimals=5)
            if "signal_exp_plus" in keys:
                pd_proc.numpy_intensity_plus = numpy.round(ddict_diffrn["signal_exp_plus"][0], decimals=5)
                pd_proc.numpy_intensity_plus_sigma = numpy.round(ddict_diffrn["signal_exp_plus"][1], decimals=5)
                pd_proc.numpy_intensity_minus = numpy.round(ddict_diffrn["signal_exp_minus"][0], decimals=5)
                pd_proc.numpy_intensity_minus_sigma = numpy.round(ddict_diffrn["signal_exp_minus"][1], decimals=5)
            else:
                pd_proc.numpy_intensity = numpy.round(ddict_diffrn["signal_exp"][0], decimals=5)
                pd_proc.numpy_intensity_sigma = numpy.round(ddict_diffrn["signal_exp"][1], decimals=5)
            pd_proc.numpy_intensity_bkg_calc = numpy.round(ddict_diffrn["signal_background"], decimals=5)
            pd_proc.numpy_excluded = ddict_diffrn["excluded_points"]
            pd_proc.numpy_to_items()
            self.pd_proc = pd_proc

        l_pd_peak = []
        l_refln = []
        for item_phase in self.phase.items:
            s_label = f"dict_in_out_{item_phase.label:}"
            if s_label in keys:
                dict_crystal = ddict_diffrn[s_label]
                dict_crystal_keys = dict_crystal.keys()
                if ("index_hkl" in dict_crystal_keys):
                    index_hkl = dict_crystal["index_hkl"]
                    pd_peak = PdPeakL(loop_name = item_phase.label)
                    pd_peak.numpy_index_h = index_hkl[0]
                    pd_peak.numpy_index_k = index_hkl[1]
                    pd_peak.numpy_index_l = index_hkl[2]

                    if "iint_plus_with_factors" in dict_crystal_keys:
                        pd_peak.numpy_intensity_plus = dict_crystal["iint_plus_with_factors"]
                    if "iint_minus_with_factors" in dict_crystal_keys:
                        pd_peak.numpy_intensity_minus = dict_crystal["iint_minus_with_factors"]
                    if "multiplicity_hkl" in dict_crystal_keys:
                        pd_peak.numpy_index_mult = dict_crystal["multiplicity_hkl"]
                    if "ttheta_hkl" in dict_crystal_keys:
                        pd_peak.numpy_ttheta = numpy.round(dict_crystal["ttheta_hkl"]* 180./numpy.pi, decimals=3)
                    pd_peak.numpy_to_items()
                    l_pd_peak.append(pd_peak)                 

                    flag_f_squared = False
                    if (("iint_plus" in dict_crystal_keys) and ("iint_minus" in dict_crystal_keys)):
                        f_squared_calc = dict_crystal["iint_plus"] + dict_crystal["iint_minus"]
                        f_polarized_squared_calc = dict_crystal["iint_plus"] - dict_crystal["iint_minus"]
                        flag_f_squared = True

                    if ((("f_nucl" in dict_crystal_keys) or ("f_charge" in dict_crystal_keys)) or flag_f_squared):
                        refln = ReflnL(loop_name = item_phase.label)
                        refln.numpy_index_h = index_hkl[0]
                        refln.numpy_index_k = index_hkl[1]
                        refln.numpy_index_l = index_hkl[2]

                        if "f_charge" in dict_crystal_keys:
                            refln.numpy_a_calc = dict_crystal["f_charge"].real
                            refln.numpy_b_calc = dict_crystal["f_charge"].imag
                        elif "f_nucl" in dict_crystal_keys:
                            refln.numpy_a_calc = dict_crystal["f_nucl"].real
                            refln.numpy_b_calc = dict_crystal["f_nucl"].imag

                        if flag_f_squared:
                            refln.numpy_f_squared_calc = f_squared_calc
                            refln.numpy_f_polarized_squared_calc = f_polarized_squared_calc

                        refln.numpy_to_items()
                        l_refln.append(refln)




        if len(l_pd_peak) > 0:
            self.add_items(l_pd_peak)
        if len(l_refln) > 0:
            self.add_items(l_refln)

    def define_background_points(self, step_ttheta: float = 1.):
        if self.is_attribute("pd_background"):
            pd_background = self.pd_background
        else:
            pd_background = PdBackgroundL()
            self.add_items(pd_background)
        pd_meas = self.pd_meas
        pd_background.define_points(pd_meas, step_ttheta)
        return 

    def estimate_background(self):
        pd_background = self.pd_background
        intensity_refinement = numpy.copy(numpy.array(
            pd_background.intensity_refinement, dtype=bool))
        pd_background.numpy_intensity_refinement = numpy.ones_like(intensity_refinement, dtype=bool)
        pd_background.numpy_to_items()
        setup=self.setup
        pd_proc = self.pd_proc
        res = pd_proc.estimate_background(
            pd_background,
            offset_ttheta=setup.offset_ttheta)
        pd_background.numpy_intensity_refinement = intensity_refinement
        pd_background.numpy_to_items()
        return res
