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
from cryspy.C_item_loop_classes.cl_1_tof_background_by_points import TOFBackgroundPointL
from cryspy.C_item_loop_classes.cl_1_tof_intensity_incident import \
    TOFIntensityIncident
from cryspy.C_item_loop_classes.cl_1_tof_meas import TOFMeasL
from cryspy.C_item_loop_classes.cl_1_tof_parameters import TOFParameters
from cryspy.C_item_loop_classes.cl_1_tof_peak import TOFPeakL
from cryspy.C_item_loop_classes.cl_1_tof_proc import TOFProcL
from cryspy.C_item_loop_classes.cl_1_tof_profile import TOFProfile

from cryspy.E_data_classes.cl_1_crystal import Crystal
# from cryspy.E_data_classes.cl_1_mag_crystal import MagCrystal


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

    CLASSES_MANDATORY = (TOFParameters, TOFProfile, PhaseL, TOFMeasL)
    CLASSES_OPTIONAL = (TOFBackground, TOFBackgroundPointL, DiffrnRadiation, Chi2, Range,
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
                    ax_hkl.plot(np_time, 0.*np_time+y_min_s -
                                y_step_s, "|", label=item.loop_name)
                    y_step_s += 0.05*y_dist_s
            res = []
            ax_s.legend(loc='upper right')
            res.append((fig_s, ax_s))
            return res
        elif self.is_attribute("tof_meas"):
            return self.tof_meas.plots()

    def get_dictionary(self):
        dict_tof = super(TOF, self).get_dictionary()

        tof_meas, range_, exclude = None, None, None

        l_obj = take_items_by_class(self, (Range, ))
        if len(l_obj) > 0:
            range_ = l_obj[0]

        l_obj = take_items_by_class(self, (ExcludeL, ))
        if len(l_obj) > 0:
            exclude = l_obj[0]

        l_obj = take_items_by_class(self, (TOFMeasL, ))
        if len(l_obj) > 0:
            tof_meas = l_obj[0]

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
                int_plus = numpy.array(
                    tof_meas.intensity_plus, dtype=float)[flag_in]
                s_int_plus = numpy.array(
                    tof_meas.intensity_plus_sigma, dtype=float)[flag_in]
                int_minus = numpy.array(
                    tof_meas.intensity_minus, dtype=float)[flag_in]
                s_int_minus = numpy.array(
                    tof_meas.intensity_minus_sigma, dtype=float)[flag_in]
                dict_tof["signal_exp_plus"] = numpy.stack(
                    [int_plus, s_int_plus], axis=0)
                dict_tof["signal_exp_minus"] = numpy.stack(
                    [int_minus, s_int_minus], axis=0)
            else:
                int_sum = numpy.array(tof_meas.intensity, dtype=float)[flag_in]
                s_int_sum = numpy.array(
                    tof_meas.intensity_sigma, dtype=float)[flag_in]
                dict_tof["signal_exp"] = numpy.stack(
                    [int_sum, s_int_sum], axis=0)

            time_in_range = time[flag_in]
            flag_exclude = numpy.zeros(time_in_range.shape, dtype=bool)
            if exclude is not None:
                for item_e in exclude.items:
                    flag_in_1 = numpy.logical_and(
                        time_in_range >= item_e.time_min,
                        time_in_range <= item_e.time_max)
                    flag_exclude = numpy.logical_or(flag_exclude, flag_in_1)
            dict_tof["excluded_points"] = flag_exclude

        return dict_tof

    def take_parameters_from_dictionary(self, ddict_diffrn, l_parameter_name: list = None, l_sigma: list = None):
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

        if "background_intensity" in keys:
            hh = ddict_diffrn["background_intensity"]
            for i_item, item in enumerate(self.tof_backgroundpoint.items):
                item.intensity = float(hh[i_item])

        if "profile_peak_shape" in keys:
            profile_peak_shape = ddict_diffrn["profile_peak_shape"]
            if profile_peak_shape == "pseudo-Voigt":
                profile_sigmas = ddict_diffrn["profile_sigmas"]
                profile_alphas = ddict_diffrn["profile_alphas"]
                profile_betas = ddict_diffrn["profile_betas"]
                profile_gammas = ddict_diffrn["profile_gammas"]
                tof_profile = self.tof_profile
                setattr(tof_profile, "alpha0", profile_alphas[0])
                setattr(tof_profile, "alpha1", profile_alphas[1])
                setattr(tof_profile, "beta0", profile_betas[0])
                setattr(tof_profile, "beta1", profile_betas[1])
                setattr(tof_profile, "sigma0", profile_sigmas[0])
                setattr(tof_profile, "sigma1", profile_sigmas[1])
                setattr(tof_profile, "sigma2", profile_sigmas[2])
                setattr(tof_profile, "gamma0", profile_gammas[0])
                setattr(tof_profile, "gamma1", profile_gammas[1])
                setattr(tof_profile, "gamma2", profile_gammas[2])
            elif profile_peak_shape == "Gauss":
                profile_sigmas = ddict_diffrn["profile_sigmas"]
                profile_alphas = ddict_diffrn["profile_alphas"]
                profile_betas = ddict_diffrn["profile_betas"]
                tof_profile = self.tof_profile
                setattr(tof_profile, "alpha0", profile_alphas[0])
                setattr(tof_profile, "alpha1", profile_alphas[1])
                setattr(tof_profile, "beta0", profile_betas[0])
                setattr(tof_profile, "beta1", profile_betas[1])
                setattr(tof_profile, "sigma0", profile_sigmas[0])
                setattr(tof_profile, "sigma1", profile_sigmas[1])
                setattr(tof_profile, "sigma2", profile_sigmas[2])
            elif profile_peak_shape == "type0m":
                profile_alphas = ddict_diffrn["profile_alphas"]
                profile_betas = ddict_diffrn["profile_betas"]
                profile_sigmas = ddict_diffrn["profile_sigmas"]
                profile_gammas = ddict_diffrn["profile_gammas"]
                profile_rs = ddict_diffrn["profile_rs"]
                tof_profile = self.tof_profile
                setattr(tof_profile, "alpha1", profile_alphas[0])
                setattr(tof_profile, "alpha2", profile_alphas[1])
                setattr(tof_profile, "beta00", profile_betas[0])
                setattr(tof_profile, "beta01", profile_betas[1])
                setattr(tof_profile, "beta10", profile_betas[2])
                setattr(tof_profile, "sigma0", profile_sigmas[0])
                setattr(tof_profile, "sigma1", profile_sigmas[1])
                setattr(tof_profile, "sigma2", profile_sigmas[2])
                setattr(tof_profile, "gamma0", profile_gammas[0])
                setattr(tof_profile, "gamma1", profile_gammas[1])
                setattr(tof_profile, "gamma2", profile_gammas[2])
                setattr(tof_profile, "r01", profile_rs[0])
                setattr(tof_profile, "r02", profile_rs[1])
                setattr(tof_profile, "r03", profile_rs[2])

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
                if name[0] == "profile_sigmas":
                    if name[1] == 0:
                        self.tof_profile.sigma0_sigma = float(sigma)
                    elif name[1] == 1:
                        self.tof_profile.sigma1_sigma = float(sigma)
                    elif name[1] == 2:
                        self.tof_profile.sigma2_sigma = float(sigma)
                if name[0] == "profile_gammas":
                    if name[1] == 0:
                        self.tof_profile.gamma0_sigma = float(sigma)
                    elif name[1] == 1:
                        self.tof_profile.gamma1_sigma = float(sigma)
                    elif name[1] == 2:
                        self.tof_profile.gamma2_sigma = float(sigma)

        if (("signal_plus" in keys) and ("signal_minus" in keys)):
            tof_proc = TOFProcL()
            tof_proc.numpy_time = numpy.round(ddict_diffrn["time"], decimals=5)
            tof_proc.numpy_intensity_plus_net = numpy.round(
                ddict_diffrn["signal_plus"], decimals=5)
            tof_proc.numpy_intensity_minus_net = numpy.round(
                ddict_diffrn["signal_minus"], decimals=5)
            if "signal_exp_plus" in keys:
                tof_proc.numpy_intensity_plus = numpy.round(
                    ddict_diffrn["signal_exp_plus"][0], decimals=5)
                tof_proc.numpy_intensity_plus_sigma = numpy.round(
                    ddict_diffrn["signal_exp_plus"][1], decimals=5)
                tof_proc.numpy_intensity_minus = numpy.round(
                    ddict_diffrn["signal_exp_minus"][0], decimals=5)
                tof_proc.numpy_intensity_minus_sigma = numpy.round(
                    ddict_diffrn["signal_exp_minus"][1], decimals=5)
            else:
                tof_proc.numpy_intensity = numpy.round(
                    ddict_diffrn["signal_exp"][0], decimals=5)
                tof_proc.numpy_intensity_sigma = numpy.round(
                    ddict_diffrn["signal_exp"][1], decimals=5)
            tof_proc.numpy_intensity_bkg_calc = numpy.round(
                ddict_diffrn["signal_background"], decimals=5)
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
                if (("index_hkl" in dict_crystal_keys) and ("time_hkl" in dict_crystal_keys)):
                    index_hkl = dict_crystal["index_hkl"]
                    time_hkl = dict_crystal["time_hkl"]

                    tof_peak = TOFPeakL(loop_name=item_phase.label)
                    int_plus_max, int_minus_max = None, None
                    if "iint_plus_with_factors" in dict_crystal_keys:
                        # int_plus_max = dict_crystal["iint_plus_with_factors"]
                        int_plus_max = dict_crystal["iint_plus"]

                    if "iint_minus_with_factors" in dict_crystal_keys:
                        # int_minus_max = dict_crystal["iint_minus_with_factors"]
                        int_minus_max = dict_crystal["iint_minus"]

                    if "f_nucl":
                        refln = ReflnL(loop_name=item_phase.label)
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
