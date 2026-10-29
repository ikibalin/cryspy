"""Description of Pd2d class."""
__author__ = 'ikibalin'
__version__ = "2020_08_19"
import numpy
from typing import NoReturn


from cryspy.A_functions_base.powder_diffraction_const_wavelength import \
    calc_ttheta_phi_by_gamma_nu

from cryspy.B_parent_classes.cl_3_data import DataN
from cryspy.D_functions_item_loop.function_1_report_powder import \
    report_powder_experiment
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


na = numpy.newaxis

class Pd2d(DataN):
    """
    Powder diffraction experiment with polarized or unpolarized neutrons (2d).

    Data items in the DIFFRN category record details about
    2d powder diffraction measurements.

    Methods
    -------
        - estimate_background
        - plots
        - report


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

    def apply_constraints(self):
        """Apply constraints."""
        pass

    def report(self):
        """Agreement with the model, phases and refinement summary."""
        return report_powder_experiment(self)

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
        elif self.is_attribute("pd2d_meas"):
            return self.pd2d_meas.plots()
        return []

    def get_dictionary(self):
        """Form dictionary. See documentation moduel CrysPy using Jupyter notebook.
        """
        ddict = super(Pd2d, self).get_dictionary()
        pd2d_meas, range_, exclude = None, None, None

        l_obj = take_items_by_class(self, (Pd2dMeas, ))
        if len(l_obj) > 0:
            pd2d_meas = l_obj[0]


        l_obj = take_items_by_class(self, (Range, ))
        if len(l_obj) > 0:
            range_ = l_obj[0]

        l_obj = take_items_by_class(self, (ExcludeL, ))
        if len(l_obj) > 0:
            exclude = l_obj[0]


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
            self.setup.offset_gamma = float(ddict_diffrn["offset_gamma"].squeeze()) * 180./numpy.pi

        if "offset_nu" in keys:
            self.setup.offset_nu = float(ddict_diffrn["offset_nu"].squeeze()) * 180./numpy.pi

        if "wavelength" in keys:
            self.setup.wavelength = float(ddict_diffrn["wavelength"].squeeze())

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
