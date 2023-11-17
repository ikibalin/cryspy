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
# from cryspy.E_data_classes.cl_1_mag_crystal import MagCrystal


class Pd(DataN):
    """
    Powder diffraction experiment with polarized or unpolarized neutrons (1d).

    Data items in the DIFFRN category record details about
    powder diffraction measurements.

    Methods
    -------
        - calc_profile
        - calc_chi_sq
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
        ddict = super(Pd, self).get_dictionary()
        pd_meas, range_, exclude = None, None, None

        l_obj = take_items_by_class(self, (PdMeasL, ))
        if len(l_obj) > 0:
            pd_meas = l_obj[0]

        l_obj = take_items_by_class(self, (Range, ))
        if len(l_obj) > 0:
            range_ = l_obj[0]

        l_obj = take_items_by_class(self, (ExcludeL, ))
        if len(l_obj) > 0:
            exclude = l_obj[0]

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

        if "texture_g1" in keys:
            pd_texture_g1 = ddict_diffrn["texture_g1"]
            pd_texture_g2 = ddict_diffrn["texture_g2"]
            pd_texture_axis = ddict_diffrn["texture_axis"]
            for i_item, item in enumerate(self.texture.items):
                item.g_1 = float(pd_texture_g1[i_item])
                item.g_2 = float(pd_texture_g2[i_item])
                item.h_ax = float(pd_texture_axis[0, i_item])
                item.k_ax = float(pd_texture_axis[1, i_item])
                item.l_ax = float(pd_texture_axis[2, i_item])

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
                if name[0] == "texture_g1":
                    self.texture.items[name[1][0]].g_1_sigma = float(sigma)
                if name[0] == "texture_g2":
                    self.texture.items[name[1][0]].g_2_sigma = float(sigma)
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
