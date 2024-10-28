"""Description of Diffrn class."""
__author__ = 'ikibalin'
__version__ = "2020_08_19"

from warnings import warn
from typing import NoReturn, List, Union
import numpy
import scipy
import scipy.misc

from cryspy.A_functions_base.matrix_operations import calc_vector_product_v1_v2_v1

from cryspy.A_functions_base.extinction import calc_extinction_sphere
from cryspy.A_functions_base.flip_ratio import calc_flip_ratio_by_structure_factors
from cryspy.A_functions_base.unit_cell import calc_sthovl_by_unit_cell_parameters, calc_volume_uc_by_unit_cell_parameters, \
    calc_eq_ccs_by_unit_cell_parameters
from cryspy.A_functions_base.structure_factor import calc_f_m_perp_by_sft

from cryspy.A_functions_base.function_1_roots import calc_roots

from cryspy.B_parent_classes.cl_2_loop import LoopN
from cryspy.B_parent_classes.cl_3_data import DataN
from cryspy.B_parent_classes.preocedures import take_items_by_class

from cryspy.C_item_loop_classes.cl_1_chi2 import Chi2
from cryspy.C_item_loop_classes.cl_1_setup import Setup
from cryspy.C_item_loop_classes.cl_1_diffrn_radiation import \
    DiffrnRadiation
from cryspy.C_item_loop_classes.cl_2_diffrn_orient_matrix import \
    DiffrnOrientMatrix
from cryspy.C_item_loop_classes.cl_1_diffrn_refln import DiffrnReflnL
from cryspy.C_item_loop_classes.cl_1_extinction import Extinction
from cryspy.C_item_loop_classes.cl_1_phase import Phase
from cryspy.C_item_loop_classes.cl_1_refln import ReflnL
from cryspy.C_item_loop_classes.cl_1_refln_susceptibility import \
    ReflnSusceptibilityL
from cryspy.C_item_loop_classes.cl_1_refine_ls import RefineLs

from cryspy.E_data_classes.cl_1_crystal import Crystal
# from cryspy.E_data_classes.cl_1_mag_crystal import MagCrystal


class Diffrn(DataN):
    """
    Single diffraction experiment with polarized neutrons.

    Data items in the DIFFRN category record details about
    single diffraction measurements.

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
        - diffrn_radiation (mandatory)
        - diffrn_orient_matrix (mandatory)
        - diffrn_refln (mandatory)
        - extinction
        - phase
        - refln_#phase_name
        - refln_susceptibility_#phase_name
        - refine_ls
    """

    CLASSES_MANDATORY = (Setup, DiffrnOrientMatrix,
                         DiffrnReflnL)
    CLASSES_OPTIONAL = (DiffrnRadiation, Extinction, Phase, ReflnL, ReflnSusceptibilityL,
                        RefineLs, Chi2)
    # CLASSES_INTERNAL = ()

    CLASSES = CLASSES_MANDATORY + CLASSES_OPTIONAL

    PREFIX = "diffrn"

    # default values for the parameters
    D_DEFAULT = {}

    def __init__(self, data_name=None, **kwargs) -> NoReturn:
        super(Diffrn, self).__init__()

        self.__dict__["items"] = []
        self.__dict__["data_name"] = data_name

        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)

    def estimate_FM(self, crystal, max_abs_f_mag: float = 5,
                    delta_f_mag: float = 0.001, flag_lambdaover2_correction=False):
        """
        Estimate magnetic structure factor.

        The method calculates magnetic structure factors (only real part) from
        flipping ratios based on give crystal structure.
        The crystal structure should be centrosymmetric
        """
        setup = self.setup
        wavelength = setup.wavelength
        field_norm = setup.field

        diffrn_refln = self.diffrn_refln
        index_h = diffrn_refln.numpy_index_h
        index_k = diffrn_refln.numpy_index_k
        index_l = diffrn_refln.numpy_index_l
        index_hkl = numpy.stack([index_h, index_k, index_l], axis=0)
        fr_exp = diffrn_refln.numpy_fr
        fr_sigma = diffrn_refln.numpy_fr_sigma

        diffrn_radiation = self.diffrn_radiation
        polarization = diffrn_radiation.polarization
        flipper_efficiency = diffrn_radiation.efficiency

        if self.is_attribute("extinction"):
            extinction = self.extinction
            model_extinction = extinction.model
            radius = extinction.radius
            mosaicity = extinction.mosaicity
        else:
            model_extinction = ""
            radius = None
            mosaicity = None

        diffrn_orient_matrix = self.diffrn_orient_matrix
        u_matrix = diffrn_orient_matrix.u
        magnetic_field = field_norm * u_matrix[2, :]

        refln = crystal.calc_refln(index_hkl)
        f_nucl = numpy.array(refln.f_calc, dtype=complex)

        cell = crystal.cell
        unit_cell_parameters = cell.get_unit_cell_parameters()
        sthovl = cell.calc_sthovl(index_hkl[0], index_hkl[1], index_hkl[2])
        eq_ccs, dder_eq_ccs = calc_eq_ccs_by_unit_cell_parameters(index_hkl, unit_cell_parameters, flag_unit_cell_parameters=False)


        if (setup.is_attribute("ratio_lambdaover2") &
                flag_lambdaover2_correction):
            ratio_lambdaover2 = setup.ratio_lambdaover2
            flag_lambdaover2 = True
        else:
            ratio_lambdaover2 = None
            flag_lambdaover2 = False

        if (flag_lambdaover2 & flag_lambdaover2_correction):
            index_2hkl = 2 * index_hkl
            f_nucl_2hkl = crystal.calc_f_nucl(index_2hkl)
        else:
            f_nucl_2hkl = None

        if (flag_lambdaover2 & flag_lambdaover2_correction):
            sft_2hkl, dder_sft_2hkl = crystal.calc_structure_factor_tensor_ccs(index_2hkl)

            fm_perp_loc_2hkl, dder_f_m_perp = calc_f_m_perp_by_sft(
                sft_2hkl, magnetic_field, eq_ccs,
                flag_sft_ccs=False,
                flag_magnetic_field=False,
                flag_eq_ccs=False)

        else:
            fm_perp_loc_2hkl = None

        def func_extinction(f_sq, flag_f_sq: bool = False):
            model = self.extinction.model
            radius = self.extinction.radius
            mosaicity = self.extinction.mosaicity
            volume_unit_cell = calc_volume_uc_by_unit_cell_parameters(unit_cell_parameters)[0]
            wavelength = self.setup.wavelength
            cos_2theta = 1.-2*numpy.square(sthovl * wavelength)
            flag_radius = True
            flag_mosaicity = True
            flag_volume_unit_cell = True
            flag_cos_2theta = True
            flag_wavelength = True
            return calc_extinction_sphere(
                f_sq, radius, mosaicity, volume_unit_cell, cos_2theta, wavelength,
                model, flag_f_sq=flag_f_sq, flag_radius=flag_radius,
                flag_mosaicity=flag_mosaicity,
                flag_volume_unit_cell=flag_volume_unit_cell,
                flag_cos_2theta=flag_cos_2theta,
                flag_wavelength=flag_wavelength)

        def function_to_find_fm(f_mag):
            """Calc flip ratio.

            Parameters
            ----------
            f_mag : float
                in muB.

            Returns
            -------
            flip_ratio : float
                flip ratio.
            """
            mag_loc_1 = u_matrix[2, 0]*f_mag*0.2695
            mag_loc_2 = u_matrix[2, 1]*f_mag*0.2695
            mag_loc_3 = u_matrix[2, 2]*f_mag*0.2695
            mag_loc = numpy.stack([mag_loc_1, mag_loc_2, mag_loc_3], axis=0)
            f_m_perp, dder_fm_perp_loc = calc_vector_product_v1_v2_v1(eq_ccs, mag_loc, flag_v1=False, flag_v2=False)

            axis_z = u_matrix[2, :]
            f_r, dder_f_r = calc_flip_ratio_by_structure_factors(
                       polarization, flipper_efficiency, f_nucl, f_m_perp, axis_z, func_extinction=func_extinction,
                       c_lambda2=ratio_lambdaover2, f_n_2h=f_nucl_2hkl, f_m_perp_2h=fm_perp_loc_2hkl,
                       flag_polarization=False, flag_flipper=False,
                       flag_f_n=False, flag_f_m_perp=False,
                       flag_c_lambda2=False,
                       flag_f_n_2h=False, flag_f_m_perp_2h=False)

            return f_r

        lFMans, lFMansMin = [], []
        lFMansMin_sigma = []
        for _ind in range(fr_exp.size):
            def f(x):
                return (function_to_find_fm(x)[_ind] - fr_exp[_ind])
            FMans = calc_roots(f, -1*max_abs_f_mag-delta_f_mag, max_abs_f_mag+delta_f_mag)
            lFMans.append(FMans)
            if len(FMans) == 0:
                lFMansMin.append(None)
                lFMansMin_sigma.append(None)
            else:
                FMans_min = min(FMans)
                lFMansMin.append(FMans_min)
                der1 = scipy.misc.derivative(f, FMans_min, 0.1, 1)
                f_mag_sigma = abs(fr_sigma[_ind] / der1)
                lFMansMin_sigma.append(f_mag_sigma)

        ls_out = [f"loop_{crystal.data_name:}"]
        ls_out.append("_estimate_index_h\n_estimate_index_k")
        ls_out.append("_estimate_index_l\n_estimate_F_M\n_estimate_F_M_sigma")
        for ind_h, ind_k, ind_l, f_mag, f_sig in \
                zip(index_h, index_k, index_l, lFMansMin, lFMansMin_sigma):
            ls_out.append(f"{ind_h:} {ind_k:} {ind_l:} {f_mag:}  {f_sig:}")

        res = LoopN.from_cif("\n".join(ls_out))

        return res

    def apply_constraints(self):
        """Apply constraints."""
        pass

    def plots(self):
        l_res = []
        if self.is_attribute("diffrn_refln"):
            diffrn_refln = self.diffrn_refln
            flag_asymmetry = False
            if self.is_attribute("chi2"):
                chi2 = self.chi2
                if chi2.is_attribute("asymmetry"):
                    flag_asymmetry = chi2.asymmetry

            if flag_asymmetry:
                fig_ax = diffrn_refln.plot_asymmetry_vs_asymmetry_calc()
            elif diffrn_refln.is_attribute("fr"):
                fig_ax = diffrn_refln.plot_fr_vs_fr_calc()
            else:
                fig_ax = diffrn_refln.plot_intensity_vs_intensity_calc()

            if fig_ax is not None:
                fig, ax = fig_ax
                ax.set_title(self.data_name + " - "+ax.title.get_text())
                l_res.append((fig, ax))
        return l_res

    def report(self):
        if self.is_attribute("diffrn_refln"):
            diffrn_refln = self.diffrn_refln
            rep_1 = diffrn_refln.report_agreement_factor_exp()
            rep_2 = ""
            if self.is_attribute("diffrn_orient_matrix"):
                ub_matrix = self.diffrn_orient_matrix
                if ub_matrix.is_defined():
                    cell = ub_matrix.cell
                    rep_2 = diffrn_refln.report_chi_sq_exp(cell=cell)
            if rep_2 == "":
                rep_2 = diffrn_refln.report_chi_sq_exp()
            return rep_1 + "\n" + rep_2
        return  ""

    def get_dictionary(self):
        """Form dictionary. See documentation moduel CrysPy using Jupyter notebook.
        """
        ddict = super(Diffrn, self).get_dictionary()

        keys = ddict.keys()
        if (("magnetic_field" in keys) and ("matrix_u" in keys)):
            ddict["magnetic_field"] = ddict["magnetic_field"]*ddict["matrix_u"][6:]

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
        if "extinction_mosaicity" in keys:
            self.extinction.mosaicity = ddict_diffrn["extinction_mosaicity"]
        if "extinction_radius" in keys:
            self.extinction.radius = ddict_diffrn["extinction_radius"]
        if "phase_scale" in keys:
            self.phase.scale = ddict_diffrn["phase_scale"]


        if len(parameter_label) > 0:
            for name, sigma in zip(l_parameter_name, l_sigma):
                if name[0] == "wavelength":
                    self.setup.wavelength_sigma = float(sigma)
                if name[0] == "beam_polarization":
                    self.diffrn_radiation.polarization_sigma = float(sigma)
                if name[0] == "flipper_efficiency":
                    self.diffrn_radiation.efficiency_sigma = float(sigma)
                if name[0] == "extinction_mosaicity":
                    self.extinction.mosaicity_sigma = float(sigma)
                if name[0] == "extinction_radius":
                    self.extinction.radius_sigma = float(sigma)
                if name[0] == "phase_scale":
                    self.phase.scale_sigma = float(sigma)

        if (("iint_plus" in keys) and ("iint_minus" in keys)):
            diffrn_refln = self.diffrn_refln
            if diffrn_refln.is_attribute("fr"):
                diffrn_refln.numpy_fr_calc = ddict_diffrn["iint_plus"]/ddict_diffrn["iint_minus"]
            else:
                diffrn_refln.numpy_intensity_calc = ddict_diffrn["intensity_calc"]
            diffrn_refln.numpy_to_items()

        l_dict_diffrn_out_keys = [hh for hh in keys if hh.startswith("dict_in_out_hkl_")]
        for dict_diffrn_out_keys in l_dict_diffrn_out_keys:
            phase_name = dict_diffrn_out_keys[len("dict_in_out_hkl_"):]
            dict_in_out_hkl = ddict_diffrn[dict_diffrn_out_keys]
            dict_in_out_hkl_keys = dict_in_out_hkl.keys()
            l_refln = []
            if (("f_nucl" in dict_in_out_hkl_keys) and ("index_hkl" in dict_in_out_hkl_keys)):
                index_hkl = dict_in_out_hkl["index_hkl"]
                refln = ReflnL(loop_name = phase_name)
                refln.numpy_index_h = index_hkl[0]
                refln.numpy_index_k = index_hkl[1]
                refln.numpy_index_l = index_hkl[2]
                refln.numpy_a_calc = dict_in_out_hkl["f_nucl"].real
                refln.numpy_b_calc = dict_in_out_hkl["f_nucl"].imag
                refln.numpy_sintlambda = dict_in_out_hkl["sthovl"]
                refln.numpy_to_items()
                l_refln.append(refln)
            self.add_items(l_refln)