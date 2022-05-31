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
from cryspy.E_data_classes.cl_1_mag_crystal import MagCrystal


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

    CLASSES_MANDATORY = (Setup, DiffrnRadiation, DiffrnOrientMatrix,
                         DiffrnReflnL)
    CLASSES_OPTIONAL = (Extinction, Phase, ReflnL, ReflnSusceptibilityL,
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

    def calc_iint_u_d_flip_ratio(
            self, index_hkl, l_crystal: List[Crystal], # or list of MagCrystal
            dict_in_out: dict = None,
            flag_internal: bool = True):
        """
        Calculate intensity for the given diffraction angle.

        Keyword Arguments
        -----------------
            h, k, l: 1D numpy array of Miller indexes
            l_crystal: list of crystal structures

        Output
        ------
            iint_u: intensity up (without scale factor)
            iint_d: intensity down (without scale factor)
            flip_ratio: flip ratio (intensity_plus/intensity_minus)

        """
        if dict_in_out is None:
            flag_dict = False
            dict_in_out_keys = []
        else:
            flag_dict = True
            dict_in_out_keys = dict_in_out.keys()

        l_label = [_.data_name for _ in l_crystal]
        if self.phase.label in l_label:
            ind = l_label.index(self.phase.label)
            crystal = l_crystal[ind]
        elif self.data_name in l_label:
            self.phase = Phase(label=self.data_name)
            ind = l_label.index(self.data_name)
            crystal = l_crystal[ind]
        else:
            warn("Crystal    not found. The first one is taken.")
            crystal = l_crystal[0]
            self.phase.label = crystal.label

        cell = crystal.cell
        unit_cell_parameters = cell.get_unit_cell_parameters()
        volume_unit_cell = cell.volume

        sthovl, dder_sthovl = calc_sthovl_by_unit_cell_parameters(index_hkl, unit_cell_parameters, flag_unit_cell_parameters=False)
        eq_ccs, dder_eq_ccs  = calc_eq_ccs_by_unit_cell_parameters(index_hkl, unit_cell_parameters, flag_unit_cell_parameters=False)

        setup = self.setup
        wavelength = setup.wavelength
        field_norm = setup.field
        try:
            ratio_lambdaover2 = setup.ratio_lambdaover2
            flag_lambdaover2 = True
        except AttributeError:
            ratio_lambdaover2 = None
            flag_lambdaover2 = False

        diffrn_orient_matrix = self.diffrn_orient_matrix
        u_matrix = diffrn_orient_matrix.u
        magnetic_field = field_norm*u_matrix[2, :]
        matrix_u = numpy.array([
            u_matrix[0,0], u_matrix[0,1], u_matrix[0,2],
            u_matrix[1,0], u_matrix[1,1], u_matrix[1,2],
            u_matrix[2,0], u_matrix[2,1], u_matrix[2,2]], dtype=float)
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
        if flag_lambdaover2:
            index_2hkl = 2 * index_hkl

        if flag_internal:
            refln = crystal.calc_refln(index_hkl)
            f_nucl = numpy.array(refln.f_calc, dtype=complex)
        else:
            f_nucl = crystal.calc_f_nucl(index_hkl)

        if flag_dict:
            dict_in_out["phase_name"] = crystal.get_name()
            dict_in_out["index_hkl"] = index_hkl
            dict_in_out["f_nucl"] = f_nucl
        
        if flag_lambdaover2:
            f_nucl_2hkl = crystal.calc_f_nucl(index_2hkl)
            if flag_dict:
                dict_in_out["f_nucl_2hkl"] = f_nucl_2hkl
        else:
            f_nucl_2hkl = None
        
        refln_s = None
        if isinstance(crystal, Crystal):
            if flag_internal:
                refln_s = crystal.calc_refln_susceptibility(index_hkl)
                sft_ccs = numpy.stack([
                    refln_s.numpy_chi_11_calc, refln_s.numpy_chi_12_calc, refln_s.numpy_chi_13_calc,
                    refln_s.numpy_chi_21_calc, refln_s.numpy_chi_22_calc, refln_s.numpy_chi_23_calc,
                    refln_s.numpy_chi_31_calc, refln_s.numpy_chi_32_calc, refln_s.numpy_chi_33_calc],
                    axis=0)
                fm_perp_loc, dder_fm_perp_loc = calc_f_m_perp_by_sft(
                    sft_ccs, magnetic_field, eq_ccs, flag_magnetic_field=False, flag_eq_ccs=False, flag_sft_ccs=False) 
            else:
                fm_perp_loc = crystal.calc_f_m_perp(index_hkl, magnetic_field, dict_in_out=dict_in_out)
                if flag_dict:
                    sft_ccs = dict_in_out["sft_ccs"]
            if flag_dict:
                dict_in_out["f_m_perp"] = fm_perp_loc
        elif isinstance(crystal, MagCrystal):
            fm_perp_loc = crystal.calc_f_mag_perp(index_hkl)

        if flag_lambdaover2:
            if isinstance(crystal, Crystal):
                dict_in_out_2 = {"phase_name": crystal.get_name()}
                fm_perp_loc_2hkl = crystal.calc_f_m_perp(index_2hkl, magnetic_field, dict_in_out=dict_in_out_2)
                if flag_dict:
                    dict_in_out["sft_ccs_2hkl"] = dict_in_out_2.pop("sft_ccs")
            elif isinstance(crystal, MagCrystal):
                fm_perp_loc_2hkl = crystal.calc_f_mag_perp(index_2hkl)
                if flag_dict:
                    dict_in_out["f_m_perp_2hkl"] = fm_perp_loc_2hkl
        else:
            fm_perp_loc_2hkl = None


        def func_extinction(f_sq, flag_f_sq: bool = False):
            cos_2theta = 1.-2*numpy.square(sthovl * wavelength)
            flag_radius = False
            flag_mosaicity = False
            flag_volume_unit_cell = False
            flag_cos_2theta = False
            flag_wavelength = False
            return calc_extinction_sphere(
                f_sq, radius, mosaicity, volume_unit_cell, cos_2theta, wavelength,
                model_extinction, flag_f_sq=flag_f_sq, flag_radius=flag_radius,
                flag_mosaicity=flag_mosaicity,
                flag_volume_unit_cell=flag_volume_unit_cell,
                flag_cos_2theta=flag_cos_2theta,
                flag_wavelength=flag_wavelength)

        flag_polarization = False
        flag_flipper = False
        flag_f_n = False
        flag_f_m_perp = False
        flag_c_lambda2 = False
        flag_f_n_2h = False 
        flag_f_m_perp_2h = False
        flip_ratio, dder_fr = calc_flip_ratio_by_structure_factors(
            polarization, flipper_efficiency, f_nucl, fm_perp_loc, matrix_u,
            func_extinction=func_extinction,
            c_lambda2=ratio_lambdaover2, f_nucl_2hkl=f_nucl_2hkl, f_m_perp_2hkl=fm_perp_loc_2hkl,
            flag_beam_polarization=flag_polarization, flag_flipper_efficiency=flag_flipper,
            flag_f_nucl=flag_f_n, flag_f_m_perp=flag_f_m_perp,
            flag_c_lambda2=flag_c_lambda2,
            flag_f_nucl_2hkl=flag_f_n_2h, flag_f_m_perp_2hkl=flag_f_m_perp_2h,
            dict_in_out=dict_in_out, flag_use_precalculated_data=flag_internal)


        dder = {}

        if flag_internal:
            refln.loop_name = crystal.data_name
            self.add_items([refln, ])
            if refln_s is not None:
                refln_s.loop_name = crystal.data_name
                self.add_items([refln_s, ])


        return flip_ratio, dder


    def calc_fr(self, cell, f_nucl, f_m_perp, delta_f_nucl=None,
                delta_f_m_perp=None):
        """
        Calculate flip ratio.

        By given nuclear structure factor and
        perpendicular component of magnetic structure factor
        (in 10**-12 cm)
        cell (only for extinction correction)
        f_nucl [hkl]
        f_m_perp [3, hkl]
        delta_f_nucl [hkl, parameters]
        delta_f_m_perp [hkl, parameters]
        """


        diffrn_radiation = self.diffrn_radiation
        polarization = diffrn_radiation.polarization
        flipper = diffrn_radiation.efficiency
        m_u = self.diffrn_orient_matrix.u

        unit_cell_parameters = cell.get_unit_cell_parameters()
        index_hkl = numpy.array([
            self.diffrn_refln.index_h,
            self.diffrn_refln.index_k,
            self.diffrn_refln.index_l], dtype=int)
        sthovl = calc_sthovl_by_unit_cell_parameters(index_hkl, unit_cell_parameters)[0]

        flag_polarization = False
        flag_flipper = False
        flag_f_n = False
        flag_f_m_perp = True

        c_lambda2 = None
        f_n_2h = None
        f_m_perp_2h = None
        flag_c_lambda2 = False
        flag_f_n_2h = False
        flag_f_m_perp_2h = False

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

        f_r_m, dder = calc_flip_ratio_by_structure_factors(
            polarization, flipper, f_nucl, f_m_perp, m_u, func_extinction=func_extinction,
            c_lambda2=c_lambda2, f_n_2h=f_n_2h, f_m_perp_2h=f_m_perp_2h,
            flag_polarization=flag_polarization, flag_flipper=flag_flipper,
            flag_f_n=flag_f_n, flag_f_m_perp=flag_f_m_perp,
            flag_c_lambda2=flag_c_lambda2,
            flag_f_n_2h=flag_f_n_2h, flag_f_m_perp_2h=flag_f_m_perp_2h)

        delta_f_r_m = (
            numpy.expand_dims(dder["f_m_perp_real"], 2)*
            delta_f_m_perp.real +
            numpy.expand_dims(dder["f_m_perp_imag"], 2)*
            delta_f_m_perp.imag).sum(axis=0)[:3]

        return f_r_m, delta_f_r_m


    def calc_chi_sq(
            self, l_crystal: List[Crystal], # or MagCrystal
            flag_internal=True, dict_in_out: dict = None):
        """
        Calculate chi square.

        Keyword Arguments
        -----------------
            - l_crystal: a list of Crystal objects of cryspy library
            - flag_internal: a flag to calculate internal objects
              (default is True)

        Output arguments
        ----------------
            - chi_sq_val: chi square of flip ratio
              (Sum_i ((y_e_i - y_m_i) / sigma_i)**2)
            - n: number of measured reflections
        """
        diffrn_refln = self.diffrn_refln
        index_h = diffrn_refln.numpy_index_h
        index_k = diffrn_refln.numpy_index_k
        index_l = diffrn_refln.numpy_index_l
        index_hkl = numpy.stack([index_h, index_k, index_l], axis=0)
        fr_exp = diffrn_refln.numpy_fr
        fr_sigma = diffrn_refln.numpy_fr_sigma
        fr_mod, dder = self.calc_iint_u_d_flip_ratio(index_hkl, l_crystal, flag_internal=flag_internal, dict_in_out=dict_in_out)

        if flag_internal:
            diffrn_refln.numpy_fr_calc = fr_mod
            # diffrn_refln.numpy_intensity_plus_calc = int_u_mod
            # diffrn_refln.numpy_intensity_minus_calc = int_d_mod
            diffrn_refln.numpy_to_items()

        flag_in = numpy.logical_not(diffrn_refln.numpy_excluded)

        chi_sq = ((fr_mod[flag_in]-fr_exp[flag_in])/fr_sigma[flag_in])**2
        chi_sq_val = (chi_sq[numpy.logical_not(numpy.isnan(chi_sq))]).sum()
        n = numpy.logical_not(numpy.isnan(chi_sq)).sum()

        if flag_internal:
            refine_ls = RefineLs(number_reflns=n,
                                 goodness_of_fit_all=chi_sq_val/float(n),
                                 weighting_scheme="sigma")
            self.refine_ls = refine_ls
        return chi_sq_val, n

    def params_to_cif(self, separator="_") -> str:
        """Save parameters in cif format."""
        ls_out = []
        l_cls = (Setup, DiffrnRadiation, DiffrnOrientMatrix, Extinction, Phase)
        l_obj = [item for item in self.items if type(item) in l_cls]
        ls_out.extend([_.to_cif(separator=separator)+"\n" for _ in l_obj])
        return "\n".join(ls_out)

    def data_to_cif(self, separator="_") -> str:
        """Save data in cif format."""
        ls_out = []
        l_cls = (DiffrnReflnL, )
        l_obj = [item for item in self.items if type(item) in l_cls]
        ls_out.extend([_.to_cif(separator=separator)+"\n" for _ in l_obj])
        return "\n".join(ls_out)

    def calc_to_cif(self, separator="_", flag=False) -> str:
        """Save calculations in cif format."""
        ls_out = []
        l_cls = (RefineLs, ReflnL, ReflnSusceptibilityL)
        l_obj = [item for item in self.items if type(item) in l_cls]
        ls_out.extend([_.to_cif(separator=separator)+"\n" for _ in l_obj])
        return "\n".join(ls_out)

    def estimate_FM(self, crystal, max_abs_f_mag: float = 5,
                    delta_f_mag: float = 0.001, flag_lambdaover2_correction=False):
        """
        Estimate magnetic structure factor.

        The method calculates magnetic structure factors (only real part) from
        flipping ratios based on give crystal structure.
        The crystal structure should be centrosymmetrical
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


            f_r, dder_f_r = calc_flip_ratio_by_structure_factors(
                       polarization, flipper_efficiency, f_nucl, f_m_perp, u_matrix, func_extinction=func_extinction,
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
            else:
                fig_ax = diffrn_refln.plot_fr_vs_fr_calc()

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
        
        self.form_object()
        ddict = {}
        setup, diffrn_refln, diffrn_orient_matrix = None, None, None
        extinction, diffrn_radiation = None, None
        phase, chi2 = None, None

        l_obj = take_items_by_class(self, (Setup, ))
        if len(l_obj) > 0:
            setup = l_obj[0]

        l_obj = take_items_by_class(self, (DiffrnReflnL, ))
        if len(l_obj) > 0:
            diffrn_refln = l_obj[0]

        l_obj = take_items_by_class(self, (DiffrnOrientMatrix, ))
        if len(l_obj) > 0:
            diffrn_orient_matrix = l_obj[0]

        l_obj = take_items_by_class(self, (Extinction, ))
        if len(l_obj) > 0:
            extinction = l_obj[0]

        l_obj = take_items_by_class(self, (DiffrnRadiation, ))
        if len(l_obj) > 0:
            diffrn_radiation = l_obj[0]

        l_obj = take_items_by_class(self, (Phase, ))
        if len(l_obj) > 0:
            phase = l_obj[0]

        l_obj = take_items_by_class(self, (Chi2, ))
        if len(l_obj) > 0:
            chi2 = l_obj[0]

        ddict["name"] = self.data_name
        ddict["type_name"] = self.get_name()
        if setup is not None:
            field = setup.field
            wavelength = setup.wavelength
            ddict["wavelength"] = numpy.array([wavelength], dtype=float)
            ddict["flags_wavelength"] = numpy.array([setup.wavelength_refinement], dtype=bool)
            if setup.is_attribute("ratio_lambdaover2"):
                ddict["c_lambda2"] = numpy.array([setup.ratio_lambdaover2], dtype=float)
                ddict["flags_c_lambda2"] = numpy.array([setup.ratio_lambdaover2_refinement], dtype=bool)
            if setup.is_attribute("temperature"):
                ddict["temperature"] = numpy.array([setup.temperature], dtype=float)
                
        if diffrn_refln is not None:
            index_hkl = numpy.array([
                diffrn_refln.index_h, diffrn_refln.index_k,
                diffrn_refln.index_l], dtype=float)
            refln_fr_es = numpy.array([
                diffrn_refln.fr, diffrn_refln.fr_sigma], dtype=float)
            refln_fr_excl = numpy.array(diffrn_refln.excluded, dtype=bool)

            ddict["index_hkl"] = index_hkl
            ddict["flip_ratio_es"] = refln_fr_es
            ddict["flip_ratio_excluded"] = refln_fr_excl

        if diffrn_orient_matrix is not None:
            u_matrix = diffrn_orient_matrix.u
            e_up = numpy.array([u_matrix[2,0], u_matrix[2,1], u_matrix[2,2]], dtype=float)
            
            ddict["matrix_u"] = numpy.array([
                diffrn_orient_matrix.u_11, diffrn_orient_matrix.u_12,
                diffrn_orient_matrix.u_13, diffrn_orient_matrix.u_21,
                diffrn_orient_matrix.u_22, diffrn_orient_matrix.u_23,
                diffrn_orient_matrix.u_31, diffrn_orient_matrix.u_32,
                diffrn_orient_matrix.u_33], dtype = float)
            ddict["matrix_ub"] = numpy.array([
                diffrn_orient_matrix.ub_11, diffrn_orient_matrix.ub_12,
                diffrn_orient_matrix.ub_13, diffrn_orient_matrix.ub_21,
                diffrn_orient_matrix.ub_22, diffrn_orient_matrix.ub_23,
                diffrn_orient_matrix.ub_31, diffrn_orient_matrix.ub_32,
                diffrn_orient_matrix.ub_33], dtype = float)

        if extinction is not None:
            model_extinction = extinction.model
            radius = extinction.radius
            mosaicity = extinction.mosaicity
            ddict["extinction_model"] = model_extinction
            ddict["extinction_radius"] = numpy.array([radius], dtype=float)
            ddict["extinction_mosaicity"] = numpy.array([mosaicity], dtype=float)
            ddict["flags_extinction_radius"] = numpy.array([extinction.radius_refinement], dtype=bool)
            ddict["flags_extinction_mosaicity"] = numpy.array([extinction.mosaicity_refinement], dtype=bool)

        if diffrn_radiation is not None:
            beam_polarization = diffrn_radiation.polarization
            flipper_efficiency = diffrn_radiation.efficiency
            ddict["beam_polarization"] = numpy.array([beam_polarization], dtype=float)
            ddict["flipper_efficiency"] = numpy.array([flipper_efficiency], dtype=float)
            ddict["flags_beam_polarization"] = numpy.array([diffrn_radiation.polarization_refinement], dtype=bool)
            ddict["flags_flipper_efficiency"] = numpy.array([diffrn_radiation.efficiency_refinement], dtype=bool)

        if phase is not None:
            ddict["phase_label"] = numpy.array([phase.label], dtype=str)
        else:
            ddict["phase_label"] = numpy.array([""], dtype=str)

        if chi2 is not None:
            if chi2.is_attribute("asymmetry"):
                ddict["flag_asymmetry"] = numpy.array([chi2.asymmetry], dtype=bool)
            else:
                ddict["flag_asymmetry"] = numpy.array([False], dtype=bool)
        else:
            ddict["flag_asymmetry"] = numpy.array([False], dtype=bool)
            
        keys = ddict.keys()
        if ((setup is not None) and ("matrix_u" in keys)):
            ddict["magnetic_field"] = field*ddict["matrix_u"][6:]
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

        if (("iint_plus" in keys) and ("iint_minus" in keys)):
            diffrn_refln = self.diffrn_refln
            diffrn_refln.numpy_fr_calc = ddict_diffrn["iint_plus"]/ddict_diffrn["iint_minus"]
            diffrn_refln.numpy_to_items()

