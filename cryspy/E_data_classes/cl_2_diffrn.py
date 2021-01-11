"""Description of Diffrn class."""
__author__ = 'ikibalin'
__version__ = "2020_08_19"

from warnings import warn
from typing import NoReturn, List
import numpy
import scipy
import scipy.misc

from cryspy.A_functions_base.function_1_matrices import scalar_product
from cryspy.A_functions_base.function_1_roots import calc_roots

from cryspy.B_parent_classes.cl_2_loop import LoopN
from cryspy.B_parent_classes.cl_3_data import DataN

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

from cryspy.D_functions_item_loop.function_1_flip_ratio import \
    calc_flip_ratio, calc_fm_perp_loc, calc_e_up_loc

from cryspy.E_data_classes.cl_1_crystal import Crystal


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
                        RefineLs)
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

    def calc_iint_u_d_flip_ratio(self, index_h, index_k, index_l, l_crystal,
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
            flip_ratio: flip ratio (intensity_up/intensity_down)

        """
        l_label = [_.data_name for _ in l_crystal]
        if self.phase.label in l_label:
            ind = l_label.index(self.phase.label)
            crystal = l_crystal[ind]
        elif self.data_name in l_label:
            self.phase = Phase(label=self.data_name)
            ind = l_label.index(self.data_name)
            crystal = l_crystal[ind]
        else:
            warn("Crystal not found. The first one is taken.")
            crystal = l_crystal[0]
            self.phase.label = crystal.label

        cell = crystal.cell
        volume_unit_cell = cell.volume
        sthovl = cell.calc_sthovl(index_h, index_k, index_l)
        k_loc_i = cell.calc_k_loc(index_h, index_k, index_l)

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
        e_up_loc = calc_e_up_loc(0., 0., 0., u_matrix)

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
            index_2h = 2 * index_h
            index_2k = 2 * index_k
            index_2l = 2 * index_l

        if flag_internal:
            refln = crystal.calc_refln(index_h, index_k, index_l)
            f_nucl = numpy.array(refln.f_calc, dtype=complex)
        else:
            f_nucl = crystal.calc_f_nucl(index_h, index_k, index_l)

        if flag_lambdaover2:
            f_nucl_2hkl = crystal.calc_f_nucl(index_2h, index_2k, index_2l)
        else:
            f_nucl_2hkl = None

        if flag_internal:
            refln_s = crystal.calc_refln_susceptibility(index_h, index_k,
                                                        index_l)
            sft_ij = (refln_s.numpy_chi_11_calc,
                      refln_s.numpy_chi_12_calc,
                      refln_s.numpy_chi_13_calc,
                      refln_s.numpy_chi_21_calc,
                      refln_s.numpy_chi_22_calc,
                      refln_s.numpy_chi_23_calc,
                      refln_s.numpy_chi_31_calc,
                      refln_s.numpy_chi_32_calc,
                      refln_s.numpy_chi_33_calc)

            sftm_ij = (refln_s.numpy_moment_11_calc,
                       refln_s.numpy_moment_12_calc,
                       refln_s.numpy_moment_13_calc,
                       refln_s.numpy_moment_21_calc,
                       refln_s.numpy_moment_22_calc,
                       refln_s.numpy_moment_23_calc,
                       refln_s.numpy_moment_31_calc,
                       refln_s.numpy_moment_32_calc,
                       refln_s.numpy_moment_33_calc)
        else:
            chi_m = crystal.calc_susceptibility_moment_tensor(
                index_h, index_k, index_l)
            sft_ij = chi_m[:9]
            sftm_ij = chi_m[9:]

        fm_perp_loc = calc_fm_perp_loc(e_up_loc, field_norm, k_loc_i, sft_ij,
                                       sftm_ij)

        if flag_lambdaover2:
            chi_m_2hkl = crystal.calc_susceptibility_moment_tensor(
                index_2h, index_2k, index_2l)
            sft_ij_2hkl = chi_m_2hkl[:9]
            sftm_ij_2hkl = chi_m_2hkl[9:]
            fm_perp_loc_2hkl = calc_fm_perp_loc(e_up_loc, field_norm, k_loc_i,
                                                sft_ij_2hkl, sftm_ij_2hkl)
        else:
            fm_perp_loc_2hkl = None

        flip_ratio, dder = calc_flip_ratio(
            sthovl, wavelength, field_norm, u_matrix, polarization,
            flipper_efficiency, volume_unit_cell, model_extinction, radius,
            mosaicity, f_nucl, fm_perp_loc,
            ratio_lambdaover2=ratio_lambdaover2, f_nucl_2hkl=f_nucl_2hkl,
            fm_perp_loc_2hkl=fm_perp_loc_2hkl)

        if flag_internal:
            refln.loop_name = crystal.data_name
            refln_s.loop_name = crystal.data_name
            self.add_items([refln, refln_s])

        return flip_ratio, dder

    # FIXME: merge with calc_iint_u_d_flip_ratio
    def calc_fr(self, cell, f_nucl, f_m_perp, delta_f_nucl=None,
                delta_f_m_perp=None):
        """
        Calculate flip ratio.

        By given nuclear structure factor and
        perpendicular component of magnetic structure factor
        (in 10**-12 cm)
        cell (only for extinction correction)
        f_nucl [hkl]
        f_m_perp [hkl]
        delta_f_nucl [hkl, parameters]
        delta_f_m_perp [hkl, parameters]
        """
        f_n_sq = (f_nucl * f_nucl.conjugate()).real
        f_m_p_x, f_m_p_y, f_m_p_z = f_m_perp
        delta_f_m_p_x, delta_f_m_p_y, delta_f_m_p_z = delta_f_m_perp

        f_m_p_sq = (f_m_p_x*f_m_p_x.conjugate() +
                    f_m_p_y*f_m_p_y.conjugate() +
                    f_m_p_z*f_m_p_z.conjugate()).real

        diffrn_radiation = self.diffrn_radiation
        diffrn_orient_matrix = self.diffrn_orient_matrix
        e_up = diffrn_orient_matrix.calc_e_up()
        p_u = float(diffrn_radiation.polarization)
        p_d = (2. * float(diffrn_radiation.efficiency) - 1.) * p_u

        f_m_p_vert = scalar_product(f_m_perp, e_up)

        if delta_f_m_perp is not None:
            delta_f_m_p_vert = scalar_product(delta_f_m_perp, e_up)

        two_f_nucl_f_m_p_vert = (f_nucl * f_m_p_vert.conjugate() +
                                 f_m_p_vert * f_nucl.conjugate()).real

        delta_two_f_nucl_f_m_p_vert = (
            f_nucl[:, numpy.newaxis]*delta_f_m_p_vert.conjugate() +
            delta_f_m_p_vert*f_nucl.conjugate()[:, numpy.newaxis]).real

        extinction = self.extinction
        diffrn_refln = self.diffrn_refln
        index_h = diffrn_refln.numpy_index_h
        index_k = diffrn_refln.numpy_index_k
        index_l = diffrn_refln.numpy_index_l

        setup = self.setup
        wavelength = setup.wavelength

        delta_f_m_p_sq = 2.*(
         numpy.real(f_m_p_x)[:, numpy.newaxis] * numpy.real(delta_f_m_p_x) +
         numpy.imag(f_m_p_x)[:, numpy.newaxis] * numpy.imag(delta_f_m_p_x) +
         numpy.real(f_m_p_y)[:, numpy.newaxis] * numpy.real(delta_f_m_p_y) +
         numpy.imag(f_m_p_y)[:, numpy.newaxis] * numpy.imag(delta_f_m_p_y) +
         numpy.real(f_m_p_z)[:, numpy.newaxis] * numpy.real(delta_f_m_p_z) +
         numpy.imag(f_m_p_z)[:, numpy.newaxis] * numpy.imag(delta_f_m_p_z))

        # extinction correction
        f_m_p_vert_sq = (f_m_p_vert * f_m_p_vert.conjugate()).real
        delta_f_m_p_vert_sq = 2.*(
            numpy.real(f_m_p_vert)[:, numpy.newaxis] *
            numpy.real(delta_f_m_p_vert) +
            numpy.imag(f_m_p_vert)[:, numpy.newaxis] *
            numpy.imag(delta_f_m_p_vert))

        fnp = two_f_nucl_f_m_p_vert
        fp_sq = f_n_sq + f_m_p_sq + two_f_nucl_f_m_p_vert
        fm_sq = f_n_sq + f_m_p_sq - two_f_nucl_f_m_p_vert
        fpm_sq = f_m_p_sq - f_m_p_vert_sq

        delta_fnp = delta_two_f_nucl_f_m_p_vert
        delta_fp_sq = delta_f_m_p_sq + delta_two_f_nucl_f_m_p_vert
        delta_fm_sq = delta_f_m_p_sq - delta_two_f_nucl_f_m_p_vert
        delta_fpm_sq = delta_f_m_p_sq - delta_f_m_p_vert_sq

        if extinction is not None:
            yp, dderp = extinction.calc_extinction(
                cell, index_h, index_k, index_l,
                fp_sq, wavelength, flag_derivatives=True)
            delta_yp = dderp["der_yext__f_sq"]
            ym, dderm = extinction.calc_extinction(
                cell, index_h, index_k, index_l,
                fm_sq, wavelength, flag_derivatives=True)
            delta_ym = dderm["der_yext__f_sq"]
            ypm, dderpm = extinction.calc_extinction(
                cell, index_h, index_k, index_l,
                fpm_sq, wavelength, flag_derivatives=True)
            delta_ypm = dderpm["der_yext__f_sq"]
        else:
            yp = numpy.ones(shape=fp_sq.shape, dtype=float),
            delta_yp = numpy.zeros(shape=fp_sq.shape, dtype=float)
            ym = numpy.ones(shape=fp_sq.shape, dtype=float)
            delta_ym = numpy.zeros(shape=fp_sq.shape, dtype=float)
            ypm = numpy.ones(shape=fp_sq.shape, dtype=float)
            delta_ypm = numpy.zeros(shape=fp_sq.shape, dtype=float)

        pppl = 0.5 * ((1 + p_u) * yp + (1 - p_u) * ym)
        ppmin = 0.5 * ((1 - p_d) * yp + (1 + p_d) * ym)
        pmpl = 0.5 * ((1 + p_u) * yp - (1 - p_u) * ym)
        pmmin = 0.5 * ((1 - p_d) * yp - (1 + p_d) * ym)

        delta_pppl = 0.5 * ((1 + p_u) * delta_yp[:, numpy.newaxis] *
                            delta_fp_sq + (1 - p_u) *
                            delta_ym[:, numpy.newaxis] * delta_fm_sq)
        delta_ppmin = 0.5 * ((1 - p_d) * delta_yp[:, numpy.newaxis] *
                             delta_fp_sq + (1 + p_d) *
                             delta_ym[:, numpy.newaxis] * delta_fm_sq)
        delta_pmpl = 0.5 * ((1 + p_u) * delta_yp[:, numpy.newaxis] *
                            delta_fp_sq - (1 - p_u) *
                            delta_ym[:, numpy.newaxis] * delta_fm_sq)
        delta_pmmin = 0.5 * ((1 - p_d) * delta_yp[:, numpy.newaxis] *
                             delta_fp_sq - (1 + p_d) *
                             delta_ym[:, numpy.newaxis] * delta_fm_sq)

        iint_u = (f_n_sq+f_m_p_vert_sq)*pppl+pmpl*fnp+ypm*fpm_sq
        iint_d = (f_n_sq+f_m_p_vert_sq)*ppmin+pmmin*fnp+ypm*fpm_sq

        delta_iint_u = (((f_n_sq + f_m_p_vert_sq)[:, numpy.newaxis] *
                         delta_pppl + delta_pmpl * fnp[:, numpy.newaxis] +
                         (delta_ypm * fpm_sq)[:, numpy.newaxis] * delta_fpm_sq)
                        + (delta_f_m_p_vert_sq) * pppl[:, numpy.newaxis] +
                        pmpl[:, numpy.newaxis] * delta_fnp +
                        ypm[:, numpy.newaxis] * delta_fpm_sq)

        delta_iint_d = (((f_n_sq + f_m_p_vert_sq)[:, numpy.newaxis] *
                         delta_ppmin + delta_pmmin * fnp[:, numpy.newaxis] +
                         (delta_ypm * fpm_sq)[:, numpy.newaxis] * delta_fpm_sq)
                        + (delta_f_m_p_vert_sq) * ppmin[:, numpy.newaxis] +
                        pmmin[:, numpy.newaxis] * delta_fnp +
                        ypm[:, numpy.newaxis] * delta_fpm_sq)

        f_r_m = iint_u/iint_d

        delta_f_r_m = (iint_d[:, numpy.newaxis] * delta_iint_u -
                       iint_u[:, numpy.newaxis] * delta_iint_d) / \
            (iint_d**2)[:, numpy.newaxis]

        return f_r_m, delta_f_r_m

    def calc_chi_sq(self, l_crystal: List[Crystal], flag_internal=True):
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
        np_h = diffrn_refln.numpy_index_h
        np_k = diffrn_refln.numpy_index_k
        np_l = diffrn_refln.numpy_index_l
        fr_exp = diffrn_refln.numpy_fr
        fr_sigma = diffrn_refln.numpy_fr_sigma
        fr_mod, dder = self.calc_iint_u_d_flip_ratio(
            np_h, np_k, np_l, l_crystal, flag_internal=flag_internal)

        if flag_internal:
            diffrn_refln.numpy_fr_calc = fr_mod
            # diffrn_refln.numpy_intensity_up_calc = int_u_mod
            # diffrn_refln.numpy_intensity_down_calc = int_d_mod
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

    def estimate_FM(self, crystal, maxAbsFM: float = 5,
                    deltaFM: float = 0.001, flag_lambdaover2_correction=False):
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
        e_up_loc = calc_e_up_loc(0., 0., 0., u_matrix)

        refln = crystal.calc_refln(index_h, index_k, index_l)
        f_nucl = numpy.array(refln.f_calc, dtype=complex)

        cell = crystal.cell
        volume_unit_cell = cell.volume
        sthovl = cell.calc_sthovl(index_h, index_k, index_l)
        k_loc_i = cell.calc_k_loc(index_h, index_k, index_l)
        k_loc_1, k_loc_2, k_loc_3 = k_loc_i

        if (setup.is_attribute("ratio_lambdaover2") &
                flag_lambdaover2_correction):
            ratio_lambdaover2 = setup.ratio_lambdaover2
            flag_lambdaover2 = True
        else:
            ratio_lambdaover2 = None
            flag_lambdaover2 = False

        if (flag_lambdaover2 & flag_lambdaover2_correction):
            index_2h = 2 * index_h
            index_2k = 2 * index_k
            index_2l = 2 * index_l
            f_nucl_2hkl = crystal.calc_f_nucl(index_2h, index_2k, index_2l)
        else:
            f_nucl_2hkl = None

        if (flag_lambdaover2 & flag_lambdaover2_correction):
            chi_m_2hkl = crystal.calc_susceptibility_moment_tensor(
                index_2h, index_2k, index_2l)
            sft_ij_2hkl = chi_m_2hkl[:9]
            sftm_ij_2hkl = chi_m_2hkl[9:]
            fm_perp_loc_2hkl = calc_fm_perp_loc(e_up_loc, field_norm, k_loc_i,
                                                sft_ij_2hkl, sftm_ij_2hkl)
        else:
            fm_perp_loc_2hkl = None


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
            mag_loc_1 = e_up_loc[0]*f_mag*0.2695
            mag_loc_2 = e_up_loc[1]*f_mag*0.2695
            mag_loc_3 = e_up_loc[2]*f_mag*0.2695

            # vector product k x (0, 0, mag_3) x k
            mag_p_1 = (k_loc_3*mag_loc_1 - k_loc_1*mag_loc_3)*k_loc_3 - \
                (k_loc_1*mag_loc_2 - k_loc_2*mag_loc_1)*k_loc_2
            mag_p_2 = (k_loc_1*mag_loc_2 - k_loc_2*mag_loc_1)*k_loc_1 - \
                (k_loc_2*mag_loc_3 - k_loc_3*mag_loc_2)*k_loc_3
            mag_p_3 = (k_loc_2*mag_loc_3 - k_loc_3*mag_loc_2)*k_loc_2 - \
                (k_loc_3*mag_loc_1 - k_loc_1*mag_loc_3)*k_loc_1

            fm_perp_loc = (mag_p_1, mag_p_2, mag_p_3)
            flip_ratio, dder = calc_flip_ratio(
                sthovl, wavelength, field_norm, u_matrix, polarization,
                flipper_efficiency, volume_unit_cell, model_extinction, radius,
                mosaicity, f_nucl, fm_perp_loc,
                ratio_lambdaover2=ratio_lambdaover2,
                f_nucl_2hkl=f_nucl_2hkl,
                fm_perp_loc_2hkl=fm_perp_loc_2hkl)
            return flip_ratio

        lFMans, lFMansMin = [], []
        lFMansMin_sigma = []
        for _ind in range(fr_exp.size):
            def f(x):
                return (function_to_find_fm(x)[_ind] - fr_exp[_ind])
            FMans = calc_roots(f, -1*maxAbsFM-deltaFM, maxAbsFM+deltaFM)
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
            fr_ax = diffrn_refln.plot_fr_vs_fr_calc()
            if fr_ax is not None:
                fr, ax = fr_ax
                ax.set_title(self.data_name + " - "+ax.title.get_text())
                l_res.append((fr, ax))        
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

# s_cont = """
#  data_mono
#  _setup_wavelength     0.840
#  _setup_field          1.000

#  _diffrn_radiation_polarization 1.0
#  _diffrn_radiation_efficiency   1.0

#  _extinction_mosaicity 100.0
#  _extinction_radius    50.0
#  _extinction_model     gauss

#  _diffrn_orient_matrix_UB_11 6.59783
#  _diffrn_orient_matrix_UB_12 -6.99807
#  _diffrn_orient_matrix_UB_13 3.3663
#  _diffrn_orient_matrix_UB_21 2.18396
#  _diffrn_orient_matrix_UB_22 -2.60871
#  _diffrn_orient_matrix_UB_23 -9.5302
#  _diffrn_orient_matrix_UB_31 7.4657
#  _diffrn_orient_matrix_UB_32 6.94702
#  _diffrn_orient_matrix_UB_33 -0.18685

#  _phase_label  Fe3O4

#  loop_
#  _diffrn_refln_index_h
#  _diffrn_refln_index_k
#  _diffrn_refln_index_l
#  _diffrn_refln_fr
#  _diffrn_refln_fr_sigma
#  0 0 8 0.64545 0.01329
#  2 0 6 1.75682 0.04540
#  0 2 6 1.67974 0.03711

# """

# obj = Diffrn.from_cif(s_cont)
# print(obj)
# for var_name in obj.get_variable_names():
#     print(var_name)
