__author__ = 'ikibalin'
__version__ = "2020_01_02"
import os
import numpy
from pycifstar import Global

import warnings
from typing import List, Tuple
from cryspy.common.cl_item_constr import ItemConstr
from cryspy.common.cl_loop_constr import LoopConstr
from cryspy.common.cl_data_constr import DataConstr
from cryspy.common.cl_fitable import Fitable

from cryspy.corecif.cl_refln import Refln, ReflnL
from cryspy.corecif.cl_diffrn_orient_matrix import DiffrnOrientMatrix
from cryspy.corecif.cl_refine_ls import RefineLs

from cryspy.magneticcif.cl_refln_susceptibility import ReflnSusceptibility, ReflnSusceptibilityL

from .cl_diffrn_refln import DiffrnRefln, DiffrnReflnL
from .cl_diffrn_radiation import DiffrnRadiation
from .cl_extinction import Extinction
from .cl_setup import Setup
from .cl_phase import Phase




class Diffrn(DataConstr):
    """
Class to describe information about single diffraction measurements

Description in cif file::

 data_mono
 _setup_wavelength     0.840
 _setup_field          1.000
 
 _diffrn_radiation_polarization 1.0
 _diffrn_radiation_efficiency   1.0

 _extinction_mosaicity 100.0
 _extinction_radius    50.0
 _extinction_model     gauss

 _diffrn_orient_matrix_UB_11 6.59783
 _diffrn_orient_matrix_UB_12 -6.99807
 _diffrn_orient_matrix_UB_13 3.3663
 _diffrn_orient_matrix_UB_21 2.18396
 _diffrn_orient_matrix_UB_22 -2.60871
 _diffrn_orient_matrix_UB_23 -9.5302
 _diffrn_orient_matrix_UB_31 7.4657
 _diffrn_orient_matrix_UB_32 6.94702
 _diffrn_orient_matrix_UB_33 -0.18685

 _phase_label  Fe3O4

 loop_
 _diffrn_refln_index_h
 _diffrn_refln_index_k
 _diffrn_refln_index_l
 _diffrn_refln_fr
 _diffrn_refln_fr_sigma
 0 0 8 0.64545 0.01329 
 2 0 6 1.75682 0.04540 
 0 2 6 1.67974 0.03711 
    """
    MANDATORY_CLASSES = (Setup, DiffrnRadiation, DiffrnOrientMatrix, DiffrnReflnL)
    OPTIONAL_CLASSES = (Extinction, Phase)
    INTERNAL_CLASSES = (RefineLs, ReflnL, ReflnSusceptibilityL)
    def __init__(self, setup=None, diffrn_radiation=None, diffrn_orient_matrix=None, 
                 diffrn_refln=None, extinction=None, phase=None,
                 data_name=""):

        super(Diffrn, self).__init__(mandatory_classes=self.MANDATORY_CLASSES,
                                     optional_classes=self.OPTIONAL_CLASSES,
                                     internal_classes=self.INTERNAL_CLASSES)

        self.data_name = data_name
        self.setup = setup
        self.diffrn_radiation = diffrn_radiation
        self.diffrn_orient_matrix = diffrn_orient_matrix
        self.diffrn_refln = diffrn_refln
        self.extinction = extinction
        self.phase = phase

        #FIXME: internal attributes:
        #self.reflns = None
        #self.refln_ss = None

        if self.is_defined:
            self.form_object
        

    @property
    def setup(self):
        l_res = self[Setup]
        if len(l_res) >= 1:
            return l_res[0]
        else:
            return None
    @setup.setter
    def setup(self, x):
        if x is None:
            pass
        elif isinstance(x, Setup):
            l_ind = []
            for _i, _obj in enumerate(self.mandatory_objs):
                if isinstance(_obj, Setup):
                    l_ind.append(_i)
            if len(l_ind) == 0:
                self.mandatory_objs.append(x)
            else:
                self.mandatory_objs[l_ind[0]] = x
            if len(l_ind) > 1:
                for _ind in l_ind.reverse():
                    self.mandatory_objs.pop(_ind)

    @property
    def diffrn_radiation(self):
        l_res = self[DiffrnRadiation]
        if len(l_res) >= 1:
            return l_res[0]
        else:
            return None
    @diffrn_radiation.setter
    def diffrn_radiation(self, x):
        if x is None:
            pass
        elif isinstance(x, DiffrnRadiation):
            l_ind = []
            for _i, _obj in enumerate(self.mandatory_objs):
                if isinstance(_obj, DiffrnRadiation):
                    l_ind.append(_i)
            if len(l_ind) == 0:
                self.mandatory_objs.append(x)
            else:
                self.mandatory_objs[l_ind[0]] = x
            if len(l_ind) > 1:
                for _ind in l_ind.reverse():
                    self.mandatory_objs.pop(_ind)

    @property
    def diffrn_orient_matrix(self):
        l_res = self[DiffrnOrientMatrix]
        if len(l_res) >= 1:
            return l_res[0]
        else:
            return None
    @diffrn_orient_matrix.setter
    def diffrn_orient_matrix(self, x):
        if x is None:
            pass
        elif isinstance(x, DiffrnOrientMatrix):
            l_ind = []
            for _i, _obj in enumerate(self.mandatory_objs):
                if isinstance(_obj, DiffrnOrientMatrix):
                    l_ind.append(_i)
            if len(l_ind) == 0:
                self.mandatory_objs.append(x)
            else:
                self.mandatory_objs[l_ind[0]] = x
            if len(l_ind) > 1:
                for _ind in l_ind.reverse():
                    self.mandatory_objs.pop(_ind)

    @property
    def diffrn_refln(self):
        l_res = self[DiffrnReflnL]
        if len(l_res) >= 1:
            return l_res[0]
        else:
            return None
    @diffrn_refln.setter
    def diffrn_refln(self, x):
        if x is None:
            pass
        elif isinstance(x, DiffrnReflnL):
            l_ind = []
            for _i, _obj in enumerate(self.mandatory_objs):
                if isinstance(_obj, DiffrnReflnL):
                    l_ind.append(_i)
            if len(l_ind) == 0:
                self.mandatory_objs.append(x)
            else:
                self.mandatory_objs[l_ind[0]] = x
            if len(l_ind) > 1:
                for _ind in l_ind.reverse():
                    self.mandatory_objs.pop(_ind)


    @property
    def extinction(self):
        l_res = self[Extinction]
        if len(l_res) >= 1:
            return l_res[0]
        else:
            return None
    @extinction.setter
    def extinction(self, x):
        if x is None:
            pass
        elif isinstance(x, Extinction):
            l_ind = []
            for _i, _obj in enumerate(self.optional_objs):
                if isinstance(_obj, Extinction):
                    l_ind.append(_i)
            if len(l_ind) == 0:
                self.optional_objs.append(x)
            else:
                self.optional_objs[l_ind[0]] = x
            if len(l_ind) > 1:
                for _ind in l_ind.reverse():
                    self.optional_objs.pop(_ind)


    @property
    def phase(self):
        l_res = self[Phase]
        if len(l_res) >= 1:
            return l_res[0]
        else:
            return None
    @phase.setter
    def phase(self, x):
        if x is None:
            pass
        elif isinstance(x, Phase):
            l_ind = []
            for _i, _obj in enumerate(self.optional_objs):
                if isinstance(_obj, Phase):
                    l_ind.append(_i)
            if len(l_ind) == 0:
                self.optional_objs.append(x)
            else:
                self.optional_objs[l_ind[0]] = x
            if len(l_ind) > 1:
                for _ind in l_ind.reverse():
                    self.optional_objs.pop(_ind)


    @property
    def refine_ls(self):
        l_res = self[RefineLs]
        if len(l_res) >= 1:
            return l_res[0]
        else:
            return None

    @property
    def refln(self):
        l_res = self[ReflnL]
        return l_res

    @property
    def refln_susceptibility(self):
        l_res = self[ReflnSusceptibilityL]
        return l_res




    def __repr__(self):
        ls_out = ["Diffrn:"]
        ls_out.append(f"{str(self):}")
        return "\n".join(ls_out)

    def _show_message(self, s_out: str):
        warnings.warn("***  Error ***", s_out, UserWarning, stacklevel=2)

    @property
    def is_variable(self) -> bool:
        """
Output: True if there is any refined parameter
        """
        res = any([_obj.is_variable for _obj in self.mandatory_objs] +
                  [_obj.is_variable for _obj in self.optional_objs])
        return res
        
    def get_variables(self) -> List:
        """
Output: the list of the refined parameters
        """
        l_variable = []
        for _obj in self.mandatory_objs:
            l_variable.extend(_obj.get_variables())
        for _obj in self.optional_objs:
            l_variable.extend(_obj.get_variables())
        return l_variable



    def calc_iint_u_d_flip_ratio(self, h, k, l, l_crystal, flag_internal=True):
        """
Calculate intensity for the given diffraction angle

Keyword arguments:

    h, k, l: 1D numpy array of Miller indexes

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
            self._show_message("Crystal not found.\nThe first one is taken.")
            crystal = l_crystal[0]
            self.phase.label = crystal.label
        wavelength = float(self.setup.wavelength)
        field_z = float(self.setup.field)
        diffrn_radiation = self.diffrn_radiation
        extinction = self.extinction
        orient_matrix = self.diffrn_orient_matrix

        orientation = orient_matrix.u       
        field_vec = numpy.array([0., 0., field_z], dtype=float)

        phi_d, chi_d, omega_d = 0., 0., 0.
        m_phi_d = numpy.array([[ numpy.cos(phi_d), numpy.sin(phi_d), 0.],
                               [-numpy.sin(phi_d), numpy.cos(phi_d), 0.],
                               [               0.,               0., 1.]], dtype=float)

        m_omega_d = numpy.array([[ numpy.cos(omega_d), numpy.sin(omega_d), 0.],
                                 [-numpy.sin(omega_d), numpy.cos(omega_d), 0.],
                                 [                 0.,                 0., 1.]], dtype=float)

        m_chi_d = numpy.array([[ numpy.cos(chi_d), 0., numpy.sin(chi_d)],
                               [               0., 1.,               0.],
                               [-numpy.sin(chi_d), 0., numpy.cos(chi_d)]], dtype=float)
        
        m_u_d = numpy.matmul(m_omega_d, numpy.matmul(m_chi_d, numpy.matmul(m_phi_d, 
                             orientation)))

        field_loc = numpy.matmul(m_u_d.transpose(), field_vec)
        
        
        field_norm = ((field_loc**2).sum())**0.5
        
        p_u = float(diffrn_radiation.polarization)
        p_d = (2.*float(diffrn_radiation.efficiency)-1.)*p_u
        
        e_u_loc = field_loc/field_norm
        
        if flag_internal:
            refln = crystal.calc_refln(h, k, l)
            f_nucl  =  numpy.array(refln.f_calc, dtype=complex)
        else:
            f_nucl = crystal.calc_f_nucl(h, k, l) 
        

        if flag_internal:
            refln_s = crystal.calc_refln_susceptibility(h, k, l)
            sft_11, sft_12, sft_13 = numpy.array(refln_s.chi_11_calc, dtype=complex), numpy.array(refln_s.chi_12_calc, dtype=complex), numpy.array(refln_s.chi_13_calc, dtype=complex)
            sft_21, sft_22, sft_23 = numpy.array(refln_s.chi_21_calc, dtype=complex), numpy.array(refln_s.chi_22_calc, dtype=complex), numpy.array(refln_s.chi_23_calc, dtype=complex)
            sft_31, sft_32, sft_33 = numpy.array(refln_s.chi_31_calc, dtype=complex), numpy.array(refln_s.chi_32_calc, dtype=complex), numpy.array(refln_s.chi_33_calc, dtype=complex)
            sftm_11, sftm_12, sftm_13 = numpy.array(refln_s.moment_11_calc, dtype=complex), numpy.array(refln_s.moment_12_calc, dtype=complex), numpy.array(refln_s.moment_13_calc, dtype=complex)
            sftm_21, sftm_22, sftm_23 = numpy.array(refln_s.moment_21_calc, dtype=complex), numpy.array(refln_s.moment_22_calc, dtype=complex), numpy.array(refln_s.moment_23_calc, dtype=complex)
            sftm_31, sftm_32, sftm_33 = numpy.array(refln_s.moment_31_calc, dtype=complex), numpy.array(refln_s.moment_32_calc, dtype=complex), numpy.array(refln_s.moment_33_calc, dtype=complex)
        else:
            chi_m = crystal.calc_susceptibility_moment_tensor(h, k, l) 
            sft_11, sft_12, sft_13, sft_21, sft_22, sft_23, sft_31, sft_32, sft_33 = chi_m[:9]
            sftm_11, sftm_12, sftm_13, sftm_21, sftm_22, sftm_23, sftm_31, sftm_32, sftm_33 = chi_m[9:] 

        cell = crystal.cell
        k_1, k_2, k_3 = cell.calc_k_loc(h, k, l)
        
        
        #not sure about e_u_loc at field < 0 
        mag_1 = sft_11*field_loc[0] + sft_12*field_loc[1] + sft_13*field_loc[2] + sftm_11*e_u_loc[0] + sftm_12*e_u_loc[1] + sftm_13*e_u_loc[2]
        mag_2 = sft_21*field_loc[0] + sft_22*field_loc[1] + sft_23*field_loc[2] + sftm_21*e_u_loc[0] + sftm_22*e_u_loc[1] + sftm_23*e_u_loc[2]
        mag_3 = sft_31*field_loc[0] + sft_32*field_loc[1] + sft_33*field_loc[2] + sftm_31*e_u_loc[0] + sftm_32*e_u_loc[1] + sftm_33*e_u_loc[2]
        
        #vector product k x mag x k        
        mag_p_1 = (k_3*mag_1 - k_1*mag_3)*k_3 - (k_1*mag_2 - k_2*mag_1)*k_2
        mag_p_2 = (k_1*mag_2 - k_2*mag_1)*k_1 - (k_2*mag_3 - k_3*mag_2)*k_3
        mag_p_3 = (k_2*mag_3 - k_3*mag_2)*k_2 - (k_3*mag_1 - k_1*mag_3)*k_1
        
        mag_p_sq = abs(mag_p_1*mag_p_1.conjugate() + mag_p_2*mag_p_2.conjugate() + 
                       mag_p_3*mag_p_3.conjugate())
        
        mag_p_e_u = mag_p_1*e_u_loc[0]+mag_p_2*e_u_loc[1]+mag_p_3*e_u_loc[2]
        
        f_nucl_sq = abs(f_nucl)**2
        mag_p_e_u_sq = abs(mag_p_e_u*mag_p_e_u.conjugate())
        fnp = (mag_p_e_u*f_nucl.conjugate()+mag_p_e_u.conjugate()*f_nucl).real
        fp_sq = f_nucl_sq + mag_p_sq + fnp
        fm_sq = f_nucl_sq + mag_p_sq - fnp
        fpm_sq = mag_p_sq - mag_p_e_u_sq

        #extinction correction     
        yp = extinction.calc_extinction(cell, h, k, l, fp_sq, wavelength)
        ym = extinction.calc_extinction(cell, h, k, l, fm_sq, wavelength)
        ypm = extinction.calc_extinction(cell, h, k, l, fpm_sq, wavelength)

        pppl = 0.5*((1+p_u)*yp+(1-p_u)*ym)
        ppmin= 0.5*((1-p_d)*yp+(1+p_d)*ym)
        pmpl = 0.5*((1+p_u)*yp-(1-p_u)*ym)
        pmmin= 0.5*((1-p_d)*yp-(1+p_d)*ym)
        """
        print("   h   k   l  f_nucl f_m_p_sq   f_np  fpm_sq")
        for h1, k1, l1, f_n_sq, f_m_sq, f_np, f_pm_sq in zip(h, k, l, f_nucl_sq, mag_p_e_u_sq, fnp, fpm_sq):
            print(" {:3} {:3} {:3} {:7.3f} {:7.3f} {:7.3f} {:7.3f}".format(
                    h1, k1, l1, f_n_sq, f_m_sq, f_np, f_pm_sq))

        print("   h   k   l      yp      ym     ypm      pppl   ppmin    pmpl   pmmin")
        for h1, k1, l1, y_p, y_m, y_pm, p_ppl, p_pmin, p_mpl, p_mmin in zip(
                h, k, l, yp, ym, ypm, pppl, ppmin, pmpl, pmmin):
            print(" {:3} {:3} {:3} {:7.3f} {:7.3f} {:7.3f}   {:7.3f} {:7.3f} {:7.3f} {:7.3f}".format(
                    h1, k1, l1, y_p, y_m, y_pm, p_ppl, p_pmin, p_mpl, p_mmin))
        """
        #integral intensities and flipping ratios
        iint_u = (f_nucl_sq+mag_p_e_u_sq)*pppl + pmpl*fnp + ypm*fpm_sq
        iint_d = (f_nucl_sq+mag_p_e_u_sq)*ppmin + pmmin*fnp + ypm*fpm_sq

        flip_ratio = iint_u/iint_d
        
        """
        d_info_out = {"iint_u": iint_u, "iint_d": iint_d, 
                      "flip_ratio": flip_ratio}
        d_info_out.update(d_info_sf)      
        print("   h   k   l  iint_u  iint_d flip_ratio")
        for h1, k1, l1, i_u, i_d, f_r in zip(h, k, l, iint_u, iint_d, flip_ratio):
            print("{:3} {:3} {:3} {:7.3f} {:7.3f} {:7.3f}".format(
                    h1, k1, l1, i_u, i_d, f_r))
        """
        
        if flag_internal:
            refln.loop_name = crystal.data_name
            refln_s.loop_name = crystal.data_name
            l_internal_objs = [refln, refln_s]
            setattr(self, "__internal_objs", l_internal_objs)

        return iint_u, iint_d, flip_ratio
    
    def calc_chi_sq(self, l_crystal: list, flag_internal=True):
        """
Calculate chi square

Keyword arguments:

    l_crystal: a list of Crystal objects of cryspy library
    flag_internal: a flag to calculate internal objects (default is True)

Output arguments:

    chi_sq_val: chi square of flip ratio (Sum_i ((y_e_i - y_m_i) / sigma_i)**2)
    n: number of measured reflections
        """
        diffrn_refln = self.diffrn_refln
        np_h = numpy.array(diffrn_refln.index_h, dtype=int)
        np_k = numpy.array(diffrn_refln.index_k, dtype=int)
        np_l = numpy.array(diffrn_refln.index_l, dtype=int)
        fr_exp = diffrn_refln.fr
        fr_sigma = diffrn_refln.fr_sigma
        int_u_mod, int_d_mod, fr_mod = self.calc_iint_u_d_flip_ratio(
                                               np_h, np_k, np_l, l_crystal, flag_internal=flag_internal)
        
        if flag_internal:
            for _fr_mod, _int_u_mod, _int_d_mod, _item in zip(fr_mod, int_u_mod, int_d_mod, diffrn_refln.item):
                _item.fr_calc = _fr_mod 
                _item._intensity_up_calc = _int_u_mod 
                _item._intensity_down_calc = _int_d_mod 
            for _obj in self.internal_objs:
                if isinstance(_obj, ReflnL):
                    for _item, _1, _2, _3 in zip(_obj.item, fr_mod, int_u_mod, int_d_mod):
                        _item.fr_calc = _1
                        _item.intensity_up_calc = _2
                        _item.intensity_down_calc = _3

        #self.refln = refln
        #self.refln_s = refln_s

        chi_sq = ((fr_mod-fr_exp)/fr_sigma)**2
        chi_sq_val = (chi_sq[numpy.logical_not(numpy.isnan(chi_sq))]).sum()
        n = numpy.logical_not(numpy.isnan(chi_sq)).sum()

        if flag_internal:
            refine_ls = RefineLs(number_reflns=n, goodness_of_fit_all=chi_sq_val/float(n), weighting_scheme="sigma")
            self.internal_objs.append(refine_ls)
        return chi_sq_val, n

    def params_to_cif(self, separator="_", flag=False) -> str: 
        ls_out = []
        l_cls = (Setup, DiffrnRadiation, DiffrnOrientMatrix, Extinction, Phase)
        l_obj = [_obj for _obj in (self.mandatory_objs + self.optional_objs) if type(_obj) in l_cls]
        ls_out.extend([_.to_cif(separator=separator, flag=flag)+"\n" for _ in l_obj])
        return "\n".join(ls_out)

    def data_to_cif(self, separator="_", flag=False) -> str: 
        ls_out = []
        l_cls = (DiffrnReflnL, )
        l_obj = [_obj for _obj in (self.mandatory_objs + self.optional_objs) if type(_obj) in l_cls]
        ls_out.extend([_.to_cif(separator=separator, flag=flag)+"\n" for _ in l_obj])
            
        return "\n".join(ls_out)

    def calc_to_cif(self, separator="_", flag=False) -> str: 
        ls_out = []
        l_cls = (RefineLs, ReflnL, ReflnSusceptibilityL)
        for _cls in l_cls:
            l_obj = [_obj for _obj in (self.internal_objs) if isinstance(_obj, _cls)]
            ls_out.extend([_.to_cif(separator=separator, flag=flag)+"\n" for _ in l_obj])

        #if self.reflns is not None:
        #    ls_out.extend([_.to_cif(separator=separator, flag=flag)+"\n" for _ in self.reflns])
        #if self.refln_ss is not None:
        #    ls_out.extend([_.to_cif(separator=separator, flag=flag)+"\n" for _ in self.refln_ss])
        return "\n".join(ls_out)
