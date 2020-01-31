__author__ = 'ikibalin'
__version__ = "2020_01_03"
import os
import math
import numpy
from pycifstar import Global

import warnings
from typing import List, Tuple
from cryspy.common.cl_item_constr import ItemConstr
from cryspy.common.cl_loop_constr import LoopConstr
from cryspy.common.cl_data_constr import DataConstr
from cryspy.common.cl_fitable import Fitable

from cryspy.corecif.cl_refln import Refln, ReflnL
from cryspy.corecif.cl_refine_ls import RefineLs
from cryspy.magneticcif.cl_refln_susceptibility import ReflnSusceptibility, ReflnSusceptibilityL

from cryspy.cif_like.cl_exclude import Exclude, ExcludeL
import cryspy.cif_like.CONSTANTS_AND_FUNCTIONS as CONSTANTS_AND_FUNCTIONS
from cryspy.cif_like.cl_phase import Phase, PhaseL
from cryspy.cif_like.cl_diffrn_radiation import DiffrnRadiation
from cryspy.cif_like.cl_extinction import Extinction
from cryspy.cif_like.cl_setup import Setup
from cryspy.cif_like.cl_range import Range
from cryspy.cif_like.cl_chi2 import Chi2
from cryspy.cif_like.cl_texture import Texture



from cryspy.pd2dcif_like.cl_pd2d_background import Pd2dBackground
from cryspy.pd2dcif_like.cl_pd2d_instr_reflex_asymmetry import Pd2dInstrReflexAsymmetry
from cryspy.pd2dcif_like.cl_pd2d_instr_resolution import Pd2dInstrResolution
from cryspy.pd2dcif_like.cl_pd2d_meas import Pd2dMeas 
from cryspy.pd2dcif_like.cl_pd2d_peak import Pd2dPeak, Pd2dPeakL
from cryspy.pd2dcif_like.cl_pd2d_proc import Pd2dProc


def calc_cos_ang(cell, h_1, k_1, l_1, h_2, k_2, l_2):
    q_1_x, q_1_y, q_1_z = cell.calc_k_loc(h_1,k_1,l_1)
    q_2_x, q_2_y, q_2_z = cell.calc_k_loc(h_2,k_2,l_2)
    q_1_sq = q_1_x*q_1_x + q_1_y*q_1_y + q_1_z*q_1_z
    q_2_sq = q_2_x*q_2_x + q_2_y*q_2_y + q_2_z*q_2_z
    q_12 = q_1_x*q_2_x + q_1_y*q_2_y + q_1_z*q_2_z
    res = q_12/(q_1_sq*q_2_sq)**0.5
    res[res>1.] = 1.
    return res


class Pd2d(DataConstr):
    """
Class to describe information about powder diffraction measurements

Description in cif file::

 data_powder2d
 _setup_wavelength     0.840
 _setup_field          1.000
 _setup_offset_2theta -0.385
 
 _chi2_sum  True
 _chi2_diff False
 _chi2_up   False
 _chi2_down False
 
 _range_2theta_min     4.000
 _range_2theta_max    80.000
 _range_phi_min       -2.000
 _range_phi_max       40.000
 
 _pd2d_background_intensity
 ;
 ;
 
 loop_
 _exclude_2theta_min
 _exclude_2theta_max
 _exclude_phi_min
 _exclude_phi_max
 0.0 1.0 0. 1.0
 
 _pd2d_instr_reflex_asymmetry_p1 0.0
 _pd2d_instr_reflex_asymmetry_p2 0.0
 _pd2d_instr_reflex_asymmetry_p3 0.0
 _pd2d_instr_reflex_asymmetry_p4 0.0
 
 _diffrn_radiation_polarization 0.0
 _diffrn_radiation_efficiency   1.0
 
 _pd2d_instr_resolution_u 16.9776
 _pd2d_instr_resolution_v -2.8357
 _pd2d_instr_resolution_w  0.5763
 _pd2d_instr_resolution_x  0.0
 _pd2d_instr_resolution_y  0.0
 
 loop_
 _phase_label
 _phase_scale
 _phase_igsize
 Fe3O4 0.02381 0.0
 
 _pd2d_meas_intensity_up
 ;
 ;

 _pd2d_meas_intensity_up_sigma
 ;
 ;

 _pd2d_meas_intensity_down
 ;
 ;

 _pd2d_meas_intensity_down_sigma
 ;
 ;
    """
    MANDATORY_CLASSES = (Pd2dBackground, Pd2dInstrResolution, Pd2dMeas, PhaseL, DiffrnRadiation,
                         Setup, Range, Chi2)
    OPTIONAL_CLASSES = (Extinction, ExcludeL, Pd2dInstrReflexAsymmetry, Texture)
    INTERNAL_CLASSES = (RefineLs, ReflnL, ReflnSusceptibilityL, Pd2dPeakL, Pd2dProc)
    def __init__(self, background=None, resolution=None, meas=None, 
                 phase=None, diffrn_radiation=None,
                 setup=None,  range=None, chi2=None,
                 extinction=None, texture=None, 
                 exclude=None, asymmetry=None,
                 data_name=""):
        super(Pd2d, self).__init__(mandatory_classes=self.MANDATORY_CLASSES,
                                   optional_classes=self.OPTIONAL_CLASSES,
                                   internal_classes=self.INTERNAL_CLASSES)

        self.data_name = data_name
        self.background = background
        self.resolution = resolution
        self.meas = meas
        self.phase = phase
        self.diffrn_radiation = diffrn_radiation
        self.setup = setup
        self.range = range
        self.chi2 = chi2
        self.extinction = extinction
        self.exclude = exclude
        self.asymmetry = asymmetry
        self.texture = texture

        #FIXME: internal attributes
        #self.peaks = None
        #self.reflns = None
        #self.refln_ss = None
        self.__dd = None

        if self.is_defined:
            self.form_object


    @property
    def background(self):
        """
        """
        l_res = self[Pd2dBackground]
        if len(l_res) >= 1:
            return l_res[0]
        else:
            return None
    @background.setter
    def background(self, x):
        if x is None:
            pass
        elif isinstance(x, Pd2dBackground):
            l_ind = []
            for _i, _obj in enumerate(self.mandatory_objs):
                if isinstance(_obj, Pd2dBackground):
                    l_ind.append(_i)
            if len(l_ind) == 0:
                self.mandatory_objs.append(x)
            else:
                self.mandatory_objs[l_ind[0]] = x
            if len(l_ind) > 1:
                for _ind in l_ind.reverse():
                    self.mandatory_objs.pop(_ind)

    @property
    def resolution(self):
        """
        """
        l_res = self[Pd2dInstrResolution]
        if len(l_res) >= 1:
            return l_res[0]
        else:
            return None
    @resolution.setter
    def resolution(self, x):
        if x is None:
            pass
        elif isinstance(x, Pd2dInstrResolution):
            l_ind = []
            for _i, _obj in enumerate(self.mandatory_objs):
                if isinstance(_obj, Pd2dInstrResolution):
                    l_ind.append(_i)
            if len(l_ind) == 0:
                self.mandatory_objs.append(x)
            else:
                self.mandatory_objs[l_ind[0]] = x
            if len(l_ind) > 1:
                for _ind in l_ind.reverse():
                    self.mandatory_objs.pop(_ind)

    @property
    def meas(self):
        """
        """
        l_res = self[Pd2dMeas]
        if len(l_res) >= 1:
            return l_res[0]
        else:
            return None
    @meas.setter
    def meas(self, x):
        if x is None:
            pass
        elif isinstance(x, Pd2dMeas):
            l_ind = []
            for _i, _obj in enumerate(self.mandatory_objs):
                if isinstance(_obj, Pd2dMeas):
                    l_ind.append(_i)
            if len(l_ind) == 0:
                self.mandatory_objs.append(x)
            else:
                self.mandatory_objs[l_ind[0]] = x
            if len(l_ind) > 1:
                for _ind in l_ind.reverse():
                    self.mandatory_objs.pop(_ind)

    @property
    def phase(self):
        """
        """
        l_res = self[PhaseL]
        if len(l_res) >= 1:
            return l_res[0]
        else:
            return None
    @phase.setter
    def phase(self, x):
        if x is None:
            pass
        elif isinstance(x, PhaseL):
            l_ind = []
            for _i, _obj in enumerate(self.mandatory_objs):
                if isinstance(_obj, PhaseL):
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
        """
        """
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
    def setup(self):
        """
        """
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
    def range(self):
        """
        """
        l_res = self[Range]
        if len(l_res) >= 1:
            return l_res[0]
        else:
            return None
    @range.setter
    def range(self, x):
        if x is None:
            pass
        elif isinstance(x, Range):
            l_ind = []
            for _i, _obj in enumerate(self.mandatory_objs):
                if isinstance(_obj, Range):
                    l_ind.append(_i)
            if len(l_ind) == 0:
                self.mandatory_objs.append(x)
            else:
                self.mandatory_objs[l_ind[0]] = x
            if len(l_ind) > 1:
                for _ind in l_ind.reverse():
                    self.mandatory_objs.pop(_ind)

    @property
    def chi2(self):
        """
        """
        l_res = self[Chi2]
        if len(l_res) >= 1:
            return l_res[0]
        else:
            return None
    @chi2.setter
    def chi2(self, x):
        if x is None:
            pass
        elif isinstance(x, Chi2):
            l_ind = []
            for _i, _obj in enumerate(self.mandatory_objs):
                if isinstance(_obj, Chi2):
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
        """
        """
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
    def exclude(self):
        """
        """
        l_res = self[ExcludeL]
        if len(l_res) >= 1:
            return l_res[0]
        else:
            return None
    @exclude.setter
    def exclude(self, x):
        if x is None:
            pass
        elif isinstance(x, ExcludeL):
            l_ind = []
            for _i, _obj in enumerate(self.optional_objs):
                if isinstance(_obj, ExcludeL):
                    l_ind.append(_i)
            if len(l_ind) == 0:
                self.optional_objs.append(x)
            else:
                self.optional_objs[l_ind[0]] = x
            if len(l_ind) > 1:
                for _ind in l_ind.reverse():
                    self.optional_objs.pop(_ind)

    @property
    def asymmetry(self):
        """
        """
        l_res = self[Pd2dInstrReflexAsymmetry]
        if len(l_res) >= 1:
            return l_res[0]
        else:
            return None
    @asymmetry.setter
    def asymmetry(self, x):
        if x is None:
            pass
        elif isinstance(x, Pd2dInstrReflexAsymmetry):
            l_ind = []
            for _i, _obj in enumerate(self.optional_objs):
                if isinstance(_obj, Pd2dInstrReflexAsymmetry):
                    l_ind.append(_i)
            if len(l_ind) == 0:
                self.optional_objs.append(x)
            else:
                self.optional_objs[l_ind[0]] = x
            if len(l_ind) > 1:
                for _ind in l_ind.reverse():
                    self.optional_objs.pop(_ind)

    @property
    def texture(self):
        """
        """
        l_res = self[Texture]
        if len(l_res) >= 1:
            return l_res[0]
        else:
            return None
    @texture.setter
    def texture(self, x):
        if x is None:
            pass
        elif isinstance(x, Texture):
            l_ind = []
            for _i, _obj in enumerate(self.optional_objs):
                if isinstance(_obj, Texture):
                    l_ind.append(_i)
            if len(l_ind) == 0:
                self.optional_objs.append(x)
            else:
                self.optional_objs[l_ind[0]] = x
            if len(l_ind) > 1:
                for _ind in l_ind.reverse():
                    self.optional_objs.pop(_ind)


    @property
    def peak(self):
        l_res = self[Pd2dPeakL]
        return l_res

    @property
    def proc(self):
        l_res = self[Pd2dProc]
        if len(l_res) >= 1:
            return l_res[0]
        else:
            return None

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
        ls_out = ["Pd2d:"]
        ls_out.append(f"{str(self):}")
        return "\n".join(ls_out)

    @property
    def is_texture(self) -> bool:
        return (self.texture is not None)

    @property
    def is_variable_offset(self):
        l_bool = []
        l_bool.append(self.setup.offset_ttheta.refinement)
        l_bool.append(self.setup.offset_phi.refinement)
        res = any(l_bool)
        return res


    def calc_profile(self, tth, phi, l_crystal, l_peak_in=None, l_refln_in=None, l_refln_susceptibility_in=None, l_dd_in=None):
        """
        calculate intensity for the given diffraction angle
        """
        proc = Pd2dProc() #it is output
        setattr(proc, "__ttheta", tth)
        setattr(proc, "__phi", phi)
        
        background = self.background
        int_bkgd = background.interpolate_by_points(tth, phi)
        setattr(proc, "__intensity_bkg_calc", int_bkgd)

        tth_rad = tth*numpy.pi/180.
        phi_rad = phi*numpy.pi/180.
        cos_theta_1d = numpy.cos(0.5*tth_rad)
        sin_phi_1d = numpy.sin(phi_rad)

        setup = self.setup
        wavelength = float(setup.wavelength)
        diffrn_radiation = self.diffrn_radiation
        if setup.offset_phi is not None:
            phi_0 = float(setup.offset_phi)
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

        res_u_2d = numpy.zeros((tth.shape[0],phi.shape[0]), dtype=float)
        res_d_2d = numpy.zeros((tth.shape[0],phi.shape[0]), dtype=float)

        texture = self.texture 
        if texture is not None:
            h_ax, k_ax, l_ax = float(texture.h_ax), float(texture.k_ax), float(texture.l_ax)
            g_1, g_2 = float(texture.g_1), float(texture.g_2)

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
        for phase_label, phase_scale, phase_igsize, peak_in, refln_in, refln_susceptibility_in, dd_in in zip(
            phase.label, phase.scale, phase.igsize, l_peak_in, l_refln_in, l_refln_susceptibility_in, l_dd_in):
            dd_out = {}
            for i_crystal, crystal in enumerate(l_crystal):
                if crystal.data_name == phase_label:
                    ind_cry = i_crystal
                    break
            if ind_cry is None:
                self._show_message(f"Crystal with name '{phase_label:}' is not found.")
                return
            crystal = l_crystal[ind_cry]
            scale = float(phase_scale)
            i_g = float(phase_igsize)


            cell = crystal.cell
            space_group = crystal.space_group
            
            if peak_in is not None:
                h = numpy.array(peak_in.index_h, dtype=int)
                k = numpy.array(peak_in.index_k, dtype=int)
                l = numpy.array(peak_in.index_l, dtype=int)
                mult = numpy.array(peak_in.index_mult, dtype=int)
            else:
                if texture is None:
                    h, k, l, mult = cell.calc_hkl(space_group, sthovl_min, sthovl_max)
                else:
                    h, k, l, mult = cell.calc_hkl_in_range(sthovl_min, sthovl_max)

            peak = Pd2dPeakL(item=[Pd2dPeak(index_h=_h, index_k=_k, index_l=_l, index_mult=_m) for _h, _k, _l, _m in zip(h, k, l, mult)])
                
            cond_1 = not(crystal.is_variable)
            cond_2 = (peak_in is not None) & (refln_in is not None)
            
            if (cond_1 & cond_2):
                f_nucl_sq, f_m_p_sin_sq = peak_in.f_nucl_sq, peak_in.f_m_p_sin_sq
                f_m_p_cos_sq, cross_sin = peak_in.f_m_p_cos_sq, peak_in.cross_sin 
                refln = refln_in
                refln_s = refln_susceptibility_in
            else:
                f_nucl_sq, f_m_p_sin_sq, f_m_p_cos_sq, cross_sin, refln, refln_s = self.calc_for_iint(h, k, l, crystal)
                refln.loop_name = phase_label
                refln_s.loop_name = phase_label
            l_refln.append(refln)
            l_refln_s.append(refln_s)

            for _item, _1, _2, _3, _4 in zip(peak.item, f_nucl_sq, f_m_p_sin_sq, f_m_p_cos_sq, cross_sin):
                setattr(_item, "f_nucl_sq", _1) 
                setattr(_item, "f_m_p_sin_sq", _2) 
                setattr(_item, "f_m_p_cos_sq", _3) 
                setattr(_item, "cross_sin", _4)

            cond_1 = dd_in is not None
            cond_2 = ((not(crystal.is_variable)) & 
                      (not(diffrn_radiation.polarization.refinement)) & 
                      (not(diffrn_radiation.efficiency.refinement)))
            if cond_1 & cond_2:
                iint_u_3d, iint_d_3d = dd_in["iint_u_3d"], dd_in["iint_d_3d"]
                cos_theta_3d, sin_phi_3d = dd_in["cos_theta_3d"], dd_in["sin_phi_3d"]
            else:                      
                cos_theta_3d, sin_phi_3d, mult_f_n_3d = numpy.meshgrid(cos_theta_1d, sin_phi_1d, mult*f_nucl_sq, indexing="ij")

                mult_f_m_c_3d = numpy.meshgrid(tth_rad, phi_rad, mult*f_m_p_cos_sq, indexing="ij")[2]

                hh_u_s_3d = numpy.meshgrid(tth_rad, phi_rad, mult*(f_m_p_sin_sq+p_u*cross_sin), indexing="ij")[2]
                hh_d_s_3d = numpy.meshgrid(tth_rad, phi_rad, mult*(f_m_p_sin_sq-p_d*cross_sin), indexing="ij")[2]

                c_a_sq_3d = (cos_theta_3d * sin_phi_3d)**2
                s_a_sq_3d = 1.-c_a_sq_3d
                iint_u_3d = (mult_f_n_3d + hh_u_s_3d*s_a_sq_3d + mult_f_m_c_3d*c_a_sq_3d)
                iint_d_3d = (mult_f_n_3d + hh_d_s_3d*s_a_sq_3d + mult_f_m_c_3d*c_a_sq_3d)
            dd_out["cos_theta_3d"] = cos_theta_3d
            dd_out["sin_phi_3d"] = sin_phi_3d
            dd_out["iint_u_3d"] = iint_u_3d
            dd_out["iint_d_3d"] = iint_d_3d

            sthovl_hkl = cell.calc_sthovl(h, k, l)
            tth_hkl_rad = numpy.where(sthovl_hkl*wavelength < 1., 2.*numpy.arcsin(sthovl_hkl*wavelength), numpy.pi)
            tth_hkl = tth_hkl_rad*180./numpy.pi

            cond_1 = dd_in is not None
            cond_2 = ((not(phase_igsize.refinement)) & 
                      (not(self.resolution.is_variable)) & 
                      (not(self.is_variable_offset)) &
                      (not(setup.wavelength.refinement)))
            if cond_1 & cond_2:
                profile_3d, tth_zs, h_pv = dd_in["profile_3d"], dd_in["tth_zs"], dd_in["h_pv"]
            else:
                profile_3d, tth_zs, h_pv = self.calc_shape_profile(tth, phi, tth_hkl, i_g)
            dd_out["profile_3d"] = profile_3d
            dd_out["tth_zs"] = tth_zs
            dd_out["h_pv"] = h_pv

            for _item, _1, _2 in zip(peak.item, tth_hkl, h_pv):
                setattr(_item, "ttheta", _1)
                setattr(_item, "width_ttheta", _2)

            #texture
            if texture is not None:
                cond_1 = dd_in is not None
                cond_2 = ((not(setup.offset_phi.refinement)))
                if cond_1 & cond_2:
                    cos_alpha_ang_3d, sin_alpha_ang_3d = dd_in["cos_alpha_ang_3d"], dd_in["sin_alpha_ang_3d"]
                else:
                    cos_alpha_ang_3d = cos_theta_3d * sin_phi_3d  
                    sin_alpha_ang_3d = numpy.sqrt(1.-cos_alpha_ang_3d**2)
                dd_out["cos_alpha_ang_3d"], dd_out["sin_alpha_ang_3d"] = cos_alpha_ang_3d, sin_alpha_ang_3d

                cond_2 = (cond_2 & (not(texture.is_variable)) & (not(cell.is_variable)))
                if cond_1 & cond_2:
                    texture_3d = dd_in["texture_3d"]
                else:
                    cos_alpha_ax = calc_cos_ang(cell, h_ax, k_ax, l_ax, h, k, l)
                    c_help = 1.-cos_alpha_ax**2
                    c_help[c_help<0.] = 0.
                    sin_alpha_ax = numpy.sqrt(c_help)
                    cos_alpha_ax_3d = numpy.meshgrid(tth, phi, cos_alpha_ax, indexing="ij")[2]
                    sin_alpha_ax_3d = numpy.meshgrid(tth, phi, sin_alpha_ax, indexing="ij")[2]
                    cos_alpha_3d = cos_alpha_ax_3d*cos_alpha_ang_3d+sin_alpha_ax_3d*sin_alpha_ang_3d
                    texture_3d = g_2 + (1.-g_2)*(1./g_1 +
                                       (g_1**2-1./g_1)*cos_alpha_3d**2)**(-1.5)
                dd_out["texture_3d"] = texture_3d

                profile_3d = profile_3d*texture_3d



            res_u_3d = profile_3d*iint_u_3d 
            res_d_3d = profile_3d*iint_d_3d 

            res_u_2d += scale*res_u_3d.sum(axis=2) 
            res_d_2d += scale*res_d_3d.sum(axis=2) 
            l_peak.append(peak)
            l_dd_out.append(dd_out)


        proc.ttheta_corrected = tth_zs
        setattr(proc, "__intensity_up_net", res_u_2d)
        setattr(proc, "__intensity_down_net", res_d_2d)
        setattr(proc, "__intensity_up_total", res_u_2d+int_bkgd)
        setattr(proc, "__intensity_down_total", res_d_2d+int_bkgd)
        l_internal_objs = l_refln + l_refln_s + l_peak
        l_internal_objs.append(proc)
        setattr(self, "__internal_objs", l_internal_objs)

        return proc, l_peak, l_refln, l_dd_out
    #def calc_y_mod(self, l_crystal, d_prof_in={}):
    #    """
    #    calculate model diffraction profiles up and down if observed data is defined
    #    """
    #    observed_data = self._p_observed_data
    #    tth = observed_data.get_val('tth')
    #
    #    wavelength = observed_data.get_val('wavelength')
    #    setup = self._p_setup
    #    setup.set_val(wavelength=wavelength)
    #
    #    field = observed_data.get_val('field')
    #    for calculated_data in self._list_calculated_data:
    #        calculated_data.set_val(field=field)
    #
    #    int_u_mod, int_d_mod, d_exp_prof_out = self.calc_profile(tth, l_crystal, d_prof_in)
    #    return int_u_mod, int_d_mod, d_exp_prof_out

    def calc_chi_sq(self, l_crystal):
        """
Calculate chi square
        """
        meas = self.meas

        tth = meas.ttheta
        phi = meas.phi
        int_u_exp = meas.intensity_up
        sint_u_exp = meas.intensity_up_sigma
        int_d_exp = meas.intensity_down
        sint_d_exp = meas.intensity_down_sigma

        #if ((self.peaks is not None) & (self.reflns is not None)):
        #    l_peak_in = self.peaks
        #    l_refln_in = self.reflns
        #    l_dd_in = self.__dd
        if ((len(self.peak)!=0) & (len(self.refln) != 0) & (len(self.refln_susceptibility) != 0)):
            l_peak_in = self.peak
            l_refln_in = self.refln
            l_refln_susceptibility_in = self.refln_susceptibility
            l_dd_in = self.__dd
        else:
            l_peak_in, l_refln_in, l_dd_in = [], [], [] 
            l_refln_susceptibility_in = []

        range_ = self.range
        cond_tth_in = numpy.ones(tth.size, dtype=bool)
        cond_tth_in = numpy.logical_and(cond_tth_in, tth >= range_.ttheta_min)
        cond_tth_in = numpy.logical_and(cond_tth_in, tth <= range_.ttheta_max)

        cond_phi_in = numpy.ones(phi.size, dtype=bool)
        cond_phi_in = numpy.logical_and(cond_phi_in, phi >= range_.phi_min)
        cond_phi_in = numpy.logical_and(cond_phi_in, phi <= range_.phi_max)

        #cond_1_in, cond_2_in = numpy.meshgrid(cond_tth_in, cond_phi_in, indexing="ij")
        #cond_in = numpy.logical_and(cond_1_in, cond_2_in)
        tth_in = tth[cond_tth_in]
        phi_in = phi[cond_phi_in]
        int_u_exp_in = int_u_exp[cond_tth_in, :][:, cond_phi_in]
        sint_u_exp_in = sint_u_exp[cond_tth_in, :][:, cond_phi_in]
        int_d_exp_in = int_d_exp[cond_tth_in, :][:, cond_phi_in]
        sint_d_exp_in = sint_d_exp[cond_tth_in, :][:, cond_phi_in]


        proc, l_peak, l_refln, l_dd_out = self.calc_profile(tth_in, phi_in, l_crystal, 
                                                     l_peak_in=l_peak_in, l_refln_in=l_refln_in, 
                                                     l_refln_susceptibility_in=l_refln_susceptibility_in, l_dd_in=l_dd_in)
        setattr(self.proc, "__intensity_up", int_u_exp_in)
        setattr(self.proc, "__intensity_up_sigma", sint_u_exp_in)
        setattr(self.proc, "__intensity_down", int_d_exp_in)
        setattr(self.proc, "__intensity_down_sigma", sint_d_exp_in)

        #self.proc = proc
        #self.peaks = l_peak
        #self.reflns = l_refln
        self.__dd = l_dd_out

        int_u_mod = proc.intensity_up_total
        int_d_mod = proc.intensity_down_total

        sint_sum_exp_in = (sint_u_exp_in**2 + sint_d_exp_in**2)**0.5

        chi_sq_u = ((int_u_mod-int_u_exp_in)/sint_u_exp_in)**2
        chi_sq_d = ((int_d_mod-int_d_exp_in)/sint_d_exp_in)**2
        chi_sq_sum = ((int_u_mod+int_d_mod-int_u_exp_in-int_d_exp_in)/sint_sum_exp_in)**2
        chi_sq_dif = ((int_u_mod-int_d_mod-int_u_exp_in+int_d_exp_in)/sint_sum_exp_in)**2


        cond_u = numpy.logical_not(numpy.isnan(chi_sq_u))
        cond_d = numpy.logical_not(numpy.isnan(chi_sq_d))
        cond_sum = numpy.logical_not(numpy.isnan(chi_sq_sum))
        cond_dif = numpy.logical_not(numpy.isnan(chi_sq_dif))

        #exclude region
        exclude = self.exclude
        if exclude is not None:
            l_excl_tth_min = exclude.ttheta_min
            l_excl_tth_max = exclude.ttheta_max
            l_excl_phi_min = exclude.phi_min
            l_excl_phi_max = exclude.phi_max
            for excl_tth_min, excl_tth_max, excl_phi_min, excl_phi_max in zip(l_excl_tth_min, l_excl_tth_max, l_excl_phi_min, l_excl_phi_max):
                cond_1 = numpy.logical_or(tth_in < 1.*excl_tth_min, tth_in > 1.*excl_tth_max)
                cond_2 = numpy.logical_or(phi_in < 1.*excl_phi_min, phi_in > 1.*excl_phi_max)
                cond_11, cond_22 = numpy.meshgrid(cond_1, cond_2, indexing="ij")
                cond_12 = numpy.logical_or(cond_11, cond_22)
                cond_u = numpy.logical_and(cond_u, cond_12)
                cond_d = numpy.logical_and(cond_d, cond_12)
                cond_sum = numpy.logical_and(cond_sum, cond_12)


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
                  int(flag_sum)*chi_sq_sum_val + int(flag_dif)*chi_sq_dif_val)
        n = (int(flag_u)*n_u + int(flag_d)*n_d + int(flag_sum)*n_sum + 
             int(flag_dif)*n_dif)
        #print(f"chi_sq_val/n: {chi_sq_val/n:.2f}           chi_sq_val: {chi_sq_val: .2f}")
        #d_exp_out = {"chi_sq_val": chi_sq_val, "n": n}
        #d_exp_out.update(d_exp_prof_out)

        refine_ls = RefineLs(number_reflns=n, goodness_of_fit_all=chi_sq_val/float(n), weighting_scheme="sigma")
        self.internal_objs.append(refine_ls)  
        return chi_sq_val, n


    def calc_for_iint(self, h, k, l, crystal):
        """
Calculate the integral intensity for h, k, l reflections
        """
        setup = self.setup
        field = float(setup.field)
        diffrn_radiation = self.diffrn_radiation
        p_u = float(diffrn_radiation.polarization)
        p_d = (2.*float(diffrn_radiation.efficiency)-1)*p_u

        refln = crystal.calc_refln(h, k, l)

        refln_s = crystal.calc_refln_susceptibility(h, k, l)

        f_nucl = numpy.array(refln.f_calc, dtype=complex)
        sft_11, sft_12, sft_13 = numpy.array(refln_s.chi_11_calc, dtype=complex), numpy.array(refln_s.chi_12_calc, dtype=complex), numpy.array(refln_s.chi_13_calc, dtype=complex)
        sft_21, sft_22, sft_23 = numpy.array(refln_s.chi_21_calc, dtype=complex), numpy.array(refln_s.chi_22_calc, dtype=complex), numpy.array(refln_s.chi_23_calc, dtype=complex)
        sft_31, sft_32, sft_33 = numpy.array(refln_s.chi_31_calc, dtype=complex), numpy.array(refln_s.chi_32_calc, dtype=complex), numpy.array(refln_s.chi_33_calc, dtype=complex)

        sftm_11, sftm_12, sftm_13 = numpy.array(refln_s.moment_11_calc, dtype=complex), numpy.array(refln_s.moment_12_calc, dtype=complex), numpy.array(refln_s.moment_13_calc, dtype=complex)
        sftm_21, sftm_22, sftm_23 = numpy.array(refln_s.moment_21_calc, dtype=complex), numpy.array(refln_s.moment_22_calc, dtype=complex), numpy.array(refln_s.moment_23_calc, dtype=complex)
        sftm_31, sftm_32, sftm_33 = numpy.array(refln_s.moment_31_calc, dtype=complex), numpy.array(refln_s.moment_32_calc, dtype=complex), numpy.array(refln_s.moment_33_calc, dtype=complex)

        _11, _12, _13 = sftm_11+field*sft_11, sftm_12+field*sft_12, sftm_13+field*sft_13
        _21, _22, _23 = sftm_21+field*sft_21, sftm_22+field*sft_22, sftm_23+field*sft_23
        _31, _32, _33 = sftm_31+field*sft_31, sftm_32+field*sft_32, sftm_33+field*sft_33

        cell = crystal.cell
        #k_loc = cell.calc_k_loc(h, k, l)
        t_11, t_12, t_13, t_21, t_22, t_23, t_31, t_32, t_33 = cell.calc_m_t(h, k, l)
        
        th_11, th_12, th_13, th_21, th_22, th_23, th_31, th_32, th_33 = CONSTANTS_AND_FUNCTIONS.calc_mRmCmRT(
                t_11, t_21, t_31, t_12, t_22, t_32, t_13, t_23, t_33,
                _11, _12, _13, _21, _22, _23, _31, _32, _33)

        #f_m_p_sin_sq = (field**2)*abs(0.5*(th_11*th_11.conjugate()+th_22*th_22.conjugate())+th_12*th_12.conjugate())
        #f_m_p_cos_sq = (field**2)*abs(th_13*th_13.conjugate()+th_23*th_23.conjugate())
        #f_m_p_field = 0.5*field*(th_11+th_22) 

        f_nucl_sq = abs(f_nucl*f_nucl.conjugate())
        f_m_p_sin_sq = abs(0.5*(th_11*th_11.conjugate()+th_22*th_22.conjugate())+th_12*th_12.conjugate())
        f_m_p_cos_sq = abs(th_13*th_13.conjugate()+th_23*th_23.conjugate())
        f_m_p_field = 0.5*(th_11+th_22) 
        cross_sin = 2.*(f_nucl.real*f_m_p_field.real+f_nucl.imag*f_m_p_field.imag)
        
        return f_nucl_sq, f_m_p_sin_sq, f_m_p_cos_sq, cross_sin, refln, refln_s

    def _gauss_pd(self, tth_2d):
        """
        one dimensional gauss powder diffraction
        """
        ag, bg = self.__ag, self.__bg
        val_1 = bg*tth_2d**2
        val_2 = numpy.where(val_1 < 5., numpy.exp(-val_1), 0.)
        self.__gauss_pd = ag*val_2

    def _lor_pd(self, tth_2d):
        """
        one dimensional lorentz powder diffraction
        """
        al, bl = self.__al, self.__bl
        self.__lor_pd = al*1./(1.+bl*tth_2d**2)
    

    def calc_shape_profile(self, tth, phi, tth_hkl, i_g=0.):
        """
        calculate profile in the range ttheta for reflections placed on 
        ttheta_hkl with i_g parameter by default equal to zero
        
        tth, phi, tth_hkl in degrees

        """
        setup = self.setup
        zero_shift = float(setup.offset_ttheta)
        tth_zs = tth-zero_shift

        resolution = self.resolution
        
        h_pv, eta, h_g, h_l, a_g, b_g, a_l, b_l = resolution.calc_resolution(
                                                  tth_hkl, i_g=i_g)


        tth_3d, phi_3d, tth_hkl_3d = numpy.meshgrid(tth_zs, phi, tth_hkl, indexing="ij")


        self.__ag = numpy.meshgrid(tth_zs, phi, a_g, indexing="ij")[2]
        self.__bg = numpy.meshgrid(tth_zs, phi, b_g, indexing="ij")[2]
        self.__al = numpy.meshgrid(tth_zs, phi, a_l, indexing="ij")[2]
        self.__bl = numpy.meshgrid(tth_zs, phi, b_l, indexing="ij")[2]
        eta_3d = numpy.meshgrid(tth_zs, phi, eta, indexing="ij")[2]
        self.__eta = eta_3d 

        self._gauss_pd(tth_3d-tth_hkl_3d)
        self._lor_pd(tth_3d-tth_hkl_3d)
        g_pd2d_3d = self.__gauss_pd 
        l_pd2d_3d = self.__lor_pd
        
        np_shape_3d = eta_3d * l_pd2d_3d + (1.-eta_3d) * g_pd2d_3d

        asymmetry = self.asymmetry
        np_ass_2d = asymmetry.calc_asymmetry(tth_zs, tth_hkl)
        np_ass_3d = np_ass_2d[:, numpy.newaxis,:]*numpy.ones(phi.size, dtype=float)[numpy.newaxis, :, numpy.newaxis]
       
        #Lorentz factor
        tth_rad = tth_zs*numpy.pi/180.
        np_lor_1d = 1./(numpy.sin(tth_rad)*numpy.sin(0.5*tth_rad))
        np_lor_3d = numpy.meshgrid(np_lor_1d, phi, tth_hkl, indexing="ij")[0]
        
        
        profile_3d = np_shape_3d*np_ass_3d*np_lor_3d
        
        return profile_3d, tth_zs, h_pv


    def params_to_cif(self, separator="_", flag=False) -> str: 
        ls_out = []
        l_cls = (Pd2dBackground, Pd2dInstrResolution, PhaseL, DiffrnRadiation, Setup, Range, Chi2,
                 Extinction, ExcludeL, Pd2dInstrReflexAsymmetry, Pd2dPeakL, Texture)
        l_obj = [_obj for _obj in (self.mandatory_objs + self.optional_objs) if type(_obj) in l_cls]
        ls_out.extend([_.to_cif(separator=separator, flag=flag)+"\n" for _ in l_obj])
        return "\n".join(ls_out)

    def data_to_cif(self, separator="_", flag=False) -> str: 
        ls_out = []
        l_cls = (Pd2dMeas, )
        l_obj = [_obj for _obj in (self.mandatory_objs + self.optional_objs) if type(_obj) in l_cls]
        ls_out.extend([_.to_cif(separator=separator, flag=flag)+"\n" for _ in l_obj])
        return "\n".join(ls_out)

    def calc_to_cif(self, separator="_", flag=False) -> str: 

        ls_out = []
        l_cls = (RefineLs, ReflnL, ReflnSusceptibilityL, Pd2dPeakL, Pd2dProc)
        for _cls in l_cls:
            l_obj = [_obj for _obj in (self.optional_objs+self.internal_objs) if isinstance(_obj, _cls)]
            ls_out.extend([_.to_cif(separator=separator, flag=flag)+"\n" for _ in l_obj])


        #ls_out = []
        #l_cls = (Pd2dProc, )
        #l_obj = [_obj for _obj in (self.mandatory_objs + self.optional_objs) if type(_obj) in l_cls]
        #ls_out.extend([_.to_cif(separator=separator, flag=flag)+"\n" for _ in l_obj])
        #if self.reflns is not None:
        #    ls_out.extend([_.to_cif(separator=separator, flag=flag)+"\n" for _ in self.reflns])
        #if self.refln_ss is not None:
        #    ls_out.extend([_.to_cif(separator=separator, flag=flag)+"\n" for _ in self.refln_ss])

        return "\n".join(ls_out)

    @property
    def remove_internal_objs(self):
        l_objs = self.internal_objs
        for _objs in l_objs:
            l_objs.remove(_objs)
        self.__dd = None            