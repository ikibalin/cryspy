__author__ = 'ikibalin'
__version__ = "2019_12_11"
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

from cryspy.pd1dcif_like.cl_pd_background import PdBackground, PdBackgroundL
from cryspy.pd1dcif_like.cl_pd_exclude import PdExclude, PdExcludeL
from cryspy.pd1dcif_like.cl_pd_instr_reflex_asymmetry import PdInstrReflexAsymmetry
from cryspy.pd1dcif_like.cl_pd_instr_resolution import PdInstrResolution
from cryspy.pd1dcif_like.cl_pd_meas import PdMeas, PdMeasL
from cryspy.pd1dcif_like.cl_pd_peak import PdPeak, PdPeakL
from cryspy.pd1dcif_like.cl_pd_proc import PdProc, PdProcL

from .cl_phase import Phase, PhaseL
from .cl_diffrn_radiation import DiffrnRadiation
from .cl_extinction import Extinction
from .cl_setup import Setup
from .cl_range import Range
from .cl_chi2 import Chi2

class Pd(DataConstr):
    """
Class to describe information about powder diffraction measurements

Description in cif file::

 data_pnd
 _setup_wavelength     0.840
 _setup_field          1.000
 _setup_offset_ttheta -0.385
 
 _chi2_sum  True
 _chi2_diff False
 _chi2_up   False
 _chi2_down False
 
 _range_ttheta_min     4.000
 _range_ttheta_max    80.000
 
 loop_
 _pd_background_ttheta
 _pd_background_intensity
  4.5 256.0
 40.0 158.0
 80.0  65.0
 
 loop_
 _pd_exclude_ttheta_min
 _pd_exclude_ttheta_max
 0.0 1.0
 
 _pd_instr_reflex_asymmetry_p1 0.0
 _pd_instr_reflex_asymmetry_p2 0.0
 _pd_instr_reflex_asymmetry_p3 0.0
 _pd_instr_reflex_asymmetry_p4 0.0
 
 _diffrn_radiation_polarization 0.0
 _diffrn_radiation_efficiency   1.0
 
 _pd_instr_resolution_u 16.9776
 _pd_instr_resolution_v -2.8357
 _pd_instr_resolution_w  0.5763
 _pd_instr_resolution_x  0.0
 _pd_instr_resolution_y  0.0
 
 loop_
 _phase_label
 _phase_scale
 _phase_igsize
 Fe3O4 0.02381 0.0
 
 loop_
 _pd_meas_ttheta
 _pd_meas_intensity_up
 _pd_meas_intensity_up_sigma
 _pd_meas_intensity_down
 _pd_meas_intensity_down_sigma
 4.0 465.80 128.97 301.88 129.30
 4.2 323.78 118.22 206.06 120.00
 4.4 307.14 115.90 230.47 116.53
  
 loop_
 _pd_peak_index_h
 _pd_peak_index_k
 _pd_peak_index_l
 _pd_peak_mult
 _pd_peak_ttheta
 _pd_peak_intensity_up
 _pd_peak_intensity_down
 _pd_peak_width_2theta
 1 1 1  8  9.748 128.15576 128.15576 0.677
 2 0 0  6 11.260   0.00000   0.00000 0.680
 2 2 0 12 15.950  94.21107  94.21107 0.716
 
 loop_
 _pd_proc_ttheta
 _pd_proc_ttheta_corrected
 _pd_proc_intensity_up_net
 _pd_proc_intensity_down_net
 _pd_proc_intensity_up_total
 _pd_proc_intensity_down_total
 _pd_proc_intensity_bkg_calc
 _pd_proc_intensity_up
 _pd_proc_intensity_up_sigma
 _pd_proc_intensity_down
 _pd_proc_intensity_down_sigma
 4.000 4.385 0.00000 0.00000 256.00000 256.00000 256.00000 465.80000 128.97000 301.88000 129.30000
 4.200 4.585 0.00000 0.00000 256.00000 256.00000 256.00000 323.78000 118.22000 206.06000 120.00000
 4.400 4.785 0.00000 0.00000 256.00000 256.00000 256.00000 307.14000 115.90000 230.47000 116.53000
     """
    MANDATORY_CLASSES = (PdBackgroundL, PdInstrResolution, PdMeasL, PhaseL, DiffrnRadiation,
                         Setup, Range, Chi2)
    OPTIONAL_CLASSES = (Extinction, PdExcludeL, PdInstrReflexAsymmetry, PdPeakL, PdProcL)
    INTERNAL_CLASSES = ()
    def __init__(self, background=None, resolution=None, meas=None, 
                 phase=None, diffrn_radiation=None,
                 setup=None,  range=None, chi2=None,
                 extinction=None,
                 exclude=None, asymmetry=None, peak=None, proc=None,
                 data_name=""):
        super(Pd, self).__init__(mandatory_classes=self.MANDATORY_CLASSES,
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
        self.peak = peak
        self.proc = proc

        if self.is_defined:
            self.form_object

    @property
    def background(self):
        """
        """
        l_res = self[PdBackgroundL]
        if len(l_res) >= 1:
            return l_res[0]
        else:
            return None
    @background.setter
    def background(self, x):
        if x is None:
            pass
        elif isinstance(x, PdBackgroundL):
            l_ind = []
            for _i, _obj in enumerate(self.mandatory_objs):
                if isinstance(_obj, PdBackgroundL):
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
        l_res = self[PdInstrResolution]
        if len(l_res) >= 1:
            return l_res[0]
        else:
            return None
    @resolution.setter
    def resolution(self, x):
        if x is None:
            pass
        elif isinstance(x, PdInstrResolution):
            l_ind = []
            for _i, _obj in enumerate(self.mandatory_objs):
                if isinstance(_obj, PdInstrResolution):
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
        l_res = self[PdMeasL]
        if len(l_res) >= 1:
            return l_res[0]
        else:
            return None
    @meas.setter
    def meas(self, x):
        if x is None:
            pass
        elif isinstance(x, PdMeasL):
            l_ind = []
            for _i, _obj in enumerate(self.mandatory_objs):
                if isinstance(_obj, PdMeasL):
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
        l_res = self[PdExcludeL]
        if len(l_res) >= 1:
            return l_res[0]
        else:
            return None
    @exclude.setter
    def exclude(self, x):
        if x is None:
            pass
        elif isinstance(x, PdExcludeL):
            l_ind = []
            for _i, _obj in enumerate(self.optional_objs):
                if isinstance(_obj, PdExcludeL):
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
        l_res = self[PdInstrReflexAsymmetry]
        if len(l_res) >= 1:
            return l_res[0]
        else:
            return None
    @asymmetry.setter
    def asymmetry(self, x):
        if x is None:
            pass
        elif isinstance(x, PdInstrReflexAsymmetry):
            l_ind = []
            for _i, _obj in enumerate(self.optional_objs):
                if isinstance(_obj, PdInstrReflexAsymmetry):
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
        """
        """
        l_res = self[PdPeakL]
        if len(l_res) >= 1:
            return l_res[0]
        else:
            return None
    @peak.setter
    def peak(self, x):
        if x is None:
            pass
        elif isinstance(x, PdPeakL):
            l_ind = []
            for _i, _obj in enumerate(self.optional_objs):
                if isinstance(_obj, PdPeakL):
                    l_ind.append(_i)
            if len(l_ind) == 0:
                self.optional_objs.append(x)
            else:
                self.optional_objs[l_ind[0]] = x
            if len(l_ind) > 1:
                for _ind in l_ind.reverse():
                    self.optional_objs.pop(_ind)

    @property
    def proc(self):
        """
        """
        l_res = self[PdProcL]
        if len(l_res) >= 1:
            return l_res[0]
        else:
            return None
    @proc.setter
    def proc(self, x):
        if x is None:
            pass
        elif isinstance(x, PdProcL):
            l_ind = []
            for _i, _obj in enumerate(self.optional_objs):
                if isinstance(_obj, PdProcL):
                    l_ind.append(_i)
            if len(l_ind) == 0:
                self.optional_objs.append(x)
            else:
                self.optional_objs[l_ind[0]] = x
            if len(l_ind) > 1:
                for _ind in l_ind.reverse():
                    self.optional_objs.pop(_ind)

    def __repr__(self):
        ls_out = ["Pd:"]
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


    def calc_profile(self, tth, l_crystal, l_peak_in=[], l_refln_in=[]):
        """
        calculate intensity for the given diffraction angle
        """
        proc = PdProc() #it is output
        proc.ttheta = tth
        
        background = self.background
        int_bkgd = background.interpolate_by_points(tth)
        proc.bkg_calc = int_bkgd


        wavelength = float(self.wavelength)
        beam_polarization = self.beam_polarization
        
        tth_min = tth.min()
        tth_max = tth.max()+3. 
        if tth_max > 180.:
            tth_max = 180.
        sthovl_min = numpy.sin(0.5*tth_min*numpy.pi/180.)/wavelength
        sthovl_max = numpy.sin(0.5*tth_max*numpy.pi/180.)/wavelength

        res_u_1d = numpy.zeros(tth.shape[0], dtype=float)
        res_d_1d = numpy.zeros(tth.shape[0], dtype=float)

        phase = self.phase
        if len(phase.label) != len(l_peak_in):
            l_peak_in = len(phase.label)*[None]
        if len(phase.label) != len(l_refln_in):
            l_refln_in = len(phase.label)*[None]
        l_peak, l_refln = [], []
        for phase_label, phase_scale, phase_igsize, peak_in, refln_in in zip(
                                 phase.label, phase.scale, phase.igsize, l_peak_in, l_refln_in):
            for i_crystal, crystal in enumerate(l_crystal):
                if crystal.label == phase_label:
                    ind_cry = i_crystal
                    break
            if ind_cry is None:
                print("Crystal with name '{:}' is not found.".format(
                        phase_label))
                return
            crystal = l_crystal[ind_cry]
            scale = float(phase_scale)
            i_g = float(phase_igsize)


            cell = crystal.cell
            space_group = crystal.space_group
            
            peak = PdPeak()
            if peak_in is not None:
                h, k, l, mult = peak_in.h, peak_in.k, peak_in.l, peak_in.mult 
            else:
                h, k, l, mult = cell.calc_hkl(space_group, sthovl_min, sthovl_max)
            peak.h, peak.k, peak.l, peak.mult = h, k, l, mult

            cond_1 = not(crystal.is_variable)
            cond_2 = (peak_in is not None) & (refln_in is not None)
            if cond_1 & cond_2:
                np_iint_u, np_iint_d, refln = peak_in.up, peak_in.down, refln_in
            else:   
                np_iint_u, np_iint_d, refln = self.calc_iint(h, k, l, crystal)

            l_refln.append(refln)
            peak.up, peak.down = np_iint_u, np_iint_d

            sthovl_hkl = cell.calc_sthovl(h, k, l)
            
            tth_hkl_rad = 2.*numpy.arcsin(sthovl_hkl*wavelength)
            tth_hkl = tth_hkl_rad*180./numpy.pi
            peak.ttheta = tth_hkl

            profile_2d, tth_zs, h_pv = self.calc_shape_profile(tth, tth_hkl, i_g)
            peak.width_2theta = h_pv
            
            np_iint_u_2d = numpy.meshgrid(tth, np_iint_u*mult, indexing="ij")[1]
            np_iint_d_2d = numpy.meshgrid(tth, np_iint_d*mult, indexing="ij")[1]
            
            res_u_2d = profile_2d*np_iint_u_2d 
            res_d_2d = profile_2d*np_iint_d_2d 

            res_u_1d += scale*res_u_2d.sum(axis=1) 
            res_d_1d += scale*res_d_2d.sum(axis=1) 
            l_peak.append(peak)
        proc.ttheta_corrected = tth_zs
        proc.up_net = res_u_1d
        proc.down_net = res_d_1d
        proc.up_total = res_u_1d+int_bkgd
        proc.down_total = res_d_1d+int_bkgd
        return proc, l_peak, l_refln
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
        calculate chi square
        """
        meas = self.meas

        tth = meas.ttheta
        int_u_exp = meas.up
        sint_u_exp = meas.up_sigma
        int_d_exp = meas.down
        sint_d_exp = meas.down_sigma

        if ((self.peaks is not None) & (self.reflns is not None)):
            l_peak_in = self.peaks
            l_refln_in = self.reflns
        else:
            l_peak_in, l_refln_in = [], [] 


        cond_in = numpy.ones(tth.shape, dtype=bool)
        cond_in = numpy.logical_and(cond_in, tth >= self.range_min)
        cond_in = numpy.logical_and(cond_in, tth <= self.range_max)

        tth_in = tth[cond_in]
        int_u_exp_in = int_u_exp[cond_in]
        sint_u_exp_in = sint_u_exp[cond_in]
        int_d_exp_in = int_d_exp[cond_in]
        sint_d_exp_in = sint_d_exp[cond_in]
        
        proc, l_peak, l_refln = self.calc_profile(tth_in, l_crystal, l_peak_in, l_refln_in)
        
        proc.up = int_u_exp_in
        proc.up_sigma = sint_u_exp_in
        proc.down = int_d_exp_in
        proc.down_sigma = sint_d_exp_in
        self.proc = proc
        self.peaks = l_peak
        self.reflns = l_refln

        int_u_mod = proc.up_total
        int_d_mod = proc.down_total

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
        exclude_2theta = self.exclude_2theta
        if exclude_2theta is not None:
            l_excl_tth_min = exclude_2theta.min
            l_excl_tth_max = exclude_2theta.max
            for excl_tth_min, excl_tth_max in zip(l_excl_tth_min, l_excl_tth_max):
                cond_1 = numpy.logical_or(tth_in < 1.*excl_tth_min, tth_in > 1.*excl_tth_max)
                cond_u = numpy.logical_and(cond_u, cond_1)
                cond_d = numpy.logical_and(cond_d, cond_1)
                cond_sum = numpy.logical_and(cond_sum, cond_1)
            

        chi_sq_u_val = (chi_sq_u[cond_u]).sum()
        n_u = cond_u.sum()
        
        chi_sq_d_val = (chi_sq_d[cond_d]).sum()
        n_d = cond_d.sum()

        chi_sq_sum_val = (chi_sq_sum[cond_sum]).sum()
        n_sum = cond_sum.sum()

        chi_sq_dif_val = (chi_sq_dif[cond_dif]).sum()
        n_dif = cond_dif.sum()

        flag_u = self.chi2_up
        flag_d = self.chi2_down
        flag_sum = self.chi2_sum
        flag_dif = self.chi2_diff
        
        chi_sq_val = (int(flag_u)*chi_sq_u_val + int(flag_d)*chi_sq_d_val + 
                  int(flag_sum)*chi_sq_sum_val + int(flag_dif)*chi_sq_dif_val)
        n = (int(flag_u)*n_u + int(flag_d)*n_d + int(flag_sum)*n_sum + 
             int(flag_dif)*n_dif)
        #d_exp_out = {"chi_sq_val": chi_sq_val, "n": n}
        #d_exp_out.update(d_exp_prof_out)
        return chi_sq_val, n


    def calc_iint(self, h, k, l, crystal):
        """
        calculate the integral intensity for h, k, l reflections
        """
        field = float(self.field)
        beam_polarization = self.beam_polarization
        p_u = float(beam_polarization.polarization)
        p_d = (2.*float(beam_polarization.efficiency)-1)*p_u

        refln = crystal.calc_sf(h, k, l)

        f_nucl, sft_11, sft_12, sft_13 = refln.f_nucl, refln.sft_11, refln.sft_12, refln.sft_13
        sft_21, sft_22, sft_23 = refln.sft_21, refln.sft_22, refln.sft_23
        sft_31, sft_32, sft_33 = refln.sft_31, refln.sft_32, refln.sft_33

        sftm_11, sftm_12, sftm_13 = refln.sftm_11, refln.sftm_12, refln.sftm_13
        sftm_21, sftm_22, sftm_23 = refln.sftm_21, refln.sftm_22, refln.sftm_23
        sftm_31, sftm_32, sftm_33 = refln.sftm_31, refln.sftm_32, refln.sftm_33

        _11, _12, _13 = sftm_11+field*sft_11, sftm_12+field*sft_12, sftm_13+field*sft_13
        _21, _22, _23 = sftm_21+field*sft_21, sftm_22+field*sft_22, sftm_23+field*sft_23
        _31, _32, _33 = sftm_31+field*sft_31, sftm_32+field*sft_32, sftm_33+field*sft_33

        cell = crystal.cell
        
        #k_loc = cell.calc_k_loc(h, k, l)
        t_11, t_12, t_13, t_21, t_22, t_23, t_31, t_32, t_33 = cell.calc_m_t(h, k, l)
        
        th_11, th_12, th_13, th_21, th_22, th_23, th_31, th_32, th_33 = calc_mRmCmRT(
                t_11, t_21, t_31, t_12, t_22, t_32, t_13, t_23, t_33,
                _11, _12, _13, _21, _22, _23, _31, _32, _33)

        #fm_p_sq = (field**2)*abs(0.5*(th_11*th_11.conjugate()+th_22*th_22.conjugate())+th_12*th_12.conjugate())
        #fm_p_field = field*0.5*(th_11+th_22) 
        fm_p_sq = abs(0.5*(th_11*th_11.conjugate()+th_22*th_22.conjugate())+th_12*th_12.conjugate())
        fm_p_field = 0.5*(th_11+th_22) 
        cross = 2.*(f_nucl.real*fm_p_field.real+f_nucl.imag*fm_p_field.imag)

        iint_u = abs(f_nucl*f_nucl.conjugate()) + fm_p_sq + p_u*cross
        iint_d = abs(f_nucl*f_nucl.conjugate()) + fm_p_sq - p_d*cross

        #d_info_out = {"iint_u": iint_u, "iint_d": iint_d}   
        #d_info_out.update(d_info_cry)
        
        return iint_u, iint_d, refln

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
    

    def calc_shape_profile(self, tth, tth_hkl, i_g=0.):
        """
        calculate profile in the range ttheta for reflections placed on 
        ttheta_hkl with i_g parameter by default equal to zero
        
        tth, tth_hkl in degrees

        """
        zero_shift = float(self.offset)
        tth_zs = tth-zero_shift

        resolution = self.resolution
        
        h_pv, eta, h_g, h_l, a_g, b_g, a_l, b_l = resolution.calc_resolution(
                                                  tth_hkl, i_g=i_g)
    
        tth_2d, tth_hkl_2d = numpy.meshgrid(tth_zs, tth_hkl, indexing="ij")

        self.__ag = numpy.meshgrid(tth_zs, a_g, indexing="ij")[1]
        self.__bg = numpy.meshgrid(tth_zs, b_g, indexing="ij")[1]
        self.__al = numpy.meshgrid(tth_zs, a_l, indexing="ij")[1]
        self.__bl = numpy.meshgrid(tth_zs, b_l, indexing="ij")[1]
        eta_2d = numpy.meshgrid(tth_zs, eta, indexing="ij")[1]
        self.__eta = eta_2d 

        self._gauss_pd(tth_2d-tth_hkl_2d)
        self._lor_pd(tth_2d-tth_hkl_2d)
        g_pd_2d = self.__gauss_pd 
        l_pd_2d = self.__lor_pd
        
        np_shape_2d = eta_2d * l_pd_2d + (1.-eta_2d) * g_pd_2d

        asymmetry = self.asymmetry
        np_ass_2d = asymmetry.calc_asymmetry(tth_zs, tth_hkl)

        #Lorentz factor
        tth_rad = tth_zs*numpy.pi/180.
        np_lor_1d = 1./(numpy.sin(tth_rad)*numpy.sin(0.5*tth_rad))
        
        np_lor_2d = numpy.meshgrid(np_lor_1d, tth_hkl, indexing="ij")[0]
        
        profile_2d = np_shape_2d*np_ass_2d*np_lor_2d
        
        return profile_2d, tth_zs, h_pv


