
"""
define classe Pd which describes the 1d powder diffraction experiment
"""
__author__ = 'ikibalin'
__version__ = "2019_09_09"
import os
import numpy
from pycifstar import Data, Loop

from cryspy.f_common.cl_fitable import Fitable
from cryspy.f_experiment.cl_beam_polarization import BeamPolarization
from .cl_pd_background import PdBackground
from .cl_pd_exclude_2theta import PdExclude2Theta
from .cl_pd_instr_reflex_asymmetry import PdInstrReflexAsymmetry
from .cl_pd_instr_resolution import PdInstrResolution
from .cl_pd_meas import PdMeas
from .cl_pd_phase import PdPhase
from .cl_pd_proc import PdProc
from .cl_pd_peak import PdPeak


from cryspy.f_crystal.cl_magnetism import calc_mRmCmRT

class Pd(object):
    """
    Class to describe information about single diffraction measurements

    Example:

    
    _diffrn_radiation_wavelength 1.40 
    _diffrn_ambient_field  1.000

    _pd_chi2_sum True
    _pd_chi2_diff False
    _pd_chi2_up False
    _pd_chi2_down False

    _pd_calib_2theta_offset -0.385404

    _pd_2theta_range_min 4.0
    _pd_2theta_range_max 80.0
    """
    def __init__(self, label="powder1", 
                 background = None, exclude_2theta=PdExclude2Theta(min=[0.], max=[1.0]), asymmetry=PdInstrReflexAsymmetry(),
                 beam_polarization = BeamPolarization(), resolution = PdInstrResolution(), meas=PdMeas(), phase=PdPhase(),
                 wavelength=1.4, field=1.0,
                 chi2_sum = True, chi2_diff = True, chi2_up = False, chi2_down = False,
                 offset = Fitable(0.), range_min = 0.1, range_max = 179.9):
        super(Pd, self).__init__()
        self.__label = None
        self.__background = None
        self.__exclude_2theta = None
        self.__asymmetry = None
        self.__beam_polarization = None
        self.__resolution = None
        self.__meas = None 
        self.__phase = None 
        self.__proc = None 
        self.__diffrn_radiation_wavelength = None 
        self.__diffrn_ambient_field = None 
        self.__pd_chi2_sum = None 
        self.__pd_chi2_diff = None 
        self.__pd_chi2_up = None 
        self.__pd_chi2_down = None 
        self.__pd_calib_2theta_offset = None 
        self.__pd_2theta_range_min = None 
        self.__pd_2theta_range_max = None 

        self.label = label
        self.background = background
        self.exclude_2theta = exclude_2theta
        self.asymmetry = asymmetry
        self.beam_polarization = beam_polarization
        self.resolution = resolution
        self.meas = meas
        self.phase = phase
        self.wavelength = wavelength
        self.field = field
        self.chi2_sum = chi2_sum
        self.chi2_diff = chi2_diff
        self.chi2_up = chi2_up
        self.chi2_down = chi2_down
        self.offset = offset
        self.range_max = range_max
        self.range_min = range_min

        self.__proc = None
        self.__peaks = None
        self.__reflns = None

    @property
    def erase_precalc(self):
        """
        Erase precalculated values
        """
        self.__proc = None
        self.__peaks = None
        self.__reflns = None

    @property
    def label(self):
        return self.__label
    @label.setter
    def label(self, x: str):
        self.__label = str(x)

    @property
    def background(self):
        return self.__background
    @background.setter
    def background(self, x):
        if isinstance(x, PdBackground):
            x_in = x
        elif isinstance(x, str):
            x_in = PdBackground()
            flag = x_in.from_cif(x)
            if not(flag):
                x_in = None
                self._show_message("A induced string can not be converted to PdBackground")
        elif x is None:
            x_in = None
        else:
            x_in = None
            self._show_message("Input type for PdBackground is not recognized")
        self.__background = x_in

    @property
    def exclude_2theta(self):
        return self.__exclude_2theta
    @exclude_2theta.setter
    def exclude_2theta(self, x):
        if isinstance(x, PdExclude2Theta):
            x_in = x
        elif isinstance(x, str):
            x_in = PdExclude2Theta()
            flag = x_in.from_cif(x)
            if not(flag):
                x_in = None
                self._show_message("A induced string can not be converted to PdExclude2Theta")
        elif x is None:
            x_in = None
        else:
            x_in = None
            self._show_message("Input type for PdExclude2Theta is not recognized")
        self.__exclude_2theta = x_in
        self.erase_precalc

    @property
    def asymmetry(self):
        return self.__asymmetry
    @asymmetry.setter
    def asymmetry(self, x):
        if isinstance(x, PdInstrReflexAsymmetry):
            x_in = x
        elif isinstance(x, str):
            x_in = PdInstrReflexAsymmetry()
            flag = x_in.from_cif(x)
            if not(flag):
                x_in = None
                self._show_message("A induced string can not be converted to PdAsymmetry")
        elif x is None:
            x_in = None
        else:
            x_in = None
            self._show_message("Input type for PdAsymmetry is not recognized")
        self.__asymmetry = x_in

    @property
    def beam_polarization(self):
        return self.__beam_polarization
    @beam_polarization.setter
    def beam_polarization(self, x):
        if isinstance(x, BeamPolarization):
            x_in = x
        elif isinstance(x, str):
            x_in = BeamPolarization()
            flag = x_in.from_cif(x)
            if not(flag):
                x_in = None
                self._show_message("A induced string can not be converted to BeamPolarization")
        elif x is None:
            x_in = None
        else:
            x_in = None
            self._show_message("Input type for BeamPolarization is not recognized")
        self.__beam_polarization = x_in

    @property
    def resolution(self):
        return self.__resolution
    @resolution.setter
    def resolution(self, x):
        if isinstance(x, PdInstrResolution):
            x_in = x
        elif isinstance(x, str):
            x_in = PdInstrResolution()
            flag = x_in.from_cif(x)
            if not(flag):
                x_in = None
                self._show_message("A induced string can not be converted to PdResolution")
        elif x is None:
            x_in = None
        else:
            x_in = None
            self._show_message("Input type for PdResolution is not recognized")
        self.__resolution = x_in

    @property
    def meas(self):
        return self.__meas
    @meas.setter
    def meas(self, x):
        if isinstance(x, PdMeas):
            x_in = x
        elif isinstance(x, str):
            x_in = PdMeas()
            flag = x_in.from_cif(x)
            if not(flag):
                x_in = None
                self._show_message("A induced string can not be converted to PdMeas")
        elif x is None:
            x_in = None
        else:
            x_in = None
            self._show_message("Input type for PdMeas is not recognized")
        self.__meas = x_in
        self.erase_precalc

    @property
    def phase(self):
        return self.__phase
    @phase.setter
    def phase(self, x):
        if isinstance(x, PdPhase):
            x_in = x
        elif isinstance(x, str):
            x_in = PdPhase()
            flag = x_in.from_cif(x)
            if not(flag):
                x_in = None
                self._show_message("A induced string can not be converted to PdPhase")
        elif x is None:
            x_in = None
        else:
            x_in = None
            self._show_message("Input type for PdPhase is not recognized")
        self.__phase = x_in
        self.erase_precalc


    @property
    def proc(self):
        return self.__proc
    @proc.setter
    def proc(self, x):
        if isinstance(x, PdProc):
            x_in = x
        elif isinstance(x, str):
            x_in = PdProc()
            flag = x_in.from_cif(x)
            if not(flag):
                x_in = None
                self._show_message("A induced string can not be converted to PdProc")
        elif x is None:
            x_in = None
        else:
            x_in = None
            self._show_message("Input type for PdProc is not recognized")
        self.__proc = x_in


    @property
    def wavelength(self):
        return self.__diffrn_radiation_wavelength
    @wavelength.setter
    def wavelength(self, x):
        if isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        self.__diffrn_radiation_wavelength = x_in
        self.erase_precalc

    @property
    def field(self):
        return self.__diffrn_ambient_field
    @field.setter
    def field(self, x):
        if isinstance(x, float):
            x_in = x
        elif x is None:
            x_in = None
        else:
            x_in = float(x)
        self.__diffrn_ambient_field = x_in
        self.erase_precalc

    @property
    def range_min(self):
        return self.__pd_2theta_range_min
    @range_min.setter
    def range_min(self, x):
        if isinstance(x, float):
            x_in = x
        elif x is None:
            x_in = None
        else:
            x_in = float(x)
        self.__pd_2theta_range_min = x_in
        self.erase_precalc

    @property
    def range_max(self):
        return self.__pd_2theta_range_max
    @range_max.setter
    def range_max(self, x):
        if isinstance(x, float):
            x_in = x
        elif x is None:
            x_in = None
        else:
            x_in = float(x)
        self.__pd_2theta_range_max = x_in
        self.erase_precalc

    @property
    def chi2_sum(self):
        return self.__pd_chi2_sum
    @chi2_sum.setter
    def chi2_sum(self, x):
        if isinstance(x, bool):
            x_in = x
        elif x is None:
            x_in = None
        else:
            x_in = bool(x)
        self.__pd_chi2_sum = x_in

    @property
    def chi2_diff(self):
        return self.__pd_chi2_diff
    @chi2_diff.setter
    def chi2_diff(self, x):
        if isinstance(x, bool):
            x_in = x
        elif x is None:
            x_in = None
        else:
            x_in = bool(x)
        self.__pd_chi2_diff = x_in

    @property
    def chi2_up(self):
        return self.__pd_chi2_up
    @chi2_up.setter
    def chi2_up(self, x):
        if isinstance(x, bool):
            x_in = x
        elif x is None:
            x_in = None
        else:
            x_in = bool(x)
        self.__pd_chi2_up = x_in

    @property
    def chi2_down(self):
        return self.__pd_chi2_down
    @chi2_down.setter
    def chi2_down(self, x):
        if isinstance(x, bool):
            x_in = x
        elif x is None:
            x_in = None
        else:
            x_in = bool(x)
        self.__pd_chi2_down = x_in

    @property
    def offset(self):
        return self.__pd_calib_2theta_offset
    @offset.setter
    def offset(self, x):
        if isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        self.__pd_calib_2theta_offset = x_in


    @property
    def peaks(self):
        return self.__peaks
    @peaks.setter
    def peaks(self, x):
        self.__peaks = x

    @property
    def reflns(self):
        return self.__reflns
    @reflns.setter
    def reflns(self, x):
        self.__reflns = x

    def _show_message(self, s_out: str):
        print("***  Error ***")
        print(s_out)

    def __repr__(self):
        ls_out = ["Pd:"]
        if self.label is not None:
            ls_out.append(" label: {:}".format(self.label))
        if self.wavelength is not None:
            ls_out.append(" wavelength: {:.3f}".format(float(self.wavelength)))
        if self.field is not None:
            ls_out.append(" field: {:.3f}".format(float(self.field)))
        if self.range_min is not None:
            ls_out.append(" range_min: {:.3f}".format(float(self.range_min)))
        if self.range_max is not None:
            ls_out.append(" range_max: {:.3f}".format(float(self.range_max)))
        if self.offset is not None:
            ls_out.append(" offset: {:}".format(self.offset.print_with_sigma))
        if self.chi2_up is not None:
            ls_out.append(" chi2_up: {:}".format(self.chi2_up))
        if self.chi2_down is not None:
            ls_out.append(" chi2_down: {:}".format(self.chi2_down))
        if self.chi2_sum is not None:
            ls_out.append(" chi2_sum: {:}".format(self.chi2_sum))
        if self.chi2_diff is not None:
            ls_out.append(" chi2_diff: {:}".format(self.chi2_diff))

        if self.background is not None:
            ls_out.append("\n"+str(self.background))
        if self.exclude_2theta is not None:
            ls_out.append("\n"+str(self.exclude_2theta))
        if self.asymmetry is not None:
            ls_out.append("\n"+str(self.asymmetry))
        if self.beam_polarization is not None:
            ls_out.append("\n"+str(self.beam_polarization))
        if self.resolution is not None:
            ls_out.append("\n"+str(self.resolution))
        if self.phase is not None:
            ls_out.append("\n"+str(self.phase))
        if self.reflns is not None:
            ls_out.extend(["\n"+str(_) for _ in self.reflns])
        if self.peaks is not None:
            ls_out.extend(["\n"+str(_) for _ in self.peaks])
        if self.proc is not None:
            ls_out.append("\n"+str(self.proc))
        if self.meas is not None:
            ls_out.append("\n"+str(self.meas))
        return "\n".join(ls_out)

    
    @property
    def is_variable(self):
        l_bool = []
        if self.background is not None:
            l_bool.append(self.background.is_variable)
        if self.exclude_2theta is not None:
            l_bool.append(self.exclude_2theta.is_variable)
        if self.asymmetry is not None:
            l_bool.append(self.asymmetry.is_variable)
        if self.beam_polarization is not None:
            l_bool.append(self.beam_polarization.is_variable)
        if self.resolution is not None:
            l_bool.append(self.resolution.is_variable)
        if self.phase is not None:
            l_bool.append(self.phase.is_variable)
        l_bool.append(self.offset.refinement)
        if self.wavelength.refinement:
            l_bool.append(self.wavelength)
        res = any(l_bool)
        return res

    def get_variables(self):
        l_variable = []
        if self.background is not None:
            l_variable.extend(self.background.get_variables())
        if self.exclude_2theta is not None:
            l_variable.extend(self.exclude_2theta.get_variables())
        if self.asymmetry is not None:
            l_variable.extend(self.asymmetry.get_variables())
        if self.beam_polarization is not None:
            l_variable.extend(self.beam_polarization.get_variables())
        if self.resolution is not None:
            l_variable.extend(self.resolution.get_variables())
        if self.phase is not None:
            l_variable.extend(self.phase.get_variables())
        if self.offset.refinement: l_variable.append(self.offset)
        if self.wavelength.refinement: l_variable.append(self.wavelength)
        return l_variable


    @property
    def to_cif(self):
        ls_out = []
        ls_out.append("data_{:}\n".format(self.label))
        str_1 = self.params_to_cif
        if str_1 != "":
            ls_out.append(str_1)
        str_2 = self.data_to_cif
        if str_2 != "":
            ls_out.append(str_2)
        str_3 = self.calc_to_cif
        if str_3 != "":
            ls_out.append(str_3)
        return "\n".join(ls_out)
    
    @property
    def params_to_cif(self):
        ls_out = []
        if self.wavelength is not None:
            ls_out.append("_diffrn_radiation_wavelength {:}".format(self.wavelength.print_with_sigma))
        if self.field is not None:
            ls_out.append("_diffrn_ambient_field {:}".format(self.field))
        if self.chi2_sum is not None:
            if self.chi2_sum:
                s_line = "True" 
            else: 
                s_line = "False"
            ls_out.append("\n_pd_chi2_sum {:}".format(s_line))
        if self.chi2_diff is not None:
            if self.chi2_diff:
                s_line = "True" 
            else: 
                s_line = "False"
            ls_out.append("_pd_chi2_diff {:}".format(s_line))
        if self.chi2_up is not None:
            if self.chi2_up:
                s_line = "True" 
            else: 
                s_line = "False"
            ls_out.append("_pd_chi2_up {:}".format(s_line))
        if self.chi2_down is not None:
            if self.chi2_down:
                s_line = "True" 
            else: 
                s_line = "False"
            ls_out.append("_pd_chi2_down {:}".format(s_line))
        if self.range_min is not None:
            ls_out.append("\n_pd_2theta_range_min {:.3f}".format( self.range_min))
        if self.range_max is not None:
            ls_out.append("_pd_2theta_range_max {:.3f}".format( self.range_max))
        if self.offset is not None:
            ls_out.append("_pd_calib_2theta_offset {:}".format(self.offset.print_with_sigma))
        if self.background is not None:
            ls_out.append("\n"+self.background.to_cif)
        if self.exclude_2theta is not None:
            ls_out.append("\n"+self.exclude_2theta.to_cif)
        if self.asymmetry is not None:
            ls_out.append("\n"+self.asymmetry.to_cif)
        if self.beam_polarization is not None:
            ls_out.append("\n"+self.beam_polarization.to_cif)
        if self.resolution is not None:
            ls_out.append("\n"+self.resolution.to_cif)
        if self.phase is not None:
            ls_out.append("\n"+self.phase.to_cif)
        return "\n".join(ls_out)

    @property
    def data_to_cif(self):
        ls_out = []
        if self.meas is not None:
            ls_out.append("\n"+self.meas.to_cif)
        return "\n".join(ls_out)

    @property
    def calc_to_cif(self):
        ls_out = []
        if self.reflns is not None:
            ls_out.extend(["\n"+_.to_cif for _ in self.reflns])
        if self.peaks is not None:
            ls_out.extend(["\n"+_.to_cif for _ in self.peaks])
        if self.proc is not None:
            ls_out.append("\n"+self.proc.to_cif)
        return "\n".join(ls_out)

    def from_cif(self, string: str):
        cif_data = Data()
        flag = cif_data.take_from_string(string)
        if not flag:
            return False
        self.label = cif_data.name
        cif_values = cif_data.items
        if cif_values is not None:
            if cif_values.is_prefix("_diffrn_radiation_wavelength"):
                self.wavelength = cif_values["_diffrn_radiation_wavelength"]
            if cif_values.is_prefix("_diffrn_ambient_field"):
                self.field = float(cif_values["_diffrn_ambient_field"])
            if cif_values.is_prefix("_pd_chi2_sum"):
                self.chi2_sum = (cif_values["_pd_chi2_sum"].value).strip().lower()=="true"
            if cif_values.is_prefix("_pd_chi2_diff"):
                self.chi2_diff = (cif_values["_pd_chi2_diff"].value).strip().lower()=="true"
            if cif_values.is_prefix("_pd_chi2_sum"):
                self.chi2_up = (cif_values["_pd_chi2_up"].value).strip().lower()=="true"
            if cif_values.is_prefix("_pd_chi2_down"):
                self.chi2_down = (cif_values["_pd_chi2_down"].value).strip().lower()=="true"
            if cif_values.is_prefix("_pd_calib_2theta_offset"):
                self.offset = cif_values["_pd_calib_2theta_offset"]
            if cif_values.is_prefix("_pd_2theta_range_min"):
                self.range_min = float(cif_values["_pd_2theta_range_min"])
            if cif_values.is_prefix("_pd_2theta_range_max"):
                self.range_max = float(cif_values["_pd_2theta_range_max"])
            if cif_values.is_prefix("_diffrn_radiation"):
                self.beam_polarization = str(cif_values)
            if cif_values.is_prefix("_pd_instr_reflex_asymmetry"): 
                self.asymmetry = "\n".join([str(_) for _ in cif_data["_pd_instr_reflex_asymmetry"]])
            if cif_values.is_prefix("_pd_instr_resolution"): 
                self.resolution = "\n".join([str(_) for _ in cif_data["_pd_instr_resolution"]])

        if cif_data.is_prefix("_pd_exclude_2theta"): self.exclude_2theta = str(cif_data["_pd_exclude_2theta"])
        if cif_data.is_prefix("_pd_phase"): self.phase = str(cif_data["_pd_phase"])
        if cif_data.is_prefix("_pd_background"): self.background = str(cif_data["_pd_background"])
        if cif_data.is_prefix("_pd_meas"): self.meas = str(cif_data["_pd_meas"])

        return True

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

        cell = crystal.cell
        
        #k_loc = cell.calc_k_loc(h, k, l)
        t_11, t_12, t_13, t_21, t_22, t_23, t_31, t_32, t_33 = cell.calc_m_t(h, k, l)
        
        th_11, th_12, th_13, th_21, th_22, th_23, th_31, th_32, th_33 = calc_mRmCmRT(
                t_11, t_21, t_31, t_12, t_22, t_32, t_13, t_23, t_33,
                sft_11, sft_12, sft_13, sft_21, sft_22, sft_23, sft_31, sft_32, 
                sft_33)

        fm_p_sq = (field**2)*abs(0.5*(th_11*th_11.conjugate()+th_22*th_22.conjugate())+th_12*th_12.conjugate())
        fm_p_field = field*0.5*(th_11+th_22) 
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


