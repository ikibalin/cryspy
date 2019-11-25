
"""
define classe Pd2d which describes the 2d powder diffraction experiment
"""
__author__ = 'ikibalin'
__version__ = "2019_09_11"
import os
import numpy
from pycifstar import Data, Loop

from cryspy.f_common.cl_fitable import Fitable
from cryspy.f_experiment.cl_beam_polarization import BeamPolarization
from .cl_pd2d_background import Pd2dBackground
from .cl_pd2d_exclude import Pd2dExclude
from .cl_pd2d_instr_reflex_asymmetry import Pd2dInstrReflexAsymmetry
from .cl_pd2d_instr_resolution import Pd2dInstrResolution
from .cl_pd2d_meas import Pd2dMeas
from .cl_pd2d_phase import Pd2dPhase
from .cl_pd2d_proc import Pd2dProc
from .cl_pd2d_peak import Pd2dPeak


from cryspy.f_crystal.cl_magnetism import calc_mRmCmRT

class Pd2d(object):
    """
    Class to describe information about single diffraction measurements

    Example:

    
    _diffrn_radiation_wavelength 1.40 
    _diffrn_ambient_field  1.000

    _pd2d_chi2_sum True
    _pd2d_chi2_diff False
    _pd2d_chi2_up False
    _pd2d_chi2_down False

    _pd2d_calib_2theta_offset -0.385404

    _pd2d_2theta_range_min 4.0
    _pd2d_2theta_range_max 80.0
    _pd2d_phi_range_min -2.0
    _pd2d_phi_range_max 40.0
    """
    def __init__(self, label="powder2", 
                 background = None, exclude=Pd2dExclude(ttheta_min=[0.], ttheta_max=[1.], phi_min=[0.0], phi_max=[0.1]), 
                 asymmetry=Pd2dInstrReflexAsymmetry(),
                 beam_polarization = BeamPolarization(), resolution = Pd2dInstrResolution(), meas=Pd2dMeas(), phase=Pd2dPhase(),
                 wavelength=1.4, field=1.0,
                 chi2_sum = True, chi2_diff = True, chi2_up = False, chi2_down = False,
                 ttheta_offset = Fitable(0.), ttheta_min = 0.1, ttheta_max = 179.9, phi_min=0., phi_max=40.):
        super(Pd2d, self).__init__()
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
        self.__pd2d_chi2_sum = None 
        self.__pd2d_chi2_diff = None 
        self.__pd2d_chi2_up = None 
        self.__pd2d_chi2_down = None 
        self.__pd2d_calib_2theta_offset = None 
        self.__pd2d_2theta_range_min = None 
        self.__pd2d_2theta_range_max = None 
        self.__pd2d_phi_range_min = None 
        self.__pd2d_phi_range_max = None 

        self.__dd = None 

        self.label = label
        self.background = background
        self.exclude = exclude
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
        self.ttheta_offset = ttheta_offset
        self.ttheta_max = ttheta_max
        self.ttheta_min = ttheta_min
        self.phi_max = phi_max
        self.phi_min = phi_min

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
        if isinstance(x, Pd2dBackground):
            x_in = x
        elif isinstance(x, str):
            x_in = Pd2dBackground()
            flag = x_in.from_cif(x)
            if not(flag):
                x_in = None
                self._show_message("A induced string can not be converted to Pd2dBackground")
        elif x is None:
            x_in = None
        else:
            x_in = None
            self._show_message("Input type for Pd2dBackground is not recognized")
        self.__background = x_in

    @property
    def exclude(self):
        return self.__exclude_2theta
    @exclude.setter
    def exclude(self, x):
        if isinstance(x, Pd2dExclude):
            x_in = x
        elif isinstance(x, str):
            x_in = Pd2dExclude()
            flag = x_in.from_cif(x)
            if not(flag):
                x_in = None
                self._show_message("A induced string can not be converted to Pd2dExclude")
        elif x is None:
            x_in = None
        else:
            x_in = None
            self._show_message("Input type for Pd2dExclude is not recognized")
        self.__exclude_2theta = x_in
        self.erase_precalc

    @property
    def asymmetry(self):
        return self.__asymmetry
    @asymmetry.setter
    def asymmetry(self, x):
        if isinstance(x, Pd2dInstrReflexAsymmetry):
            x_in = x
        elif isinstance(x, str):
            x_in = Pd2dInstrReflexAsymmetry()
            flag = x_in.from_cif(x)
            if not(flag):
                x_in = None
                self._show_message("A induced string can not be converted to Pd2dAsymmetry")
        elif x is None:
            x_in = None
        else:
            x_in = None
            self._show_message("Input type for Pd2dAsymmetry is not recognized")
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
        if isinstance(x, Pd2dInstrResolution):
            x_in = x
        elif isinstance(x, str):
            x_in = Pd2dInstrResolution()
            flag = x_in.from_cif(x)
            if not(flag):
                x_in = None
                self._show_message("A induced string can not be converted to Pd2dResolution")
        elif x is None:
            x_in = None
        else:
            x_in = None
            self._show_message("Input type for Pd2dResolution is not recognized")
        self.__resolution = x_in

    @property
    def meas(self):
        return self.__meas
    @meas.setter
    def meas(self, x):
        if isinstance(x, Pd2dMeas):
            x_in = x
        elif isinstance(x, str):
            x_in = Pd2dMeas()
            flag = x_in.from_cif(x)
            if not(flag):
                x_in = None
                self._show_message("A induced string can not be converted to Pd2dMeas")
        elif x is None:
            x_in = None
        else:
            x_in = None
            self._show_message("Input type for Pd2dMeas is not recognized")
        self.__meas = x_in
        self.erase_precalc


    @property
    def phase(self):
        return self.__phase
    @phase.setter
    def phase(self, x):
        if isinstance(x, Pd2dPhase):
            x_in = x
        elif isinstance(x, str):
            x_in = Pd2dPhase()
            flag = x_in.from_cif(x)
            if not(flag):
                x_in = None
                self._show_message("A induced string can not be converted to Pd2dPhase")
        elif x is None:
            x_in = None
        else:
            x_in = None
            self._show_message("Input type for Pd2dPhase is not recognized")
        self.__phase = x_in
        self.erase_precalc



    @property
    def proc(self):
        return self.__proc
    @proc.setter
    def proc(self, x):
        if isinstance(x, Pd2dProc):
            x_in = x
        elif isinstance(x, str):
            x_in = Pd2dProc()
            flag = x_in.from_cif(x)
            if not(flag):
                x_in = None
                self._show_message("A induced string can not be converted to Pd2dProc")
        elif x is None:
            x_in = None
        else:
            x_in = None
            self._show_message("Input type for Pd2dProc is not recognized")
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
    def ttheta_min(self):
        return self.__pd2d_2theta_range_min
    @ttheta_min.setter
    def ttheta_min(self, x):
        if isinstance(x, float):
            x_in = x
        elif x is None:
            x_in = None
        else:
            x_in = float(x)
        self.__pd2d_2theta_range_min = x_in
        self.erase_precalc


    @property
    def phi_min(self):
        return self.__pd2d_phi_range_min
    @phi_min.setter
    def phi_min(self, x):
        if isinstance(x, float):
            x_in = x
        elif x is None:
            x_in = None
        else:
            x_in = float(x)
        self.__pd2d_phi_range_min = x_in
        self.erase_precalc



    @property
    def ttheta_max(self):
        return self.__pd2d_2theta_range_max
    @ttheta_max.setter
    def ttheta_max(self, x):
        if isinstance(x, float):
            x_in = x
        elif x is None:
            x_in = None
        else:
            x_in = float(x)
        self.__pd2d_2theta_range_max = x_in
        self.erase_precalc



    @property
    def phi_max(self):
        return self.__pd2d_phi_range_max
    @phi_max.setter
    def phi_max(self, x):
        if isinstance(x, float):
            x_in = x
        elif x is None:
            x_in = None
        else:
            x_in = float(x)
        self.__pd2d_phi_range_max = x_in
        self.erase_precalc


    @property
    def chi2_sum(self):
        return self.__pd2d_chi2_sum
    @chi2_sum.setter
    def chi2_sum(self, x):
        if isinstance(x, bool):
            x_in = x
        elif x is None:
            x_in = None
        else:
            x_in = bool(x)
        self.__pd2d_chi2_sum = x_in

    @property
    def chi2_diff(self):
        return self.__pd2d_chi2_diff
    @chi2_diff.setter
    def chi2_diff(self, x):
        if isinstance(x, bool):
            x_in = x
        elif x is None:
            x_in = None
        else:
            x_in = bool(x)
        self.__pd2d_chi2_diff = x_in

    @property
    def chi2_up(self):
        return self.__pd2d_chi2_up
    @chi2_up.setter
    def chi2_up(self, x):
        if isinstance(x, bool):
            x_in = x
        elif x is None:
            x_in = None
        else:
            x_in = bool(x)
        self.__pd2d_chi2_up = x_in

    @property
    def chi2_down(self):
        return self.__pd2d_chi2_down
    @chi2_down.setter
    def chi2_down(self, x):
        if isinstance(x, bool):
            x_in = x
        elif x is None:
            x_in = None
        else:
            x_in = bool(x)
        self.__pd2d_chi2_down = x_in

    @property
    def ttheta_offset(self):
        return self.__pd2d_calib_2theta_offset
    @ttheta_offset.setter
    def ttheta_offset(self, x):
        if isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        self.__pd2d_calib_2theta_offset = x_in

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
        ls_out = ["Pd2d:"]
        if self.label is not None:
            ls_out.append(" label: {:}".format(self.label))
        if self.wavelength is not None:
            ls_out.append(" wavelength: {:.3f}".format(float(self.wavelength)))
        if self.field is not None:
            ls_out.append(" field: {:.3f}".format(float(self.field)))
        if self.ttheta_min is not None:
            ls_out.append(" ttheta_min: {:.3f}".format(float(self.ttheta_min)))
        if self.ttheta_max is not None:
            ls_out.append(" ttheta_max: {:.3f}".format(float(self.ttheta_max)))
        if self.phi_min is not None:
            ls_out.append(" phi_min: {:.3f}".format(float(self.phi_min)))
        if self.phi_max is not None:
            ls_out.append(" phi_max: {:.3f}".format(float(self.phi_max)))
        if self.ttheta_offset is not None:
            ls_out.append(" ttheta_offset: {:}".format(self.ttheta_offset.print_with_sigma))
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
        if self.exclude is not None:
            ls_out.append("\n"+str(self.exclude))
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
        if self.exclude is not None:
            l_bool.append(self.exclude.is_variable)
        if self.asymmetry is not None:
            l_bool.append(self.asymmetry.is_variable)
        if self.beam_polarization is not None:
            l_bool.append(self.beam_polarization.is_variable)
        if self.resolution is not None:
            l_bool.append(self.resolution.is_variable)
        if self.phase is not None:
            l_bool.append(self.phase.is_variable)
        l_bool.append(self.is_variable_offset)
        if self.wavelength.refinement:
            l_bool.append(self.wavelength)
        res = any(l_bool)
        return res

    @property
    def is_variable_offset(self):
        l_bool = []
        l_bool.append(self.ttheta_offset.refinement)
        res = any(l_bool)
        return res

    def get_variables(self):
        l_variable = []
        if self.background is not None:
            l_variable.extend(self.background.get_variables())
        if self.exclude is not None:
            l_variable.extend(self.exclude.get_variables())
        if self.asymmetry is not None:
            l_variable.extend(self.asymmetry.get_variables())
        if self.beam_polarization is not None:
            l_variable.extend(self.beam_polarization.get_variables())
        if self.resolution is not None:
            l_variable.extend(self.resolution.get_variables())
        if self.phase is not None:
            l_variable.extend(self.phase.get_variables())
        if self.ttheta_offset.refinement: l_variable.append(self.ttheta_offset)
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
            ls_out.append("\n_pd2d_chi2_sum {:}".format(s_line))
        if self.chi2_diff is not None:
            if self.chi2_diff:
                s_line = "True" 
            else: 
                s_line = "False"
            ls_out.append("_pd2d_chi2_diff {:}".format(s_line))
        if self.chi2_up is not None:
            if self.chi2_up:
                s_line = "True" 
            else: 
                s_line = "False"
            ls_out.append("_pd2d_chi2_up {:}".format(s_line))
        if self.chi2_down is not None:
            if self.chi2_down:
                s_line = "True" 
            else: 
                s_line = "False"
            ls_out.append("_pd2d_chi2_down {:}".format(s_line))
        if self.ttheta_min is not None:
            ls_out.append("\n_pd2d_2theta_range_min {:.3f}".format( self.ttheta_min))
        if self.ttheta_max is not None:
            ls_out.append("_pd2d_2theta_range_max {:.3f}".format( self.ttheta_max))
        if self.phi_min is not None:
            ls_out.append("\n_pd2d_phi_range_min {:.3f}".format( self.phi_min))
        if self.ttheta_max is not None:
            ls_out.append("_pd2d_phi_range_max {:.3f}".format( self.phi_max))
        if self.ttheta_offset is not None:
            ls_out.append("_pd2d_calib_2theta_offset {:}".format(self.ttheta_offset.print_with_sigma))
        if self.background is not None:
            ls_out.append("\n"+self.background.to_cif)
        if self.exclude is not None:
            ls_out.append("\n"+self.exclude.to_cif)
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
            if cif_values.is_prefix("_pd2d_chi2_sum"):
                self.chi2_sum = (cif_values["_pd2d_chi2_sum"].value).strip().lower()=="true"
            if cif_values.is_prefix("_pd2d_chi2_diff"):
                self.chi2_diff = (cif_values["_pd2d_chi2_diff"].value).strip().lower()=="true"
            if cif_values.is_prefix("_pd2d_chi2_sum"):
                self.chi2_up = (cif_values["_pd2d_chi2_up"].value).strip().lower()=="true"
            if cif_values.is_prefix("_pd2d_chi2_down"):
                self.chi2_down = (cif_values["_pd2d_chi2_down"].value).strip().lower()=="true"
            if cif_values.is_prefix("_pd2d_calib_2theta_offset"):
                self.ttheta_offset = cif_values["_pd2d_calib_2theta_offset"]
            if cif_values.is_prefix("_pd2d_2theta_range_min"):
                self.ttheta_min = float(cif_values["_pd2d_2theta_range_min"])
            if cif_values.is_prefix("_pd2d_2theta_range_max"):
                self.ttheta_max = float(cif_values["_pd2d_2theta_range_max"])
            if cif_values.is_prefix("_pd2d_phi_range_min"):
                self.phi_min = float(cif_values["_pd2d_phi_range_min"])
            if cif_values.is_prefix("_pd2d_phi_range_max"):
                self.phi_max = float(cif_values["_pd2d_phi_range_max"])
            if cif_values.is_prefix("_diffrn_radiation"):
                self.beam_polarization = str(cif_values)
            if cif_values.is_prefix("_pd2d_background_2theta_phi_intensity"): 
                self.background = str(cif_values["_pd2d_background_2theta_phi_intensity"])
            if cif_values.is_prefix("_pd2d_instr_reflex_asymmetry"): 
                self.asymmetry = "\n".join([str(_) for _ in cif_values["_pd2d_instr_reflex_asymmetry"]])
            if cif_values.is_prefix("_pd2d_instr_resolution"): 
                self.resolution = "\n".join([str(_) for _ in cif_values["_pd2d_instr_resolution"]])
            if cif_values.is_prefix("_pd2d_meas_2theta_phi_intensity"): 
                self.meas = "\n".join([str(_) for _ in cif_values["_pd2d_meas_2theta_phi_intensity"]])
        if cif_data.is_prefix("_pd2d_exclude"): self.exclude = str(cif_data["_pd2d_exclude"])
        if cif_data.is_prefix("_pd2d_phase"): self.phase = str(cif_data["_pd2d_phase"])

        return True

    def calc_profile(self, tth, phi, l_crystal, l_peak_in=[], l_refln_in=[], l_dd_in = []):
        """
        calculate intensity for the given diffraction angle
        """
        proc = Pd2dProc() #it is output
        proc.ttheta = tth
        proc.phi = phi
        
        background = self.background
        int_bkgd = background.interpolate_by_points(tth, phi)
        proc.bkg_calc = int_bkgd
        tth_rad = tth*numpy.pi/180.
        phi_rad = phi*numpy.pi/180.
        cos_theta_1d = numpy.cos(0.5*tth_rad)
        sin_phi_1d = numpy.sin(phi_rad)

        wavelength = float(self.wavelength)
        beam_polarization = self.beam_polarization

        p_u = float(beam_polarization.polarization)
        p_d = (2.*float(beam_polarization.efficiency)-1.)*p_u
        

        tth_min = tth.min()
        tth_max = tth.max()+3. 
        if tth_max > 180.:
            tth_max = 180.
        sthovl_min = numpy.sin(0.5*tth_min*numpy.pi/180.)/wavelength
        sthovl_max = numpy.sin(0.5*tth_max*numpy.pi/180.)/wavelength

        res_u_2d = numpy.zeros((tth.shape[0],phi.shape[0]), dtype=float)
        res_d_2d = numpy.zeros((tth.shape[0],phi.shape[0]), dtype=float)


        phase = self.phase
        if len(phase.label) != len(l_peak_in):
            l_peak_in = len(phase.label) * [None] 
        if len(phase.label) != len(l_refln_in):
            l_refln_in = len(phase.label) * [None] 
        if len(phase.label) != len(l_dd_in):
            l_dd_in = len(phase.label) * [None] 
        


        l_peak, l_refln, l_dd_out = [], [], []
        for phase_label, phase_scale, phase_igsize, peak_in, refln_in, dd_in in zip(
            phase.label, phase.scale, phase.igsize, l_peak_in, l_refln_in, l_dd_in):
            dd_out = {}
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
            
            peak = Pd2dPeak()
            if peak_in is not None:
                h, k, l, mult = peak_in.h, peak_in.k, peak_in.l, peak_in.mult
            else:
                h, k, l, mult = cell.calc_hkl(space_group, sthovl_min, sthovl_max)
            peak.h, peak.k, peak.l, peak.mult = h, k, l, mult

                
            cond_1 = not(crystal.is_variable)
            cond_2 = (peak_in is not None) & (refln_in is not None)
            if cond_1 & cond_2:
                f_nucl_sq, f_m_p_sin_sq = peak_in.f_nucl_sq, peak_in.f_m_p_sin_sq
                f_m_p_cos_sq, cross_sin = peak_in.f_m_p_cos_sq, peak_in.cross_sin 
                refln = refln_in
            else:
                f_nucl_sq, f_m_p_sin_sq, f_m_p_cos_sq, cross_sin, refln = self.calc_for_iint(h, k, l, crystal)
            l_refln.append(refln)

            peak.f_nucl_sq, peak.f_m_p_sin_sq = f_nucl_sq, f_m_p_sin_sq
            peak.f_m_p_cos_sq, peak.cross_sin = f_m_p_cos_sq, cross_sin

            cond_1 = dd_in is not None
            cond_2 = ((not(crystal.is_variable)) & 
                      (not(beam_polarization.polarization.refinement)) & 
                      (not(beam_polarization.efficiency.refinement)))
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
            
            tth_hkl_rad = 2.*numpy.arcsin(sthovl_hkl*wavelength)
            tth_hkl = tth_hkl_rad*180./numpy.pi
            peak.ttheta = tth_hkl


            cond_1 = dd_in is not None
            cond_2 = ((not(phase_igsize.refinement)) & 
                      (not(self.resolution.is_variable)) & 
                      (not(self.is_variable_offset)) &
                      (not(self.wavelength.refinement)))
            if cond_1 & cond_2:
                profile_3d, tth_zs, h_pv = dd_in["profile_3d"], dd_in["tth_zs"], dd_in["h_pv"]
            else:
                profile_3d, tth_zs, h_pv = self.calc_shape_profile(tth, phi, tth_hkl, i_g)
            dd_out["profile_3d"] = profile_3d
            dd_out["tth_zs"] = tth_zs
            dd_out["h_pv"] = h_pv
            peak.width_2theta = h_pv
            
            
            res_u_3d = profile_3d*iint_u_3d 
            res_d_3d = profile_3d*iint_d_3d 

            res_u_2d += scale*res_u_3d.sum(axis=2) 
            res_d_2d += scale*res_d_3d.sum(axis=2) 
            l_peak.append(peak)
            l_dd_out.append(dd_out)
        proc.ttheta_corrected = tth_zs
        proc.up_net = res_u_2d
        proc.down_net = res_d_2d
        proc.up_total = res_u_2d+int_bkgd
        proc.down_total = res_d_2d+int_bkgd
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
        calculate chi square
        """
        meas = self.meas

        tth = meas.ttheta
        phi = meas.phi
        int_u_exp = meas.up
        sint_u_exp = meas.up_sigma
        int_d_exp = meas.down
        sint_d_exp = meas.down_sigma

        if ((self.peaks is not None) & (self.reflns is not None)):
            l_peak_in = self.peaks
            l_refln_in = self.reflns
            l_dd_in = self.__dd
        else:
            l_peak_in, l_refln_in, l_dd_in = [], [], [] 

        cond_tth_in = numpy.ones(tth.size, dtype=bool)
        cond_tth_in = numpy.logical_and(cond_tth_in, tth >= self.ttheta_min)
        cond_tth_in = numpy.logical_and(cond_tth_in, tth <= self.ttheta_max)

        cond_phi_in = numpy.ones(phi.size, dtype=bool)
        cond_phi_in = numpy.logical_and(cond_phi_in, phi >= self.phi_min)
        cond_phi_in = numpy.logical_and(cond_phi_in, phi <= self.phi_max)

        #cond_1_in, cond_2_in = numpy.meshgrid(cond_tth_in, cond_phi_in, indexing="ij")
        #cond_in = numpy.logical_and(cond_1_in, cond_2_in)
        tth_in = tth[cond_tth_in]
        phi_in = phi[cond_phi_in]
        int_u_exp_in = int_u_exp[cond_tth_in, :][:, cond_phi_in]
        sint_u_exp_in = sint_u_exp[cond_tth_in, :][:, cond_phi_in]
        int_d_exp_in = int_d_exp[cond_tth_in, :][:, cond_phi_in]
        sint_d_exp_in = sint_d_exp[cond_tth_in, :][:, cond_phi_in]



        proc, l_peak, l_refln, l_dd_out = self.calc_profile(tth_in, phi_in, l_crystal, l_peak_in, l_refln_in, l_dd_in)
        proc.up = int_u_exp_in
        proc.up_sigma = sint_u_exp_in
        proc.down = int_d_exp_in
        proc.down_sigma = sint_d_exp_in
        self.proc = proc
        self.peaks = l_peak
        self.reflns = l_refln
        self.__dd = l_dd_out

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


    def calc_for_iint(self, h, k, l, crystal):
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

        #f_m_p_sin_sq = (field**2)*abs(0.5*(th_11*th_11.conjugate()+th_22*th_22.conjugate())+th_12*th_12.conjugate())
        #f_m_p_cos_sq = (field**2)*abs(th_13*th_13.conjugate()+th_23*th_23.conjugate())
        #f_m_p_field = 0.5*field*(th_11+th_22) 

        f_nucl_sq = abs(f_nucl*f_nucl.conjugate())
        f_m_p_sin_sq = abs(0.5*(th_11*th_11.conjugate()+th_22*th_22.conjugate())+th_12*th_12.conjugate())
        f_m_p_cos_sq = abs(th_13*th_13.conjugate()+th_23*th_23.conjugate())
        f_m_p_field = 0.5*(th_11+th_22) 
        cross_sin = 2.*(f_nucl.real*f_m_p_field.real+f_nucl.imag*f_m_p_field.imag)
        
        return f_nucl_sq, f_m_p_sin_sq, f_m_p_cos_sq, cross_sin, refln

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
        zero_shift = float(self.ttheta_offset)
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


