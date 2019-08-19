"""
define classes to calculated data
"""
__author__ = 'ikibalin'
__version__ = "2019_04_06"
import os
import numpy

from neupy.f_common.cl_variable import Variable
from neupy.f_experiment.f_powder_1d.cl_resolution_powder_1d import ResolutionPowder1D
from neupy.f_experiment.f_powder_1d.cl_factor_lorentz_powder_1d import FactorLorentzPowder1D
from neupy.f_experiment.f_powder_1d.cl_asymmetry_powder_1d import AsymmetryPowder1D
from neupy.f_experiment.f_powder_1d.cl_background_powder_1d import BackgroundPowder1D
    
from neupy.f_experiment.cl_beam_polarization import BeamPolarization

#Description of setup class


class SetupPowder1D(dict):
    """
    Class to describe characteristics of powder diffractometer
    """
    def __init__(self, label="exp", wave_length=1.4, zero_shift=0., 
                 resolution=ResolutionPowder1D(), factor_lorentz=FactorLorentzPowder1D(), 
                 asymmetry=AsymmetryPowder1D(), beam_polarization=BeamPolarization(),
                 background=BackgroundPowder1D()):
        super(SetupPowder1D, self).__init__()
        self._p_label = None
        self._p_wave_length = None
        self._p_zero_shift = None
        self._p_resolution = None
        self._p_factor_lorentz = None
        self._p_asymmetry = None
        self._p_beam_polarization = None
        self._p_background = None
        
        self._refresh(label, wave_length, zero_shift, resolution, 
                      factor_lorentz, asymmetry, beam_polarization, background)

    def __repr__(self):
        lsout = """SetupPowder1D:\n label: {:}\n wave_length: {:}
 zero_shift: {:}\n{:}\n{:}\n{:}\n{:}\n{:}""".format(self._p_label, 
 self._p_wave_length, self._p_zero_shift, self._p_resolution, 
 self._p_factor_lorentz, self._p_asymmetry, self._p_beam_polarization, 
 self._p_background)
        return lsout

    def _refresh(self, label, wave_length, zero_shift, resolution, 
                 factor_lorentz, asymmetry, beam_polarization, background):
        if not(isinstance(label, type(None))):
            self._p_label = label
        if not(isinstance(wave_length, type(None))):
            self._p_wave_length = wave_length
        if not(isinstance(zero_shift, type(None))):
            self._p_zero_shift = zero_shift
        if not(isinstance(resolution, type(None))):
            self._p_resolution = resolution
        if not(isinstance(factor_lorentz, type(None))):
            self._p_factor_lorentz = factor_lorentz
        if not(isinstance(asymmetry, type(None))):
            self._p_asymmetry = asymmetry
        if not(isinstance(beam_polarization, type(None))):
            self._p_beam_polarization = beam_polarization
        if not(isinstance(background, type(None))):
            self._p_background = background
            
    def set_val(self, label=None, wave_length=None, zero_shift=None, 
                resolution=None, factor_lorentz=None, asymmetry=None, 
                beam_polarization=None, background=None):
        
        self._refresh(label, wave_length, zero_shift, resolution, 
                      factor_lorentz, asymmetry, beam_polarization, background)
        
    def get_val(self, label):
        lab = "_p_"+label
        
        if lab in self.__dict__.keys():
            val = self.__dict__[lab]
            if isinstance(val, type(None)):
                self.set_val()
                val = self.__dict__[lab]
        else:
            print("The value '{:}' is not found".format(lab))
            val = None
        return val

    def list_vals(self):
        """
        give a list of parameters with small descripition
        """
        lsout = """
Parameters:
label is just to make labe
wave_length is to describe wave_length in angstrems
zero_shift is to describe zeroshift in degrees

resolution is a class to describe resolution of powder diffractometer
factor_lorentz is a class to describe factor Lorentz
asymmetry is a class to descibe the asymmetry
beam_polarization is a class to describe beam polarization
background  is Background class
        """
        print(lsout)

    def _gauss_pd(self, tth_2d):
        """
        one dimensional gauss powder diffraction
        """
        ag, bg = self._p_ag, self._p_bg
        val_1 = bg*tth_2d**2
        val_2 = numpy.where(val_1 < 5., numpy.exp(-val_1), 0.)
        self._p_gauss_pd = ag*val_2
        
    def _lor_pd(self, tth_2d):
        """
        one dimensional lorentz powder diffraction
        """
        al, bl = self._p_al, self._p_bl
        self._p_lor_pd = al*1./(1.+bl*tth_2d**2)
    
    def calc_shape_profile(self, tth, tth_hkl, i_g=0.):
        """
        calculate profile in the range ttheta for reflections placed on 
        ttheta_hkl with i_g parameter by defoult equal to zero
        
        tth, tth_hkl in degrees

        """
        
        resolution = self._p_resolution
        
        h_pv, eta, h_g, h_l, a_g, b_g, a_l, b_l = resolution.calc_resolution(
                                                  tth_hkl, i_g=i_g)

        
        tth_2d, tth_hkl_2d = numpy.meshgrid(tth, tth_hkl, indexing="ij")

        self._p_ag = numpy.meshgrid(tth, a_g, indexing="ij")[1]
        self._p_bg = numpy.meshgrid(tth, b_g, indexing="ij")[1]
        self._p_al = numpy.meshgrid(tth, a_l, indexing="ij")[1]
        self._p_bl = numpy.meshgrid(tth, b_l, indexing="ij")[1]
        eta_2d = numpy.meshgrid(tth, eta, indexing="ij")[1]
        self._p_eta = eta_2d 

        self._gauss_pd(tth_2d-tth_hkl_2d)
        self._lor_pd(tth_2d-tth_hkl_2d)
        g_pd_2d = self._p_gauss_pd 
        l_pd_2d = self._p_lor_pd
        
        profile_2d = eta_2d * l_pd_2d + (1.-eta_2d) * g_pd_2d
        
        return profile_2d
    
    def calc_profile(self, tth, tth_hkl, i_g):
        """
        tth and tth_hkl in degrees
        """
        zero_shift = 1*self._p_zero_shift
        tth_zs = tth-zero_shift
        np_shape_2d = self.calc_shape_profile(tth_zs, tth_hkl, i_g=i_g)
        asymmetry = self.get_val("asymmetry")
        np_ass_2d = asymmetry.calc_asymmetry(tth_zs, tth_hkl)
        factor_lorentz = self.get_val("factor_lorentz")
        np_lor_1d = factor_lorentz.calc_f_lorentz(tth_zs)
        
        
        np_lor_2d = numpy.meshgrid(np_lor_1d, tth_hkl, indexing="ij")[0]
        
        
        profile_2d = np_shape_2d*np_ass_2d*np_lor_2d
        return profile_2d 
    

    def calc_background(self, tth):
        """
        estimates background points on the given ttheta positions
        """
        background = self._p_background
        int_bkgd = background.interpolate_by_points(tth)
        return int_bkgd 
                    
    def is_variable(self):
        """
        without extinction
        """
        beam_polarization = self.get_val("beam_polarization")
        resolution = self.get_val("resolution")
        asymmetry = self.get_val("asymmetry")
        res = any([isinstance(self._p_zero_shift, Variable), 
                   beam_polarization.is_variable(), 
                   resolution.is_variable(), 
                   asymmetry.is_variable()])
        return res   

    def get_variables(self):
        l_variable = []
        if isinstance(self._p_zero_shift, Variable):
            l_variable.append(self._p_zero_shift)
        
        background = self.get_val("background")
        l_var = background.get_variables()
        l_variable.extend(l_var)
        
        beam_polarization = self.get_val("beam_polarization")
        l_var = beam_polarization.get_variables()
        l_variable.extend(l_var)
        
        resolution = self.get_val("resolution")
        l_var = resolution.get_variables()
        l_variable.extend(l_var)
        
        asymmetry = self.get_val("asymmetry")
        l_var = asymmetry.get_variables()
        l_variable.extend(l_var)
        
        return l_variable

