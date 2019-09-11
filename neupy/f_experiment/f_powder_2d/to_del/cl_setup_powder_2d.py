"""
define classes to calculated data
"""
__author__ = 'ikibalin'
__version__ = "2019_04_06"
import os
import numpy
import scipy.interpolate

#class BeamPolarization
from neupy.f_experiment.f_powder_2d.cl_asymmetry_powder_2d import AsymmetryPowder2D
from neupy.f_experiment.f_powder_2d.cl_background_powder_2d import BackgroundPowder2D
from neupy.f_experiment.f_powder_2d.cl_factor_lorentz_powder_2d import FactorLorentzPowder2D
from neupy.f_experiment.f_powder_2d.cl_resolution_powder_2d import ResolutionPowder2D


from neupy.f_experiment.cl_beam_polarization import BeamPolarization
from neupy.f_common.cl_variable import Variable
#Description of setup class


class SetupPowder2D(dict):
    """
    Class to describe characteristics of powder diffractometer
    """
    def __init__(self, label="exp", wave_length=1.4, zero_shift=0., 
                 resolution=ResolutionPowder2D(), factor_lorentz=FactorLorentzPowder2D(), 
                 asymmetry=AsymmetryPowder2D(), beam_polarization=BeamPolarization(),
                 background=BackgroundPowder2D()):
        super(SetupPowder2D, self).__init__()
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
        lsout = """SetupPowder2D:\n label: {:}\n wave_length: {:}
 zero_shift: {:}\n{:}\n{:}\n{:}\n{:}\n{:}""".format(self._p_label, 
 self._p_wave_length, self._p_zero_shift, self._p_resolution, 
 self._p_factor_lorentz, self._p_asymmetry, self._p_beam_polarization, 
 self._p_background)
        return lsout

    def _refresh(self, label, wave_length, zero_shift, resolution, 
                 factor_lorentz, asymmetry, beam_polarization, background):
        if label is not None:
            self._p_label = label
        if wave_length is not None:
            self._p_wave_length = wave_length
        if zero_shift is not None:
            self._p_zero_shift = zero_shift
        if resolution is not None:
            self._p_resolution = resolution
        if factor_lorentz is not None:
            self._p_factor_lorentz = factor_lorentz
        if asymmetry is not None:
            self._p_asymmetry = asymmetry
        if beam_polarization is not None:
            self._p_beam_polarization = beam_polarization
        if background is not None:
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

    def _gauss_pd(self, tth_3d):
        """
        one dimensional gauss powder diffraction
        """
        ag, bg = self._p_ag, self._p_bg
        val_1 = bg*tth_3d**2
        val_2 = numpy.where(val_1 < 5., numpy.exp(-val_1), 0.)
        self._p_gauss_pd = ag*val_2
        
    def _lor_pd(self, tth_3d):
        """
        one dimensional lorentz powder diffraction
        """
        al, bl = self._p_al, self._p_bl
        self._p_lor_pd = al*1./(1.+bl*tth_3d**2)
    
    def calc_shape_profile(self, tth, phi, tth_hkl, i_g=0.):
        """
        calculate profile in the range ttheta for reflections placed on 
        ttheta_hkl with i_g parameter by defoult equal to zero
        
        tth, phi, tth_hkl in degrees

        """
        
        resolution = self._p_resolution
        
        h_pv, eta, h_g, h_l, a_g, b_g, a_l, b_l = resolution.calc_resolution(
                                                  tth_hkl, i_g=i_g)

        
        tth_3d, phi_3d, tth_hkl_3d = numpy.meshgrid(tth, phi, tth_hkl, indexing="ij")

        self._p_ag = numpy.meshgrid(tth, phi, a_g, indexing="ij")[2]
        self._p_bg = numpy.meshgrid(tth, phi, b_g, indexing="ij")[2]
        self._p_al = numpy.meshgrid(tth, phi, a_l, indexing="ij")[2]
        self._p_bl = numpy.meshgrid(tth, phi, b_l, indexing="ij")[2]
        eta_3d = numpy.meshgrid(tth, phi, eta, indexing="ij")[2]
        self._p_eta = eta_3d 

        self._gauss_pd(tth_3d-tth_hkl_3d)
        self._lor_pd(tth_3d-tth_hkl_3d)
        g_pd_3d = self._p_gauss_pd 
        l_pd_3d = self._p_lor_pd
        
        profile_3d = eta_3d * l_pd_3d + (1.-eta_3d) * g_pd_3d
        
        return profile_3d
    
    def calc_profile(self, tth, phi, tth_hkl, i_g, d_setup_in={}):
        """
        tth and tth_hkl in degrees
        """
        d_setup_out = {}

        zero_shift = 1*self._p_zero_shift
        
        tth_zs = tth-zero_shift

        np_shape_3d = self.calc_shape_profile(tth_zs, phi, tth_hkl, i_g=i_g)
        d_setup_out["np_shape_3d"] = np_shape_3d

        asymmetry = self.get_val("asymmetry")
        #dimension (tth, hkl)
        np_ass_2d = asymmetry.calc_asymmetry(tth_zs, tth_hkl)
        d_setup_out["np_ass_2d"] = np_ass_2d

        #dimension (tth, phi, hkl)
        np_ass_3d = np_ass_2d[:, numpy.newaxis,:]*numpy.ones(phi.size, dtype=float)[numpy.newaxis, :, numpy.newaxis]
        
        factor_lorentz = self.get_val("factor_lorentz")
        np_lor_1d = factor_lorentz.calc_f_lorentz(tth_zs)
        
        
        np_lor_3d = numpy.meshgrid(np_lor_1d, phi, tth_hkl, indexing="ij")[0]
        
        
        profile_3d = np_shape_3d*np_ass_3d*np_lor_3d
        return profile_3d, d_setup_out
    

    def calc_background(self, tth, phi):
        """
        estimates background points on the given ttheta positions
        """                   
        background = self._p_background
        int_bkgd = background.interpolate_by_points(tth, phi)
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
