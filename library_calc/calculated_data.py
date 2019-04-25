"""
define classes to describe calculated data
"""

__author__ = 'ikibalin'
__version__ = "2019_04_16"
import os
import numpy

from crystal import *

class CalculatedDataPowder1D(dict):
    """
    Calculate the model data for 1D powder diffraction experiment
    """
    def __init__(self, scale=1., field=1., crystal=Crystal()):
        super(CalculatedDataPowder1D, self).__init__()
        self._p_scale = None
        self._p_field = None
        self._p_crystal = None
        self._refresh(scale, field, crystal)

    def __repr__(self):
        lsout = """Calculated data 1D:\n scale {:}\n field {:}\n{:}""".format(
                self._p_scale, self._p_field, self._p_crystal)
        return lsout

    def _refresh(self, scale, field, crystal):
        if not(isinstance(scale, type(None))):
            self._p_scale = scale
        if not(isinstance(field, type(None))):
            self._p_field = field
        if not(isinstance(crystal, type(None))):
            self._p_crystal = crystal
            
    def set_val(self, scale=None, field=None, crystal=None):
        self._refresh(scale, field, crystal)
        
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
scale is the scale factor for crystal
field is the value of magnetic field applied along vertical direction in Tesla
crystal is the definition of crystal 
        """
        print(lsout)
    
    def calc_iint(self, h, k, l, beam_polarization):
        """
        calculate the integral intensity for h, k, l reflections
        """
        crystal = self._p_crystal
        field = self._p_field
        p_u = beam_polarization.get_val("p_u")
        p_d = beam_polarization.get_val("p_d")

        f_nucl, sft_11, sft_12, sft_13, sft_21, sft_22, sft_23, sft_31, sft_32, sft_33 = crystal.calc_sf(h, k, l)
        
        cell = crystal.get_val("cell")
        
        #k_loc = cell.calc_k_loc(h, k, l)
        t_11, t_12, t_13, t_21, t_22, t_23, t_31, t_32, t_33 = cell.calc_m_t(h, k, l)
        
        th_11, th_12, th_13, th_21, th_22, th_23, th_31, th_32, th_33 = calc_mRmCmRT(
                t_11, t_21, t_31, t_12, t_22, t_32, t_13, t_23, t_33,
                sft_11, sft_12, sft_13, sft_21, sft_22, sft_23, sft_31, sft_32, 
                sft_33)
        fm_p_sq = (field**2)*abs(0.5*(th_11*th_11.conjugate()+th_22*th_22.conjugate())+th_12*th_12.conjugate())
        fm_p_field = field*0.5*(th_11+th_22) 
        cross = 2.*(f_nucl.real*fm_p_field.real+f_nucl.imag*fm_p_field.imag)
        #lkloc=[cfunc.calck(hkl,mB) for hkl in lhkl]

        iint_u = self._p_scale * (abs(f_nucl*f_nucl.conjugate()) +
                 fm_p_sq + p_u*cross)

        iint_d = self._p_scale * (abs(f_nucl*f_nucl.conjugate()) +
                 fm_p_sq - p_d*cross)
        
        #I_p, I_m = self.calc_extinc_powder(h, k, l, fn, fm_perp_eup, fm_p_sq,ext, p_up, p_down, ucp, wavelength)
        
        print("   h   k   l fn_real fn_imag")
        for h1, k1, l1, f, c11, c12, c13, c21, c22, c23, c31, c32, c33 in zip(h, k, l, f_nucl, sft_11, sft_12, sft_13, sft_21, sft_22, sft_23, sft_31, sft_32, sft_33):
            print(""" {:3} {:3} {:3} {:7.3f} {:7.3f}
            {:7.3f}+{:7.3f}i {:7.3f}+{:7.3f}i {:7.3f}+{:7.3f}i
            {:7.3f}+{:7.3f}i {:7.3f}+{:7.3f}i {:7.3f}+{:7.3f}i
            {:7.3f}+{:7.3f}i {:7.3f}+{:7.3f}i {:7.3f}+{:7.3f}i""".format(
                    h1, k1, l1, f.real, f.imag, c11.real, c11.imag, c12.real, 
                    c12.imag, c13.real, c13.imag, c21.real, c21.imag, c22.real, 
                    c22.imag, c23.real, c23.imag, c31.real, c31.imag, c32.real, 
                    c32.imag, c33.real, c33.imag))
        return iint_u, iint_d
    

class CalculatedDataPowder2D(dict):
    """
    Calculate the model data for 2D powder diffraction experiment
    """
    def __init__(self, scale=1., field=1., crystal=Crystal()):
        super(CalculatedDataPowder2D, self).__init__()
        self._p_scale = None
        self._p_field = None
        self._p_crystal = None
        self._refresh(scale, field, crystal)

    def __repr__(self):
        lsout = """Calculated data 2D:\n scale {:}\n field {:}\n{:}""".format(
                self._p_scale, self._p_field, self._p_crystal)
        return lsout

    def _refresh(self, scale, field, crystal):
        if not(isinstance(scale, type(None))):
            self._p_scale = scale
        if not(isinstance(field, type(None))):
            self._p_field = field
        if not(isinstance(crystal, type(None))):
            self._p_crystal = crystal
            
    def set_val(self, scale=None, field=None, crystal=None):
        self._refresh(scale, field, crystal)
        
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
scale is the scale factor for crystal
field is the value of magnetic field applied along vertical direction in Tesla
crystal is the definition of crystal 
        """
        print(lsout)
    
    def calc_for_iint(self, h, k, l):
        """
        calculate the integral intensity for h, k, l reflections
        """
        crystal = self._p_crystal
        field = self._p_field
        p_u = beam_polarization.get_val("p_u")
        p_d = beam_polarization.get_val("p_d")

        f_nucl, sft_11, sft_12, sft_13, sft_21, sft_22, sft_23, sft_31, sft_32, sft_33 = crystal.calc_sf(h, k, l)
        
        cell = crystal.get_val("cell")
        
        #k_loc = cell.calc_k_loc(h, k, l)
        t_11, t_12, t_13, t_21, t_22, t_23, t_31, t_32, t_33 = cell.calc_m_t(h, k, l)
        
        th_11, th_12, th_13, th_21, th_22, th_23, th_31, th_32, th_33 = calc_mRmCmRT(
                t_11, t_21, t_31, t_12, t_22, t_32, t_13, t_23, t_33,
                sft_11, sft_12, sft_13, sft_21, sft_22, sft_23, sft_31, sft_32, 
                sft_33)
        fm_p_sq = (field**2)*abs(0.5*(th_11*th_11.conjugate()+th_22*th_22.conjugate())+th_12*th_12.conjugate())
        fm_p_field = field*0.5*(th_11+th_22) 
        cross = 2.*(f_nucl.real*fm_p_field.real+f_nucl.imag*fm_p_field.imag)
        return fm_p_sq, fm_p_field, cross
        
if (__name__ == "__main__"):
  pass

