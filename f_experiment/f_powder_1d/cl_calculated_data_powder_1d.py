"""
define classes to describe calculated data
"""

__author__ = 'ikibalin'
__version__ = "2019_04_16"
import os
import numpy

from f_crystal.cl_crystal import *
from f_common.cl_variable import *


class CalculatedDataPowder1D(dict):
    """
    Calculate the model data for 1D powder diffraction experiment
    """
    def __init__(self, name=None, scale=1., field=1.):
        super(CalculatedDataPowder1D, self).__init__()
        self._p_name = None
        self._p_scale = None
        self._p_field = None
        self._refresh(name, scale, field)

    def __repr__(self):
        lsout = """CalculatedDataPowder1D: \n name: {:}\n scale: {:}
 field: {:}""".format(self._p_name, self._p_scale, self._p_field)
        return lsout

    def _refresh(self, name, scale, field):
        if name is not None:
            self._p_name = name
        if scale is not None:
            self._p_scale = scale
        if field is not None:
            self._p_field = field
            
    def set_val(self, name=None, scale=None, field=None):
        self._refresh(name, scale, field)
        
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
name is the name of CalculatedData1D
scale is the scale factor for crystal
field is the value of magnetic field applied along vertical direction in Tesla
crystal is the definition of crystal 
        """
        print(lsout)
    
    def calc_iint(self, h, k, l, beam_polarization, crystal):
        """
        calculate the integral intensity for h, k, l reflections
        """
        field = 1.*self._p_field
        p_u = 1.*beam_polarization.get_val("p_u")
        p_d = (2.*beam_polarization.get_val("flipper_efficiency")-1)*p_u
        

        f_nucl, sft_11, sft_12, sft_13, sft_21, sft_22, sft_23, sft_31, sft_32, sft_33, d_info_cry = crystal.calc_sf(h, k, l)
        
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

        iint_u = abs(f_nucl*f_nucl.conjugate()) + fm_p_sq + p_u*cross
                 

        iint_d = abs(f_nucl*f_nucl.conjugate()) + fm_p_sq - p_d*cross
        d_info_out = {"iint_u": iint_u, "iint_d": iint_d}   
        d_info_out.update(d_info_cry)
        
        #I_p, I_m = self.calc_extinc_powder(h, k, l, fn, fm_perp_eup, fm_p_sq,ext, p_up, p_down, ucp, wave_length)
        #
        #print("   h   k   l fn_real fn_imag")
        #for h1, k1, l1, f, c11, c12, c13, c21, c22, c23, c31, c32, c33 in zip(h, k, l, f_nucl, sft_11, sft_12, sft_13, sft_21, sft_22, sft_23, sft_31, sft_32, sft_33):
        #    print(""" {:3} {:3} {:3} {:7.3f} {:7.3f}
        #    {:7.3f}+{:7.3f}i {:7.3f}+{:7.3f}i {:7.3f}+{:7.3f}i
        #    {:7.3f}+{:7.3f}i {:7.3f}+{:7.3f}i {:7.3f}+{:7.3f}i
        #    {:7.3f}+{:7.3f}i {:7.3f}+{:7.3f}i {:7.3f}+{:7.3f}i""".format(
        #            h1, k1, l1, f.real, f.imag, c11.real, c11.imag, c12.real, 
        #            c12.imag, c13.real, c13.imag, c21.real, c21.imag, c22.real, 
        #            c22.imag, c23.real, c23.imag, c31.real, c31.imag, c32.real, 
        #            c32.imag, c33.real, c33.imag))
        return iint_u, iint_d, d_info_out
    
    def is_variable(self):
        """
        without extinction
        """
        res = any([isinstance(self._p_scale, Variable)])
        return res        

    def get_variables(self):
        l_variable = []
        if isinstance(self._p_scale, Variable):
            l_variable.append(self._p_scale)
        return l_variable

    
if (__name__ == "__main__"):
  pass

