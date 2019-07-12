"""
define classes to describe calculated data
"""

__author__ = 'ikibalin'
__version__ = "2019_04_16"
import os
import numpy

from f_crystal.cl_crystal import *
from f_common.cl_variable import *


    
class CalculatedDataPowder2D(dict):
    """
    Calculate the model data for 2D powder diffraction experiment
    """
    def __init__(self, name=None, scale=1., field=1.):
        super(CalculatedDataPowder2D, self).__init__()
        self._p_name = None
        self._p_scale = None
        self._p_field = None
        self._refresh(name, scale, field)

    def __repr__(self):
        lsout = """CalculatedDataPowder2D:\n name: {:}\n scale: {:}
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
name is the name of CalculatedDataPowder2D
scale is the scale factor for crystal
field is the value of magnetic field applied along vertical direction in Tesla
        """
        print(lsout)
    
    def calc_for_iint(self, h, k, l, crystal, d_map={}):
        """
        calculate the integral intensity for h, k, l reflections
        Output: f_nucl_sq, f_m_p_sin_sq, f_m_p_cos_sq, cross_sin
        """
        if d_map == {}:
            d_map.update(self.plot_map())
        #if not(d_map["flag"]|(d_map["out"] is None)):
        #    f_nucl_sq, f_m_p_sin_sq, f_m_p_cos_sq, cross_sin = d_map["out"]
        #    return f_nucl_sq, f_m_p_sin_sq, f_m_p_cos_sq, cross_sin
        
        
        field = 1.*self._p_field

        #d_sf = d_map["sf"]
        #if not(d_sf["flag"]|(d_sf["out"] is None)):
        #    f_nucl, sft_11, sft_12, sft_13, sft_21, sft_22, sft_23, sft_31, sft_32, sft_33 = d_sf["out"]
        #else:
        f_nucl, sft_11, sft_12, sft_13, sft_21, sft_22, sft_23, sft_31, sft_32, sft_33, d_info_cry = crystal.calc_sf(h, k, l)
        
        cell = crystal.get_val("cell")
        
        #k_loc = cell.calc_k_loc(h, k, l)
        t_11, t_12, t_13, t_21, t_22, t_23, t_31, t_32, t_33 = cell.calc_m_t(h, k, l)
        
        th_11, th_12, th_13, th_21, th_22, th_23, th_31, th_32, th_33 = calc_mRmCmRT(
                t_11, t_21, t_31, t_12, t_22, t_32, t_13, t_23, t_33,
                sft_11, sft_12, sft_13, sft_21, sft_22, sft_23, sft_31, sft_32, 
                sft_33)
        f_nucl_sq = abs(f_nucl*f_nucl.conjugate())
        f_m_p_sin_sq = (field**2)*abs(0.5*(th_11*th_11.conjugate()+th_22*th_22.conjugate())+th_12*th_12.conjugate())
        f_m_p_cos_sq = (field**2)*abs(th_13*th_13.conjugate()+th_23*th_23.conjugate())
        f_m_p_field = 0.5*field*(th_11+th_22) 
        cross_sin = 2.*(f_nucl.real*f_m_p_field.real+f_nucl.imag*f_m_p_field.imag)
        #print("   h   k   l f_nucl_sq f_m_p_sin_sq f_m_p_cos_sq cross_sin")
        #for h_1, k_1, l_1, hh_1, hh_2, hh_3, hh_4  in zip(h, k, l, f_nucl_sq, f_m_p_sin_sq, f_m_p_cos_sq, cross_sin):
        #    print(""" {:3} {:3} {:3} {:9.3f} {:9.3f} {:9.3f} {:9.3f}""".format(
        #            h_1, k_1, l_1, hh_1, hh_2, hh_3, hh_4))        
        d_info_out = {"f_nucl_sq": f_nucl_sq, "f_m_p_sin_sq": f_m_p_sin_sq, 
                      "f_m_p_cos_sq": f_m_p_cos_sq, "cross_sin": cross_sin}   
        d_info_out.update(d_info_cry)        
        return f_nucl_sq, f_m_p_sin_sq, f_m_p_cos_sq, cross_sin, d_info_out 

    def plot_map(self):
        b_variable = self.is_variable()       
        d_map = {"flag": b_variable, "out":None}
        return d_map

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

