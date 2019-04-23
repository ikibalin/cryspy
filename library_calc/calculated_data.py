"""
define classes to describe calculated data
"""

__author__ = 'ikibalin'
__version__ = "2019_04_16"
import os
import numpy

from crystal import *

class CalculatedData1DPD(dict):
    """
    Calculate the model data for 1D powder diffraction experiment
    """
    def __init__(self, scale=1., crystal=Crystal()):
        super(CalculatedData1DPD, self).__init__()
        self._p_scale = None
        self._p_crystal = None
        self._refresh(scale, crystal)

    def __repr__(self):
        lsout = """Calculated data:\n scale {:}\n crystal: {:}""".format(
                self._p_scale, self._p_crystal)
        return lsout

    def _refresh(self, scale, crystal):
        if not(isinstance(scale, type(None))):
            self._p_scale = scale
        if not(isinstance(crystal, type(None))):
            self._p_crystal = crystal
            
    def set_val(self, scale=None, crystal=None):
        self._refresh(scale, crystal)
        
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
crystal is the definition of crystal 
        """
        print(lsout)
    
    def calc_iint(self, h, k,l):
        """
        calculate the integral intensity for h, k, l reflections
        """
        crystal = self._p_crystal
        f_nucl = crystal.calc_fn(h, k, l)
        iint = self._p_scale * abs(f_nucl*f_nucl.conjugate())
        print("   h   k   l       Iint")
        for h1, k1, l1, iint1 in zip(h, k, l, iint):
            print(" {:3} {:3} {:3} {:10.3f}".format(h1, k1, l1, iint1))
        return iint
    

        
if (__name__ == "__main__"):
  pass

