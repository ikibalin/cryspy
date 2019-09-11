"""
define classes to calculated data
"""
__author__ = 'ikibalin'
__version__ = "2019_04_06"
import os
import numpy
import scipy.interpolate

#class BeamPolarization
from neupy.f_experiment.cl_beam_polarization import BeamPolarization
from neupy.f_common.cl_variable import Variable
#Description of setup class

class FactorLorentzPowder2D(dict):
    """
    Lorentz Factor for one dimensional powder diffraction
    """
    def __init__(self):
        super(FactorLorentzPowder2D, self).__init__()
        dd= {}
        self.update(dd)
        
    def __repr__(self):
        lsout = """FactorLorentzPowder2D:"""
        return lsout

    def _refresh(self):
        print("'_refresh' is not introduced for FactorLorentzPD")
            
    def set_val(self):
        print("'set_val' is not introduced for FactorLorentzPD")
        
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
no parameters
        """
        print(lsout)
        
    
    def calc_f_lorentz(self, tth):
        """
        Lorentz factor
        tth should be in degrees
        """
        tth_rad = tth*numpy.pi/180.
        factor_lorentz = 1./(numpy.sin(tth_rad)*numpy.sin(0.5*tth_rad))
        return factor_lorentz 


