"""
define classes to calculated data
"""
__author__ = 'ikibalin'
__version__ = "2019_04_06"
import os
import numpy

from f_experiment.f_powder_1d.cl_setup_powder_1d import BeamPolarization
from f_common.cl_variable import *
    

#Description of setup class

class SetupSingle(dict):
    """
    Class to describe characteristics of powder diffractometer
    """
    def __init__(self, label="exp", wave_length=1.4, beam_polarization=BeamPolarization()):
        super(SetupSingle, self).__init__()
        self._p_label = None
        self._p_wave_length = None
        self._p_beam_polarization = None
        
        self._refresh(label, wave_length, beam_polarization)

    def __repr__(self):
        lsout = """SetupSingle:\n label: {:}\n wave_length: {:}
{:}""".format(self._p_label, self._p_wave_length, 
 self._p_beam_polarization)
        return lsout

    def _refresh(self, label, wave_length, beam_polarization):
        if label is not None:
            self._p_label = label
        if wave_length  is not None:
            self._p_wave_length = wave_length
        if beam_polarization  is not None:
            self._p_beam_polarization = beam_polarization
            
    def set_val(self, label=None, wave_length=None, beam_polarization=None):
        
        self._refresh(label, wave_length, beam_polarization)
        
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
Parameters for setup_single:
label is just to make labe
wave_length is to describe wave_length in angstrems

beam_polarization is a class to describe beam polarization
        """
        print(lsout)

    def is_variable(self):
        """
        without extinction
        """
        beam_polarization = self.get_val("beam_polarization")
        res = any([beam_polarization.is_variable()])
        return res   
    
    def get_variables(self):
        l_variable = []
        beam_polarization = self.get_val("beam_polarization")
        l_var = beam_polarization.get_variables()
        l_variable.extend(l_var)
        return l_variable    