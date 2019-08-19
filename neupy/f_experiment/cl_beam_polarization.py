"""
define classes to calculated data
"""
__author__ = 'ikibalin'
__version__ = "2019_04_06"
import os
import numpy

from neupy.f_common.cl_variable import Variable
    
#Description of setup class

class BeamPolarization(dict):
    """
    Describe the polarisation of the beam
    """
    def __init__(self, p_u = 1.0, flipper_efficiency = 1.0):
        super(BeamPolarization, self).__init__()
        self._p_p_u = None
        self._p_flipper_efficiency = None
        
        self._refresh(p_u, flipper_efficiency)
        
    def __repr__(self):
        lsout = """BeamPolarization: \n p_u: {:}\n flipper_efficiency: {:}""".format(
                self.get_val("p_u"), self.get_val("flipper_efficiency"))
        return lsout

    def _refresh(self, p_u, flipper_efficiency):
        if p_u is not None:
            self._p_p_u = p_u
        if flipper_efficiency is not None:
            self._p_flipper_efficiency = flipper_efficiency
            
    def set_val(self, p_u=None, flipper_efficiency=None):
        self._refresh(p_u, flipper_efficiency)
        
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
p_u, flipper_efficiency is describe the polarization of the incident beam and the flipper efficiency p_down = (2*eff-1)*p_up: (p_up = (n_up - n_down) / (n_up + n_down))
         and down
        """
        print(lsout)

    def is_variable(self):
        """
        without extinction
        """
        res = any([isinstance(self._p_p_u, Variable), 
                   isinstance(self._p_flipper_efficiency, Variable)])
        return res        
    
    def get_variables(self):
        l_variable = []
        if isinstance(self._p_p_u, Variable):
            l_variable.append(self._p_p_u)
        if isinstance(self._p_flipper_efficiency, Variable):
            l_variable.append(self._p_flipper_efficiency)
        return l_variable
        
