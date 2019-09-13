"""
define classes to calculated data
"""
__author__ = 'ikibalin'
__version__ = "2019_09_06"
import os
import numpy


from pycifstar import Data
from cryspy.f_common.cl_fitable import Fitable


#Description of setup class

class BeamPolarization(object):
    """
    Describe the beam polarisation

    Example:
    
    _diffrn_radiation_efficiency 1.0
    _diffrn_radiation_polarization -0.87

    """
    def __init__(self, polarization = Fitable(1.0), efficiency = Fitable(1.0)):
        super(BeamPolarization, self).__init__()
        self.__polarization = None
        self.__efficiency = None

        self.polarization = polarization
        self.efficiency = efficiency
        
    @property
    def polarization(self):
        """
        The polarization of the incident beam. 

        The permitted range is -1.0 -> 1.0
        """
        return self.__polarization
    @polarization.setter
    def polarization(self, x):
        if isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        self.__polarization = x_in

    @property
    def efficiency(self):
        """
        The efficiency of the efficiency. 

        The permitted range is -1.0 -> 1.0
        """
        return self.__efficiency
    @efficiency.setter
    def efficiency(self, x):
        if isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        self.__efficiency = x_in


    def __repr__(self):
        ls_out = []
        ls_out.append("BeamPolarization: ")
        ls_out.append(" polarization: {}".format(self.polarization.print_with_sigma))
        ls_out.append(" efficiency: {}".format(self.efficiency.print_with_sigma))
        return "\n".join(ls_out)


    @property
    def is_variable(self):
        """
        Output: True if there is any refined parameter
        """
        res = any([self.polarization.refinement, 
                   self.efficiency.refinement])
        return res        
    
    def get_variables(self):
        """
        Output: the list of the refined parameters
        """
        l_variable = []
        if self.polarization.refinement:
            l_variable.append(self.polarization)
        if self.efficiency.refinement:
            l_variable.append(self.efficiency)
        return l_variable

    def _show_message(self, s_out: str):
        print("***  Error ***")
        print(s_out)
        
    @property
    def to_cif(self):
        ls_out = ["_diffrn_radiation_polarization {:}".format(self.polarization.print_with_sigma)]
        ls_out.append("_diffrn_radiation_efficiency {:}".format(self.efficiency.print_with_sigma))
        return "\n".join(ls_out)


    def from_cif(self, string: str):
        cif_data = Data()
        flag = cif_data.take_from_string(string)
        if not flag:
            return False
        flag = False
        if cif_data.is_value("_diffrn_radiation_polarization"):
            self.polarization = cif_data["_diffrn_radiation_polarization"] # CIFvalue
        if cif_data.is_value("_diffrn_radiation_efficiency"):
            self.efficiency = cif_data["_diffrn_radiation_efficiency"] # CIFvalue
        return True

