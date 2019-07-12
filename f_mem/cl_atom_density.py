
"""
define classe to describe information about atom in cell
"""
__author__ = 'ikibalin'
__version__ = "2019_07_09"
import os
import numpy

import f_common.cl_variable 

class AtomDensity(dict):
    """
    Class to describe all information concerning the density in atom
    """
    def __init__(self, name=None):
        super(AtomDensity, self).__init__()
        self._p_name = None
        self._p_x = None
        self._p_y = None
        self._p_z = None
        self_p_structure_factor_tensor = None
        self._p_density = None
        
    def __repr__(self):
        ls_out = """CellDensity:\nname: {:}
            """.format(self._p_name)
        return ls_out


    def _refresh(self, name):
        if name is not None:
            self._p_name = name

    def set_val(self, name=None):
        self._refresh(name)

      
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
name is the name of mem
        """
        print(lsout)
        
    def calc_fourier_transform(self):
        entropy = None
        return entropy
        
        
    def calc_magnetic_structure_factor(self):
        flip_ratio = None
        return flip_ratio
        
    def calc_structure_factor_tensor(self):
        res = None
        return res
        