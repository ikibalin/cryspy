
"""
define classe to describe information about MEM
"""
__author__ = 'ikibalin'
__version__ = "2019_07_09"
import os
import numpy

import f_crystal.cl_crystal


import f_mem.cl_cell_density
import f_mem.cl_observed_data_mem

import f_common.cl_variable 

class MemReconstruction(dict):
    """
    Class to describe all information concerning the mem reconstruction
    """
    def __init__(self, name=None, file_dir=None, file_out=None,
                 cell_density=None):
        super(MemReconstruction, self).__init__()
        self._p_name = None
        self._p_file_dir = None
        self._p_file_out = None

        self._p_cell_density = None
        self._list_crystal = []
        self._list_observed_data_mem = []

        self._refresh(name, file_dir, file_out, cell_density)

    def __repr__(self):
        ls_out = """MemReconstruction:\nname: {:} \n file_dir: {:}
 file_out: {:}\n{:}""".format(self._p_name, self._p_file_dir, self._p_file_out, 
 self._p_cell_density)
        for crystal in self._list_crystal:
            ls_out += "{:}".format(crystal)
        for observed_data_mem in self._list_observed_data_mem:
            ls_out += "{:}".format(observed_data_mem)
        return ls_out


    def _refresh(self, name, file_dir, file_out, cell_density):
        if name is not None:
            self._p_name = name
        if file_dir is not None:
            self._p_file_dir = file_dir
        if file_out is not None:
            self._p_file_out = file_out
        if cell_density is not None:
            self._p_cell_density = cell_density



    def set_val(self, name=None, file_dir=None, file_out=None, cell_density=None):
        self._refresh(name, file_dir, file_out, cell_density)

      
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
        
    def calc_entropy(self):
        entropy = None
        return entropy
        
    def calc_chi_sq(self):
        chi_sq = None
        return chi_sq
        
    def calc_minimizer(self):
        minimizer = None
        return minimizer
        
    def calc_flip_ratio(self):
        flip_ratio = None
        return flip_ratio
        
    def minimize(self):
        res = None
        return res
        
    def add_crystal(self, crystal):
        self._list_crystal.append(crystal)

    def add_observed_data_mem(self, observed_data_mem):
        self._list_observed_data_mem.append(observed_data_mem)
