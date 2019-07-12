
"""
define classe to describe information about density in cell
"""
__author__ = 'ikibalin'
__version__ = "2019_07_09"
import os
import numpy

import f_mem.cl_atom_density

import f_common.cl_variable 

class CellDensity(dict):
    """
    Class to describe all information concerning the density in cell
    """
    def __init__(self, name=None, points_number_a=None, 
                 points_number_b=None, points_number_c=None,
                 file_dir=None, file_name=None):
        super(CellDensity, self).__init__()
        self._p_name = None
        self._p_points_number_a = None
        self._p_points_number_b = None
        self._p_points_number_c = None
        self._p_file_dir = None
        self._p_file_name = None
        self._list_atom_density = []
        self._refresh(name, points_number_a, 
                 points_number_b, points_number_c,
                 file_dir, file_name)
        
    def __repr__(self):
        ls_out = """CellDensity:\nname: {:}\n points_number_a: {:}
 points_number_b: {:}\n points_number_c: {:}\n file_dir: {:}
 file_name: {:}""".format(self._p_name, self._p_points_number_a, 
 self._p_points_number_b, self._p_points_number_c, self._p_file_dir, self._p_file_name)
        for atom_density in self._list_atom_density:
            ls_out += "{:}".format(atom_density)
        return ls_out


    def _refresh(self, name, points_number_a, points_number_b, 
                 points_number_c, file_dir, file_name):
        if name is not None:
            self._p_name = name
        if points_number_a is not None:
            self._p_points_number_a = int(points_number_a)
        if points_number_b is not None:
            self._p_points_number_b = int(points_number_b)
        if points_number_c is not None:
            self._p_points_number_c = int(points_number_c)
        if file_dir is not None:
            self._p_file_dir = file_dir
        if file_name is not None:
            self._p_file_name = file_name

    def set_val(self, name=None, points_number_a=None, 
                 points_number_b=None, points_number_c=None,
                 file_dir=None, file_name=None):
        self._refresh(name, points_number_a, points_number_b, 
                 points_number_c, file_dir, file_name)

      
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

    def create_density(self):
        points_number_a = 1*self._p_points_number_a
        points_number_b = 1*self._p_points_number_b
        points_number_c = 1*self._p_points_number_c
        #np_frac_x = numpy.linspace(0., 1., points_number_a, endpoint=False)
        #np_frac_y = numpy.linspace(0., 1., points_number_b, endpoint=False)
        #np_frac_z = numpy.linspace(0., 1., points_number_c, endpoint=False)

        #np_frac_x_3d, np_frac_y_3d, np_frac_z_3d = numpy.meshgrid([np_frac_x, np_frac_y, np_frac_z], indexing ="ij")

        val = 1./float(points_number_a*points_number_b*points_number_c)
        propability = val*numpy.ones(shape=(points_number_a, points_number_b, points_number_c), dtype=float, order='C')

        
    

    def calc_fourier_transform(self):
        entropy = None
        return entropy
        
    def write_density(self, f_name):
        chi_sq = None
        return chi_sq
        
    def read_density(self, f_name):
        minimizer = None
        return minimizer
        
    def calc_magnetic_structure_factor(self):
        flip_ratio = None
        return flip_ratio
        
    def calc_structure_factor_tensor(self):
        res = None
        return res
        