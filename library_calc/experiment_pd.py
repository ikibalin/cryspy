
"""
define classe to describe experiment
"""
__author__ = 'ikibalin'
__version__ = "2019_04_06"
import os
import numpy

from observed_data import *
from calculated_data import *
from setup_1d_pd import *

class Experiment1DPD(dict):
    """
    Class to describe all information concerning to the experiment
    """
    def __init__(self, setup = Setup1DPD(), 
                 list_calculated_data = [], 
                 observed_data = ObservedData1DPD()):
        super(Experiment1DPD, self).__init__()
        self._p_setup = None
        self._list_calculated_data = []
        self._p_observed_data = None
        
        self._refresh(setup, observed_data)

    def __repr__(self):
        lsout = """Experiment:\n setup {:}\n observed data: {:} """.format(
                self._p_setup, self._p_observed_data)
        return lsout

    def _refresh(self, setup, observed_data):
        if not(isinstance(setup, type(None))):
            self._p_setup = setup
        if not(isinstance(observed_data, type(None))):
            self._p_observed_data = observed_data

            
    def set_val(self, setup=None, observed_data=None):
        self._refresh(setup, observed_data)
        
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
setup is to describe parameters of diffractometer 
observed_data is the experimental data
        """
        print(lsout)
    def add_calculated_data(self, observed_data):
        self._list_calculated_data.append(observed_data)

    def del_calculated_data(self, ind):
        self._list_calculated_data.pop(ind)        

    def replace_calculated_data(self, ind, observed_data):
        self._list_calculated_data.pop(ind)
        self._list_calculated_data.insert(ind, observed_data)
    
    def calc_profile(self, tth):
        """
        calculate intensity for the given diffraction angle
        """
        setup = self._p_setup
        background = setup.calc_background(tth)
        wavelength = setup.get_val("wavelength")
        tth_min = tth.min()
        tth_max = tth.max()
        sthovl_min = numpy.sin(0.5*tth_min*numpy.pi/180.)/wavelength
        sthovl_max = numpy.sin(0.5*tth_max*numpy.pi/180.)/wavelength

        res_1d = numpy.zeros(tth.shape[0], dtype=float)

        for calculated_data in self._list_calculated_data:
            
            i_g = calculated_data.get_val("crystal").get_val("i_g")
            cell = calculated_data.get_val("crystal").get_val("cell")
            space_groupe = calculated_data.get_val("crystal").get_val("space_groupe")
            h, k, l, mult = setup.calc_hkl(cell, space_groupe, sthovl_min, sthovl_max)
            
            np_iint = calculated_data.calc_iint(h, k, l)
            sthovl_hkl = cell.calc_sthovl(h, k, l)
            
            tth_hkl_rad = 2.*numpy.arcsin(sthovl_hkl*wavelength)
            tth_hkl = tth_hkl_rad*180./numpy.pi
            profile_2d = setup.calc_profile(tth, tth_hkl, i_g)
            
            
            np_iint_2d = numpy.meshgrid(tth, np_iint*mult, indexing="ij")[1]
            
            res_2d = profile_2d*np_iint_2d 

            res_1d += res_2d.sum(axis=1) 
            
        return res_1d+background
    
    def calc_chi_sq(self):
        """
        calculate chi square
        """
        observed_data = self._p_observed_data
        tth = observed_data.get_val('tth')
        int_u_exp = observed_data.get_val('int_u')
        sint_u_exp = observed_data.get_val('sint_u')
        int_u_mod = self.calc_profile(tth)
        chi_sq = ((int_u_mod-int_u_exp)/sint_u_exp)**2
        return chi_sq.sum()
    
if (__name__ == "__main__"):
  pass
