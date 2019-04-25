
"""
define classe to describe experiment
"""
__author__ = 'ikibalin'
__version__ = "2019_04_06"
import os
import numpy

from observed_data import *
from calculated_data import *
from setup_powder_1d import *

class ExperimentPowder1D(dict):
    """
    Class to describe all information concerning to the experiment
    """
    def __init__(self, setup = SetupPowder1D(), 
                 list_calculated_data = [], 
                 observed_data = ObservedDataPowder1D()):
        super(ExperimentPowder1D, self).__init__()
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
        beam_polarization = setup.get_val("beam_polarization")
        
        tth_min = tth.min()
        tth_max = tth.max()
        sthovl_min = numpy.sin(0.5*tth_min*numpy.pi/180.)/wavelength
        sthovl_max = numpy.sin(0.5*tth_max*numpy.pi/180.)/wavelength

        res_u_1d = numpy.zeros(tth.shape[0], dtype=float)
        res_d_1d = numpy.zeros(tth.shape[0], dtype=float)

        for calculated_data in self._list_calculated_data:
            
            i_g = calculated_data.get_val("crystal").get_val("i_g")
            cell = calculated_data.get_val("crystal").get_val("cell")
            space_groupe = calculated_data.get_val("crystal").get_val("space_groupe")
            h, k, l, mult = setup.calc_hkl(cell, space_groupe, sthovl_min, sthovl_max)
            
            np_iint_u, np_iint_d = calculated_data.calc_iint(h, k, l, beam_polarization)
            sthovl_hkl = cell.calc_sthovl(h, k, l)
            
            tth_hkl_rad = 2.*numpy.arcsin(sthovl_hkl*wavelength)
            tth_hkl = tth_hkl_rad*180./numpy.pi
            profile_2d = setup.calc_profile(tth, tth_hkl, i_g)
            
            
            np_iint_u_2d = numpy.meshgrid(tth, np_iint_u*mult, indexing="ij")[1]

            np_iint_d_2d = numpy.meshgrid(tth, np_iint_d*mult, indexing="ij")[1]
            
            res_u_2d = profile_2d*np_iint_u_2d 
            res_d_2d = profile_2d*np_iint_d_2d 

            res_u_1d += res_u_2d.sum(axis=1) 
            res_d_1d += res_d_2d.sum(axis=1) 
            
        return res_u_1d+background, res_d_1d+background
    
    def calc_chi_sq(self):
        """
        calculate chi square
        """
        observed_data = self._p_observed_data

        tth = observed_data.get_val('tth')
        int_u_exp = observed_data.get_val('int_u')
        sint_u_exp = observed_data.get_val('sint_u')
        int_d_exp = observed_data.get_val('int_d')
        sint_d_exp = observed_data.get_val('sint_d')

        int_u_mod, int_d_mod = self.calc_profile(tth)
        sint_sum_exp = (sint_u_exp**2 + sint_d_exp**2)**0.5

        chi_sq_u = ((int_u_mod-int_u_exp)/sint_u_exp)**2
        chi_sq_d = ((int_d_mod-int_d_exp)/sint_d_exp)**2

        chi_sq_sum = ((int_u_mod+int_d_mod-int_u_exp-int_d_exp)/sint_sum_exp)**2
        chi_sq_dif = ((int_u_mod-int_d_mod-int_u_exp+int_d_exp)/sint_sum_exp)**2
        return chi_sq_u.sum(), chi_sq_d.sum(), chi_sq_sum.sum(), chi_sq_dif.sum()


class ExperimentPowder2D(dict):
    """
    Class to describe all information concerning to the experiment
    """
    def __init__(self, setup = SetupPowder2D(), 
                 list_calculated_data = [], 
                 observed_data = ObservedDataPowder2D()):
        super(ExperimentPowder2D, self).__init__()
        self._p_setup = None
        self._list_calculated_data = []
        self._p_observed_data = None
        
        self._refresh(setup, observed_data)

    def __repr__(self):
        lsout = """Experiment 2D:\n{:}\n{:} """.format(
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
    
    def calc_profile(self, tth, phi):
        """
        calculate intensity for the given diffraction angle
        """
        setup = self._p_setup
        background = setup.calc_background(tth, phi)
        wavelength = setup.get_val("wavelength")
        beam_polarization = setup.get_val("beam_polarization")
        
        tth_min = tth.min()
        tth_max = tth.max()
        sthovl_min = numpy.sin(0.5*tth_min*numpy.pi/180.)/wavelength
        sthovl_max = numpy.sin(0.5*tth_max*numpy.pi/180.)/wavelength

        res_u_1d = numpy.zeros(tth.shape[0], dtype=float)
        res_d_1d = numpy.zeros(tth.shape[0], dtype=float)

        for calculated_data in self._list_calculated_data:
            
            i_g = calculated_data.get_val("crystal").get_val("i_g")
            cell = calculated_data.get_val("crystal").get_val("cell")
            space_groupe = calculated_data.get_val("crystal").get_val("space_groupe")
            h, k, l, mult = setup.calc_hkl(cell, space_groupe, sthovl_min, sthovl_max)
            
            n_sq, m_sq, cross = calculated_data.calc_for_iint(h, k, l)
             
            np_iint_u, np_iint_d
            sthovl_hkl = cell.calc_sthovl(h, k, l)
            
            tth_hkl_rad = 2.*numpy.arcsin(sthovl_hkl*wavelength)
            tth_hkl = tth_hkl_rad*180./numpy.pi
            profile_2d = setup.calc_profile(tth, phi, tth_hkl, i_g)
            
            
            np_iint_u_2d = numpy.meshgrid(tth, np_iint_u*mult, indexing="ij")[1]

            np_iint_d_2d = numpy.meshgrid(tth, np_iint_d*mult, indexing="ij")[1]
            
            res_u_2d = profile_2d*np_iint_u_2d 
            res_d_2d = profile_2d*np_iint_d_2d 

            res_u_1d += res_u_2d.sum(axis=1) 
            res_d_1d += res_d_2d.sum(axis=1) 
            
        return res_u_1d+background, res_d_1d+background
    
    def calc_chi_sq(self):
        """
        calculate chi square
        """
        observed_data = self._p_observed_data

        tth = observed_data.get_val('tth')
        int_u_exp = observed_data.get_val('int_u')
        sint_u_exp = observed_data.get_val('sint_u')
        int_d_exp = observed_data.get_val('int_d')
        sint_d_exp = observed_data.get_val('sint_d')

        int_u_mod, int_d_mod = self.calc_profile(tth)
        sint_sum_exp = (sint_u_exp**2 + sint_d_exp**2)**0.5

        chi_sq_u = ((int_u_mod-int_u_exp)/sint_u_exp)**2
        chi_sq_d = ((int_d_mod-int_d_exp)/sint_d_exp)**2

        chi_sq_sum = ((int_u_mod+int_d_mod-int_u_exp-int_d_exp)/sint_sum_exp)**2
        chi_sq_dif = ((int_u_mod-int_d_mod-int_u_exp+int_d_exp)/sint_sum_exp)**2
        return chi_sq_u.sum(), chi_sq_d.sum(), chi_sq_sum.sum(), chi_sq_dif.sum()
    
    
if (__name__ == "__main__"):
  pass
