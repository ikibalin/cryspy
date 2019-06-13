
"""
define classe to describe experiment
"""
__author__ = 'ikibalin'
__version__ = "2019_06_07"
import os
import numpy
import sys

import matplotlib
import matplotlib.pyplot  



from cl_observed_data_single_domain import *
from cl_calculated_data_single import *
from cl_setup_single import *

from cl_variable import *




class ExperimentSingleDomain(dict):
    """
    Class to describe all information concerning to the experiment for monocrystal
    """
    def __init__(self, name=None, setup = SetupSingle(), 
                 list_calculated_data = [], 
                 observed_data = ObservedDataSingleDomain(), file_out=None, file_dir=None, scale_domain=None):
        super(ExperimentSingleDomain, self).__init__()
        self._p_name = None
        self._p_setup = None
        self._list_calculated_data = []
        self._p_observed_data = None
        self._p_file_out = None
        self._p_file_dir = None
        self._p_scale_domain = []
        self._refresh(name, setup, observed_data, file_out, file_dir, scale_domain)

    def __repr__(self):
        ls_out = """ExperimentSingleDomain:\n name: {:}\n file_out: {:}\n{:}\n{:}""".format(
                self._p_name, self._p_file_out, self._p_setup, self._p_observed_data)
        
        ls_calculated_data = []
        for calculated_data in self._list_calculated_data:
            ls_calculated_data.append("{:}".format(calculated_data))

        ls_out += "\n\n\nCalculatedData:\n\n"+"\n\n".join(ls_calculated_data)
        return ls_out

    def _refresh(self, name, setup, observed_data, file_out, file_dir, scale_domain):
        if name is not None:
            self._p_name = name
        if setup is not None:
            self._p_setup = setup
        if observed_data is not None:
            self._p_observed_data = observed_data
        if file_out is not None:
            self._p_file_out = file_out
        if file_dir is not None:
            self._p_file_dir = file_dir
        if scale_domain is not None:
            self._p_scale_domain = [hh for hh in scale_domain]
            
    def set_val(self, name=None, setup=None, observed_data=None, file_out=None, file_dir=None, scale_domain=None):
        self._refresh(name, setup, observed_data, file_out, file_dir, scale_domain)
        
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
name is the name of experiment (should be unique)
setup is to describe parameters of diffractometer 
observed_data is the experimental data
file_out is the file name of model data
        """
        print(lsout)
    def add_calculated_data(self, observed_data):
        self._list_calculated_data.append(observed_data)

    def del_calculated_data(self, ind):
        self._list_calculated_data.pop(ind)        

    def replace_calculated_data(self, ind, observed_data):
        self._list_calculated_data.pop(ind)
        self._list_calculated_data.insert(ind, observed_data)
    
    def calc_iint_u_d_flip_ratio(self, h, k, l, l_crystal):
        """
        calculate intensity for the given diffraction angle
        """
        setup = self._p_setup
        wave_length = setup.get_val("wave_length")
        beam_polarization = setup.get_val("beam_polarization")
        
        #calculations only for first crystal phase
        for calculated_data in self._list_calculated_data[0:1]:
            ind_cry = None
            observed_data_name = calculated_data.get_val("name")
            for i_crystal, crystal in enumerate(l_crystal):
                if crystal.get_val("name") == observed_data_name:
                    ind_cry = i_crystal
                    break
            if ind_cry is None:
                print("Crystal with name '{:}' is not found.".format(
                    observed_data_name))
                return
            crystal = l_crystal[ind_cry]
            
            iint_u, iint_d, flip_ratio, d_info_cd = calculated_data.calc_iint_u_d_flip_ratio(
                              h, k, l, beam_polarization, wave_length, crystal)

        d_info_out = {}
        d_info_out.update(d_info_cd)        
            
        return iint_u, iint_d, flip_ratio, d_info_out
    
    def calc_chi_sq(self, l_crystal, d_map={}):
        """
        calculate chi square
        """
        if d_map == {}:
            d_map.update(self.plot_map())
        #if not(d_map["flag"]|(d_map["out"] is None)):
        #    chi_sq_val, n = d_map["out"]
        #    return chi_sq_val, n        
        observed_data = self._p_observed_data

        np_h = observed_data.get_val('h')
        np_k = observed_data.get_val('k')
        np_l = observed_data.get_val('l')
        flip_ratio_exp = observed_data.get_val('flip_ratio')
        sflip_ratio_exp = observed_data.get_val('sflip_ratio')

        wave_length = observed_data.get_val('wave_length')
        setup = self._p_setup
        setup.set_val(wave_length=wave_length)

        field = observed_data.get_val('field')
        np_orientation = observed_data.get_val('orientation')


        n_domain = np_h.shape[-1]
        for i_domain in range(n_domain):
            h, k, l = np_h[:, i_domain], np_k[:, i_domain], np_l[:, i_domain]
            orientation = np_orientation[:, :, i_domain]
            
            for calculated_data in self._list_calculated_data:
                calculated_data.set_val(field=field, orientation=orientation)
            
            
            iint_u_mod, iint_d_mod, flip_ratio_mod, d_info_cd = self.calc_iint_u_d_flip_ratio(
                                               h, k, l, l_crystal)
            
            if i_domain == 0:
                np_iint_u = iint_u_mod[:, numpy.newaxis]
                np_iint_d = iint_d_mod[:, numpy.newaxis]
                np_flip_ratio = flip_ratio_mod[:, numpy.newaxis]
            else:
                np_iint_u = numpy.concatenate((np_iint_u, iint_u_mod[:, numpy.newaxis]), axis=1)
                np_iint_d = numpy.concatenate((np_iint_d, iint_d_mod[:, numpy.newaxis]), axis=1)
                np_flip_ratio = numpy.concatenate((np_flip_ratio, flip_ratio_mod[:, numpy.newaxis]), axis=1)
        
        
        l_scale_domain = [1.*hh for hh in self._p_scale_domain]
        if len(l_scale_domain) < n_domain:
            l_scale_domain.extend([1. for hh in range(n_domain-len(l_scale_domain))])
        np_scale_domain = numpy.array(l_scale_domain, dtype=float)
        l_flip_ratio_mod = []
        for iint_u_domain, iint_d_domain in zip(np_iint_u, np_iint_d):
            flag = numpy.logical_not(numpy.isnan(iint_u_domain))
            iint_u = (np_scale_domain[flag]*iint_u_domain[flag]).sum()
            iint_d = (np_scale_domain[flag]*iint_d_domain[flag]).sum()
            if all(list(numpy.logical_not(flag))):
                l_flip_ratio_mod.append(1.)
            else:
                l_flip_ratio_mod.append(iint_u/iint_d)

        flip_ratio_mod = numpy.array(l_flip_ratio_mod)
        chi_sq = ((flip_ratio_mod-flip_ratio_exp)/sflip_ratio_exp)**2
        
        chi_sq_val = (chi_sq[numpy.logical_not(numpy.isnan(chi_sq))]).sum()
        
        n = numpy.logical_not(numpy.isnan(chi_sq)).sum()

        d_info_out = {"chi_sq_val": chi_sq_val, "n": n}
        d_info_out.update(d_info_cd)        
        return chi_sq_val, n, d_info_out
    
    def plot_map(self):
        b_variable = self.is_variable()
        d_map = {"flag": b_variable, "out":None}
        return d_map
    
    def is_variable_iint(self):
        lres = []
        for calculated_data in self._list_calculated_data:
            lres.append(calculated_data.is_variable())
        return lres

    def is_variable_setup(self):
        setup = self.get_val("setup")
        res = setup.is_variable()
        return res

    def is_variable(self):
        """
        without extinction
        """
        lres = self.is_variable_iint() 
        lres.append(self.is_variable_setup())
        res = any(lres) 
        return res 
    
    def get_variables(self):
        l_variable = []
        setup = self.get_val("setup")
        l_var = setup.get_variables()
        l_variable.extend(l_var)
        for calculated_data in self._list_calculated_data:
            l_var = calculated_data.get_variables()
            l_variable.extend(l_var)
        return l_variable

    def save_exp_mod_data(self, l_crystal):
        observed_data = self.get_val("observed_data")
        np_h = observed_data.get_val("h")
        np_k = observed_data.get_val("k")
        np_l = observed_data.get_val("l")
        fr_exp = observed_data.get_val("flip_ratio")
        sfr_exp = observed_data.get_val("sflip_ratio")
        
        h, k, l = np_h[:,0], np_k[:,0], np_l[:,0]
        iint_u_mod, iint_d_mod, fr_mod, d_info = self.calc_iint_u_d_flip_ratio(h, k, l, l_crystal)
        
        s_1 = "   h   k   l       FR_exp        sigma       FR_mod  iint_up_mod iint_down_mod\n"
        l_s_2 = ["{:4}{:4}{:4} {:12.5f} {:12.5f} {:12.5f} {:12.5f} {:12.5f}".format(
                hh_1, hh_2, hh_3, hh_4, hh_5, hh_6, hh_7, hh_8) for 
                hh_1, hh_2, hh_3, hh_4, hh_5, hh_6, hh_7, hh_8 in 
                zip(h, k, l, fr_exp, sfr_exp, fr_mod, iint_u_mod, iint_d_mod)]

        s_int = s_1 + "\n".join(l_s_2)
        #hkl should be added
        s_out = s_int
        
        file_out, file_dir = self._p_file_out, self._p_file_dir
        if ((file_out is None) | (file_dir is None)):
            print("File to save model data is not defined.")
            return
        f_name = os.path.join(file_dir, file_out)
        fid = open(f_name, "w")
        fid.write(s_out)
        fid.close()

    def print_report(self, l_crystal):
        s_out = "{:}".format(self)
        
        chi_sq_val, n, d_info = self.calc_chi_sq(l_crystal)
        
        ls_out = []
        ls_out.append("\n\n   h   k   l   fn_real   fn_imag   iint_up iint_down flip_ratio")

        h, k, l, f_nucl = d_info["h"], d_info["k"], d_info["l"], d_info["f_nucl"]
        flip_ratio = d_info["flip_ratio"]
        iint_u, iint_d = d_info["iint_u"], d_info["iint_d"]
        sft_11, sft_12, sft_13 = d_info["sft_11"], d_info["sft_12"], d_info["sft_13"]
        sft_21, sft_22, sft_23 = d_info["sft_21"], d_info["sft_22"], d_info["sft_23"]
        sft_31, sft_32, sft_33 = d_info["sft_31"], d_info["sft_32"], d_info["sft_33"]

        for h1, k1, l1, f, i_u, i_d, f_r in zip(
            h, k, l, f_nucl, iint_u, iint_d, flip_ratio):
                ls_out.append(""" {:3} {:3} {:3}  {:8.3f}  {:8.3f}  {:8.1f}  {:8.1f} {:8.5f}""".format(
     h1, k1, l1, f.real, f.imag, i_u, i_d, f_r))
        ls_out.append("\n\nStructure factor tensor in a*, b*, c* (real part, mu_B) ")
        ls_out.append("\n   h   k   l   sft_11   sft_22   sft_33    sft_12   sft_13   sft_23")
        for h1, k1, l1, c11, c22,  c33, c12, c13, c23 in zip(
            h, k, l, sft_11, sft_22, sft_33, sft_12, sft_13, sft_23):
                ls_out.append(""" {:3} {:3} {:3} {:8.3f} {:8.3f} {:8.3f}  {:8.3f} {:8.3f} {:8.3f}""".format(
                        h1, k1, l1, c11.real, c22.real, c33.real, c12.real, c13.real, c23.real))
                
        return s_out+"\n".join(ls_out)


    
if (__name__ == "__main__"):
  pass
