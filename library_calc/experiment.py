
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
from setup_powder_2d import *

from variable import *


def save_2d(np_tth, np_phi, int_2d):
    line = "       {:7}    ".format(len(np_phi)) + " ".join(["{:10.2f}".format(hh) for hh in np_tth])
    lsout_u = [line]
    for int_1d, phi in zip(int_2d.transpose(), np_phi):
        line_u = "    {:10.5f}    ".format(phi) + " ".join(["{:10.2f}".format(hh) for hh in int_1d])
        lsout_u.append(line_u)

    return "\n".join(lsout_u)

class ExperimentSingle(dict):
    """
    Class to describe all information concerning to the experiment for monocrystal
    """
    def __init__(self, name=None, setup = SetupSingle(), 
                 list_calculated_data = [], 
                 observed_data = ObservedDataSingle(), f_out=None):
        super(ExperimentSingle, self).__init__()
        self._p_name = None
        self._p_setup = None
        self._list_calculated_data = []
        self._p_observed_data = None
        self._p_f_out = None
        self._refresh(name, setup, observed_data, f_out)

    def __repr__(self):
        ls_out = """ExperimentSingle:\n name: {:}\n f_out: {:}\n{:}\n{:}""".format(
                self._p_name, self._p_setup, self._p_observed_data, 
                self._p_f_out)
        
        ls_calculated_data = []
        for calculated_data in self._list_calculated_data:
            ls_calculated_data.append("{:}".format(calculated_data))

        ls_out += "\n\n\nCalculatedData:\n\n"+"\n\n".join(ls_calculated_data)
        return ls_out

    def _refresh(self, name, setup, observed_data, f_out):
        if name is not None:
            self._p_name = name
        if setup is not None:
            self._p_setup = setup
        if observed_data is not None:
            self._p_observed_data = observed_data
        if f_out is not None:
            self._p_f_out = f_out
            
    def set_val(self, name=None, setup=None, observed_data=None, f_out=None):
        self._refresh(name, setup, observed_data, f_out)
        
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
f_out is the file name of model data
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
            
            iint_u, iint_d, flip_ratio = calculated_data.calc_iint_u_d_flip_ratio(
                              h, k, l, beam_polarization, wave_length, crystal)
            
        return iint_u, iint_d, flip_ratio 
    
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

        h = observed_data.get_val('h')
        k = observed_data.get_val('k')
        l = observed_data.get_val('l')
        flip_ratio_exp = observed_data.get_val('flip_ratio')
        sflip_ratio_exp = observed_data.get_val('sflip_ratio')

        wave_length = observed_data.get_val('wave_length')
        setup = self._p_setup
        setup.set_val(wave_length=wave_length)

        field = observed_data.get_val('field')
        orientation = observed_data.get_val('orientation')
        for calculated_data in self._list_calculated_data:
            calculated_data.set_val(field=field, orientation=orientation)

        int_u_mod, int_d_mod, flip_ratio_mod = self.calc_iint_u_d_flip_ratio(
                                               h, k, l, l_crystal)
        

        chi_sq = ((flip_ratio_mod-flip_ratio_exp)/sflip_ratio_exp)**2
        chi_sq_val = (chi_sq[numpy.logical_not(numpy.isnan(chi_sq))]).sum()
        n = numpy.logical_not(numpy.isnan(chi_sq)).sum()
        d_map["out"] = chi_sq_val, n
        return chi_sq_val, n
    
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
        h = observed_data.get_val("h")
        k = observed_data.get_val("k")
        l = observed_data.get_val("l")
        fr_exp = observed_data.get_val("flip_ratio")
        sfr_exp = observed_data.get_val("sflip_ratio")
        
        
        iint_u_mod, iint_d_mod, fr_mod = self.calc_iint_u_d_flip_ratio(h, k, l, l_crystal)
        
        s_1 = "   h   k   l       FR_exp        sigma       FR_mod  iint_up_mod iint_down_mod\n"
        l_s_2 = ["{:4}{:4}{:4} {:12.5f} {:12.5f} {:12.5f} {:12.5f} {:12.5f}".format(
                hh_1, hh_2, hh_3, hh_4, hh_5, hh_6, hh_7, hh_8) for 
                hh_1, hh_2, hh_3, hh_4, hh_5, hh_6, hh_7, hh_8 in 
                zip(h, k, l, fr_exp, sfr_exp, fr_mod, iint_u_mod, iint_d_mod)]

        s_int = s_1 + "\n".join(l_s_2)
        #hkl should be added
        s_out = s_int
        
        f_out = self._p_f_out
        if f_out is None:
            print("File to save model data is not defined.")
            return
        fid = open(f_out, "w")
        fid.write(s_out)
        fid.close()



class ExperimentPowder1D(dict):
    """
    Class to describe all information concerning to the experiment for powder 1d
    """
    def __init__(self, name=None, setup=SetupPowder1D(), 
                 list_calculated_data=[], 
                 observed_data=ObservedDataPowder1D(), mode_chi_sq="up, down", 
                 f_out=None):
        super(ExperimentPowder1D, self).__init__()
        self._p_name = None
        self._p_setup = None
        self._list_calculated_data = []
        self._p_observed_data = None
        self._p_mode_chi_sq = None
        self._p_f_out = None
        
        self._refresh(name, setup, observed_data, mode_chi_sq, f_out)

    def __repr__(self):
        ls_out = """ExperimentPowder1D:\nname: {:}
 mode_chi_sq: {:}\n f_out: {:}\n{:}\n{:}""".format(self._p_name, 
 self._p_mode_chi_sq, self._p_f_out, self._p_setup, self._p_observed_data)
        
        ls_calculated_data = []
        for calculated_data in self._list_calculated_data:
            ls_calculated_data.append("{:}".format(calculated_data))

        ls_out += "\n\n\nCalculatedData:\n\n"+"\n\n".join(ls_calculated_data)
        return ls_out

    def _refresh(self, name, setup, observed_data, mode_chi_sq, f_out):
        if name is not None:
            self._p_name = name
        if setup is not None:
            self._p_setup = setup
        if observed_data is not None:
            self._p_observed_data = observed_data
        if mode_chi_sq is not None:
            self._p_mode_chi_sq = mode_chi_sq
        if f_out is not None:
            self._p_f_out = f_out
            
    def set_val(self, name=None, setup=None, observed_data=None, 
                mode_chi_sq=None, f_out=None):
        self._refresh(name, setup, observed_data, mode_chi_sq, f_out)
        
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
name is the name of experiment
setup is to describe parameters of diffractometer 
observed_data is the experimental data
mode_chi_sq is define the refinement: "up down sum diff"
f_out is the file name of experimental and model data
        """
        print(lsout)
    def add_calculated_data(self, observed_data):
        self._list_calculated_data.append(observed_data)

    def del_calculated_data(self, ind):
        self._list_calculated_data.pop(ind)        

    def replace_calculated_data(self, ind, observed_data):
        self._list_calculated_data.pop(ind)
        self._list_calculated_data.insert(ind, observed_data)
    
    def calc_profile(self, tth, l_crystal):
        """
        calculate intensity for the given diffraction angle
        """
        setup = self._p_setup
        background = setup.get_val("background")
        int_bkgd = background.interpolate_by_points(tth)
        wave_length = setup.get_val("wave_length")
        beam_polarization = setup.get_val("beam_polarization")
        
        tth_min = tth.min()
        tth_max = tth.max()
        sthovl_min = numpy.sin(0.5*tth_min*numpy.pi/180.)/wave_length
        sthovl_max = numpy.sin(0.5*tth_max*numpy.pi/180.)/wave_length

        res_u_1d = numpy.zeros(tth.shape[0], dtype=float)
        res_d_1d = numpy.zeros(tth.shape[0], dtype=float)
        
        for calculated_data in self._list_calculated_data:
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

            scale = 1.*calculated_data.get_val("scale")
            
            i_g = 1.*crystal.get_val("i_g")
            cell = crystal.get_val("cell")
            space_groupe = crystal.get_val("space_groupe")
            h, k, l, mult = setup.calc_hkl(cell, space_groupe, sthovl_min, sthovl_max)
            
            np_iint_u, np_iint_d = calculated_data.calc_iint(h, k, l, 
                                                    beam_polarization, crystal)
            sthovl_hkl = cell.calc_sthovl(h, k, l)
            
            tth_hkl_rad = 2.*numpy.arcsin(sthovl_hkl*wave_length)
            tth_hkl = tth_hkl_rad*180./numpy.pi
            profile_2d = setup.calc_profile(tth, tth_hkl, i_g)
            
            
            np_iint_u_2d = numpy.meshgrid(tth, np_iint_u*mult, indexing="ij")[1]

            np_iint_d_2d = numpy.meshgrid(tth, np_iint_d*mult, indexing="ij")[1]
            
            res_u_2d = profile_2d*np_iint_u_2d 
            res_d_2d = profile_2d*np_iint_d_2d 

            res_u_1d += scale*res_u_2d.sum(axis=1) 
            res_d_1d += scale*res_d_2d.sum(axis=1) 
            
        return res_u_1d+int_bkgd, res_d_1d+int_bkgd
    
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

        tth = observed_data.get_val('tth')
        int_u_exp = observed_data.get_val('int_u')
        sint_u_exp = observed_data.get_val('sint_u')
        int_d_exp = observed_data.get_val('int_d')
        sint_d_exp = observed_data.get_val('sint_d')
        
        wave_length = observed_data.get_val('wave_length')
        setup = self._p_setup
        setup.set_val(wave_length=wave_length)

        field = observed_data.get_val('field')
        for calculated_data in self._list_calculated_data:
            calculated_data.set_val(field=field)

        int_u_mod, int_d_mod = self.calc_profile(tth, l_crystal)
        sint_sum_exp = (sint_u_exp**2 + sint_d_exp**2)**0.5

        chi_sq_u = ((int_u_mod-int_u_exp)/sint_u_exp)**2
        chi_sq_d = ((int_d_mod-int_d_exp)/sint_d_exp)**2

        chi_sq_sum = ((int_u_mod+int_d_mod-int_u_exp-int_d_exp)/sint_sum_exp)**2
        chi_sq_dif = ((int_u_mod-int_d_mod-int_u_exp+int_d_exp)/sint_sum_exp)**2

        
        chi_sq_u_val = (chi_sq_u[numpy.logical_not(numpy.isnan(chi_sq_u))]).sum()
        n_u = numpy.logical_not(numpy.isnan(chi_sq_u)).sum()
        
        chi_sq_d_val = (chi_sq_d[numpy.logical_not(numpy.isnan(chi_sq_d))]).sum()
        n_d = numpy.logical_not(numpy.isnan(chi_sq_d)).sum()

        chi_sq_sum_val = (chi_sq_sum[numpy.logical_not(numpy.isnan(chi_sq_sum))]).sum()
        n_sum = numpy.logical_not(numpy.isnan(chi_sq_sum)).sum()

        chi_sq_dif_val = (chi_sq_dif[numpy.logical_not(numpy.isnan(chi_sq_dif))]).sum()
        n_sum = numpy.logical_not(numpy.isnan(chi_sq_dif)).sum()

        mode_chi_sq = self._p_mode_chi_sq
        flag_u = "up" in mode_chi_sq 
        flag_d = "down" in mode_chi_sq 
        flag_sum = "sum" in mode_chi_sq 
        flag_dif = "up" in mode_chi_sq 
        
        chi_sq_val = (int(flag_u)*chi_sq_u_val + int(flag_d)*chi_sq_d_val + 
                  int(flag_sum)*chi_sq_sum_val + int(flag_dif)*chi_sq_dif_val)
        n = (int(flag_u)*n_u + int(flag_d)*n_d + int(flag_sum)*n_sum + 
             int(flag_dif)*n_sum)
        d_map["out"] = chi_sq_val, n
        return chi_sq_val, n
    
    def plot_map(self):
        b_variable = self.is_variable()
        d_map = {"flag": b_variable, "out":None}
        return d_map
    
    def is_variable_iint(self):
        lres = []
        for calculated_data in self._list_calculated_data:
            lres.append(calculated_data.is_variable())
        return lres

    def is_variable_profile_s(self):
        setup = self.get_val("setup")
        res_s = setup.is_variable()
        lres = [res_s]
        return lres

    def is_variable_profile(self):
        lres = self.is_variable_iint() 
        lres.extend(self.is_variable_profile_s())
        res = any(lres) 
        return res

    def is_variable(self):
        """
        without extinction
        """
        res = any([self.is_variable_profile()])
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
        tth = observed_data.get_val("tth")
        int_u_exp = observed_data.get_val("int_u")
        sint_u_exp = observed_data.get_val("sint_u")
        int_d_exp = observed_data.get_val("int_d")
        sint_d_exp = observed_data.get_val("sint_d")
        
        int_u_mod, int_d_mod = self.calc_profile(tth, l_crystal) 
        
        s_1 = "    ttheta       exp_up        sigma       mod_up     exp_down        sigma     mod_down\n"
        l_s_2 = ["{:10.2f} {:12.5f} {:12.5f} {:12.5f} {:12.5f} {:12.5f} {:12.5f}".format(
                hh_1, hh_2, hh_3, hh_4, hh_5, hh_6, hh_7) for 
                hh_1, hh_2, hh_3, hh_4, hh_5, hh_6, hh_7 in 
                zip(tth, int_u_exp, sint_u_exp, int_u_mod, 
                         int_d_exp, sint_d_exp, int_d_mod)]

        s_int = s_1 + "\n".join(l_s_2)
        #hkl should be added
        s_out = s_int
        
        f_out = self._p_f_out
        if f_out is None:
            print("File to save model data is not defined.")
            return
        fid = open(f_out, "w")
        fid.write(s_out)
        fid.close()


class ExperimentPowder2D(dict):
    """
    Class to describe all information concerning to the experiment for powder 2d
    """
    def __init__(self, name=None, setup=SetupPowder2D(), 
                 list_calculated_data=[], 
                 observed_data=ObservedDataPowder2D(), mode_chi_sq="up, down", 
                 f_out=None):
        super(ExperimentPowder2D, self).__init__()
        self._p_name = None
        self._p_setup = None
        self._list_calculated_data = []
        self._p_observed_data = None
        self._p_mode_chi_sq = None
        self._p_f_out = None
        
        self._refresh(name, setup, observed_data, mode_chi_sq, f_out)

    def __repr__(self):
        ls_out = """ExperimentPowder2D:\n name: {:}
 mode_chi_sq: {:}\n f_out: {:}\n{:}\n{:}""".format(self._p_name, 
     self._p_mode_chi_sq, self._p_f_out, self._p_setup, self._p_observed_data)
        
        ls_calculated_data = []
        for calculated_data in self._list_calculated_data:
            ls_calculated_data.append("{:}".format(calculated_data))

        ls_out += "\n\n\nCalculatedData:\n\n"+"\n\n".join(ls_calculated_data)
            
        return ls_out

    def _refresh(self, name, setup, observed_data, mode_chi_sq, f_out):
        if name is not None:
            self._p_name = name
        if setup is not None:
            self._p_setup = setup
        if observed_data is not None:
            self._p_observed_data = observed_data
        if mode_chi_sq is not None:
            self._p_mode_chi_sq = mode_chi_sq
        if f_out is not None:
            self._p_f_out = f_out

            
    def set_val(self, name=None, setup=None, observed_data=None, 
                mode_chi_sq=None, f_out=None):
        self._refresh(name, setup, observed_data, mode_chi_sq, f_out)
        
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
name is the name of an experiment (should be unique)
setup is to describe parameters of diffractometer 
observed_data is the experimental data
mode_chi_sq is define the refinement: "up down sum diff"
f_out is the file name to save model profiles
        """
        print(lsout)
    def add_calculated_data(self, observed_data):
        self._list_calculated_data.append(observed_data)

    def del_calculated_data(self, ind):
        self._list_calculated_data.pop(ind)        

    def replace_calculated_data(self, ind, observed_data):
        self._list_calculated_data.pop(ind)
        self._list_calculated_data.insert(ind, observed_data)
    
    def calc_profile(self, tth, phi, l_crystal, d_map={}):
        """
        calculate intensity for the given diffraction angle
        
        tth and phi is 1D data
        """
        if d_map == {}:
            d_map.update(self.plot_map())
        d_profile = d_map["profile"]
        #if not(d_profile["flag"]|(d_profile["out"] is None)):
        #    res_u_2d_int_bkgd, res_d_2d_int_bkgd = d_profile["out"]
        #    return res_u_2d_int_bkgd, res_d_2d_int_bkgd
        
        tth_rad = tth*numpy.pi/180.
        phi_rad = phi*numpy.pi/180.
        
        setup = self._p_setup
        background = setup.get_val("background")
        int_bkgd = background.interpolate_by_points(tth, phi)
        
        wave_length = setup.get_val("wave_length")
        beam_polarization = setup.get_val("beam_polarization")

        p_u = beam_polarization.get_val("p_u")
        p_d = beam_polarization.get_val("p_d")
        
        tth_min = tth.min()
        tth_max = tth.max()
        sthovl_min = numpy.sin(0.5*tth_min*numpy.pi/180.)/wave_length
        sthovl_max = numpy.sin(0.5*tth_max*numpy.pi/180.)/wave_length

        res_u_2d = numpy.zeros((tth.shape[0],phi.shape[0]), dtype=float)
        res_d_2d = numpy.zeros((tth.shape[0],phi.shape[0]), dtype=float)

        for calculated_data in self._list_calculated_data:
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

            scale = 1.*calculated_data.get_val("scale")
            i_g = 1.*crystal.get_val("i_g")
            cell = crystal.get_val("cell")
            space_groupe = crystal.get_val("space_groupe")
            
            
            d_hkl = d_map[("hkl", observed_data_name)]
            if not(d_hkl["flag"]|(d_hkl["out"] is None)):
                h, k, l, mult = d_hkl["out"]
            else:
                h, k, l, mult = setup.calc_hkl(cell, space_groupe, sthovl_min, sthovl_max)
                d_hkl["out"] = (h, k, l, mult)
            mult_3d = numpy.meshgrid(tth, phi, mult, indexing="ij")[2]
            
            
            #d_for_iint = d_map[("for_iint", observed_data_name)]
            #if not(d_for_iint["flag"]|(d_for_iint["out"] is None)):
            #    f_nucl_sq, f_m_p_sin_sq, f_m_p_cos_sq, cross_sin = d_for_iint["out"]
            #else:
            #    #dimension: (tth_hkl)
            f_nucl_sq, f_m_p_sin_sq, f_m_p_cos_sq, cross_sin = calculated_data.calc_for_iint(
                        h, k, l, crystal)#d_for_iint
            
            #print("   h   k   l mult   f_nucl_sq f_m_p_sin_sq f_m_p_cos_sq cross_sin")
            #for h_1, k_1, l_1, hh_1, hh_2, hh_3, hh_4, hh_5  in zip(h, k, l, mult, f_nucl_sq, f_m_p_sin_sq, f_m_p_cos_sq, cross_sin):
            #    print(""" {:3} {:3} {:3} {:4} {:11.3f} {:11.3f} {:11.3f} {:11.3f}""".format(
            #            h_1, k_1, l_1, hh_1, hh_2, hh_3, hh_4, hh_5))   
             
            tth_r_3d, phi_r_3d, f_n_3d = numpy.meshgrid(tth_rad, phi_rad, f_nucl_sq, indexing="ij")
            f_m_s_3d = numpy.meshgrid(tth, phi, f_m_p_sin_sq, indexing="ij")[2]
            f_m_c_3d = numpy.meshgrid(tth, phi, f_m_p_cos_sq, indexing="ij")[2]
            f_c_s_3d = numpy.meshgrid(tth, phi, cross_sin, indexing="ij")[2]
            c_a_sq_3d = (numpy.cos(0.5*tth_r_3d)*numpy.sin(phi_r_3d))**2
            s_a_sq_3d = 1.-c_a_sq_3d

            i_u_3d = mult_3d*(f_n_3d +(f_m_s_3d+p_u*f_c_s_3d)*s_a_sq_3d+f_m_c_3d*c_a_sq_3d)
            i_d_3d = mult_3d*(f_n_3d +(f_m_s_3d-p_d*f_c_s_3d)*s_a_sq_3d+f_m_c_3d*c_a_sq_3d)

            sthovl_hkl = cell.calc_sthovl(h, k, l)
            tth_hkl_rad = 2.*numpy.arcsin(sthovl_hkl*wave_length)
            tth_hkl = tth_hkl_rad*180./numpy.pi
            

            d_profile_s = d_map[("profile", observed_data_name)]
            if not(d_profile_s["flag"]|(d_profile_s["out"] is None)):
                profile_3d = d_profile_s["out"]
            else:
                #dimension: (tth, phi, tth_hkl)
                profile_3d = setup.calc_profile(tth, phi, tth_hkl, i_g)
                d_profile_s["out"] = profile_3d
            
            
            res_u_3d = profile_3d*i_u_3d
            res_d_3d = profile_3d*i_d_3d

            res_u_2d += scale*res_u_3d.sum(axis=2) 
            res_d_2d += scale*res_d_3d.sum(axis=2) 
        d_profile["out"] = res_u_2d+int_bkgd, res_d_2d+int_bkgd
        return res_u_2d+int_bkgd, res_d_2d+int_bkgd
    
    def calc_chi_sq(self, l_crystal, d_map={}):
        """
        calculate chi square
        """
        if d_map == {}:
            d_map.update(self.plot_map())
        if not(d_map["flag"]|(d_map["out"] is None)):
            chi_sq_val, n = d_map["out"]
            return chi_sq_val, n
        
        observed_data = self._p_observed_data

        tth = observed_data.get_val('tth')
        phi = observed_data.get_val('phi')
        int_u_exp = observed_data.get_val('int_u')
        sint_u_exp = observed_data.get_val('sint_u')
        int_d_exp = observed_data.get_val('int_d')
        sint_d_exp = observed_data.get_val('sint_d')
        
        wave_length = observed_data.get_val('wave_length')
        setup = self._p_setup
        setup.set_val(wave_length=wave_length)

        field = observed_data.get_val('field')
        for calculated_data in self._list_calculated_data:
            calculated_data.set_val(field=field)

        
            
            
        int_u_mod, int_d_mod = self.calc_profile(tth, phi, l_crystal, d_map)
        sint_sum_exp = (sint_u_exp**2 + sint_d_exp**2)**0.5

        chi_sq_u = ((int_u_mod-int_u_exp)/sint_u_exp)**2
        chi_sq_d = ((int_d_mod-int_d_exp)/sint_d_exp)**2

        chi_sq_sum = ((int_u_mod+int_d_mod-int_u_exp-int_d_exp)/sint_sum_exp)**2
        chi_sq_dif = ((int_u_mod-int_d_mod-int_u_exp+int_d_exp)/sint_sum_exp)**2

        chi_sq_u_val = (chi_sq_u[numpy.logical_not(numpy.isnan(chi_sq_u))]).sum()
        n_u = numpy.logical_not(numpy.isnan(chi_sq_u)).sum()
        
        chi_sq_d_val = (chi_sq_d[numpy.logical_not(numpy.isnan(chi_sq_d))]).sum()
        n_d = numpy.logical_not(numpy.isnan(chi_sq_d)).sum()

        chi_sq_sum_val = (chi_sq_sum[numpy.logical_not(numpy.isnan(chi_sq_sum))]).sum()
        n_sum = numpy.logical_not(numpy.isnan(chi_sq_sum)).sum()

        chi_sq_dif_val = (chi_sq_dif[numpy.logical_not(numpy.isnan(chi_sq_dif))]).sum()
        n_sum = numpy.logical_not(numpy.isnan(chi_sq_dif)).sum()
        
        mode_chi_sq = self._p_mode_chi_sq
        flag_u = "up" in mode_chi_sq 
        flag_d = "down" in mode_chi_sq 
        flag_sum = "sum" in mode_chi_sq 
        flag_dif = "diff" in mode_chi_sq 
        
        chi_sq_val = (int(flag_u)*chi_sq_u_val + int(flag_d)*chi_sq_d_val + 
                  int(flag_sum)*chi_sq_sum_val + int(flag_dif)*chi_sq_dif_val)
        n = (int(flag_u)*n_u + int(flag_d)*n_d + int(flag_sum)*n_sum + 
             int(flag_dif)*n_sum)
        
        d_map["out"] = (chi_sq_val, n)
        return chi_sq_val, n

    def plot_map(self):
        b_variable = self.is_variable()
        d_map = {"flag": b_variable, "out":None}
        d_background = {"flag": False, "out":None}

        b_variable_profile = self.is_variable_profile()
        d_profile = {"flag": b_variable_profile, "out":None}
        d_profile.update({"background":d_background})
        d_map.update({"profile":d_profile})
        
        b_variable_profile_s = self.is_variable_profile_s()
        for calculated_data, b_1 in zip(self._list_calculated_data, 
                                        b_variable_profile_s):
            d_for_iint = calculated_data.plot_map()
            name = calculated_data.get_val("name")
            
            d_hkl = {"flag": False, "out":None}
            #for setup.calc_profile
            d_profile = {"flag": b_1, "out":None}
            d_map.update({("hkl", name): d_hkl, 
                          ("for_iint", name): d_for_iint,
                          ("profile", name): d_profile})
        return d_map
    
    def is_variable_for_iint(self):
        lres = []
        for calculated_data in self._list_calculated_data:
            lres.append(calculated_data.is_variable())
        return lres

    def is_variable_profile_s(self):
        setup = self.get_val("setup")
        res_s = setup.is_variable()
        lres = [res_s]
        return lres

    def is_variable_profile(self):
        lres = self.is_variable_for_iint() 
        lres.extend(self.is_variable_profile_s())
        res = any(lres) 
        return res
    
    def is_variable(self):
        """
        without extinction
        """
        res = any([self.is_variable_profile()])
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
        tth = observed_data.get_val("tth")
        phi = observed_data.get_val("phi")
        int_u_exp = observed_data.get_val("int_u")
        sint_u_exp = observed_data.get_val("sint_u")
        int_d_exp = observed_data.get_val("int_d")
        sint_d_exp = observed_data.get_val("sint_d")
        
        
        int_u_mod, int_d_mod = self.calc_profile(tth, phi, l_crystal) 

        s_int_u_exp = save_2d(tth, phi, int_u_exp)
        s_sint_u_exp = save_2d(tth, phi, sint_u_exp)
        s_int_d_exp = save_2d(tth, phi, int_d_exp)
        s_sint_d_exp = save_2d(tth, phi, sint_d_exp)
        s_int_u_mod = save_2d(tth, phi, int_u_mod)
        s_int_d_mod = save_2d(tth, phi, int_d_mod)
        
        #list of hkl should be added
        s_out = (s_int_u_exp+3*"\n"+s_sint_u_exp+3*"\n"+s_int_u_mod+3*"\n"+
                 s_int_d_exp+3*"\n"+s_sint_d_exp+3*"\n"+s_int_d_mod)
        
        f_out = self._p_f_out
        if f_out is None:
            print("File to save model data is not defined.")
            return
        fid = open(f_out, "w")
        fid.write(s_out)
        fid.close()
        
    
if (__name__ == "__main__"):
  pass
