
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

class ExperimentSingle(dict):
    """
    Class to describe all information concerning to the experiment for monocrystal
    """
    def __init__(self, name=None, setup = SetupSingle(), 
                 list_calculated_data = [], 
                 observed_data = ObservedDataSingle()):
        super(ExperimentSingle, self).__init__()
        self._p_name = None
        self._p_setup = None
        self._list_calculated_data = []
        self._p_observed_data = None
        
        self._refresh(name, setup, observed_data)

    def __repr__(self):
        ls_out = """ExperimentSingle:\n name: {:}\n{:}\n{:}""".format(
                self._p_name, self._p_setup, self._p_observed_data)
        
        ls_calculated_data = []
        for calculated_data in self._list_calculated_data:
            ls_calculated_data.append("{:}".format(calculated_data))

        ls_out += "\n\n\nCalculatedData:\n\n"+"\n\n".join(ls_calculated_data)
        return ls_out

    def _refresh(self, name, setup, observed_data):
        if name is not None:
            self._p_name = name
        if setup is not None:
            self._p_setup = setup
        if observed_data is not None:
            self._p_observed_data = observed_data

            
    def set_val(self, name=None, setup=None, observed_data=None):
        self._refresh(name, setup, observed_data)
        
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
        """
        print(lsout)
    def add_calculated_data(self, observed_data):
        self._list_calculated_data.append(observed_data)

    def del_calculated_data(self, ind):
        self._list_calculated_data.pop(ind)        

    def replace_calculated_data(self, ind, observed_data):
        self._list_calculated_data.pop(ind)
        self._list_calculated_data.insert(ind, observed_data)
    
    def calc_iint_u_d_flip_ratio(self, h, k, l):
        """
        calculate intensity for the given diffraction angle
        """
        setup = self._p_setup
        wave_length = setup.get_val("wave_length")
        beam_polarization = setup.get_val("beam_polarization")
        
        #calculations only for first crystal phase
        for calculated_data in self._list_calculated_data[0:1]:
            iint_u, iint_d, flip_ration = calculated_data.calc_iint_u_d_flip_ratio(
                                      h, k, l, beam_polarization, wave_length)
            
        return iint_u, iint_d, flip_ration 
    
    def calc_chi_sq(self, d_map={}):
        """
        calculate chi square
        """
        if d_map == {}:
            d_map.update(self.plot_map())
        if not(d_map["flag"]|(d_map["out"] is None)):
            chi_sq_val, n = d_map["out"]
            return chi_sq_val, n        
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
                                               h, k, l)
        

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

class ExperimentPowder1D(dict):
    """
    Class to describe all information concerning to the experiment for powder 1d
    """
    def __init__(self, name=None, setup=SetupPowder1D(), 
                 list_calculated_data=[], 
                 observed_data=ObservedDataPowder1D(), mode_chi_sq="up, down"):
        super(ExperimentPowder1D, self).__init__()
        self._p_name = None
        self._p_setup = None
        self._list_calculated_data = []
        self._p_observed_data = None
        self._p_mode_chi_sq = None
        
        self._refresh(name, setup, observed_data, mode_chi_sq)

    def __repr__(self):
        ls_out = """ExperimentPowder1D:\nname: {:}\n{:}\n{:}
 mode_chi_sq: {:}""".format(self._p_name, self._p_setup, self._p_observed_data, 
                            self._p_mode_chi_sq)
        
        ls_calculated_data = []
        for calculated_data in self._list_calculated_data:
            ls_calculated_data.append("{:}".format(calculated_data))

        ls_out += "\n\n\nCalculatedData:\n\n"+"\n\n".join(ls_calculated_data)
        return ls_out

    def _refresh(self, name, setup, observed_data, mode_chi_sq):
        if name is not None:
            self._p_name = name
        if setup is not None:
            self._p_setup = setup
        if observed_data is not None:
            self._p_observed_data = observed_data
        if mode_chi_sq is not None:
            self._p_mode_chi_sq = mode_chi_sq

            
    def set_val(self, name=None, setup=None, observed_data=None, 
                mode_chi_sq=None):
        self._refresh(name, setup, observed_data, mode_chi_sq)
        
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
            scale = 1.*calculated_data.get_val("scale")
            
            i_g = 1.*calculated_data.get_val("crystal").get_val("i_g")
            cell = calculated_data.get_val("crystal").get_val("cell")
            space_groupe = calculated_data.get_val("crystal").get_val("space_groupe")
            h, k, l, mult = setup.calc_hkl(cell, space_groupe, sthovl_min, sthovl_max)
            
            np_iint_u, np_iint_d = calculated_data.calc_iint(h, k, l, beam_polarization)
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
    
    def calc_chi_sq(self, d_map={}):
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

        int_u_mod, int_d_mod = self.calc_profile(tth)
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
        lres = []
        for calculated_data in self._list_calculated_data:
            crystal = calculated_data.get_val("crystal")
            cell = crystal.get_val("cell")
            res = (cell.is_variable()|res_s)
            lres.append(res)
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
   

class ExperimentPowder2D(dict):
    """
    Class to describe all information concerning to the experiment for powder 2d
    """
    def __init__(self, name=None, setup=SetupPowder2D(), 
                 list_calculated_data=[], 
                 observed_data=ObservedDataPowder2D(), mode_chi_sq="up, down"):
        super(ExperimentPowder2D, self).__init__()
        self._p_name = None
        self._p_setup = None
        self._list_calculated_data = []
        self._p_observed_data = None
        self._p_mode_chi_sq = None
        
        self._refresh(name, setup, observed_data, mode_chi_sq)

    def __repr__(self):
        ls_out = """ExperimentPowder2D:\n name: {:}\n{:}\n{:}
 mode_chi_sq: {:}""".format(self._p_name, self._p_setup, self._p_observed_data, 
                            self._p_mode_chi_sq)
        
        ls_calculated_data = []
        for calculated_data in self._list_calculated_data:
            ls_calculated_data.append("{:}".format(calculated_data))

        ls_out += "\n\n\nCalculatedData:\n\n"+"\n\n".join(ls_calculated_data)
            
        return ls_out

    def _refresh(self, name, setup, observed_data, mode_chi_sq):
        if name is not None:
            self._p_name = name
        if setup is not None:
            self._p_setup = setup
        if observed_data is not None:
            self._p_observed_data = observed_data
        if mode_chi_sq is not None:
            self._p_mode_chi_sq = mode_chi_sq

            
    def set_val(self, name=None, setup=None, observed_data=None, 
                mode_chi_sq=None):
        self._refresh(name, setup, observed_data, mode_chi_sq)
        
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
        """
        print(lsout)
    def add_calculated_data(self, observed_data):
        self._list_calculated_data.append(observed_data)

    def del_calculated_data(self, ind):
        self._list_calculated_data.pop(ind)        

    def replace_calculated_data(self, ind, observed_data):
        self._list_calculated_data.pop(ind)
        self._list_calculated_data.insert(ind, observed_data)
    
    def calc_profile(self, tth, phi, d_map={}):
        """
        calculate intensity for the given diffraction angle
        
        tth and phi is 1D data
        """
        if d_map == {}:
            d_map.update(self.plot_map())
        d_profile = d_map["profile"]
        if not(d_profile["flag"]|(d_profile["out"] is None)):
            res_u_2d_int_bkgd, res_d_2d_int_bkgd = d_profile["out"]
            return res_u_2d_int_bkgd, res_d_2d_int_bkgd
        
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
            name = calculated_data.get_val("name")
            scale = 1.*calculated_data.get_val("scale")
            i_g = 1.*calculated_data.get_val("crystal").get_val("i_g")
            cell = calculated_data.get_val("crystal").get_val("cell")
            space_groupe = calculated_data.get_val("crystal").get_val("space_groupe")
            
            
            d_hkl = d_map[("hkl", name)]
            if not(d_hkl["flag"]|(d_hkl["out"] is None)):
                h, k, l, mult = d_hkl["out"]
            else:
                h, k, l, mult = setup.calc_hkl(cell, space_groupe, sthovl_min, sthovl_max)
                d_hkl["out"] = (h, k, l, mult)
            mult_3d = numpy.meshgrid(tth, phi, mult, indexing="ij")[2]
            
            
            d_for_iint = d_map[("for_iint", name)]
            if not(d_for_iint["flag"]|(d_for_iint["out"] is None)):
                f_nucl_sq, f_m_p_sin_sq, f_m_p_cos_sq, cross_sin = d_for_iint["out"]
            else:
                #dimension: (tth_hkl)
                f_nucl_sq, f_m_p_sin_sq, f_m_p_cos_sq, cross_sin = calculated_data.calc_for_iint(h, k, l, d_for_iint)
            
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
            

            d_profile_s = d_map[("profile", name)]
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
    
    def calc_chi_sq(self, d_map={}):
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

        
            
            
        int_u_mod, int_d_mod = self.calc_profile(tth, phi, d_map)
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
        lres = []
        for calculated_data in self._list_calculated_data:
            crystal = calculated_data.get_val("crystal")
            cell = crystal.get_val("cell")
            res = (cell.is_variable()|res_s)
            lres.append(res)
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

    
if (__name__ == "__main__"):
  pass
