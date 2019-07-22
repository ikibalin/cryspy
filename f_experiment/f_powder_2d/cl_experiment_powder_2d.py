
"""
define classe to describe experiment
"""
__author__ = 'ikibalin'
__version__ = "2019_04_06"
import os
import numpy

from f_experiment.f_powder_2d.cl_observed_data_powder_2d import *
from f_experiment.f_powder_2d.cl_calculated_data_powder_2d import *
from f_experiment.f_powder_2d.cl_setup_powder_2d import *

from f_common.cl_variable import *


def save_2d(np_tth, np_phi, int_2d):
    line = "       {:7}    ".format(len(np_phi)) + " ".join(["{:10.2f}".format(hh) for hh in np_tth])
    lsout_u = [line]
    for int_1d, phi in zip(int_2d.transpose(), np_phi):
        line_u = "    {:10.5f}    ".format(phi) + " ".join(["{:10.2f}".format(hh) for hh in int_1d])
        lsout_u.append(line_u)

    return "\n".join(lsout_u)


class ExperimentPowder2D(dict):
    """
    Class to describe all information concerning to the experiment for powder 2d
    """
    def __init__(self, name=None, setup=SetupPowder2D(), 
                 list_calculated_data=[], 
                 observed_data=ObservedDataPowder2D(), flag_chi2_up=None, flag_chi2_down=None, 
                 flag_chi2_sum=None, flag_chi2_diff=None, 
                 file_out=None, file_dir=None):
        super(ExperimentPowder2D, self).__init__()
        self._p_name = None
        self._p_setup = None
        self._list_calculated_data = []
        self._p_observed_data = None
        self._p_flag_chi2_up = None
        self._p_flag_chi2_down = None
        self._p_flag_chi2_sum = None
        self._p_flag_chi2_diff = None
        self._p_file_out = None
        self._p_file_dir = None
        
        self._refresh(name, setup, observed_data, flag_chi2_up, flag_chi2_down, 
        flag_chi2_sum, flag_chi2_diff, file_out, file_dir)

    def __repr__(self):
        ls_out = """ExperimentPowder2D:\n name: {:}
 file_out: {:}\n{:}\n{:}""".format(self._p_name, 
     self._p_file_out, self._p_setup, self._p_observed_data)
        
        ls_calculated_data = []
        for calculated_data in self._list_calculated_data:
            ls_calculated_data.append("{:}".format(calculated_data))

        ls_out += "\n\n\nCalculatedData:\n\n"+"\n\n".join(ls_calculated_data)
            
        return ls_out

    def _refresh(self, name, setup, observed_data, flag_chi2_up, flag_chi2_down, 
                flag_chi2_sum, flag_chi2_diff, file_out, file_dir):
        if name is not None:
            self._p_name = name
        if setup is not None:
            self._p_setup = setup
        if observed_data is not None:
            self._p_observed_data = observed_data
        if flag_chi2_up is not None:
            self._p_flag_chi2_up = flag_chi2_up
        if flag_chi2_down is not None:
            self._p_flag_chi2_down = flag_chi2_down
        if flag_chi2_sum is not None:
            self._p_flag_chi2_sum = flag_chi2_sum
        if flag_chi2_diff is not None:
            self._p_flag_chi2_diff = flag_chi2_diff
        if file_out is not None:
            self._p_file_out = file_out
        if file_dir is not None:
            self._p_file_dir = file_dir

            
    def set_val(self, name=None, setup=None, observed_data=None, 
                flag_chi2_up=None, flag_chi2_down=None, flag_chi2_sum=None, 
                flag_chi2_diff=None, file_out=None, file_dir=None):
        self._refresh(name, setup, observed_data, flag_chi2_up, flag_chi2_down, 
                flag_chi2_sum, flag_chi2_diff, file_out, file_dir)
        
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
flag_chi2_up, flag_chi2_down are flags for refinement "up", "down"
flag_chi2_sum, flag_chi2_diff are flags for refinement "sum", "diff"
file_out is the file name to save model profiles 
file_dir is the working directory
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

        p_u = 1.*beam_polarization.get_val("p_u")
        p_d = (2.*beam_polarization.get_val("flipper_efficiency")-1)*p_u
        
        tth_min = tth.min()
        tth_max = tth.max()+3.
        if tth_max > 180.:
            tth_max = 180.
        sthovl_min = numpy.sin(0.5*tth_min*numpy.pi/180.)/wave_length
        sthovl_max = numpy.sin(0.5*tth_max*numpy.pi/180.)/wave_length

        res_u_2d = numpy.zeros((tth.shape[0],phi.shape[0]), dtype=float)
        res_d_2d = numpy.zeros((tth.shape[0],phi.shape[0]), dtype=float)
        
        l_d_info = []
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
            
            
            #d_hkl = d_map[("hkl", observed_data_name)]
            #if not(d_hkl["flag"]|(d_hkl["out"] is None)):
            #    h, k, l, mult = d_hkl["out"]
            #else:
            h, k, l, mult = cell.calc_hkl(space_groupe, sthovl_min, sthovl_max)
            #d_hkl["out"] = (h, k, l, mult)
            mult_3d = numpy.meshgrid(tth, phi, mult, indexing="ij")[2]
            
            
            #d_for_iint = d_map[("for_iint", observed_data_name)]
            #if not(d_for_iint["flag"]|(d_for_iint["out"] is None)):
            #    f_nucl_sq, f_m_p_sin_sq, f_m_p_cos_sq, cross_sin = d_for_iint["out"]
            #else:
            #    #dimension: (tth_hkl)
            f_nucl_sq, f_m_p_sin_sq, f_m_p_cos_sq, cross_sin, d_info_cd = calculated_data.calc_for_iint(
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
            

            #if not(d_profile_s["flag"]|(d_profile_s["out"] is None)):
            #    profile_3d = d_profile_s["out"]
            #else:
            #    #dimension: (tth, phi, tth_hkl)
            profile_3d = setup.calc_profile(tth, phi, tth_hkl, i_g)
            #d_profile_s["out"] = profile_3d
            
            
            res_u_3d = profile_3d*i_u_3d
            res_d_3d = profile_3d*i_d_3d

            res_u_2d += scale*res_u_3d.sum(axis=2) 
            res_d_2d += scale*res_d_3d.sum(axis=2) 
            
            
            d_info = {}
            d_info.update(d_info_cd)
            l_d_info.append(d_info)
        
        d_info_out = {"crystal": l_d_info}            
        #d_profile["out"] = res_u_2d+int_bkgd, res_d_2d+int_bkgd
        return res_u_2d+int_bkgd, res_d_2d+int_bkgd, d_info_out
    
    def calc_chi_sq(self, l_crystal, d_map={}):
        """
        calculate chi square
        """
        #if d_map == {}:
        #    d_map.update(self.plot_map())
        #if not(d_map["flag"]|(d_map["out"] is None)):
        #    chi_sq_val, n = d_map["out"]
        #    return chi_sq_val, n
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

        int_u_mod, int_d_mod, d_info_prof = self.calc_profile(tth, phi, l_crystal)
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
        
        flag_u = self._p_flag_chi2_up
        flag_d = self._p_flag_chi2_down
        flag_sum = self._p_flag_chi2_sum
        flag_dif = self._p_flag_chi2_diff
        
        chi_sq_val = (int(flag_u)*chi_sq_u_val + int(flag_d)*chi_sq_d_val + 
                  int(flag_sum)*chi_sq_sum_val + int(flag_dif)*chi_sq_dif_val)
        n = (int(flag_u)*n_u + int(flag_d)*n_d + int(flag_sum)*n_sum + 
             int(flag_dif)*n_sum)
        
        d_info_out = {"chi_sq_val": chi_sq_val, "n": n}
        d_info_out.update(d_info_prof)
        return chi_sq_val, n, d_info_out        

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
        
        d_info_in = {}
        int_u_mod, int_d_mod, d_info_out = self.calc_profile(tth, phi, l_crystal, d_info_in) 

        s_int_sum_exp = save_2d(tth, phi, int_u_exp+int_d_exp)
        s_int_sum_mod = save_2d(tth, phi, int_u_mod+int_d_mod)
        s_int_dif_exp = save_2d(tth, phi, int_u_exp-int_d_exp)
        s_int_dif_mod = save_2d(tth, phi, int_u_mod-int_d_mod)
        
        s_int_sum_dif = save_2d(tth, phi, int_u_exp+int_d_exp-int_u_mod-int_d_mod)
        s_int_dif_dif = save_2d(tth, phi, int_u_exp-int_d_exp-int_u_mod+int_d_mod)

        
        #list of hkl should be added
        s_out = (s_int_sum_exp+3*"\n"+s_int_sum_mod+3*"\n"+s_int_sum_dif+3*"\n"+
                 s_int_dif_exp+3*"\n"+s_int_dif_mod+3*"\n"+s_int_dif_dif)
        
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
        
        chi_sq_val, n, d_info_full = self.calc_chi_sq(l_crystal)
        
        ls_out = []
        l_d_info = d_info_full["crystal"]
        for d_info in l_d_info:
            ls_out.append("\n\n   h   k   l   fn_real   fn_imag  f_nucl_sq f_m_p_sin_sq f_m_p_cos_sq cross_sin")
            h, k, l, f_nucl = d_info["h"], d_info["k"], d_info["l"], d_info["f_nucl"]
            f_nucl_sq, f_m_p_sin_sq = d_info["f_nucl_sq"], d_info["f_m_p_sin_sq"]
            f_m_p_cos_sq, cross_sin = d_info["f_m_p_cos_sq"], d_info["cross_sin"]
            
            sft_11, sft_12, sft_13 = d_info["sft_11"], d_info["sft_12"], d_info["sft_13"]
            sft_21, sft_22, sft_23 = d_info["sft_21"], d_info["sft_22"], d_info["sft_23"]
            sft_31, sft_32, sft_33 = d_info["sft_31"], d_info["sft_32"], d_info["sft_33"]

            for h1, k1, l1, f, i_1, i_2, i_3, i_4 in zip(h, k, l, f_nucl, f_nucl_sq, f_m_p_sin_sq, f_m_p_cos_sq, cross_sin):
                ls_out.append(""" {:3} {:3} {:3}  {:8.3f}  {:8.3f}  {:8.1f}  {:8.1f}  {:8.1f}  {:8.1f}""".format(
                        h1, k1, l1, f.real, f.imag, i_1, i_2, i_3, i_4))
            ls_out.append("\n\nStructure factor tensor in a*, b*, c* (real part, mu_B) ")
            ls_out.append("\n   h   k   l   sft_11   sft_22   sft_33    sft_12   sft_13   sft_23")
            for h1, k1, l1, c11, c22,  c33, c12, c13, c23 in zip(
                    h, k, l, sft_11, sft_22, sft_33, sft_12, sft_13, sft_23):
                ls_out.append(""" {:3} {:3} {:3} {:8.3f} {:8.3f} {:8.3f}  {:8.3f} {:8.3f} {:8.3f}""".format(
                        h1, k1, l1, c11.real, c22.real, c33.real, c12.real, c13.real, c23.real))
        
        return s_out+"\n".join(ls_out)

    
if (__name__ == "__main__"):
  pass
