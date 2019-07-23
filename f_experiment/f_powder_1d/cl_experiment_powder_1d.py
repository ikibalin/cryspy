
"""
define classe to describe experiment
"""
__author__ = 'ikibalin'
__version__ = "2019_04_06"
import os
import numpy

from f_experiment.f_powder_1d.cl_observed_data_powder_1d import *
from f_experiment.f_powder_1d.cl_calculated_data_powder_1d import *
from f_experiment.f_powder_1d.cl_setup_powder_1d import *

from f_common.cl_variable import *

class ExperimentPowder1D(dict):
    """
    Class to describe all information concerning to the experiment for powder 1d
    """
    def __init__(self, name=None, setup=SetupPowder1D(), 
                 list_calculated_data=[], 
                 observed_data=ObservedDataPowder1D(), flag_chi2_up=False, 
                 flag_chi2_down=False, flag_chi2_sum=False, flag_chi2_diff=False, 
                 file_out=None, file_dir=None, excl_tth_min=[], excl_tth_max=[]):
        super(ExperimentPowder1D, self).__init__()
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
        self._p_excl_tth_min = None
        self._p_excl_tth_max = None

        self._refresh(name, setup, observed_data, flag_chi2_up, flag_chi2_down, 
        flag_chi2_sum, flag_chi2_diff, file_out, file_dir, excl_tth_min, excl_tth_max)

    def __repr__(self):
        ls_out = """ExperimentPowder1D:\nname: {:}
 file_out: {:}\n{:}\n{:}""".format(self._p_name, 
  self._p_file_out, self._p_setup, self._p_observed_data)
        
        ls_calculated_data = []
        for calculated_data in self._list_calculated_data:
            ls_calculated_data.append("{:}".format(calculated_data))

        ls_out += "\n\n\nCalculatedData:\n\n"+"\n\n".join(ls_calculated_data)
        return ls_out

    def _refresh(self, name, setup, observed_data, flag_chi2_up, flag_chi2_down, 
                flag_chi2_sum, flag_chi2_diff, file_out, file_dir, excl_tth_min, excl_tth_max):
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
        if excl_tth_min is not None:
            self._p_excl_tth_min = excl_tth_min
        if excl_tth_max is not None:
            self._p_excl_tth_max = excl_tth_max
            
    def set_val(self, name=None, setup=None, observed_data=None, 
                flag_chi2_up=None, flag_chi2_down=None, 
                flag_chi2_sum=None, flag_chi2_diff=None, file_out=None, file_dir=None, 
                excl_tth_min=None, excl_tth_max=None):
        self._refresh(name, setup, observed_data, flag_chi2_up, flag_chi2_down, 
                flag_chi2_sum, flag_chi2_diff, file_out, file_dir, excl_tth_min, excl_tth_max)
        
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
flag_chi2_up, flag_chi2_down are flags for  refinement: "up down"
flag_chi2_sum, flag_chi2_diff  are flags for  refinement: "sum diff"
file_out is the file name of experimental and model data (basename)
file_dir is the working directory
excl_tth_min is the list of excluded ttheta from up, down and sum (difference is taken into account), minimal
excl_tth_max is the list of excluded ttheta from up, down and sum (difference is taken into account), maximal
        """
        print(lsout)
    def add_calculated_data(self, observed_data):
        self._list_calculated_data.append(observed_data)

    def del_calculated_data(self, ind):
        self._list_calculated_data.pop(ind)        

    def replace_calculated_data(self, ind, observed_data):
        self._list_calculated_data.pop(ind)
        self._list_calculated_data.insert(ind, observed_data)
    
    def calc_profile(self, tth, l_crystal, d_in={}):
        """
        calculate intensity for the given diffraction angle
        """
        setup = self._p_setup
        background = setup.get_val("background")
        int_bkgd = background.interpolate_by_points(tth)
        wave_length = setup.get_val("wave_length")
        beam_polarization = setup.get_val("beam_polarization")
        
        tth_min = tth.min()
        tth_max = tth.max()+3. 
        if tth_max > 180.:
            tth_max = 180.
        sthovl_min = numpy.sin(0.5*tth_min*numpy.pi/180.)/wave_length
        sthovl_max = numpy.sin(0.5*tth_max*numpy.pi/180.)/wave_length

        res_u_1d = numpy.zeros(tth.shape[0], dtype=float)
        res_d_1d = numpy.zeros(tth.shape[0], dtype=float)
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
            h, k, l, mult = cell.calc_hkl(space_groupe, sthovl_min, sthovl_max)
            
            np_iint_u, np_iint_d, d_info_cd = calculated_data.calc_iint(h, k, l, 
                                                    beam_polarization, crystal)
            sthovl_hkl = cell.calc_sthovl(h, k, l)
            
            tth_hkl_rad = 2.*numpy.arcsin(sthovl_hkl*wave_length)
            tth_hkl = tth_hkl_rad*180./numpy.pi
            d_info_cd.update({"tth_hkl": tth_hkl})
            profile_2d = setup.calc_profile(tth, tth_hkl, i_g)
            
            
            np_iint_u_2d = numpy.meshgrid(tth, np_iint_u*mult, indexing="ij")[1]

            np_iint_d_2d = numpy.meshgrid(tth, np_iint_d*mult, indexing="ij")[1]
            
            res_u_2d = profile_2d*np_iint_u_2d 
            res_d_2d = profile_2d*np_iint_d_2d 

            res_u_1d += scale*res_u_2d.sum(axis=1) 
            res_d_1d += scale*res_d_2d.sum(axis=1) 
            
            d_info = {}
            d_info.update(d_info_cd)
            l_d_info.append(d_info)
        
        d_info_out = {"crystal": l_d_info}
            
        return res_u_1d+int_bkgd, res_d_1d+int_bkgd, d_info_out
    
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

        int_u_mod, int_d_mod, d_info_prof = self.calc_profile(tth, l_crystal)
        sint_sum_exp = (sint_u_exp**2 + sint_d_exp**2)**0.5

        chi_sq_u = ((int_u_mod-int_u_exp)/sint_u_exp)**2
        chi_sq_d = ((int_d_mod-int_d_exp)/sint_d_exp)**2

        chi_sq_sum = ((int_u_mod+int_d_mod-int_u_exp-int_d_exp)/sint_sum_exp)**2
        chi_sq_dif = ((int_u_mod-int_d_mod-int_u_exp+int_d_exp)/sint_sum_exp)**2

        cond_u = numpy.logical_not(numpy.isnan(chi_sq_u))
        cond_d = numpy.logical_not(numpy.isnan(chi_sq_d))
        cond_sum = numpy.logical_not(numpy.isnan(chi_sq_sum))
        cond_dif = numpy.logical_not(numpy.isnan(chi_sq_dif))
        
        #exclude region
        l_excl_tth_min = self._p_excl_tth_min
        l_excl_tth_max = self._p_excl_tth_max
        for excl_tth_min, excl_tth_max in zip(l_excl_tth_min, l_excl_tth_max):
            cond_1 = numpy.logical_or(tth < 1.*excl_tth_min, tth > 1.*excl_tth_max)
            cond_u = numpy.logical_and(cond_u, cond_1)
            cond_d = numpy.logical_and(cond_d, cond_1)
            cond_sum = numpy.logical_and(cond_sum, cond_1)

        chi_sq_u_val = (chi_sq_u[cond_u]).sum()
        n_u = cond_u.sum()
        
        chi_sq_d_val = (chi_sq_d[cond_d]).sum()
        n_d = cond_d.sum()

        chi_sq_sum_val = (chi_sq_sum[cond_sum]).sum()
        n_sum = cond_sum.sum()

        chi_sq_dif_val = (chi_sq_dif[cond_dif]).sum()
        n_sum = cond_dif.sum()

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
        setup = self.get_val("setup")
        zero_shift = setup.get_val("zero_shift")
        int_u_mod, int_d_mod, d_info = self.calc_profile(tth, l_crystal) 
        
        s_1 = "    ttheta       exp_up     s_exp_up       mod_up     exp_down   s_exp_down     mod_down\n"
        l_s_2 = ["{:10.2f} {:12.5f} {:12.5f} {:12.5f} {:12.5f} {:12.5f} {:12.5f}".format(
                hh_1, hh_2, hh_3, hh_4, hh_5, hh_6, hh_7) for 
                hh_1, hh_2, hh_3, hh_4, hh_5, hh_6, hh_7 in 
                zip(tth, int_u_exp, sint_u_exp, int_u_mod, 
                         int_d_exp, sint_d_exp, int_d_mod)]

        s_int = s_1 + "\n".join(l_s_2)
        s_out = s_int
        #hkl should be added
        ls_int = []
        for i_phase, d_info_cd in enumerate(d_info["crystal"]):
            ls_int.append("\n\n\n\n#phase {:}".format(i_phase+1))
            tth_hkl = d_info_cd["tth_hkl"] + 1.*zero_shift 

            np_h = d_info_cd["h"]
            np_k = d_info_cd["k"]
            np_l = d_info_cd["l"]
            for h, k, l, tth in zip(np_h, np_k, np_l, tth_hkl):
                ls_int.append(" {:} {:} {:} {:}".format(h, k, l, tth))
        s_out += "\n".join(ls_int)
        
        
        
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
            ls_out.append("\n\n   h   k   l   fn_real   fn_imag   iint_up iint_down")
            h, k, l, f_nucl = d_info["h"], d_info["k"], d_info["l"], d_info["f_nucl"]
            iint_u, iint_d = d_info["iint_u"], d_info["iint_d"]
            sft_11, sft_12, sft_13 = d_info["sft_11"], d_info["sft_12"], d_info["sft_13"]
            sft_21, sft_22, sft_23 = d_info["sft_21"], d_info["sft_22"], d_info["sft_23"]
            sft_31, sft_32, sft_33 = d_info["sft_31"], d_info["sft_32"], d_info["sft_33"]

            for h1, k1, l1, f, i_u, i_d in zip(h, k, l, f_nucl, iint_u, iint_d):
                ls_out.append(""" {:3} {:3} {:3}  {:8.3f}  {:8.3f}  {:8.1f}  {:8.1f}""".format(
                        h1, k1, l1, f.real, f.imag, i_u, i_d))
            ls_out.append("\n\nStructure factor tensor in a*, b*, c* (real part, mu_B) ")
            ls_out.append("\n   h   k   l   sft_11   sft_22   sft_33    sft_12   sft_13   sft_23")
            for h1, k1, l1, c11, c22,  c33, c12, c13, c23 in zip(
                    h, k, l, sft_11, sft_22, sft_33, sft_12, sft_13, sft_23):
                ls_out.append(""" {:3} {:3} {:3} {:8.3f} {:8.3f} {:8.3f}  {:8.3f} {:8.3f} {:8.3f}""".format(
                        h1, k1, l1, c11.real, c22.real, c33.real, c12.real, c13.real, c23.real))
        
        return s_out+"\n".join(ls_out)

    
if (__name__ == "__main__"):
  pass
