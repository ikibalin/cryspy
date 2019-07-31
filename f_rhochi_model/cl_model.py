"""
define classe to describe model of experiments
"""
__author__ = 'ikibalin'
__version__ = "2019_04_06"
import os
import sys
import numpy
import scipy.optimize
import time

import matplotlib.pyplot

from f_experiment.f_powder_1d.cl_experiment_powder_1d import *
from f_experiment.f_powder_2d.cl_experiment_powder_2d import *
from f_experiment.f_single.cl_experiment_single import *
from f_common.cl_variable import *
from f_common.error_simplex import *

class Model(dict):
    """
    Class to describe model
    """
    def __init__(self, l_experiment=None, l_variable=None, l_crystal=None, name="", file_out=None, file_dir=None):
        super(Model, self).__init__()
        self._list_experiment = []
        self._list_crystal = []
        self._list_variable = []
        self._p_name = None
        self._p_file_out = None
        self._p_file_dir = None
        self._refresh(l_experiment, l_crystal, l_variable, name, file_out, file_dir)

    def __repr__(self):
        ls_out = """Model\n name: {:}\n file_out: {:}\n""".format(self._p_name, 
                                  self._p_file_out)
        ls_exp = []
        for experiment in self._list_experiment:
            ls_exp.append("{:}".format(experiment))

        ls_cry = []
        for crystal in self._list_crystal:
            ls_cry.append("{:}".format(crystal))
            
        
        ls_out += "\n\n\nExperiment:\n\n"+"\n\n".join(ls_exp)
        ls_out += "\n\n\nCrystal:\n\n"+"\n\n".join(ls_cry)
        return ls_out

    def _refresh(self, l_experiment=None, l_crystal=None, l_variable=None, 
                 name=None, file_out=None, file_dir=None):
        if l_experiment is not None:
            self._list_experiment = l_experiment
        if l_crystal is not None:
            self._list_crystal = l_crystal
        if l_variable is not None:
            self._list_variable = l_variable
        if name is not None:
            self._p_name = name
        if file_out is not None:
            self._p_file_out = file_out
        if file_dir is not None:
            self._p_file_dir = file_dir
            
    def set_val(self, l_experiment=None, l_crystal=None, l_variable=None, 
                name=None, file_out=None, file_dir=None):
        self._refresh(l_experiment, l_crystal, l_variable, name, file_out, file_dir)
        
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

name is the name of model parameter
file_out is the file name for listing, full path should be given
    
        """
        print(lsout)


    def add_experiment(self, experiment):
        self._list_experiment.append(experiment)

    def del_experiment(self, ind):
        self._list_experiment.pop(ind)        

    def replace_experiment(self, ind, experiment):
        self._list_experiment.pop(ind)
        self._list_experiment.insert(ind, experiment)

    def add_crystal(self, crystal):
        self._list_crystal.append(crystal)

    def del_crystal(self, ind):
        self._list_crystal.pop(ind)        

    def replace_crystal(self, ind, crystal):
        self._list_crystal.pop(ind)
        self._list_crystal.insert(ind, crystal)

    def add_variable(self, variable):
        self._list_variable.append(variable)

    def del_variable(self, ind):
        self._list_variable.pop(ind)        

    def replace_variable(self, ind, variable):
        self._list_variable.pop(ind)
        self._list_variable.insert(ind, variable)
    
    def calc_chi_sq(self, d_info_in={}):
        """
        calculate the integral intensity for h, k, l reflections
        """
        self.apply_constraint()

        #if not(d_map["flag"]|(d_map["out"] is None)):
        #    chi_sq, n_res = d_map["out"]
        #    return chi_sq, n_res 
            
        l_crystal = self._list_crystal

        chi_sq_res, n_res = 0., 0.
        if "l_experiment_info" in d_info_in.keys():
            l_experiment_info_in = d_info_in["l_experiment_info"]
        else:
            l_experiment_info_in = [{} for hh in self._list_experiment] 
        l_experiment_info_out = []

        for experiment, d_info_in_2 in zip(self._list_experiment, l_experiment_info_in):
            name = experiment.get_val("name")
            chi_sq, n, d_info_out_2 = experiment.calc_chi_sq(l_crystal, d_info_in_2)
            #print("experiment: {:} {:} {:}".format(name, chi_sq, n))
            chi_sq_res += chi_sq
            n_res += n
            l_experiment_info_out.append(d_info_out_2)
        d_info_out = {}
        d_info_out["l_experiment_info"] = l_experiment_info_out
        return chi_sq_res, n_res, d_info_out 
    
    def apply_constraint(self):
        for crystal in self._list_crystal:
            crystal.apply_constraint()
    
    def refine_model(self, d_info_in={}):
        """
        optimization
        """
            
        l_variable = self.get_variables()
        if l_variable == []:
            print("No variables found")
            return None
        
        #to load d_map
        chi_sq, n, d_info_out = self.calc_chi_sq(d_info_in)
        l_param_0 = [variable["val"] for variable in l_variable]

        def tempfunc(l_param):
            for variable, param in zip(l_variable, l_param):
                variable["val"] = param
            chi_sq, n, d_hh = self.calc_chi_sq(d_info_out)
            return chi_sq

        #if self.ref[0].refin:
        print("starting chi_sq/n: {:.2f} (n = {:}).".format(chi_sq*1./n, int(n)))
        print("\nrefinement started for parameters:")
        ls_out = " ".join(["{:12}".format(var[3].rjust(12)) if len(var[3])<=12 
                           else "{:12}".format(var[3][-12:]) for var in l_variable]) + "       chi_sq"
        print(ls_out)
        aa = time.time()
        """
        res, m_error, infodict, errmsg, ier = \
            scipy.optimize.leastsq(tempfunc, l_param_0, full_output=1)

        """
        res = scipy.optimize.minimize(tempfunc, l_param_0, method='Nelder-Mead', 
                                      callback=self._f_callback, options = {"fatol": 0.01*n})
        bb = time.time()
        print("refinement complete, time {:.2f} sec.\n\nfinal chi_sq/n: {:.2f}\nstarted chi_sq/n: {:.2f}".format(bb-aa, res.fun*1./n, chi_sq*1./n))
        
        #it is not checked
        m_error, dist_hh = error_estimation_simplex(res["final_simplex"][0], res["final_simplex"][1], tempfunc)
        
        l_sigma = []
        for i, val_2 in zip(range(m_error.shape[0]), dist_hh):
            #slightly change definition, instead of (n-k) here is n
            error = (abs(m_error[i,i])*1./n)**0.5
            if m_error[i,i] < 0.:
                print(50*"*"+"\nError is incorrect\n(negative diagonal elements of Hessian)\n"+50*"*")
            if val_2 > error:
                print(50*"*"+"\nErrors is incorrect\n(minimum is not found)\n"+50*"*")
                
            l_sigma.append(max(error,val_2))
            
        for variable, sigma in zip(l_variable, l_sigma):
            variable[4] = sigma
        return res
    
    def print_report(self):
        ls_out = []
        ls_out.append("Model: {:}\n".format(self._p_name))
        ls_out.append("Crystals:\n")
        l_crystal = self._list_crystal
        for crystal in l_crystal:
            s_out_1 = crystal.print_report()
            ls_out.append(s_out_1)
        ls_out.append("\n\nExperiments:\n")
        for experiment in self._list_experiment:
            s_out_1 = experiment.print_report(l_crystal)
            ls_out.append(s_out_1)
        return "\n".join(ls_out)
    
    def _f_callback(self, *arg):
        res_x = arg[0]
        ls_out = " ".join(["{:12.5f}".format(hh) for hh in res_x])
        if len(arg) > 1:
            res_fun = arg[1]
            ls_out += " {:12.1f}".format(res_fun.fun)
        print(ls_out)
    
    def plot_map(self):
        d_map = {"flag": self.is_variable(), "out":None}
        for experiment in self._list_experiment:
            d_exp = experiment.plot_map()
            name = experiment.get_val("name")
            d_map.update({("chi_sq",name):d_exp})
        for crystal in self._list_crystal:
            d_cry = crystal.plot_map()
            name = crystal.get_val("name")
            #???
            d_map.update({("crystal",name):d_cry})
        return d_map

    
    def is_variable(self):
        lres = []
        for experiment in self._list_experiment:
            lres.append(experiment.is_variable())
        for crystal in self._list_crystal:
            lres.append(crystal.is_variable())
        res = any(lres)
        return res

    def get_variables(self):
        l_variable = []
        for experiment in self._list_experiment:
            l_var = experiment.get_variables()
            l_variable.extend(l_var)
        for crystal in self._list_crystal:
            l_var = crystal.get_variables()
            l_variable.extend(l_var)
        return l_variable
    
    def save_exp_mod_data(self):
        l_crystal = self._list_crystal
        for experiment in self._list_experiment:
            experiment.save_exp_mod_data(l_crystal)
    


if (__name__ == "__main__"):
    pass
