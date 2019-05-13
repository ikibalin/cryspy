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

from experiment import *
from variable import *

class Model(dict):
    """
    Class to describe model
    """
    def __init__(self, l_experiment=None, l_variable=None, l_crystal=None):
        super(Model, self).__init__()
        self._list_experiment = []
        self._list_crystal = []
        self._list_variable = []
        self._refresh(l_experiment, l_crystal, l_variable)

    def __repr__(self):
        ls_out = """Model\n """.format()
        ls_exp = []
        for epxeriment in self._list_experiment:
            ls_exp.append("{:}".format(epxeriment))

        ls_cry = []
        for crystal in self._list_crystal:
            ls_cry.append("{:}".format(crystal))
            
        
        ls_out += "\n\n\nExperiment:\n\n"+"\n\n".join(ls_exp)
        ls_out += "\n\n\nCrystal:\n\n"+"\n\n".join(ls_cry)
        return ls_out

    def _refresh(self, l_experiment=None, l_crystal=None, l_variable=None):
        if l_experiment is not None:
            self._list_experiment = l_experiment
        if l_crystal is not None:
            self._list_crystal = l_crystal
        if l_variable is not None:
            self._list_variable = l_variable
            
    def set_val(self, l_experiment=None, l_crystal=None, l_variable=None):
        self._refresh(l_experiment, l_crystal, l_variable)
        
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

    def print_vals(self):
        """
        give a list of parameters with small descripition
        """
        lsout = """
Parameters:

None
    
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
        self._list_crystal.insert(ind, experiment)

    def add_variable(self, variable):
        self._list_variable.append(variable)

    def del_variable(self, ind):
        self._list_variable.pop(ind)        

    def replace_variable(self, ind, variable):
        self._list_variable.pop(ind)
        self._list_variable.insert(ind, variable)
    
    def calc_chi_sq(self, d_map={}):
        """
        calculate the integral intensity for h, k, l reflections
        """
        if d_map == {}:
            d_map.update(self.plot_map())
        #if not(d_map["flag"]|(d_map["out"] is None)):
        #    chi_sq, n_res = d_map["out"]
        #    return chi_sq, n_res 
            
        l_crystal = self._list_crystal

        chi_sq_res, n_res = 0., 0.

        for experiment in self._list_experiment:
            name = experiment.get_val("name")
            chi_sq, n = experiment.calc_chi_sq(l_crystal)#d_map[("chi_sq", name)]
            chi_sq_res += chi_sq
            n_res += n

        d_map["out"] = (chi_sq, n_res)
        return chi_sq, n_res 
    
    def refine_model(self, d_map={}):
        """
        optimization
        """
        if d_map == {}:
            d_map.update(self.plot_map())
            
        l_variable = self.get_variables()
        if l_variable == []:
            print("No variables found")
            return None
        
        #to load d_map
        chi_sq, n = self.calc_chi_sq()#d_map
        l_param_0 = [variable["val"] for variable in l_variable]

        def tempfunc(l_param):
            for variable, param in zip(l_variable, l_param):
                variable["val"] = param
            d_map_new = {}
            chi_sq, n = self.calc_chi_sq(d_map_new)#d_map
            return chi_sq

        #if self.ref[0].refin:
        print("starting chi_sq/n: {:.2f}".format(chi_sq*1./n))
        print("\nrefinement started for parameters:")
        ls_out = " ".join(["{:12}".format(var[3].rjust(12)) if len(var[3])<=12 
                           else "{:12}".format(var[3][-12:]) for var in l_variable])
        print(ls_out)
        aa = time.time()
        res = scipy.optimize.minimize(tempfunc, l_param_0, method='Nelder-Mead', callback=self._f_callback)
        bb = time.time()
        print("refinement complete, time {:.2f} sec.\n\nfinal chi_sq/n: {:.2f}".format(bb-aa, res.fun*1./n))
        #print("\n\nResult is:\n", res)
        
        #very rough estimation of the errorbar        
        l_sigma = [float(hh) for hh in numpy.std(res["final_simplex"][0],axis=0)]
        for variable, sigma in zip(l_variable, l_sigma):
            variable[4] = sigma
        return res
    
    def _f_callback(self, res_x):
        ls_out = " ".join(["{:12.5f}".format(hh) for hh in res_x])
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
