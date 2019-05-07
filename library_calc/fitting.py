"""
define classe to describe fitting of experiments
"""
__author__ = 'ikibalin'
__version__ = "2019_04_06"
import os
import numpy
import scipy.optimize
import time

from experiment import *
from variable import *

class Fitting(dict):
    """
    Class to describe fitting
    """
    def __init__(self, l_experiment=None, l_variable=None):
        super(Fitting, self).__init__()
        self._list_experiment = []
        self._list_variable = []
        self._refresh(l_experiment, l_variable)

    def __repr__(self):
        ls_out = """Fitting\n """.format()
        ls_exp = []
        for epxeriment in self._list_experiment:
            ls_exp.append("{:}".format(epxeriment))
            
        ls_var = []
        for variable in self._list_variable:
            ls_var.append("{:}".format(variable))
        
        ls_out += "\n\n\nExperiment:\n\n"+"\n\n".join(ls_exp)
        ls_out += "\n\n\nVariable:\n\n"+"\n\n".join(ls_var)
        return ls_out

    def _refresh(self, l_experiment=None, l_variable=None):
        if l_experiment is not None:
            self._list_experiment = l_experiment
        if l_variable is not None:
            self._list_variable = l_variable
            
    def set_val(self, l_experiment=None, l_variable=None):
        self._refresh(l_experiment, l_variable)
        
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
        if not(d_map["flag"]|(d_map["out"] is None)):
            chi_sq, n_res = d_map["out"]
            return chi_sq, n_res 
            
        chi_sq_res, n_res = 0., 0.
        for experiment in self._list_experiment:
            name = experiment.get_val("name")
            chi_sq, n = experiment.calc_chi_sq(d_map[("chi_sq", name)])
            chi_sq_res += chi_sq
            n_res += n

        d_map["out"] = (chi_sq, n_res)
        return chi_sq, n_res 
    
    def refinement(self, d_map={}):
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
        chi_sq, n = self.calc_chi_sq(d_map)
        
        l_param_0 = [variable["val"] for variable in l_variable]

        def tempfunc(l_param):
            for variable, param in zip(l_variable, l_param):
                variable["val"] = param
            chi_sq, n = self.calc_chi_sq(d_map)
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
        return d_map

    
    def is_variable(self):
        lres = []
        for experiment in self._list_experiment:
            lres.append(experiment.is_variable())
        res = any(lres)
        return res

    def get_variables(self):
        l_variable = []
        for experiment in self._list_experiment:
            l_var = experiment.get_variables()
            l_variable.extend(l_var)
        return l_variable
    
if (__name__ == "__main__"):
  pass
