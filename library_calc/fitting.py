"""
define classe to describe fitting of experiments
"""
__author__ = 'ikibalin'
__version__ = "2019_04_06"
import os
import numpy
import scipy.optimize

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
        
        ls_out += "\n\n\nExperiments\n\n"+"\n\n".join(ls_exp)
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
    
    def calc_chi_sq(self):
        """
        calculate the integral intensity for h, k, l reflections
        """
        chi_sq_res, n_res = 0., 0.
        for experiment in self._list_experiment:
            chi_sq, n = experiment.calc_chi_sq()
            chi_sq_res += chi_sq
            n_res += n
        return chi_sq, n_res 
    
    def refinement(self):
        """
        optimization
        """
        l_variable = [hh for hh in self._list_variable if hh[1]]

        l_param_0 = [variable["val"] for variable in l_variable]

        def tempfunc(l_param):
            for variable, param in zip(l_variable, l_param):
                variable["val"] = param
            chi_sq, n = self.calc_chi_sq()
            return chi_sq

        #if self.ref[0].refin:
        res = scipy.optimize.minimize(tempfunc, l_param_0, method='Nelder-Mead')
        print("refinement is finished")
        print("\n\nResult is:\n", res)
        return res




    
if (__name__ == "__main__"):
  pass
