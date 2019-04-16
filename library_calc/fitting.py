"""
define classe to describe fitting of experiments
"""
__author__ = 'ikibalin'
__version__ = "2019_04_06"
import os
import numpy

from experiment_pd import *


class Fitting(dict):
    """
    Class to describe fitting
    """
    def __init__(self):
        super(Fitting, self).__init__()
        self._list_experiment = []
        self._refresh()

    def __repr__(self):
        lsout = """Fitting\n """.format()
        return lsout

    def _refresh(self):
        pass
        #if not(isinstance(setup, type(None))):
        #    self._p_setup = setup

            
    def set_val(self):
        self._refresh()
        
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

    
        """
        print(lsout)


    def add_experiment(self, experiment):
        self._list_experiment.append(experiment)

    def del_experiment(self, ind):
        self._list_experiment.pop(ind)        

    def replace_experiment(self, ind, experiment):
        self._list_experiment.pop(ind)
        self._list_experiment.insert(ind, experiment)
    
    def calc_chi_sq(self):
        """
        calculate the integral intensity for h, k, l reflections
        """
        chi_sq = 0.
        for experiment in self._list_experiment:
            chi_sq += experiment.calc_chi_sq()
            
        return chi_sq 
    




    
if (__name__ == "__main__"):
  pass
