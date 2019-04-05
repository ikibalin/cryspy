"""
define resolution for the powder diffractometer along ttheta
"""
__author__ = 'ikibalin'
__version__ = "2019_03_29"
import numpy


class UnitCell(dict):
    """
    Peak Profile
    """
    def __init__(self, a = 1.0, b = 1.0, c = 1.0, alpha = 1.0, beta = 90.0, 
                 gamma= 90., singony = "Monoclinic"):
        super(Crystal, self).__init__()
        dd= {"a": a, "b": b, "c": c, "alpha": alpha, "beta": beta, "gamma": gamma,
             "singony": singony}
        self.update(dd)
        
    def __repr__(self):
        lsout = """Unit cell: \n a: {:}\n b: {:}\n c: {:}\n alpha: {:}
 beta: {:}\n gamma: {:}\n singony: {:}""".format(self["a"], self["b"], 
                 self["c"], self["alpha"], self["beta"], self["gamma"], self["singony"])
        return lsout
    
    def calc_model(self, a = None, b = None, c = None, alpha = None, 
                   beta = None, gamma= None, singony = None):
        if a != None:
            self["a"] = a
        if b != None:
            self["b"] = b
        if c != None:
            self["c"] = c
        if alpha != None:
            self["alpha"] = alpha
        if beta != None:
            self["beta"] = beta
        if gamma != None:
            self["gamma"] = gamma
            
        d_out = dict()
        self.update(d_out)
        
        
    def set_vals(self, d_vals, refresh = False):
        """
        Set values 
        """
        keys = d_vals.keys()
        llab = ["a", "b", "c", "alpha", "beta", "gamma", "singony"]
        llab_in = [(hh in keys) for hh in llab]
        for lab, cond_h in zip(llab, llab_in):
            if cond_h:
                self[lab] = d_vals[lab]
                
            
        if refresh:
            #tth = self["tth"]
            self.calc_model()

            
    def soft_copy(self):
        """
        Soft copy of the object with saving the links on the internal parameter of the object
        """
        obj_new = UnitCell(a=self["a"], b= self["b"], c = self["c"], 
                           alpha = self["alpha"], beta=self["beta"], 
                           gamma=self["gamma"], singony=self["singony"])
        #llab = ["fn", "sft"]
        #keys = self.keys()
        #llab_in = [(hh in keys) for hh in llab]
        #for lab, cond_h in zip(llab, llab_in):
        #    if cond_h:
        #        obj_new[lab] = self[lab]
        return obj_new



class Crystal(dict):
    """
    Peak Profile
    """
    def __init__(self, crystal_symmetry = CrystSymmetry(), 
                 unit_cell = UnitCell(), atom = Atom()):
        super(Crystal, self).__init__()
        dd= {"crystal_symmetry": crystal_symmetry, "unit_cell": unit_cell, 
             "atom": atom}
        self.update(dd)
        
    def __repr__(self):
        lsout = """Phase: \n crystal symmetry: {:}\n unit cell {:}
 atom {:}""".format(self["crystal_symmetry"], self["unit_cell"], self["atom"])
        return lsout
    
    def calc_model(self, reflection, crystal_symmetry = None, unit_cell = None,
                   atom = None):
        if crystal_symmetry != None:
            self["crystal_symmetry"] = crystal_symmetry 
        if unit_cell != None:
            self["unit_cell"] = unit_cell 
        if atom != None:
            self["atom"] = atom 
            
        np_fn_1d = None #complex 1D numpy array of nuclear structure factors over list of hkl
        np_sft_2d = None #complex 2D numpy array of tensor structure factors /11, 22, 33, 12, 13, 23/... over list of hkl 
        d_out = dict(fn = np_fn_1d, sft = np_sft_2d)
        self.update(d_out)
        
        
    def set_vals(self, d_vals, refresh = False):
        """
        Set values 
        """
        keys = d_vals.keys()
        llab = ["crystal_symmetry", "unit_cell", "atom"]
        llab_in = [(hh in keys) for hh in llab]
        for lab, cond_h in zip(llab, llab_in):
            if cond_h:
                self[lab] = d_vals[lab]
                
        #redo it
        llab = ["resolution", "u", "v", "w", "x", "y"]
        llab_in = [(hh in keys) for hh in llab]
        cond = any(llab_in)
        if cond:
            d_val_as = {}
            for lab, cond_h in zip(llab, llab_in):
                if cond_h:
                    d_val_as[lab] = d_vals[lab]
            self["crystal_symmetry"].set_vals(d_val_as, refresh)

        #redo it
        llab = ["p1", "p2", "p3", "p4"]
        llab_in = [(hh in keys) for hh in llab]
        cond = any(llab_in)
        if cond:
            d_val_as = {}
            for lab, cond_h in zip(llab, llab_in):
                if cond_h:
                    d_val_as[lab] = d_vals[lab]
            self["unit_cell"].set_vals(d_val_as, refresh)

        #redo it
        llab = ["p1", "p2", "p3", "p4"]
        llab_in = [(hh in keys) for hh in llab]
        cond = any(llab_in)
        if cond:
            d_val_as = {}
            for lab, cond_h in zip(llab, llab_in):
                if cond_h:
                    d_val_as[lab] = d_vals[lab]
            self["atom"].set_vals(d_val_as, refresh)
            
        if refresh:
            tth = self["tth"]
            self.calc_model(tth)

            
    def soft_copy(self):
        """
        Soft copy of the object with saving the links on the internal parameter of the object
        """
        obj_new = Crystal(crystal_symmetry = self["crystal_symmetry"], 
                          unit_cell = self["unit_cell"], atom = self["atom"])
        llab = ["fn", "sft"]
        keys = self.keys()
        llab_in = [(hh in keys) for hh in llab]
        for lab, cond_h in zip(llab, llab_in):
            if cond_h:
                obj_new[lab] = self[lab]
        return obj_new


             
#import Variable
#v_u = Variable.Variable(0.2)

        
if (__name__ == "__main__"):
  pass
