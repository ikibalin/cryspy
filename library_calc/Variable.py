"""
define the class variable.

It has the similar operation as a number, but also have some additional 
attributes
"""
__author__ = 'ikibalin'
__version__ = "2019_04_02$"

import numpy

class Variable(dict):
    """
    general class for Variable
    """
    value = 0.
    refinement = False
    constraint = ""
    def __init__(self, val=0., ref=False, name="", constr=""):
        super(Variable, self).__init__()
        self[0] = val 
        self[1] = ref
        self[2] = constr
        self[3] = name
    def __pos__(self):
        """
        output is float
        """
        return self.value
    def __neg__(self):
        """
        output is float
        """
        return -1*self.value
    def __abs__(self):
        """
        output is float
        """
        return abs(self.value)
    def __round__(self, n):
        """
        output is float
        """
        return round(self.value, n)
    def __add__(self, var2):
        """
        output is float
        """
        if type(var2) is numpy.ndarray:
            res_1 = self.value+var2
            res = res_1.astype(var2.dtype)
        else:
            res = self.value+var2
        return res
    def __radd__(self, var2):
        """
        output is float
        """
        if type(var2) is numpy.ndarray:
            res_1 = var2+self.value
            res = res_1.astype(var2.dtype)
        else:
            res = var2+self.value
        return res
    def __sub__(self, var2):
        """
        output is float
        """
        if type(var2) is numpy.ndarray:
            res_1 = self.value-var2
            res = res_1.astype(var2.dtype)
        else:
            res = self.value-var2
        return res
    def __rsub__(self, var2):
        """
        output is float
        """
        if type(var2) is numpy.ndarray:
            res_1 = var2-self.value
            res = res_1.astype(var2.dtype)
        else:
            res = var2-self.value
        return res
    def __mul__(self, var2):
        """
        output is float
        """
        if type(var2) is numpy.ndarray:
            res_1 = self.value*var2
            res = res_1.astype(var2.dtype)
        else:
            res = self.value*var2
        return res 
    def __rmul__(self, var2):
        """
        output is float
        """
        if type(var2) is numpy.ndarray:
            res_1 = self.value*var2
            res = res_1.astype(var2.dtype)
        else:
            res = self.value*var2
        return res
    def __div__(self, var2):
        """
        output is float
        """
        if type(var2) is numpy.ndarray:
            res_1 = self.value*1./var2
            res = res_1.astype(var2.dtype)
        else:
            res = self.value*1./var2
        return res 
    def __rdiv__(self, var2):
        """
        output is float
        """
        return var2*1./self.value
    def __pow__(self, var2):
        """
        output is float
        """
        return self.value**var2
    def __rpow__(self, var2):
        """
        output is float
        """
        return var2**self.value
    def __mod__(self, var2):
        """
        output is float
        """
        return s-elf.value%var2
    def __rmod__(self, var2):
        """
        output is float
        """
        return var2%self.value
    def __floordiv__(self, var2):
        """
        output is float
        """
        return self.value//var2
    def __rfloordiv__(self, var2):
        """
        output is float
        """
        return var2//self.value
    def __ne__(self, var2):
        """
        output is bool
        """
        return  self.value != var2    
    def __rne__(self, var2):
        """
        output is bool
        """
        return  var2 != self.value  
    def __gt__(self, var2):
        """
        output is bool
        """
        return self.value > var2
    def __ge__(self, var2):
        """
        output is bool
        """
        return self.value >= var2
    def __lt__(self, var2):
        """
        output is bool
        """
        return  self.value  < var2
    def __le__(self, var2):
        """
        output is bool
        """
        return  self.value  <= var2
    def __and__(self, var2):
        """
        output is bool
        """
        return  (self.refinement & var2)
    def __rand__(self, var2):
        """
        output is bool
        """
        return  (var2 & self.refinement)
    def __or__(self, var2):
        """
        output is bool
        """
        return  (self.refinement | var2)
    def __ror__(self, var2):
        """
        output is bool
        """
        return  (var2 | self.refinement)
    def __getitem__(self, i):
        if i==3:
            return self.name
        elif i==2:
            return self.constraint
        elif i == 1:
            return self.refinement
        else:
            return self.value
    def __setitem__(self, i, v):
        #self.check(v)
        if i==3:
            cond = isinstance(v, str)            
            if cond:
                self.name = v 
            else:
                print ("Name should be given as a string")
        elif i==2:
            cond = isinstance(v, str)            
            if cond:
                self.constraint = v 
            else:
                print ("Constraints should be given as a string")
        elif i == 1:
            cond = isinstance(v, bool)
            if cond:
                self.refinement = v 
            else:
                print ("Refinement should have a bool type variable")
        else:
            cond1 = isinstance(v, float)
            cond2 = isinstance(v, int)
            if (cond1 | cond2):
                self.value = v 
            else:
                print ("Variable should be integer of float type")
        return
    def __repr__(self):
        ls_out = " {:} (refinement: {:}".format(self.value, self.refinement)
        if self.name != "":
            ls_out += ", name: {:})".format(self.name)
        else:
            ls_out += ")"
        #id(self)
        #if self.constraint != "":
        #    lsout_add = "    constraint:{:}".format(self.constraint)
        #else:
        #    lsout_add = "    constraint: None"
        # "".join([lsout, lsout_add])
        return ls_out
        
    