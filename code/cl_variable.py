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
    sigma = 0.
    refinement = False
    constraint = ""
    def __init__(self, val=0., ref=False, name="", constr="", sigma=0.):
        super(Variable, self).__init__()
        self[0] = val 
        self[1] = ref
        self[2] = constr
        self[3] = name
        self[4] = sigma
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
        return self.value%var2
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
        if i==4:
            return self.sigma
        elif i==3:
            return self.name
        elif i==2:
            return self.constraint
        elif i == 1:
            return self.refinement
        else:
            return self.value
    def __setitem__(self, i, v):
        #self.check(v)
        if i==4:
            cond1 = isinstance(v, float)
            cond2 = isinstance(v, int)
            if (cond1 | cond2):
                self.sigma = v 
            else:
                print ("Sigma should be integer or float type")
        elif i==3:
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
                print ("Variable should be integer or float type")
        return
    
    def __repr__(self):
        ls_out = self.print_with_sigma()
        
        #if self.name != "":
        #    ls_out += ", {:}".format(self.name)
            
        #id(self)
        #if self.constraint != "":
        #    lsout_add = "    constraint:{:}".format(self.constraint)
        #else:
        #    lsout_add = "    constraint: None"
        # "".join([lsout, lsout_add])
        return ls_out
    
    def print_with_sigma(self):
        if self.sigma != 0.:
            val_hh = numpy.log10(self.sigma)
            n_power = int(val_hh)
            if val_hh <= 0:
                val_1 = numpy.round(self.value, decimals = -1*int(n_power)+2)
                val_2 = numpy.round(self.sigma, decimals = -1*int(n_power)+2)
                i_val_11 = int(val_1)
                s_sign = ""
                if val_1 < 0.:
                    s_sign = "-"
                s_val_12 = ("{:}".format(int(abs(val_1)%1.*10**(-n_power+2)))).rjust(int(-n_power+2),"0")
                val_2=int(val_2*10**(-n_power+2))
                ls_out = " {:}{:}.{:}({:})".format(s_sign, abs(i_val_11), s_val_12, int(val_2))
            else:
                val_1 = numpy.round(self.value)
                val_2 = numpy.round(self.sigma)
                ls_out = " {:}({:})".format(int(val_1), int(val_2))
        else:
            ls_out = "{:}".format(self.value)
        return ls_out 


