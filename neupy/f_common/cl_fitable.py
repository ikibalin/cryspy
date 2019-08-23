"""
define the class Fitable.
"""
__author__ = 'ikibalin'
__version__ = "2019_08_22"

import numpy

class Fitable(object):
    """
    Fitable class to describe variable parameter
    """
    def __init__(self, value=0., sigma=0., refinement=False, 
                 name=None, constraint=None, minimal=None, maximal=None):
        super(Fitable, self).__init__()
        self.value = value
        self.sigma = sigma
        self.refinement = refinement
        self.name = name
        self.constraint = constraint
        self.minimal = minimal
        self.maximal = maximal
    def __setattr__(self, name, value):
        name_s = name.lower()
        if name_s in ["value", "sigma", "minimal", "maximal"]:
            func = float
        elif name_s in ["refinement"]:
            func = bool
        elif name_s in ["name", "constraint"]:
            func = str
        else:
            func = lambda x: x
        if value is not None:
            self.__dict__[name_s] = func(value)
        else:
            self.__dict__[name_s] = None
    def __array__(self):
        return array(self.value)

    def __float__(self):
        return self.value
    def __repr__(self):
        ls_out = ["Fitable object:"]
        if self.name is not None:
            ls_out.append("name: {:}".format(self.name))
        if self.value is not None:
            ls_out.append("value: {:.3f}".format(self.value))
        if self.sigma is not None:
            ls_out.append("sigma: {:.3f}".format(self.sigma))
        if self.minimal is not None:
            ls_out.append("minimal: {:.3f}".format(self.minimal))
        if self.maximal is not None:
            ls_out.append("maximal: {:.3f}".format(self.maximal))
        if self.refinement is not None:
            ls_out.append("refinement: {:}".format(self.refinement))
        if self.constraint is not None:
            ls_out.append("constraint: {:.3f}".format(self.constraint))
        return "\n".join(ls_out)
    def __pos__(self):
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
        elif isinstance(var2, Fitable):
            res = Fitable(value=self.value+var2.value,
                          refinement=(self.refinement|var2.refinement),
                          sigma=(self.sigma**2+var2.sigma**2)**0.5
                          )
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
        elif isinstance(var2, Fitable):
            res = Fitable(value=self.value+var2.value,
                          refinement=(self.refinement|var2.refinement),
                          sigma=(self.sigma**2+var2.sigma**2)**0.5
                          )
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
        elif isinstance(var2, Fitable):
            res = Fitable(value=self.value-var2.value,
                          refinement=(self.refinement|var2.refinement),
                          sigma=(self.sigma**2+var2.sigma**2)**0.5
                          )
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
        elif isinstance(var2, Fitable):
            res = Fitable(value=var2.value-self.value,
                          refinement=(self.refinement|var2.refinement),
                          sigma=(self.sigma**2+var2.sigma**2)**0.5
                          )
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
        elif isinstance(var2, Fitable):
            res = Fitable(value=var2.value*self.value,
                          refinement=(self.refinement|var2.refinement),
                          sigma=((self.sigma*var2.value)**2+(var2.sigma*self.value)**2)**0.5
                          )
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
        elif isinstance(var2, Fitable):
            res = Fitable(value=var2.value*self.value,
                          refinement=(self.refinement|var2.refinement),
                          sigma=((self.sigma*var2.value)**2+(var2.sigma*self.value)**2)**0.5
                          )
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
        elif isinstance(var2, Fitable):
            hh = (((var2.value*self.sigma)**2+(self.value*var2.sigma)**2)**0.5)/(var2.value)**2
            res = Fitable(value=self.value*1./var2.value, 
                          refinement=(self.refinement | var2.refinement),
                          sigma=hh)
        else:
            res = self.value*1./var2
        return res 
    def __rdiv__(self, var2):
        """
        output is float
        """
        if type(var2) is numpy.ndarray:
            res_1 = var2*1./self.value
            res = res_1.astype(var2.dtype)
        elif isinstance(var2, Fitable):
            hh = (((var2.value*self.sigma)**2+(self.value*var2.sigma)**2)**0.5)/(self.value)**2
            res = Fitable(value=var2.value*1./self.value, 
                          refinement=(self.refinement | var2.refinement),
                          sigma=hh)
        else:
            res = var2*1./self.value
        return res 
        return var2*1./self.value
    def __truediv__(self, var2):
        """
        output is float
        """
        if type(var2) is numpy.ndarray:
            res_1 = self.value*1./var2
            res = res_1.astype(var2.dtype)
        elif isinstance(var2, Fitable):
            hh = (((var2.value*self.sigma)**2+(self.value*var2.sigma)**2)**0.5)/(var2.value)**2
            res = Fitable(value=self.value*1./var2.value, 
                          refinement=(self.refinement | var2.refinement),
                          sigma=hh)
        else:
            res = self.value*1./var2
        return res 
    def __rtruediv__(self, var2):
        """
        output is float
        """
        if type(var2) is numpy.ndarray:
            res_1 = var2*1./self.value
            res = res_1.astype(var2.dtype)
        elif isinstance(var2, Fitable):
            hh = (((var2.value*self.sigma)**2+(self.value*var2.sigma)**2)**0.5)/(self.value)**2
            res = Fitable(value=var2.value*1./self.value, 
                          refinement=(self.refinement | var2.refinement),
                          sigma=hh)
        else:
            res = var2*1./self.value
        return res 
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

    def __float__(self):
        return float(self[0])

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


