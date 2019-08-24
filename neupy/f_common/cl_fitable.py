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
    @property
    def value(self):
        return self.__value
    @value.setter
    def value(self, x):
        try:
            self.__value = float(x)
        except:
            self.__value = None
    @property
    def sigma(self):
        return self.__sigma
    @sigma.setter
    def sigma(self, x):
        try:
            self.__sigma = float(x)
        except:
            self.__sigma = None
    @property
    def refinement(self):
        return self.__refinement
    @refinement.setter
    def refinement(self, x):
        try:
            self.__refinement = bool(x)
        except:
            self.__refinement = False
    @property
    def name(self):
        return self.__name
    @name.setter
    def name(self, x):
        if x is not None:
            try:
                self.__name = str(x)
            except:
                self.__name = None
        else:
            self.__name = None
    @property
    def constraint(self):
        return self.__constraint
    @constraint.setter
    def constraint(self, x):
        if x is not None:
            try:
                self.__constraint = str(x)
            except:
                self.__constraint = None
        else:
            self.__constraint = None
    @property
    def minimal(self):
        return self.__minimal
    @minimal.setter
    def minimal(self, x):
        try:
            self.__minimal = float(x)
        except:
            self.__minimal = None
    @property
    def maximal(self):
        return self.__maximal
    @maximal.setter
    def maximal(self, x):
        try:
            self.__maximal = float(x)
        except:
            self.__maximal = None
    def __array__(self):
        return numpy.array(self.value)
    def __float__(self):
        return self.value
    def __bool__(self):
        return self.refinement
    def __repr__(self):
        ls_out = [self.print_with_sigma()]
        if self.name is not None:
            ls_out.append("name: {:}".format(self.name))
        if self.minimal is not None:
            ls_out.append("minimal: {:.3f}".format(self.minimal))
        if self.maximal is not None:
            ls_out.append("maximal: {:.3f}".format(self.maximal))
        if self.refinement:
            ls_out.append("refinement: {:}".format(self.refinement))
        if self.constraint is not None:
            ls_out.append("constraint: {:}".format(self.constraint))
        res = ls_out[0]
        if len(ls_out) > 1:
            res = "{:} ({:})".format(res, ", ".join(ls_out[1:]))            
        return res
    def print_with_sigma(self):
        if not((self.sigma == 0.)|(self.sigma is None)):
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
    def __add__(self, var2):
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
    def __truediv__(self, var2):
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
    def __pos__(self):
        return self.value
    def __neg__(self):
        return -1*self.value
    def __abs__(self):
        return abs(self.value)
    def __round__(self, n):
        return round(self.value, n)
    def __pow__(self, var2):
        return self.value**var2
    def __rpow__(self, var2):
        return var2**self.value
    def __mod__(self, var2):
        return self.value%var2
    def __rmod__(self, var2):
        return var2%self.value
    def __floordiv__(self, var2):
        return self.value//var2
    def __rfloordiv__(self, var2):
        return var2//self.value
    def __ne__(self, var2):
        return  self.value != var2    
    def __rne__(self, var2):
        return  var2 != self.value  
    def __gt__(self, var2):
        return self.value > var2
    def __ge__(self, var2):
        return self.value >= var2
    def __lt__(self, var2):
        return  self.value  < var2
    def __le__(self, var2):
        return  self.value  <= var2
    def __and__(self, var2):
        return  (self.refinement & var2)
    def __rand__(self, var2):
        return  (var2 & self.refinement)
    def __or__(self, var2):
        return  (self.refinement | var2)
    def __ror__(self, var2):
        return  (var2 | self.refinement)

