"""
define resolution for the powder diffractometer along ttheta
"""
__author__ = 'ikibalin'
__version__ = "2019_03_29$"
import numpy

class Variable():
    """
    general class for Variable
    """
    value = 0.
    refinement = False
    constraint = ""
    def __init__(self, val = 0., ref = False, constr = ""):
        super(Variable, self).__init__()
        self[0] = val 
        self[1] = ref
        self[2] = constr
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
        return self.value+var2
    def __radd__(self, var2):
        """
        output is float
        """
        return self.value+var2
    def __sub__(self, var2):
        """
        output is float
        """
        return self.value-var2
    def __rsub__(self, var2):
        """
        output is float
        """
        return var2-self.value
    def __mul__(self, var2):
        """
        output is float
        """
        return self.value*var2
    def __rmul__(self, var2):
        """
        output is float
        """
        return self.value*var2
    def __div__(self, var2):
        """
        output is float
        """
        return self.value*1./var2
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
        if i==2:
            return self.constraint
        elif i == 1:
            return self.refinement
        else:
            return self.value
    def __setitem__(self, i, v):
        #self.check(v)
        if i==2:
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
        lsout = """For variable with id {:}:\n    value: {:}
    refinement: {:}\n""".format(id(self), self.value, self.refinement)
        if self.constraint != "":
            lsout_add = "    constraint:{:}".format(self.constraint)
        else:
            lsout_add = "    constraint: None"
        return "".join([lsout, lsout_add])
        
    

class Resolution():
    """
    Resoulution of the diffractometer
    """
    U = 0.
    V = 0.
    W = 0.
    Ig = 0.
    X = 0.
    Y = 0.
    def __init__(self, U = 0, V = 0, W = 0, Ig = 0, X = 0, Y = 0):
        super(Resolution, self).__init__()
        dd= {"U":U, "V":V, "W":W, "Ig":Ig, "X":X, "Y":Y}
        self.__dict__() = dd
    
    def calc_hg(self, ttheta):
        """
        ttheta in radians, could be array
        """
        U = self["U"]
        V = self["V"]
        W = self["W"]
        Ig = self["Ig"]
        tan_ttheta = numpy.tan(ttheta)
        tan_ttheta_sq = tan_ttheta**2
        cos_ttheta = numpy.cos(ttheta)
        res_sq = U*tan_ttheta_sq+V*tan_ttheta+W+Ig/cos_ttheta
        return numpy.sqrt(res_sq)
    
bb = Variable("sf","qsf","")
cc = Variable(12,True,"sqdf sfd")
cc[0]=86
print(cc)
cc&True

bb-cc
cc-6
6-cc

bb.print()
bb = Resolution()    
bb["U"]
print (bb.calc_hg(41))
    
if (__name__ == "__main__"):
  pass