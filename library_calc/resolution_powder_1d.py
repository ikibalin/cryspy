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
        
    

class Resolution(dict):
    """
    Resoulution of the diffractometer
    """
    def __init__(self, U = 0, V = 0, W = 0, Ig = 0, X = 0, Y = 0):
        super(Resolution, self).__init__()
        self._p_tan_tth = None
        self._p_tan_tth_sq = None
        self._p_cos_tth = None
        self._p_hg = None
        self._p_hl = None
        self._p_hpv = None
        self._p_eta = None
        self._p_ag = None
        self._p_bg = None
        self._p_al = None
        self._p_bl = None

        dd= {"U":U, "V":V, "W":W, "Ig":Ig, "X":X, "Y":Y}
        self.update(dd)
        
    def __repr__(self):
        lsout = """Resolution: 
 U {:}\n V {:}\n W {:}""".format(self["U"],  self["V"],  self["W"])
        return lsout
    
    def _calc_tancos(self, tth):
        self._p_t_tth = numpy.tan(tth)
        self._p_t_tth_sq = self._p_tan_tth**2
        res = numpy.cos(tth)
        self._p_c_tth = res
        self._p_ic_tth = 1./res
        
    def _calc_hg(self, tth):
        """
        ttheta in radians, could be array
        gauss size
        """
        res_sq = (self["U"]*self._p_t_tth_sq + self["V"]*self._p_t_tth + 
                  self["W"]*self._p_c_tth + self["Ig"]*self._p_ic_tth**2)
        self._p_hg = numpy.sqrt(res_sq)
        
    def _calc_hl(self, tth):
        """
        ttheta in radians, could be array
        lorentz site
        """
        self._p_hl = self["X"]*self._p_t_tth + self["Y"]*self._p_ic_tth


    def _calc_hpveta(self):
        """
        ttheta in radians, could be array
        pseudo-Voight function
        """
        hg = self._p_hg
        hl = self._p_hl
        hg_2, hl_2 = hg**2, hl**2
        hg_3, hl_3 = hg_2*hg, hl_2*hl
        hg_4, hl_4 = hg_3*hg, hl_3*hl
        hg_5, hl_5 = hg_4*hg, hl_4*hl
        c_2, c_3, c_4, c_5 = 2.69269, 2.42843, 4.47163, 0.07842
        hpv = (hg_5 + c_2*hg_4*hl + c_3*hg_3*hl_2 + 
                       c_4*hg_2*hl_3 + c_5*hg*hl_4 + hl_5)**0.2
        hh = hl*1./hpv
        self._p_hpv = hpv 
        self._p_eta = 1.36603*hh - 0.47719*hh**2 + 0.11116*hh**3

    def _calc_agbg(self):
        hpv = self._p_hpv
        self._p_ag = (2./hpv)*(math.log(2.)/math.pi)**0.5
        self._p_bg = 4*math.log(2)/(hpv**2)
        
    def _calc_albl(self):
        hpv = self._p_hpv
        self._p_al = 2./(math.pi*hpv )
        self._p_bl = 4./(hpv**2)
    
    def calc_param(self, tth, U = None, V = None, W = None, Ig = None, 
                    X = None, Y = None):
        """
        Calculate parameters for tth
        """
        if U != None:
            self["U"] = U
        if V != None:
            self["V"] = V
        if W != None:
            self["W"] = W
        if Ig != None:
            self["Ig"] = Ig
        if X != None:
            self["X"] = X
        if Y != None:
            self["Y"] = Y
            
        self._calc_tancos(tth)
        self._calc_hg(tth)
        self._calc_hl(tth)
        self._calc_hpveta()
        self._calc_agbg()
        self._calc_albl()
        d_out = dict(ag = self._p_ag, bg = self._p_bg, 
                     al = self._p_al, bl = self._p_bl,
                     hg = self._p_hg, hl = self._p_hl,
                     hpv = self._p_hpv, eta = self._p_eta)
        return d_out


class PeakProfile(dict):
    """
    Peak Profile
    """
    def __init__(self, ag = 0, bg = 0, al = 0, bl = 0):
        super(PeakProfile, self).__init__()
        dd= {"ag":ag, "bg":bg, "al":al, "bl":bl}
        self.update(dd)
        
    def __repr__(self):
        lsout = """Profile: \n ag {:}\n bg {:}\n al {:}
 bl {:}""".format(self["ag"],  self["bg"],  self["al"],  self["bl"])
        return lsout
        
    def gauss_pd(self, tth):
        """
        one dimensional gauss powder diffraction
        """
        ag, bg = self["ag"], self["bg"]
        return ag*numpy.exp(-bg*tth**2)
        
    def lor_pd(self, tth):
        """
        one dimensional lorentz powder diffraction
        """
        al, bl = self["al"], self["bl"]
        return al*1./(1.+bl*tth**2)



class Assymetry(dict):
    """
    Assymetry of the diffractometer
    """
    def __init__(self, p1 = 0, p2 = 0, p3 = 0, p4 = 0):
        super(Assymetry, self).__init__()
        dd= {"p1":p1, "p2":p2, "p3":p3, "p4":p4}
        self.update(dd)
        
    def __repr__(self):
        lsout = """Assymetry: \n p1 {:}\n p2 {:}\n p3 {:}
 p4 {:}""".format(self["p1"],  self["p2"],  self["p3"],  self["p4"])
        return lsout
        
    def _func_fa(self, tth):
        """
        for assymmetry correction
        """ 
        return 2*tth*numpy.exp(-tth**2)
        
    def _func_fb(self, tth):
        """
        for assymmetry correction
        """ 
        return 2.*(2.*tth**2-3.)* self._func_fa(tth)
        
    def calc_assym(self, tth, tth_hkl):
        fa = self._func_fa(tth)
        fb = self._func_fb(tth)
        res = 1.+(p1*fa+p2*fb)*1./numpy.tanh(0.5*tth_hkl)+(p3*fa+p4*fb)*1./numpy.tanh(tth_hkl)
        return res

    

cc = Variable(2,True,"sqdf sfd")
bb = Resolution(W=2,U=cc,V=3)    
cc[0]=74
bb["U"]
bb["U"]
print (bb.calc_hg(41))
    
if (__name__ == "__main__"):
  pass
