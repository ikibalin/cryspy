"""
define classes to calculated data
"""
__author__ = 'ikibalin'
__version__ = "2019_04_06"
import os
import numpy

from neupy.f_common.cl_variable import Variable
    
#Description of setup class

class ResolutionPowder1D(dict):
    """
    Resoulution of the diffractometer
    """
    def __init__(self, u = 0, v = 0, w = 0.01, x = 0, y = 0):
        super(ResolutionPowder1D, self).__init__()
        self._p_u = None
        self._p_v = None
        self._p_w = None
        self._p_x = None
        self._p_y = None

        self._p_tan_th = None
        self._p_tan_th_sq = None
        self._p_cos_th = None
        self._p_hg = None
        self._p_hl = None
        self._p_hpv = None
        self._p_eta = None
        self._p_ag = None
        self._p_bg = None
        self._p_al = None
        self._p_bl = None
        
        self._refresh(u, v, w, x, y)
        
    def __repr__(self):
        lsout = """Resolution: 
 u: {:}\n v: {:}\n w: {:}\n x: {:}\n y: {:}""".format(self.get_val("u"),  
 self.get_val("v"), self.get_val("w"), self.get_val("x"), self.get_val("y"))
        return lsout

    def _refresh(self, u, v, w, x, y):
        if not(isinstance(u, type(None))):
            self._p_u = u
        if not(isinstance(v, type(None))):
            self._p_v = v
        if not(isinstance(w, type(None))):
            self._p_w = w
        if not(isinstance(x, type(None))):
            self._p_x = x
        if not(isinstance(y, type(None))):
            self._p_y = y
            
    def set_val(self, u=None, v=None, w=None, x=None, y=None):
        self._refresh(u, v, w, x, y)
        
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
u, v, w are coefficients to describe the Gauss part of peak shape

x, y are coefficients to describe the Lorentz part of peak shape
        """
        print(lsout)

    
    def _calc_tancos(self, th_hkl):
        """
        tth_hkl in radianas
        calculate tangenth (theta)
        """
        self._p_t_th = numpy.tan(th_hkl)
        self._p_t_th_sq = self._p_t_th**2
        res = numpy.cos(th_hkl)

        self._p_c_th = res
        self._p_ic_th = 1./res
        
    def _calc_hg(self, i_g = 0.):
        """
        ttheta in radians, could be array
        gauss size
        """
        u, v, w = 1.*self._p_u, 1.*self._p_v, 1.*self._p_w
        res_sq = (u*self._p_t_th_sq + v*self._p_t_th + w + 
                  i_g*self._p_ic_th**2)
        self._p_hg = numpy.sqrt(res_sq)
        
    def _calc_hl(self):
        """
        ttheta in radians, could be array
        lorentz site
        """
        x, y = self._p_x, self._p_y
        self._p_hl = x*self._p_t_th + y*self._p_ic_th


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

        self._p_ag = (2./hpv)*(numpy.log(2.)/numpy.pi)**0.5
        self._p_bg = 4*numpy.log(2)/(hpv**2)
        
    def _calc_albl(self):
        hpv = self._p_hpv
        self._p_al = 2./(numpy.pi*hpv )
        self._p_bl = 4./(hpv**2)
    
    def calc_resolution(self, tth_hkl, i_g = 0.):
        """
        Calculate parameters for tth
        tth_hkl in degrees
        """
        self._calc_tancos(0.5*tth_hkl*numpy.pi/180.)
        self._calc_hg(i_g = i_g)
        self._calc_hl()
        self._calc_hpveta()
        self._calc_agbg()
        self._calc_albl()

        a_g = self._p_ag  
        b_g = self._p_bg  
        a_l = self._p_al  
        b_l = self._p_bl 
        h_g = self._p_hg  
        h_l = self._p_hl 
        h_pv = self._p_hpv 
        eta = self._p_eta
        
        return h_pv, eta, h_g, h_l, a_g, b_g, a_l, b_l
    
    def is_variable(self):
        """
        without extinction
        """
        res = any([isinstance(self._p_u, Variable), 
                   isinstance(self._p_v, Variable), 
                   isinstance(self._p_w, Variable), 
                   isinstance(self._p_x, Variable), 
                   isinstance(self._p_y, Variable)])
        return res        

    def get_variables(self):
        l_variable = []
        if isinstance(self._p_u, Variable):
            l_variable.append(self._p_u)
        if isinstance(self._p_v, Variable):
            l_variable.append(self._p_v)
        if isinstance(self._p_w, Variable):
            l_variable.append(self._p_w)
        if isinstance(self._p_x, Variable):
            l_variable.append(self._p_x)
        if isinstance(self._p_y, Variable):
            l_variable.append(self._p_y)
        return l_variable


