"""
define classes to describe crystal 
"""
__author__ = 'ikibalin'
__version__ = "2019_04_06"
import os
import numpy

from neupy.f_common.cl_variable import Variable
from neupy.f_interface.cl_abstract_extinction import AbstractExtinction


class Extinction(AbstractExtinction):
    """
    Extinction
    """
    def __init__(self, domain_radius=0., mosaicity=0., model="gauss"):
        super(Extinction, self).__init__()
        self._p_domain_radius = None
        self._p_mosaicity = None
        self._p_model = None
        self._refresh(domain_radius, mosaicity, model)

    def __repr__(self):
        lsout = """Extinction: \n model: {:}\n domain_radius: {:}
 mosaicity: {:}""".format(self._p_model, self._p_domain_radius, 
 self._p_mosaicity)
        return lsout


    def _refresh(self, domain_radius, mosaicity, model):
        
        if domain_radius is not None:
            self._p_domain_radius = domain_radius
        if mosaicity is not None:
            self._p_mosaicity = mosaicity
        if model is not None:
            self._p_model = model

    def set_val(self, domain_radius=None, mosaicity=None, model=None):
        self._refresh(domain_radius, mosaicity, model)
        
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
model can be "gauss" of "lorentz" model
domain_radius is the averaged radius of domains in Angstrems (???)
mosaicity is the mosaicity of domains in Minutes (???)
        """
        print(lsout)
        
    def calc_extinction(self, cell, h, k, l, f_sq, wave_length):
        """
        f_sq in 10-12cm
        extinction for spherical model
        """
        r = 1.*self._p_domain_radius
        g = 1.*self._p_mosaicity
        model = self._p_model
        kk = 1.
        vol = cell.get_val("vol")
        sthovl = cell.calc_sthovl(h, k, l)
        stheta = sthovl * wave_length

        s2theta = 2. * stheta * (1. - stheta**2)**0.5
        c2theta = 1. - 2. * stheta**2
    
        q = (f_sq*kk/vol**2)*(wave_length**3)*1./s2theta
    
        t = 1.5*r
        alpha = 1.5*r*s2theta*1./wave_length
        x = 2./3*q*alpha*t
    
        A = 0.20 + 0.45 * c2theta
        B = 0.22 - 0.12 * (0.5-c2theta)**2
        yp = (1.+2.*x+(A*x**2)*1./(1.+B*x))**(-0.5)
        
        ag = numpy.zeros(h.shape, dtype=float)
        al = numpy.zeros(h.shape, dtype=float)
        
        flag = alpha != 0.
        ag[flag] = alpha[flag]*g*(g**2+0.5*alpha[flag]**2)**(-0.5)
        al[flag] = alpha[flag]*g*1./(g+alpha[flag]*2./3.)
        
        if model == "gauss":
            xs = 2./3.*q*ag*t
            A = 0.58 + 0.48 * c2theta + 0.24 * c2theta**2
            B = 0.02 - 0.025 * c2theta
            #print("A, B", A, B)
            ys = (1+2.12*xs+(A*xs**2)*1./(1+B*xs))**(-0.5)
        elif model == "lorentz":
            xs = 2./3.*q*al*t
            A = 0.025 + 0.285 * c2theta
            B = -0.45 * c2theta
            flag = c2theta>0
            B[flag] = 0.15 - 0.2 * (0.75-c2theta[flag])**2
            ys = (1+2*xs+(A*xs**2)*1./(1+B*xs))**(-0.5)
        else:
            ys = 1.
        #print("ys", ys)
        yext = yp * ys
        return yext
            
    def is_variable(self):
        """
        without extinction
        """
        res = any([isinstance(self._p_domain_radius, Variable), 
                   isinstance(self._p_mosaicity, Variable)])
        return res

    def get_variables(self):
        l_variable = []
        if isinstance(self._p_domain_radius, Variable):
            l_variable.append(self._p_domain_radius)
        if isinstance(self._p_mosaicity, Variable):
            l_variable.append(self._p_mosaicity)
        return l_variable

