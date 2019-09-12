"""
define classes to describe crystal 
"""
__author__ = 'ikibalin'
__version__ = "2019_09_02"
import os
import numpy

from pycifstar import Global
from neupy.f_common.cl_fitable import Fitable

class Extinction(object):
    """
    Extinction
    """
    def __init__(self, radius=Fitable(0.), 
                       mosaicity=Fitable(0.), 
                       model="gauss"):
        super(Extinction, self).__init__()
        self.__refine_ls_extinction_coef_radius = None
        self.__refine_ls_extinction_coef_mosaicity = None
        self.__refine_ls_extinction_model = None
        self.model = model
        self.radius = radius
        self.mosaicity = mosaicity

    def __repr__(self):
        ls_out = []
        ls_out.append("Extinction:")
        ls_out.append("model: {:}".format(self.model))
        ls_out.append("radius: {:}".format(self.radius))
        ls_out.append("mosaicity: {:}".format(self.mosaicity))
        return "\n".join(ls_out)

    @property
    def model(self):
        """
        The model of extinction

        """
        return self.__refine_ls_extinction_model
    @model.setter
    def model(self, x):
        l_model = ["gauss", "lorentz"]
        if x is None:
            x_in = "gauss"
        else:
            x_in = str(x).strip().lower()
            if not(x_in in l_model):
                x_in = "gauss"
                self._show_message("Can not recognize the introduced model '{:}'.\nThe supported models: {}.\nThe gauss model was choosen".format(
                                    x, " ".join(l_model)))
        self.__refine_ls_extinction_model = x_in

    @property
    def radius(self):
        """
        The domain radius in angstrem
        """
        return self.__refine_ls_extinction_coef_radius
    @radius.setter
    def radius(self, x):
        if isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        self.__refine_ls_extinction_coef_radius = x_in

    @property
    def mosaicity(self):
        """
        The domains' mosaicity in angstrem
        """
        return self.__refine_ls_extinction_coef_mosaicity
    @mosaicity.setter
    def mosaicity(self, x):
        if isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        self.__refine_ls_extinction_coef_mosaicity = x_in

    def calc_extinction(self, cell, h, k, l, f_sq, wavelength):
        """
        f_sq in 10-12cm
        extinction for spherical model
        """
        r = 1.*self.radius
        g = 1.*self.mosaicity
        model = self.model
        kk = 1.
        vol = cell.volume
        sthovl = cell.calc_sthovl(h=h, k=k, l=l)
        stheta = sthovl * wavelength

        s2theta = 2. * stheta * (1. - stheta**2)**0.5
        c2theta = 1. - 2. * stheta**2
    
        q = (f_sq*kk/vol**2)*(wavelength**3)*1./s2theta
    
        t = 1.5*r
        alpha = 1.5*r*s2theta*1./wavelength
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
            
    @property
    def is_variable(self):
        """
        without extinction
        """
        res = any([self.radius.refinement, 
                   self.mosaicity.refinement])
        return res

    def get_variables(self):
        l_variable = []
        if self.radius.refinement:
            l_variable.append(self.radius)
        if self.mosaicity.refinement:
            l_variable.append(self.mosaicity)
        return l_variable

    def _show_message(self, s_out: str):
        print("***  Error ***")
        print(s_out)

    @property
    def to_cif(self):
        ls_out = ["_refine_ls_extinction_coef_mosaicity {:}".format(self.mosaicity.print_with_sigma)]
        ls_out.append("_refine_ls_extinction_coef_radius {:}".format(self.radius.print_with_sigma))
        return "\n".join(ls_out)


    def from_cif(self, string: str):
        cif_data = Global()
        flag = cif_data.take_from_string(string)
        if not flag:
            return False
        flag = False
        if cif_data.is_value("_refine_ls_extinction_coef_mosaicity"):
            self.mosaicity = cif_data["_refine_ls_extinction_coef_mosaicity"] # CIFvalue
        if cif_data.is_value("_refine_ls_extinction_coef_radius"):
            self.radius = cif_data["_refine_ls_extinction_coef_radius"] # CIFvalue
        return True

