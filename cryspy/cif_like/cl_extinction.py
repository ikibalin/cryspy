__author__ = 'ikibalin'
__version__ = "2019_12_10"
import os
import numpy

import warnings
from typing import List, Tuple


from cryspy.common.cl_item_constr import ItemConstr
from cryspy.common.cl_loop_constr import LoopConstr
from cryspy.common.cl_fitable import Fitable

class Extinction(ItemConstr):
    """
Describe the extinction.

Description in cif file::

 _extinction_mosaicity 100.0
 _extinction_radius    50.0
 _extinction_model     gauss
    """
    MANDATORY_ATTRIBUTE = ("mosaicity", "radius", "model")
    OPTIONAL_ATTRIBUTE = ()
    INTERNAL_ATTRIBUTE = ()
    ACCESIBLE_MODEL = ("gauss", "lorentz")
    PREFIX = "extinction"
    def __init__(self, mosaicity=None, radius=None, model=None):
        super(Extinction, self).__init__(mandatory_attribute=self.MANDATORY_ATTRIBUTE, 
                                         optional_attribute=self.OPTIONAL_ATTRIBUTE, 
                                         internal_attribute=self.INTERNAL_ATTRIBUTE,
                                         prefix=self.PREFIX)

        self.model = model
        self.radius = radius
        self.mosaicity = mosaicity

        if self.is_defined:
            self.form_object

    @property
    def radius(self):
        """
The domain radius in angstrem
        """
        return getattr(self, "__radius")
    @radius.setter
    def radius(self, x):
        if x is None:
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__radius", x_in)

    @property
    def mosaicity(self):
        """
The domains' mosaicity in angstrem
        """
        return getattr(self, "__mosaicity")
    @mosaicity.setter
    def mosaicity(self, x):
        if x is None:
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__mosaicity", x_in)


    @property
    def model(self):
        """
The model of extinction.

:Accesible paramters: - gauss
                      - lorentz
        """
        return getattr(self, "__model")
    @model.setter
    def model(self, x):
        if x is None:
            x_in = None
        else:
            x_in = str(x)
            if not(x_in in self.ACCESIBLE_MODEL):
                warnings.warn(f"model '{x_in:}' is not supported", UserWarning, stacklevel=2)
                x_in = None
        setattr(self, "__model", x_in)

    def __repr__(self):
        ls_out = []
        ls_out.append(f"Extinction:\n{str(self):}")
        return "\n".join(ls_out)

    def _show_message(self, s_out: str):
        warnings.warn("***  Error ***\n"+s_out, UserWarning, stacklevel=2)
        
    @property
    def is_variable(self) -> bool:
        """
Output: True if there is any refined parameter
        """
        res = any([self.radius.refinement, 
                   self.mosaicity.refinement])
        return res

    def get_variables(self) -> List:
        """
Output: the list of the refined parameters
        """
        l_variable = []
        if self.radius.refinement: l_variable.append(self.radius)
        if self.mosaicity.refinement: l_variable.append(self.mosaicity)
        return l_variable



    def calc_extinction(self, cell, h, k, l, f_sq, wavelength, flag_derivative_f_sq=False):
        """
        f_sq in 10-12cm
        extinction for spherical model
        """
        r = float(self.radius)
        g = float(self.mosaicity)
        model = self.model
        kk = 1.
        vol = cell.volume
        sthovl = cell.calc_sthovl(h=h, k=k, l=l)
        stheta = sthovl * wavelength

        s2theta = 2. * stheta * (1. - stheta**2)**0.5
        c2theta = 1. - 2. * stheta**2
    
        q = (f_sq*kk/vol**2)*(wavelength**3)*1./s2theta
        delta_f_sq_q = (kk/vol**2)*(wavelength**3)*1./s2theta

        t = 1.5*r
        alpha = 1.5*r*s2theta*1./wavelength
        x = 2./3*q*alpha*t
        delta_f_sq_x = 2./3*delta_f_sq_q*alpha*t

        A = 0.20 + 0.45 * c2theta
        B = 0.22 - 0.12 * (0.5-c2theta)**2
        yp = (1.+2.*x+(A*x**2)*1./(1.+B*x))**(-0.5)
        delta_f_sq_yp = -0.5*(yp**3)*(2.*delta_f_sq_x+
                          2.*A*delta_f_sq_x*x*1./(1.+B*x)-
                          B*delta_f_sq_x*(A*x**2)*1./((1.+B*x)**2)) 
        
        
        ag = numpy.zeros(h.shape, dtype=float)
        al = numpy.zeros(h.shape, dtype=float)
        
        flag = alpha != 0.
        ag[flag] = alpha[flag]*g*(g**2+0.5*alpha[flag]**2)**(-0.5)
        al[flag] = alpha[flag]*g*1./(g+alpha[flag]*2./3.)
        
        if model == "gauss":
            xs = 2./3.*q*ag*t
            delta_f_sq_xs = 2./3.*delta_f_sq_q*ag*t
            A = 0.58 + 0.48 * c2theta + 0.24 * c2theta**2
            B = 0.02 - 0.025 * c2theta
            #print("A, B", A, B)
            ys = (1+2.12*xs+(A*xs**2)*1./(1+B*xs))**(-0.5)
            delta_f_sq_ys = -0.5*(ys**3)*(2.12*delta_f_sq_xs+2.*(A*xs*delta_f_sq_xs)*1./(1+B*xs)-B*delta_f_sq_xs*(A*xs**2)*1./((1+B*xs)**2))
        elif model == "lorentz":
            xs = 2./3.*q*al*t
            delta_f_sq_xs = 2./3.*delta_f_sq_q*al*t

            A = 0.025 + 0.285 * c2theta
            B = -0.45 * c2theta
            flag = c2theta>0
            B[flag] = 0.15 - 0.2 * (0.75-c2theta[flag])**2
            ys = (1.+2.*xs+(A*xs**2)*1./(1+B*xs))**(-0.5)
            delta_f_sq_ys = -0.5*(ys**3)*(2.*delta_f_sq_xs+2.*(A*xs*delta_f_sq_xs)*1./(1+B*xs)-B*delta_f_sq_xs*(A*xs**2)*1./((1.+B*xs)**2))
        else:
            ys = 1.
            delta_f_sq_ys = 0.
        #print("ys", ys)
        yext = yp * ys
        delta_f_sq_yext = delta_f_sq_yp * ys + yp * delta_f_sq_ys
        if flag_derivative_f_sq:
            return yext, delta_f_sq_yext 
        else:  
            return yext
    