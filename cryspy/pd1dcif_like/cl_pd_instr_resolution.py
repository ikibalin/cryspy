__author__ = 'ikibalin'
__version__ = "2019_12_10"
import os
import numpy

import warnings
from typing import List, Tuple

from cryspy.common.cl_item_constr import ItemConstr
from cryspy.common.cl_loop_constr import LoopConstr
from cryspy.common.cl_fitable import Fitable

class PdInstrResolution(ItemConstr):
    """
PdInstrResolution describes resolution of 1d powder diffractometer

Description in cif file::

 _pd_instr_resolution_u 16.9776
 _pd_instr_resolution_v -2.8357
 _pd_instr_resolution_w  0.5763
 _pd_instr_resolution_x  0.0000
 _pd_instr_resolution_y  0.0000
    """
    MANDATORY_ATTRIBUTE = ("u", "v", "w", "x", "y")
    OPTIONAL_ATTRIBUTE = ()
    INTERNAL_ATTRIBUTE = ()
    PREFIX = "pd_instr_resolution"
    def __init__(self, u=None, v=None, w=None, x=None, y=None):
        super(PdInstrResolution, self).__init__(mandatory_attribute=self.MANDATORY_ATTRIBUTE, 
                                                optional_attribute=self.OPTIONAL_ATTRIBUTE, 
                                                internal_attribute=self.INTERNAL_ATTRIBUTE,
                                                prefix=self.PREFIX)

        self.u = u
        self.v = v
        self.w = w
        self.x = x
        self.y = y

        if self.is_defined:
            self.form_object

        # self.__tan_th = None
        # self.__tan_th_sq = None
        # self.__cos_th = None
        # self.__hg = None
        # self.__hl = None
        # self.__hpv = None
        # self.__eta = None
        # self.__ag = None
        # self.__bg = None
        # self.__al = None
        # self.__bl = None

    @property
    def u(self):
        return getattr(self, "__u")
    @u.setter
    def u(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__u", x_in)

    @property
    def v(self):
        return getattr(self, "__v")
    @v.setter
    def v(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__v", x_in)

    @property
    def w(self):
        return getattr(self, "__w")
    @w.setter
    def w(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__w", x_in)

    @property
    def x(self):
        return getattr(self, "__x")
    @x.setter
    def x(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__x", x_in)

    @property
    def y(self):
        return getattr(self, "__y")
    @y.setter
    def y(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__y", x_in)


    def __repr__(self):
        ls_out = []
        ls_out.append(f"PdInstrResolution:\n{str(self):}")
        return "\n".join(ls_out)


    def _show_message(self, s_out: str):
        warnings.warn("***  Error ***\n"+s_out, UserWarning, stacklevel=2)

    @property
    def is_variable(self) -> bool:
        """
Output: True if there is any refined parameter
        """
        res = any([self.u.refinement, 
                   self.v.refinement,
                   self.w.refinement,
                   self.x.refinement,
                   self.y.refinement])
        return res        
    
    def get_variables(self) -> List:
        """
Output: the list of the refined parameters
        """
        l_variable = []
        if self.u.refinement: l_variable.append(self.u)
        if self.v.refinement: l_variable.append(self.v)
        if self.w.refinement: l_variable.append(self.w)
        if self.x.refinement: l_variable.append(self.x)
        if self.y.refinement: l_variable.append(self.y)
        return l_variable

    
    def _calc_tancos(self, th_hkl):
        """
tth_hkl in radianas
calculate tangenth (theta)
        """
        self.__t_th = numpy.tan(th_hkl)
        self.__t_th_sq = self.__t_th**2
        res = numpy.cos(th_hkl)

        self.__c_th = res
        self.__ic_th = 1./res
        
    def _calc_hg(self, i_g = 0.):
        """
        ttheta in radians, could be array
        gauss size
        """
        u, v, w = float(self.u), float(self.v), float(self.w)
        res_sq = (u*self.__t_th_sq + v*self.__t_th + w + 
                  i_g*self.__ic_th**2)
        self.__hg = numpy.sqrt(res_sq)
        
    def _calc_hl(self):
        """
        ttheta in radians, could be array
        lorentz site
        """
        x, y = float(self.x), float(self.y)
        self.__hl = x*self.__t_th + y*self.__ic_th


    def _calc_hpveta(self):
        """
        ttheta in radians, could be array
        pseudo-Voight function
        """
        hg = self.__hg
        hl = self.__hl

        hg_2, hl_2 = hg**2, hl**2
        hg_3, hl_3 = hg_2*hg, hl_2*hl
        hg_4, hl_4 = hg_3*hg, hl_3*hl
        hg_5, hl_5 = hg_4*hg, hl_4*hl
        c_2, c_3, c_4, c_5 = 2.69269, 2.42843, 4.47163, 0.07842
        hpv = (hg_5 + c_2*hg_4*hl + c_3*hg_3*hl_2 + 
                       c_4*hg_2*hl_3 + c_5*hg*hl_4 + hl_5)**0.2
        hh = hl*1./hpv
        self.__hpv = hpv 
        self.__eta = 1.36603*hh - 0.47719*hh**2 + 0.11116*hh**3

    def _calc_agbg(self):
        hpv = self.__hpv

        self.__ag = (2./hpv)*(numpy.log(2.)/numpy.pi)**0.5
        self.__bg = 4*numpy.log(2)/(hpv**2)
        
    def _calc_albl(self):
        hpv = self.__hpv
        self.__al = 2./(numpy.pi*hpv )
        self.__bl = 4./(hpv**2)
    
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

        a_g = self.__ag  
        b_g = self.__bg  
        a_l = self.__al  
        b_l = self.__bl 
        h_g = self.__hg  
        h_l = self.__hl 
        h_pv = self.__hpv 
        eta = self.__eta
        
        return h_pv, eta, h_g, h_l, a_g, b_g, a_l, b_l
    
