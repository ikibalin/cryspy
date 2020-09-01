__author__ = 'ikibalin'
__version__ = "2019_12_10"
import os
import numpy

import warnings
from typing import List, Tuple

from cryspy.common.cl_item_constr import ItemConstr
from cryspy.common.cl_loop_constr import LoopConstr
from cryspy.common.cl_fitable import Fitable

class PdInstrReflexAsymmetry(ItemConstr):
    """
PdInstrReflexAsymmetry describes asymmetry of Bragg reflections for 1d powder diffractometer

Description in cif file::

 _pd_instr_reflex_asymmetry_p1 0.0
 _pd_instr_reflex_asymmetry_p2 0.0
 _pd_instr_reflex_asymmetry_p3 0.0
 _pd_instr_reflex_asymmetry_p4 0.0
    """
    MANDATORY_ATTRIBUTE = ("p1", "p2", "p3", "p4")
    OPTIONAL_ATTRIBUTE = ()
    INTERNAL_ATTRIBUTE = ()
    PREFIX = "pd_instr_reflex_asymmetry"
    def __init__(self, p1=None, p2=None, p3=None, p4=None):
        super(PdInstrReflexAsymmetry, self).__init__(mandatory_attribute=self.MANDATORY_ATTRIBUTE, 
                                                optional_attribute=self.OPTIONAL_ATTRIBUTE, 
                                                internal_attribute=self.INTERNAL_ATTRIBUTE,
                                                prefix=self.PREFIX)
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        self.p4 = p4

        if self.is_defined:
            self.form_object

    @property
    def p1(self):
        return getattr(self, "__p1")
    @p1.setter
    def p1(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__p1", x_in)

    @property
    def p2(self):
        return getattr(self, "__p2")
    @p2.setter
    def p2(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__p2", x_in)

    @property
    def p3(self):
        return getattr(self, "__p3")
    @p3.setter
    def p3(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__p3", x_in)

    @property
    def p4(self):
        return getattr(self, "__p4")
    @p4.setter
    def p4(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__p4", x_in)

    @property
    def is_variable(self) -> bool:
        """
Output: True if there is any refined parameter
        """
        res = any([self.p1.refinement, 
                   self.p2.refinement,
                   self.p3.refinement,
                   self.p4.refinement])
        return res        
    
    def get_variables(self) -> List:
        """
Output: the list of the refined parameters
        """
        l_variable = []
        if self.p1.refinement: l_variable.append(self.p1)
        if self.p2.refinement: l_variable.append(self.p2)
        if self.p3.refinement: l_variable.append(self.p3)
        if self.p4.refinement: l_variable.append(self.p4)
        return l_variable
        
    def _func_fa(self, tth):
        """
        For assymmetry correction F_a(z)
        """ 
        return 2*tth*numpy.exp(-tth**2)
        
    def _func_fb(self, tth):
        """
        For assymmetry correction F_b(z)
        """ 
        return 2.*(2.*tth**2-3.)* self._func_fa(tth)
        
    def calc_asymmetry(self, tth, tth_hkl, fwhm):
        """
Calculate asymmetry coefficients for  on the given list ttheta for 
bragg reflections flaced on the position ttheta_hkl
tth and tth_hkl in degrees

:IMPORTANT: THERE IS MISTAKE (look page 54 in FullProf Manual)
        """
        tth_2d, tth_hkl_2d = numpy.meshgrid(tth, tth_hkl, indexing="ij")
        np_zero = numpy.zeros(tth_2d.shape, dtype = float)
        np_one = numpy.ones(tth_2d.shape, dtype = float)
        val_1, val_2 = np_zero, np_zero
        
        z_2d = (tth_2d - tth_hkl_2d)/fwhm[numpy.newaxis, :]
        
        p1, p2 = float(self.p1), float(self.p2)
        p3, p4 = float(self.p3), float(self.p4)
        flag_1, flag_2 = False, False
        if ((p1!= 0.)|(p3!= 0.)):
            flag_1 = True
            fa = self._func_fa(z_2d)
        if ((p2!= 0.)|(p4!= 0.)):
            flag_2 = True
            fb = self._func_fb(z_2d)
            
        flag_3, flag_4 = False, False
        if ((p1!= 0.)|(p2!= 0.)):
            if flag_1:
                val_1 += p1*fa
                flag_3 = True
            if flag_2:
                val_1 += p2*fb
                flag_3 = True
            if flag_3:
                c1 = 1./numpy.tanh(0.5*tth_hkl)
                val_1 *= c1[numpy.newaxis, :]

        if ((p3!= 0.)|(p4!= 0.)):
            if flag_1:
                val_2 += p3*fa
                flag_4 = True
            if flag_2:
                val_2 += p4*fb
                flag_4 = True
            if flag_4:
                c2 = 1./numpy.tanh(tth_hkl)
                val_2 *= c2[numpy.newaxis, :]

        asymmetry_2d = np_one+val_1+val_2
        return asymmetry_2d
    
