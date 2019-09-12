"""
define class _pd2d_instr_reflex_asymmetry to describe asymmetry of Bragg reflections for 1d powder diffractometer
"""

__author__ = 'ikibalin'
__version__ = "2019_09_10"

import os
import numpy


from pystar import Data
from neupy.f_common.cl_fitable import Fitable



class Pd2dInstrReflexAsymmetry(dict):
    """
    Pd2dInstrReflexAsymmetry describes asymmetry of Bragg reflections for 2d powder diffractometer

    Example:

    _pd2d_instr_reflex_asymmetry_p1 0.0
    _pd2d_instr_reflex_asymmetry_p2 0.0
    _pd2d_instr_reflex_asymmetry_p3 0.0
    _pd2d_instr_reflex_asymmetry_p4 0.0
    """
    def __init__(self, p1 = Fitable(0.), p2 = Fitable(0.), 
                       p3 = Fitable(0.), p4 = Fitable(0.)):
        super(Pd2dInstrReflexAsymmetry, self).__init__()
        self.__pd2d_instr_reflex_asymmetry_p1 = None
        self.__pd2d_instr_reflex_asymmetry_p2 = None
        self.__pd2d_instr_reflex_asymmetry_p3 = None
        self.__pd2d_instr_reflex_asymmetry_p4 = None

        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        self.p4 = p4

    @property
    def p1(self):
        return self.__pd2d_instr_reflex_asymmetry_p1
    @p1.setter
    def p1(self, x):
        if isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        self.__pd2d_instr_reflex_asymmetry_p1 = x_in

    @property
    def p2(self):
        return self.__pd2d_instr_reflex_asymmetry_p2
    @p2.setter
    def p2(self, x):
        if isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        self.__pd2d_instr_reflex_asymmetry_p2 = x_in

    @property
    def p3(self):
        return self.__pd2d_instr_reflex_asymmetry_p3
    @p3.setter
    def p3(self, x):
        if isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        self.__pd2d_instr_reflex_asymmetry_p3 = x_in

    @property
    def p4(self):
        return self.__pd2d_instr_reflex_asymmetry_p4
    @p4.setter
    def p4(self, x):
        if isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        self.__pd2d_instr_reflex_asymmetry_p4 = x_in


    def __repr__(self):
        ls_out = ["Pd2dInstrResolution:"]
        ls_out.append(" p1: {:}".format(self.p1.print_with_sigma))
        ls_out.append(" p2: {:}".format(self.p2.print_with_sigma))
        ls_out.append(" p3: {:}".format(self.p3.print_with_sigma))
        ls_out.append(" p4: {:}".format(self.p4.print_with_sigma))
        return "\n".join(ls_out)

    @property
    def is_variable(self):
        """
        Output: True if there is any refined parameter
        """
        res = any([self.p1.refinement, 
                   self.p2.refinement,
                   self.p3.refinement,
                   self.p4.refinement])
        return res        
    
    def get_variables(self):
        """
        Output: the list of the refined parameters
        """
        l_variable = []
        if self.p1.refinement:
            l_variable.append(self.p1)
        if self.p2.refinement:
            l_variable.append(self.p2)
        if self.p3.refinement:
            l_variable.append(self.p3)
        if self.p4.refinement:
            l_variable.append(self.p4)
        return l_variable

    def _show_message(self, s_out: str):
        print("***  Error ***")
        print(s_out)
        
    @property
    def to_cif(self):
        ls_out = []
        ls_out.append("_pd2d_instr_reflex_asymmetry_p1 {:}".format(self.p1.print_with_sigma))
        ls_out.append("_pd2d_instr_reflex_asymmetry_p2 {:}".format(self.p2.print_with_sigma))
        ls_out.append("_pd2d_instr_reflex_asymmetry_p3 {:}".format(self.p3.print_with_sigma))
        ls_out.append("_pd2d_instr_reflex_asymmetry_p4 {:}".format(self.p4.print_with_sigma))
        return "\n".join(ls_out)


    def from_cif(self, string: str):
        cif_data = Data()
        flag = cif_data.take_from_string(string)
        if not flag:
            return False
        flag = False
        if cif_data.is_value("_pd2d_instr_reflex_asymmetry_p1"):
            self.p1 = cif_data["_pd2d_instr_reflex_asymmetry_p1"] # CIFvalue
        if cif_data.is_value("_pd2d_instr_reflex_asymmetry_p2"):
            self.p2 = cif_data["_pd2d_instr_reflex_asymmetry_p2"] # CIFvalue
        if cif_data.is_value("_pd2d_instr_reflex_asymmetry_p3"):
            self.p3 = cif_data["_pd2d_instr_reflex_asymmetry_p3"] # CIFvalue
        if cif_data.is_value("_pd2d_instr_reflex_asymmetry_p4"):
            self.p4 = cif_data["_pd2d_instr_reflex_asymmetry_p4"] # CIFvalue
        return True

        
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
        
    def calc_asymmetry(self, tth, tth_hkl):
        """
        Calculate asymmetry coefficients for  on the given list ttheta for 
        bragg reflections flaced on the position ttheta_hkl
        tth and tth_hkl in degrees
        
        IMPORTANT: THERE IS MISTAKE (look page 54 in FullProf Manual)
        """
        tth_2d, tth_hkl_2d = numpy.meshgrid(tth, tth_hkl, indexing="ij")
        np_zero = numpy.zeros(tth_2d.shape, dtype = float)
        np_one = numpy.ones(tth_2d.shape, dtype = float)
        val_1, val_2 = np_zero, np_zero
        
        
        p1, p2 = float(self.p1), float(self.p2)
        p3, p4 = float(self.p3), float(self.p4)
        flag_1, flag_2 = False, False
        if ((p1!= 0.)|(p3!= 0.)):
            flag_1 = True
            fa = self._func_fa(tth)
        if ((p2!= 0.)|(p4!= 0.)):
            flag_2 = True
            fb = self._func_fb(tth)
            
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
                c1_2d = numpy.meshgrid(tth, c1, indexing="ij")[1]
                val_1 *= c1_2d

        if ((p3!= 0.)|(p4!= 0.)):
            if flag_1:
                val_2 += p3*fa
                flag_4 = True
            if flag_2:
                val_2 += p4*fb
                flag_4 = True
            if flag_4:
                c2 = 1./numpy.tanh(tth_hkl)
                c2_2d = numpy.meshgrid(tth, c2, indexing="ij")[1]
                val_2 *= c2_2d

        asymmetry_2d = np_one+val_1+val_2
        return asymmetry_2d
    
