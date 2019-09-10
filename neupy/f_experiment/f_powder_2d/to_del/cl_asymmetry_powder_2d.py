"""
define classes to calculated data
"""
__author__ = 'ikibalin'
__version__ = "2019_04_06"
import os
import numpy
import scipy.interpolate

#class BeamPolarization
from neupy.f_experiment.cl_beam_polarization import BeamPolarization
from neupy.f_common.cl_variable import Variable
#Description of setup class


class AsymmetryPowder2D(dict):
    """
    Asymmetry of the diffractometer
    """
    def __init__(self, p1 = 0., p2 = 0., p3 = 0., p4 = 0.):
        super(AsymmetryPowder2D, self).__init__()
        self._p_p1 = None
        self._p_p2 = None
        self._p_p3 = None
        self._p_p4 = None
        self._refresh(p1, p2, p3, p4)
        
    def __repr__(self):
        lsout = """AsymmetryPowder2D:\n p1: {:}\n p2: {:}\n p3: {:}
 p4: {:}""".format(self.get_val("p1"),  self.get_val("p2"),  
 self.get_val("p3"),  self.get_val("p4"))
        return lsout

    def _refresh(self, p1, p2, p3, p4):
        if not(isinstance(p1, type(None))):
            self._p_p1 = p1
        if not(isinstance(p2, type(None))):
            self._p_p2 = p2
        if not(isinstance(p3, type(None))):
            self._p_p3 = p3
        if not(isinstance(p4, type(None))):
            self._p_p4 = p4
            
    def set_val(self, p1=None, p2=None, p3=None, p4=None):
        self._refresh(p1, p2, p3, p4)
        
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
p1, p2, p3, p4 are coefficients to describe the assymetry shape for all 
               reflections like in FullProf
        """
        print(lsout)

        
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
        
        
        p1, p2 = self.get_val("p1"), self.get_val("p2")
        p3, p4 = self.get_val("p3"), self.get_val("p4")
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
    
    def is_variable(self):
        """
        without extinction
        """
        res = any([isinstance(self._p_p1, Variable), 
                   isinstance(self._p_p2, Variable), 
                   isinstance(self._p_p3, Variable), 
                   isinstance(self._p_p4, Variable)])
        return res        
    
    def get_variables(self):
        l_variable = []
        if isinstance(self._p_p1, Variable):
            l_variable.append(self._p_p1)
        if isinstance(self._p_p2, Variable):
            l_variable.append(self._p_p2)
        if isinstance(self._p_p3, Variable):
            l_variable.append(self._p_p3)
        if isinstance(self._p_p4, Variable):
            l_variable.append(self._p_p4)
        return l_variable
