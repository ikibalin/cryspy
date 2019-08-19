"""
define classes to describe crystal 
"""
__author__ = 'ikibalin'
__version__ = "2019_04_06"
import os
import numpy

from neupy.f_common.cl_variable import Variable
from neupy.f_interface.cl_abstract_adp import AbstractADP


class ADP(AbstractADP):
    """
    ADP
    """
    def __init__(self, u_11 = 0., u_22 = 0., u_33 = 0., 
                 u_12 = 0., u_13 = 0., u_23 = 0., b_iso = 0.):
        super(ADP, self).__init__()
        self._p_u_11 = None
        self._p_u_22 = None
        self._p_u_33 = None
        self._p_u_12 = None
        self._p_u_13 = None
        self._p_u_23 = None
        self._p_b_iso = None
        self._refresh(u_11, u_22, u_33, u_12, u_13, u_23, b_iso)

    def __repr__(self):
        lsout = """Debye Waller: \n u_11: {:}, u_22: {:}, u_33: {:}
 u_12: {:}, u_13: {:}, u_23: {:}\n b_iso: {:}""".format(
 self._p_u_11, self._p_u_22, self._p_u_33, self._p_u_12, 
 self._p_u_13, self._p_u_23, self._p_b_iso)
        return lsout


    def _refresh(self, u_11, u_22, u_33, u_12, u_13, u_23, b_iso):
        
        if not(isinstance(u_11, type(None))):
            self._p_u_11 = u_11
        if not(isinstance(u_22, type(None))):
            self._p_u_22 = u_22
        if not(isinstance(u_33, type(None))):
            self._p_u_33 = u_33
        if not(isinstance(u_12, type(None))):
            self._p_u_12 = u_12
        if not(isinstance(u_13, type(None))):
            self._p_u_13 = u_13
        if not(isinstance(u_23, type(None))):
            self._p_u_23 = u_23
        if not(isinstance(b_iso, type(None))):
            self._p_b_iso = b_iso

    def set_val(self, u_11=None, u_22=None, u_33=None, u_12=None, 
                u_13=None, u_23=None, b_iso=None):
        self._refresh(u_11, u_22, u_33, u_12, u_13, u_23, b_iso)
        
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
beta_ij are Debye-Waller factor
b_iso is the isotropical Debye-Waller factor
        """
        print(lsout)
        
    def _calc_power_dwf_iso(self, sthovl):
        """
        isotropic harmonic Debye-Waller factor
        """
        b_iso = 1.*self._p_b_iso
        sthovl_sq = sthovl**2
        b_iso_2d, sthovl_sq_2d = numpy.meshgrid(sthovl_sq, b_iso, indexing="ij")
        
        power_dwf_iso_2d = b_iso_2d*sthovl_sq_2d
        return power_dwf_iso_2d

    def calc_power_dwf_aniso(self, space_group, cell, h, k, l):
        """
        anisotropic harmonic Debye-Waller factor
        
        h,k,l is 1D (temporary solution)
        """
        r_11, r_12 = space_group.get_val("r_11"), space_group.get_val("r_12")
        r_13, r_21 = space_group.get_val("r_13"), space_group.get_val("r_21")
        r_22, r_23 = space_group.get_val("r_22"), space_group.get_val("r_23")
        r_31, r_32 = space_group.get_val("r_31"), space_group.get_val("r_32")
        r_33 = space_group.get_val("r_33")
  
        b_11, b_22, b_33, b_12, b_13, b_23 = self.calc_beta(cell)

        np_h, np_b_11, np_r_11 = numpy.meshgrid(h, b_11, r_11, indexing="ij")
        np_k, np_b_22, np_r_22 = numpy.meshgrid(k, b_22, r_22, indexing="ij")
        np_l, np_b_33, np_r_33 = numpy.meshgrid(l, b_33, r_33, indexing="ij")
        np_h, np_b_12, np_r_12 = numpy.meshgrid(h, b_12, r_12, indexing="ij")
        np_h, np_b_13, np_r_13 = numpy.meshgrid(h, b_13, r_13, indexing="ij")
        np_h, np_b_23, np_r_23 = numpy.meshgrid(h, b_23, r_23, indexing="ij")
        np_r_21 = numpy.meshgrid(h, b_23, r_21, indexing="ij")[2]
        np_r_31 = numpy.meshgrid(h, b_23, r_31, indexing="ij")[2]
        np_r_32 = numpy.meshgrid(h, b_23, r_32, indexing="ij")[2]
        
        np_h_s = np_h*np_r_11 + np_k*np_r_21 + np_l*np_r_31
        np_k_s = np_h*np_r_12 + np_k*np_r_22 + np_l*np_r_32
        np_l_s = np_h*np_r_13 + np_k*np_r_23 + np_l*np_r_33
        
        power_dwf_aniso = (np_b_11*np_h_s**2 + np_b_22*np_k_s**2 + 
                       np_b_33*np_l_s**2 + 2.*np_b_12*np_h_s*np_k_s + 
                       2.*np_b_13*np_h_s*np_l_s + 2.*np_b_23*np_k_s*np_l_s)
        
        return power_dwf_aniso 
        
    def calc_dwf(self, space_group, cell, h, k, l):
        """
        calculate Debye-Waller factor
        """
        sthovl = cell.calc_sthovl(h, k, l)
        #dimensions (hkl, atoms in assymmetric unit cell)
        power_iso_2d = self._calc_power_dwf_iso(sthovl)
        #dimensions (hkl, atoms in assymmetric unit cell, el.symmetry)
        power_aniso_3d = self.calc_power_dwf_aniso(space_group, cell, h, k, l)
        power_3d = power_iso_2d[:, :, numpy.newaxis] + power_aniso_3d
        dwf_3d = numpy.exp(-power_3d)
        return dwf_3d
        
    def calc_beta(self, cell):
        """
        calculate beta_ij from U_ij
        """
        ia, ib, ic = 1.*cell.get_val("ia"), 1.*cell.get_val("ib"), 1.*cell.get_val("ic")
        u_11, u_22, u_33 = 1.*self._p_u_11, 1.*self._p_u_22, 1.*self._p_u_33
        u_12, u_13, u_23 = 1.*self._p_u_12, 1.*self._p_u_13, 1.*self._p_u_23
        
        beta_11 = 2.*numpy.pi**2 * u_11 *ia**2
        beta_22 = 2.*numpy.pi**2 * u_22 *ib**2
        beta_33 = 2.*numpy.pi**2 * u_33 *ic**2
        beta_12 = 2.*numpy.pi**2 * u_12 *ia*ib
        beta_13 = 2.*numpy.pi**2 * u_13 *ia*ic
        beta_23 = 2.*numpy.pi**2 * u_23 *ib*ic
        return beta_11, beta_22, beta_33, beta_12, beta_13, beta_23

