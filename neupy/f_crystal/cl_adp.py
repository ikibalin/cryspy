"""
internal class to calculate Debye-Waller factor
"""
__author__ = 'ikibalin'
__version__ = "2019_08_29"
import os
import numpy

from neupy.f_common.cl_fitable import Fitable


class ADP(object):
    """
    ADP
    """
    def __init__(self, u_11=numpy.array([0], dtype=float), 
                       u_22=numpy.array([0], dtype=float), 
                       u_33=numpy.array([0], dtype=float), 
                       u_12=numpy.array([0], dtype=float), 
                       u_13=numpy.array([0], dtype=float), 
                       u_23=numpy.array([0], dtype=float), 
                       b_iso=numpy.array([0], dtype=float)):
        super(ADP, self).__init__()
        self.__atom_site_aniso_u_11 = None
        self.__atom_site_aniso_u_22 = None
        self.__atom_site_aniso_u_33 = None
        self.__atom_site_aniso_u_12 = None
        self.__atom_site_aniso_u_13 = None
        self.__atom_site_aniso_u_23 = None
        self.__atom_site_b_iso_or_equiv = None

        self.u_11 = u_11
        self.u_22 = u_22
        self.u_33 = u_33
        self.u_12 = u_12
        self.u_13 = u_13
        self.u_23 = u_23
        self.b_iso = b_iso

    def _trans_to_float_array(self, x):
        if isinstance(x, numpy.ndarray):
            x_out = x.astype(float)
        else:
            x_out = numpy.array([x], dtype=float)
        return x_out
        
    @property
    def u_11(self):
        """

        reference:
        """
        return self.__atom_site_aniso_u_11
    @u_11.setter
    def u_11(self, x):
        self.__atom_site_aniso_u_11 = self._trans_to_float_array(x)

    @property
    def u_22(self):
        """

        reference:
        """
        return self.__atom_site_aniso_u_22
    @u_22.setter
    def u_22(self, x):
        self.__atom_site_aniso_u_22 = self._trans_to_float_array(x)

    @property
    def u_33(self):
        """

        reference:
        """
        return self.__atom_site_aniso_u_33
    @u_33.setter
    def u_33(self, x):
        self.__atom_site_aniso_u_33 = self._trans_to_float_array(x)

    @property
    def u_12(self):
        """

        reference:
        """
        return self.__atom_site_aniso_u_12
    @u_12.setter
    def u_12(self, x):
        self.__atom_site_aniso_u_12 = self._trans_to_float_array(x)

    @property
    def u_13(self):
        """

        reference:
        """
        return self.__atom_site_aniso_u_13
    @u_13.setter
    def u_13(self, x):
        self.__atom_site_aniso_u_13 = self._trans_to_float_array(x)

    @property
    def u_23(self):
        """

        reference:
        """
        return self.__atom_site_aniso_u_23
    @u_23.setter
    def u_23(self, x):
        self.__atom_site_aniso_u_23 = self._trans_to_float_array(x)


    @property
    def b_iso(self):
        """

        reference:
        """
        return self.__atom_site_b_iso_or_equiv
    @b_iso.setter
    def b_iso(self, x):
        self.__atom_site_b_iso_or_equiv = self._trans_to_float_array(x)



    def __repr__(self):
        ls_out = ["Debye Waller:\n   U_iso    U_11     U_22     U_33    U_12     U_13     U_23"]
        ls_out.extend(["{:8.5f} {:8.5f} {:8.5f}{:8.5f} {:8.5f} {:8.5f} {:8.5f}".format(
            hh_1, hh_2, hh_3, hh_4, hh_5, hh_6, hh_7) for hh_1, hh_2, hh_3, hh_4, hh_5, hh_6, hh_7 
            in zip(self.b_iso, self.u_11, self.u_22, self.u_33, self.u_12, self.u_13, self.u_23)])
        return "\n".join(ls_out)

    def _calc_power_dwf_iso(self, sthovl):
        """
        isotropic harmonic Debye-Waller factor
        """
        b_iso = 1.*self.b_iso
        sthovl_sq = sthovl**2
        b_iso_2d, sthovl_sq_2d = numpy.meshgrid(sthovl_sq, b_iso, indexing="ij")
        
        power_dwf_iso_2d = b_iso_2d*sthovl_sq_2d
        return power_dwf_iso_2d

    def calc_power_dwf_aniso(self, space_group, cell, h, k, l):
        """
        anisotropic harmonic Debye-Waller factor
        
        h,k,l is 1D (temporary solution)
        """
        r_11, r_12 = space_group.r_11, space_group.r_12
        r_13, r_21 = space_group.r_13, space_group.r_21
        r_22, r_23 = space_group.r_22, space_group.r_23
        r_31, r_32 = space_group.r_31, space_group.r_32
        r_33 = space_group.r_33
  
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
        ia, ib, ic = 1.*cell.ia, 1.*cell.ib, 1.*cell.ic
        u_11, u_22, u_33 = 1.*self.u_11, 1.*self.u_22, 1.*self.u_33
        u_12, u_13, u_23 = 1.*self.u_12, 1.*self.u_13, 1.*self.u_23
        
        beta_11 = 2.*numpy.pi**2 * u_11 *ia**2
        beta_22 = 2.*numpy.pi**2 * u_22 *ib**2
        beta_33 = 2.*numpy.pi**2 * u_33 *ic**2
        beta_12 = 2.*numpy.pi**2 * u_12 *ia*ib
        beta_13 = 2.*numpy.pi**2 * u_13 *ia*ic
        beta_23 = 2.*numpy.pi**2 * u_23 *ib*ic
        return beta_11, beta_22, beta_33, beta_12, beta_13, beta_23

