"""
define classes to describe _orient_matrix_UB and  _orient_matrix_U
"""
__author__ = 'ikibalin'
__version__ = "2019_09_02"
import os
import numpy
import copy

from pycifstar import Global
from neupy.f_common.cl_fitable import Fitable
from neupy.f_crystal.cl_cell import Cell


class OrientMatrix(object):
    """
    describe orientation matrix
    """
    def __init__(self,  u=numpy.array([[1., 0., 0.],[0., 1., 0.],[0., 0., 1.]], dtype=float),
                        cell=Cell(), ub_from=None, ub_ccsl=None, gamma_0=None, nu_0=None, wavelength=None):
        super(OrientMatrix, self).__init__()
        self.__orient_matrix_ub_11 = None
        self.__orient_matrix_ub_12 = None
        self.__orient_matrix_ub_13 = None
        self.__orient_matrix_ub_21 = None
        self.__orient_matrix_ub_22 = None
        self.__orient_matrix_ub_23 = None
        self.__orient_matrix_ub_31 = None
        self.__orient_matrix_ub_32 = None
        self.__orient_matrix_ub_33 = None

        self.__orient_matrix_u_11 = None
        self.__orient_matrix_u_12 = None
        self.__orient_matrix_u_13 = None
        self.__orient_matrix_u_21 = None
        self.__orient_matrix_u_22 = None
        self.__orient_matrix_u_23 = None
        self.__orient_matrix_u_31 = None
        self.__orient_matrix_u_32 = None
        self.__orient_matrix_u_33 = None

        self.__gamma_0 = None
        self.__nu_0 = None
        self.__wavelength = None
        
        flag_ub_from = ub_from is not None
        flag_ub_ccsl = ub_ccsl is not None

        if flag_ub_from:
            self.ub_from = ub_from
        elif flag_ub_ccsl:
            self.ub_ccsl = ub_ccsl
        else:
            self.cell = cell
            self.u = u
        self.gamma_0 = gamma_0
        self.nu_0 = nu_0
        self.wavelength = wavelength

    @property
    def gamma_0(self):
        return self.__gamma_0
    @gamma_0.setter
    def gamma_0(self, x):
        if x is None:
            x_in = 0.
        else:
            x_in = float(x)
        self.__gamma_0 = x_in

    @property
    def nu_0(self):
        return self.__nu_0
    @nu_0.setter
    def nu_0(self, x):
        if x is None:
            x_in = 0.
        else:
            x_in = float(x)
        self.__nu_0 = x_in

    @property
    def wavelength(self):
        return self.__wavelength
    @wavelength.setter
    def wavelength(self, x):
        if x is None:
            x_in = 1.4
        else:
            x_in = float(x)
        self.__wavelength = x_in



    @property
    def u_11(self):
        return self.__orient_matrix_u_11
    @property
    def u_12(self):
        return self.__orient_matrix_u_12
    @property
    def u_13(self):
        return self.__orient_matrix_u_13

    @property
    def u_21(self):
        return self.__orient_matrix_u_21
    @property
    def u_22(self):
        return self.__orient_matrix_u_22
    @property
    def u_23(self):
        return self.__orient_matrix_u_23

    @property
    def u_31(self):
        return self.__orient_matrix_u_31
    @property
    def u_32(self):
        return self.__orient_matrix_u_32
    @property
    def u_33(self):
        return self.__orient_matrix_u_33

    @property
    def ub_11(self):
        return self.__orient_matrix_ub_11
    @property
    def ub_12(self):
        return self.__orient_matrix_ub_12
    @property
    def ub_13(self):
        return self.__orient_matrix_ub_13

    @property
    def ub_21(self):
        return self.__orient_matrix_ub_21
    @property
    def ub_22(self):
        return self.__orient_matrix_ub_22
    @property
    def ub_23(self):
        return self.__orient_matrix_ub_23

    @property
    def ub_31(self):
        return self.__orient_matrix_ub_31
    @property
    def ub_32(self):
        return self.__orient_matrix_ub_32
    @property
    def ub_33(self):
        return self.__orient_matrix_ub_33


    @property
    def ub_ccsl(self):
        """
        The UB matrix (CCSL)
        """
        ub = numpy.array([[self.ub_11, self.ub_12, self.ub_13], 
                          [self.ub_21, self.ub_22, self.ub_23], 
                          [self.ub_31, self.ub_32, self.ub_33]], dtype=float)
        return ub
    @ub_ccsl.setter
    def ub_ccsl(self, x):
        if isinstance(x, numpy.ndarray):
            x_in = x
        else:
            try:
                x_in = numpy.array(x, dtype=float)
            except:
                self._show_message("The UB matrix should be numpy.ndarray with size(3x3)")
                return
        flag = x_in.shape == (3,3)
        if flag:
            self.__orient_matrix_ub_11 = float(x_in[0,0])
            self.__orient_matrix_ub_12 = float(x_in[0,1])
            self.__orient_matrix_ub_13 = float(x_in[0,2])
            self.__orient_matrix_ub_21 = float(x_in[1,0])
            self.__orient_matrix_ub_22 = float(x_in[1,1])
            self.__orient_matrix_ub_23 = float(x_in[1,2])
            self.__orient_matrix_ub_31 = float(x_in[2,0])
            self.__orient_matrix_ub_32 = float(x_in[2,1])
            self.__orient_matrix_ub_33 = float(x_in[2,2])
            self.__recalc_u_cell
        else:
                self._show_message("The UB matrix should be numpy.ndarray with size(3x3)")
                return

    @property
    def cell(self):
        """
        The cell
        """
        return self.__cell
    @cell.setter
    def cell(self, x):
        if isinstance(x, Cell):
            x_in = x
        else:
            self._show_message("A type of induced element is not recognized to convert it into Cell")
            x_in=Cell()
        self.__cell = x_in

    @property
    def ub_from(self):
        """
        The UB matrix (from) 
        """
        #not sure it should be checked
        ub = numpy.array([[-1.*self.ub_21, -1.*self.ub_22, -1.*self.ub_23], 
                          [self.ub_11, self.ub_12, self.ub_13], 
                          [self.ub_31, self.ub_32, self.ub_33]], dtype=float)
        return ub
    @ub_from.setter
    def ub_from(self, x):
        self.__from_ub_from_to_ub_ccsl(x)

    @property
    def u(self):
        """
        The U matrix (orientation matrix) 
        """
        #not sure it should be checked
        u = numpy.array([[self.u_11, self.u_12, self.u_13], 
                         [self.u_21, self.u_22, self.u_23], 
                         [self.u_31, self.u_32, self.u_33]], dtype=float)
        return u
    @u.setter
    def u(self, x):
        if isinstance(x, numpy.ndarray):
            x_in = x
        else:
            try:
                x_in = numpy.array(x, dtype=float)
            except:
                self._show_message("The U matrix should be numpy.ndarray with size(3x3)")
                return
        flag = x_in.shape == (3,3)
        if flag:
            self.__orient_matrix_u_11 = float(x_in[0,0])
            self.__orient_matrix_u_12 = float(x_in[0,1])
            self.__orient_matrix_u_13 = float(x_in[0,2])
            self.__orient_matrix_u_21 = float(x_in[1,0])
            self.__orient_matrix_u_22 = float(x_in[1,1])
            self.__orient_matrix_u_23 = float(x_in[1,2])
            self.__orient_matrix_u_31 = float(x_in[2,0])
            self.__orient_matrix_u_32 = float(x_in[2,1])
            self.__orient_matrix_u_33 = float(x_in[2,2])
            self.__recalc_ub
        else:
                self._show_message("The U matrix should be numpy.ndarray with size(3x3)")
                return



    def __from_ub_from_to_ub_ccsl(self, ub_from):
        ub_ccsl = numpy.array([[ub_from[1, 0], ub_from[1, 1], ub_from[1, 2]],
                               [-1.*ub_from[0, 0], -1.*ub_from[0, 1], -1.*ub_from[0, 2]],
                               [ub_from[2, 0], ub_from[2, 1], ub_from[2, 2]]], dtype=float)
        self.ub_ccsl = ub_ccsl


    @property
    def __recalc_u_cell(self):
        ub_from = self.ub_from
    
        v_b_1 = ub_from[:, 0]
        v_b_2 = ub_from[:, 1]
        v_b_3 = ub_from[:, 2]

        b_1 = float(numpy.sqrt(((v_b_1*v_b_1).sum())))
        b_2 = float(numpy.sqrt(((v_b_2*v_b_2).sum())))
        b_3 = float(numpy.sqrt(((v_b_3*v_b_3).sum())))

        c_i_g = (v_b_1*v_b_2).sum()/(b_1*b_2)
        c_i_b = (v_b_1*v_b_3).sum()/(b_1*b_3)
        c_i_a = (v_b_2*v_b_3).sum()/(b_2*b_3)

        i_a, i_b, i_g = numpy.arccos(c_i_a), numpy.arccos(c_i_b), numpy.arccos(c_i_g)

        s_i_a, s_i_b, s_i_g = numpy.sin(i_a), numpy.sin(i_b), numpy.sin(i_g)
    
        i_vol=b_1*b_2*b_3*(1-c_i_a**2-c_i_b**2-c_i_g**2+2.*c_i_a*c_i_b*c_i_g)**0.5

        c_a = (c_i_b*c_i_g-c_i_a)/(s_i_b*s_i_g)
        c_b = (c_i_a*c_i_g-c_i_b)/(s_i_a*s_i_g)
        c_g = (c_i_a*c_i_b-c_i_g)/(s_i_a*s_i_b)

        a_1 = b_2*b_3*s_i_a/i_vol
        a_2 = b_1*b_3*s_i_b/i_vol
        a_3 = b_1*b_2*s_i_g/i_vol
        alpha, beta, gamma = numpy.arccos(c_a), numpy.arccos(c_b), numpy.arccos(c_g)

        cell = self.cell
        cell.a, cell.b, cell.c = a_1, a_2, a_3
        cell.alpha, cell.beta, cell.gamma = float(alpha*180./numpy.pi), float(beta*180./numpy.pi), float(gamma*180./numpy.pi),

        m_ib = cell.m_ib
        u = numpy.dot(ub_from, m_ib)
        
        self.__orient_matrix_u_11 = u[0, 0]
        self.__orient_matrix_u_12 = u[0, 1]
        self.__orient_matrix_u_13 = u[0, 2]
        self.__orient_matrix_u_21 = u[1, 0]
        self.__orient_matrix_u_22 = u[1, 1]
        self.__orient_matrix_u_23 = u[1, 2]
        self.__orient_matrix_u_31 = u[2, 0]
        self.__orient_matrix_u_32 = u[2, 1]
        self.__orient_matrix_u_33 = u[2, 2]

    @property
    def __recalc_ub(self):
        u = self.u
        cell = self.cell
        m_ib = cell.m_ib
        ub_from = numpy.dot(u, m_ib)
        ub_ccsl = numpy.array([[ub_from[1, 0], ub_from[1, 1], ub_from[1, 2]],
                       [-1.*ub_from[0, 0], -1.*ub_from[0, 1], -1.*ub_from[0, 2]],
                       [ub_from[2, 0], ub_from[2, 1], ub_from[2, 2]]], dtype=float)

        self.__orient_matrix_ub_11 = ub_ccsl[0, 0]
        self.__orient_matrix_ub_12 = ub_ccsl[0, 1]
        self.__orient_matrix_ub_13 = ub_ccsl[0, 2]
        self.__orient_matrix_ub_21 = ub_ccsl[1, 0]
        self.__orient_matrix_ub_22 = ub_ccsl[1, 1]
        self.__orient_matrix_ub_23 = ub_ccsl[1, 2]
        self.__orient_matrix_ub_31 = ub_ccsl[2, 0]
        self.__orient_matrix_ub_32 = ub_ccsl[2, 1]
        self.__orient_matrix_ub_33 = ub_ccsl[2, 2]


    def __repr__(self):
        u = self.u
        ub_ccsl = self.ub_ccsl
        ub_from = self.ub_from
        cell = self.cell
        wavelength = self.wavelength
        gamma_0 = self.gamma_0
        nu_0 = self.nu_0
        ls_out = ["OrientMatrix:"]
        if wavelength is not None:
            ls_out.append("wavelength:   {:.3f} Angstrem".format(wavelength))
        if gamma_0 is not None:
            ls_out.append("offset gamma: {:.3f} deg.".format(gamma_0))
        if nu_0 is not None:
            ls_out.append("offset nu:    {:.3f} deg.".format(nu_0))

        ls_out.extend(["\nOrientation matrix U:\n{:12.5f}{:12.5f}{:12.5f}\n{:12.5f}{:12.5f}{:12.5f}\n{:12.5f}{:12.5f}{:12.5f}".format(
                        u[0, 0], u[0, 1], u[0, 2], 
                        u[1, 0], u[1, 1], u[1, 2], 
                        u[2, 0], u[2, 1], u[2, 2])])
        ls_out.extend(["\nThe UB matrix (CCSL):\n{:12.5f}{:12.5f}{:12.5f}\n{:12.5f}{:12.5f}{:12.5f}\n{:12.5f}{:12.5f}{:12.5f}".format(
                        ub_ccsl[0, 0], ub_ccsl[0, 1], ub_ccsl[0, 2], 
                        ub_ccsl[1, 0], ub_ccsl[1, 1], ub_ccsl[1, 2], 
                        ub_ccsl[2, 0], ub_ccsl[2, 1], ub_ccsl[2, 2])])
        ls_out.extend(["\nThe UB matrix (from):\n{:12.5f}{:12.5f}{:12.5f}\n{:12.5f}{:12.5f}{:12.5f}\n{:12.5f}{:12.5f}{:12.5f}".format(
                        ub_from[0, 0], ub_from[0, 1], ub_from[0, 2], 
                        ub_from[1, 0], ub_from[1, 1], ub_from[1, 2], 
                        ub_from[2, 0], ub_from[2, 1], ub_from[2, 2])])
        ls_out.append(str(cell))
        return "\n".join(ls_out)

    def _show_message(self, s_out: str):
        print("***  Error ***")
        print(s_out)

    @property
    def to_cif(self):
        ls_out = []
        ls_out.append("_diffrn_orient_matrix_UB_11 {:}".format(self.ub_11))
        ls_out.append("_diffrn_orient_matrix_UB_12 {:}".format(self.ub_12))
        ls_out.append("_diffrn_orient_matrix_UB_13 {:}".format(self.ub_13))
        ls_out.append("_diffrn_orient_matrix_UB_21 {:}".format(self.ub_21))
        ls_out.append("_diffrn_orient_matrix_UB_22 {:}".format(self.ub_22))
        ls_out.append("_diffrn_orient_matrix_UB_23 {:}".format(self.ub_23))
        ls_out.append("_diffrn_orient_matrix_UB_31 {:}".format(self.ub_31))
        ls_out.append("_diffrn_orient_matrix_UB_32 {:}".format(self.ub_32))
        ls_out.append("_diffrn_orient_matrix_UB_33 {:}".format(self.ub_33))
        return "\n".join(ls_out)


    def from_cif(self, string: str):
        cif_global = Global()
        flag = cif_global.take_from_string(string)
        if not flag:
            return False
        flag = False
        ub_ccsl = numpy.zeros((3, 3), dtype=float)
        if cif_global.is_value("_diffrn_orient_matrix_UB_11"):
            ub_ccsl[0, 0] = float(cif_global["_diffrn_orient_matrix_UB_11"]) 
        if cif_global.is_value("_diffrn_orient_matrix_UB_12"):
            ub_ccsl[0, 1] = float(cif_global["_diffrn_orient_matrix_UB_12"]) 
        if cif_global.is_value("_diffrn_orient_matrix_UB_13"):
            ub_ccsl[0, 2] = float(cif_global["_diffrn_orient_matrix_UB_13"])

        if cif_global.is_value("_diffrn_orient_matrix_UB_21"):
            ub_ccsl[1, 0] = float(cif_global["_diffrn_orient_matrix_UB_21"]) 
        if cif_global.is_value("_diffrn_orient_matrix_UB_22"):
            ub_ccsl[1, 1] = float(cif_global["_diffrn_orient_matrix_UB_22"]) 
        if cif_global.is_value("_diffrn_orient_matrix_UB_23"):
            ub_ccsl[1, 2] = float(cif_global["_diffrn_orient_matrix_UB_23"])

        if cif_global.is_value("_diffrn_orient_matrix_UB_31"):
            ub_ccsl[2, 0] = float(cif_global["_diffrn_orient_matrix_UB_31"]) 
        if cif_global.is_value("_diffrn_orient_matrix_UB_32"):
            ub_ccsl[2, 1] = float(cif_global["_diffrn_orient_matrix_UB_32"]) 
        if cif_global.is_value("_diffrn_orient_matrix_UB_33"):
            ub_ccsl[2, 2] = float(cif_global["_diffrn_orient_matrix_UB_33"])
        self.ub_ccsl = ub_ccsl
        return True


    def read_ubfrom(self, f_name):
        fid = open(f_name)
        l_cont = fid.readlines()
        fid.close()
        wavelength = float(l_cont[0])
        ub_from = numpy.array([[float(hh_2) for hh_2 in hh_1.strip().split()] for hh_1 in l_cont[1:4]], dtype=float)
        gamma_0, nu_0 = [float(hh) for hh in l_cont[4].strip().split()]
        self.__gamma_0 = gamma_0
        self.__nu_0 = nu_0
        self.__wavelength = wavelength
        
        self.ub_from = ub_from

    def calc_angle(self, h,k,l):
        ub_from = self.ub_from
        q_ub = numpy.dot(ub_from, [h, k, l])
        q = numpy.linalg.norm(q_ub)
        wavelength = self.wavelength
        gamma_0 = self.gamma_0
        nu_0 = self.nu_0

        q_final = [- q**2 * wavelength / 2.0, 0, q_ub[2]]
        q_y2 = q**2 - q_final[0]**2 - q_final[2]**2
        flag_hh_2 = True

        if q_y2<0:
            print("Angles are not found.")
            flag_hh_2 = False
        if flag_hh_2:
            q_final[1] = numpy.sqrt(q_y2)
            k_f = [q_final[0] + 1.0/wavelength, q_final[1], q_final[2]]
            nu = numpy.arcsin(k_f[2]*wavelength)/numpy.pi*180. - nu_0
        
            gamma = numpy.arctan2(k_f[1],k_f[0])/numpy.pi*180. - gamma_0
            phi = (numpy.arctan2(q_final[1], q_final[0]) -numpy.arctan2(q_ub[1], q_ub[0]))/numpy.pi*180
            [phi, gamma, nu] = [phi if phi > 0. else phi + 360., gamma if gamma > 0. else gamma + 360., nu]
            print("gamma is {:7.3f}   nu is {:7.3f}   phi is {:7.3f}".format(gamma, nu, phi))
