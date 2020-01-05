__author__ = 'ikibalin'
__version__ = "2019_12_06"

import os
import numpy
import copy
import warnings

from pycifstar import Global


from typing import List, Tuple
from cryspy.common.cl_item_constr import ItemConstr
from cryspy.common.cl_loop_constr import LoopConstr
from cryspy.common.cl_fitable import Fitable


from cryspy.corecif.cl_cell import Cell


class DiffrnOrientMatrix(ItemConstr):
    """
Data items in the DIFFRN_ORIENT_MATRIX category record details
about the orientation matrix used in the measurement of the
diffraction intensities.

Description in cif file::

 _diffrn_orient_matrix_UB_11           -0.04170
 _diffrn_orient_matrix_UB_12           -0.01429
 _diffrn_orient_matrix_UB_13           -0.02226
 _diffrn_orient_matrix_UB_21           -0.00380
 _diffrn_orient_matrix_UB_22           -0.05578
 _diffrn_orient_matrix_UB_23           -0.05048
 _diffrn_orient_matrix_UB_31            0.00587
 _diffrn_orient_matrix_UB_32           -0.13766
 _diffrn_orient_matrix_UB_33            0.02277

    """
    MANDATORY_ATTRIBUTE = ("ub_11", "ub_12", "ub_13", "ub_21", "ub_22", "ub_23",
                           "ub_31", "ub_32", "ub_33")
    OPTIONAL_ATTRIBUTE = ("id", "type")
    INTERNAL_ATTRIBUTE = ("u_11", "u_12", "u_13", "u_21", "u_22", "u_23",
                          "u_31", "u_32", "u_33", "cell")
    PREFIX = "diffrn_orient_matrix"
    def __init__(self,  ub_11=None, ub_12=None, ub_13=None, ub_21=None, ub_22=None, ub_23=None, ub_31=None, ub_32=None, ub_33=None,
                        id=None):
        super(DiffrnOrientMatrix, self).__init__(mandatory_attribute=self.MANDATORY_ATTRIBUTE, 
                                                 optional_attribute=self.OPTIONAL_ATTRIBUTE, 
                                                 internal_attribute=self.INTERNAL_ATTRIBUTE,
                                                 prefix=self.PREFIX)
        self.ub_11 = ub_11
        self.ub_12 = ub_12
        self.ub_13 = ub_13
        self.ub_21 = ub_21
        self.ub_22 = ub_22
        self.ub_23 = ub_23
        self.ub_31 = ub_31
        self.ub_32 = ub_32
        self.ub_33 = ub_33

        if self.is_defined:
            self.form_object


    @property
    def type(self):
        """
A description of the orientation matrix type and how it should
be applied to define the orientation of the crystal precisely
with respect to the diffractometer axes.

Type: char        
        """
        return getattr(self, "__type")
    @type.setter
    def type(self, x):
        if x is None:
            x_in = None
        else:
            x_in = str(x)
        setattr(self, "__type", x_in)


    @property
    def ub_11(self):
        """
The elements of the diffractometer orientation matrix. These
define the dimensions of the reciprocal cell and its orientation
to the local diffractometer axes.

Type: numb
        """
        return getattr(self, "__ub_11")
    @ub_11.setter
    def ub_11(self, x):
        if x is None:
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__ub_11", x_in)


    @property
    def ub_12(self):
        """
See definition for **_diffrn_orient_matrix.ub_11**
        """
        return getattr(self, "__ub_12")
    @ub_12.setter
    def ub_12(self, x):
        if x is None:
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__ub_12", x_in)

    @property
    def ub_13(self):
        """
See definition for **_diffrn_orient_matrix.ub_11**
        """
        return getattr(self, "__ub_13")
    @ub_13.setter
    def ub_13(self, x):
        if x is None:
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__ub_13", x_in)

    @property
    def ub_21(self):
        """
See definition for **_diffrn_orient_matrix.ub_11**
        """
        return getattr(self, "__ub_21")
    @ub_21.setter
    def ub_21(self, x):
        if x is None:
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__ub_21", x_in)

    @property
    def ub_22(self):
        """
See definition for **_diffrn_orient_matrix.ub_11**
        """
        return getattr(self, "__ub_22")
    @ub_22.setter
    def ub_22(self, x):
        if x is None:
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__ub_22", x_in)

    @property
    def ub_23(self):
        """
See definition for **_diffrn_orient_matrix.ub_11**
        """
        return getattr(self, "__ub_23")
    @ub_23.setter
    def ub_23(self, x):
        if x is None:
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__ub_23", x_in)

    @property
    def ub_31(self):
        """
See definition for **_diffrn_orient_matrix.ub_11**
        """
        return getattr(self, "__ub_31")
    @ub_31.setter
    def ub_31(self, x):
        if x is None:
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__ub_31", x_in)

    @property
    def ub_32(self):
        """
See definition for **_diffrn_orient_matrix.ub_11**
        """
        return getattr(self, "__ub_32")
    @ub_32.setter
    def ub_32(self, x):
        if x is None:
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__ub_32", x_in)

    @property
    def ub_33(self):
        """
See definition for **_diffrn_orient_matrix.ub_11**
        """
        return getattr(self, "__ub_33")
    @ub_33.setter
    def ub_33(self, x):
        if x is None:
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__ub_33", x_in)

    @property
    def u_11(self):
        return getattr(self, "__u_11")
    @property
    def u_12(self):
        return getattr(self, "__u_12")
    @property
    def u_13(self):
        return getattr(self, "__u_13")
    @property
    def u_21(self):
        return getattr(self, "__u_21")
    @property
    def u_22(self):
        return getattr(self, "__u_22")
    @property
    def u_23(self):
        return getattr(self, "__u_23")
    @property
    def u_31(self):
        return getattr(self, "__u_31")
    @property
    def u_32(self):
        return getattr(self, "__u_32")
    @property
    def u_33(self):
        return getattr(self, "__u_33")
    @property
    def cell(self):
        return getattr(self, "__cell")



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
            self.ub_11 = float(x_in[0,0])
            self.ub_12 = float(x_in[0,1])
            self.ub_13 = float(x_in[0,2])
            self.ub_21 = float(x_in[1,0])
            self.ub_22 = float(x_in[1,1])
            self.ub_23 = float(x_in[1,2])
            self.ub_31 = float(x_in[2,0])
            self.ub_32 = float(x_in[2,1])
            self.ub_33 = float(x_in[2,2])
            self.__recalc_u_cell
        else:
                self._show_message("The UB matrix should be numpy.ndarray with size(3x3)")
                return

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
                warnings.warn("The U matrix should be numpy.ndarray with size(3x3)", UserWarning, stacklevel=2)
                return
        flag = x_in.shape == (3,3)
        if flag:
            setattr(self, "__u_11", float(x_in[0,0]))
            setattr(self, "__u_12", float(x_in[0,1]))
            setattr(self, "__u_13", float(x_in[0,2]))
            setattr(self, "__u_21", float(x_in[1,0]))
            setattr(self, "__u_22", float(x_in[1,1]))
            setattr(self, "__u_23", float(x_in[1,2]))
            setattr(self, "__u_31", float(x_in[2,0]))
            setattr(self, "__u_32", float(x_in[2,1]))
            setattr(self, "__u_33", float(x_in[2,2]))
            self.__recalc_ub
        else:
                warnings.warn("The U matrix should be numpy.ndarray with size(3x3)", UserWarning, stacklevel=2)
                return



    def __from_ub_from_to_ub_ccsl(self, ub_from):
        ub_ccsl = numpy.array([[ub_from[1, 0], ub_from[1, 1], ub_from[1, 2]],
                               [-1.*ub_from[0, 0], -1.*ub_from[0, 1], -1.*ub_from[0, 2]],
                               [ub_from[2, 0], ub_from[2, 1], ub_from[2, 2]]], dtype=float)
        self.ub_ccsl = ub_ccsl

    @property
    def form_object(self)->bool:
        flag = False
        if self.is_defined:
            self.__recalc_u_cell
            flag = True
        return flag

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

        if self.cell is None:
            cell = Cell(length_a=a_1, length_b=a_2, length_c=a_3,
                        angle_alpha=float(alpha*180./numpy.pi),
                        angle_beta=float(beta*180./numpy.pi),
                        angle_gamma=float(gamma*180./numpy.pi))
            setattr(self, "__cell", cell)
        else:    
            cell = self.cell
            cell.length_a, cell.length_b, cell.length_c = a_1, a_2, a_3
            cell.angle_alpha, cell.angle_beta, cell.angle_gamma = float(alpha*180./numpy.pi), float(beta*180./numpy.pi), float(gamma*180./numpy.pi),
        cell.form_object
        m_ib = cell.m_ib
        u = numpy.dot(ub_from,  m_ib)

        setattr(self, "__u_11", float(u[0,0]))
        setattr(self, "__u_12", float(u[0,1]))
        setattr(self, "__u_13", float(u[0,2]))
        setattr(self, "__u_21", float(u[1,0]))
        setattr(self, "__u_22", float(u[1,1]))
        setattr(self, "__u_23", float(u[1,2]))
        setattr(self, "__u_31", float(u[2,0]))
        setattr(self, "__u_32", float(u[2,1]))
        setattr(self, "__u_33", float(u[2,2]))

    @property
    def __recalc_ub(self):
        u = self.u
        cell = self.cell
        m_b = cell.m_b
        ub_from = numpy.dot(u, m_b)
        ub_ccsl = numpy.array([[ub_from[1, 0], ub_from[1, 1], ub_from[1, 2]],
                       [-1.*ub_from[0, 0], -1.*ub_from[0, 1], -1.*ub_from[0, 2]],
                       [ub_from[2, 0], ub_from[2, 1], ub_from[2, 2]]], dtype=float)

        self.ub_11 = ub_ccsl[0, 0]
        self.ub_12 = ub_ccsl[0, 1]
        self.ub_13 = ub_ccsl[0, 2]
        self.ub_21 = ub_ccsl[1, 0]
        self.ub_22 = ub_ccsl[1, 1]
        self.ub_23 = ub_ccsl[1, 2]
        self.ub_31 = ub_ccsl[2, 0]
        self.ub_32 = ub_ccsl[2, 1]
        self.ub_33 = ub_ccsl[2, 2]


    def __repr__(self):
        ls_out = ["OrientMatrix:"]
        ls_out.append(str(self))
        return "\n".join(ls_out)

    def _show_message(self, s_out: str):
        warnings.warn("***  Error ***\n"+s_out, UserWarning, stacklevel=2)


    #def read_ubfrom(self, f_name="ubfrom.raf"):
    #    fid = open(f_name)
    #    l_cont = fid.readlines()
    #    fid.close()
    #    wavelength = float(l_cont[0])
    #    ub_from = numpy.array([[float(hh_2) for hh_2 in hh_1.strip().split()] for hh_1 in l_cont[1:4]], dtype=float)
    #    gamma_0, nu_0 = [float(hh) for hh in l_cont[4].strip().split()]
    #    self.__gamma_0 = gamma_0
    #    self.__nu_0 = nu_0
    #    self.__wavelength = wavelength
    #    
    #    self.ub_from = ub_from

    def calc_angle(self, h, k, l, wavelength=1.4):
        ub_from = self.ub_from
        q_ub = numpy.dot(ub_from, [h, k, l])
        q = numpy.linalg.norm(q_ub)
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
