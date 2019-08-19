"""
define classes to describe crystal 
"""
__author__ = 'ikibalin'
__version__ = "2019_04_06"
import os
import numpy

from neupy.f_common.cl_variable import Variable
from neupy.f_interface.cl_abstract_cell import AbstractCell

class Cell(AbstractCell):
    """
    Cell parameters
    """
    def __init__(self, a = 1.0, b = 1.0, c = 1.0, alpha = 90.0, beta = 90.0, 
                 gamma= 90., singony = "Triclinic"):
        super(Cell, self).__init__()
        self._p_a = None
        self._p_b = None
        self._p_c = None
        self._p_alpha = None
        self._p_beta = None
        self._p_gamma = None
        self._p_singony = None
        
        self._p_cos_a = None
        self._p_cos_b = None
        self._p_cos_g = None
        self._p_cos_a_sq = None
        self._p_cos_b_sq = None
        self._p_cos_g_sq = None
        self._p_sin_a = None
        self._p_sin_b = None
        self._p_sin_g = None
        self._p_sin_a_sq = None
        self._p_sin_b_sq = None
        self._p_sin_g_sq = None
        
        self._p_ia = None
        self._p_ib = None
        self._p_ic = None
        self._p_ialpha = None
        self._p_ibeta = None
        self._p_igamma = None        

        self._p_cos_ia = None
        self._p_cos_ib = None
        self._p_cos_ig = None
        self._p_cos_ia_sq = None
        self._p_cos_ib_sq = None
        self._p_cos_ig_sq = None
        self._p_sin_ia = None
        self._p_sin_ib = None
        self._p_sin_ig = None
        self._p_sin_ia_sq = None
        self._p_sin_ib_sq = None
        self._p_sin_ig_sq = None
        
        self._p_vol = None
        self._p_ivol = None
        self._p_m_b = None
        self._p_m_ib = None

        self._refresh(a, b, c, alpha, beta, gamma, singony)
        self.set_val()
        
    def __repr__(self):
        lsout = """Cell: \n a: {:}\n b: {:}\n c: {:}\n alpha: {:}
 beta: {:}\n gamma: {:}\n Bravais lattice: {:}\n""".format(self._p_a, self._p_b, 
                 self._p_c, self._p_alpha, self._p_beta, self._p_gamma, 
                 self._p_singony)
        if self._p_m_b is not None:
             lsout += """ B matrix is:\n {:9.5f} {:9.5f} {:9.5f}
 {:9.5f} {:9.5f} {:9.5f}\n {:9.5f} {:9.5f} {:9.5f}\n""".format(self._p_m_b[0, 0],
 self._p_m_b[0, 1], self._p_m_b[0, 2], self._p_m_b[1, 0], self._p_m_b[1, 1], 
 self._p_m_b[1, 2], self._p_m_b[2, 0], self._p_m_b[2, 1], self._p_m_b[2, 2])
        if self._p_m_ib is not None:
             lsout += """ inversed B matrix is:\n {:9.5f} {:9.5f} {:9.5f}
 {:9.5f} {:9.5f} {:9.5f}\n {:9.5f} {:9.5f} {:9.5f}\n""".format(self._p_m_ib[0, 0],
 self._p_m_ib[0, 1], self._p_m_ib[0, 2], self._p_m_ib[1, 0], self._p_m_ib[1, 1], 
 self._p_m_ib[1, 2], self._p_m_ib[2, 0], self._p_m_ib[2, 1], self._p_m_ib[2, 2])
        return lsout
    
    def _refresh(self, a, b, c, alpha, beta, gamma, singony):
        """
        refresh variables
        """
        if a is not None:
            self._p_a = a
        if b is not None:
            self._p_b = b
        if c is not None:
            self._p_c = c
        if alpha is not None:
            self._p_alpha = alpha
        if beta is not None:
            self._p_beta = beta
        if gamma is not None:
            self._p_gamma = gamma
        if singony is not None:
            self._p_singony = singony

        cond = any([hh is not None for hh in [a, b, c, alpha, beta, gamma, 
                                          singony]])
        if cond:
            self.apply_constraint()
    
    def _constr_singony(self):
        singony = self._p_singony
        if singony == "Cubic":
            self._p_b = 1.*self._p_a
            self._p_c = 1.*self._p_a
            self._p_alpha = 90.
            self._p_beta = 90.
            self._p_gamma = 90.
        elif singony == "Hexagonal":
            self._p_b = 1.*self._p_a
            self._p_alpha = 90.
            self._p_beta = 90.
            self._p_gamma = 120.        
        elif singony == "Rhombohedral":
            self._p_b = 1.*self._p_a
            self._p_c = 1.*self._p_a
            self._p_beta = 1.*self._p_alpha
            self._p_gamma = 1.*self._p_alpha
        elif singony == "Trigonal":
            self._p_b = 1.*self._p_a
            self._p_c = 1.*self._p_a
        elif singony == "Tetragonal":
            self._p_b = 1.*self._p_a
            self._p_alpha = 90.
            self._p_beta = 90.
            self._p_gamma = 90.
        elif singony == "Orthorhombic":
            self._p_alpha = 90.
            self._p_beta = 90.
            self._p_gamma = 90.
        elif singony == "Monoclinic":
            self._p_alpha = 90.
            self._p_gamma = 90.

    def set_val(self, a = None, b = None, c = None, alpha = None, 
                   beta = None, gamma= None, singony = None):
        self._refresh(a, b, c, alpha, beta, gamma, singony)


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
a, b, c - unit cell parameters in Angstrems
alpha, beta, gamma - angles in degrees

singony - singony: "Triclinic", "Monoclinic", "Orthorhombic", "Tetragonal", 
                   "Trigonal", "Hexagonal" or "Cubic"

ia, ib, ic - inverse unit cell parameters in Angstrems**-1
ialpha, ibeta, igamma - angles of inverse cell in degrees

vol - volume of unit cell in Angstrems**3
ivol - volume of inverse cell in Angstrems**3

m_b - matrix B (x is along ia, y is in (ia,ib) plane, z is vector product x, y)
m_ib - inverse B matrix 
        """
        print(lsout)
            
    def _calc_cos_abc(self):
        rad=numpy.pi/180.
        self._p_cos_a = numpy.cos(self._p_alpha*rad)
        self._p_cos_b = numpy.cos(self._p_beta*rad)
        self._p_cos_g = numpy.cos(self._p_gamma*rad)
        
        self._p_sin_a = numpy.sin(self._p_alpha*rad)
        self._p_sin_b = numpy.sin(self._p_beta*rad)
        self._p_sin_g = numpy.sin(self._p_gamma*rad)
        
        self._p_cos_a_sq = self._p_cos_a**2
        self._p_cos_b_sq = self._p_cos_b**2
        self._p_cos_g_sq = self._p_cos_g**2

        self._p_sin_a_sq = 1.-self._p_cos_a_sq
        self._p_sin_b_sq = 1.-self._p_cos_b_sq
        self._p_sin_g_sq = 1.-self._p_cos_g_sq
        
    def _calc_cos_iabc(self):
        rad=numpy.pi/180.
        self._p_cos_ia = numpy.cos(self._p_ialpha*rad)
        self._p_cos_ib = numpy.cos(self._p_ibeta*rad)
        self._p_cos_ig = numpy.cos(self._p_igamma*rad)
        
        self._p_sin_ia = numpy.sin(self._p_ialpha*rad)
        self._p_sin_ib = numpy.sin(self._p_ibeta*rad)
        self._p_sin_ig = numpy.sin(self._p_igamma*rad)
        
        self._p_cos_ia_sq = self._p_cos_ia**2
        self._p_cos_ib_sq = self._p_cos_ib**2
        self._p_cos_ig_sq = self._p_cos_ig**2

        self._p_sin_a_sq = 1.-self._p_cos_a_sq
        self._p_sin_b_sq = 1.-self._p_cos_b_sq
        self._p_sin_g_sq = 1.-self._p_cos_g_sq

    def _calc_volume(self):
        a = 1.*self._p_a
        b = 1.*self._p_b
        c = 1.*self._p_c
        c_a = self._p_cos_a
        c_b = self._p_cos_b
        c_g = self._p_cos_g
        c_a_sq = self._p_cos_a_sq
        c_b_sq = self._p_cos_b_sq
        c_g_sq = self._p_cos_g_sq
        vol = a*b*c*(1.-c_a_sq-c_b_sq-c_g_sq+2.*c_a*c_b*c_g)**0.5
        self._p_vol = vol
        
    
    def _calc_iucp(self):
        """
        calculate inverse unit cell
        """
        irad = 180./numpy.pi

        a = 1.*self._p_a
        b = 1.*self._p_b
        c = 1.*self._p_c
        c_a = self._p_cos_a
        c_b = self._p_cos_b
        c_g = self._p_cos_g
        s_a = self._p_sin_a
        s_b = self._p_sin_b
        s_g = self._p_sin_g
        vol = self._p_vol
        
        self._p_ialpha = numpy.arccos((c_b*c_g-c_a)/(s_b*s_g))*irad
        self._p_ibeta = numpy.arccos((c_g*c_a-c_b)/(s_g*s_a))*irad
        self._p_igamma = numpy.arccos((c_a*c_b-c_g)/(s_a*s_b))*irad

        self._p_ia = b*c*s_a/vol
        self._p_ib = c*a*s_b/vol
        self._p_ic = a*b*s_g/vol


    def _calc_m_b(self):
        """
        calculate matrix B 
        """
        c = 1.*self._p_c

        ia = self._p_ia 
        ib = self._p_ib 
        ic = self._p_ic 
        
        c_a = self._p_cos_a
        
        #ic_a = self._p_cos_ia 
        ic_b = self._p_cos_ib 
        ic_g = self._p_cos_ig 
        #is_a = self._p_sin_ia 
        is_b = self._p_sin_ib 
        is_g = self._p_sin_ig 
        
        self._p_m_b = numpy.array([[ia,  ib*ic_g,  ic*ic_b],
            [0.,  ib*is_g, -ic*is_b*c_a],
            [0.,       0.,  1./c]], dtype = float)

    def _calc_m_ib(self):
        """
        calculate inverse B matrix 
        """
        x1 = self._p_m_b[0,0]
        x2 = self._p_m_b[1,1]
        x3 = self._p_m_b[2,2]
        x4 = self._p_m_b[0,1]
        x5 = self._p_m_b[0,2]
        x6 = self._p_m_b[1,2]
        #B=[[x1,x4,x5],
        #   [0.,x2,x6],
        #   [0.,0.,x3]]
        #it shuld be checked
        #iB=numpy.linalg.inv(B)
        y1 = 1./x1
        y2 = 1./x2
        y3 = 1./x3
        y4 = -1*x4*1./(x1*x2)
        y6 = -1*x6*1./(x2*x3)
        y5 = (x4*x6-x2*x5)*1./(x1*x2*x3)
        
        self._p_m_ib = numpy.array([[y1,y4,y5],[0.,y2,y6],[0.,0.,y3]], 
                                   dtype = float)
            
                
        
    
    def calc_sthovl(self, h, k, l):
        """
        calculate sin(theta)/lambda for list of hkl reflections
        """
            
        a = 1.*self._p_a
        b = 1.*self._p_b
        c = 1.*self._p_c
        c_a = self._p_cos_a
        c_b = self._p_cos_b
        c_g = self._p_cos_g
        c_a_sq = self._p_cos_a_sq
        c_b_sq = self._p_cos_b_sq
        c_g_sq = self._p_cos_g_sq
        s_a_sq = self._p_sin_a_sq
        s_b_sq = self._p_sin_b_sq
        s_g_sq = self._p_sin_g_sq

        A=( 1. - c_a_sq - c_b_sq - c_g_sq + 2.*c_a*c_b*c_g)
        B1 = (s_a_sq*(h*1./a)**2+s_b_sq*(k*1./b)**2+s_g_sq*(l*1./c)**2)
        B2 = 2.*(k*l*c_a)/(b*c)+2.*(h*l*c_b)/(a*c)+2.*(h*k*c_g)/(a*b)
        #it should be checked, I am not sure
        B = B1-B2
        inv_d = (B*1./A)**0.5
        return 0.5*inv_d

    def calc_k_loc(self, h, k, l):
        """
        calculate unity scattering vector
        """
        m_b = self.get_val("m_b")
        k_x = m_b[0, 0]*h + m_b[0, 1]*k +m_b[0, 2]*l
        k_y = m_b[1, 0]*h + m_b[1, 1]*k +m_b[1, 2]*l
        k_z = m_b[2, 0]*h + m_b[2, 1]*k +m_b[2, 2]*l
        
        k_norm = (k_x**2 + k_y**2 + k_z**2)**0.5
        if not((type(h) is float)|(type(h) is int)|(type(h) is numpy.float64)):
            k_norm[k_norm == 0.] = 1.
        elif k_norm == 0.:
            k_norm = 1.
        
        k_x = k_x/k_norm
        k_y = k_y/k_norm
        k_z = k_z/k_norm
        
        return k_x, k_y, k_z
        
    def calc_m_t(self, h, k, l):
        """define rotation matrix to have new z axis along kloc
        Rotation matrix is defined by Euler angles
        """
        m_b = self.get_val("m_b")
        k_x = m_b[0, 0]*h + m_b[0, 1]*k +m_b[0, 2]*l
        k_y = m_b[1, 0]*h + m_b[1, 1]*k +m_b[1, 2]*l
        k_z = m_b[2, 0]*h + m_b[2, 1]*k +m_b[2, 2]*l
        
        k_norm = (k_x**2 + k_y**2 + k_z**2)**0.5
        k_norm[k_norm == 0.] = 1.
        
        k_x = k_x/k_norm
        k_y = k_y/k_norm
        k_z = k_z/k_norm
        

        al = numpy.zeros(k_x.shape, dtype=float)
        
        be = numpy.arccos(k_z)
        sb = numpy.sin(be)
        flag = (sb != 0.)
        
        sa1 = k_x[flag]*1./sb[flag]
        ca2 = -1*k_y[flag]*1./sb[flag]
        sa1[sa1>1] = 1.
        sa1[sa1<-1] = -1.
            
        ca2[ca2>1] = 1.
        ca2[ca2<-1] = -1.

        al1 = numpy.arcsin(sa1)
        al2 = numpy.arccos(ca2)
        
        al_sh = numpy.copy(al1)
        al_sh[sa1 > 0.] = al2[sa1 > 0.]
        al_sh[sa1 <= 0.] = 2.*numpy.pi-al2[sa1 <= 0.]
        al_sh[numpy.abs(al2-al1)<0.00001] = al1[numpy.abs(al2-al1)<0.00001]

        al[flag] = al_sh
            
        ga=0.
        ca, cb, cg = numpy.cos(al), numpy.cos(be), numpy.cos(ga)
        sa, sb, sg = numpy.sin(al), numpy.sin(be), numpy.sin(ga)
        t_11, t_12, t_13 = ca*cg-sa*cb*sg, -ca*sg-sa*cb*cg,  sa*sb
        t_21, t_22, t_23 = sa*cg+ca*cb*sg, -sa*sg+ca*cb*cg, -ca*sb
        t_31, t_32, t_33 =          sb*sg,           sb*cg,     cb
        
        flag = (((sa*sb-k_x)**2+(-ca*sb-k_y)**2+(cb-k_z)**2)>0.0001)
        if any(flag):
            print("Mistake with k_loc")
            print("Program is stopped")
            quit()
        return t_11, t_12, t_13, t_21, t_22, t_23, t_31, t_32, t_33 
    
    def is_variable(self):
        """
        without extinction
        """
        res = any([isinstance(self._p_a, Variable), 
                   isinstance(self._p_b, Variable),
                   isinstance(self._p_c, Variable),
                   isinstance(self._p_alpha, Variable),
                   isinstance(self._p_beta, Variable),
                   isinstance(self._p_gamma, Variable)])
        return res
        
    def get_variables(self):
        l_variable = []
        if isinstance(self._p_a, Variable):
            l_variable.append(self._p_a)
        if isinstance(self._p_b, Variable):
            l_variable.append(self._p_b)
        if isinstance(self._p_c, Variable):
            l_variable.append(self._p_c)
        if isinstance(self._p_alpha, Variable):
            l_variable.append(self._p_alpha)
        if isinstance(self._p_beta, Variable):
            l_variable.append(self._p_beta)
        if isinstance(self._p_gamma, Variable):
            l_variable.append(self._p_gamma)
        return l_variable
    
    def apply_constraint(self):
        self._constr_singony()
        self._calc_cos_abc()
        self._calc_volume()
        self._calc_iucp()
        self._calc_cos_iabc()
        self._calc_m_b()
        self._calc_m_ib()

    def calc_hkl(self, space_group, sthovl_min, sthovl_max):
        """
        give a list of reflections hkl for cell in the range sthovl_min, sthovl_max
        

        """
        lhkl,lmult=[],[]
        lhklres=[]

        hmax = int(2.*self.get_val('a')*sthovl_max)
        kmax = int(2.*self.get_val('b')*sthovl_max)
        lmax = int(2.*self.get_val('c')*sthovl_max)
        hmin, kmin, lmin = -1*hmax, -1*kmax, -1*lmax

        hmin=0
        
        lorig = space_group.get_val("orig")
        lsymm = space_group.get_val("el_symm")
        for h in range(hmin,hmax+1,1):
            for k in range(kmin,kmax+1,1):
                for l in range(lmin,lmax+1,1):
                    flag=(abs(sum([numpy.exp(2.*numpy.pi*1j*(orig[0]*h+orig[1]*k+orig[2]*l)) for orig in lorig]))>0.00001)
                    #flag=True
                    if (flag):
                        lhkls=[(h*symm[1]+k*symm[5]+l*symm[9], h*symm[2]+k*symm[6]+l*symm[10], h*symm[3]+k*symm[7]+l*symm[11]) for symm in lsymm]
                        lhkls.extend([(-hkl[0],-hkl[1],-hkl[2]) for hkl in lhkls])
                        lhkls.sort(key=lambda x:10000*x[0]+100*x[1]+x[2])
                        if (not(lhkls[-1] in lhkl)):
                            lhkl.append(lhkls[-1])
                            lmult.append(len(set(lhkls)))
                            
        lhklsthovl=[(hkl, self.calc_sthovl(hkl[0], hkl[1], hkl[2]), mult) for hkl, mult in zip(lhkl, lmult)]
        lhklsthovl.sort(key=lambda x: x[1])
        lhklres = [hklsthovl[0] for hklsthovl in lhklsthovl if ((hklsthovl[1]>sthovl_min) & (hklsthovl[1]<sthovl_max))]
        lmultres = [hklsthovl[2] for hklsthovl in lhklsthovl if ((hklsthovl[1]>sthovl_min) & (hklsthovl[1]<sthovl_max))]

        h = numpy.array([hh[0] for hh in lhklres], dtype=int)
        k = numpy.array([hh[1] for hh in lhklres], dtype=int)
        l = numpy.array([hh[2] for hh in lhklres], dtype=int)
        mult = numpy.array(lmultres, dtype=int)
        return h, k, l, mult


    def calc_hkl_in_range(self, sthovl_min, sthovl_max):
        h_max = int(2.*self.get_val('a')*sthovl_max)
        k_max = int(2.*self.get_val('b')*sthovl_max)
        l_max = int(2.*self.get_val('c')*sthovl_max)
        h_min, k_min, l_min = -1*h_max, -1*k_max, -1*l_max

        np_h = numpy.array(range(h_min, h_max+1, 1), dtype=int)
        np_k = numpy.array(range(k_min, k_max+1, 1), dtype=int)
        np_l = numpy.array(range(l_min, l_max+1, 1), dtype=int)
        h_3d, k_3d, l_3d = numpy.meshgrid(np_h, np_k, np_l, indexing="ij")
        
        sthovl_3d = self.calc_sthovl(h_3d, k_3d, l_3d)
        flag_1 = sthovl_3d >= sthovl_min
        flag_2 = sthovl_3d <= sthovl_max
        flag_12 = numpy.logical_and(flag_1, flag_2)
        
        h = h_3d[flag_12]
        k = k_3d[flag_12]
        l = l_3d[flag_12]
        mult = numpy.ones(h.size, dtype=int)
        sthovl = sthovl_3d[flag_12]
        arg_sort = numpy.argsort(sthovl)
        return h[arg_sort], k[arg_sort], l[arg_sort], mult[arg_sort] 

        
if (__name__ == "__main__"):
  pass

