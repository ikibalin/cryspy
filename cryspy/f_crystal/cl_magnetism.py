"""
internal class to calculate magnetic form factor
"""
__author__ = 'ikibalin'
__version__ = "2019_09_02"
import os
import numpy


def calc_mRmCmRT(r11, r12, r13, r21, r22, r23, r31, r32, r33,
                 c11, c12, c13, c21, c22, c23, c31, c32, c33):
    """
    calculate matrix multiplication R*C*RT, when matrices are expressed through 
    its component and can be expressed as nD-array
    """
    rc_11, rc_12 = r11*c11+r12*c21+r13*c31, r11*c12+r12*c22+r13*c32
    rc_13 = r11*c13+r12*c23+r13*c33
    rc_21, rc_22 = r21*c11+r22*c21+r23*c31, r21*c12+r22*c22+r23*c32
    rc_23 = r21*c13+r22*c23+r23*c33
    rc_31, rc_32 = r31*c11+r32*c21+r33*c31, r31*c12+r32*c22+r33*c32
    rc_33 = r31*c13+r32*c23+r33*c33

    #dimension (atoms, symmetry)
    rcrt_11 = (rc_11*r11+rc_12*r12+rc_13*r13)
    rcrt_12 = (rc_11*r21+rc_12*r22+rc_13*r23)
    rcrt_13 = (rc_11*r31+rc_12*r32+rc_13*r33)

    rcrt_21 = (rc_21*r11+rc_22*r12+rc_23*r13)
    rcrt_22 = (rc_21*r21+rc_22*r22+rc_23*r23)
    rcrt_23 = (rc_21*r31+rc_22*r32+rc_23*r33)

    rcrt_31 = (rc_31*r11+rc_32*r12+rc_33*r13)
    rcrt_32 = (rc_31*r21+rc_32*r22+rc_33*r23)
    rcrt_33 = (rc_31*r31+rc_32*r32+rc_33*r33)
    return rcrt_11, rcrt_12, rcrt_13, rcrt_21, rcrt_22, rcrt_23, rcrt_31, rcrt_32, rcrt_33



class Magnetism(object):
    """
    Magnetism
    """
    def __init__(self, factor_lande=numpy.array([], dtype=float),
                       kappa=numpy.array([], dtype=float),
                       chi_11=numpy.array([], dtype=float),
                       chi_12=numpy.array([], dtype=float),
                       chi_13=numpy.array([], dtype=float),
                       chi_22=numpy.array([], dtype=float),
                       chi_23=numpy.array([], dtype=float),
                       chi_33=numpy.array([], dtype=float),
                       j0_A=numpy.array([], dtype=float),
                       j0_a=numpy.array([], dtype=float),
                       j0_B=numpy.array([], dtype=float),
                       j0_b=numpy.array([], dtype=float),
                       j0_C=numpy.array([], dtype=float),
                       j0_c=numpy.array([], dtype=float),
                       j0_D=numpy.array([], dtype=float),
                       j2_A=numpy.array([], dtype=float),
                       j2_a=numpy.array([], dtype=float),
                       j2_B=numpy.array([], dtype=float),
                       j2_b=numpy.array([], dtype=float),
                       j2_C=numpy.array([], dtype=float),
                       j2_c=numpy.array([], dtype=float),
                       j2_D=numpy.array([], dtype=float)):
        super(Magnetism, self).__init__()
        self.__atom_site_aniso_magnetism_chi_11 = None
        self.__atom_site_aniso_magnetism_chi_12 = None
        self.__atom_site_aniso_magnetism_chi_13 = None
        self.__atom_site_aniso_magnetism_chi_22 = None
        self.__atom_site_aniso_magnetism_chi_23 = None
        self.__atom_site_aniso_magnetism_chi_33 = None
        self.__atom_site_magnetism_kappa = None
        self.__atom_site_magnetism_factor_lande = None
        self.__j0_A = None
        self.__j0_a = None
        self.__j0_B = None
        self.__j0_b = None
        self.__j0_C = None
        self.__j0_c = None
        self.__j0_D = None
        self.__j2_A = None
        self.__j2_a = None
        self.__j2_B = None
        self.__j2_b = None
        self.__j2_C = None
        self.__j2_c = None
        self.__j2_D = None

        self.__matrix_chi_loc = None

        self.factor_lande = factor_lande
        self.kappa = kappa
        self.chi_11 = chi_11
        self.chi_12 = chi_12
        self.chi_13 = chi_13
        self.chi_22 = chi_22
        self.chi_23 = chi_23
        self.chi_33 = chi_33
        self.j0_A = j0_A
        self.j0_a = j0_a
        self.j0_B = j0_B
        self.j0_b = j0_b
        self.j0_C = j0_C
        self.j0_c = j0_c
        self.j0_D = j0_D
        self.j2_A = j2_A
        self.j2_a = j2_a
        self.j2_B = j2_B
        self.j2_b = j2_b
        self.j2_C = j2_C
        self.j2_c = j2_c
        self.j2_D = j2_D

    def _trans_to_float_array(self, x):
        if isinstance(x, numpy.ndarray):
            x_out = x.astype(float)
        else:
            x_out = numpy.array([x], dtype=float)
        return x_out

    @property
    def factor_lande(self):
        return self.__atom_site_magnetism_factor_lande
    @factor_lande.setter
    def factor_lande(self, x):
        self.__atom_site_magnetism_factor_lande = self._trans_to_float_array(x)
    @property
    def kappa(self):
        return self.__atom_site_magnetism_kappa
    @kappa.setter
    def kappa(self, x):
        self.__atom_site_magnetism_kappa = self._trans_to_float_array(x)        
    @property
    def chi_11(self):
        return self.__atom_site_aniso_magnetism_chi_11
    @chi_11.setter
    def chi_11(self, x):
        self.__atom_site_aniso_magnetism_chi_11 = self._trans_to_float_array(x)    
    @property
    def chi_22(self):
        return self.__atom_site_aniso_magnetism_chi_22
    @chi_22.setter
    def chi_22(self, x):
        self.__atom_site_aniso_magnetism_chi_22 = self._trans_to_float_array(x)    
    @property
    def chi_33(self):
        return self.__atom_site_aniso_magnetism_chi_33
    @chi_33.setter
    def chi_33(self, x):
        self.__atom_site_aniso_magnetism_chi_33 = self._trans_to_float_array(x)    
    @property
    def chi_12(self):
        return self.__atom_site_aniso_magnetism_chi_12
    @chi_12.setter
    def chi_12(self, x):
        self.__atom_site_aniso_magnetism_chi_12 = self._trans_to_float_array(x)    
    @property
    def chi_13(self):
        return self.__atom_site_aniso_magnetism_chi_13
    @chi_13.setter
    def chi_13(self, x):
        self.__atom_site_aniso_magnetism_chi_13 = self._trans_to_float_array(x)    
    @property
    def chi_23(self):
        return self.__atom_site_aniso_magnetism_chi_23
    @chi_23.setter
    def chi_23(self, x):
        self.__atom_site_aniso_magnetism_chi_23 = self._trans_to_float_array(x)    

    @property
    def j0_A(self):
        return self.__j0_A
    @j0_A.setter
    def j0_A(self, x):
        self.__j0_A = self._trans_to_float_array(x)    
    @property
    def j0_a(self):
        return self.__j0_a
    @j0_a.setter
    def j0_a(self, x):
        self.__j0_a = self._trans_to_float_array(x)    
    @property
    def j0_B(self):
        return self.__j0_B
    @j0_B.setter
    def j0_B(self, x):
        self.__j0_B = self._trans_to_float_array(x)    
    @property
    def j0_b(self):
        return self.__j0_b
    @j0_b.setter
    def j0_b(self, x):
        self.__j0_b = self._trans_to_float_array(x)    
    @property
    def j0_C(self):
        return self.__j0_C
    @j0_C.setter
    def j0_C(self, x):
        self.__j0_C = self._trans_to_float_array(x)    
    @property
    def j0_c(self):
        return self.__j0_c
    @j0_c.setter
    def j0_c(self, x):
        self.__j0_c = self._trans_to_float_array(x)    
    @property
    def j0_D(self):
        return self.__j0_D
    @j0_D.setter
    def j0_D(self, x):
        self.__j0_D = self._trans_to_float_array(x)    


    @property
    def j2_A(self):
        return self.__j2_A
    @j2_A.setter
    def j2_A(self, x):
        self.__j2_A = self._trans_to_float_array(x)    
    @property
    def j2_a(self):
        return self.__j2_a
    @j2_a.setter
    def j2_a(self, x):
        self.__j2_a = self._trans_to_float_array(x)    
    @property
    def j2_B(self):
        return self.__j2_B
    @j2_B.setter
    def j2_B(self, x):
        self.__j2_B = self._trans_to_float_array(x)    
    @property
    def j2_b(self):
        return self.__j2_b
    @j2_b.setter
    def j2_b(self, x):
        self.__j2_b = self._trans_to_float_array(x)    
    @property
    def j2_C(self):
        return self.__j2_C
    @j2_C.setter
    def j2_C(self, x):
        self.__j2_C = self._trans_to_float_array(x)    
    @property
    def j2_c(self):
        return self.__j2_c
    @j2_c.setter
    def j2_c(self, x):
        self.__j2_c = self._trans_to_float_array(x)    
    @property
    def j2_D(self):
        return self.__j2_D
    @j2_D.setter
    def j2_D(self, x):
        self.__j2_D = self._trans_to_float_array(x)    

    @property
    def matrix_chi_loc(self):
        return self.__matrix_chi_loc

    def __repr__(self):
        ls_out = ["Magnetism susceptibility:\n  chi_11   chi_22   chi_33  chi_12   chi_13   chi_23"]
        ls_out.extend(["{:8.5f} {:8.5f} {:8.5f}{:8.5f} {:8.5f} {:8.5f}".format(
            hh_1, hh_2, hh_3, hh_4, hh_5, hh_6) for hh_1, hh_2, hh_3, hh_4, hh_5, hh_6 
            in zip(self.chi_11, self.chi_22, self.chi_33, self.chi_12, self.chi_13, self.chi_23)])
        return "\n".join(ls_out)

    
    def _calc_chi_loc(ia, ib, ic, matrix_ib):
        """
        IT SHOULD BE DELETED


        representation of chi in crystallographic coordinate system defined as x||a*, z||c, y= [z x] (right handed)
        expressions are taken from international tables
        matrix_ib is inversed matrix B
        ia, ib, ic is inversed unit cell parameters (it can be estimated from matrix matrix_ib)

        X = B x, x = iB X
        xT*CHI*x = XT iBT CHI iB X
    
        output chiLOC = iBT CHI iB
        """
        matrix_chi = numpy.array(
                [[self.chi_11, self.chi_12, self.chi_13],
                 [self.chi_12, self.chi_22, self.chi_23],
                 [self.chi_13, self.chi_23, self.chi_33]], 
                 dtype = float)
        #mchi=[[chi[0],chi[3],chi[4]],[chi[3],chi[1],chi[5]],[chi[4],chi[5],chi[2]]]
        #[a,b,c,alpha,beta,gamma]=ucp
        y1 = matrix_ib[0,0]
        y2 = matrix_ib[1,1]
        y3 = matrix_ib[2,2]
        y4 = matrix_ib[0,1]
        y5 = matrix_ib[0,2]
        y6 = matrix_ib[1,2]
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
        matrix_ib_norm = matrix_ib
        #!!!! to check it !!! before it was [:, 0]
        matrix_ib_norm[0, :] *= ia
        matrix_ib_norm[1, :] *= ib
        matrix_ib_norm[2, :] *= ic
        print("There is mistake!!!!!")
        matrix_ibt_norm = matrix_ib_norm.transpose()
        #it is not compatible with case, vhen chi_ij is 1D array 
        ibt_chi = numpy.matmul(matrix_ibt_norm, matrix_chi)
        matrix_chi_loc = numpy.matmul(ibt_chi, matrix_ib_norm)
        self.__matrix_chi_loc = matrix_chi_loc
    
    def calc_form_factor_tensor(self, space_group, cell, h, k, l):
        """
        give components of form factor tensor:
            fft_11, fft_12, fft_13
            fft_21, fft_22, fft_23
            fft_31, fft_32, fft_33
            
        in 3 dimension (hkl, atoms, symmetry elements)            
        """
        sthovl = cell.calc_sthovl(h, k, l)
        #dimension (hkl, atoms)
        form_factor_2d = self._calc_form_factor(sthovl)
        
        r_11, r_12 = space_group.r_11, space_group.r_12
        r_13, r_21 = space_group.r_13, space_group.r_21
        r_22, r_23 = space_group.r_22, space_group.r_23
        r_31, r_32 = space_group.r_31, space_group.r_32
        r_33 = space_group.r_33

        chi_11, chi_22 = self.chi_11, self.chi_22 
        chi_33, chi_12 = self.chi_33, self.chi_12
        chi_13, chi_23 = self.chi_13, self.chi_23
        chi_21, chi_31, chi_32 = chi_12, chi_13, chi_23 
        
        c11, r11 = numpy.meshgrid(chi_11, r_11, indexing="ij")
        c22, r22 = numpy.meshgrid(chi_22, r_22, indexing="ij")
        c33, r33 = numpy.meshgrid(chi_33, r_33, indexing="ij")
        c12, r12 = numpy.meshgrid(chi_12, r_12, indexing="ij")
        c13, r13 = numpy.meshgrid(chi_13, r_13, indexing="ij")
        c23, r23 = numpy.meshgrid(chi_23, r_23, indexing="ij")
        c21, r21 = numpy.meshgrid(chi_21, r_21, indexing="ij")
        c31, r31 = numpy.meshgrid(chi_31, r_31, indexing="ij")
        c32, r32 = numpy.meshgrid(chi_32, r_32, indexing="ij")
        
        
        rcrt_11, rcrt_12, rcrt_13, rcrt_21, rcrt_22, rcrt_23, rcrt_31, rcrt_32, rcrt_33 = calc_mRmCmRT(
                r11, r12, r13, r21, r22, r23, r31, r32, r33,
                c11, c12, c13, c21, c22, c23, c31, c32, c33)        
        
        #dimension (hkl, atoms, symmetry)
        fft_11 = form_factor_2d[:,:,numpy.newaxis]*rcrt_11[numpy.newaxis, :,:]
        fft_12 = form_factor_2d[:,:,numpy.newaxis]*rcrt_12[numpy.newaxis, :,:]
        fft_13 = form_factor_2d[:,:,numpy.newaxis]*rcrt_13[numpy.newaxis, :,:]
        fft_21 = form_factor_2d[:,:,numpy.newaxis]*rcrt_21[numpy.newaxis, :,:]
        fft_22 = form_factor_2d[:,:,numpy.newaxis]*rcrt_22[numpy.newaxis, :,:]
        fft_23 = form_factor_2d[:,:,numpy.newaxis]*rcrt_23[numpy.newaxis, :,:]
        fft_31 = form_factor_2d[:,:,numpy.newaxis]*rcrt_31[numpy.newaxis, :,:]
        fft_32 = form_factor_2d[:,:,numpy.newaxis]*rcrt_32[numpy.newaxis, :,:]
        fft_33 = form_factor_2d[:,:,numpy.newaxis]*rcrt_33[numpy.newaxis, :,:]
        
        #ortogonalization should be done
        
        return fft_11, fft_12, fft_13, fft_21, fft_22, fft_23, fft_31, fft_32, fft_33
    
    def calc_chi_rot(self, matrix_chi, elsymm):
        """
        calculate R*chi*RT
        rotation of chi by element of symmetry
        """
        [b1,r11,r12,r13,b2,r21,r22,r23,b3,r31,r32,r33]=elsymm
        matrix_r = numpy.array([[r11, r12, r13], [r21, r22, r23], 
                               [r31, r32, r33]], dtype=float)
        matrix_rt = matrix_r.transpose()
        r_chi = numpy.matmul(matrix_r, matrix_chi)
        
        matrix_chi_rot = numpy.matmul(r_chi, matrix_rt)
        return matrix_chi_rot 
    
    def _calc_form_factor(self, sthovl):
        """Calculate magnetic form factor in frame of Spherical model (Int.Tabl.C.p.592)\n
        LFactor is Lande factor\n
        coeff0 is a list [A,a,B,b,C,c,D] at n=0\n
        coeff2 is a list [A,a,B,b,C,c,D] at n=2\n
        lsthovl is list sin(theta)/lambda in Angstrems**-1

        Calculation of magnetic form factor <j0>,<j2>,<j4>,<j6>\n
        coeff is a list [A,a,B,b,C,c,D] at n=0,2,4,6
        For help see International Table Vol.C p.460
        """
        #not sure about kappa, it is here just for test, by default it is 1.0
        kappa = self.kappa
        factor_lande = self.factor_lande
        j0_A = self.j0_A
        j0_a = self.j0_a
        j0_B = self.j0_B
        j0_b = self.j0_b
        j0_C = self.j0_C
        j0_c = self.j0_c
        j0_D = self.j0_D
        j2_A = self.j2_A
        j2_a = self.j2_a
        j2_B = self.j2_B
        j2_b = self.j2_b
        j2_C = self.j2_C
        j2_c = self.j2_c
        j2_D = self.j2_D     
        
        np_sthovl, np_kappa = numpy.meshgrid(sthovl, kappa, indexing ="ij")
        np_factor_lande = numpy.meshgrid(sthovl, factor_lande, indexing ="ij")[1]
        np_j0_A = numpy.meshgrid(sthovl, j0_A, indexing ="ij")[1]
        np_j0_a = numpy.meshgrid(sthovl, j0_a, indexing ="ij")[1]
        np_j0_B = numpy.meshgrid(sthovl, j0_B, indexing ="ij")[1]
        np_j0_b = numpy.meshgrid(sthovl, j0_b, indexing ="ij")[1]
        np_j0_C = numpy.meshgrid(sthovl, j0_C, indexing ="ij")[1]
        np_j0_c = numpy.meshgrid(sthovl, j0_c, indexing ="ij")[1]
        np_j0_D = numpy.meshgrid(sthovl, j0_D, indexing ="ij")[1]
        np_j2_A = numpy.meshgrid(sthovl, j2_A, indexing ="ij")[1]
        np_j2_a = numpy.meshgrid(sthovl, j2_a, indexing ="ij")[1]
        np_j2_B = numpy.meshgrid(sthovl, j2_B, indexing ="ij")[1]
        np_j2_b = numpy.meshgrid(sthovl, j2_b, indexing ="ij")[1]
        np_j2_C = numpy.meshgrid(sthovl, j2_C, indexing ="ij")[1]
        np_j2_c = numpy.meshgrid(sthovl, j2_c, indexing ="ij")[1]
        np_j2_D = numpy.meshgrid(sthovl, j2_D, indexing ="ij")[1]

        np_h = (np_sthovl*1./np_kappa)**2
        j0_av = (np_j0_A*numpy.exp(-np_j0_a*np_h)+
                 np_j0_B*numpy.exp(-np_j0_b*np_h)+
                 np_j0_C*numpy.exp(-np_j0_c*np_h)+np_j0_D)
        j2_av = (np_j2_A*numpy.exp(-np_j2_a*np_h)+
                 np_j2_B*numpy.exp(-np_j2_b*np_h)+
                 np_j2_C*numpy.exp(-np_j2_c*np_h)+np_j2_D)*np_h
        #according https://journals.aps.org/prb/pdf/10.1103/PhysRevB.79.140405
        #mismatch with international tables where (1.0-2.0/np_factor_lande)
        form_factor = j0_av+(2.0/np_factor_lande-1.0)*j2_av
        return form_factor

