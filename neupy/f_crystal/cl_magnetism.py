"""
define classes to describe crystal 
"""
__author__ = 'ikibalin'
__version__ = "2019_04_06"
import os
import numpy

from neupy.f_common.cl_variable import Variable
from neupy.f_interface.cl_abstract_magnetism import AbstractMagnetism


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



class Magnetism(AbstractMagnetism):
    """
    Magnetism
    """
    def __init__(self, kappa=1.0, factor_lande=2.0, chi_11=0., chi_22=0., 
                 chi_33=0., chi_12=0., chi_13=0., chi_23=0., j0_A=0., j0_a=0., 
                 j0_B=0., j0_b=0., j0_C=0., j0_c=0., j0_D=0., j2_A=0., j2_a=0., 
                 j2_B=0., j2_b=0., j2_C=0., j2_c=0., j2_D=0.):
        super(Magnetism, self).__init__()
        self._p_chi_11 = None
        self._p_chi_22 = None
        self._p_chi_33 = None
        self._p_chi_12 = None
        self._p_chi_13 = None
        self._p_chi_23 = None
        self._p_kappa = None
        self._p_factor_lande = None
        self._p_j0_A = None
        self._p_j0_a = None
        self._p_j0_B = None
        self._p_j0_b = None
        self._p_j0_C = None
        self._p_j0_c = None
        self._p_j0_D = None
        self._p_j2_A = None
        self._p_j2_a = None
        self._p_j2_B = None
        self._p_j2_b = None
        self._p_j2_C = None
        self._p_j2_c = None
        self._p_j2_D = None
        self._refresh(chi_11, chi_22, chi_33, chi_12, chi_13, chi_23,kappa, 
                      factor_lande, j0_A, j0_a, j0_B, j0_b, j0_C, j0_c, j0_D, 
                      j2_A, j2_a, j2_B, j2_b, j2_C, j2_c, j2_D)
        
    def __repr__(self):
        lsout = """Magnetism: \n chi_11: {:}\n chi_22: {:}\n chi_33: {:}
 chi_12: {:}\n chi_13: {:}\n chi_23: {:}\n kappa: {:}
 factor_lande: {:}""".format(
 self._p_chi_11, self._p_chi_22, self._p_chi_33, self._p_chi_12, 
 self._p_chi_13, self._p_chi_23, self._p_kappa, self._p_factor_lande)
        return lsout


    def _refresh(self, chi_11, chi_22, chi_33, chi_12, chi_13, chi_23, kappa, 
                      factor_lande, j0_A, j0_a, j0_B, j0_b, j0_C, j0_c, j0_D, 
                      j2_A, j2_a, j2_B, j2_b, j2_C, j2_c, j2_D):
        
        if not(isinstance(chi_11, type(None))):
            self._p_chi_11 = chi_11
        if not(isinstance(chi_22, type(None))):
            self._p_chi_22 = chi_22
        if not(isinstance(chi_33, type(None))):
            self._p_chi_33 = chi_33
        if not(isinstance(chi_12, type(None))):
            self._p_chi_12 = chi_12
        if not(isinstance(chi_13, type(None))):
            self._p_chi_13 = chi_13
        if not(isinstance(chi_23, type(None))):
            self._p_chi_23 = chi_23
        if not(isinstance(kappa, type(None))):
            self._p_kappa = kappa 
        if not(isinstance(factor_lande, type(None))):
            self._p_factor_lande = factor_lande 
        if not(isinstance(j0_A, type(None))):
            self._p_j0_A = j0_A 
        if not(isinstance(j0_a, type(None))):
            self._p_j0_a = j0_a 
        if not(isinstance(j0_B, type(None))):
            self._p_j0_B = j0_B 
        if not(isinstance(j0_b, type(None))):
            self._p_j0_b = j0_b 
        if not(isinstance(j0_C, type(None))):
            self._p_j0_C = j0_C 
        if not(isinstance(j0_c, type(None))):
            self._p_j0_c = j0_c 
        if not(isinstance(j0_D, type(None))):
            self._p_j0_D = j0_D 
        if not(isinstance(j2_A, type(None))):
            self._p_j2_A = j2_A 
        if not(isinstance(j2_a, type(None))):
            self._p_j2_a = j2_a 
        if not(isinstance(j2_B, type(None))):
            self._p_j2_B = j2_B 
        if not(isinstance(j2_b, type(None))):
            self._p_j2_b = j2_b 
        if not(isinstance(j2_C, type(None))):
            self._p_j2_C = j2_C 
        if not(isinstance(j2_c, type(None))):
            self._p_j2_c = j2_c 
        if not(isinstance(j2_D, type(None))):
            self._p_j2_D = j2_D 
            

    def set_val(self, chi_11=None, chi_22=None, chi_33=None, chi_12=None, 
                chi_13=None, chi_23=None, kappa=None, factor_lande=None, 
                j0_A=None, j0_a=None, j0_B=None, j0_b=None, j0_C=None, 
                j0_c=None, j0_D=None, j2_A=None, j2_a=None, j2_B=None, 
                j2_b=None, j2_C=None, j2_c=None, j2_D=None):
        self._refresh(chi_11, chi_22, chi_33, chi_12, chi_13, chi_23, kappa, 
                      factor_lande, j0_A, j0_a, j0_B, j0_b, j0_C, j0_c, j0_D, 
                      j2_A, j2_a, j2_B, j2_b, j2_C, j2_c, j2_D)
        
        
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

chi_ij are the susceptibility vector
kappa define the size of the radial function (equals 1. by default)
factor_lande is the factor Lande (equals 2. by default)
        """
        print(lsout)
    
    def _calc_chi_loc(ia, ib, ic, matrix_ib):
        """
        representation of chi in crystallographic coordinate system defined as x||a*, z||c, y= [z x] (right handed)
        expressions are taken from international tables
        matrix_ib is inversed matrix B
        ia, ib, ic is inversed unit cell parameters (it can be estimated from matrix matrix_ib)

        X = B x, x = iB X
        xT*CHI*x = XT iBT CHI iB X
    
        output chiLOC = iBT CHI iB
        """
        matrix_chi = numpy.array(
                [[self["chi_11"], self["chi_12"], self["chi_13"]],
                 [self["chi_12"], self["chi_22"], self["chi_23"]],
                 [self["chi_13"], self["chi_23"], self["chi_33"]]], 
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
        matrix_ib_norm[:,0] *= ia
        matrix_ib_norm[:,1] *= ib
        matrix_ib_norm[:,2] *= ic
        
        matrix_ibt_norm = matrix_ib_norm.transpose()
        #it is not compatible with case, vhen chi_ij is 1D array 
        ibt_chi = numpy.matmul(matrix_ibt_norm, matrix_chi)
        matrix_chi_loc = numpy.matmul(ibt_chi, matrix_ib_norm)
        d_out = dict(matrix_chi_loc = matrix_chi_loc)
        self.update(d_out)
    
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
        
        r_11, r_12 = space_group.get_val("r_11"), space_group.get_val("r_12")
        r_13, r_21 = space_group.get_val("r_13"), space_group.get_val("r_21")
        r_22, r_23 = space_group.get_val("r_22"), space_group.get_val("r_23")
        r_31, r_32 = space_group.get_val("r_31"), space_group.get_val("r_32")
        r_33 = space_group.get_val("r_33")

        chi_11, chi_22 = self._p_chi_11, self._p_chi_22 
        chi_33, chi_12 = self._p_chi_33, self._p_chi_12
        chi_13, chi_23 = self._p_chi_13, self._p_chi_23
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
    
    def calc_chi_rot(matrix_chi, elsymm):
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
        kappa = self._p_kappa
        factor_lande = self._p_factor_lande
        j0_A = self._p_j0_A
        j0_a = self._p_j0_a
        j0_B = self._p_j0_B
        j0_b = self._p_j0_b
        j0_C = self._p_j0_C
        j0_c = self._p_j0_c
        j0_D = self._p_j0_D
        j2_A = self._p_j2_A
        j2_a = self._p_j2_a
        j2_B = self._p_j2_B
        j2_b = self._p_j2_b
        j2_C = self._p_j2_C
        j2_c = self._p_j2_c
        j2_D = self._p_j2_D     
        
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

