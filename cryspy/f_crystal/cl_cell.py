__author__ = 'ikibalin'
__version__ = "2019_08_25"
import os
import numpy
from pycifstar import Global

from cryspy.f_common.cl_fitable import Fitable
from cryspy.f_crystal.cl_space_group import SpaceGroup

class Cell(object):
    """
    Data items in the Cell class record details about the
    crystallographic cell parameters and their measurement.

    Description in cif file:

        _cell_length_a                     5.959(1)
        _cell_length_b                     14.956(1)
        _cell_length_c                     19.737(3)
        _cell_angle_alpha                  90
        _cell_angle_beta                   90
        _cell_angle_gamma                  90
    """
    def __init__(self, a = Fitable(1.0), 
                       b = Fitable(1.0), 
                       c = Fitable(1.0), 
                       alpha = Fitable(90.0), 
                       beta = Fitable(90.0), 
                       gamma= Fitable(90.), bravais_lattice = "triclinic"):
        super(Cell, self).__init__()
        self.__cell_length_a = None
        self.__cell_length_b = None
        self.__cell_length_c = None
        self.__cell_angle_alpha = None
        self.__cell_angle_beta = None
        self.__cell_angle_gamma = None
        self.__cell_setting = None

        self.__cos_a = None
        self.__cos_b = None
        self.__cos_g = None
        self.__cos_a_sq = None
        self.__cos_b_sq = None
        self.__cos_g_sq = None
        self.__sin_a = None
        self.__sin_b = None
        self.__sin_g = None
        self.__sin_a_sq = None
        self.__sin_b_sq = None
        self.__sin_g_sq = None
        
        self.__cell_reciprocal_length_a = None
        self.__cell_reciprocal_length_b = None
        self.__cell_reciprocal_length_c = None
        self.__cell_reciprocal_angle_alpha = None
        self.__cell_reciprocal_angle_beta = None
        self.__cell_reciprocal_angle_gamma = None        

        self.__cos_ia = None
        self.__cos_ib = None
        self.__cos_ig = None
        self.__cos_ia_sq = None
        self.__cos_ib_sq = None
        self.__cos_ig_sq = None
        self.__sin_ia = None
        self.__sin_ib = None
        self.__sin_ig = None
        self.__sin_ia_sq = None
        self.__sin_ib_sq = None
        self.__sin_ig_sq = None
        
        self.__cell_volume = None
        self.__cell_ivolume = None
        self.__m_b = None
        self.__m_ib = None

        self.a = a # type: Fitable
        self.b = b # type: Fitable
        self.c = c # type: Fitable
        self.alpha = alpha # type: Fitable
        self.beta = beta # type: Fitable
        self.gamma = gamma # type: Fitable
        self.bravais_lattice = bravais_lattice # type: str
        
    @property
    def a(self):
        """
        Unit-cell lengths in angstroms corresponding to the structure
        reported. 
        The permitted range is 0.0 -> infinity

        reference: https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Icell_length_.html
        """
        return self.__cell_length_a
    @a.setter
    def a(self, x):
        if isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        self.__cell_length_a = x_in
        self.apply_constraint()        
    @property
    def b(self):
        """
        Unit-cell lengths in angstroms corresponding to the structure
        reported. 
        The permitted range is 0.0 -> infinity
        
        reference: https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Icell_length_.html
        """
        return self.__cell_length_b
    @b.setter
    def b(self, x):
        if isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        self.__cell_length_b = x_in
        self.apply_constraint()        
    @property
    def c(self):
        """
        Unit-cell lengths in angstroms corresponding to the structure
        reported. 
        The permitted range is 0.0 -> infinity
        
        reference: https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Icell_length_.html
        """
        return self.__cell_length_c
    @c.setter
    def c(self, x):
        if isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        self.__cell_length_c = x_in
        self.apply_constraint()        
    @property
    def alpha(self):
        """
        Unit-cell angles of the reported structure in degrees.
        The permitted range is 0.0 -> 180.0 
        Enumeration default: 90.0
        
        reference: https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Icell_angle_.html
        """
        return self.__cell_angle_alpha
    @alpha.setter
    def alpha(self, x):
        if isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        self.__cell_angle_alpha = x_in
        self.apply_constraint()        
    @property
    def beta(self):
        """
        Unit-cell angles of the reported structure in degrees.
        The permitted range is 0.0 -> 180.0 
        Enumeration default: 90.0
        
        reference: https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Icell_angle_.html
        """
        return self.__cell_angle_beta
    @beta.setter
    def beta(self, x):
        if isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        self.__cell_angle_beta = x_in
        self.apply_constraint()        
    @property
    def gamma(self):
        """
        Unit-cell angles of the reported structure in degrees.
        The permitted range is 0.0 -> 180.0 
        Enumeration default: 90.0
        
        reference: https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Icell_angle_.html
        """
        return self.__cell_angle_gamma
    @gamma.setter
    def gamma(self, x):
        if isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        self.__cell_angle_gamma = x_in
        self.apply_constraint()        
    @property
    def bravais_lattice(self):
        """
        The cell settings for this space-group symmetry.

        The data value must be one of the following:

        triclinic	
        monoclinic	
        orthorhombic	
        tetragonal	
        rhombohedral	
        trigonal	
        hexagonal	
        cubic

        reference: https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Isymmetry_cell_setting.html
        """
        return self.__cell_setting
    @bravais_lattice.setter
    def bravais_lattice(self, x):
        try:
            x_in = str(x).lower()
        except:
            x_in = "triclinic"
        l_bravais_lattice = ["cubic", "hexagonal", "rhombohedral", "trigonal", "tetragonal", 
                             "orthorhombic", "monoclinic", "triclinic"]
        if x_in not in l_bravais_lattice:
            print("Introduced bravais_lattice is not found.")
            print("Please try one of them: {:}.".format(", ".join(l_bravais_lattice)))
            x_in = "Triclinic"
        self.__cell_setting = x_in
        self.apply_constraint()        
    @property
    def volume(self):
        """
        Cell volume V in angstroms cubed.
        
        V = a b c [1 - cos^2^(alpha) - cos^2^(beta) - cos^2^(gamma)
               + 2 cos(alpha) cos(beta) cos(gamma) ] ^1/2^

        The permitted range is 0.0 -> infinity

        reference: https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Icell_volume.html
        """
        return self.__cell_volume
    @property
    def ivolume(self):
        """
        The reciprocal-cell volume V in angstroms cubed.
        (not realized)
        """
        return self.__cell_ivolume
    @property
    def ia(self):
        """
        The reciprocal-cell lengths in inverse angstroms.  These are
        related to the real cell by:

        recip-a = b*c*sin(alpha)/V
        recip-b = c*a*sin(beta)/V
        recip-c = a*b*sin(gamma)/V

        where V is the cell volume.

        Ref: Buerger, M. J. (1942). X-ray Crystallography, p. 360.
           New York: John Wiley & Sons Inc.

        The permitted range is 0.0 -> infinity

        reference: https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Icell_reciprocal_length_.html
        """
        return self.__cell_reciprocal_length_a
    @property
    def ib(self):
        """
        The reciprocal-cell lengths in inverse angstroms.  These are
        related to the real cell by:

        recip-a = b*c*sin(alpha)/V
        recip-b = c*a*sin(beta)/V
        recip-c = a*b*sin(gamma)/V

        where V is the cell volume.

        Ref: Buerger, M. J. (1942). X-ray Crystallography, p. 360.
           New York: John Wiley & Sons Inc.

        The permitted range is 0.0 -> infinity
        
        reference: https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Icell_reciprocal_length_.html
        """
        return self.__cell_reciprocal_length_b
    @property
    def ic(self):
        """
        The reciprocal-cell lengths in inverse angstroms.  These are
        related to the real cell by:

        recip-a = b*c*sin(alpha)/V
        recip-b = c*a*sin(beta)/V
        recip-c = a*b*sin(gamma)/V

        where V is the cell volume.

        Ref: Buerger, M. J. (1942). X-ray Crystallography, p. 360.
           New York: John Wiley & Sons Inc.

        The permitted range is 0.0 -> infinity
        
        reference: https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Icell_reciprocal_length_.html
        """
        return self.__cell_reciprocal_length_c
    @property
    def ialpha(self):
        """
        The angles defining the reciprocal cell in degrees. These
        are related to those in the real cell by:

        cos(recip-alpha)
          = [cos(beta)*cos(gamma) - cos(alpha)]/[sin(beta)*sin(gamma)]

        cos(recip-beta)
          = [cos(gamma)*cos(alpha) - cos(beta)]/[sin(gamma)*sin(alpha)]

        cos(recip-gamma)
         = [cos(alpha)*cos(beta) - cos(gamma)]/[sin(alpha)*sin(beta)]

        Ref: Buerger, M. J. (1942). X-ray Crystallography, p. 360.
           New York: John Wiley & Sons Inc.

        The permitted range is 0.0 -> 180.0 
        Enumeration default: 90.0

        reference: https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Icell_reciprocal_angle_.html
        """
        return self.__cell_reciprocal_angle_alpha
    @property
    def ibeta(self):
        """
        The angles defining the reciprocal cell in degrees. These
        are related to those in the real cell by:

        cos(recip-alpha)
          = [cos(beta)*cos(gamma) - cos(alpha)]/[sin(beta)*sin(gamma)]

        cos(recip-beta)
          = [cos(gamma)*cos(alpha) - cos(beta)]/[sin(gamma)*sin(alpha)]

        cos(recip-gamma)
         = [cos(alpha)*cos(beta) - cos(gamma)]/[sin(alpha)*sin(beta)]

        Ref: Buerger, M. J. (1942). X-ray Crystallography, p. 360.
           New York: John Wiley & Sons Inc.

        The permitted range is 0.0 -> 180.0 
        Enumeration default: 90.0

        reference: https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Icell_reciprocal_angle_.html
        """
        return self.__cell_reciprocal_angle_beta
    @property
    def igamma(self):
        """
        The angles defining the reciprocal cell in degrees. These
        are related to those in the real cell by:

        cos(recip-alpha)
          = [cos(beta)*cos(gamma) - cos(alpha)]/[sin(beta)*sin(gamma)]

        cos(recip-beta)
          = [cos(gamma)*cos(alpha) - cos(beta)]/[sin(gamma)*sin(alpha)]

        cos(recip-gamma)
         = [cos(alpha)*cos(beta) - cos(gamma)]/[sin(alpha)*sin(beta)]

        Ref: Buerger, M. J. (1942). X-ray Crystallography, p. 360.
           New York: John Wiley & Sons Inc.

        The permitted range is 0.0 -> 180.0 
        Enumeration default: 90.0

        reference: https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Icell_reciprocal_angle_.html
        """
        return self.__cell_reciprocal_angle_gamma
    @property
    def m_b(self):
        """
        B matrix
        """
        return self.__m_b
    @property
    def m_ib(self):
        """
        Inversed B matrix
        """
        return self.__m_ib
        
    @property
    def cos_a(self):
        """
        cos(alpha)
        """
        return self.__cos_a
    @property
    def cos_b(self):
        """
        cos(beta)
        """
        return self.__cos_b
    @property
    def cos_g(self):
        """
        cos(gamma)
        """
        return self.__cos_g
    @property
    def cos_ib(self):
        """
        cos(beta*)
        """
        return self.__cos_ib
    @property
    def sin_ib(self):
        """
        sin(beta*)
        """
        return self.__sin_ib
    @property
    def cos_ig(self):
        """
        cos(gamma*)
        """
        return self.__cos_ig
    @property
    def sin_ig(self):
        """
        sin(gamma*)
        """
        return self.__sin_ig
        
    def __repr__(self):
        ls_out = ["""Cell: \n a: {:}\n b: {:}\n c: {:}\n alpha: {:}
 beta: {:}\n gamma: {:}\n bravais_lattice: {:}""".format(self.a, self.b, 
                 self.c, self.alpha, self.beta, self.gamma, 
                 self.bravais_lattice)]
        if self.__cell_volume is not None:
            ls_out.append(" volume: {:.3f}".format(self.__cell_volume))

        if self.__m_b is not None:
             ls_out.append(""" B matrix is:\n {:9.5f} {:9.5f} {:9.5f}
 {:9.5f} {:9.5f} {:9.5f}\n {:9.5f} {:9.5f} {:9.5f}""".format(self.__m_b[0, 0],
 self.__m_b[0, 1], self.__m_b[0, 2], self.__m_b[1, 0], self.__m_b[1, 1], 
 self.__m_b[1, 2], self.__m_b[2, 0], self.__m_b[2, 1], self.__m_b[2, 2]))
        if self.__m_ib is not None:
             ls_out.append(""" inversed B matrix is:\n {:9.5f} {:9.5f} {:9.5f}
 {:9.5f} {:9.5f} {:9.5f}\n {:9.5f} {:9.5f} {:9.5f}""".format(self.__m_ib[0, 0],
 self.__m_ib[0, 1], self.__m_ib[0, 2], self.__m_ib[1, 0], self.__m_ib[1, 1], 
 self.__m_ib[1, 2], self.__m_ib[2, 0], self.__m_ib[2, 1], self.__m_ib[2, 2]))
        return "\n".join(ls_out)
    
    
    def _constr_bravais_lattice(self):
        bravais_lattice = self.bravais_lattice
        if bravais_lattice == "cubic":
            self.a.constraint_flag = False
            self.__cell_length_b = Fitable(self.a.value, self.a.sigma, False, "_cell_length_b", True)
            self.__cell_length_c = Fitable(self.a.value, self.a.sigma, False, "_cell_length_c", True)
            self.__cell_angle_alpha = Fitable(90., None, False, "_cell_angle_alpha", True)
            self.__cell_angle_beta = Fitable(90., None, False, "_cell_angle_beta", True)
            self.__cell_angle_gamma = Fitable(90., None, False, "_cell_angle_gamma", True)
        elif bravais_lattice == "hexagonal":
            self.a.constraint_flag = False
            self.c.constraint_flag = False
            self.__cell_length_b = Fitable(self.a.value, self.a.sigma, False, "_cell_length_b", True)
            self.__cell_angle_alpha = Fitable(90., None, False, "_cell_angle_alpha", True)
            self.__cell_angle_beta = Fitable(90., None, False, "_cell_angle_beta", True)
            self.__cell_angle_gamma = Fitable(120., None, False, "_cell_angle_gamma", True)
        elif bravais_lattice == "rhombohedral":
            self.a.constraint_flag = False
            self.alpha.constraint_flag = False
            self.__cell_length_b = Fitable(self.a.value, self.a.sigma, False, "_cell_length_b", True)
            self.__cell_length_c = Fitable(self.a.value, self.a.sigma, False, "_cell_length_c", True)
            self.__cell_angle_beta = Fitable(self.alpha.value, None, False, "_cell_angle_beta", True)
            self.__cell_angle_gamma = Fitable(self.alpha.value, None, False, "_cell_angle_gamma", True)
        elif bravais_lattice == "trigonal":
            self.a.constraint_flag = False
            self.alpha.constraint_flag = False
            self.beta.constraint_flag = False
            self.gamma.constraint_flag = False
            self.__cell_length_b = Fitable(self.a.value, self.a.sigma, False, "_cell_length_b", True)
            self.__cell_length_c = Fitable(self.a.value, self.a.sigma, False, "_cell_length_c", True)
        elif bravais_lattice == "tetragonal":
            self.a.constraint_flag = False
            self.c.constraint_flag = False
            self.__cell_length_b = Fitable(self.a.value, self.a.sigma, False, "_cell_length_b", True)
            self.__cell_angle_alpha = Fitable(90., None, False, "_cell_angle_alpha", True)
            self.__cell_angle_beta = Fitable(90., None, False, "_cell_angle_beta", True)
            self.__cell_angle_gamma = Fitable(90., None, False, "_cell_angle_gamma", True)
        elif bravais_lattice == "orthorhombic":
            self.a.constraint_flag = False
            self.b.constraint_flag = False
            self.c.constraint_flag = False
            self.__cell_angle_alpha = Fitable(90., None, False, "_cell_angle_alpha", True)
            self.__cell_angle_beta = Fitable(90., None, False, "_cell_angle_beta", True)
            self.__cell_angle_gamma = Fitable(90., None, False, "_cell_angle_gamma", True)
        elif bravais_lattice == "monoclinic":
            self.a.constraint_flag = False
            self.b.constraint_flag = False
            self.c.constraint_flag = False
            self.beta.constraint_flag = False
            self.__cell_angle_alpha = Fitable(90., None, False, "_cell_angle_alpha", True)
            self.__cell_angle_gamma = Fitable(90., None, False, "_cell_angle_gamma", True)

    def _calc_cos_abc(self):
        rad=numpy.pi/180.
        self.__cos_a = numpy.cos(self.alpha*rad)
        self.__cos_b = numpy.cos(self.beta*rad)
        self.__cos_g = numpy.cos(self.gamma*rad)
        
        self.__sin_a = numpy.sin(self.alpha*rad)
        self.__sin_b = numpy.sin(self.beta*rad)
        self.__sin_g = numpy.sin(self.gamma*rad)
        
        self.__cos_a_sq = self.__cos_a**2
        self.__cos_b_sq = self.__cos_b**2
        self.__cos_g_sq = self.__cos_g**2

        self.__sin_a_sq = 1.-self.__cos_a_sq
        self.__sin_b_sq = 1.-self.__cos_b_sq
        self.__sin_g_sq = 1.-self.__cos_g_sq
        
    def _calc_cos_iabc(self):
        rad=numpy.pi/180.
        self.__cos_ia = numpy.cos(self.__cell_reciprocal_angle_alpha*rad)
        self.__cos_ib = numpy.cos(self.__cell_reciprocal_angle_beta*rad)
        self.__cos_ig = numpy.cos(self.__cell_reciprocal_angle_gamma*rad)
        
        self.__sin_ia = numpy.sin(self.__cell_reciprocal_angle_alpha*rad)
        self.__sin_ib = numpy.sin(self.__cell_reciprocal_angle_beta*rad)
        self.__sin_ig = numpy.sin(self.__cell_reciprocal_angle_gamma*rad)
        
        self.__cos_ia_sq = self.__cos_ia**2
        self.__cos_ib_sq = self.__cos_ib**2
        self.__cos_ig_sq = self.__cos_ig**2

        self.__sin_a_sq = 1.-self.__cos_a_sq
        self.__sin_b_sq = 1.-self.__cos_b_sq
        self.__sin_g_sq = 1.-self.__cos_g_sq

    def _calc_volume(self):
        a = 1.*self.a
        b = 1.*self.b
        c = 1.*self.c
        c_a = self.__cos_a
        c_b = self.__cos_b
        c_g = self.__cos_g
        c_a_sq = self.__cos_a_sq
        c_b_sq = self.__cos_b_sq
        c_g_sq = self.__cos_g_sq
        vol = a*b*c*(1.-c_a_sq-c_b_sq-c_g_sq+2.*c_a*c_b*c_g)**0.5
        self.__cell_volume = vol
        
    
    def _calc_iucp(self):
        """
        calculate inverse unit cell
        """
        irad = 180./numpy.pi

        a = 1.*self.a
        b = 1.*self.b
        c = 1.*self.c
        c_a = self.__cos_a
        c_b = self.__cos_b
        c_g = self.__cos_g
        s_a = self.__sin_a
        s_b = self.__sin_b
        s_g = self.__sin_g
        vol = self.__cell_volume
        
        self.__cell_reciprocal_angle_alpha = numpy.arccos((c_b*c_g-c_a)/(s_b*s_g))*irad
        self.__cell_reciprocal_angle_beta = numpy.arccos((c_g*c_a-c_b)/(s_g*s_a))*irad
        self.__cell_reciprocal_angle_gamma = numpy.arccos((c_a*c_b-c_g)/(s_a*s_b))*irad

        self.__cell_reciprocal_length_a = b*c*s_a/vol
        self.__cell_reciprocal_length_b = c*a*s_b/vol
        self.__cell_reciprocal_length_c = a*b*s_g/vol


    def _calc_m_b(self):
        """
        calculate matrix B 
        """
        c = 1.*self.c

        ia = self.__cell_reciprocal_length_a 
        ib = self.__cell_reciprocal_length_b 
        ic = self.__cell_reciprocal_length_c 
        
        c_a = self.__cos_a
        
        #ic_a = self._p_cos_ia 
        ic_b = self.__cos_ib 
        ic_g = self.__cos_ig 
        #is_a = self._p_sin_ia 
        is_b = self.__sin_ib 
        is_g = self.__sin_ig 
        
        self.__m_b = numpy.array([[ia,  ib*ic_g,  ic*ic_b],
            [0.,  ib*is_g, -ic*is_b*c_a],
            [0.,       0.,  1./c]], dtype = float)

    def _calc_m_ib(self):
        """
        calculate inverse B matrix 
        """
        x1 = self.__m_b[0,0]
        x2 = self.__m_b[1,1]
        x3 = self.__m_b[2,2]
        x4 = self.__m_b[0,1]
        x5 = self.__m_b[0,2]
        x6 = self.__m_b[1,2]
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
        
        self.__m_ib = numpy.array([[y1,y4,y5],[0.,y2,y6],[0.,0.,y3]], 
                                   dtype = float)
            
                
    def _show_message(self, s_out: str):
        print("***  Error ***")
        print(s_out)

    
    def calc_sthovl(self, h=None, k=None, l=None, hkl=None, l_hkl=None, f_print=False):
        """
        Calculate sin(theta)/lambda for list of hkl reflections.

        Keyword arguments:
        h, k, l --- Miller indices (example: h=1, k=0, l=2)
        or
        hkl --- a turple of Miller indices (example: hkl=(1, 0, 2))
        or
        l_hkl --- a list of turples of Miller indices (example: l_hkl=[(1, 0, 2), (1, 1, 3)])
        f_print --- a flag to print output information in terminal
        """
        cond_1 = all([h is not None, k is not None, l is not None])
        cond_2 = hkl is not None
        cond_3 = l_hkl is not None
        if cond_1:
            np_h = h # type: numpy.array or float
            np_k = k # type: numpy.array or float
            np_l = l # type: numpy.array or float
        elif cond_2:
            np_h = hkl[0] # type: float
            np_k = hkl[1] # type: float
            np_l = hkl[2] # type: float
        elif cond_3:
            np_h = numpy.array([hh[0] for hh in l_hkl], dtype=float) # type: numpy.array
            np_k = numpy.array([hh[1] for hh in l_hkl], dtype=float) # type: numpy.array
            np_l = numpy.array([hh[2] for hh in l_hkl], dtype=float) # type: numpy.array
        else: 
            self._show_message("Did not found correct input. Expected h, k, l or hkl or l_hkl")
            return
        a = 1.*self.a
        b = 1.*self.b
        c = 1.*self.c
        c_a = self.__cos_a
        c_b = self.__cos_b
        c_g = self.__cos_g
        c_a_sq = self.__cos_a_sq
        c_b_sq = self.__cos_b_sq
        c_g_sq = self.__cos_g_sq
        s_a_sq = self.__sin_a_sq
        s_b_sq = self.__sin_b_sq
        s_g_sq = self.__sin_g_sq

        A=( 1. - c_a_sq - c_b_sq - c_g_sq + 2.*c_a*c_b*c_g)
        B1 = (s_a_sq*(np_h*1./a)**2+s_b_sq*(np_k*1./b)**2+s_g_sq*(np_l*1./c)**2)
        B2 = 2.*(np_k*np_l*c_a)/(b*c)+2.*(np_h*np_l*c_b)/(a*c)+2.*(np_h*np_k*c_g)/(a*b)
        #it should be checked, I am not sure
        B = B1-B2
        inv_d = (B*1./A)**0.5
        res = 0.5*inv_d
        if f_print:
            ls_out = ["    h     k     l     sthovl"]
            try:
                _ = (hh for hh in res)
                f_iter = True # type: bool
            except TypeError:
                f_iter = False # type: bool
            if f_iter:
                for hh_1, hh_2, hh_3, hh_4 in zip(np_h, np_k, np_l, res):
                    ls_out.append("{:5.1f} {:5.1f} {:5.1f} {:10.5f}".format(hh_1, hh_2, hh_3, hh_4))
            else:
                ls_out.append("{:5.1f} {:5.1f} {:5.1f} {:10.5f}".format(np_h, np_k, np_l, res))
            print("\n".join(ls_out))
        return res

    def calc_k_loc(self, h, k, l):
        """
        Calculate unity scattering vector.

        Keyword arguments:
        h, k, l -- Miller indices
        """
        m_b = self.__m_b
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
        """
        define rotation matrix to have new z axis along kloc
        Rotation matrix is defined by Euler angles
        """
        m_b = self.__m_b
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

    @property
    def is_variable(self):
        """
        Output: True if there is any refined parameter
        """
        res = any([self.a.refinement,
                   self.b.refinement,
                   self.c.refinement,
                   self.alpha.refinement,
                   self.beta.refinement,
                   self.gamma.refinement])
        return res
        
    def get_variables(self):
        """
        Output: the list of the refined parameters
        """
        l_variable = []
        if self.a.refinement:
            l_variable.append(self.a)
        if self.b.refinement:
            l_variable.append(self.b)
        if self.c.refinement:
            l_variable.append(self.c)
        if self.alpha.refinement:
            l_variable.append(self.alpha)
        if self.beta.refinement:
            l_variable.append(self.beta)
        if self.gamma.refinement:
            l_variable.append(self.gamma)
        return l_variable
    
    def apply_constraint(self):
        if self.is_defined:
            self._constr_bravais_lattice()
            self._calc_cos_abc()
            self._calc_volume()
            self._calc_iucp()
            self._calc_cos_iabc()
            self._calc_m_b()
            self._calc_m_ib()

    @property
    def is_defined(self):
        """
        Output: True if all started parameters are given
        """

        cond = any([self.a is None, self.b is None, self.c is None,
                    self.alpha is None, self.beta is None, self.gamma is None,
                    self.bravais_lattice is None])
        return not(cond)

    def calc_hkl(self, space_group, sthovl_min, sthovl_max):
        """
        a list of reflections hkl for cell in the range sthovl_min, sthovl_max
        taking into account the space group
        """
        if not(self.is_defined):
            print("Object 'Cell' is not fully defined for calculations")
            return None
        lhkl,lmult=[],[]
        l_hklres=[]

        hmax = int(2.*self.a*sthovl_max)
        kmax = int(2.*self.b*sthovl_max)
        lmax = int(2.*self.c*sthovl_max)
        hmin, kmin, lmin = -1*hmax, -1*kmax, -1*lmax

        hmin=0
        
        l_orig = space_group.orig
        l_symm = space_group.el_symm
        for h in range(hmin,hmax+1,1):
            for k in range(kmin,kmax+1,1):
                for l in range(lmin,lmax+1,1):
                    flag=(abs(sum([numpy.exp(2.*numpy.pi*1j*(orig[0]*h+orig[1]*k+orig[2]*l)) for orig in l_orig]))>0.00001)
                    #flag=True
                    if (flag):
                        lhkls=[(h*symm[1]+k*symm[5]+l*symm[9], h*symm[2]+k*symm[6]+l*symm[10], h*symm[3]+k*symm[7]+l*symm[11]) for symm in l_symm]
                        lhkls.extend([(-hkl[0],-hkl[1],-hkl[2]) for hkl in lhkls])
                        lhkls.sort(key=lambda x:10000*x[0]+100*x[1]+x[2])
                        if (not(lhkls[-1] in lhkl)):
                            lhkl.append(lhkls[-1])
                            lmult.append(len(set(lhkls)))
                            
        l_hklsthovl=[(hkl, self.calc_sthovl(hkl[0], hkl[1], hkl[2]), mult) for hkl, mult in zip(lhkl, lmult)]
        l_hklsthovl.sort(key=lambda x: x[1])
        l_hklres = [hklsthovl[0] for hklsthovl in l_hklsthovl if ((hklsthovl[1]>sthovl_min) & (hklsthovl[1]<sthovl_max))]
        l_multres = [hklsthovl[2] for hklsthovl in l_hklsthovl if ((hklsthovl[1]>sthovl_min) & (hklsthovl[1]<sthovl_max))]

        h = numpy.array([hh[0] for hh in l_hklres], dtype=int)
        k = numpy.array([hh[1] for hh in l_hklres], dtype=int)
        l = numpy.array([hh[2] for hh in l_hklres], dtype=int)
        mult = numpy.array(l_multres, dtype=int)
        return h, k, l, mult


    def calc_hkl_in_range(self, sthovl_min, sthovl_max):
        """
        give a list of reflections hkl for cell in the range sthovl_min, sthovl_max
        """
        if not(self.is_defined):
            print("Object 'Cell' is not fully defined for calculations")
            return None
        h_max = int(2.*self.a*sthovl_max)
        k_max = int(2.*self.b*sthovl_max)
        l_max = int(2.*self.c*sthovl_max)
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

    @property
    def to_cif(self):
        ls_out = ["_cell_length_a {:}".format(self.a.print_with_sigma)]
        ls_out.append("_cell_length_b {:}".format(self.b.print_with_sigma))
        ls_out.append("_cell_length_c {:}".format(self.c.print_with_sigma))
        ls_out.append("_cell_angle_alpha {:}".format(self.alpha.print_with_sigma))
        ls_out.append("_cell_angle_beta {:}".format(self.beta.print_with_sigma))
        ls_out.append("_cell_angle_gamma {:}".format(self.beta.print_with_sigma))
        return "\n".join(ls_out)
    
    def from_cif(self, string: str):
        cif_global = Global()
        flag = cif_global.take_from_string(string)
        if not flag:
            return False
        flag = False
        if cif_global.is_value("cell_length_a"):
            self.a = cif_global["cell_length_a"] # CIFvalue
        if cif_global.is_value("cell_length_b"):
            self.b = cif_global["cell_length_b"] # CIFvalue
        if cif_global.is_value("cell_length_c"):
            self.c = cif_global["cell_length_c"] # CIFvalue
        if cif_global.is_value("cell_angle_alpha"):
            self.alpha = cif_global["cell_angle_alpha"] # CIFvalue
        if cif_global.is_value("cell_angle_beta"):
            self.beta = cif_global["cell_angle_beta"] # CIFvalue
        if cif_global.is_value("cell_angle_gamma"):
            self.gamma = cif_global["cell_angle_gamma"] # CIFvalue
        if (cif_global.is_value("_space_group_name_H-M_alt") & cif_global.is_value("_space_group_it_coordinate_system_code")):
            space_groupe = SpaceGroup(spgr_given_name=cif_global["_space_group_name_H-M_alt"].value, 
                                      spgr_choice=cif_global["_space_group_it_coordinate_system_code"].value)
            self.bravais_lattice =space_groupe.bravais_lattice

        return True

