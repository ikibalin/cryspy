"""
define resolution for the powder diffractometer along ttheta
"""
__author__ = 'ikibalin'
__version__ = "2019_04_06"
import os
import numpy


class Cell(dict):
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
        lsout = """Unit cell: \n a: {:}\n b: {:}\n c: {:}\n alpha: {:}
 beta: {:}\n gamma: {:}\n singony: {:}""".format(self._p_a, self._p_b, 
                 self._p_c, self._p_alpha, self._p_beta, self._p_gamma, 
                 self._p_singony)
        return lsout
    
    def _refresh(self, a, b, c, alpha, beta, gamma, singony):
        """
        refresh variables
        """
        if a != None:
            self._p_a = a
        if b != None:
            self._p_b = b
        if c != None:
            self._p_c = c
        if alpha != None:
            self._p_alpha = alpha
        if beta != None:
            self._p_beta = beta
        if gamma != None:
            self._p_gamma = gamma
        if singony != None:
            self._p_singony = singony

        cond = any([hh != None for hh in [a, b, c, alpha, beta, gamma, 
                                          singony]])
        if cond:
            self._constr_singony()
    
    def _constr_singony(self):
        singony = self._p_singony
        if singony == "Cubic":
            self._p_b = self._p_a
            self._p_c = self._p_a
            self._p_alpha = 90.
            self._p_beta = 90.
            self._p_gamma = 90.
        elif singony == "Hexagonal":
            self._p_b = self._p_a
            self._p_alpha = 90.
            self._p_beta = 90.
            self._p_gamma = 120.        
        elif singony == "Trigonal":
            self._p_b = self._p_a
            self._p_c = self._p_a
        elif singony == "Tetragonal":
            self._p_b = self._p_a
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
        self._calc_cos_abc()
        self._calc_volume()
        self._calc_iucp()
        self._calc_cos_iabc()
        self._calc_m_b()
        self._calc_m_ib()

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
        a = self._p_a
        b = self._p_b
        c = self._p_c
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

        a = self._p_a
        b = self._p_b
        c = self._p_c
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
        c = self._p_c

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
            
                
        
    
    def calc_sthovl(self, h, k, l, a = None, b = None, c = None, alpha = None, 
                   beta = None, gamma= None, singony = None):
        """
        calculate sin(theta)/lambda for list of hkl reflections
        """
        cond = any([hh != None for hh in [a, b, c, alpha, beta, gamma, 
                                          singony]])
        if cond:
            self.set_val(a=a, b=b, c=c, alpha=alpha, beta=beta, gamma=gamma, 
                         singony=singony)
            
        a = self._p_a
        b = self._p_b
        c = self._p_c
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


    def calc_q_hkl(self, h, k, l, a = None, b = None, c = None, alpha = None, 
                   beta = None, gamma= None, singony = None):
        """
        calculate sin(theta)/lambda for list of hkl reflections
        """
        cond_1 = any([hh != None for hh in [a, b, c, alpha, beta, gamma, singony]])
        
        cond = cond_1
        if cond:
            self.set_val(a=a, b=b, c=c, alpha=alpha, beta=beta, gamma=gamma, singony=singony)

        matrix_ib = self.get_val("m_b")
        v_ia = matrix_ib[:,0]
        v_ib = matrix_ib[:,1]
        v_ic = matrix_ib[:,2]
        #mesh grid
        q_hkl = h*v_ia+k*v_ib+l*v_ic
        return q_hkl 
        
        

class SpaceGroupe(dict):
    """
    Space Groupe
    """
    def __init__(self, spgr_given_name = "P1", spgr_choice = "1",
                 f_dir_prog = os.getcwd()):
        super(SpaceGroupe, self).__init__()
        
        if isinstance(spgr_choice, float):
            spgr_choice = "{:}".format(int(spgr_choice))
        
        self._p_spgr_given_name = None
        self._p_spgr_choice = None
        self._p_f_dir_prog = None
        self._p_spgr_table = None

        self._p_centr = None
        self._p_el_symm = None
        self._p_orig = None
        self._p_p_centr = None
        self._p_spgr_name = None
        self._p_spgr_number = None

        self._p_r_11 = None
        self._p_r_12 = None
        self._p_r_13 = None
        self._p_r_21 = None
        self._p_r_22 = None
        self._p_r_23 = None
        self._p_r_31 = None
        self._p_r_32 = None
        self._p_r_33 = None

        self._p_b_1 = None
        self._p_b_2 = None
        self._p_b_3 = None

        f_itables=os.path.join(f_dir_prog,"itables.txt")
        self._read_el_cards(f_itables)        
        self._refresh(spgr_given_name, spgr_choice, f_dir_prog)
        self.set_val()
        
    def __repr__(self):
        lsout = """Space group: \n name: {:}\n choiсe: {:}
 directory: '{:}'""".format(self._p_spgr_given_name, self._p_spgr_choice, 
                            self._p_f_dir_prog)
        return lsout

    def _refresh(self, spgr_given_name, spgr_choice, f_dir_prog):
        if spgr_given_name != None:
            self._p_spgr_given_name = spgr_given_name
        if spgr_choice != None:
            self._p_spgr_choice = spgr_choice
        if f_dir_prog != None:
            self._p_f_dir_prog = f_dir_prog
            

    def set_val(self, spgr_given_name = None, spgr_choice = None,
                   f_dir_prog = None):
        self._refresh(spgr_given_name, spgr_choice, f_dir_prog)
        
        self._get_symm()
        self._calc_rotation_matrix_anb_b()
        
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
spgr_given_name is number or name of the space groupe
spgr_choice is choise of origin, 1, 2, "abc", "bac"
f_dir_prog is directory where the file "itables.txt" it is 

centr is inversion center
p_centr is position of inversin center
el_symm is element of symmetry
orig is packing
spgr_name is name of space groupe
spgr_number is number of space groupe

r_11, r_12, r_13  
r_21, r_22, r_23    element of symmetry in form of element of rotation matrix
r_31, r_32, r_33 

b_1, b_2,  b_3 is translation vecto for symmetry elements
        """
        print(lsout)
        
    def _read_el_cards(self, f_itables):
        """
        reading information about space grooupe from file fitables to list of cards ldcard
        Info in file fitables:
        
        1 P1               Triclinic
        choice: 1
        centr: false
        pcentr: 0, 0, 0
        symmetry: X,Y,Z
        
        2 P-1              Triclinic
        ...
        """
        fid = open(f_itables, "r")
        lcontent = fid.readlines()
        fid.close()
    
        lcontent = [hh.strip() for hh in lcontent if hh.strip() != ""]
        ldcard = []
        dcard = None
        for hh in lcontent:
            lhelp = hh.split()
            if lhelp[0].isdigit():
                if dcard != None:
                    ldcard.append(dcard)
                dcard = {"number":lhelp[0], "name": lhelp[1], "syngony": lhelp[2]}
            else:
                lhelp = hh.split(":")
                if (lhelp[0].strip() in dcard.keys()):
                    dcard[lhelp[0].strip()].append(lhelp[1].strip())
                else:
                    dcard[lhelp[0].strip()] = [lhelp[1].strip()]
        self._p_spgr_table = ldcard
        

    def _get_symm(self):
        """
        get symmetry from space group
        """
        
        spgr_choice = self._p_spgr_choice
        
        spgr_given_name = self._p_spgr_given_name
        

        if spgr_given_name.isdigit():
            spgr_n = spgr_given_name
            spgr_name = ""
        else:
            spgr_n = ""
            spgr_name = spgr_given_name
        
        spgr_table = self._p_spgr_table

        for dcard in spgr_table:
            if (((dcard["number"] == spgr_n)|(dcard["name"] == spgr_name))&(dcard["choice"][0] == spgr_choice)):
                flag = True
                break
        if (not flag):
            print("Space groupe is not found")
            return
        
        flag = False
            
        lelsymm = []
        for ssymm in dcard["symmetry"]:
            lelsymm.append(self._trans_str_to_el_symm(ssymm))
        centr = dcard["centr"][0]=="true"
        pcentr = [float(hh) for hh in dcard["pcentr"][0].split(",")]
        fletter = dcard["name"][0]
        spgr = dcard["name"]
        number = int(dcard["number"])
        if (fletter == "P"):
            lorig = [(0, 0, 0)]
        elif fletter == "C":
            lorig = [(0, 0, 0), (0.5, 0.5, 0)]
        elif fletter == "I":
            lorig = [(0, 0, 0), (0.5, 0.5, 0.5)]
        elif fletter == "F":
            lorig = [(0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5)]
        elif (fletter == "R"):
            if spgr_choice == "1":
                lorig = [(0, 0, 0), (0.66667, 0.33333, 0.33333), (0.33334, 0.66666, 0.66666)]
            else:
                lorig = [(0, 0, 0)]
        else:
            print("Undefined syngony")

            
        self._p_centr = centr
        self._p_el_symm = lelsymm
        self._p_orig = lorig
        self._p_p_centr = pcentr
        self._p_spgr_name = spgr
        self._p_spgr_number = number
        
    def _calc_rotation_matrix_anb_b(self):
        """
        give representation for rotation matrix: r_11, r_22, r_33, r_12, r_13, r_23 and vector b_1, b_2, b_3
        """
        lel_symm = self._p_el_symm
        b_1 = numpy.array([hh[0] for hh in lel_symm], dtype = float)
        r_11 = numpy.array([hh[1] for hh in lel_symm], dtype = int)
        r_12 = numpy.array([hh[2] for hh in lel_symm], dtype = int)
        r_13 = numpy.array([hh[3] for hh in lel_symm], dtype = int)

        b_2 = numpy.array([hh[4] for hh in lel_symm], dtype = float)
        r_21 = numpy.array([hh[5] for hh in lel_symm], dtype = int)
        r_22 = numpy.array([hh[6] for hh in lel_symm], dtype = int)
        r_23 = numpy.array([hh[7] for hh in lel_symm], dtype = int)

        b_3 = numpy.array([hh[8] for hh in lel_symm], dtype = float)
        r_31 = numpy.array([hh[9] for hh in lel_symm], dtype = int)
        r_32 = numpy.array([hh[10] for hh in lel_symm], dtype = int)
        r_33 = numpy.array([hh[11] for hh in lel_symm], dtype = int)
        
        self._p_r_11 = r_11
        self._p_r_12 = r_12
        self._p_r_13 = r_13

        self._p_r_21 = r_21
        self._p_r_22 = r_22
        self._p_r_23 = r_23

        self._p_r_31 = r_31
        self._p_r_32 = r_32
        self._p_r_33 = r_33

        self._p_b_1 = b_1
        self._p_b_2 = b_2
        self._p_b_3 = b_3

        
        

    def calc_hkl_equiv(self, h, k, l):
        """
        give equivalent reflections of hkl and its multiplicity
        """
        r_11 = self._p_r_11
        r_12 = self._p_r_12
        r_13 = self._p_r_13
        r_21 = self._p_r_21
        r_22 = self._p_r_22
        r_23 = self._p_r_23
        r_31 = self._p_r_31
        r_32 = self._p_r_32
        r_33 = self._p_r_33

        h_s = r_11*h + r_21*k + r_31*l 
        k_s = r_12*h + r_22*k + r_32*l 
        l_s = r_13*h + r_23*k + r_33*l 
        
        hkl_s = numpy.vstack([h_s, k_s, l_s])
        hkl_s = numpy.hstack([hkl_s,-1*hkl_s])
        hkl_s_un = numpy.unique(hkl_s, axis=1)
        multiplicity = int(round(hkl_s.shape[1]*1./hkl_s_un.shape[1]))
        h_s, k_s, l_s = hkl_s_un[0, :], hkl_s_un[1, :], hkl_s_un[2, :]
        return h_s, k_s, l_s, multiplicity

    def calc_xyz_mult(self, x, y, z):
        """
        give unique x,y,z elements and calculate multiplicit for given x,y,z fract
        """
        r_11 = self._p_r_11
        r_12 = self._p_r_12
        r_13 = self._p_r_13
        r_21 = self._p_r_21
        r_22 = self._p_r_22
        r_23 = self._p_r_23
        r_31 = self._p_r_31
        r_32 = self._p_r_32
        r_33 = self._p_r_33
        b_1 = self._p_b_1
        b_2 = self._p_b_2
        b_3 = self._p_b_3
        
        lorig = self._p_orig
        centr = self._p_centr
        p_centr = self._p_p_centr

        x,y,z=0.125,0.125,0.125

        x_s = numpy.round(numpy.mod(r_11*x + r_12*y + r_13*z + b_1, 1), 5)
        y_s = numpy.round(numpy.mod(r_21*x + r_22*y + r_23*z + b_2, 1), 5)
        z_s = numpy.round(numpy.mod(r_31*x + r_32*y + r_33*z + b_3, 1), 5)

        x_o = [orig[0] for orig in lorig]
        y_o = [orig[1] for orig in lorig]
        z_o = [orig[2] for orig in lorig]
        
        x_s_2d, x_o_2d = numpy.meshgrid(x_s, x_o)
        y_s_2d, y_o_2d = numpy.meshgrid(y_s, y_o)
        z_s_2d, z_o_2d = numpy.meshgrid(z_s, z_o)
        
        x_s_2d = numpy.round(numpy.mod(x_s_2d+x_o_2d, 1), 5)
        y_s_2d = numpy.round(numpy.mod(y_s_2d+y_o_2d, 1), 5)
        z_s_2d = numpy.round(numpy.mod(z_s_2d+z_o_2d, 1), 5)

        x_s = x_s_2d.flatten()
        y_s = y_s_2d.flatten()
        z_s = z_s_2d.flatten()

        if centr:
            x_s_h = numpy.round(numpy.mod(2.*p_centr[0]-1.*x_s, 1), 5)
            y_s_h = numpy.round(numpy.mod(2.*p_centr[1]-1.*y_s, 1), 5)
            z_s_h = numpy.round(numpy.mod(2.*p_centr[2]-1.*z_s, 1), 5)
            x_s =numpy.hstack([x_s, x_s_h])
            y_s =numpy.hstack([y_s, y_s_h])
            z_s =numpy.hstack([z_s, z_s_h])
                        
        xyz_s = numpy.vstack([x_s, y_s, z_s])
        
        xyz_s_un = numpy.unique(xyz_s, axis=1)
        multiplicity = int(round(xyz_s.shape[1]*1./xyz_s_un.shape[1]))
        x_s, y_s, z_s = xyz_s_un[0, :], xyz_s_un[1, :], xyz_s_un[2, :]
        return x_s, y_s, z_s, multiplicity

    
    def _trans_str_to_el_symm(self, str1):
        """
        transform string to element of symmetry: (x,y,-z) -> 0.0 1 0 0  0.0 0 1 0  0.0 0 0 -1
        """
        str2="".join(str1.split(" "))
        lhelp1,lhelp2,lhelp3=[],[],[]
        lhelp1=[hh for hh in str2.split('(') if hh!=""]
        [lhelp2.extend(hh.split(')')) for hh in lhelp1 if hh!=""]
        [lhelp3.extend(hh.split(',')) for hh in lhelp2 if hh!=""]
        lAx=['x','y','z']
        lelsymm=[]
        for hh in lhelp3:
            elsymmh=[0.0,0,0,0]
            strh=hh
            for inum,Ax in enumerate(lAx):
                if (strh.find(Ax)!=-1):
                    if (strh.find("+"+Ax)!=-1):
                        elsymmh[inum+1]=1
                        strh="".join(strh.split("+"+Ax))
                    elif (strh.find("-"+Ax)!=-1):
                        elsymmh[inum+1]=-1
                        strh="".join(strh.split("-"+Ax))
                    else:
                        elsymmh[inum+1]=1
                        strh="".join(strh.split(Ax))
            if (strh==""):
                pass
            elif (strh.find("/")!=-1):
                lhelp1=strh.split("/")
                elsymmh[0]=float(lhelp1[0])/float(lhelp1[1])
            else:
                elsymmh[0]=float(strh)
            lelsymm.append(elsymmh)
        elsymm=[]
        [elsymm.extend(hh) for hh in lelsymm]
        return elsymm
    



class Fract(dict):
    """
    Fract of atom_site(s) in unit cell.
    """
    def __init__(self, x = 0., y = 0., z = 0.):
        super(Fract, self).__init__()
        dd= {"x": x, "y": y, "z": z}
        self.update(dd)


    def __repr__(self):
        lsout = """Fract: \n xyz: {:} {:} {:}""".format(
                self["x"], self["y"], self["z"])
        return lsout


    def _calc_phase(self, h, k, l, space_groupe, x = None, y = None, 
                    z = None):
        """
        calculate phase: exp(-2 pi i * (h*x+k*y+l*z))
        r_11, r_22, r_33, r_12, r_13, r_23 are element of symmetry 
        """
        if x != None:
            self["x"] = x
        if y != None:
            self["y"] = y
        if z != None:
            self["z"] = z
        x, y, z = self["x"], self["y"], self["z"]

        r_11, r_12 = space_groupe["r_11"], space_groupe["r_12"]
        r_13, r_21 = space_groupe["r_13"], space_groupe["r_21"]
        r_22, r_23 = space_groupe["r_22"], space_groupe["r_23"]
        r_31, r_32 = space_groupe["r_31"], space_groupe["r_32"]
        r_33 = space_groupe["r_33"]
        
        np_h, np_x, np_r_11 = numpy.meshgrid(h, x, r_11, indexing="ij")
        np_k, np_y, np_r_22 = numpy.meshgrid(k, y, r_22, indexing="ij")
        np_l, np_z, np_r_33 = numpy.meshgrid(l, z, r_33, indexing="ij")
        np_r_12 = numpy.meshgrid(h, x, r_12, indexing="ij")[2]
        np_r_13 = numpy.meshgrid(k, y, r_13, indexing="ij")[2]
        np_r_23 = numpy.meshgrid(l, z, r_23, indexing="ij")[2]
        np_r_21 = numpy.meshgrid(h, x, r_21, indexing="ij")[2]
        np_r_31 = numpy.meshgrid(k, y, r_31, indexing="ij")[2]
        np_r_32 = numpy.meshgrid(l, z, r_32, indexing="ij")[2]
        
        np_h_s = np_h*np_r_11 + np_k*np_r_12 + np_l*np_r_13
        np_k_s = np_h*np_r_21 + np_k*np_r_22 + np_l*np_r_23
        np_l_s = np_h*np_r_31 + np_k*np_r_32 + np_l*np_r_33
        
        phase = numpy.exp(2*numpy.pi*1j*(np_x*np_h_s + np_y*np_k_s+ np_z*np_l_s))
        
        d_out = dict(phase = phase)
        self.update(d_out)
        
    def _calc_multiplicity(self, r_11, r_12, r_13, r_21, r_22, r_23, r_31, 
                    r_32, r_33, b_1, b_2, b_3, x = None, y = None, z = None):
        """
        calculate atom multiplicity
        """
        if x != None:
            self["x"] = x
        if y != None:
            self["y"] = y
        if z != None:
            self["z"] = z
        x, y, z = self["x"], self["y"], self["z"]
        b_1, b_2, b_3 = space_groupe["b_1"], space_groupe["b_2"], space_groupe["b_3"]

        r_11, r_12 = space_groupe["r_11"], space_groupe["r_12"]
        r_13, r_21 = space_groupe["r_13"], space_groupe["r_21"]
        r_22, r_23 = space_groupe["r_22"], space_groupe["r_23"]
        r_31, r_32 = space_groupe["r_31"], space_groupe["r_32"]
        r_33 = space_groupe["r_33"]

        lorig = space_groupe["orig"]
        centr = space_groupe["centr"]
        pcentr = space_groupe["pcentr"]
        
        

        
    def els4pos(self, space_groupe):
        """
        give the lelements of symmetry which transfer atom to the same atom
        """
        
        lelsymm = space_groupe["elsymm"]
        lorig = space_groupe["orig"]
        centr = space_groupe["centr"]
        pcentr = space_groupe["pcentr"]
    
        lelsat = []
        lelsuniqat, lcoorduniqat = [], []
        [x, y, z] = xyz
        x, y, z = x%1, y%1, z%1
        for els in lelsymm:
            for orig in lorig:
                xat = (els[0] + els[1]*x + els[ 2]*y + els[ 3]*z+orig[0])%1
                yat = (els[4] + els[5]*x + els[ 6]*y + els[ 7]*z+orig[1])%1
                zat = (els[8] + els[9]*x + els[10]*y + els[11]*z+orig[2])%1
                elsn = [els[0]+orig[0],els[1],els[2],els[3],els[4]+orig[1],els[5],els[6],els[7],
                els[8]+orig[2],els[9],els[10],els[11]]
                if ((abs(xat-x)<10**-5)&(abs(yat-y)<10**-5)&(abs(zat-z)<10**-5)): lelsat.append(elsn)
                xyzatu = (round(xat,4),round(yat,4),round(zat,4))
                if (not(xyzatu in lcoorduniqat)):
                    lcoorduniqat.append(xyzatu)
                    lelsuniqat.append(elsn)
                if (centr):
                    elsn=[2*pcentr[0]-els[0]-orig[0],-1*els[1],-1*els[2],-1*els[3],
                        2*pcentr[1]-els[4]-orig[1],-1*els[5],-1*els[6],-1*els[7],
                        2*pcentr[2]-els[8]-orig[2],-1*els[9],-1*els[10],-1*els[11]]
                    xat,yat,zat=(2*pcentr[0]-xat)%1,(2*pcentr[1]-yat)%1,(2*pcentr[2]-zat)%1
                    if ((abs(xat-x)<10**-5)&(abs(yat-y)<10**-5)&(abs(zat-z)<10**-5)): lelsat.append(elsn)
                    xyzatu=(round(xat,4),round(yat,4),round(zat,4))
                    if (not(xyzatu in lcoorduniqat)):
                        lcoorduniqat.append(xyzatu)
                        lelsuniqat.append(elsn)
        return lelsat,lelsuniqat
    
    
    def calc_model(self, h, k, l, space_groupe, x = None, y = None, z = None):
        
        if x != None:
            self["x"] = x
        if y != None:
            self["y"] = y
        if z != None:
            self["z"] = z
        x, y, z = self["x"], self["y"], self["z"]

        
        self._calc_phase(self, h, k, l, space_groupe)
        self._calc_multiplicity(self, h, k, l, space_groupe)
        
        d_out = dict(h=h, k=k, l=l, space_groupe=space_groupe)
        self.update(d_out)
        
        
    def set_vals(self, d_vals, refresh = False):
        """
        Set values 
        """
        keys = d_vals.keys()
        llab = ["x", "y", "z"]
        llab_in = [(hh in keys) for hh in llab]
        for lab, cond_h in zip(llab, llab_in):
            if cond_h:
                self[lab] = d_vals[lab]
            
        if refresh:
            h = self["h"] 
            k = self["k"] 
            l = self["l"] 
            
            r_11 = self["r_11"] 
            r_12 = self["r_12"] 
            r_13 = self["r_13"] 
            r_21 = self["r_21"] 
            r_22 = self["r_22"] 
            r_23 = self["r_23"] 
            r_31 = self["r_31"] 
            r_32 = self["r_32"] 
            r_33 = self["r_33"]
            
            self.calc_model(h, k, l, r_11, r_12, r_13, r_21, r_22, r_23, r_31, 
                    r_32, r_33)

            
    def soft_copy(self):
        """
        Soft copy of the object with saving the links on the internal parameter of the object
        """
        obj_new = Fract(x=self["x"], y=self["y"], z=self["z"])
        llab = ["h", "k", "l", "r_11", "r_22", "r_33", "r_12", "r_13", "r_23", 
                "phase"]
        keys = self.keys()
        llab_in = [(hh in keys) for hh in llab]
        for lab, cond_h in zip(llab, llab_in):
            if cond_h:
                obj_new[lab] = self[lab]
        return obj_new


    
class ADP(dict):
    """
    ADP
    """
    def __init__(self, beta_11 = 0., beta_22 = 0., beta_33 = 0., 
                 beta_12 = 0., beta_13 = 0., beta_23 = 0., b_iso = 0.):
        super(ADP, self).__init__()
        dd= {"beta_11": beta_11, "beta_22": beta_22, "beta_33": beta_33,
             "beta_12": beta_12, "beta_13": beta_13, "beta_23": beta_23,
             "b_iso": b_iso}
        self.update(dd)

    def __repr__(self):
        lsout = """Debye Waller: \n beta_11: {:}\n beta_22: {:}\n beta_33: {:}
 beta_12: {:}\n beta_13: {:}\n beta_23: {:}""".format(
 self["beta_11"], self["beta_22"], self["beta_33"], self["beta_12"], 
 self["beta_13"], self["beta_23"])
        return lsout


    def _calc_dwf_iso(self, sthovl):
        """
        isotropic harmonic Debye-Waller factor
        """
        b_iso = self["b_iso"]
        sthovl_sq = sthovl**2
        np_biso, np_sthovl_sq = numpy.meshgrid(sthovl_sq, b_iso, indexing="ij")
        
        dwf_iso = numpy.exp(-np_biso*np_sthovl_sq)
        d_out = dict(dwf_iso = dwf_iso)#2 dimensional
        self.update(d_out)


    def calc_dwf_aniso(self, h, k, l, space_groupe = None):
        """
        anisotropic harmonic Debye-Waller factor
        
        h,k,l is 1D (temporary solution)
        """
        beta_11, beta_22 = self["beta_11"], self["beta_22"] 
        beta_33, beta_12 = self["beta_33"], self["beta_12"]
        beta_13, beta_23 = self["beta_13"], self["beta_23"]
        

        np_h, np_beta_11 = numpy.meshgrid(h, beta_11, indexing="ij")
        np_k, np_beta_22 = numpy.meshgrid(k, beta_22, indexing="ij")
        np_l, np_beta_33 = numpy.meshgrid(l, beta_33, indexing="ij")
        np_h, np_beta_12 = numpy.meshgrid(h, beta_12, indexing="ij")
        np_h, np_beta_13 = numpy.meshgrid(h, beta_13, indexing="ij")
        np_h, np_beta_23 = numpy.meshgrid(h, beta_23, indexing="ij")
        
        dwf_aniso = numpy.exp(-1*(np_beta_11*np_h**2 + np_beta_22*np_k**2 + 
                           np_beta_33*np_l**2 + 2.*np_beta_12*np_h*np_k + 
                    2.*np_beta_13*np_h*np_l + 2.*np_beta_23*np_k*np_l))
        
        return dwf_aniso 

        
    def calc_model(self, h, k, l, sthovl, r_11, r_12, r_13, r_21, r_22, r_23, 
                        r_31, r_32, r_33,
                   beta_11 = None, beta_22 = None, beta_33 = None, 
                   beta_12 = None, beta_13 = None, beta_23 = None, b_iso = None):
        
        if beta_11 != None:
            self["beta_11"] = beta_11 
        if beta_22 != None:
            self["beta_22"] = beta_22 
        if beta_33 != None:
            self["beta_33"] = beta_33
        if beta_12 != None:
            self["beta_12"] = beta_12 
        if beta_13 != None:
            self["beta_13"] = beta_13 
        if beta_23 != None:
            self["beta_23"] = beta_23 
        
        self._calc_dwf_iso(self, sthovl)
        self._calc_dwf_aniso(self, h, k, l)
        
        d_out = dict(h=h, k=k, l=l, r_11=r_11, r_22=r_22, r_33=r_33, r_12=r_12, 
                     r_13=r_13, r_23=r_23, r_21=r_21, r_31=r_31, r_32=r_32, 
                     sthovl=sthovl)
        self.update(d_out)
        
        
    def set_vals(self, d_vals, refresh = False):
        """
        Set values 
        """
        keys = d_vals.keys()
        llab = ["beta_11", "beta_22", "beta_33", "beta_12", "beta_13", 
                "beta_23"]
        llab_in = [(hh in keys) for hh in llab]
        for lab, cond_h in zip(llab, llab_in):
            if cond_h:
                self[lab] = d_vals[lab]
            
        if refresh:
            h = self["h"] 
            k = self["k"] 
            l = self["l"] 
            sthovl = self["sthovl"]
            r_11 = self["r_11"]
            r_22 = self["r_22"]
            r_33 = self["r_33"]
            r_12 = self["r_12"]
            r_13 = self["r_13"]
            r_23 = self["r_23"]
            r_21 = self["r_21"]
            r_31 = self["r_31"]
            r_32 = self["r_32"]
            
            self.calc_model(h, k, l, sthovl, r_11, r_12, r_13, r_21, r_22, 
                            r_23, r_31, r_32, r_33)

            
    def soft_copy(self):
        """
        Soft copy of the object with saving the links on the internal parameter of the object
        """
        obj_new = ADP(beta_11 = self["beta_11"], 
                              beta_22 = self["beta_22"], 
                   beta_33 = self["beta_33"], beta_12 = self["beta_12"], 
                   beta_13 = self["beta_13"], beta_23 = self["beta_23"], 
                   b_iso = self["b_iso"])
        llab = ["dwf_iso", "dwf_aniso", "h", "k", "l", "sthovl", "r_11", 
                "r_22", "r_33", "r_12", "r_13", "r_23", "r_21", "r_31", "r_32"]
        keys = self.keys()
        llab_in = [(hh in keys) for hh in llab]
        for lab, cond_h in zip(llab, llab_in):
            if cond_h:
                obj_new[lab] = self[lab]
        return obj_new



class Magnetism(dict):
    """
    Magnetism
    """
    def __init__(self, kappa = 1.0, factor_lande = 2.0, type = "",
                 chi_11 = 0., chi_22 = 0., chi_33 = 0., 
                 chi_12 = 0., chi_13 = 0., chi_23 = 0.):
        super(Magnetism, self).__init__()
        dd= {"kappa": kappa, "factor_lande": factor_lande , 
             "chi_11": chi_11, "chi_22": chi_22, "chi_33": chi_33,
             "chi_12": chi_12, "chi_13": chi_13, "chi_23": chi_23}
        self.update(dd)
        
    def __repr__(self):
        lsout = """Magnetism: \n chi_11: {:}\n chi_22: {:}\n chi_33: {:}
 chi_12: {:}\n chi_13: {:}\n chi_23: {:}\n kappa: {:}
 factor_lande: {:}""".format(
 self["chi_11"], self["chi_22"], self["chi_33"], self["chi_12"], 
 self["chi_13"], self["chi_23"], self["kappa"], self["factor_lande"])
        return lsout
    
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
    
    def calc_model(self, ia, ib, ic, matrix_ib, chi_11=None, chi_22=None, 
                   chi_33=None, chi_12=None, chi_13=None, chi_23=None):
        
        if chi_11 != None:
            self["chi_11"] = chi_11 
        if chi_22 != None:
            self["chi_22"] = chi_22 
        if chi_33 != None:
            self["chi_33"] = chi_33
        if chi_12 != None:
            self["chi_12"] = chi_12 
        if chi_13 != None:
            self["chi_13"] = chi_13 
        if chi_23 != None:
            self["chi_23"] = chi_23 
        
        self._calc_chi_loc(ia, ib, ic, matrix_ib)
        d_out = dict(ia=ia, ib=ib, ic=ic, matrix_ib=matrix_ib)
        self.update(d_out)

    
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

    
    def set_vals(self, d_vals, refresh = False):
        """
        Set values 
        """
        keys = d_vals.keys()
        llab = ["chi_11", "chi_22", "chi_33", "chi_12", "chi_13", "chi_23"]
                
        llab_in = [(hh in keys) for hh in llab]
        for lab, cond_h in zip(llab, llab_in):
            if cond_h:
                self[lab] = d_vals[lab]
        self._constr_singony()
            
        if refresh:
            ia = self["ia"] 
            ib = self["ib"] 
            ic = self["ic"] 
            matrix_ib = self["matrix_ib"]
            self.calc_model(ia, ib, ic, matrix_ib)

            
    def soft_copy(self):
        """
        Soft copy of the object with saving the links on the internal parameter of the object
        """
        obj_new = Magnetism(chi_11=self["chi_11"], chi_22=self["chi_22"], 
                            chi_33=self["chi_33"], chi_12=self["chi_12"], 
                            chi_13=self["chi_13"], chi_23=self["chi_23"])
        llab = ["matrix_chi_loc", "ia", "ib", "ic", "matrix_ib"]
        keys = self.keys()
        llab_in = [(hh in keys) for hh in llab]
        for lab, cond_h in zip(llab, llab_in):
            if cond_h:
                obj_new[lab] = self[lab]
        return obj_new



    
    
class AtomSite(dict):
    """
    AtomSite
    """
    def __init__(self, atom_name="", atom_nucl_type="", b_scat=0., 
                 occupation =1.,
                 fract=Fract(), adp=ADP(),
                 magnetism=Magnetism()):
        super(AtomSite, self).__init__()
        dd= {"atom_name": atom_name, "atom_nucl_type": atom_nucl_type, 
             "b_scat": b_scat, "occupation": occupation, "fract": fract, 
             "adp": adp, "magnetism": magnetism,  }
        self.update(dd)
        
    def __repr__(self):
        lsout = """AtomSite: \n name: {:}\n nucl_type: {:}, b_scat: {} 
 occupation {:}\n fract {:}
 adp: {}\n magnetism: {}""".format(self["atom_name"], 
 self["atom_nucl_type"],  self["b_scat"],  self["occupation"], self["fract"], 
 self["adp"], self["magnetism"])
        return lsout
    
    def calc_fn(self, h, k, l, space_groupe, cell, fract=None, adp=None):
        """
        calculate nuclear structure factor
        """
        if adp != None:
            self["adp"] = adp
        if fract != None:
            self["fract"] = fract

        fract.calc_mult
        space_groupe.calc_xyz_mult()

        multiplicity = atom_site["multiplicity"] 
        occupation = atom_site["occupation"]
        b_scat = atom_site["b_scat"]
        
        atom_site["fract"].calc_model(h, k, l, r_11, r_12, r_13, r_21, r_22, 
            r_23, r_31, r_32, r_33)
        atom_site["adp"].calc_model(h, k, l, sthovl, r_11, r_12, r_13, 
            r_21, r_22, r_23, r_31, r_32, r_33)
        
        phase = atom_site["fract"]["phase"]
        
        dwf_aniso = atom_site["adp"]["dwf_aniso"]
        
        hh_1 = phase*dwf_aniso
        hh_2 = hh_1.sum(axes = 2)
        
        h_1 = multiplicity*occupation*b_scat
        
        hh_1 = numpy.meshgrid(h, h_1, indexing = "ij")[1]
        
        # b_scat * occ * mult * sum_el.symm.(phase)
        hh_3 = hh_1*hh_2
        f_hkl_as = hh_3.sum(axes=1)*1./len(lelsymm)#1 dimensional array
        
        np_h, np_orig_x = numpy.meshgrid(h, orig_x, indexing = "ij")
        np_k, np_orig_y = numpy.meshgrid(k, orig_y, indexing = "ij")
        np_l, np_orig_z = numpy.meshgrid(l, orig_z, indexing = "ij")
        
        
        np_f_hkl_as = numpy.exp(2*numpy.pi*1j*(np_h*np_orig_x+np_k*np_orig_y+np_l*np_orig_z))
        f_hkl_as = np_f_hkl_as.sum(axes=1)*1./len(lorig)

        if (centr):
            orig = crystal_symmetry["p_centr"]
            f_nucl = (f_hkl_as+f_hkl_as.conjugate()*numpy.exp(2.*2.*numpy.pi*1j* (h*orig[0]+k*orig[1]+l*orig[2])))*0.5
        else:
            f_nucl = Fhklas
    
    
    
class Crystal(dict):
    """
    Crystal
    """
    def __init__(self, crystal_symmetry = SpaceGroupe(), 
                 cell = Cell(), atom_site = AtomSite()):
        super(Crystal, self).__init__()
        dd= {"crystal_symmetry": crystal_symmetry, "cell": cell, 
             "atom_site": AtomSite}
        self.update(dd)
        
    def __repr__(self):
        lsout = """Phase: \n crystal symmetry: {:}\n cell {:}
 atom_site {:}""".format(self["crystal_symmetry"], self["cell"], self["atom_site"])
        return lsout

    
    def calc_fn(self, h, k, l, crystal_symmetry = None, cell = None, 
                atom_site = None):
        """
        calculate nuclear structure factor
        """
        if crystal_symmetry != None:
            self["crystal_symmetry"] = crystal_symmetry
        if cell != None:
            self["cell"] = cell 
        if atom_site != None:
            self["atom_site"] = atom_site 

        crystal_symmetry = self["crystal_symmetry"]
        cell = self["cell"]
        atom_site = self["atom_site"] 
        
        sthovl = cell.calc_sthovl(h, k, l)
        
        crystal_symmetry.calc_model()

        lelsymm = crystal_symmetry["elsymm"]
        lorig = crystal_symmetry["orig"]
        centr = crystal_symmetry["centr"]
        
        r_11 = crystal_symmetry["r_11"]
        r_12 = crystal_symmetry["r_12"]
        r_13 = crystal_symmetry["r_13"]

        r_21 = crystal_symmetry["r_21"]
        r_22 = crystal_symmetry["r_22"]
        r_23 = crystal_symmetry["r_23"]

        r_31 = crystal_symmetry["r_31"]
        r_32 = crystal_symmetry["r_32"]
        r_33 = crystal_symmetry["r_33"]
        
        multiplicity = atom_site["multiplicity"] 
        occupation = atom_site["occupation"]
        b_scat = atom_site["b_scat"]
        
        atom_site["fract"].calc_model(h, k, l, r_11, r_12, r_13, r_21, r_22, 
            r_23, r_31, r_32, r_33)
        atom_site["adp"].calc_model(h, k, l, sthovl, r_11, r_12, r_13, 
            r_21, r_22, r_23, r_31, r_32, r_33)
        
        phase = atom_site["fract"]["phase"]
        
        dwf_aniso = atom_site["adp"]["dwf_aniso"]
        
        hh_1 = phase*dwf_aniso
        hh_2 = hh_1.sum(axes = 2)
        
        h_1 = multiplicity*occupation*b_scat
        
        hh_1 = numpy.meshgrid(h, h_1, indexing = "ij")[1]
        
        # b_scat * occ * mult * sum_el.symm.(phase)
        hh_3 = hh_1*hh_2
        f_hkl_as = hh_3.sum(axes=1)*1./len(lelsymm)#1 dimensional array
        
        np_h, np_orig_x = numpy.meshgrid(h, orig_x, indexing = "ij")
        np_k, np_orig_y = numpy.meshgrid(k, orig_y, indexing = "ij")
        np_l, np_orig_z = numpy.meshgrid(l, orig_z, indexing = "ij")
        
        
        np_f_hkl_as = numpy.exp(2*numpy.pi*1j*(np_h*np_orig_x+np_k*np_orig_y+np_l*np_orig_z))
        f_hkl_as = np_f_hkl_as.sum(axes=1)*1./len(lorig)

        if (centr):
            orig = crystal_symmetry["p_centr"]
            f_nucl = (f_hkl_as+f_hkl_as.conjugate()*numpy.exp(2.*2.*numpy.pi*1j* (h*orig[0]+k*orig[1]+l*orig[2])))*0.5
        else:
            f_nucl = Fhklas
        return f_nucl


    def calc_sft(self):
        """
        calculate structure factor tensor
        """
        pass
    
    def calc_model(self, reflection, crystal_symmetry = None, cell = None,
                   atom_site = None):
        if crystal_symmetry != None:
            self["crystal_symmetry"] = crystal_symmetry 
        if cell != None:
            self["cell"] = cell 
        if atom_site != None:
            self["atom_site"] = atom_site 
            
        np_fn_1d = self.calc_fn(h, k, l)#complex 1D numpy array of nuclear structure factors over list of hkl
        np_sft_2d = None #complex 2D numpy array of tensor structure factors /11, 22, 33, 12, 13, 23/... over list of hkl 
        d_out = dict(fn = np_fn_1d, sft = np_sft_2d)
        self.update(d_out)
        
        
    def set_vals(self, d_vals, refresh = False):
        """
        Set values 
        """
        keys = d_vals.keys()
        llab = ["crystal_symmetry", "cell", "atom_site"]
        llab_in = [(hh in keys) for hh in llab]
        for lab, cond_h in zip(llab, llab_in):
            if cond_h:
                self[lab] = d_vals[lab]
                
        #redo it
        llab = ["resolution", "u", "v", "w", "x", "y"]
        llab_in = [(hh in keys) for hh in llab]
        cond = any(llab_in)
        if cond:
            d_val_as = {}
            for lab, cond_h in zip(llab, llab_in):
                if cond_h:
                    d_val_as[lab] = d_vals[lab]
            self["crystal_symmetry"].set_vals(d_val_as, refresh)

        #redo it
        llab = ["p1", "p2", "p3", "p4"]
        llab_in = [(hh in keys) for hh in llab]
        cond = any(llab_in)
        if cond:
            d_val_as = {}
            for lab, cond_h in zip(llab, llab_in):
                if cond_h:
                    d_val_as[lab] = d_vals[lab]
            self["cell"].set_vals(d_val_as, refresh)

        #redo it
        llab = ["p1", "p2", "p3", "p4"]
        llab_in = [(hh in keys) for hh in llab]
        cond = any(llab_in)
        if cond:
            d_val_as = {}
            for lab, cond_h in zip(llab, llab_in):
                if cond_h:
                    d_val_as[lab] = d_vals[lab]
            self["atom_site"].set_vals(d_val_as, refresh)
            
        if refresh:
            tth = self["tth"]
            self.calc_model(tth)

            
    def soft_copy(self):
        """
        Soft copy of the object with saving the links on the internal parameter of the object
        """
        obj_new = Crystal(crystal_symmetry = self["crystal_symmetry"], 
                          cell = self["cell"], atom_site = self["atom_site"])
        llab = ["fn", "sft"]
        keys = self.keys()
        llab_in = [(hh in keys) for hh in llab]
        for lab, cond_h in zip(llab, llab_in):
            if cond_h:
                obj_new[lab] = self[lab]
        return obj_new


cell = Cell()
cell.keys()
cell.calc_model()
cell.keys()

h = numpy.array([1, 0, 0, 1], dtype=int)
k = numpy.array([0, 1, 0, 1], dtype=int)
l = numpy.array([0, 0, 1, 1], dtype=int)

x = numpy.array([0, 0.1, 0.2], dtype=int)
y = numpy.array([0, 0, 0.2], dtype=int)
z = numpy.array([0, 0, 0], dtype=int)



cell.calc_sthovl(h, k, l)
cell["sthovl"]

space_groupe = SpaceGroupe(spgr_given_name = "Fd-3m")
space_groupe.calc_model()
space_groupe.keys()
space_groupe["spgr_choice"]
cell.calc_model()
crystal = Crystal()             
#import Variable
#v_u = Variable.Variable(0.2)

        
if (__name__ == "__main__"):
  pass
