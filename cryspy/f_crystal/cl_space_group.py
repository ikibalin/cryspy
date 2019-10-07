"""
define classes to describe space group
"""
__author__ = 'ikibalin'
__version__ = "2019_08_27"
import os
import numpy
from pycifstar import Global

from cryspy.f_common.cl_fitable import Fitable

class SpaceGroup(object):
    """
    Space Group
    """
    def __init__(self, spgr_given_name = "P1", spgr_choice = "1",
                 f_dir_prog = os.path.dirname(__file__)):
        #         f_dir_prog = os.getcwd()):
        super(SpaceGroup, self).__init__()
        
        self.__spgr_given_name = None
        self.__spgr_choice = None
        self.__f_dir_prog = None
        self.__spgr_table = None


        self.__singony = "Triclinic"

        self.__centr = None
        self.__el_symm = None
        self.__orig = None
        self.__p_centr = None
        self.__spgr_name = None
        self.__spgr_number = None

        self.__r_11 = None
        self.__r_12 = None
        self.__r_13 = None
        self.__r_21 = None
        self.__r_22 = None
        self.__r_23 = None
        self.__r_31 = None
        self.__r_32 = None
        self.__r_33 = None

        self.__b_1 = None
        self.__b_2 = None
        self.__b_3 = None

        self.f_dir_prog = f_dir_prog
        self.spgr_name = spgr_given_name
        self.spgr_choice = spgr_choice


    @property
    def spgr_given_name(self):
        """

        reference:
        """
        return self.__spgr_given_name
    @spgr_given_name.setter
    def spgr_given_name(self, x):
        self.__spgr_given_name = "".join(x.split()).strip("\"").strip("'")
        if self.is_defined:
            _ = self._form_object()

    @property
    def spgr_choice(self):
        """

        reference:
        """
        return self.__spgr_choice
    @spgr_choice.setter
    def spgr_choice(self, x):
        if isinstance(x, float):
            x_in = "{:}".format(int(x))
        else:
            x_in = str(x)
        self.__spgr_choice = x_in.strip()
        if self.is_defined:
            _ = self._form_object()

    @property
    def f_dir_prog(self):
        """

        reference:
        """
        return self.__f_dir_prog
    @f_dir_prog.setter
    def f_dir_prog(self, x):
        f_itables = os.path.join(x, "tables", "itables.txt")
        self._read_el_cards(f_itables)    
        self.__f_dir_prog = x
        if self.is_defined:
            _ = self._form_object()



    @property
    def spgr_name(self):
        return self.__spgr_name
    @spgr_name.setter
    def spgr_name(self, x):
        self.__spgr_given_name = "".join(str(x).split()) # __spgr_given_name is not mistake
        if self.is_defined:
            _ = self._form_object()

    @property
    def spgr_number(self):
        return self.__spgr_number
    @property
    def singony(self):
        return self.__singony
    @property
    def bravais_lattice(self):
        return self.__singony
    @property
    def centr(self):
        return self.__centr
    @property
    def el_symm(self):
        return self.__el_symm
    @property
    def orig(self):
        return self.__orig
    @property
    def p_centr(self):
        return self.__p_centr
    @property
    def r_11(self):
        return self.__r_11
    @property
    def r_22(self):
        return self.__r_22
    @property
    def r_33(self):
        return self.__r_33
    @property
    def r_12(self):
        return self.__r_12
    @property
    def r_13(self):
        return self.__r_13
    @property
    def r_23(self):
        return self.__r_23
    @property
    def r_21(self):
        return self.__r_21
    @property
    def r_31(self):
        return self.__r_31
    @property
    def r_32(self):
        return self.__r_32
    @property
    def b_1(self):
        return self.__b_1
    @property
    def b_2(self):
        return self.__b_2
    @property
    def b_3(self):
        return self.__b_3



    def _show_message(self, s_out: str):
        print("***  Error ***")
        print(s_out)

    def __repr__(self):
        ls_out = ["SpaceGroup:\n given name: {:}\n choice: {:}".format(self.spgr_given_name, self.spgr_choice)]
        if self.spgr_name is not None:
            ls_out.append(" name: {:}".format(self.spgr_name))
        if self.spgr_number is not None:
            ls_out.append(" number: {:}".format(self.spgr_number))
        #ls_out.append(" {:}\n directory: '{:}'".format(self._trans_el_symm_to_str(), self.f_dir_prog))
        return "\n".join(ls_out)

    @property
    def is_defined(self):
        """
        Output: True if all started parameters are given
        """
        cond = any([self.spgr_given_name is None, 
                    self.spgr_choice is None, 
                    self.f_dir_prog is None])
        return not(cond)
        
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
                dcard = {"number":lhelp[0], "name": lhelp[1], "singony": lhelp[2]}
            else:
                lhelp = hh.split(":")
                if (lhelp[0].strip() in dcard.keys()):
                    dcard[lhelp[0].strip()].append(lhelp[1].strip())
                else:
                    dcard[lhelp[0].strip()] = [lhelp[1].strip()]
        ldcard.append(dcard)
        self._p_spgr_table = ldcard
        

    def _get_symm(self):
        """
        get symmetry from space group
        """
        
        spgr_choice = self.spgr_choice
        spgr_given_name = self.spgr_given_name
        

        if spgr_given_name.isdigit():
            spgr_n = spgr_given_name
            spgr_name = ""
        else:
            spgr_n = ""
            spgr_name = spgr_given_name
        
        spgr_table = self._p_spgr_table
        flag = False
        for dcard in spgr_table:
            if (((dcard["number"] == spgr_n)|(dcard["name"] == spgr_name))&(dcard["choice"][0] == spgr_choice)):
                flag = True
                break
        if (not flag):
            print("Space group is not found: {:} {:} {:}".format(spgr_n, spgr_name, spgr_choice))
            return
        
        flag = False
            
        l_el_symm = []
        for ssymm in dcard["symmetry"]:
            l_el_symm.append(self._trans_str_to_el_symm(ssymm))
        centr = dcard["centr"][0]=="true"
        pcentr = [float(hh) for hh in dcard["pcentr"][0].split(",")]
        fletter = dcard["name"][0]
        spgr = dcard["name"]
        number = int(dcard["number"])
        singony = dcard["singony"]
        if (fletter == "P"):
            l_orig = [(0, 0, 0)]
        elif fletter == "C":
            l_orig = [(0, 0, 0), (0.5, 0.5, 0)]
        elif fletter == "I":
            l_orig = [(0, 0, 0), (0.5, 0.5, 0.5)]
        elif fletter == "F":
            l_orig = [(0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5)]
        elif (fletter == "R"):
            if spgr_choice == "1":
                l_orig = [(0, 0, 0), (0.66667, 0.33333, 0.33333), (0.33334, 0.66666, 0.66666)]
            else:
                l_orig = [(0, 0, 0)]
        else:
            print("Undefined syngony")

            
        self.__centr = centr
        self.__el_symm = l_el_symm
        self.__orig = l_orig
        self.__p_centr = pcentr
        self.__spgr_name = spgr
        self.__spgr_number = number
        self.__singony = singony
        
    def _calc_rotation_matrix_anb_b(self):
        """
        give representation for rotation matrix: r_11, r_22, r_33, r_12, r_13, r_23 and vector b_1, b_2, b_3
        """
        lel_symm = self.el_symm
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
        
        self.__r_11 = r_11
        self.__r_12 = r_12
        self.__r_13 = r_13

        self.__r_21 = r_21
        self.__r_22 = r_22
        self.__r_23 = r_23

        self.__r_31 = r_31
        self.__r_32 = r_32
        self.__r_33 = r_33

        self.__b_1 = b_1
        self.__b_2 = b_2
        self.__b_3 = b_3

    def _form_object(self):
        if self.is_defined:
            self._get_symm()
            self._calc_rotation_matrix_anb_b()
        else:
            self._show_message("Object is not properly defined")
            return False
        return True

    def calc_hkl_equiv(self, h, k, l):
        """
        give equivalent reflections of hkl and its multiplicity
        """

        r_11 = self.r_11
        r_12 = self.r_12
        r_13 = self.r_13
        r_21 = self.r_21
        r_22 = self.r_22
        r_23 = self.r_23
        r_31 = self.r_31
        r_32 = self.r_32
        r_33 = self.r_33

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
        r_11 = self.r_11
        r_12 = self.r_12
        r_13 = self.r_13
        r_21 = self.r_21
        r_22 = self.r_22
        r_23 = self.r_23
        r_31 = self.r_31
        r_32 = self.r_32
        r_33 = self.r_33
        b_1 = self.b_1
        b_2 = self.b_2
        b_3 = self.b_3
        
        l_orig = self.orig
        centr = self.centr
        p_centr = self.p_centr

        x_s = numpy.round(numpy.mod(r_11*x + r_12*y + r_13*z + b_1, 1), 6)
        y_s = numpy.round(numpy.mod(r_21*x + r_22*y + r_23*z + b_2, 1), 6)
        z_s = numpy.round(numpy.mod(r_31*x + r_32*y + r_33*z + b_3, 1), 6)

        x_o = [orig[0] for orig in l_orig]
        y_o = [orig[1] for orig in l_orig]
        z_o = [orig[2] for orig in l_orig]
        
        x_s_2d, x_o_2d = numpy.meshgrid(x_s, x_o)
        y_s_2d, y_o_2d = numpy.meshgrid(y_s, y_o)
        z_s_2d, z_o_2d = numpy.meshgrid(z_s, z_o)
        
        x_s_2d = numpy.round(numpy.mod(x_s_2d+x_o_2d, 1), 6)
        y_s_2d = numpy.round(numpy.mod(y_s_2d+y_o_2d, 1), 6)
        z_s_2d = numpy.round(numpy.mod(z_s_2d+z_o_2d, 1), 6)

        x_s = x_s_2d.flatten()
        y_s = y_s_2d.flatten()
        z_s = z_s_2d.flatten()

        if centr:
            x_s_h = numpy.round(numpy.mod(2.*p_centr[0]-1.*x_s, 1), 6)
            y_s_h = numpy.round(numpy.mod(2.*p_centr[1]-1.*y_s, 1), 6)
            z_s_h = numpy.round(numpy.mod(2.*p_centr[2]-1.*z_s, 1), 6)
            x_s =numpy.hstack([x_s, x_s_h])
            y_s =numpy.hstack([y_s, y_s_h])
            z_s =numpy.hstack([z_s, z_s_h])
                        
        xyz_s = numpy.vstack([x_s, y_s, z_s])

        xyz_s_un = numpy.unique(xyz_s, axis=1)
        n_atom = int(round(xyz_s.shape[1]*1./xyz_s_un.shape[1]))
        x_s, y_s, z_s = xyz_s_un[0, :], xyz_s_un[1, :], xyz_s_un[2, :]
        
        return x_s, y_s, z_s, n_atom
    
    def calc_el_symm_for_xyz(self, x_in, y_in, z_in):
        x ,y, z = x_in%1., y_in%1., z_in%1.
        r_11 = self.r_11
        r_12 = self.r_12
        r_13 = self.r_13
        r_21 = self.r_21
        r_22 = self.r_22
        r_23 = self.r_23
        r_31 = self.r_31
        r_32 = self.r_32
        r_33 = self.r_33
        b_1 = self.b_1
        b_2 = self.b_2
        b_3 = self.b_3
        
        l_orig = self.orig
        centr = self.centr
        p_centr = self.p_centr


        x_o = [orig[0] for orig in l_orig]
        y_o = [orig[1] for orig in l_orig]
        z_o = [orig[2] for orig in l_orig]
        
        r_11_2d, x_o_2d = numpy.meshgrid(r_11, x_o, indexing="ij")
        r_12_2d, y_o_2d = numpy.meshgrid(r_12, y_o, indexing="ij")
        r_13_2d, z_o_2d = numpy.meshgrid(r_13, z_o, indexing="ij")

        r_21_2d = numpy.meshgrid(r_21, x_o, indexing="ij")[0]
        r_22_2d = numpy.meshgrid(r_22, y_o, indexing="ij")[0]
        r_23_2d = numpy.meshgrid(r_23, z_o, indexing="ij")[0]
        
        r_31_2d = numpy.meshgrid(r_31, x_o, indexing="ij")[0]
        r_32_2d = numpy.meshgrid(r_32, y_o, indexing="ij")[0]
        r_33_2d = numpy.meshgrid(r_33, z_o, indexing="ij")[0]
        
        b_1_2d = numpy.meshgrid(b_1, x_o, indexing="ij")[0]
        b_2_2d = numpy.meshgrid(b_2, y_o, indexing="ij")[0]
        b_3_2d = numpy.meshgrid(b_3, z_o, indexing="ij")[0]
        
        b_1_2d = b_1_2d + x_o_2d
        b_2_2d = b_2_2d + x_o_2d
        b_3_2d = b_3_2d + x_o_2d

        e_11 = r_11_2d.flatten()
        e_12 = r_12_2d.flatten()
        e_13 = r_13_2d.flatten()

        e_21 = r_21_2d.flatten()
        e_22 = r_22_2d.flatten()
        e_23 = r_23_2d.flatten()

        e_31 = r_31_2d.flatten()
        e_32 = r_32_2d.flatten()
        e_33 = r_33_2d.flatten()

        e_1 = b_1_2d.flatten()
        e_2 = b_2_2d.flatten()
        e_3 = b_3_2d.flatten()

        if centr:
            me_11, me_12, me_13 = -1*e_11, -1*e_12, -1*e_13
            me_21, me_22, me_23 = -1*e_21, -1*e_22, -1*e_23
            me_31, me_32, me_33 = -1*e_31, -1*e_32, -1*e_33
            me_1 = 2.*p_centr[0]-1.*e_1
            me_2 = 2.*p_centr[1]-1.*e_2
            me_3 = 2.*p_centr[2]-1.*e_3 
            
            e_11 = numpy.hstack([e_11, me_11])
            e_12 = numpy.hstack([e_12, me_12])
            e_13 = numpy.hstack([e_13, me_13])
                                          
            e_21 = numpy.hstack([e_21, me_21])
            e_22 = numpy.hstack([e_22, me_22])
            e_23 = numpy.hstack([e_23, me_23])
                                          
            e_31 = numpy.hstack([e_31, me_31])
            e_32 = numpy.hstack([e_32, me_32])
            e_33 = numpy.hstack([e_33, me_33])

            e_1 = numpy.hstack([e_1, me_1])
            e_2 = numpy.hstack([e_2, me_2])
            e_3 = numpy.hstack([e_3, me_3])

            
        x_s = numpy.round(numpy.mod(e_11*x + e_12*y + e_13*z + e_1, 1), 5)
        y_s = numpy.round(numpy.mod(e_21*x + e_22*y + e_23*z + e_2, 1), 5)
        z_s = numpy.round(numpy.mod(e_31*x + e_32*y + e_33*z + e_3, 1), 5)
            
        xyz_s = numpy.vstack([x_s, y_s, z_s])
        
        xyz_s_un, unique_inverse = numpy.unique(xyz_s, return_inverse=True, axis=1)
        x_s, y_s, z_s = xyz_s_un[0, :], xyz_s_un[1, :], xyz_s_un[2, :]
        ind = (numpy.where((x-x_s)**2+(y-y_s)**2+(z-z_s)**2 < 0.00001))[0][0]

        flag = unique_inverse == ind
        o_11, o_12, o_13 = e_11[flag], e_12[flag], e_13[flag]
        o_21, o_22, o_23 = e_21[flag], e_22[flag], e_23[flag]
        o_31, o_32, o_33 = e_31[flag], e_32[flag], e_33[flag]
        o_3, o_2, o_3 = e_3[flag], e_2[flag], e_3[flag]
        return o_11, o_12, o_13, o_21, o_22, o_23, o_31, o_32, o_33, o_3, o_2, o_3


    
    def calc_atom_mult(self, np_x, np_y, np_z):
        """
        calculate atom multiplicity
        """
        lmult=[]
        for x, y, z in zip(np_x, np_y, np_z):
            np_x_s = self.calc_xyz_mult(x, y, z)[0]
            lmult.append(np_x_s.shape[0])
        np_multiplicity = numpy.array(lmult, dtype=int)
        
        return np_multiplicity
    
    def _trans_str_to_el_symm(self, str1):
        """
        transform string to element of symmetry: (x,y,-z) -> 0.0 1 0 0  0.0 0 1 0  0.0 0 0 -1
        """
        str2 = "".join(str1.split(" "))
        l_help1, l_help2, l_help3 = [], [], []
        l_help1 = [hh for hh in str2.split('(') if hh != ""]
        [l_help2.extend(hh.split(')')) for hh in l_help1 if hh != ""]
        [l_help3.extend(hh.split(',')) for hh in l_help2 if hh != ""]
        l_Ax = ['x', 'y', 'z']
        l_el_symm = []
        for hh in l_help3:
            el_symm_h = [0.0, 0, 0, 0]
            str_h = hh
            for i_num, Ax in enumerate(l_Ax):
                if (str_h.find(Ax) != -1):
                    if (str_h.find("+"+Ax) != -1):
                        el_symm_h[i_num+1] = 1
                        str_h = "".join(str_h.split("+"+Ax))
                    elif (str_h.find("-"+Ax) != -1):
                        el_symm_h[i_num+1] = -1
                        str_h = "".join(str_h.split("-"+Ax))
                    else:
                        el_symm_h[i_num+1] = 1
                        str_h = "".join(str_h.split(Ax))
            if (str_h==""):
                pass
            elif (str_h.find("/") != -1):
                l_help1 = str_h.split("/")
                el_symm_h[0] = float(l_help1[0])/float(l_help1[1])
            else:
                el_symm_h[0] = float(str_h)
            l_el_symm.append(el_symm_h)
        el_symm = []
        [el_symm.extend(hh) for hh in l_el_symm]
        return el_symm

    def _trans_el_symm_to_str(self):
        ls_out = []
        l_el_symm = self.el_symm
        centr = self.centr
        if centr:
            p_centr = self.p_centr
            ls_out.append("inversion center at ({:.3f}, {:.3f}, {:.3f})".format(
                    p_centr[0], p_centr[1], p_centr[2]))
        if l_el_symm == []:
            return ""
        for el_symm in l_el_symm:
            s_x = ""
            if el_symm[0] != 0.: s_x+="{:.3f}".format(el_symm[0])
            if el_symm[1] == 1: s_x+="+x"
            if el_symm[1] == -1: s_x+="-x"
            if el_symm[2] == 1: s_x+="+y"
            if el_symm[2] == -1: s_x+="-y"
            if el_symm[3] == 1: s_x+="+z"
            if el_symm[3] == -1: s_x+="-z"
            if s_x.startswith("+"): s_x = s_x[1:]
            
            s_y = ""
            if el_symm[4] != 0.: s_y+="{:.3f}".format(el_symm[4])
            if el_symm[5] == 1: s_y+="+x"
            if el_symm[5] == -1: s_y+="-x"
            if el_symm[6] == 1: s_y+="+y"
            if el_symm[6] == -1: s_y+="-y"
            if el_symm[7] == 1: s_y+="+z"
            if el_symm[7] == -1: s_y+="-z"
            if s_y.startswith("+"): s_y = s_y[1:]

            s_z = ""
            if el_symm[8] != 0.: s_z+="{:.3f}".format(el_symm[8])
            if el_symm[9]==1: s_z+="+x"
            if el_symm[9] == -1: s_z+="-x"
            if el_symm[10] == 1: s_z+="+y"
            if el_symm[10] == -1: s_z+="-y"
            if el_symm[11] == 1: s_z+="+z"
            if el_symm[11] == -1: s_z+="-z"
            if s_z.startswith("+"): s_z = s_z[1:]

            line=" {:}, {:}, {:}".format(s_x, s_y, s_z)
            ls_out.append(line)
        return "\n".join(ls_out)

    def calc_asymmetric_cell(self, n_a, n_b, n_c):
        """
        give the numbers in asymmetric cell

        n_a is the number of points along a axis
        n_b is the numper of points along b axis
        n_c is the numper of points along c axis

        na, n_b, nc should be divided on 24: 8 and 3

        output: 
        l_coord is a list of coordinates in asymmetric cell (frac_x = n_x/n_a and so on)
        l_symm contains a list of symmetry given as (n_symm, centr, n_orig)
        """

        n_a_new = int(round(n_a/24))*24
        n_b_new = int(round(n_b/24))*24
        n_c_new = int(round(n_c/24))*24

        #print("na: {:}, n_b: {:}, n_c: {:}".format(n_a_new, n_b_new, n_c_new))

        l_el_symm = self.el_symm
        f_centr = self.centr
        p_centr = self.p_centr
        l_orig = self.orig
        l_coord = []
        
        spgr_choice = self.spgr_choice 
        spgr_name = self.spgr_name 
        spgr_number = self.spgr_number
        
        if (spgr_number==227) & (spgr_choice=="2"):
            n_a_new = int(round(n_a/8))*8
            for n_x in range(-n_a_new//8, 3*n_a_new//8+1):
                for n_y in range(-n_a_new//8, 0+1):
                    for n_z in range(-n_a_new//4, 0+1):
                        cond_1 = (n_y < min([n_a_new//4-n_x, n_x]))
                        cond_2 = (n_z >= -n_y-n_a_new//4)
                        cond_3 = (n_z <= n_y)
                        if (cond_1 & cond_2 & cond_3):
                            coord_x, coord_y = float(n_x)/float(n_a_new), float(n_y)/float(n_a_new)
                            coord_z = float(n_z)/float(n_a_new)
                            #print(" {:3} {:3} {:3}".format(n_x, n_y, n_z), " ", " {:6.3f} {:6.3f} {:6.3f}".format(coord_x, coord_y, coord_z))
                            l_coord.append((coord_x, coord_y, coord_z))
        return l_coord

    def from_cif(self, string: str):
        cif_global = Global()
        flag = cif_global.take_from_string(string)
        if not flag:
            return False
        flag = False
        if cif_global.is_value("_space_group_name_H-M_alt"):
            self.spgr_name = cif_global["_space_group_name_H-M_alt"].value
        if cif_global.is_value("_space_group_it_coordinate_system_code"):
            self.spgr_choice = cif_global["_space_group_it_coordinate_system_code"].value
        return True

    @property
    def to_cif(self):
        ls_out = []
        ls_out.append("_space_group_name_H-M_alt {:}".format(self.spgr_name))
        ls_out.append("_space_group_it_coordinate_system_code {:}".format(self.spgr_choice))
        return "\n".join(ls_out)