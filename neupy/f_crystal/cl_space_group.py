"""
define classes to describe crystal 
"""
__author__ = 'ikibalin'
__version__ = "2019_04_06"
import os
import numpy

from neupy.f_common.cl_variable import Variable
from neupy.f_interface.cl_abstract_space_group import AbstractSpaceGroup

class SpaceGroup(AbstractSpaceGroup):
    """
    Space Group
    """
    def __init__(self, spgr_given_name = "P1", spgr_choice = "1",
                 f_dir_prog = os.path.dirname(__file__)):
        #         f_dir_prog = os.getcwd()):
        super(SpaceGroup, self).__init__()
        
        self._p_spgr_given_name = None
        self._p_spgr_choice = None
        self._p_f_dir_prog = None
        self._p_spgr_table = None

        self._p_singony = "Triclinic"

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

        self._refresh(spgr_given_name, spgr_choice, f_dir_prog)
        self.set_val()
        
    def __repr__(self):
        ls_out = ["SpaceGroup:\n given name: {:}\n choice: {:}".format(self._p_spgr_given_name, self._p_spgr_choice)]
        if self._p_spgr_name is not None:
            ls_out.append(" name: {:}".format(self._p_spgr_name))
        if self._p_spgr_number is not None:
            ls_out.append(" number: {:}".format(self._p_spgr_number))
        ls_out.append(" {:}\n directory: '{:}'".format(self._trans_el_symm_to_str(), self._p_f_dir_prog))
        return "\n".join(ls_out)

    def _refresh(self, spgr_given_name, spgr_choice, f_dir_prog):
        
        if f_dir_prog is not None:
            f_itables = os.path.join(f_dir_prog, "tables", "itables.txt")
            self._read_el_cards(f_itables)        
            self._p_f_dir_prog = f_dir_prog
        if spgr_given_name is not None:
            hh = "".join(spgr_given_name.split())
            self._p_spgr_given_name = hh
        if spgr_choice is not None:
            if isinstance(spgr_choice, float):
                spgr_choice = "{:}".format(int(spgr_choice))
            self._p_spgr_choice = spgr_choice = spgr_choice.strip()
            

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
spgr_given_name is number or name of the space group
spgr_choice is choise of origin, 1, 2, "abc", "bac"
f_dir_prog is directory where the file "itables.txt" it is 

centr is inversion center
p_centr is position of inversin center
el_symm is element of symmetry
orig is packing
spgr_name is name of space group
spgr_number is number of space group
singony is singony of the space group

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
                dcard = {"number":lhelp[0], "name": lhelp[1], "singony": lhelp[2]}
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
        flag = False
        for dcard in spgr_table:
            if (((dcard["number"] == spgr_n)|(dcard["name"] == spgr_name))&(dcard["choice"][0] == spgr_choice)):
                flag = True
                break
        if (not flag):
            print("Space group is not found: {:} {:} {:}".format(spgr_n, spgr_name, spgr_choice))
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
        singony = dcard["singony"]
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
        self._p_singony = singony
        
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

        x_s = numpy.round(numpy.mod(r_11*x + r_12*y + r_13*z + b_1, 1), 4)
        y_s = numpy.round(numpy.mod(r_21*x + r_22*y + r_23*z + b_2, 1), 4)
        z_s = numpy.round(numpy.mod(r_31*x + r_32*y + r_33*z + b_3, 1), 4)

        x_o = [orig[0] for orig in lorig]
        y_o = [orig[1] for orig in lorig]
        z_o = [orig[2] for orig in lorig]
        
        x_s_2d, x_o_2d = numpy.meshgrid(x_s, x_o)
        y_s_2d, y_o_2d = numpy.meshgrid(y_s, y_o)
        z_s_2d, z_o_2d = numpy.meshgrid(z_s, z_o)
        
        x_s_2d = numpy.round(numpy.mod(x_s_2d+x_o_2d, 1), 4)
        y_s_2d = numpy.round(numpy.mod(y_s_2d+y_o_2d, 1), 4)
        z_s_2d = numpy.round(numpy.mod(z_s_2d+z_o_2d, 1), 4)

        x_s = x_s_2d.flatten()
        y_s = y_s_2d.flatten()
        z_s = z_s_2d.flatten()

        if centr:
            x_s_h = numpy.round(numpy.mod(2.*p_centr[0]-1.*x_s, 1), 4)
            y_s_h = numpy.round(numpy.mod(2.*p_centr[1]-1.*y_s, 1), 4)
            z_s_h = numpy.round(numpy.mod(2.*p_centr[2]-1.*z_s, 1), 4)
            x_s =numpy.hstack([x_s, x_s_h])
            y_s =numpy.hstack([y_s, y_s_h])
            z_s =numpy.hstack([z_s, z_s_h])
                        
        xyz_s = numpy.vstack([x_s, y_s, z_s])
        
        xyz_s_un = numpy.unique(xyz_s, axis=1)
        n_atom = int(round(xyz_s.shape[1]*1./xyz_s_un.shape[1]))
        x_s, y_s, z_s = xyz_s_un[0, :], xyz_s_un[1, :], xyz_s_un[2, :]
        
        
        
        return x_s, y_s, z_s, n_atom
    
    def calc_el_symm_for_xyz(self, x, y, z):
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


        x_o = [orig[0] for orig in lorig]
        y_o = [orig[1] for orig in lorig]
        z_o = [orig[2] for orig in lorig]
        
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

    def _trans_el_symm_to_str(self):
        ls_out = []
        l_el_symm = self._p_el_symm
        centr = self._p_centr
        if centr:
            p_centr = self._p_p_centr
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

    def calc_assymmetric_cell(self, n_a, n_b, n_c):
        """
        give the numbers in assymmetric cell

        n_a is the number of points along a axis
        n_b is the numper of points along b axis
        n_c is the numper of points along c axis

        na, n_b, nc should be devided on 24: 8 and 3

        output: 
        l_coord is a list of coordinates in assymmetric cell (frac_x = n_x/n_a and so on)
        l_symm contains a list of symmetry given as (n_symm, centr, n_orig)
        """

        n_a_new = int(round(n_a/24))*24
        n_b_new = int(round(n_b/24))*24
        n_c_new = int(round(n_c/24))*24

        #print("na: {:}, n_b: {:}, n_c: {:}".format(n_a_new, n_b_new, n_c_new))

        l_el_symm = self._p_el_symm
        f_centr = self._p_centr
        p_centr = self._p_p_centr
        l_orig = self._p_orig

        l_coord = []
        
        spgr_choice = self._p_spgr_choice 
        spgr_name = self._p_spgr_name 
        spgr_number = self._p_spgr_number
        
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

