"""
define classes to describe AtomType
"""
__author__ = 'ikibalin'
__version__ = "2019_08_26"
import os
import numpy


from neupy.f_common.cl_fitable import Fitable

class AtomType(object):
    """
    Data items in the ATOM_TYPE category record details about
    properties of the atoms that occupy the atom sites, such as the
    atomic scattering factors.
    
    Description in cif file:

    loop_
    _atom_type_symbol
    _atom_type_oxidation_number
    _atom_type_number_in_cell
    _atom_type_scat_dispersion_real
    _atom_type_scat_dispersion_imag
    _atom_type_scat_source
      C 0 72  .017  .009  International_Tables_Vol_IV_Table_2.2B
      H 0 100  0     0    International_Tables_Vol_IV_Table_2.2B
      O 0 12  .047  .032  International_Tables_Vol_IV_Table_2.2B
      N 0 4   .029  .018  International_Tables_Vol_IV_Table_2.2B

    reference: https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Catom_type.html
    
    NOTE: Now is much more than this but it should be reduced in future
    """    
    def __init__(self, name="H", 
                 type_n="H", type_m="Fe3", #todo replace by one atom_type_symbol
                 flag_m=False, 
                 x=0., y=0., z=0., 
                 adp_type="Uiso", b_iso=0., 
                 u_11=0., u_22=0., u_33=0., 
                 u_12=0., u_13=0., u_23=0., 
                 chi_type='Cani', 
                 chi_11=0., chi_22=0., chi_33=0., chi_12=0., 
                 chi_13=0., chi_23=0., 
                 kappa=1., factor_lande=2.,
                 occupation = 1., 
                 f_dir_prog = os.path.dirname(__file__)):
        super(AtomType, self).__init__()

        self.__atom_tyme_name = None
        self.__atom_tyme_type_n = None
        self.__atom_tyme_type_m = None
        self.__atom_tyme_flag_m = None
        
        self._p_b_occupation = None

        self._p_x = None
        self._p_y = None
        self._p_z = None
        self._p_b_scat = None

        self._p_adp_type = None
        self._p_b_iso = None
        self._p_u_11 = None 
        self._p_u_22 = None 
        self._p_u_33 = None 
        self._p_u_12 = None 
        self._p_u_13 = None 
        self._p_u_23 = None 
        
        self._p_chi_type = None 
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
        self._p_f_dir_prog = None
        self._handbook_nucl = []
        self._handbook_mag = []
        
        self._refresh(name, type_n, type_m, flag_m, x, y, z, adp_type, b_iso, 
                      u_11, u_22, u_33, u_12, u_13, u_23, chi_type, chi_11, 
                      chi_22, chi_33, chi_12, chi_13, chi_23, kappa, 
                      factor_lande, occupation, f_dir_prog)

    def __repr__(self):
        ls_out = """AtomType:
 name: {}            
 type_n: {:} (b_scat: {:} cm**-12)
 occupation: {:}
 fract  x: {:}  y: {:}  z: {:}\n  adp_type: {:}\n""".format(
 self._p_name, self._p_type_n, self._p_b_scat, self._p_occupation, self._p_x, 
 self._p_y, self._p_z, self._p_adp_type)
        if self._p_adp_type == "uiso":
             ls_out += """ b_iso: {:} Ang**2""".format(self._p_b_iso)
        elif self._p_adp_type == "uani":
             ls_out += """ u_11: {:}, u_22: {:}, u_33: {:} (in Ang**2)
 u_12: {:}, u_13: {:}, u_23: {:} (in Ang**2)""".format(self._p_u_11, 
 self._p_u_22, self._p_u_33, self._p_u_12, self._p_u_13, self._p_u_23)
        if self._p_flag_m:
            ls_out += """\n\n type_m: {:}
 kappa: {:}
 lande factor: {:}
 chi_type: {}
 chi_11: {:}, chi_22: {:}, chi_33: {:} (in mu_B)
 chi_12: {:}, chi_13: {:}, chi_23: {:} (in mu_B) 
 j0   A: {:}, a: {:}, B: {:}, b: {:}
      C: {:}, c: {:}, D: {:}
 j2   A: {:}, a: {:}, B: {:}, b: {:}
      C: {:}, c: {:}, D: {:}
""".format(self._p_type_m, self._p_kappa, self._p_factor_lande, 
self._p_chi_type, self._p_chi_11, self._p_chi_22, 
self._p_chi_33, self._p_chi_12, self._p_chi_13, self._p_chi_23, self._p_j0_A,
self._p_j0_a, self._p_j0_B, self._p_j0_b, self._p_j0_C, self._p_j0_c, 
self._p_j0_D, self._p_j2_A, self._p_j2_a, self._p_j2_B, self._p_j2_b, 
self._p_j2_C, self._p_j2_c, self._p_j2_D)
        return ls_out
    
    def _refresh(self, name, type_n, type_m, flag_m, x, y, z, adp_type, b_iso, 
                 u_11, u_22, u_33, u_12, u_13, u_23, chi_type, chi_11, chi_22, 
                 chi_33, chi_12, chi_13, chi_23, kappa, factor_lande, 
                 occupation, f_dir_prog):
        
        if f_dir_prog is not None:
            self._p_f_dir_prog = f_dir_prog
            self._load_handbook_n()
            self._load_handbook_m()

        if name is not None:
            self._p_name = name
        if type_n is not None:
            self._p_type_n = type_n
            self._get_b_scat(type_n)
        if type_m is not None:
            self._p_type_m = type_m
            self._get_j0j2(type_m)
        if flag_m is not None:
            self._p_flag_m = flag_m

        if x is not None:
            #it is not good solution, but may be it is not needed at all
            if isinstance(x, Variable):
                x[0] = numpy.mod(x[0], 1.)
                self._p_x = x
            else:
                self._p_x = numpy.mod(x, 1.)
        if y is not None:
            if isinstance(y, Variable):
                y[0] = numpy.mod(y[0], 1.)
                self._p_y = y
            else:
                self._p_y = numpy.mod(y, 1.)
        if z is not None:
            if isinstance(z, Variable):
                z[0] = numpy.mod(z[0], 1.)
                self._p_z = z
            else:
                self._p_z = numpy.mod(z, 1.)
    
        if occupation is not None:
            self._p_occupation = occupation

        if adp_type is not None:
            if adp_type.lower().startswith("uani"):
                self._p_adp_type = "uani"
            else:
                self._p_adp_type = "uiso"
        if b_iso is not None:
            self._p_b_iso = b_iso
        if u_11 is not None:
            self._p_u_11 = u_11 
        if u_22 is not None:
            self._p_u_22 = u_22 
        if u_33 is not None:
            self._p_u_33 = u_33 
        if u_12 is not None:
            self._p_u_12 = u_12 
        if u_13 is not None:
            self._p_u_13 = u_13 
        if u_23 is not None:
            self._p_u_23 = u_23 
    
            
        if chi_type is not None:
            if chi_type.lower().startswith("cani"):
                self._p_chi_type = "cani"
            else:
                self._p_chi_type = "ciso"
        if chi_11 is not None:
            self._p_chi_11 = chi_11 
        if chi_22 is not None:
            self._p_chi_22 = chi_22 
        if chi_33 is not None:
            self._p_chi_33 = chi_33 
        if chi_12 is not None:
            self._p_chi_12 = chi_12 
        if chi_13 is not None:
            self._p_chi_13 = chi_13 
        if chi_23 is not None:
            self._p_chi_23 = chi_23 
            
        if kappa is not None:
            self._p_kappa = kappa 
        if factor_lande is not None:
            self._p_factor_lande = factor_lande 

            
    def set_val(self, name=None, type_n=None, type_m=None, flag_m=None, x=None, 
                y=None, z=None, adp_type=None, b_iso=None, u_11=None, 
                u_22=None, u_33=None, u_12=None, u_13=None, u_23=None, 
                chi_type=None, chi_11=None, chi_22=None, chi_33=None, 
                chi_12=None, chi_13=None, chi_23=None, kappa=None, 
                factor_lande=None, occupation=None, f_dir_prog=None):
        self._refresh(name, type_n, type_m, flag_m, x, y, z, adp_type, b_iso, 
                      u_11, u_22, u_33, u_12, u_13, u_23, chi_type, chi_11, 
                      chi_22, chi_33, chi_12, chi_13, chi_23, kappa, 
                      factor_lande, occupation, f_dir_prog)
        
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
Parameters of AtomType:
name is the name/label of atom
type_n is nuclear type of atom to defince b_scat
type_m is magnetic type of atom to define magetic atom
flag_m is True if atom magnetic and False if not.

x, y, z is fraction of atom in the crystal
b_scat is scattering amplitude
occupation is occupation factor
b_iso is isotropical atomic vibrations

adp_type
u_11, u_22, u_33 
u_12, u_13, u_23       is anisotropical atomic vibrations in Ang**2

chi_type
chi_11, chi_22, chi_33
chi_12, chi_13, chi_23          is susceptibility in Bohr magneton

kappa is expansion/contraction coefficient (by default 1.)
factor_lande is factor lande (by default 2.)

j0_A, j0_a, j0_B, j0_b, j0_C, j0_c, j0_D is coefficient to calculate <j0>

j2_A, j2_a, j2_B, j2_b, j2_C, j2_c, j2_D is coefficient to calculate <j2>

f_dir_prog is directory with file 'bscat.tab', 'formmag.tab'
        """
        print(lsout)
        
    def _load_handbook_n(self):
        f_name = os.path.join(self._p_f_dir_prog, "tables", "bscat.tab")
        fid = open(f_name, 'r')
        lcont = fid.readlines()
        fid.close()
        lcont = [line for line in lcont if not(line.startswith("#"))]
        ldcard = []
        for line in lcont:
            lhelp = line.strip().split()
            
            sline = lhelp[2].replace("i","j")
            sline = sline.split("(")[0]
            try:
                if sline.rfind("j") != -1:
                    b_scat = 0.1*complex(sline)
                else:
                    b_scat = 0.1*float(sline)
            except:
                b_scat = 0.
            dcard = {"type_n": lhelp[0], "b_scat": b_scat}
            ldcard.append(dcard)
        self._handbook_nucl = ldcard
    
    def _load_handbook_m(self):
        f_name = os.path.join(self._p_f_dir_prog, "tables", "formmag.tab")
        fid = open(f_name, 'r')
        lcont = fid.readlines()
        fid.close()
        lcont = [line for line in lcont if line.startswith("F")]
        ldcard = []
        for line in lcont:
            lhelp = line.strip().split()
            dcard = {"type_m": lhelp[1], "order": int(lhelp[2]),
                     "A": float(lhelp[3]),"a": float(lhelp[4]),
                     "B": float(lhelp[5]),"b": float(lhelp[6]),
                     "C": float(lhelp[7]),"c": float(lhelp[8]),
                     "D": float(lhelp[9])}
            ldcard.append(dcard)
        self._handbook_mag = ldcard

    def _get_b_scat(self, type_n):
        """
        Take b_scat
        """

        ldcard = self._handbook_nucl 
        flag = False
        for dcard in ldcard:
            if (dcard["type_n"] == type_n):
                self._p_b_scat = dcard["b_scat"]
                flag = True
            elif flag:
                break
        if not(flag):
            print("Can not find b_scat for '{:}'".format(type_n))
            
    def _get_j0j2(self, type_m):
        """
        Take coefficients for <j0> and <j2>
        """
        ldcard = self._handbook_mag 
        flag_0, flag_2 = False, False
        for dcard in ldcard:
            if ((dcard["type_m"] == type_m)&(dcard["order"] == 0)):
                self._p_j0_A = dcard["A"]
                self._p_j0_a = dcard["a"]
                self._p_j0_B = dcard["B"]
                self._p_j0_b = dcard["b"]
                self._p_j0_C = dcard["C"]
                self._p_j0_c = dcard["c"]
                self._p_j0_D = dcard["D"]
                flag_0 = True
            elif ((dcard["type_m"] == type_m)&(dcard["order"] == 2)):
                self._p_j2_A = dcard["A"]
                self._p_j2_a = dcard["a"]
                self._p_j2_B = dcard["B"]
                self._p_j2_b = dcard["b"]
                self._p_j2_C = dcard["C"]
                self._p_j2_c = dcard["c"]
                self._p_j2_D = dcard["D"]
                flag_2 = True
            elif (flag_0 & flag_2):
                break
        if not(flag_0):
            print("Can not find coefficients <j0> for '{:}'".format(type_m))
        if not(flag_2):
            print("Can not find coefficients <j2> for '{:}'".format(type_m))

    def is_variable_phase(self):
        res = any([isinstance(self._p_x, Variable), 
                   isinstance(self._p_y, Variable),
                   isinstance(self._p_z, Variable)])
        return res

    def is_variable_adp(self):
        if self._p_adp_type == "uiso":
            res = isinstance(self._p_b_iso, Variable)
        else:
            res = any([isinstance(self._p_u_11, Variable), 
                   isinstance(self._p_u_22, Variable),
                   isinstance(self._p_u_33, Variable),
                   isinstance(self._p_u_12, Variable),
                   isinstance(self._p_u_13, Variable),
                   isinstance(self._p_u_23, Variable)])
        return res

    def is_variable_magnetism(self):
        res_1 = any([isinstance(self._p_kappa, Variable), 
                   isinstance(self._p_factor_lande, Variable), 
                   isinstance(self._p_chi_11, Variable)])
        res_2 = False
        if self._p_chi_type == "cani":
            res_2 = any([isinstance(self._p_chi_22, Variable),
                   isinstance(self._p_chi_33, Variable),
                   isinstance(self._p_chi_12, Variable),
                   isinstance(self._p_chi_13, Variable),
                   isinstance(self._p_chi_23, Variable)])
        res = res_1 | res_2
        return res
    
    def is_variable(self):
        res = any([self.is_variable_phase(), 
                   self.is_variable_adp(),
                   self.is_variable_magnetism(),
                   isinstance(self._p_occupation, Variable)])
        return res

    def get_variables(self):
        l_variable = []
        if isinstance(self._p_x, Variable):
            l_variable.append(self._p_x)
        if isinstance(self._p_y, Variable):
            l_variable.append(self._p_y)
        if isinstance(self._p_z, Variable):
            l_variable.append(self._p_z)

        if isinstance(self._p_occupation, Variable):
            l_variable.append(self._p_occupation)
        if isinstance(self._p_b_scat, Variable):
            l_variable.append(self._p_b_scat)
            
        if self._p_adp_type == "uiso":
            if isinstance(self._p_b_iso, Variable):
                l_variable.append(self._p_b_iso)
        else:
            if isinstance(self._p_u_11, Variable):
                l_variable.append(self._p_u_11)
            if isinstance(self._p_u_22, Variable):
                l_variable.append(self._p_u_22)
            if isinstance(self._p_u_33, Variable):
                l_variable.append(self._p_u_33)
            if isinstance(self._p_u_12, Variable):
                l_variable.append(self._p_u_12)
            if isinstance(self._p_u_13, Variable):
                l_variable.append(self._p_u_13)
            if isinstance(self._p_u_23, Variable):
                l_variable.append(self._p_u_23)
            
        if isinstance(self._p_kappa, Variable):
            l_variable.append(self._p_kappa)
        if isinstance(self._p_factor_lande, Variable):
            l_variable.append(self._p_factor_lande)

        if isinstance(self._p_chi_11, Variable):
            l_variable.append(self._p_chi_11)
        if self._p_chi_type == "cani":  
            if isinstance(self._p_chi_22, Variable):
                l_variable.append(self._p_chi_22)
            if isinstance(self._p_chi_33, Variable):
                l_variable.append(self._p_chi_33)
            if isinstance(self._p_chi_12, Variable):
                l_variable.append(self._p_chi_12)
            if isinstance(self._p_chi_13, Variable):
                l_variable.append(self._p_chi_13)
            if isinstance(self._p_chi_23, Variable):
                l_variable.append(self._p_chi_23)

        return l_variable

    def calc_beta(self, cell):
        """
        calculate beta_ij from U_ij
        """
        ia, ib, ic = 1.*cell.get_val("ia"), 1.*cell.get_val("ib"), 1.*cell.get_val("ic")
        u_11, u_22, u_33 = 1.*self._p_u_11, 1.*self._p_u_22, 1.*self._p_u_33
        u_12, u_13, u_23 = 1.*self._p_u_12, 1.*self._p_u_13, 1.*self._p_u_23
        beta_11 = 2.*numpy.pi**2*u_11*ia**2
        beta_22 = 2.*numpy.pi**2*u_22*ib**2
        beta_33 = 2.*numpy.pi**2*u_33*ic**2
        beta_12 = 2.*numpy.pi**2*u_12*ia*ib
        beta_13 = 2.*numpy.pi**2*u_13*ia*ic
        beta_23 = 2.*numpy.pi**2*u_23*ib*ic
        return beta_11, beta_22, beta_33, beta_12, beta_13, beta_23


    def calc_chi_iso(self, cell):
        c_a = cell.get_val("cos_a")
        s_ib = cell.get_val("sin_ib")
        s_ig = cell.get_val("sin_ig")
        c_ib = cell.get_val("cos_ib")
        c_ig = cell.get_val("cos_ig")
        #not sure, it is better to check
        chi_11 = 1.*self._p_chi_11
        chi_22 = chi_11
        chi_33 = chi_11
        chi_12 = chi_11*c_ig
        chi_13 = chi_11*c_ib
        chi_23 = chi_11*(c_ib*c_ig-s_ib*s_ig*c_a)
        return chi_11, chi_22, chi_33, chi_12, chi_13, chi_23
    
    def apply_chi_iso(self, cell):
        chi_11, chi_22, chi_33, chi_12, chi_13, chi_23 = self.calc_chi_iso(cell)
        self._p_chi_22 = chi_22
        self._p_chi_33 = chi_33
        self._p_chi_12 = chi_12
        self._p_chi_13 = chi_13
        self._p_chi_23 = chi_23

        
    def apply_constraint(self, space_group, cell):
        chi_type = self.get_val("chi_type")
        if chi_type == "ciso":
            self.apply_chi_iso(cell)
        self.apply_space_group_constraint(space_group)

    def calc_constr_number(self, space_group):
        """
        according to table 1 in Peterse, Palm, Acta Cryst.(1966), 20, 147
        """
        x, y, z = 1.*self._p_x, 1.*self._p_y, 1.*self._p_z
        o_11, o_12, o_13, o_21, o_22, o_23, o_31, o_32, o_33, o_3, o_2, o_3 = space_group.calc_el_symm_for_xyz(x,y,z)
        
        b_11, b_22, b_33, b_12, b_13, b_23 = 107, 181, 41, 7, 19, 1
        
        i_11 = (o_11*o_11*b_11 + o_12*o_12*b_22 + o_13*o_13*b_33 + 
                o_11*o_12*b_12 + o_11*o_13*b_13 + o_12*o_13*b_23 + 
                o_12*o_11*b_12 + o_13*o_11*b_13 + o_13*o_12*b_23)
                
        i_22 = (o_21*o_21*b_11 + o_22*o_22*b_22 + o_23*o_23*b_33 + 
                o_21*o_22*b_12 + o_21*o_23*b_13 + o_22*o_23*b_23 + 
                o_22*o_21*b_12 + o_23*o_21*b_13 + o_23*o_22*b_23)
                
        i_33 = (o_31*o_31*b_11 + o_32*o_32*b_22 + o_33*o_33*b_33 + 
                o_31*o_32*b_12 + o_31*o_33*b_13 + o_32*o_33*b_23 + 
                o_32*o_31*b_12 + o_33*o_31*b_13 + o_33*o_32*b_23) 
                
        i_12 = (o_11*o_21*b_11 + o_12*o_22*b_22 + o_13*o_23*b_33 + 
                o_11*o_22*b_12 + o_11*o_23*b_13 + o_12*o_23*b_23 + 
                o_12*o_21*b_12 + o_13*o_21*b_13 + o_13*o_22*b_23) 
                
        i_13 = (o_11*o_31*b_11 + o_12*o_32*b_22 + o_13*o_33*b_33 + 
                o_11*o_32*b_12 + o_11*o_33*b_13 + o_12*o_33*b_23 + 
                o_12*o_31*b_12 + o_13*o_31*b_13 + o_13*o_32*b_23)
                
        i_23 = (o_21*o_31*b_11 + o_22*o_32*b_22 + o_23*o_33*b_33 + 
                o_21*o_32*b_12 + o_21*o_33*b_13 + o_22*o_33*b_23 + 
                o_22*o_31*b_12 + o_23*o_31*b_13 + o_23*o_32*b_23) 
        
        r_11, r_22, r_33, r_12, r_13, r_23 = i_11.sum(), i_22.sum(), i_33.sum(), i_12.sum(), i_13.sum(), i_23.sum()
        
        if r_13 == 0:
            if r_23 == 0:
                if r_12 == 0:
                    if r_11 == r_22:
                        if r_22 == r_33:
                            numb = 17
                        else:
                            numb = 8
                    elif r_22 == r_33:
                        numb = 12
                    else:
                        numb = 4
                elif r_11 == r_22:
                    if r_22 == 2*r_12:
                        numb = 16
                    else:
                        numb = 5
                elif r_22 == 2*r_12:
                    numb = 14
                else:
                    numb = 2
            elif r_22 == r_33:
                numb = 9
            else:
                numb = 3
        elif r_23 == 0:
            if r_22 == 2*r_12:
                numb = 13
            else:
                numb = 1
        elif r_23 == r_13:
            if r_22 == r_33:
                numb = 18
            else:
                numb = 6
        elif r_12 == r_13:
            numb = 10
        elif r_11 == 22:
            numb = 7
        elif r_22 == r_33:
            numb = 11
        elif r_22 == 2*r_12:
            numb = 15
        else:
            numb = 0 #no constraint
        return numb
    
    def apply_space_group_constraint(self, space_group):
        """
        according to table 1 in Peterse, Palm, Acta Cryst.(1966), 20, 147
        """
        numb = self.calc_constr_number(space_group)
        if numb == 1:
            self._p_u_12 = 0.
            self._p_u_23 = 0.

            self._p_chi_12 = 0.
            self._p_chi_23 = 0.
        elif numb == 2:
            self._p_u_23 = 0.
            self._p_u_13 = 0.

            self._p_chi_23 = 0.
            self._p_chi_13 = 0.
        elif numb == 3:
            self._p_u_12 = 0.
            self._p_u_13 = 0.

            self._p_chi_12 = 0.
            self._p_chi_13 = 0.
        elif numb == 4:
            self._p_u_12 = 0.
            self._p_u_23 = 0.
            self._p_u_13 = 0.

            self._p_chi_12 = 0.
            self._p_chi_23 = 0.
            self._p_chi_13 = 0.
        elif numb == 5:
            self._p_u_22 = 1.*self._p_u_11
            self._p_u_23 = 0.
            self._p_u_13 = 0.

            self._p_chi_22 = 1.*self._p_chi_11
            self._p_chi_23 = 0.
            self._p_chi_13 = 0.
        elif numb == 6:
            self._p_u_22 = 1.*self._p_u_11
            self._p_u_23 = 1.*self._p_u_13

            self._p_chi_22 = 1.*self._p_chi_11
            self._p_chi_23 = 1.*self._p_chi_13
        elif numb == 7:
            self._p_u_22 = 1.*self._p_u_11
            self._p_u_23 = -1.*self._p_u_13

            self._p_chi_22 = 1.*self._p_chi_11
            self._p_chi_23 = -1.*self._p_chi_13
        elif numb == 8:
            self._p_u_22 = 1.*self._p_u_11
            self._p_u_12 = 0.
            self._p_u_23 = 0.
            self._p_u_13 = 0.

            self._p_chi_22 = 1.*self._p_chi_11
            self._p_chi_12 = 0.
            self._p_chi_23 = 0.
            self._p_chi_13 = 0.
        elif numb == 9:
            self._p_u_33 = 1.*self._p_u_22
            self._p_u_12 = 0.
            self._p_u_13 = 0.

            self._p_chi_33 = 1.*self._p_chi_22
            self._p_chi_12 = 0.
            self._p_chi_13 = 0.
        elif numb == 10:
            self._p_u_33 = 1.*self._p_u_22
            self._p_u_13 = 1.*self._p_u_12 

            self._p_chi_33 = 1.*self._p_chi_22
            self._p_chi_13 = 1.*self._p_chi_12 
        elif numb == 11:
            self._p_u_33 = 1.*self._p_u_22
            self._p_u_13 = -1.*self._p_u_12 

            self._p_chi_33 = 1.*self._p_chi_22
            self._p_chi_13 = -1.*self._p_chi_12 
        elif numb == 12:
            self._p_u_33 = 1.*self._p_u_22
            self._p_u_12 = 0.
            self._p_u_23 = 0.
            self._p_u_13 = 0.

            self._p_chi_33 = 1.*self._p_chi_22
            self._p_chi_12 = 0.
            self._p_chi_23 = 0.
            self._p_chi_13 = 0.
        elif numb == 13:
            self._p_u_12 = 0.5*self._p_u_22
            self._p_u_23 = 0.

            self._p_chi_12 = 0.5*self._p_chi_22
            self._p_chi_23 = 0.
        elif numb == 14:
            self._p_u_12 = 0.5*self._p_u_22
            self._p_u_23 = 0.
            self._p_u_13 = 0.

            self._p_chi_12 = 0.5*self._p_chi_22
            self._p_chi_23 = 0.
            self._p_chi_13 = 0.
        elif numb == 15:
            self._p_u_12 = 0.5*self._p_u_22
            self._p_u_23 = 2.0*self._p_u_13

            self._p_chi_12 = 0.5*self._p_chi_22
            self._p_chi_23 = 2.0*self._p_chi_13
        elif numb == 16:
            self._p_u_22 = 1.0*self._p_u_11
            self._p_u_12 = 0.5*self._p_u_11
            self._p_u_23 = 0.
            self._p_u_13 = 0.

            self._p_chi_22 = 1.0*self._p_chi_11
            self._p_chi_12 = 0.5*self._p_chi_11
            self._p_chi_23 = 0.
            self._p_chi_13 = 0.
        elif numb == 17:
            self._p_u_22 = 1.0*self._p_u_11
            self._p_u_33 = 1.0*self._p_u_11
            self._p_u_12 = 0.
            self._p_u_23 = 0.
            self._p_u_13 = 0.

            self._p_chi_22 = 1.0*self._p_chi_11
            self._p_chi_33 = 1.0*self._p_chi_11
            self._p_chi_12 = 0.
            self._p_chi_23 = 0.
            self._p_chi_13 = 0.
        elif numb == 18:
            self._p_u_22 = 1.0*self._p_u_11
            self._p_u_33 = 1.0*self._p_u_11
            self._p_u_23 = 1.0*self._p_u_12
            self._p_u_13 = 1.0*self._p_u_12

            self._p_chi_22 = 1.0*self._p_chi_11
            self._p_chi_33 = 1.0*self._p_chi_11
            self._p_chi_23 = 1.0*self._p_chi_12
            self._p_chi_13 = 1.0*self._p_chi_12
            
  