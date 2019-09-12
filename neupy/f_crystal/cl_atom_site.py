"""
define classes to describe AtomSite
"""
__author__ = 'ikibalin'
__version__ = "2019_08_30"
import os
import numpy


from pystar import Global
from neupy.f_common.cl_fitable import Fitable
from neupy.f_crystal.cl_fract import Fract

class AtomSite(object):
    """
    Data items in the ATOM_SITE category record details about
    the atom sites in a crystal structure, such as the positional
    coordinates.

    
    Description in cif file:

    loop_                                     
    _atom_site_label          
    _atom_site_type_symbol   
    _atom_site_fract_x       
    _atom_site_fract_y       
    _atom_site_fract_z       
    _atom_site_adp_type       
    _atom_site_B_iso_or_equiv
    _atom_site_occupancy     
     Fe3A   Fe  0.12500 0.12500 0.12500  uani   0.0   1.0
     Fe3B   Fe  0.50000 0.50000 0.50000  uani   0.0   1.0
     O1     O   0.25521 0.25521 0.25521  uiso   0.0   1.0

    
    reference: https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Catom_site.html
    """    
    def __init__(self, label=[],  type_symbol=[], x=[], y=[], z=[], 
                       adp_type=[], b_iso=[], occupancy=[],
                       f_dir_prog = os.path.dirname(__file__)):
        super(AtomSite, self).__init__()

        self.__atom_site_label = []
        self.__atom_site_type_symbol = []
        self.__atom_site_fract_x = []
        self.__atom_site_fract_y = []
        self.__atom_site_fract_z = []
        self.__atom_site_adp_type = []
        self.__atom_site_b_iso_or_equiv = []
        self.__atom_site_occupancy = []
        self.__f_dir_prog = None
        #internal data
        self.__handbook_nucl = None
        self.__handbook_mag = None
        self.__scat_length_neutron = None
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



        self.label = label
        self.type_symbol = type_symbol
        self.x = x
        self.y = y
        self.z = z
        self.adp_type = adp_type
        self.b_iso = b_iso
        self.occupancy = occupancy
        self.f_dir_prog = f_dir_prog
        
    def __repr__(self):
        ls_out = ["AtomSite:"]
        ls_out.append(" label type_symbol x        y        z        adp_type b_iso    occupancy")
        for hh_1, hh_2, hh_3, hh_4, hh_5, hh_6, hh_7, hh_8 in zip(
            self.label, self.type_symbol, self.x, self.y, self.z, self.adp_type, self.b_iso, self.occupancy):
            ls_out.append(" {:8} {:8} {:8} {:8} {:8} {:8} {:8} {:8}".format(hh_1, hh_2,
                    hh_3.print_with_sigma, hh_4.print_with_sigma, hh_5.print_with_sigma,
                    hh_6, hh_7.print_with_sigma, hh_8.print_with_sigma))
        l_type_symbol = self.type_symbol
        l_ind = [l_type_symbol.index(hh) for hh in set(l_type_symbol)]
        l_scat = self.scat_length_neutron
        l_j0_A, l_j0_a, l_j0_B, l_j0_b, l_j0_C, l_j0_c, l_j0_D, l_j2_A, l_j2_a, l_j2_B, l_j2_b, l_j2_C, l_j2_c, l_j2_D = self.j0j2
        l_1, l_2, l_3 = [l_type_symbol[hh] for hh in l_ind], [l_scat[hh] for hh in l_ind], [l_j0_A[hh] for hh in l_ind]
        l_4, l_5, l_6 = [l_j0_a[hh] for hh in l_ind], [l_j0_B[hh] for hh in l_ind], [l_j0_b[hh] for hh in l_ind]
        l_7, l_8, l_9 = [l_j0_C[hh] for hh in l_ind], [l_j0_c[hh] for hh in l_ind], [l_j0_D[hh] for hh in l_ind]
        l_10, l_11, l_12 = [l_j2_A[hh] for hh in l_ind], [l_j2_a[hh] for hh in l_ind], [l_j2_B[hh] for hh in l_ind]
        l_13, l_14, l_15 = [l_j2_b[hh] for hh in l_ind], [l_j2_C[hh] for hh in l_ind], [l_j2_c[hh] for hh in l_ind]
        l_16 = [l_j2_D[hh] for hh in l_ind]

        for hh_1, hh_2, hh_3, hh_4, hh_5, hh_6, hh_7, hh_8, hh_9, hh_10, hh_11, hh_12, hh_13, hh_14, hh_15, hh_16 in zip(
            l_1, l_2, l_3, l_4, l_5, l_6, l_7, l_8, l_9, l_10, l_11, l_12, l_13, l_14, l_15, l_16):
            ls_out.append("\n for type_symbol '{:}'".format(hh_1))
            ls_out.append(" scat_length_neutron: {:.4f} cm**-12".format(hh_2))
            ls_out.append(" j0_A:{:8.4f},  j0_a:{:8.4f},  j0_B:{:8.4f},  j0_b:{:8.4f}".format(hh_3, hh_4, hh_5, hh_6))
            ls_out.append(" j0_C:{:8.4f},  j0_c:{:8.4f},  j0_D:{:8.4f}".format(hh_7, hh_8, hh_9))
            ls_out.append(" j2_A:{:8.4f},  j2_a:{:8.4f},  j2_B:{:8.4f},  j2_b:{:8.4f}".format(hh_10, hh_11, hh_12, hh_13))
            ls_out.append(" j2_C:{:8.4f},  j2_c:{:8.4f},  j2_D:{:8.4f}".format(hh_14, hh_15, hh_16))
        
        return "\n".join(ls_out)

    @property
    def label(self):
        """
        The _atom_site_label is a unique identifier for a particular site
        in the crystal. 

        Type: char

        reference: https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Iatom_site_label.html
        """
        return tuple(self.__atom_site_label)
    @label.setter
    def label(self, l_x):
        l_x_in = []
        for x in l_x:
            x_in = str(x).strip()
            l_x_in.append(x_in)
        self.__atom_site_label = l_x_in
        len_x = len(l_x_in)

        len_1 = len(self.__atom_site_type_symbol)
        if len_1 > len_x:
            self.__atom_site_type_symbol = self.__atom_site_type_symbol[:len_x]
        elif len_1 < len_x:
            l_fitable = ["Fe3+" for hh in range(len_x-len_1)] #default
            self.__atom_site_type_symbol.extend(l_fitable)

        len_1 = len(self.__atom_site_fract_x)
        if len_1 > len_x:
            self.__atom_site_fract_x = self.__atom_site_fract_x[:len_x]
        elif len_1 < len_x:
            l_fitable = [Fitable(value=0., name="_atom_site_fract_x") for hh in range(len_x-len_1)]
            self.__atom_site_fract_x.extend(l_fitable)

        len_1 = len(self.__atom_site_fract_y)
        if len_1 > len_x:
            self.__atom_site_fract_y = self.__atom_site_fract_y[:len_x]
        elif len_1 < len_x:
            l_fitable = [Fitable(value=0., name="_atom_site_fract_y") for hh in range(len_x-len_1)]
            self.__atom_site_fract_y.extend(l_fitable)

        len_1 = len(self.__atom_site_adp_type)
        if len_1 > len_x:
            self.__atom_site_adp_type = self.__atom_site_adp_type[:len_x]
        elif len_1 < len_x:
            l_fitable = ["uiso" for hh in range(len_x-len_1)] #deafault value
            self.__atom_site_adp_type.extend(l_fitable)

        len_1 = len(self.__atom_site_b_iso_or_equiv)
        if len_1 > len_x:
            self.__atom_site_b_iso_or_equiv = self.__atom_site_b_iso_or_equiv[:len_x]
        elif len_1 < len_x:
            l_fitable = [Fitable(value=0., name="_atom_site_b_iso_or_equiv") for hh in range(len_x-len_1)]
            self.__atom_site_b_iso_or_equiv.extend(l_fitable)

        len_1 = len(self.__atom_site_occupancy)
        if len_1 > len_x:
            self.__atom_site_occupancy = self.__atom_site_occupancy[:len_x]
        elif len_1 < len_x:
            l_fitable = [Fitable(value=1., name="_atom_site_occupancy") for hh in range(len_x-len_1)]
            self.__atom_site_occupancy.extend(l_fitable)

    @property
    def type_symbol(self):
        """
        A code to identify the atom species (singular or plural)
        occupying this site.
        This code must match a corresponding _atom_type_symbol. 

        Type: char

        reference: https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Iatom_site_type_symbol.html
        """
        return tuple(self.__atom_site_type_symbol)
    @type_symbol.setter
    def type_symbol(self, l_x):
        l_x_in = []
        for x in l_x:
            x_in = str(x).strip()
            l_x_in.append(x_in)
        len_x = len(l_x_in)
        len_1 = len(self.__atom_site_type_symbol)
        if len_1 < len_x:
            l_x_in = l_x_in[:len_1]
        elif len_1 > len_x:
            l_x_in.extend(["Fe3+" for hh in range(len_1-len_x)])
        self.__atom_site_type_symbol = l_x_in
        self.__scat_length_neutron = None
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

    @property
    def x(self):
        """
        Atom-site coordinates as fractions of the _cell_length_ values.

        Default: 0.

        Type: float

        reference: https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Iatom_site_fract_.html
        """
        return tuple(self.__atom_site_fract_x)
    @x.setter
    def x(self, l_x):
        l_fitable = []
        for x in l_x:
            if isinstance(x, Fitable):
                x_in = x
            else:
                x_in = Fitable()
                flag = x_in.take_it(x)
            l_fitable.append(x_in)
        len_x = len(l_fitable)
        len_1 = len(self.__atom_site_label)
        if len_1 < len_x:
            l_fitable = l_fitable[:len_1]
        elif len_1 > len_x:
            l_fitable.extend([Fitable(value=0., name="_atom_site_fract_x") for hh in range(len_1-len_x)])
        self.__atom_site_fract_x = l_fitable


    @property
    def y(self):
        """
        see help for x parameter
        """
        return tuple(self.__atom_site_fract_y)
    @y.setter
    def y(self, l_x):
        l_fitable = []
        for x in l_x:
            if isinstance(x, Fitable):
                x_in = x
            else:
                x_in = Fitable()
                flag = x_in.take_it(x)
            l_fitable.append(x_in)
        len_x = len(l_fitable)
        len_1 = len(self.__atom_site_label)
        if len_1 < len_x:
            l_fitable = l_fitable[:len_1]
        elif len_1 > len_x:
            l_fitable.extend([Fitable(value=0., name="_atom_site_fract_y") for hh in range(len_1-len_x)])
        self.__atom_site_fract_y = l_fitable


    @property
    def z(self):
        """
        see help for x parameter
        """
        return tuple(self.__atom_site_fract_z)
    @z.setter
    def z(self, l_x):
        l_fitable = []
        for x in l_x:
            if isinstance(x, Fitable):
                x_in = x
            else:
                x_in = Fitable()
                flag = x_in.take_it(x)
            l_fitable.append(x_in)
        len_x = len(l_fitable)
        len_1 = len(self.__atom_site_label)
        if len_1 < len_x:
            l_fitable = l_fitable[:len_1]
        elif len_1 > len_x:
            l_fitable.extend([Fitable(value=0., name="_atom_site_fract_z") for hh in range(len_1-len_x)])
        self.__atom_site_fract_z = l_fitable

    @property
    def adp_type(self):
        """
        A standard code used to describe the type of atomic displacement
        parameters used for the site. 

        Type: [uani, uiso, uovl, umpe, bani, biso, bovl] #supported only [uani, uiso]

        reference: https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Iatom_site_adp_type.html
        """
        return tuple(self.__atom_site_adp_type)
    @adp_type.setter
    def adp_type(self, l_x):
        l_list = ["uani", "uiso"]
        l_x_in = []
        for x in l_x:
            x_in = str(x).strip().lower()
            if (x_in not in l_list): x_in = "uiso"
            l_x_in.append(x_in)
        len_x = len(l_x_in)
        len_1 = len(self.__atom_site_adp_type)
        if len_1 < len_x:
            l_x_in = l_x_in[:len_1]
        elif len_1 > len_x:
            l_x_in.extend(["uiso" for hh in range(len_1-len_x)])
        self.__atom_site_adp_type = l_x_in


    @property
    def b_iso(self):
        """
        Isotropic atomic displacement parameter, or equivalent isotropic
        atomic displacement parameter, B(equiv), in angstroms squared,
        calculated from anisotropic displacement components.

        B(equiv) = (1/3) sum~i~[sum~j~(B^ij^ a*~i~ a*~j~ a~i~ a~j~)]

        a     = the real-space cell lengths
        a*    = the reciprocal-space cell lengths
        B^ij^ = 8 pi^2^ U^ij^

        Ref: Fischer, R. X. & Tillmanns, E. (1988). Acta Cryst. C44,
             775-776.

        The IUCr Commission on Nomenclature recommends against the use
        of B for reporting atomic displacement parameters. U, being
        directly proportional to B, is preferred.

        Default: 0.

        reference: https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Iatom_site_B_iso_or_equiv.html
        """
        return tuple(self.__atom_site_b_iso_or_equiv)
    @b_iso.setter
    def b_iso(self, l_x):
        l_fitable = []
        for x in l_x:
            if isinstance(x, Fitable):
                x_in = x
            else:
                x_in = Fitable()
                flag = x_in.take_it(x)
            l_fitable.append(x_in)
        len_x = len(l_fitable)
        len_1 = len(self.__atom_site_label)
        if len_1 < len_x:
            l_fitable = l_fitable[:len_1]
        elif len_1 > len_x:
            l_fitable.extend([Fitable(value=0., name="_atom_site_b_iso_or_equiv") for hh in range(len_1-len_x)])
        self.__atom_site_b_iso_or_equiv = l_fitable


    @property
    def occupancy(self):
        """
        The fraction of the atom type present at this site.
        The sum of the occupancies of all the atom types at this site
        may not significantly exceed 1.0 unless it is a dummy site. The
        value must lie in the 99.97% Gaussian confidence interval
        -3u =< x =< 1 + 3u. The _enumeration_range of 0.0:1.0 is thus
        correctly interpreted as meaning (0.0 - 3u) =< x =< (1.0 + 3u).

        Default: 1.0
        
        reference: https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Iatom_site_occupancy.html
        """
        return tuple(self.__atom_site_occupancy)
    @occupancy.setter
    def occupancy(self, l_x):
        l_fitable = []
        for x in l_x:
            if isinstance(x, Fitable):
                x_in = x
            else:
                x_in = Fitable()
                flag = x_in.take_it(x)
            if x_in.value < 0.: x_in.value = 0.
            if x_in.value > 0.: x_in.value = 1.
            l_fitable.append(x_in)
        len_x = len(l_fitable)
        len_1 = len(self.__atom_site_label)
        if len_1 < len_x:
            l_fitable = l_fitable[:len_1]
        elif len_1 > len_x:
            l_fitable.extend([Fitable(value=1., name="_atom_site_occupancy") for hh in range(len_1-len_x)])
        self.__atom_site_occupancy = l_fitable

    @property
    def f_dir_prog(self):
        """

        reference:
        """
        return self.__f_dir_prog
    @f_dir_prog.setter
    def f_dir_prog(self, x):
        self.__f_dir_prog = x
        self._load_handbook_n()    
        self._load_handbook_m()

    def _show_message(self, s_out: str):
        print("***  Error ***")
        print(s_out)

    @property
    def is_defined(self):
        """
        Output: True if all started parameters are given
        """
        cond = (self.label != [])
        return cond
    @property
    def is_variable(self):
        res = (any([hh.refinement for hh in self.x]) | 
               any([hh.refinement for hh in self.y]) |
               any([hh.refinement for hh in self.z]) |
               any([hh.refinement for hh in self.b_iso]) |
               any([hh.refinement for hh in self.occupancy]))
        return res


    def get_variables(self):
        l_variable = []
        l_variable.extend([hh for hh in self.x if hh.refinement])
        l_variable.extend([hh for hh in self.y if hh.refinement])
        l_variable.extend([hh for hh in self.z if hh.refinement])
        l_variable.extend([hh for hh in self.b_iso if hh.refinement])
        l_variable.extend([hh for hh in self.occupancy if hh.refinement])
        return l_variable

    @property
    def to_cif(self):
        ls_out = []
        if self.is_defined:
            ls_out.append("loop_")
            ls_out.append("_atom_site_label")
            ls_out.append("_atom_site_type_symbol")
            ls_out.append("_atom_site_fract_x")
            ls_out.append("_atom_site_fract_y")
            ls_out.append("_atom_site_fract_z")
            ls_out.append("_atom_site_adp_type")
            ls_out.append("_atom_site_B_iso_or_equiv")
            ls_out.append("_atom_site_occupancy")
            for hh_1, hh_2, hh_3, hh_4, hh_5, hh_6, hh_7, hh_8 in zip(self.label, self.type_symbol, 
                self.x, self.y, self.z, self.adp_type, self.b_iso, self.occupancy):
                ls_out.append("{:} {:} {:} {:} {:} {:} {:} {:}".format(hh_1, hh_2, 
                        hh_3.print_with_sigma, hh_4.print_with_sigma, hh_5.print_with_sigma,
                        hh_6, hh_7.print_with_sigma, hh_8.print_with_sigma))
        return "\n".join(ls_out)

    def from_cif(self, string: str):
        cif_global = Global()
        flag = cif_global.take_from_string(string)
        if not flag:
            return False
        flag = False
        flag = cif_global.is_prefix("_atom_site")
        if flag:
            cif_loop = cif_global["_atom_site"]
            l_name = cif_loop.names
            if "_atom_site_label" in l_name:
                self.label = cif_loop["_atom_site_label"]
            if "_atom_site_type_symbol" in l_name:
                self.type_symbol = cif_loop["_atom_site_type_symbol"]
            if "_atom_site_fract_x" in l_name:
                l_val = cif_loop["_atom_site_fract_x"]
                l_fitable = []
                for val in l_val:
                    fitable = Fitable(name="_atom_site_fract_x")
                    fitable.take_it(val)
                    l_fitable.append(fitable)
                self.x = l_fitable
            if "_atom_site_fract_y" in l_name:
                l_val = cif_loop["_atom_site_fract_y"]
                l_fitable = []
                for val in l_val:
                    fitable = Fitable(name="_atom_site_fract_y")
                    fitable.take_it(val)
                    l_fitable.append(fitable)
                self.y = l_fitable
            if "_atom_site_fract_z" in l_name:
                l_val = cif_loop["_atom_site_fract_z"]
                l_fitable = []
                for val in l_val:
                    fitable = Fitable(name="_atom_site_fract_z")
                    fitable.take_it(val)
                    l_fitable.append(fitable)
                self.z = l_fitable
            if "_atom_site_adp_type" in l_name:
                self.adp_type = cif_loop["_atom_site_adp_type"]
            if "_atom_site_b_iso_or_equiv" in l_name:
                l_val = cif_loop["_atom_site_b_iso_or_equiv"]
                l_fitable = []
                for val in l_val:
                    fitable = Fitable(name="_atom_site_b_iso_or_equiv")
                    fitable.take_it(val)
                    l_fitable.append(fitable)
                self.b_iso = l_fitable
            if "_atom_site_occupancy" in l_name:
                l_val = cif_loop["_atom_site_occupancy"]
                l_fitable = []
                for val in l_val:
                    fitable = Fitable(name="_atom_site_occupancy")
                    fitable.take_it(val)
                    l_fitable.append(fitable)
                self.occupancy = l_fitable
        else:
            self.label = []
        return True

    def _load_handbook_n(self):
        f_name = os.path.join(self.__f_dir_prog, "tables", "bscat.tab")
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
        self.__handbook_nucl = ldcard

    def _load_handbook_m(self):
        f_name = os.path.join(self.__f_dir_prog, "tables", "formmag.tab")
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
        self.__handbook_mag = ldcard

    @property
    def scat_length_neutron(self):
        if self.__scat_length_neutron is None:
            l_res = [self._get_scat_length_neutron(type_symbol) for type_symbol in self.type_symbol]
            res = tuple(l_res)
            self.__scat_length_neutron = l_res 
        else:
            res = tuple(self.__scat_length_neutron)
        return res

    @property
    def b_scat(self):
        return self.scat_length_neutron

    def _get_scat_length_neutron(self, type_n):
        """
        Take scat_length_neutron
        """
        str_1 = type_n.strip().lower()
        str_1 = "".join([hh if hh.isalpha() else ' ' for hh in str_1 ]).split(" ")[0]
        l_dcard = self.__handbook_nucl 
        flag = False
        for dcard in l_dcard:
            if (dcard["type_n"].lower() == str_1):
                res = dcard["b_scat"]
                flag = True
            elif flag:
                break
        if not(flag):
            res = 0.
            self._show_message("Can not find b_scat for '{:}'.\n It is putted as 0.".format(type_n))
        return res 


    @property
    def j0j2(self):
        flag = any([self.__j0_A is None, self.__j0_a is None, self.__j0_B is None, self.__j0_b is None, 
                    self.__j0_C is None, self.__j0_c is None, self.__j0_D is None, 
                    self.__j2_A is None, self.__j2_a is None, self.__j2_B is None, self.__j2_b is None, 
                    self.__j2_C is None, self.__j2_c is None, self.__j2_D is None])

        if flag:
            l_j0_A, l_j0_a, l_j0_B, l_j0_b, l_j0_C, l_j0_c, l_j0_D = [], [], [], [], [], [], []
            l_j2_A, l_j2_a, l_j2_B, l_j2_b, l_j2_C, l_j2_c, l_j2_D = [], [], [], [], [], [], []
            for type_symbol in self.type_symbol:
                j0_A, j0_a, j0_B, j0_b, j0_C, j0_c, j0_D, j2_A, j2_a, j2_B, j2_b, j2_C, j2_c, j2_D = self._get_j0j2(type_symbol)
                l_j0_A.append(j0_A)
                l_j0_a.append(j0_a)
                l_j0_B.append(j0_B)
                l_j0_b.append(j0_b)
                l_j0_C.append(j0_C)
                l_j0_c.append(j0_c)
                l_j0_D.append(j0_D)
                l_j2_A.append(j2_A)
                l_j2_a.append(j2_a)
                l_j2_B.append(j2_B)
                l_j2_b.append(j2_b)
                l_j2_C.append(j2_C)
                l_j2_c.append(j2_c)
                l_j2_D.append(j2_D)
            self.__j0_A, self.__j0_a, self.__j0_B, self.__j0_b = tuple(l_j0_A), tuple(l_j0_a), tuple(l_j0_B), tuple(l_j0_b)
            self.__j0_C, self.__j0_c, self.__j0_D = tuple(l_j0_C), tuple(l_j0_c), tuple(l_j0_D)
            self.__j2_A, self.__j2_a, self.__j2_B, self.__j2_b = tuple(l_j2_A), tuple(l_j2_a), tuple(l_j2_B), tuple(l_j2_b)
            self.__j2_C, self.__j2_c, self.__j2_D = tuple(l_j2_C), tuple(l_j2_c), tuple(l_j2_D)
            l_res = tuple([l_j0_A, l_j0_a, l_j0_B, l_j0_b, l_j0_C, l_j0_c, l_j0_D, 
                           l_j2_A, l_j2_a, l_j2_B, l_j2_b, l_j2_C, l_j2_c, l_j2_D])
        else:
            j0_A, j0_a, j0_B, j0_b = self.__j0_A, self.__j0_a, self.__j0_B, self.__j0_b
            j0_C, j0_c, j0_D = self.__j0_C, self.__j0_c, self.__j0_D
            j2_A, j2_a, j2_B, j2_b = self.__j2_A, self.__j2_a, self.__j2_B, self.__j2_b
            j2_C, j2_c, j2_D = self.__j2_C, self.__j2_c, self.__j2_D
            l_res = tuple([j0_A, j0_a, j0_B, j0_b, j0_C, j0_c, j0_D, 
                           j2_A, j2_a, j2_B, j2_b, j2_C, j2_c, j2_D])
        return l_res

    @property
    def j0_A(self):
        if self.__j0_A is None:
            l_res = self.j0j2
            res = l_res[0]
        else:
            res = self.__j0_A
        return res
    @property
    def j0_a(self):
        if self.__j0_a is None:
            l_res = self.j0j2
            res = l_res[1]
        else:
            res = self.__j0_a
        return res
    @property
    def j0_B(self):
        if self.__j0_B is None:
            l_res = self.j0j2
            res = l_res[2]
        else:
            res = self.__j0_B
        return res
    @property
    def j0_b(self):
        if self.__j0_b is None:
            l_res = self.j0j2
            res = l_res[3]
        else:
            res = self.__j0_b
        return res
    @property
    def j0_C(self):
        if self.__j0_C is None:
            l_res = self.j0j2
            res = l_res[4]
        else:
            res = self.__j0_C
        return res
    @property
    def j0_c(self):
        if self.__j0_c is None:
            l_res = self.j0j2
            res = l_res[5]
        else:
            res = self.__j0_c
        return res
    @property
    def j0_D(self):
        if self.__j0_D is None:
            l_res = self.j0j2
            res = l_res[6]
        else:
            res = self.__j0_D
        return res
    @property
    def j2_A(self):
        if self.__j2_A is None:
            l_res = self.j0j2
            res = l_res[7]
        else:
            res = self.__j2_A
        return res
    @property
    def j2_a(self):
        if self.__j2_a is None:
            l_res = self.j0j2
            res = l_res[8]
        else:
            res = self.__j2_a
        return res
    @property
    def j2_B(self):
        if self.__j2_B is None:
            l_res = self.j0j2
            res = l_res[9]
        else:
            res = self.__j2_B
        return res
    @property
    def j2_b(self):
        if self.__j2_b is None:
            l_res = self.j0j2
            res = l_res[10]
        else:
            res = self.__j2_b
        return res
    @property
    def j2_C(self):
        if self.__j2_C is None:
            l_res = self.j0j2
            res = l_res[11]
        else:
            res = self.__j2_C
        return res
    @property
    def j2_c(self):
        if self.__j2_c is None:
            l_res = self.j0j2
            res = l_res[12]
        else:
            res = self.__j2_c
        return res
    @property
    def j2_D(self):
        if self.__j2_D is None:
            l_res = self.j0j2
            res = l_res[13]
        else:
            res = self.__j2_D
        return res


    def _get_j0j2(self, type_m):
        """
        Take coefficients for <j0> and <j2>
        """
        str_1 = type_m.strip().lower()
        str_1 = "".join([hh if hh.isalnum() else ' ' for hh in str_1 ]).split(" ")[0]

        ldcard = self.__handbook_mag 
        flag_0, flag_2 = False, False
        for dcard in ldcard:
            if ((dcard["type_m"].lower() == str_1)&(dcard["order"] == 0)):
                j0_A = dcard["A"]
                j0_a = dcard["a"]
                j0_B = dcard["B"]
                j0_b = dcard["b"]
                j0_C = dcard["C"]
                j0_c = dcard["c"]
                j0_D = dcard["D"]
                flag_0 = True
            elif ((dcard["type_m"].lower() == str_1)&(dcard["order"] == 2)):
                j2_A = dcard["A"]
                j2_a = dcard["a"]
                j2_B = dcard["B"]
                j2_b = dcard["b"]
                j2_C = dcard["C"]
                j2_c = dcard["c"]
                j2_D = dcard["D"]
                flag_2 = True
            elif (flag_0 & flag_2):
                break
        if not(flag_0):
            j0_A, j0_a, j0_B, j0_b, j0_C, j0_c, j0_D = 0., 0., 0., 0., 0., 0. ,0.
            self._show_message("Can not find coefficients <j0> for '{:}'.\n It is setted as 0.".format(type_m))
        if not(flag_2):
            j2_A, j2_a, j2_B, j2_b, j2_C, j2_c, j2_D = 0., 0., 0., 0., 0., 0. ,0.
            self._show_message("Can not find coefficients <j2> for '{:}.\n It is setted as 0.'".format(type_m))
        return j0_A, j0_a, j0_B, j0_b, j0_C, j0_c, j0_D, j2_A, j2_a, j2_B, j2_b, j2_C, j2_c, j2_D

    def _form_fract(self):
        fract = Fract(x=numpy.array(self.x, dtype=float), 
                      y=numpy.array(self.y, dtype=float), 
                      z=numpy.array(self.z, dtype=float))
        return fract

    def calc_constr_number(self, space_group):
        """
        according to table 1 in Peterse, Palm, Acta Cryst.(1966), 20, 147
        """
        l_numb = []
        for _1, _2, _3 in zip(self.x, self.y, self.z):
            x, y, z = float(_1), float(_2), float(_3)
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
            l_numb.append(numb)
        return l_numb