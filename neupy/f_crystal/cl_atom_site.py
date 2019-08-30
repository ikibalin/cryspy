"""
define classes to describe AtomSite
"""
__author__ = 'ikibalin'
__version__ = "2019_08_30"
import os
import numpy


from pystar import CIFglobal
from neupy.f_common.cl_fitable import Fitable

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
                       adp_type=[], b_iso=[], occupancy=[]):
        super(AtomSite, self).__init__()

        self.__atom_site_label = []
        self.__atom_site_type_symbol = []
        self.__atom_site_fract_x = []
        self.__atom_site_fract_y = []
        self.__atom_site_fract_z = []
        self.__atom_site_adp_type = []
        self.__atom_site_b_iso_or_equiv = []
        self.__atom_site_occupancy = []

        self.label = label
        self.type_symbol = type_symbol
        self.x = x
        self.y = y
        self.z = z
        self.adp_type = adp_type
        self.b_iso = b_iso
        self.occupancy = occupancy
        
    def __repr__(self):
        ls_out = ["AtomSite:"]

                       
        ls_out.append(" label type_symbol x        y        z        adp_type b_iso    occupancy")
        for hh_1, hh_2, hh_3, hh_4, hh_5, hh_6, hh_7, hh_8 in zip(
            self.label, self.type_symbol, self.x, self.y, self.z, self.adp_type, self.b_iso, self.occupancy):
            ls_out.append(" {:8} {:8} {:8} {:8} {:8} {:8} {:8} {:8}".format(hh_1, hh_2,
                    hh_3.print_with_sigma, hh_4.print_with_sigma, hh_5.print_with_sigma,
                    hh_6, hh_7.print_with_sigma, hh_8.print_with_sigma))
        return "\n".join(ls_out)

    @property
    def label(self):
        """
        The _atom_site_label is a unique identifier for a particular site
        in the crystal. 

        Type: char

        reference: https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Iatom_site_label.html
        """
        return self.__atom_site_label
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
        return self.__atom_site_type_symbol
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


    @property
    def x(self):
        """
        Atom-site coordinates as fractions of the _cell_length_ values.

        Default: 0.

        Type: float

        reference: https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Iatom_site_fract_.html
        """
        return self.__atom_site_fract_x
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
        return self.__atom_site_fract_y
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
        return self.__atom_site_fract_z
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
        return self.__atom_site_adp_type
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
        return self.__atom_site_b_iso_or_equiv
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
        return self.__atom_site_occupancy
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
        cif_global = CIFglobal()
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

