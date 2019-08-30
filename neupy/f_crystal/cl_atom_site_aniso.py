"""
define classes to describe AtomSiteAniso
"""
__author__ = 'ikibalin'
__version__ = "2019_08_30"
import os
import numpy


from pystar import CIFglobal
from neupy.f_common.cl_fitable import Fitable

class AtomSiteAniso(object):
    """
    Data items in the ATOM_SITE_ANISO category record details about
    the atom sites in a crystal structure, such as atomic displacement parameters.

    Data items in the ATOM_site_ANISO category record details about
    magnetic properties of the atoms that occupy the atom sites.
    
    Description in cif file:

    loop_
    _atom_site_aniso_label
    _atom_site_aniso_U_11
    _atom_site_aniso_U_22
    _atom_site_aniso_U_33
    _atom_site_aniso_U_12
    _atom_site_aniso_U_13
    _atom_site_aniso_U_23
     O1   .071(1) .076(1) .0342(9) .008(1)   .0051(9) -.0030(9) 
     C2   .060(2) .072(2) .047(1)  .002(2)   .013(1)  -.009(1)  

    
    reference: https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Catom_site.html
    """    
    def __init__(self, label=[],  
                 u_11=[], u_22=[], u_33=[],
                 u_12=[], u_13=[], u_23=[]):
        super(AtomSiteAniso, self).__init__()

        self.__atom_site_aniso_label = []
        self.__atom_site_aniso_u_11 = []
        self.__atom_site_aniso_u_12 = []
        self.__atom_site_aniso_u_13 = []
        self.__atom_site_aniso_u_22 = []
        self.__atom_site_aniso_u_23 = []
        self.__atom_site_aniso_u_33 = []

        self.label = label
        self.u_11 = u_11
        self.u_12 = u_12
        self.u_13 = u_13
        self.u_22 = u_22
        self.u_23 = u_23
        self.u_33 = u_33
        
    def __repr__(self):
        ls_out = ["AtomSiteAniso:"]
        ls_out.append(" label    u_11     u_22     u_33     u_12     u_13     u_23    ")
        for hh_1, hh_2, hh_3, hh_4, hh_5, hh_6, hh_7 in zip(
            self.label, self.u_11, self.u_22, self.u_33, self.u_12, self.u_13, self.u_23):
            ls_out.append(" {:8} {:8} {:8} {:8} {:8} {:8} {:8}".format(hh_1, 
                    hh_2.print_with_sigma, hh_3.print_with_sigma, hh_4.print_with_sigma, 
                    hh_5.print_with_sigma, hh_6.print_with_sigma, hh_7.print_with_sigma))
        return "\n".join(ls_out)

    @property
    def label(self):
        """
        Anisotropic atomic displacement parameters are looped in
        a separate list. This code must match the
        _atom_site_label of the associated atom in the atom coordinate
        list and conform with the same rules described in
        _atom_site_label.

        Type: char

        reference: https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Iatom_site_aniso_label.html
        """
        return tuple(self.__atom_site_aniso_label)
    @label.setter
    def label(self, l_x):
        l_x_in = []
        for x in l_x:
            x_in = str(x).strip()
            l_x_in.append(x_in)
        self.__atom_site_aniso_label = l_x_in
        len_x = len(l_x_in)

        len_1 = len(self.__atom_site_aniso_u_11)
        if len_1 > len_x:
            self.__atom_site_aniso_u_11 = self.__atom_site_aniso_u_11[:len_x]
        elif len_1 < len_x:
            l_fitable = [Fitable(value=0., name="_atom_site_aniso_u_11") for hh in range(len_x-len_1)]
            self.__atom_site_aniso_u_11.extend(l_fitable)

        len_1 = len(self.__atom_site_aniso_u_22)
        if len_1 > len_x:
            self.__atom_site_aniso_u_22 = self.__atom_site_aniso_u_22[:len_x]
        elif len_1 < len_x:
            l_fitable = [Fitable(value=0., name="_atom_site_aniso_u_22") for hh in range(len_x-len_1)]
            self.__atom_site_aniso_u_22.extend(l_fitable)

        len_1 = len(self.__atom_site_aniso_u_33)
        if len_1 > len_x:
            self.__atom_site_aniso_u_33 = self.__atom_site_aniso_u_33[:len_x]
        elif len_1 < len_x:
            l_fitable = [Fitable(value=0., name="_atom_site_aniso_u_33") for hh in range(len_x-len_1)]
            self.__atom_site_aniso_u_33.extend(l_fitable)

        len_1 = len(self.__atom_site_aniso_u_11)
        if len_1 > len_x:
            self.__atom_site_aniso_u_12 = self.__atom_site_aniso_u_12[:len_x]
        elif len_1 < len_x:
            l_fitable = [Fitable(value=0., name="_atom_site_aniso_u_12") for hh in range(len_x-len_1)]
            self.__atom_site_aniso_u_12.extend(l_fitable)

        len_1 = len(self.__atom_site_aniso_u_13)
        if len_1 > len_x:
            self.__atom_site_aniso_u_13 = self.__atom_site_aniso_u_13[:len_x]
        elif len_1 < len_x:
            l_fitable = [Fitable(value=0., name="_atom_site_aniso_u_13") for hh in range(len_x-len_1)]
            self.__atom_site_aniso_u_13.extend(l_fitable)

        len_1 = len(self.__atom_site_aniso_u_23)
        if len_1 > len_x:
            self.__atom_site_aniso_u_23 = self.__atom_site_aniso_u_23[:len_x]
        elif len_1 < len_x:
            l_fitable = [Fitable(value=0., name="_atom_site_aniso_u_23") for hh in range(len_x-len_1)]
            self.__atom_site_aniso_u_23.extend(l_fitable)

    @property
    def u_11(self):
        """
        These are the standard anisotropic atomic displacement
        components in angstroms squared which appear in the
        structure-factor term

        T = exp{-2pi^2^ sum~i~ [sum~j~ (U^ij^ h~i~ h~j~ a*~i~ a*~j~) ] }

        h = the Miller indices
        a* = the reciprocal-space cell lengths

        The unique elements of the real symmetric matrix are
        entered by row.

        Default: 0.

        Type: float

        reference: https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Iatom_site_aniso_U_.html
        """
        return tuple(self.__atom_site_aniso_u_11)
    @u_11.setter
    def u_11(self, l_x):
        l_fitable = []
        for x in l_x:
            if isinstance(x, Fitable):
                x_in = x
            else:
                x_in = Fitable()
                flag = x_in.take_it(x)
            l_fitable.append(x_in)
        len_x = len(l_fitable)
        len_1 = len(self.__atom_site_aniso_label)
        if len_1 < len_x:
            l_fitable = l_fitable[:len_1]
        elif len_1 > len_x:
            l_fitable.extend([Fitable(value=0., name= "_atom_site_aniso_u_11") for hh in range(len_1-len_x)])
        self.__atom_site_aniso_u_11 = l_fitable


    @property
    def u_22(self):
        """
        see help for u_11 parameter
        """
        return tuple(self.__atom_site_aniso_u_22)
    @u_22.setter
    def u_22(self, l_x):
        l_fitable = []
        for x in l_x:
            if isinstance(x, Fitable):
                x_in = x
            else:
                x_in = Fitable()
                flag = x_in.take_it(x)
            l_fitable.append(x_in)
        len_x = len(l_fitable)
        len_1 = len(self.__atom_site_aniso_label)
        if len_1 < len_x:
            l_fitable = l_fitable[:len_1]
        elif len_1 > len_x:
            l_fitable.extend([Fitable(value=0., name= "_atom_site_aniso_u_22") for hh in range(len_1-len_x)])
        self.__atom_site_aniso_u_22 = l_fitable


    @property
    def u_33(self):
        """
        see help for u_11 parameter
        """
        return tuple(self.__atom_site_aniso_u_33)
    @u_33.setter
    def u_33(self, l_x):
        l_fitable = []
        for x in l_x:
            if isinstance(x, Fitable):
                x_in = x
            else:
                x_in = Fitable()
                flag = x_in.take_it(x)
            l_fitable.append(x_in)
        len_x = len(l_fitable)
        len_1 = len(self.__atom_site_aniso_label)
        if len_1 < len_x:
            l_fitable = l_fitable[:len_1]
        elif len_1 > len_x:
            l_fitable.extend([Fitable(value=0., name= "_atom_site_aniso_u_33") for hh in range(len_1-len_x)])
        self.__atom_site_aniso_u_33 = l_fitable


    @property
    def u_12(self):
        """
        see help for u_11 parameter
        """
        return tuple(self.__atom_site_aniso_u_12)
    @u_12.setter
    def u_12(self, l_x):
        l_fitable = []
        for x in l_x:
            if isinstance(x, Fitable):
                x_in = x
            else:
                x_in = Fitable()
                flag = x_in.take_it(x)
            l_fitable.append(x_in)
        len_x = len(l_fitable)
        len_1 = len(self.__atom_site_aniso_label)
        if len_1 < len_x:
            l_fitable = l_fitable[:len_1]
        elif len_1 > len_x:
            l_fitable.extend([Fitable(value=0., name= "_atom_site_aniso_u_12") for hh in range(len_1-len_x)])
        self.__atom_site_aniso_u_12 = l_fitable


    @property
    def u_13(self):
        """
        see help for u_11 parameter
        """
        return tuple(self.__atom_site_aniso_u_13)
    @u_13.setter
    def u_13(self, l_x):
        l_fitable = []
        for x in l_x:
            if isinstance(x, Fitable):
                x_in = x
            else:
                x_in = Fitable()
                flag = x_in.take_it(x)
            l_fitable.append(x_in)
        len_x = len(l_fitable)
        len_1 = len(self.__atom_site_aniso_label)
        if len_1 < len_x:
            l_fitable = l_fitable[:len_1]
        elif len_1 > len_x:
            l_fitable.extend([Fitable(value=0., name= "_atom_site_aniso_u_13") for hh in range(len_1-len_x)])
        self.__atom_site_aniso_u_13 = l_fitable


    @property
    def u_23(self):
        """
        see help for u_11 parameter
        """
        return tuple(self.__atom_site_aniso_u_23)
    @u_23.setter
    def u_23(self, l_x):
        l_fitable = []
        for x in l_x:
            if isinstance(x, Fitable):
                x_in = x
            else:
                x_in = Fitable()
                flag = x_in.take_it(x)
            l_fitable.append(x_in)
        len_x = len(l_fitable)
        len_1 = len(self.__atom_site_aniso_label)
        if len_1 < len_x:
            l_fitable = l_fitable[:len_1]
        elif len_1 > len_x:
            l_fitable.extend([Fitable(value=0., name= "_atom_site_aniso_u_23") for hh in range(len_1-len_x)])
        self.__atom_site_aniso_u_23 = l_fitable

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
        res = (any([hh.refinement for hh in self.u_11]) | 
               any([hh.refinement for hh in self.u_22]) |
               any([hh.refinement for hh in self.u_33]) |
               any([hh.refinement for hh in self.u_12]) |
               any([hh.refinement for hh in self.u_13]) |
               any([hh.refinement for hh in self.u_23]))
        return res


    def get_variables(self):
        l_variable = [hh for hh in self.u_11 if hh.refinement]
        l_variable.extend([hh for hh in self.u_22 if hh.refinement])
        l_variable.extend([hh for hh in self.u_33 if hh.refinement])
        l_variable.extend([hh for hh in self.u_12 if hh.refinement])
        l_variable.extend([hh for hh in self.u_13 if hh.refinement])
        l_variable.extend([hh for hh in self.u_23 if hh.refinement])

        return l_variable

    @property
    def to_cif(self):
        ls_out = []
        if self.is_defined:
            ls_out.append("loop_")
            ls_out.append("_atom_site_aniso_label")
            ls_out.append("_atom_site_aniso_u_11")
            ls_out.append("_atom_site_aniso_u_22")
            ls_out.append("_atom_site_aniso_u_33")
            ls_out.append("_atom_site_aniso_u_12")
            ls_out.append("_atom_site_aniso_u_13")
            ls_out.append("_atom_site_aniso_u_23")
            for hh_1, hh_2, hh_3, hh_4, hh_5, hh_6, hh_7 in zip(self.label, 
                self.u_11, self.u_22, self.u_33, self.u_12, self.u_13, self.u_23):
                ls_out.append("{:} {:} {:} {:} {:} {:} {:}".format(hh_1, 
                        hh_2.print_with_sigma, hh_3.print_with_sigma, hh_4.print_with_sigma, 
                        hh_5.print_with_sigma, hh_6.print_with_sigma, hh_7.print_with_sigma))
        return "\n".join(ls_out)

    def from_cif(self, string: str):
        cif_global = CIFglobal()
        flag = cif_global.take_from_string(string)
        if not flag:
            return False
        flag = False
        flag = cif_global.is_prefix("_atom_site_aniso")
        if flag:
            cif_loop = cif_global["_atom_site_aniso"]
            l_name = cif_loop.names
            if "_atom_site_aniso_label" in l_name:
                self.label = cif_loop["_atom_site_aniso_label"]

            if "_atom_site_aniso_u_11" in l_name:
                l_val = cif_loop["_atom_site_aniso_u_11"]
                l_fitable = []
                for val in l_val:
                    fitable = Fitable(name="_atom_site_aniso_u_11")
                    fitable.take_it(val)
                    l_fitable.append(fitable)
                self.u_11 = l_fitable
            if "_atom_site_aniso_u_22" in l_name:
                l_val = cif_loop["_atom_site_aniso_u_22"]
                l_fitable = []
                for val in l_val:
                    fitable = Fitable(name="_atom_site_aniso_u_22")
                    fitable.take_it(val)
                    l_fitable.append(fitable)
                self.u_22 = l_fitable
            if "_atom_site_aniso_u_33" in l_name:
                l_val = cif_loop["_atom_site_aniso_u_33"]
                l_fitable = []
                for val in l_val:
                    fitable = Fitable(name="_atom_site_aniso_u_33")
                    fitable.take_it(val)
                    l_fitable.append(fitable)
                self.u_33 = l_fitable
            if "_atom_site_aniso_u_12" in l_name:
                l_val = cif_loop["_atom_site_aniso_u_12"]
                l_fitable = []
                for val in l_val:
                    fitable = Fitable(name="_atom_site_aniso_u_12")
                    fitable.take_it(val)
                    l_fitable.append(fitable)
                self.u_12 = l_fitable
            if "_atom_site_aniso_u_13" in l_name:
                l_val = cif_loop["_atom_site_aniso_u_13"]
                l_fitable = []
                for val in l_val:
                    fitable = Fitable(name="_atom_site_aniso_u_13")
                    fitable.take_it(val)
                    l_fitable.append(fitable)
                self.u_13 = l_fitable
            if "_atom_site_aniso_u_23" in l_name:
                l_val = cif_loop["_atom_site_aniso_u_23"]
                l_fitable = []
                for val in l_val:
                    fitable = Fitable(name="_atom_site_aniso_u_23")
                    fitable.take_it(val)
                    l_fitable.append(fitable)
                self.u_23 = l_fitable
        else:
            self.label = []

        return True

