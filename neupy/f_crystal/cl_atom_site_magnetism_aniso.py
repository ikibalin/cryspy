"""
define classes to describe AtomSiteMagnetismAniso
"""
__author__ = 'ikibalin'
__version__ = "2019_08_29"
import os
import numpy


from pystar import CIFglobal
from neupy.f_common.cl_fitable import Fitable

class AtomSiteMagnetismAniso(object):
    """
    Data items in the ATOM_SITE_MAGNETISM_ANISO category record details about
    magnetic properties of the atoms that occupy the atom sites.
    
    Description in cif file:

    loop_  
    _atom_site_magnetism_aniso_label
    _atom_site_magnetism_aniso_chi_type
    _atom_site_magnetism_aniso_chi_11
    _atom_site_magnetism_aniso_chi_12
    _atom_site_magnetism_aniso_chi_13
    _atom_site_magnetism_aniso_chi_22
    _atom_site_magnetism_aniso_chi_23
    _atom_site_magnetism_aniso_chi_33
     Fe3A cani -3.468(74) 0.0 0.0 -3.468 0.0 -3.468
     Fe3B cani 3.041      0.0 0.0  3.041 0.0  3.041
    
    """    
    def __init__(self, label=[], chi_type=[], 
                 chi_11=[], chi_22=[], chi_33=[],
                 chi_12=[], chi_13=[], chi_23=[]):
        super(AtomSiteMagnetismAniso, self).__init__()

        self.__atom_site_magnetism_aniso_label = []
        self.__atom_site_magnetism_aniso_chi_type = []
        self.__atom_site_magnetism_aniso_chi_11 = []
        self.__atom_site_magnetism_aniso_chi_12 = []
        self.__atom_site_magnetism_aniso_chi_13 = []
        self.__atom_site_magnetism_aniso_chi_22 = []
        self.__atom_site_magnetism_aniso_chi_23 = []
        self.__atom_site_magnetism_aniso_chi_33 = []

        self.label = label
        self.chi_type = chi_type
        self.chi_11 = chi_11
        self.chi_12 = chi_12
        self.chi_13 = chi_13
        self.chi_22 = chi_22
        self.chi_23 = chi_23
        self.chi_33 = chi_33
        
    def __repr__(self):
        ls_out = ["AtomSiteMagnetismAniso:"]
        ls_out.append(" label    chi_type chi_11   chi_22   chi_33   chi_12   chi_13   chi_23   ")
        for hh_1, hh_2, hh_3, hh_4, hh_5, hh_6, hh_7, hh_8 in zip(
            self.label, self.chi_type, self.chi_11, self.chi_22, self.chi_33, self.chi_12, self.chi_13, self.chi_23):
            ls_out.append(" {:8} {:8} {:8} {:8} {:8} {:8} {:8} {:8}".format(hh_1, hh_2, 
                    hh_3.print_with_sigma, hh_4.print_with_sigma, hh_5.print_with_sigma, 
                    hh_6.print_with_sigma, hh_7.print_with_sigma, hh_8.print_with_sigma))
        return "\n".join(ls_out)

    @property
    def label(self):
        """
        The _atom_site_magnetism_aniso_label is a unique identifier for a particular site
        in the crystal. 

        Type: char
        """
        return self.__atom_site_magnetism_aniso_label
    @label.setter
    def label(self, l_x):
        l_x_in = []
        for x in l_x:
            x_in = str(x).strip()
            l_x_in.append(x_in)
        self.__atom_site_magnetism_aniso_label = l_x_in
        len_x = len(l_x_in)

        len_1 = len(self.__atom_site_magnetism_aniso_chi_type)
        if len_1 > len_x:
            self.__atom_site_magnetism_aniso_chi_type = self.__atom_site_magnetism_aniso_chi_type[:len_x]
        elif len_1 < len_x:
            l_fitable = ["Fe3+" for hh in range(len_x-len_1)]
            self.__atom_site_magnetism_aniso_chi_type.extend(l_fitable)

        len_1 = len(self.__atom_site_magnetism_aniso_chi_11)
        if len_1 > len_x:
            self.__atom_site_magnetism_aniso_chi_11 = self.__atom_site_magnetism_aniso_chi_11[:len_x]
        elif len_1 < len_x:
            l_fitable = [Fitable(value=0., name="__atom_site_magnetism_aniso_chi_11") for hh in range(len_x-len_1)]
            self.__atom_site_magnetism_aniso_chi_11.extend(l_fitable)

        len_1 = len(self.__atom_site_magnetism_aniso_chi_22)
        if len_1 > len_x:
            self.__atom_site_magnetism_aniso_chi_22 = self.__atom_site_magnetism_aniso_chi_22[:len_x]
        elif len_1 < len_x:
            l_fitable = [Fitable(value=0., name="__atom_site_magnetism_aniso_chi_22") for hh in range(len_x-len_1)]
            self.__atom_site_magnetism_aniso_chi_22.extend(l_fitable)

        len_1 = len(self.__atom_site_magnetism_aniso_chi_33)
        if len_1 > len_x:
            self.__atom_site_magnetism_aniso_chi_33 = self.__atom_site_magnetism_aniso_chi_33[:len_x]
        elif len_1 < len_x:
            l_fitable = [Fitable(value=0., name="__atom_site_magnetism_aniso_chi_33") for hh in range(len_x-len_1)]
            self.__atom_site_magnetism_aniso_chi_33.extend(l_fitable)

        len_1 = len(self.__atom_site_magnetism_aniso_chi_11)
        if len_1 > len_x:
            self.__atom_site_magnetism_aniso_chi_12 = self.__atom_site_magnetism_aniso_chi_12[:len_x]
        elif len_1 < len_x:
            l_fitable = [Fitable(value=0., name="__atom_site_magnetism_aniso_chi_12") for hh in range(len_x-len_1)]
            self.__atom_site_magnetism_aniso_chi_12.extend(l_fitable)

        len_1 = len(self.__atom_site_magnetism_aniso_chi_13)
        if len_1 > len_x:
            self.__atom_site_magnetism_aniso_chi_13 = self.__atom_site_magnetism_aniso_chi_13[:len_x]
        elif len_1 < len_x:
            l_fitable = [Fitable(value=0., name="__atom_site_magnetism_aniso_chi_13") for hh in range(len_x-len_1)]
            self.__atom_site_magnetism_aniso_chi_13.extend(l_fitable)

        len_1 = len(self.__atom_site_magnetism_aniso_chi_23)
        if len_1 > len_x:
            self.__atom_site_magnetism_aniso_chi_23 = self.__atom_site_magnetism_aniso_chi_23[:len_x]
        elif len_1 < len_x:
            l_fitable = [Fitable(value=0., name="__atom_site_magnetism_aniso_chi_23") for hh in range(len_x-len_1)]
            self.__atom_site_magnetism_aniso_chi_23.extend(l_fitable)

    @property
    def chi_type(self):
        """
        Chi type. 

        Default: Fe3+

        Type: char
        """
        return self.__atom_site_magnetism_aniso_chi_type
    @chi_type.setter
    def chi_type(self, l_x):
        l_fitable = []
        for x in l_x:
            x_in = str(x).strip()
            l_fitable.append(x_in)
        len_x = len(l_fitable)
        len_1 = len(self.__atom_site_magnetism_aniso_label)
        if len_1 < len_x:
            l_fitable = l_fitable[:len_1]
        elif len_1 > len_x:
            l_fitable.extend(["Fe3+" for hh in range(len_1-len_x)])
        self.__atom_site_magnetism_aniso_chi_type = l_fitable

    @property
    def chi_11(self):
        """
        chi_11

        Default: 0.

        Type: float
        """
        return self.__atom_site_magnetism_aniso_chi_11
    @chi_11.setter
    def chi_11(self, l_x):
        l_fitable = []
        for x in l_x:
            x_in = Fitable()
            flag = x_in.take_it(x)
            l_fitable.append(x_in)
        len_x = len(l_fitable)
        len_1 = len(self.__atom_site_magnetism_aniso_label)
        if len_1 < len_x:
            l_fitable = l_fitable[:len_1]
        elif len_1 > len_x:
            l_fitable.extend([Fitable(value=0., name="__atom_site_magnetism_aniso_chi_11") for hh in range(len_1-len_x)])
        self.__atom_site_magnetism_aniso_chi_11 = l_fitable


    @property
    def chi_22(self):
        """
        chi_22

        Default: 0.

        Type: float
        """
        return self.__atom_site_magnetism_aniso_chi_22
    @chi_22.setter
    def chi_22(self, l_x):
        l_fitable = []
        for x in l_x:
            x_in = Fitable()
            flag = x_in.take_it(x)
            l_fitable.append(x_in)
        len_x = len(l_fitable)
        len_1 = len(self.__atom_site_magnetism_aniso_label)
        if len_1 < len_x:
            l_fitable = l_fitable[:len_1]
        elif len_1 > len_x:
            l_fitable.extend([Fitable(value=0., name="__atom_site_magnetism_aniso_chi_22") for hh in range(len_1-len_x)])
        self.__atom_site_magnetism_aniso_chi_22 = l_fitable


    @property
    def chi_33(self):
        """
        chi_33

        Default: 0.

        Type: float
        """
        return self.__atom_site_magnetism_aniso_chi_33
    @chi_33.setter
    def chi_33(self, l_x):
        l_fitable = []
        for x in l_x:
            x_in = Fitable()
            flag = x_in.take_it(x)
            l_fitable.append(x_in)
        len_x = len(l_fitable)
        len_1 = len(self.__atom_site_magnetism_aniso_label)
        if len_1 < len_x:
            l_fitable = l_fitable[:len_1]
        elif len_1 > len_x:
            l_fitable.extend([Fitable(value=0., name="__atom_site_magnetism_aniso_chi_33") for hh in range(len_1-len_x)])
        self.__atom_site_magnetism_aniso_chi_33 = l_fitable


    @property
    def chi_12(self):
        """
        chi_12

        Default: 0.

        Type: float
        """
        return self.__atom_site_magnetism_aniso_chi_12
    @chi_12.setter
    def chi_12(self, l_x):
        l_fitable = []
        for x in l_x:
            x_in = Fitable()
            flag = x_in.take_it(x)
            l_fitable.append(x_in)
        len_x = len(l_fitable)
        len_1 = len(self.__atom_site_magnetism_aniso_label)
        if len_1 < len_x:
            l_fitable = l_fitable[:len_1]
        elif len_1 > len_x:
            l_fitable.extend([Fitable(value=0., name="__atom_site_magnetism_aniso_chi_12") for hh in range(len_1-len_x)])
        self.__atom_site_magnetism_aniso_chi_12 = l_fitable


    @property
    def chi_13(self):
        """
        chi_13

        Default: 0.

        Type: float
        """
        return self.__atom_site_magnetism_aniso_chi_13
    @chi_13.setter
    def chi_13(self, l_x):
        l_fitable = []
        for x in l_x:
            x_in = Fitable()
            flag = x_in.take_it(x)
            l_fitable.append(x_in)
        len_x = len(l_fitable)
        len_1 = len(self.__atom_site_magnetism_aniso_label)
        if len_1 < len_x:
            l_fitable = l_fitable[:len_1]
        elif len_1 > len_x:
            l_fitable.extend([Fitable(value=0., name="__atom_site_magnetism_aniso_chi_13") for hh in range(len_1-len_x)])
        self.__atom_site_magnetism_aniso_chi_13 = l_fitable


    @property
    def chi_23(self):
        """
        chi_23

        Default: 0.

        Type: float
        """
        return self.__atom_site_magnetism_aniso_chi_23
    @chi_23.setter
    def chi_23(self, l_x):
        l_fitable = []
        for x in l_x:
            x_in = Fitable()
            flag = x_in.take_it(x)
            l_fitable.append(x_in)
        len_x = len(l_fitable)
        len_1 = len(self.__atom_site_magnetism_aniso_label)
        if len_1 < len_x:
            l_fitable = l_fitable[:len_1]
        elif len_1 > len_x:
            l_fitable.extend([Fitable(value=0., name="__atom_site_magnetism_aniso_chi_23") for hh in range(len_1-len_x)])
        self.__atom_site_magnetism_aniso_chi_23 = l_fitable

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
        res = (any([hh.refinement for hh in self.chi_11]) | 
               any([hh.refinement for hh in self.chi_22]) |
               any([hh.refinement for hh in self.chi_33]) |
               any([hh.refinement for hh in self.chi_12]) |
               any([hh.refinement for hh in self.chi_13]) |
               any([hh.refinement for hh in self.chi_23]))
        return res


    def get_variables(self):
        l_variable = [hh for hh in self.chi_11 if hh.refinement]
        l_variable.extend([hh for hh in self.chi_22 if hh.refinement])
        l_variable.extend([hh for hh in self.chi_33 if hh.refinement])
        l_variable.extend([hh for hh in self.chi_12 if hh.refinement])
        l_variable.extend([hh for hh in self.chi_13 if hh.refinement])
        l_variable.extend([hh for hh in self.chi_23 if hh.refinement])

        return l_variable

    @property
    def to_cif(self):
        ls_out = []
        if self.is_defined:
            ls_out.append("loop_")
            ls_out.append("_atom_site_magnetism_aniso_label")
            ls_out.append("_atom_site_magnetism_aniso_chi_type")
            ls_out.append("_atom_site_magnetism_aniso_chi_11")
            ls_out.append("_atom_site_magnetism_aniso_chi_22")
            ls_out.append("_atom_site_magnetism_aniso_chi_33")
            ls_out.append("_atom_site_magnetism_aniso_chi_12")
            ls_out.append("_atom_site_magnetism_aniso_chi_13")
            ls_out.append("_atom_site_magnetism_aniso_chi_23")
            for hh_1, hh_2, hh_3, hh_4, hh_5, hh_6, hh_7, hh_8 in zip(self.label, self.chi_type, 
                self.chi_11, self.chi_22, self.chi_33, self.chi_12, self.chi_13, self.chi_23):
                ls_out.append("{:} {:} {:} {:} {:} {:} {:} {:}".format(hh_1, hh_2, 
                        hh_3.print_with_sigma, hh_4.print_with_sigma, hh_5.print_with_sigma,
                        hh_6.print_with_sigma, hh_7.print_with_sigma, hh_8.print_with_sigma))
        return "\n".join(ls_out)

    def from_cif(self, string: str):
        cif_global = CIFglobal()
        flag = cif_global.take_from_string(string)
        if not flag:
            return False
        flag = False
        flag = cif_global.is_prefix("_atom_site_magnetism_aniso")
        if flag:
            cif_loop = cif_global["_atom_site_magnetism_aniso"]
            l_name = cif_loop.names
            if "_atom_site_magnetism_aniso_label" in l_name:
                self.label = cif_loop["_atom_site_magnetism_aniso_label"]

            if "_atom_site_magnetism_aniso_chi_type" in l_name:
                self.chi_type = cif_loop["_atom_site_magnetism_aniso_chi_type"]
            if "_atom_site_magnetism_aniso_chi_11" in l_name:
                l_val = cif_loop["_atom_site_magnetism_aniso_chi_11"]
                l_fitable = []
                for val in l_val:
                    fitable = Fitable(name="_atom_site_magnetism_aniso_chi_11")
                    fitable.take_it(val)
                    l_fitable.append(fitable)
                self.chi_11 = l_fitable
            if "_atom_site_magnetism_aniso_chi_22" in l_name:
                l_val = cif_loop["_atom_site_magnetism_aniso_chi_22"]
                l_fitable = []
                for val in l_val:
                    fitable = Fitable(name="_atom_site_magnetism_aniso_chi_22")
                    fitable.take_it(val)
                    l_fitable.append(fitable)
                self.chi_22 = l_fitable
            if "_atom_site_magnetism_aniso_chi_33" in l_name:
                l_val = cif_loop["_atom_site_magnetism_aniso_chi_33"]
                l_fitable = []
                for val in l_val:
                    fitable = Fitable(name="_atom_site_magnetism_aniso_chi_33")
                    fitable.take_it(val)
                    l_fitable.append(fitable)
                self.chi_33 = l_fitable
            if "_atom_site_magnetism_aniso_chi_12" in l_name:
                l_val = cif_loop["_atom_site_magnetism_aniso_chi_12"]
                l_fitable = []
                for val in l_val:
                    fitable = Fitable(name="_atom_site_magnetism_aniso_chi_12")
                    fitable.take_it(val)
                    l_fitable.append(fitable)
                self.chi_12 = l_fitable
            if "_atom_site_magnetism_aniso_chi_13" in l_name:
                l_val = cif_loop["_atom_site_magnetism_aniso_chi_13"]
                l_fitable = []
                for val in l_val:
                    fitable = Fitable(name="_atom_site_magnetism_aniso_chi_13")
                    fitable.take_it(val)
                    l_fitable.append(fitable)
                self.chi_13 = l_fitable
            if "_atom_site_magnetism_aniso_chi_23" in l_name:
                l_val = cif_loop["_atom_site_magnetism_aniso_chi_23"]
                l_fitable = []
                for val in l_val:
                    fitable = Fitable(name="_atom_site_magnetism_aniso_chi_23")
                    fitable.take_it(val)
                    l_fitable.append(fitable)
                self.chi_23 = l_fitable
        else:
            self.label = []

        return True

