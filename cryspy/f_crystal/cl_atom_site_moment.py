"""
define classes to describe AtomSiteMoment (from mcif)
"""
__author__ = 'ikibalin'
__version__ = "2019_11_07"
import os
import numpy


from pycifstar import Global
from cryspy.f_common.cl_fitable import Fitable

class AtomSiteMoment(object):
    """
    This category provides a loop for presenting the magnetic moments 
    of atoms in one of several coordinate systems. This is a child category 
    of the AtomSite category, so that the magnetic moments can either be 
    listed alongside the non-magnetic atom properties in the main AtomSite loop
    (not realized) or be listed in a separate loop (realized)

    Category key: atom_site_moment_label

    Description in cif file:

    loop_                                     
    _atom_site_moment_label
    _atom_site_moment_crystalaxis_x
    _atom_site_moment_crystalaxis_y
    _atom_site_moment_crystalaxis_z
    Fe3A 4.8  0.0  0.0
    Fe3B 0.0 -4.5  0.0

    
    """    
    def __init__(self, label=[], crystalaxis_x=[], crystalaxis_y=[], crystalaxis_z=[]):
        super(AtomSiteMoment, self).__init__()

        self.__atom_site_moment_label = []
        self.__atom_site_moment_crystalaxis_x = []
        self.__atom_site_moment_crystalaxis_y = []
        self.__atom_site_moment_crystalaxis_z = []
        self.label = label
        self.crystalaxis_x = crystalaxis_x
        self.crystalaxis_y = crystalaxis_y
        self.crystalaxis_z = crystalaxis_z
        
    def __str__(self):
        ls_out = []
        ls_out.append(" label    crystalaxis_x  crystalaxis_y crystalaxis_z   ")
        for _1, _2, _3, _4 in zip(self.label, self.crystalaxis_x, self.crystalaxis_y, self.crystalaxis_z):
            ls_out.append(f" {_1:8} {_2.print_with_sigma:8} {_3.print_with_sigma:8} {_4.print_with_sigma:8}")
        return "\n".join(ls_out)

    def __repr__(self):
        ls_out = ["AtomSiteMagnetism:"]
        ls_out.append(str(self))
        return "\n".join(ls_out)

    @property
    def label(self):
        """
        This label in a unique identifier for a particular site in the asymmetric unit of the crystal unit cell

        Type: char
        """
        return tuple(self.__atom_site_moment_label)
    @label.setter
    def label(self, l_x):
        l_x_in = []
        for x in l_x:
            x_in = str(x).strip()
            l_x_in.append(x_in)
        self.__atom_site_moment_label = l_x_in
        len_x = len(l_x_in)
        len_1 = len(self.__atom_site_moment_crystalaxis_x)
        if len_1 > len_x:
            self.__atom_site_moment_crystalaxis_x = self.__atom_site_moment_crystalaxis_x[:len_x]
        elif len_1 < len_x:
            l_fitable = [Fitable(value=2., name="_atom_site_moment_crystalaxis_x") for hh in range(len_x-len_1)]
            self.__atom_site_moment_crystalaxis_x.extend(l_fitable)

        len_1 = len(self.__atom_site_moment_crystalaxis_y)
        if len_1 > len_x:
            self.__atom_site_moment_crystalaxis_y = self.__atom_site_moment_crystalaxis_y[:len_x]
        elif len_1 < len_x:
            l_fitable = [Fitable(value=2., name="_atom_site_moment_crystalaxis_y") for hh in range(len_x-len_1)]
            self.__atom_site_moment_crystalaxis_y.extend(l_fitable)

            
        len_1 = len(self.__atom_site_moment_crystalaxis_z)
        if len_1 > len_x:
            self.__atom_site_moment_crystalaxis_z = self.__atom_site_moment_crystalaxis_z[:len_x]
        elif len_1 < len_x:
            l_fitable = [Fitable(value=2., name="_atom_site_moment_crystalaxis_z") for hh in range(len_x-len_1)]
            self.__atom_site_moment_crystalaxis_z.extend(l_fitable)

    @property
    def crystalaxis_x(self):
        """
        The atom-site magnetic moment vector specified as a 
        projection onto the axes of the unit cell.

        in mu_B

        Default: 0.

        Type: float in 
        """
        return tuple(self.__atom_site_moment_crystalaxis_x)
    @crystalaxis_x.setter
    def crystalaxis_x(self, l_x):
        l_fitable = []
        for x in l_x:
            if isinstance(x, Fitable):
                x_in = x
            else:
                x_in = Fitable()
                flag = x_in.take_it(x)
            l_fitable.append(x_in)
        len_x = len(l_fitable)
        len_1 = len(self.__atom_site_moment_label)
        if len_1 < len_x:
            l_fitable = l_fitable[:len_1]
        elif len_1 > len_x:
            l_fitable.extend([Fitable(value=0., name="_atom_site_moment_crystalaxis_x") for hh in range(len_1-len_x)])
        self.__atom_site_moment_crystalaxis_x = l_fitable


    @property
    def crystalaxis_y(self):
        """
        The atom-site magnetic moment vector specified as a 
        projection onto the axes of the unit cell.

        in mu_B

        Default: 0.

        Type: float in 
        """
        return tuple(self.__atom_site_moment_crystalaxis_y)
    @crystalaxis_y.setter
    def crystalaxis_y(self, l_x):
        l_fitable = []
        for x in l_x:
            if isinstance(x, Fitable):
                x_in = x
            else:
                x_in = Fitable()
                flag = x_in.take_it(x)
            l_fitable.append(x_in)
        len_x = len(l_fitable)
        len_1 = len(self.__atom_site_moment_label)
        if len_1 < len_x:
            l_fitable = l_fitable[:len_1]
        elif len_1 > len_x:
            l_fitable.extend([Fitable(value=0., name="_atom_site_moment_crystalaxis_y") for hh in range(len_1-len_x)])
        self.__atom_site_moment_crystalaxis_y = l_fitable


    @property
    def crystalaxis_z(self):
        """
        The atom-site magnetic moment vector specified as a 
        projection onto the axes of the unit cell.

        in mu_B

        Default: 0.

        Type: float in 
        """
        return tuple(self.__atom_site_moment_crystalaxis_z)
    @crystalaxis_z.setter
    def crystalaxis_z(self, l_x):
        l_fitable = []
        for x in l_x:
            if isinstance(x, Fitable):
                x_in = x
            else:
                x_in = Fitable()
                flag = x_in.take_it(x)
            l_fitable.append(x_in)
        len_x = len(l_fitable)
        len_1 = len(self.__atom_site_moment_label)
        if len_1 < len_x:
            l_fitable = l_fitable[:len_1]
        elif len_1 > len_x:
            l_fitable.extend([Fitable(value=0., name="_atom_site_moment_crystalaxis_z") for hh in range(len_1-len_x)])
        self.__atom_site_moment_crystalaxis_z = l_fitable


    @property
    def moment(self):
        np_x = numpy.array(self.crystalaxis_x, dtype=float)
        np_y = numpy.array(self.crystalaxis_y, dtype=float)
        np_z = numpy.array(self.crystalaxis_z, dtype=float)
        np_moment = numpy.sqrt(numpy.square(np_x)+numpy.square(np_y)+numpy.square(np_z))
        return np_moment

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
        res = (any([hh.refinement for hh in self.crystalaxis_x]) | any([hh.refinement for hh in self.crystalaxis_y]) |
               any([hh.refinement for hh in self.crystalaxis_z]))
        return res


    def get_variables(self):
        l_variable = [hh for hh in self.crystalaxis_x if hh.refinement]
        l_variable.extend([hh for hh in self.crystalaxis_y if hh.refinement])
        l_variable.extend([hh for hh in self.crystalaxis_z if hh.refinement])
        return l_variable

    @property
    def to_cif(self):
        ls_out = []
        if self.is_defined:
            ls_out.append("loop_")
            ls_out.append("_atom_site_moment_label")
            ls_out.append("_atom_site_moment_crystalaxis_x")
            ls_out.append("_atom_site_moment_crystalaxis_y")
            ls_out.append("_atom_site_moment_crystalaxis_z")
            for _1, _2, _3, _4 in zip(self.label, self.crystalaxis_x, self.crystalaxis_y, self.crystalaxis_z):
                ls_out.append(f"{_1:} {_2.print_with_sigma:} {_3.print_with_sigma:} {_4.print_with_sigma:}")
        return "\n".join(ls_out)

    def from_cif(self, string: str):
        cif_global = Global()
        flag = cif_global.take_from_string(string)
        if not flag:
            return False
        flag = False
        flag = cif_global.is_prefix("_atom_site_moment")
        if flag:
            cif_loop = cif_global["_atom_site_moment"]
            l_name = cif_loop.names
            if "_atom_site_moment_label" in l_name:
                self.label = cif_loop["_atom_site_moment_label"]

            if "_atom_site_moment_crystalaxis_x" in l_name:
                l_val = cif_loop["_atom_site_moment_crystalaxis_x"]
                l_fitable = []
                for val in l_val:
                    fitable = Fitable(name="_atom_site_moment_crystalaxis_x")
                    fitable.take_it(val)
                    l_fitable.append(fitable)
                self.crystalaxis_x = l_fitable
            if "_atom_site_moment_crystalaxis_y" in l_name:
                l_val = cif_loop["_atom_site_moment_crystalaxis_y"]
                l_fitable = []
                for val in l_val:
                    fitable = Fitable(name="_atom_site_moment_crystalaxis_y")
                    fitable.take_it(val)
                    l_fitable.append(fitable)
                self.crystalaxis_y = l_fitable
            if "_atom_site_moment_crystalaxis_z" in l_name:
                l_val = cif_loop["_atom_site_moment_crystalaxis_z"]
                l_fitable = []
                for val in l_val:
                    fitable = Fitable(name="_atom_site_moment_crystalaxis_z")
                    fitable.take_it(val)
                    l_fitable.append(fitable)
                self.crystalaxis_z = l_fitable
        else:
            self.label = []

        return True

    def calc_zeeman(self, field_cryst):
        """
        !!! It is valid only for cubic crystals
        moment*sin(moment^field) at angle decreesing 
        """
        h_a, h_b, h_c = field_cryst[0], field_cryst[1], field_cryst[2]
        if abs(h_a)+abs(h_b)+abs(h_c) ==0.:
            np_val = numpy.zeros(len(self.label), dtype=float)
            return np_val
        else:
            mod_h = (abs(h_a)**2+abs(h_b)**2+abs(h_c)**2)**0.5
            h_a_n, h_b_n, h_c_n = h_a/mod_h, h_b/mod_h, h_c/mod_h
            np_moment_sq = numpy.square(self.moment)
            np_x = numpy.array(self.crystalaxis_x, dtype=float)
            np_y = numpy.array(self.crystalaxis_y, dtype=float)
            np_z = numpy.array(self.crystalaxis_z, dtype=float)
            np_val = h_a_n*np_x+h_b_n*np_y+h_c_n*np_z
            np_res = numpy.sqrt(np_moment_sq-numpy.square(np_val))
        return np_res