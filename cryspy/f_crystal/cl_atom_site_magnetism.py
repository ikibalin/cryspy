"""
define classes to describe AtomSiteMagnetism
"""
__author__ = 'ikibalin'
__version__ = "2019_08_29"
import os
import numpy


from pycifstar import Global
from cryspy.f_common.cl_fitable import Fitable

class AtomSiteMagnetism(object):
    """
    Data items in the ATOM_SITE_MAGNETISM category record details about
    magnetic properties of the atoms that occupy the atom sites.
    
    Description in cif file:

    loop_                                     
    _atom_site_magnetism_label
    _atom_site_magnetism_lande
    _atom_site_magnetism_kappa
    Fe3A 2.0 1.0
    Fe3B 2.0 1.0

    
    """    
    def __init__(self, label=[], lande=[], kappa=[]):
        super(AtomSiteMagnetism, self).__init__()

        self.__atom_site_magnetism_label = []
        self.__atom_site_magnetism_lande = []
        self.__atom_site_magnetism_kappa = []
        self.label = label
        self.lande = lande
        self.kappa = kappa
        
    def __repr__(self):
        ls_out = ["AtomSiteMagnetism:"]
        ls_out.append(str(self))
        return "\n".join(ls_out)        
        
    def __str__(self):
        ls_out = []
        ls_out.append(" label    kappa    lande   ")
        for hh_1, hh_2, hh_3 in zip(self.label, self.kappa, self.lande):
            ls_out.append(" {:8} {:8} {:8}".format(hh_1, hh_2.print_with_sigma, hh_3.print_with_sigma))
        return "\n".join(ls_out)

    @property
    def label(self):
        """
        The _atom_site_magnetism_label is a unique identifier for a particular site
        in the crystal. 

        Type: char
        """
        return tuple(self.__atom_site_magnetism_label)
    @label.setter
    def label(self, l_x):
        l_x_in = []
        for x in l_x:
            x_in = str(x).strip()
            l_x_in.append(x_in)
        self.__atom_site_magnetism_label = l_x_in
        len_x = len(l_x_in)
        len_1 = len(self.__atom_site_magnetism_lande)
        if len_1 > len_x:
            self.__atom_site_magnetism_lande = self.__atom_site_magnetism_lande[:len_x]
        elif len_1 < len_x:
            l_fitable = [Fitable(value=2., name="_atom_site_magnetism_lande") for hh in range(len_x-len_1)]
            self.__atom_site_magnetism_lande.extend(l_fitable)

        len_1 = len(self.__atom_site_magnetism_kappa)
        if len_1 > len_x:
            self.__atom_site_magnetism_kappa = self.__atom_site_magnetism_kappa[:len_x]
        elif len_1 < len_x:
            l_fitable = [Fitable(value=1., name="_atom_site_magnetism_kappa") for hh in range(len_x-len_1)]
            self.__atom_site_magnetism_kappa.extend(l_fitable)

    @property
    def lande(self):
        """
        Lande factor. 

        Default: 2

        Type: float
        """
        return tuple(self.__atom_site_magnetism_lande)
    @lande.setter
    def lande(self, l_x):
        l_fitable = []
        for x in l_x:
            if isinstance(x, Fitable):
                x_in = x
            else:
                x_in = Fitable()
                flag = x_in.take_it(x)
            l_fitable.append(x_in)
        len_x = len(l_fitable)
        len_1 = len(self.__atom_site_magnetism_label)
        if len_1 < len_x:
            l_fitable = l_fitable[:len_1]
        elif len_1 > len_x:
            l_fitable.extend([Fitable(value=2., name="_atom_site_magnetism_lande") for hh in range(len_1-len_x)])
        self.__atom_site_magnetism_lande = l_fitable

    @property
    def kappa(self):
        """
        kappa 

        Default: 1

        Type: float
        """
        return tuple(self.__atom_site_magnetism_kappa)
    @kappa.setter
    def kappa(self, l_x):
        l_fitable = []
        for x in l_x:
            if isinstance(x, Fitable):
                x_in = x
            else:
                x_in = Fitable()
                flag = x_in.take_it(x)
            l_fitable.append(x_in)
        len_x = len(l_fitable)
        len_1 = len(self.__atom_site_magnetism_label)
        if len_1 < len_x:
            l_fitable = l_fitable[:len_1]
        elif len_1 > len_x:
            l_fitable.extend([Fitable(value=1., name="_atom_site_magnetism_kappa") for hh in range(len_1-len_x)])
        self.__atom_site_magnetism_kappa = l_fitable

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
        res = (any([hh.refinement for hh in self.lande]) | any([hh.refinement for hh in self.kappa]))
        return res


    def get_variables(self):
        l_variable = [hh for hh in self.lande if hh.refinement]
        l_variable.extend([hh for hh in self.kappa if hh.refinement])

        return l_variable

    @property
    def to_cif(self):
        ls_out = []
        if self.is_defined:
            ls_out.append("loop_")
            ls_out.append("_atom_site_magnetism_label")
            ls_out.append("_atom_site_magnetism_lande")
            ls_out.append("_atom_site_magnetism_kappa")
            for hh_1, hh_2, hh_3 in zip(self.label, self.lande, self.kappa):
                ls_out.append("{:} {:} {:}".format(hh_1, 
                                                   hh_2.print_with_sigma, 
                                                   hh_3.print_with_sigma))
        return "\n".join(ls_out)

    def from_cif(self, string: str):
        cif_global = Global()
        flag = cif_global.take_from_string(string)
        if not flag:
            return False
        flag = False
        flag = cif_global.is_prefix("_atom_site_magnetism")
        if flag:
            cif_loop = cif_global["_atom_site_magnetism"]
            l_name = cif_loop.names
            if "_atom_site_magnetism_label" in l_name:
                self.label = cif_loop["_atom_site_magnetism_label"]

            if "_atom_site_magnetism_lande" in l_name:
                l_val = cif_loop["_atom_site_magnetism_lande"]
                l_fitable = []
                for val in l_val:
                    fitable = Fitable(name="_atom_site_magnetism_lande")
                    fitable.take_it(val)
                    l_fitable.append(fitable)
                self.lande = l_fitable
            if "_atom_site_magnetism_kappa" in l_name:
                l_val = cif_loop["_atom_site_magnetism_kappa"]
                l_fitable = []
                for val in l_val:
                    fitable = Fitable(name="_atom_site_magnetism_kappa")
                    fitable.take_it(val)
                    l_fitable.append(fitable)
                self.kappa = l_fitable
        else:
            self.label = []

        return True

    def _form_lande_kappa(self, atom_site):
        label = numpy.array(atom_site.label, dtype=str)
        label_magnetism = numpy.array(self.label, dtype=str)
        if not(set(label_magnetism).issubset(set(label))):
            self._show_message("Unknown 'aniso_label'")
            return False

        np_index = numpy.array([int(numpy.argwhere(label==hh)[0]) for hh in label_magnetism], dtype=int)

        lande_in = 2.*numpy.ones(label.shape, dtype=float)
        kappa_in = numpy.ones(label.shape, dtype=float)

        lande = numpy.array(self.lande, dtype=float)
        kappa = numpy.array(self.kappa, dtype=float)

        lande_in[np_index], kappa_in[np_index] = lande, kappa

        return lande_in, kappa_in