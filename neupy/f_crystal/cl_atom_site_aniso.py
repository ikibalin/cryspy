"""
define classes to describe AtomSiteAniso
"""
__author__ = 'ikibalin'
__version__ = "2019_08_30"
import os
import numpy


from pystar import CIFglobal
from neupy.f_common.cl_fitable import Fitable
from neupy.f_crystal.cl_adp import ADP

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

    def _form_adp(self, atom_site):
        label = numpy.array(atom_site.label, dtype=str)
        b_iso = numpy.array(atom_site.b_iso, dtype=float)
        adp_type = numpy.array(atom_site.adp_type, dtype=str)

        b_iso_in = numpy.zeros(b_iso.shape, dtype=float)
        b_iso_in = numpy.where(adp_type=="uiso", b_iso, b_iso_in)
        label_aniso = numpy.array(self.label, dtype=str)
        if not(set(label_aniso).issubset(set(label))):
            self._show_message("Unknown 'aniso_label'")
            return False
        
        u_11 = numpy.array(self.u_11, dtype=float)
        u_22 = numpy.array(self.u_22, dtype=float)
        u_33 = numpy.array(self.u_33, dtype=float)
        u_12 = numpy.array(self.u_12, dtype=float)
        u_13 = numpy.array(self.u_13, dtype=float)
        u_23 = numpy.array(self.u_23, dtype=float)


        np_index = numpy.array([int(numpy.argwhere(label==hh)[0]) for hh in label_aniso], dtype=int)

        u_11_in = numpy.zeros(label.shape, dtype=float)
        u_22_in = numpy.zeros(label.shape, dtype=float)
        u_33_in = numpy.zeros(label.shape, dtype=float)
        u_12_in = numpy.zeros(label.shape, dtype=float)
        u_13_in = numpy.zeros(label.shape, dtype=float)
        u_23_in = numpy.zeros(label.shape, dtype=float)

        u_11_in[np_index], u_22_in[np_index], u_33_in[np_index] = u_11, u_22, u_33
        u_12_in[np_index], u_13_in[np_index], u_23_in[np_index] = u_12, u_13, u_23

        np_zeros = numpy.zeros(label.shape, dtype=float)
        u_11_in = numpy.where(adp_type=="uani", u_11_in, np_zeros)
        u_22_in = numpy.where(adp_type=="uani", u_22_in, np_zeros)
        u_33_in = numpy.where(adp_type=="uani", u_33_in, np_zeros)
        u_12_in = numpy.where(adp_type=="uani", u_12_in, np_zeros)
        u_13_in = numpy.where(adp_type=="uani", u_13_in, np_zeros)
        u_23_in = numpy.where(adp_type=="uani", u_23_in, np_zeros)

        adp = ADP(u_11=u_11_in, u_22=u_22_in, u_33=u_33_in, 
                  u_12=u_12_in, u_13=u_13_in, u_23=u_23_in, 
                  b_iso=b_iso_in)
        return adp

    def apply_space_group_constraint(self, atom_site, space_group):
        """
        according to table 1 in Peterse, Palm, Acta Cryst.(1966), 20, 147
        """
        l_numb = atom_site.calc_constr_number(space_group)
        label_aniso = self.label
        label = atom_site.label
        l_ind = [label.index(_1) for _1 in label_aniso]
        for index, u_11, u_22, u_33, u_12, u_13, u_23 in zip(l_ind, self.u_11, self.u_22, self.u_33, 
                                                                    self.u_12, self.u_13, self.u_23):
            numb = l_numb[index]
            if numb == 1:
                u_12.value = 0.
                u_12.refinement = False
                u_23.value = 0.
                u_23.refinement = False
            elif numb == 2:
                u_23.value = 0.
                u_23.refinement = False
                u_13.value = 0.
                u_13.refinement = False
            elif numb == 3:
                u_12.value = 0.
                u_12.refinement = False
                u_13.value = 0.
                u_13.refinement = False
            elif numb == 4:
                u_12.value = 0.
                u_12.refinement = False
                u_13.value = 0.
                u_13.refinement = False
                u_23.value = 0.
                u_23.refinement = False
            elif numb == 5:
                u_22.value = u_11.value
                u_22.refinement = False
                u_13.value = 0.
                u_13.refinement = False
                u_23.value = 0.
                u_23.refinement = False
            elif numb == 6:
                u_22.value = u_11.value
                u_22.refinement = False
                u_23.value = u_13.value 
                u_23.refinement = False
            elif numb == 7:
                u_22.value = u_11.value
                u_22.refinement = False
                u_23.value = -1.*u_13.value 
                u_23.refinement = False
            elif numb == 8:
                u_22.value = u_11.value
                u_22.refinement = False
                u_12.value = 0.
                u_12.refinement = False
                u_13.value = 0.
                u_13.refinement = False
                u_23.value = 0.
                u_23.refinement = False
            elif numb == 9:
                u_33.value = u_22.value
                u_33.refinement = False
                u_12.value = 0.
                u_12.refinement = False
                u_13.value = 0.
                u_13.refinement = False
            elif numb == 10:
                u_33.value = u_22.value
                u_33.refinement = False
                u_13.value = 1.*u_12.value 
                u_13.refinement = False
            elif numb == 11:
                u_33.value = u_22.value
                u_33.refinement = False
                u_13.value = -1.*u_12.value 
                u_13.refinement = False
            elif numb == 12:
                u_33.value = u_22.value
                u_33.refinement = False
                u_12.value = 0.
                u_12.refinement = False
                u_13.value = 0.
                u_13.refinement = False
                u_23.value = 0.
                u_23.refinement = False
            elif numb == 13:
                u_12.value = 0.5*u_22.value
                u_12.refinement = False
                u_23.value = 0.
                u_23.refinement = False
            elif numb == 14:
                u_12.value = 0.5*u_22.value
                u_12.refinement = False
                u_13.value = 0.
                u_13.refinement = False
                u_23.value = 0.
                u_23.refinement = False
            elif numb == 15:
                u_12.value = 0.5*u_22.value
                u_12.refinement = False
                u_23.value = 2.*u_13.value
                u_23.refinement = False
            elif numb == 16:
                u_22.value = 1.0*u_11.value
                u_22.refinement = False
                u_12.value = 0.5*u_11.value
                u_12.refinement = False
                u_13.value = 0.
                u_13.refinement = False
                u_23.value = 0.
                u_23.refinement = False
            elif numb == 17:
                u_22.value = 1.0*u_11.value
                u_22.refinement = False
                u_33.value = 1.0*u_11.value
                u_33.refinement = False
                u_12.value = 0.
                u_12.refinement = False
                u_13.value = 0.
                u_13.refinement = False
                u_23.value = 0.
                u_23.refinement = False
            elif numb == 18:
                u_22.value = 1.0*u_11.value
                u_22.refinement = False
                u_33.value = 1.0*u_11.value
                u_33.refinement = False
                u_13.value = 1.0*u_12.value
                u_13.refinement = False
                u_23.value = 1.0*u_12.value
                u_23.refinement = False


    def calc_beta(self, cell):
        """
        calculate beta_ij from U_ij
        """
        ia, ib, ic = cell.ia, cell.ib, cell.ic
        u_11, u_22, u_33 = numpy.array(self.u_11, float), numpy.array(self.u_22, float), numpy.array(self.u_33, float)
        u_12, u_13, u_23 = numpy.array(self.u_12, float), numpy.array(self.u_13, float), numpy.array(self.u_23, float)
        beta_11 = 2.*numpy.pi**2*u_11*ia**2
        beta_22 = 2.*numpy.pi**2*u_22*ib**2
        beta_33 = 2.*numpy.pi**2*u_33*ic**2
        beta_12 = 2.*numpy.pi**2*u_12*ia*ib
        beta_13 = 2.*numpy.pi**2*u_13*ia*ic
        beta_23 = 2.*numpy.pi**2*u_23*ib*ic
        return beta_11, beta_22, beta_33, beta_12, beta_13, beta_23