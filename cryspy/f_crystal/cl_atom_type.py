"""
define classes to describe AtomType
"""
__author__ = 'ikibalin'
__version__ = "2019_08_26"
import os
import numpy


from cryspy.f_common.cl_fitable import Fitable

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
    
    """    
    def __init__(self, symbol=[], scat_length_neutron=[]):
        super(AtomType, self).__init__()

        self.__atom_type_symbol = []
        self.__atom_type_scat_length_neutron = []

    @property
    def symbol(self):
        """
        """
        return tuple(self.__atom_type_symbol)
    @symbol.setter
    def symbol(self, l_x):
        l_x_in = []
        for x in l_x:
            x_in = str(x).strip()
            l_x_in.append(x_in)
        self.__atom_type_symbol = l_x_in
        len_x = len(l_x_in)
        len_1 = len(self.__atom_type_scat_length_neutron)
        if len_1 > len_x:
            self.__atom_type_scat_length_neutron = self.__atom_type_scat_length_neutron[:len_x]
        elif len_1 < len_x:
            l_fitable = [Fitable(value=0., name="_atom_type_scat_length_neutron") for hh in range(len_x-len_1)]
            self.__atom_type_scat_length_neutron.extend(l_fitable)


    @property
    def scat_length_neutron(self):
        """
        in 10**-12 cm 
        """
        return tuple(self.__atom_type_scat_length_neutron)
    @scat_length_neutron.setter
    def scat_length_neutron(self, l_x):
        l_fitable = []
        for x in l_x:
            if isinstance(x, Fitable):
                x_in = x
            else:
                x_in = Fitable()
                flag = x_in.take_it(x)
            l_fitable.append(x_in)
        len_x = len(l_fitable)
        len_1 = len(self.__atom_type_symbol)
        if len_1 < len_x:
            l_fitable = l_fitable[:len_1]
        elif len_1 > len_x:
            l_fitable.extend([Fitable(value=0., name="_atom_type_scat_length_neutron") for hh in range(len_1-len_x)])
        self.__atom_type_scat_length_neutron = l_fitable

    @property
    def b_scat(self):
        return self.scat_length_neutron
    @b_scat.setter
    def b_scat(self, l_x):
        self.scat_length_neutron = l_x
        
    def __repr__(self):
        ls_out = ["AtomType:"]
        ls_out.append(str(self))
        return "\n".join(ls_out)     

    def __str__(self):
        ls_out = []
        ls_out.append(" symbol   scat_length_neutron")
        for _1, _2 in zip(self.symbol, self.scat_length_neutron):
            ls_out.append(" {:8} {:8}".format(_1, _2.print_with_sigma))
        return "\n".join(ls_out)

    def _show_message(self, s_out: str):
        print("***  Error ***")
        print(s_out)

    @property
    def is_defined(self):
        """
        Output: True if all started parameters are given
        """
        cond = (self.symbol != [])
        return cond

    @property
    def is_variable(self):
        res = (any([hh.refinement for hh in self.scat_length_neutron]))
        return res
    def get_variables(self):
        l_variable = []
        l_variable.extend([hh for hh in self.scat_length_neutron if hh.refinement])
        return l_variable

    @property
    def to_cif(self):
        ls_out = []
        if self.is_defined:
            ls_out.append("loop_")
            ls_out.append("_atom_type_symbol")
            ls_out.append("_atom_type_scat_length_neutron")
            for _1, _2 in zip(self.symbol, self.scat_length_neutron):
                ls_out.append("{:} {:}".format(_1, _2.print_with_sigma))
        return "\n".join(ls_out)

    def from_cif(self, string: str):
        cif_global = Global()
        flag = cif_global.take_from_string(string)
        if not flag:
            return False
        flag = False
        flag = cif_global.is_prefix("_atom_type")
        if flag:
            cif_loop = cif_global["_atom_type"]
            l_name = cif_loop.names
            if "_atom_type_symbol" in l_name:
                self.symbol = cif_loop["_atom_type_symbol"]

            if "_atom_type_scat_length_neutron" in l_name:
                l_val = cif_loop["_atom_type_scat_length_neutron"]
                l_fitable = []
                for val in l_val:
                    fitable = Fitable(name="_atom_type_scat_length_neutron")
                    fitable.take_it(val)
                    l_fitable.append(fitable)
                self.scat_length_neutron = l_fitable
        else:
            self.symbol = []
        return True

    def _form_scat_length_neutron(self, atom_site):
        """
        Note: It is not correct!!! It should be changed if it is neaded
        """
        label = numpy.array(atom_site.label, dtype=str)
        label_magnetism = numpy.array(self.label, dtype=str)
        if not(set(label_magnetism).issubset(set(label))):
            self._show_message("Unknown 'aniso_label'")
            return False
        np_index = numpy.array([int(numpy.argwhere(label==hh)[0]) for hh in label_magnetism], dtype=int)
        scat_length_neutron_in = numpy.zeros(label.shape, dtype=float)
        scat_length_neutron = numpy.array(self.scat_length_neutron, dtype=float)
        scat_length_neutron_in[np_index] = scat_length_neutron
        return scat_length_neutron_in