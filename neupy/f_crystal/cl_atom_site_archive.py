"""
define classes to describe atom_site 
"""
__author__ = 'ikibalin'
__version__ = "2019_08_29"
import os
import numpy


from neupy.f_common.cl_fitable import Fitable
from neupy.f_crystal.cl_fract import Fract
from neupy.f_crystal.cl_adp import ADP
from neupy.f_crystal.cl_magnetism import Magnetism


class AtomSite(object):
    """
    Data items in the ATOM_SITE category record details about
    the atom sites in a crystal structure, such as the positional
    coordinates, atomic displacement parameters, and magnetic moments
    and directions.

    Description in cif file:

    Example 1 - based on data set TOZ of Willis, Beckwith & Tozer [Acta Cryst. (1991), C47, 2276-2277]. 
    loop_
    _atom_site_label
    _atom_site_fract_x
    _atom_site_fract_y
    _atom_site_fract_z
    _atom_site_U_iso_or_equiv
    _atom_site_adp_type
    _atom_site_calc_flag
    _atom_site_calc_attached_atom
      O1   .4154(4)  .5699(1)  .3026(0)  .060(1)  Uani  ?    ?
      C2   .5630(5)  .5087(2)  .3246(1)  .060(2)  Uani  ?    ?
      C3   .5350(5)  .4920(2)  .3997(1)  .048(1)  Uani  ?    ?
      N4   .3570(3)  .5558(1)  .4167(0)  .039(1)  Uani  ?    ?
      C5   .3000(5)  .6122(2)  .3581(1)  .045(1)  Uani  ?    ?
      O21  .6958(5)  .4738(2)  .2874(1)  .090(2)  Uani  ?    ?
      C31  .4869(6)  .3929(2)  .4143(2)  .059(2)  Uani  ?    ?
     # - - - - data truncated for brevity - - - -
      H321C  .04(1)  .318(3)   .320(2)   .14000   Uiso  ?    ?
      H322A  .25(1)  .272(4)   .475(3)   .19000   Uiso  ?    ?
      H322B  .34976  .22118    .40954    .19000   Uiso  calc C322
      H322C  .08(1)  .234(4)   .397(3)   .19000   Uiso  ?    ?

    Example 2 - based on data set TOZ of Willis, Beckwith & Tozer [Acta Cryst. (1991), C47, 2276-2277]. 
 
    loop_
    _atom_site_aniso_label
    _atom_site_aniso_U_11  #_atom_site_aniso_B_11 is not introduced
    _atom_site_aniso_U_22
    _atom_site_aniso_U_33
    _atom_site_aniso_U_12
    _atom_site_aniso_U_13
    _atom_site_aniso_U_23
    _atom_site_aniso_type_symbol
     O1   .071(1) .076(1) .0342(9) .008(1)   .0051(9) -.0030(9) O
     C2   .060(2) .072(2) .047(1)  .002(2)   .013(1)  -.009(1)  C
     C3   .038(1) .060(2) .044(1)  .007(1)   .001(1)  -.005(1)  C
     N4   .037(1) .048(1) .0325(9) .0025(9)  .0011(9) -.0011(9) N
     C5   .043(1) .060(1) .032(1)  .001(1)  -.001(1)   .001(1)  C
     # - - - - data truncated for brevity - - - -
     O21  .094(2) .109(2) .068(1)  .023(2)   .038(1)  -.010(1)  O
     C51  .048(2) .059(2) .049(1)  .002(1)  -.000(1)   .007(1)  C
     C511 .048(2) .071(2) .097(3) -.008(2)  -.003(2)   .010(2)  C
     C512 .078(2) .083(2) .075(2)  .009(2)  -.005(2)   .033(2)  C
     C513 .074(2) .055(2) .075(2)  .004(2)   .001(2)  -.010(2)  C
     # - - - - data truncated for brevity - - - -


    reference: https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Catom_site.html
    """
    def __init__(self, label=numpy.array([], dtype=str), 
                       type_symbol=numpy.array([], dtype=str),
                       frac_x=numpy.array([], dtype=Fitable),
                       frac_y=numpy.array([], dtype=Fitable),
                       frac_z=numpy.array([], dtype=Fitable),
                       b_iso=numpy.array([], dtype=Fitable),
                       adp_type=numpy.array([], dtype=str),
                       occupancy=numpy.array([], dtype=Fitable),

                       aniso_label=numpy.array([], dtype=str),
                       u_11=numpy.array([], dtype=Fitable),
                       u_22=numpy.array([], dtype=Fitable),
                       u_33=numpy.array([], dtype=Fitable),
                       u_12=numpy.array([], dtype=Fitable),
                       u_13=numpy.array([], dtype=Fitable),
                       u_23=numpy.array([], dtype=Fitable),

                       magnetism_aniso_label=numpy.array([], dtype=str),
                       factor_lande=numpy.array([], dtype=Fitable),
                       kappa=numpy.array([], dtype=Fitable),
                       magnetism_type=numpy.array([], dtype=str),
                       magnetism_type_symbol=numpy.array([], dtype=str),
                       chi_11=numpy.array([], dtype=Fitable),
                       chi_12=numpy.array([], dtype=Fitable),
                       chi_13=numpy.array([], dtype=Fitable),
                       chi_22=numpy.array([], dtype=Fitable),
                       chi_23=numpy.array([], dtype=Fitable),
                       chi_33=numpy.array([], dtype=Fitable),

                       l_atom_type = []
                       ):
        super(AtomSite, self).__init__()
        #information from first loop: fract
        self.__atom_site_label = None
        self.__atom_site_type_symbol = None
        self.__atom_site_fract_x = None
        self.__atom_site_fract_y = None
        self.__atom_site_fract_z = None
        self.__atom_site_b_iso_or_equiv = None
        self.__atom_site_adp_type = None
        self.__atom_site_occupancy = None

        #information from second loop: adp
        self.__atom_site_aniso_label = None
        self.__atom_site_aniso_u_11 = None
        self.__atom_site_aniso_u_22 = None
        self.__atom_site_aniso_u_33 = None
        self.__atom_site_aniso_u_12 = None
        self.__atom_site_aniso_u_13 = None
        self.__atom_site_aniso_u_23 = None

        #information from third loop: magnetism
        #todo: it should be separeted 
        self.__atom_site_magnetism_aniso_label = None
        self.__atom_site_magnetism_factor_lande = None
        self.__atom_site_magnetism_kappa = None
        self.__atom_site_magnetism_type = None
        self.__atom_site_magnetism_type_symbol = None
        self.__atom_site_aniso_magnetism_chi_11 = None
        self.__atom_site_aniso_magnetism_chi_12 = None
        self.__atom_site_aniso_magnetism_chi_13 = None
        self.__atom_site_aniso_magnetism_chi_22 = None
        self.__atom_site_aniso_magnetism_chi_23 = None
        self.__atom_site_aniso_magnetism_chi_33 = None

        self.__atom_type_list = []

        #internal flags
        self.__flag_refresh = False

        self.label = label
        self.type_symbol = type_symbol
        self.x = frac_x
        self.y = frac_y
        self.z = frac_z
        self.b_iso = b_iso
        self.adp_type = adp_type
        self.occupancy = occupancy

        self.aniso_label = aniso_label
        self.u_11 = u_11
        self.u_22 = u_22
        self.u_33 = u_33
        self.u_12 = u_12
        self.u_13 = u_13
        self.u_23 = u_23

        self.magnetism_aniso_label = magnetism_aniso_label
        self.factor_lande = factor_lande
        self.kappa = kappa
        self.magnetism_type = magnetism_type
        self.magnetism_type_symbol = magnetism_type_symbol
        self.chi_11 = chi_11
        self.chi_12 = chi_12
        self.chi_13 = chi_13
        self.chi_22 = chi_22
        self.chi_23 = chi_23
        self.chi_33 = chi_33

        self.atom_type = l_atom_type #atom type 
        

        #internal classes
        self.__fract = None
        self.__adp = None
        self.__magnetism = None
        #internal parameters
        self.__atom_site_scat_length_neutron = None #taken from _atom_type_scat_length_neutron
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


    def _trans_to_fitable_array(self, x, name=None):
        if isinstance(x, numpy.ndarray):
            if x.size > 0:
                if isinstance(x[0], Fitable):
                    x_out = x
                else:
                    x_out = numpy.array([Fitable(hh, None, False, name) for hh in x], dtype=Fitable)
            else:
                x_out = numpy.array([], dtype=Fitable)
        else:
            try:
                x_out = numpy.array([Fitable(x, None, False, name)], dtype=Fitable)
            except:
                x_out = None
        return x_out
    def _trans_to_str_array(self, x):
        if isinstance(x, numpy.ndarray):
            x_out = x.astype(dtype=str)
        else:
            try:
                x_out = numpy.array([str(x)], dtype=str)
            except:
                x_out = None
        return x_out

    @property
    def label(self):
        """

        reference:
        """
        return self.__atom_site_label
    @label.setter
    def label(self, x):
        self.__flag_refresh = True
        self.__atom_site_label = self._trans_to_str_array(x)
    @property
    def type_symbol(self):
        """
        A code to identify the atom species (singular or plural)
        occupying this site.
        This code must match a corresponding _atom_type_symbol. The
        specification of this code is optional if component 0 of the
        _atom_site_label is used for this purpose. See _atom_type_symbol.

        Appears in list containing _atom_site_label
        Must match data name_atom_type_symbol

        May match subsidiary data name(s) 
        _atom_site_aniso_type_symbol

        Type: char
        
        reference https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Iatom_site_type_symbol.html        """
        return self.__atom_site_type_symbol
    @type_symbol.setter
    def type_symbol(self, x):
        self.__flag_refresh = True
        self.__atom_site_type_symbol = self._trans_to_str_array(x)
    @property
    def x(self):
        """

        reference:
        """
        return self.__atom_site_fract_x
    @x.setter
    def x(self, x):
        self.__flag_refresh = True
        self.__atom_site_fract_x = self._trans_to_fitable_array(x, name="_atom_site_fract_x")
    @property
    def y(self):
        """

        reference:
        """
        return self.__atom_site_fract_y
    @x.setter
    def y(self, x):
        self.__flag_refresh = True
        self.__atom_site_fract_y = self._trans_to_fitable_array(x, name="_atom_site_fract_y")
    @property
    def z(self):
        """

        reference:
        """
        return self.__atom_site_fract_z
    @z.setter
    def z(self, x):
        self.__flag_refresh = True
        self.__atom_site_fract_z = self._trans_to_fitable_array(x, name="_atom_site_fract_z")
    @property
    def b_iso(self):
        """

        reference:
        """
        return self.__atom_site_b_iso_or_equiv
    @b_iso.setter
    def b_iso(self, x):
        self.__flag_refresh = True
        self.__atom_site_b_iso_or_equiv = self._trans_to_fitable_array(x, name="_atom_site_u_iso_or_equiv")
    @property
    def adp_type(self):
        """

        reference:
        """
        return self.__atom_site_adp_type
    @adp_type.setter
    def adp_type(self, x):
        self.__flag_refresh = True
        self.__atom_site_adp_type = self._trans_to_str_array(x)
    @property
    def occupancy(self):
        """

        reference:
        """
        return self.__atom_site_occupancy
    @occupancy.setter
    def occupancy(self, x):
        self.__flag_refresh = True
        self.__atom_site_occupancy = self._trans_to_fitable_array(x, name="_atom_site_occupancy")



    @property
    def aniso_label(self):
        """

        reference:
        """
        return self.__atom_site_aniso_label
    @aniso_label.setter
    def aniso_label(self, x):
        self.__flag_refresh = True
        self.__atom_site_aniso_label = self._trans_to_str_array(x)
    @property
    def u_11(self):
        """

        reference:
        """
        return self.__atom_site_aniso_u_11
    @u_11.setter
    def u_11(self, x):
        self.__flag_refresh = True
        self.__atom_site_aniso_u_11 = self._trans_to_fitable_array(x, name="_atom_site_aniso_u_11")
    @property
    def u_22(self):
        """

        reference:
        """
        return self.__atom_site_aniso_u_22
    @u_22.setter
    def u_22(self, x):
        self.__flag_refresh = True
        self.__atom_site_aniso_u_22 = self._trans_to_fitable_array(x, name="_atom_site_aniso_u_22")
    @property
    def u_33(self):
        """

        reference:
        """
        return self.__atom_site_aniso_u_33
    @u_33.setter
    def u_33(self, x):
        self.__flag_refresh = True
        self.__atom_site_aniso_u_33 = self._trans_to_fitable_array(x, name="_atom_site_aniso_u_33")
    @property
    def u_12(self):
        """

        reference:
        """
        return self.__atom_site_aniso_u_12
    @u_12.setter
    def u_12(self, x):
        self.__flag_refresh = True
        self.__atom_site_aniso_u_12 = self._trans_to_fitable_array(x, name="_atom_site_aniso_u_12")
    @property
    def u_13(self):
        """

        reference:
        """
        return self.__atom_site_aniso_u_13
    @u_13.setter
    def u_13(self, x):
        self.__flag_refresh = True
        self.__atom_site_aniso_u_13 = self._trans_to_fitable_array(x, name="_atom_site_aniso_u_13")
    @property
    def u_23(self):
        """

        reference:
        """
        return self.__atom_site_aniso_u_23
    @u_23.setter
    def u_23(self, x):
        self.__flag_refresh = True
        self.__atom_site_aniso_u_23 = self._trans_to_fitable_array(x, name="_atom_site_aniso_u_23")



    @property
    def magnetism_aniso_label(self):
        """

        reference:
        """
        return self.__atom_site_magnetism_aniso_label
    @magnetism_aniso_label.setter
    def magnetism_aniso_label(self, x):
        self.__flag_refresh = True
        self.__atom_site_magnetism_aniso_label = self._trans_to_str_array(x)
    @property
    def factor_lande(self):
        """

        reference:
        """
        return self.__atom_site_magnetism_factor_lande
    @factor_lande.setter
    def factor_lande(self, x):
        self.__flag_refresh = True
        self.__atom_site_magnetism_factor_lande = self._trans_to_fitable_array(x, name="_atom_site_magnetism_factor_lande")
    @property
    def kappa(self):
        """

        reference:
        """
        return self.__atom_site_magnetism_kappa
    @kappa.setter
    def kappa(self, x):
        self.__flag_refresh = True
        self.__atom_site_magnetism_kappa = self._trans_to_fitable_array(x, name="_atom_site_magnetism_kappa")
    @property
    def magnetism_type(self):
        """

        reference:
        """
        return self.__atom_site_magnetism_type
    @magnetism_type.setter
    def magnetism_type(self, x):
        self.__flag_refresh = True
        self.__atom_site_magnetism_type = self._trans_to_str_array(x)
    @property
    def magnetism_type_symbol(self):
        """

        reference:
        """
        return self.__atom_site_magnetism_type_symbol
    @magnetism_type_symbol.setter
    def magnetism_type_symbol(self, x):
        self.__flag_refresh = True
        self.__atom_site_magnetism_type_symbol = self._trans_to_str_array(x)
    @property
    def chi_11(self):
        """

        reference:
        """
        return self.__atom_site_aniso_magnetism_chi_11
    @chi_11.setter
    def chi_11(self, x):
        self.__flag_refresh = True
        self.__atom_site_aniso_magnetism_chi_11 = self._trans_to_fitable_array(x, name="_atom_site_aniso_magnetism_chi_11")
    @property
    def chi_22(self):
        """

        reference:
        """
        return self.__atom_site_aniso_magnetism_chi_22
    @chi_22.setter
    def chi_22(self, x):
        self.__flag_refresh = True
        self.__atom_site_aniso_magnetism_chi_22 = self._trans_to_fitable_array(x, name="_atom_site_aniso_magnetism_chi_22")
    @property
    def chi_33(self):
        """

        reference:
        """
        return self.__atom_site_aniso_magnetism_chi_33
    @chi_33.setter
    def chi_33(self, x):
        self.__flag_refresh = True
        self.__atom_site_aniso_magnetism_chi_33 = self._trans_to_fitable_array(x, name="_atom_site_aniso_magnetism_chi_33")
    @property
    def chi_12(self):
        """

        reference:
        """
        return self.__atom_site_aniso_magnetism_chi_12
    @chi_12.setter
    def chi_12(self, x):
        self.__flag_refresh = True
        self.__atom_site_aniso_magnetism_chi_12 = self._trans_to_fitable_array(x, name="_atom_site_aniso_magnetism_chi_12")
    @property
    def chi_13(self):
        """

        reference:
        """
        return self.__atom_site_aniso_magnetism_chi_13
    @chi_13.setter
    def chi_13(self, x):
        self.__flag_refresh = True
        self.__atom_site_aniso_magnetism_chi_13 = self._trans_to_fitable_array(x, name="_atom_site_aniso_magnetism_chi_13")
    @property
    def chi_23(self):
        """

        reference:
        """
        return self.__atom_site_aniso_magnetism_chi_23
    @chi_23.setter
    def chi_23(self, x):
        self.__flag_refresh = True
        self.__atom_site_aniso_magnetism_chi_23 = self._trans_to_fitable_array(x, name="_atom_site_aniso_magnetism_chi_23")
    @property
    def atom_type(self):
        """

        reference:
        """
        return self.__atom_type_list
    @atom_type.setter
    def atom_type(self, x):
        self.__flag_refresh = True
        self.__atom_type_list = x



    @property
    def is_defined(self):
        """
        Output: True if all started parameters are given and size are the same
        """

        cond_1 = any([self.x is None, self.y is None, self.z is None])
        return not(cond_1)

    @property
    def scat_length_neutron(self):
        """

        reference:
        """
        return self.__atom_site_scat_length_neutron
    @property
    def j0_A(self):
        """

        reference:
        """
        return self.__j0_A
    @property
    def j0_a(self):
        """

        reference:
        """
        return self.__j0_a
    @property
    def j0_B(self):
        """

        reference:
        """
        return self.__j0_B
    @property
    def j0_b(self):
        """

        reference:
        """
        return self.__j0_b
    @property
    def j0_C(self):
        """

        reference:
        """
        return self.__j0_C
    @property
    def j0_c(self):
        """

        reference:
        """
        return self.__j0_c
    @property
    def j0_D(self):
        """

        reference:
        """
        return self.__j0_D
    @property
    def j2_A(self):
        """

        reference:
        """
        return self.__j2_A
    @property
    def j2_a(self):
        """

        reference:
        """
        return self.__j2_a
    @property
    def j2_B(self):
        """

        reference:
        """
        return self.__j2_B
    @property
    def j2_b(self):
        """

        reference:
        """
        return self.__j2_b
    @property
    def j2_C(self):
        """

        reference:
        """
        return self.__j2_C
    @property
    def j2_c(self):
        """

        reference:
        """
        return self.__j2_c
    @property
    def j2_D(self):
        """

        reference:
        """
        return self.__j2_D



    def __repr__(self):
        ls_out = ["""AtomSite:\n number of atoms: """]
        return "\n".join(ls_out)
    
    def add_atom(self, atom):
        pass

    def del_atom(self, ind):
        pass
                
    def _show_message(self, s_out: str):
        print("***  Error ***")
        print(s_out)

    def _apply_chi_iso(self, cell):
        pass
        return

    def _take_from_atom_type(self):
        l_atom_type = self.atom_type
        l_type_symbol = [hh.type_symbol for hh in l_atom_type]
        
        if any([hh not in l_type_symbol for hh in self.type_symbol]):
            l_hh = [hh for hh in self.type_symbol if hh not in l_type_symbol]
            self._show_message("The following _atom_type_symbols are not defined\n {:}".format(
                " ".join(l_hh)))
            return False

        l_ind = [l_type_symbol.index(hh) for hh in self.type_symbol]
        l_scat_length_neutron = [l_atom_type[ind].scat_length_neutron for ind in l_ind]
        self.__atom_site_scat_length_neutron = numpy.array(l_scat_length_neutron, dtype=float) #parameters are not refined

        if any([hh not in l_type_symbol for hh in self.magnetism_type_symbol]):
            l_hh = [hh for hh in self.type_symbol if hh not in l_type_symbol]
            self._show_message("The following _atom_type_symbols are not defined\n {:}".format(
                " ".join(l_hh)))
            return False

        l_ind = [l_type_symbol.index(hh) for hh in self.magnetism_type_symbol]
        self.__j0_A = numpy.array([l_atom_type[ind].j0_A for ind in l_ind], dtype=float)
        self.__j0_a = numpy.array([l_atom_type[ind].j0_a for ind in l_ind], dtype=float)
        self.__j0_B = numpy.array([l_atom_type[ind].j0_B for ind in l_ind], dtype=float)
        self.__j0_b = numpy.array([l_atom_type[ind].j0_b for ind in l_ind], dtype=float)
        self.__j0_C = numpy.array([l_atom_type[ind].j0_C for ind in l_ind], dtype=float)
        self.__j0_c = numpy.array([l_atom_type[ind].j0_c for ind in l_ind], dtype=float)
        self.__j0_D = numpy.array([l_atom_type[ind].j0_D for ind in l_ind], dtype=float)
        self.__j2_A = numpy.array([l_atom_type[ind].j2_A for ind in l_ind], dtype=float)
        self.__j2_a = numpy.array([l_atom_type[ind].j2_a for ind in l_ind], dtype=float)
        self.__j2_B = numpy.array([l_atom_type[ind].j2_B for ind in l_ind], dtype=float)
        self.__j2_b = numpy.array([l_atom_type[ind].j2_b for ind in l_ind], dtype=float)
        self.__j2_C = numpy.array([l_atom_type[ind].j2_C for ind in l_ind], dtype=float)
        self.__j2_c = numpy.array([l_atom_type[ind].j2_c for ind in l_ind], dtype=float)
        self.__j2_D = numpy.array([l_atom_type[ind].j2_D for ind in l_ind], dtype=float)
        return True

    def _form_fract_array(self):
        if not(self.is_defined):
            return None
        fract = Fract(x=self.x.astype(float), 
                     y=self.y.astype(float),
                     z=self.z.astype(float))
        self.__fract = fract

    def _form_adp_array(self, cell):
        label = self.label.astype(str)
        b_iso = self.b_iso #it is not used but is should be reconverted to u_ij

        #default values
        u_11_in = numpy.zeros(label.shape, dtype=float)
        u_22_in = numpy.zeros(label.shape, dtype=float)
        u_33_in = numpy.zeros(label.shape, dtype=float)
        u_12_in = numpy.zeros(label.shape, dtype=float)
        u_13_in = numpy.zeros(label.shape, dtype=float)
        u_23_in = numpy.zeros(label.shape, dtype=float)

        aniso_label = self.aniso_label.astype(str)
        if not(set(aniso_label).issubset(set(label))):
            self._show_message("Unknown 'aniso_label'")
            return False

        u_11, u_22, u_33 = self.u_11.astype(float), self.u_22.astype(float), self.u_33.astype(float)
        u_12, u_13, u_23 = self.u_12.astype(float), self.u_13.astype(float), self.u_23.astype(float)

        np_index = numpy.array([int(numpy.argwhere(label==hh)[0]) for hh in aniso_label], dtype=int)

        u_11_in[np_index], u_22_in[np_index], u_33_in[np_index] = u_11, u_22, u_33
        u_12_in[np_index], u_13_in[np_index], u_23_in[np_index] = u_12, u_13, u_23

        adp = ADP(u_11=u_11_in, u_22=u_22_in, u_33=u_33_in,
                  u_12=u_12_in, u_13=u_13_in, u_23=u_23_in)
        self.__adp = adp
        return True

    def _form_magnetism_array(self):
        label = self.label.astype(str)

        magnetism_aniso_label = self.magnetism_aniso_label.astype(str)
        if not(set(magnetism_aniso_label).issubset(set(label))):
            self._show_message("Unknown 'magnetism_aniso_label'")
            return False

        #default values
        kappa_in = numpy.ones(label.shape, dtype=float)
        factor_lande_in = 2.*numpy.ones(label.shape, dtype=float) 
        chi_11_in = numpy.zeros(label.shape, dtype=float)
        chi_22_in = numpy.zeros(label.shape, dtype=float)
        chi_33_in = numpy.zeros(label.shape, dtype=float)
        chi_12_in = numpy.zeros(label.shape, dtype=float)
        chi_13_in = numpy.zeros(label.shape, dtype=float)
        chi_23_in = numpy.zeros(label.shape, dtype=float)
        j0_A_in = numpy.zeros(label.shape, dtype=float)
        j0_a_in = numpy.zeros(label.shape, dtype=float)
        j0_B_in = numpy.zeros(label.shape, dtype=float)
        j0_b_in = numpy.zeros(label.shape, dtype=float)
        j0_C_in = numpy.zeros(label.shape, dtype=float)
        j0_c_in = numpy.zeros(label.shape, dtype=float)
        j0_D_in = numpy.zeros(label.shape, dtype=float)
        j2_A_in = numpy.zeros(label.shape, dtype=float)
        j2_a_in = numpy.zeros(label.shape, dtype=float)
        j2_B_in = numpy.zeros(label.shape, dtype=float)
        j2_b_in = numpy.zeros(label.shape, dtype=float)
        j2_C_in = numpy.zeros(label.shape, dtype=float)
        j2_c_in = numpy.zeros(label.shape, dtype=float)
        j2_D_in = numpy.zeros(label.shape, dtype=float)

        kappa, factor_lande = self.kappa, self.factor_lande

        chi_11, chi_22, chi_33 = self.chi_11.astype(float), self.chi_22.astype(float), self.chi_33.astype(float)
        chi_12, chi_13, chi_23 = self.chi_12.astype(float), self.chi_13.astype(float), self.chi_23.astype(float)

        j0_A, j0_a, j0_B, j0_b = self.j0_A, self.j0_a, self.j0_B, self.j0_b
        j0_C, j0_c, j0_D = self.j0_C, self.j0_c, self.j0_D

        j2_A, j2_a, j2_B, j2_b = self.j2_A, self.j2_a, self.j2_B, self.j2_b
        j2_C, j2_c, j2_D = self.j2_C, self.j2_c, self.j2_D


        np_index = numpy.array([int(numpy.argwhere(label==hh)[0]) for hh in magnetism_aniso_label], dtype=int)

        kappa_in[np_index], factor_lande_in[np_index] = kappa, factor_lande

        chi_11_in[np_index], chi_22_in[np_index], chi_33_in[np_index] = chi_11, chi_22, chi_33
        chi_12_in[np_index], chi_13_in[np_index], chi_23_in[np_index] = chi_12, chi_13, chi_23

        j0_A_in[np_index], j0_a_in[np_index], j0_B_in[np_index], j0_b_in[np_index] = j0_A, j0_a, j0_B, j0_b
        j0_C_in[np_index], j0_c_in[np_index], j0_D_in[np_index] = j0_C, j0_c, j0_D

        j2_A_in[np_index], j2_a_in[np_index], j2_B_in[np_index], j2_b_in[np_index] = j2_A, j2_a, j2_B, j2_b
        j2_C_in[np_index], j2_c_in[np_index], j2_D_in[np_index] = j2_C, j2_c, j2_D

        magnetism = Magnetism(kappa=kappa_in, factor_lande=factor_lande_in, 
                              chi_11=chi_11_in, chi_22=chi_22_in, chi_33=chi_33_in,
                              chi_12=chi_12_in, chi_13=chi_13_in, chi_23=chi_23_in,
                              j0_A=j0_A_in, j0_a=j0_a_in, j0_B=j0_B_in, j0_b=j0_b_in,
                              j0_C=j0_C_in, j0_c=j0_c_in, j0_D=j0_D_in,
                              j2_A=j2_A_in, j2_a=j2_a_in, j2_B=j2_B_in, j2_b=j2_b_in,
                              j2_C=j2_C_in, j2_c=j2_c_in, j2_D=j2_D_in)

        self.__magnetism = magnetism
        return True

    def _form_arrays(self, cell):
        if not(self.is_defined):
            return None

        self._take_from_atom_type()
        self._apply_chi_iso(cell)

        self._form_fract_array()
        self._form_adp_array(cell)
        self._form_magnetism_array()


    
    def calc_sf(self, space_group, cell, h, k, l, d_in={}):
        """
        calculate nuclear structure factor
        """
        if (self.is_variable() | self.__flag_refresh):
            self._form_arrays(cell)

        occupancy = self.occupancy
        scat_length_neutron = self.__atom_site_scat_length_neutron

        fract = self.__fract
        adp = self.__adp
        magnetism = self.__magnetism


        x, y, z = fract.x, fract.y, fract.z

        atom_multiplicity = space_group.calc_atom_mult(x, y, z)
        occ_mult = occupancy*atom_multiplicity 
        

        phase_3d = fract.calc_phase(space_group, h, k, l)#3d object
        dwf_3d = adp.calc_dwf(space_group, cell, h, k, l)
        ff_11, ff_12, ff_13, ff_21, ff_22, ff_23, ff_31, ff_32, ff_33 = \
                   magnetism.calc_form_factor_tensor(space_group, cell, h, k, l)

        hh = phase_3d*dwf_3d

        phase_2d = hh.sum(axis=2)

        ft_11 = (ff_11*hh).sum(axis=2)
        ft_12 = (ff_12*hh).sum(axis=2)
        ft_13 = (ff_13*hh).sum(axis=2)
        ft_21 = (ff_21*hh).sum(axis=2)
        ft_22 = (ff_22*hh).sum(axis=2)
        ft_23 = (ff_23*hh).sum(axis=2)
        ft_31 = (ff_31*hh).sum(axis=2)
        ft_32 = (ff_32*hh).sum(axis=2)
        ft_33 = (ff_33*hh).sum(axis=2)
        
        b_scat_2d = numpy.meshgrid(h, scat_length_neutron, indexing="ij")[1]
        occ_mult_2d = numpy.meshgrid(h, occ_mult, indexing="ij")[1]
        
        l_el_symm = space_group.el_symm
        l_orig = space_group.orig
        centr = space_group.centr

        #calculation of nuclear structure factor        
        hh = phase_2d * b_scat_2d * occ_mult_2d
        f_hkl_as = hh.sum(axis=1)*1./len(l_el_symm)
        
        orig_x = [hh[0] for hh in l_orig]
        orig_y = [hh[1] for hh in l_orig]
        orig_z = [hh[2] for hh in l_orig]
        
        np_h, np_orig_x = numpy.meshgrid(h, orig_x, indexing = "ij")
        np_k, np_orig_y = numpy.meshgrid(k, orig_y, indexing = "ij")
        np_l, np_orig_z = numpy.meshgrid(l, orig_z, indexing = "ij")
        
        np_orig_as = numpy.exp(2*numpy.pi*1j*(np_h*np_orig_x+np_k*np_orig_y+np_l*np_orig_z))
        f_hkl_as = f_hkl_as*np_orig_as.sum(axis=1)*1./len(l_orig)

        if (centr):
            orig = space_group.get_val("p_centr")
            f_nucl = 0.5*(f_hkl_as+f_hkl_as.conjugate()*numpy.exp(2.*2.*numpy.pi*1j* (h*orig[0]+k*orig[1]+l*orig[2])))
        else:
            f_nucl = f_hkl_as

        #calculation of structure factor tensor
        sft_as_11 = (ft_11 * occ_mult_2d).sum(axis=1)*1./len(l_el_symm)
        sft_as_12 = (ft_12 * occ_mult_2d).sum(axis=1)*1./len(l_el_symm)
        sft_as_13 = (ft_13 * occ_mult_2d).sum(axis=1)*1./len(l_el_symm)
        sft_as_21 = (ft_21 * occ_mult_2d).sum(axis=1)*1./len(l_el_symm)
        sft_as_22 = (ft_22 * occ_mult_2d).sum(axis=1)*1./len(l_el_symm)
        sft_as_23 = (ft_23 * occ_mult_2d).sum(axis=1)*1./len(l_el_symm)
        sft_as_31 = (ft_31 * occ_mult_2d).sum(axis=1)*1./len(l_el_symm)
        sft_as_32 = (ft_32 * occ_mult_2d).sum(axis=1)*1./len(l_el_symm)
        sft_as_33 = (ft_33 * occ_mult_2d).sum(axis=1)*1./len(l_el_symm)

        sft_as_11 = sft_as_11 * np_orig_as.sum(axis=1)*1./len(l_orig)
        sft_as_12 = sft_as_12 * np_orig_as.sum(axis=1)*1./len(l_orig)
        sft_as_13 = sft_as_13 * np_orig_as.sum(axis=1)*1./len(l_orig)
        sft_as_21 = sft_as_21 * np_orig_as.sum(axis=1)*1./len(l_orig)
        sft_as_22 = sft_as_22 * np_orig_as.sum(axis=1)*1./len(l_orig)
        sft_as_23 = sft_as_23 * np_orig_as.sum(axis=1)*1./len(l_orig)
        sft_as_31 = sft_as_31 * np_orig_as.sum(axis=1)*1./len(l_orig)
        sft_as_32 = sft_as_32 * np_orig_as.sum(axis=1)*1./len(l_orig)
        sft_as_33 = sft_as_33 * np_orig_as.sum(axis=1)*1./len(l_orig)
    
        if (centr):
            orig = space_group.get_val("p_centr")
            hh = numpy.exp(2.*2.*numpy.pi*1j* (h*orig[0]+k*orig[1]+l*orig[2]))
            sft_11 = 0.5*(sft_as_11+sft_as_11.conjugate()*hh)
            sft_12 = 0.5*(sft_as_12+sft_as_12.conjugate()*hh)
            sft_13 = 0.5*(sft_as_13+sft_as_13.conjugate()*hh)
            sft_21 = 0.5*(sft_as_21+sft_as_21.conjugate()*hh)
            sft_22 = 0.5*(sft_as_22+sft_as_22.conjugate()*hh)
            sft_23 = 0.5*(sft_as_23+sft_as_23.conjugate()*hh)
            sft_31 = 0.5*(sft_as_31+sft_as_31.conjugate()*hh)
            sft_32 = 0.5*(sft_as_32+sft_as_32.conjugate()*hh)
            sft_33 = 0.5*(sft_as_33+sft_as_33.conjugate()*hh)
        else:
            sft_11, sft_12, sft_13 = sft_as_11, sft_as_12, sft_as_13
            sft_21, sft_22, sft_23 = sft_as_21, sft_as_22, sft_as_23
            sft_31, sft_32, sft_33 = sft_as_31, sft_as_32, sft_as_33            
            

        d_map["out"] = f_nucl, sft_11, sft_12, sft_13, sft_21, sft_22, sft_23, sft_31, sft_32, sft_33
        return f_nucl, sft_11, sft_12, sft_13, sft_21, sft_22, sft_23, sft_31, sft_32, sft_33

    
    def is_variable_phase(self):
        res = any([atom_type.is_variable_phase() for atom_type in 
                   self.__atom_type_list])
        return res

    def is_variable_adp(self):
        res = any([atom_type.is_variable_adp() for atom_type in 
                   self.__atom_type_list])
        return res

    def is_variable_magnetism(self):
        res = any([atom_type.is_variable_magnetism() for atom_type in 
                   self.__atom_type_list])
        return res
    
    def is_variable(self):
        res = any([atom_type.is_variable() for atom_type in 
                   self.__atom_type_list])
        return res
    
    def get_variables(self):
        l_variable = []
        for atom_type in self.__atom_type_list:
            l_var = atom_type.get_variables()    
            l_variable.extend(l_var)
        return l_variable
    
    def apply_constraint(self, space_group, cell):
        for atom_type in self.__atom_type_list:
            atom_type.apply_constraint(space_group, cell)
        self._form_arrays(cell)
