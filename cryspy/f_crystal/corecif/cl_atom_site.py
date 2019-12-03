"""
define classes to describe AtomSite
"""
__author__ = 'ikibalin'
__version__ = "2019_12_30"
import os
import numpy
from pycifstar import Global

import warnings
from typing import List, Tuple
from cryspy.f_common.cl_item_constr import ItemConstr
from cryspy.f_common.cl_loop_constr import LoopConstr

from cryspy.f_common.cl_fitable import Fitable
from cryspy.f_crystal.cl_fract import Fract

class AtomSite(ItemConstr):
    """
AtomSite:
=========
Data items in the ATOM_SITE category record details about
the atom sites in a crystal structure, such as the positional
coordinates.

Description in cif file:
-------------------------
_atom_site_label            Fe3A
_atom_site_type_symbol      Fe 
_atom_site_fract_x          0.12500
_atom_site_fract_y          0.12500
_atom_site_fract_z          0.12500
_atom_site_adp_type         Uani
_atom_site_B_iso_or_equiv   0.0 
_atom_site_occupancy        1.0

Mandatory attribute: 
---------------------
- label
- type_symbol
- fract_x
- fract_y
- fract_z

Optional attribute: 
---------------------
- occupancy
- adp_type
- wyckoff_symbol
- u_iso_or_equiv
- u_equiv_geom_mean
- b_iso_or_equiv
- cartn_x
- cartn_y
- cartn_z

Class methods:
----------------
- 

Methods:
---------
- get_symop

Reference: 
-----------
`the page <https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Catom_site.html>`_
    """    
    MANDATORY_ATTRIBUTE = ("label", "type_symbol", "fract_x", "fract_y", "fract_z")
    OPTIONAL_ATTRIBUTE = ("occupancy", "adp_type", "wyckoff_symbol", 
                          "u_iso_or_equiv", "u_equiv_geom_mean", "b_iso_or_equiv",
                          "cartn_x", "cartn_y", "cartn_z")
    INTERNAL_ATTRIBUTE = ("scat_length_neutron", "atom_type_scat", "space_group_wyckoff")
    ACCESIBLE_ADP_TYPE = ("Uani", "Uiso", "Uovl", "Umpe", "Bani", "Biso", "Bovl")
    PREFIX = "atom_site"

    def __init__(self, label=None, type_symbol=None, fract_x=None, fract_y=None, fract_z=None, 
    occupancy=None, adp_type=None, wyckoff_symbol=None, 
    u_iso_or_equiv=None, u_equiv_geom_mean=None, b_iso_or_equiv=None, 
    cartn_x=None, cartn_y=None, cartn_z=None):
        super(AtomSite, self).__init__(mandatory_attribute=self.MANDATORY_ATTRIBUTE, 
                                       optional_attribute=self.OPTIONAL_ATTRIBUTE, 
                                       internal_attribute=self.INTERNAL_ATTRIBUTE,
                                       prefix=self.PREFIX)

        self.label = label
        self.type_symbol = type_symbol
        self.fract_x = fract_x
        self.fract_y = fract_y
        self.fract_z = fract_z
        self.occupancy = occupancy
        self.adp_type = adp_type
        self.wyckoff_symbol = wyckoff_symbol
        self.u_iso_or_equiv = u_iso_or_equiv
        self.u_equiv_geom_mean = u_equiv_geom_mean
        self.b_iso_or_equiv = b_iso_or_equiv
        self.cartn_x = cartn_x
        self.cartn_y = cartn_y
        self.cartn_z = cartn_z

        if self.is_defined:
            self.form_object

    @property
    def label(self)->str:
        """
        The _atom_site_label is a unique identifier for a particular site
        in the crystal. 

        Type: char

        Reference: 
        -----------
        `link <https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Iatom_site_label.html>`_
        """
        return getattr(self, "__label")
    @label.setter
    def label(self, x):
        if x is None:
            x_in = None
        else:
            x_in = str(x)
        setattr(self, "__label", x_in)

    @property
    def type_symbol(self):
        """
        A code to identify the atom species (singular or plural)
        occupying this site.
        This code must match a corresponding _atom_type_symbol. 

        Type: char

        Reference: 
        ------------
        `link <https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Iatom_site_type_symbol.html>`_
        """
        return getattr(self, "__type_symbol")
    @type_symbol.setter
    def type_symbol(self, x):
        if x is None:
            x_in = None
        else:
            x_in = str(x)
        setattr(self, "__type_symbol", x_in)

    @property
    def fract_x(self):
        """
        Atom-site coordinates as fractions of the _cell_length_ values.

        Type: float

        Reference: 
        `link <https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Iatom_site_fract_.html>`_
        """
        return getattr(self, "__fract_x")
    @fract_x.setter
    def fract_x(self, x):
        if x is None:
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__fract_x", x_in)

    @property
    def fract_y(self):
        """
        Atom-site coordinates as fractions of the _cell_length_ values.

        Type: float

        Reference: 
        `link <https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Iatom_site_fract_.html>`_
        """
        return getattr(self, "__fract_y")
    @fract_y.setter
    def fract_y(self, x):
        if x is None:
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__fract_y", x_in)

    @property
    def fract_z(self):
        """
        Atom-site coordinates as fractions of the _cell_length_ values.

        Type: float

        Reference: 
        `link <https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Iatom_site_fract_.html>`_
        """
        return getattr(self, "__fract_z")
    @fract_z.setter
    def fract_z(self, x):
        if x is None:
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__fract_z", x_in)


    @property
    def adp_type(self):
        """
        A standard code used to describe the type of atomic displacement
        parameters used for the site. 

        The data value must be one of the following:
        ----------------------------------------------
        "Uani", "Uiso", "Uovl", "Umpe", "Bani", "Biso", "Bovl"

        Reference: 
        -----------
        `link https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Iatom_site_adp_type.html`_
        """
        return getattr(self, "__adp_type")
    @adp_type.setter
    def adp_type(self, x):
        if x is None:
            x_in = None
        else:
            if not(x in self.ACCESIBLE_ADP_TYPE):
                warnings.warn(f"adp_type '{x_in:}' is not supported", UserWarning)
                x_in = None            
            x_in = str(x)
        setattr(self, "__adp_type", x_in)

    @property
    def occupancy(self):
        """
        The fraction of the atom type present at this site.
        The sum of the occupancies of all the atom types at this site
        may not significantly exceed 1.0 unless it is a dummy site. The
        value must lie in the 99.97% Gaussian confidence interval
        -3u =< x =< 1 + 3u. The _enumeration_range of 0.0:1.0 is thus
        correctly interpreted as meaning (0.0 - 3u) =< x =< (1.0 + 3u).

        
        Reference:
        -------------
        `link https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Iatom_site_occupancy.html`_
        """
        return getattr(self, "__occupancy")
    @occupancy.setter
    def occupancy(self, x):
        if x is None:
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__occupancy", x_in)

    @property
    def wyckoff_symbol(self):
        """
        The Wyckoff symbol (letter) as listed in the space-group tables of
        International Tables for Crystallography Vol. A (2002).
        """
        return getattr(self, "__wyckoff_symbol")
    @wyckoff_symbol.setter
    def wyckoff_symbol(self, x):
        if x is None:
            x_in = None
        else: 
            x_in = str(x)
        setattr(self, "__wyckoff_symbol", x_in)

    @property
    def u_iso_or_equiv(self):
        """
        Isotropic atomic displacement parameter, or equivalent isotropic
        atomic displacement parameter, U equiv , in ångströms squared, 
        calculated from anisotropic atomic displacement parameters.


        U(equiv) = (1/3) sum~i~[sum~j~(U^ij^ a*~i~ a*~j~ a~i~ a~j~)]

        a     = the real-space cell lengths
        a*    = the reciprocal-space cell lengths

        Ref: Fischer, R. X. & Tillmanns, E. (1988). Acta Cryst. C44,
             775-776.
        """
        return getattr(self, "__u_iso_or_equiv")
    @u_iso_or_equiv.setter
    def u_iso_or_equiv(self, x):
        if x is None:
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__u_iso_or_equiv", x_in)


    @property
    def u_equiv_geom_mean(self):
        """
        Equivalent isotropic atomic displacement parameter, U equiv , in
        angströms squared, calculated as the geometric mean of the
        anisotropic atomic displacement parameters.

        U_equiv = (U_i U_j U_k)**(1/3)

        where U_n = the principal components of the orthogonalized U_ij.
        """
        return getattr(self, "__u_equiv_geom_mean")
    @u_equiv_geom_mean.setter
    def u_equiv_geom_mean(self, x):
        if x is None:
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__u_equiv_geom_mean", x_in)


    @property
    def b_iso_or_equiv(self):
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

        Reference:
        ------------
        `link <https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Iatom_site_B_iso_or_equiv.html>`_
        """
        return getattr(self, "__b_iso_or_equiv")
    @b_iso_or_equiv.setter
    def b_iso_or_equiv(self, x):
        if x is None:
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__b_iso_or_equiv", x_in)

    @property
    def cartn_x(self):
        """
        The atom-site coordinates in angstroms specified according to a set
        of orthogonal Cartesian axes related to the cell axes as specified by
        the _atom_sites_Cartn_transform_axes description.
        """
        return getattr(self, "__cartn_x")
    @cartn_x.setter
    def cartn_x(self, x):
        if x is None:
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__cartn_x", x_in)

    @property
    def cartn_y(self):
        """
        The atom-site coordinates in angstroms specified according to a set
        of orthogonal Cartesian axes related to the cell axes as specified by
        the _atom_sites_Cartn_transform_axes description.
        """
        return getattr(self, "__cartn_y")
    @cartn_y.setter
    def cartn_y(self, x):
        if x is None:
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__cartn_y", x_in)

    @property
    def cartn_z(self):
        """
        The atom-site coordinates in angstroms specified according to a set
        of orthogonal Cartesian axes related to the cell axes as specified by
        the _atom_sites_Cartn_transform_axes description.
        """
        return getattr(self, "__cartn_z")
    @cartn_z.setter
    def cartn_z(self, x):
        if x is None:
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__cartn_z", x_in)


    def _show_message(self, s_out: str):
        warnings.warn("***  Error ***", UserWarning)
        warnings.warn(s_out, UserWarning)


    @property
    def is_variable(self):
        res = any([self.frac_x.refinement,
                   self.frac_y.refinement,
                   self.frac_z.refinement,
                   self.u_iso_or_equiv.refinement,
                   self.u_equiv_geom_mean.refinement,
                   self.b_iso_or_equiv.refinement,
                   self.occupancy.refinement])
        return res


    def get_variables(self):
        l_variable = []
        if self.frac_x.refinement: l_variable.append(self.frac_x)
        if self.frac_y.refinement: l_variable.append(self.frac_y)
        if self.frac_z.refinement: l_variable.append(self.frac_z)
        if self.u_iso_or_equiv.refinement: l_variable.append(self.u_iso_or_equiv)
        if self.u_equiv_geom_mean.refinement: l_variable.append(self.u_equiv_geom_mean)
        if self.b_iso_or_equiv.refinement: l_variable.append(self.b_iso_or_equiv)
        if self.occupancy.refinement: l_variable.append(self.occupancy)
        return l_variable


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


    @property
    def scat_length_neutron(self):
        if self.__scat_length_neutron is None:
            l_res = [self._get_scat_length_neutron(type_symbol) for type_symbol in self.type_symbol]
            res = tuple(l_res)
            self.__scat_length_neutron = l_res 
        else:
            res = tuple(self.__scat_length_neutron)
        return res


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



class AtomSiteL(LoopConstr):
    """
    AtomSiteL:
    ============
    Data items in the ATOM_SITE category record details about
    the atom sites in a crystal structure, such as the positional
    coordinates.

    Description in cif file:
    ---------------------------
    loop_                                     
    _atom_site_label          
    _atom_site_type_symbol   
    _atom_site_fract_x       
    _atom_site_fract_y       
    _atom_site_fract_z       
    _atom_site_adp_type       
    _atom_site_B_iso_or_equiv
    _atom_site_occupancy     
     Fe3A   Fe  0.12500 0.12500 0.12500  Uani   0.0   1.0
     Fe3B   Fe  0.50000 0.50000 0.50000  Uani   0.0   1.0
     O1     O   0.25521 0.25521 0.25521  Uiso   0.0   1.0

    Mandatory attribute: 
    ---------------------
    - label
    - type_symbol
    - fract_x
    - fract_y
    - fract_z
    
    Optional attribute: 
    ---------------------
    - occupancy
    - adp_type
    - wyckoff_symbol
    - u_iso_or_equiv
    - u_equiv_geom_mean
    - b_iso_or_equiv
    - cartn_x
    - cartn_y
    - cartn_z
    
    Class methods:
    -----------------
    - 
    
    Methods:
    ---------
    - get_symop
    
    Reference: 
    -----------
    `the page <https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Catom_site.html>`_    
    """
    CATEGORY_KEY = ("label", )
    ITEM_CLASS = AtomSite
    def __init__(self, item=[], loop_name=""):
        super(AtomSiteL, self).__init__(category_key=self.CATEGORY_KEY, item_class=self.ITEM_CLASS, loop_name=loop_name)
        self.item = item
