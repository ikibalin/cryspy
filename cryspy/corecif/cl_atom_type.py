"""
define classes AtomType, AtomTypeL
"""

__author__ = 'ikibalin'
__version__ = "2019_12_05"
import os
import numpy
from pycifstar import Global

import warnings
from typing import List, Tuple
from cryspy.common.cl_item_constr import ItemConstr
from cryspy.common.cl_loop_constr import LoopConstr
from cryspy.common.cl_fitable import Fitable



class AtomType(ItemConstr):
    """
AtomType
============
Data items in the ATOM_TYPE category record details about
properties of the atoms that occupy the atom sites, such as the
atomic scattering factors.

Description in cif file:
:: 

 _atom_type_symbol                  C
 _atom_type_oxidation_number        0
 _atom_type_number_in_cell          72
 _atom_type_scat_dispersion_real    0.017
 _atom_type_scat_dispersion_imag    0.009
 _atom_type_scat_source             International_Tables_Vol_IV_Table_2.2B

Reference: 
--------------
`iucr.org <https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Catom_type.html>`_
FIXME:
replace attribute 'analytical_mass' by 'analytical_mass_%' as in CoreCIF
    """
    MANDATORY_ATTRIBUTE = ("symbol")
    OPTIONAL_ATTRIBUTE = ("analytical_mass", "description", "number_in_cell", "oxidation_number", "radius_bond", "radius_contact", 
                          "cromer_mann_a1", "cromer_mann_a2", "cromer_mann_a3", "cromer_mann_a4", 
                          "cromer_mann_b1", "cromer_mann_b2", "cromer_mann_b3", "cromer_mann_b4", "cromer_mann_c", 
                          "scat_dispersion_real", "scat_dispersion_imag", "dispersion_source", "scat_length_neutron",
                          "scat_source", "scat_versus_stol_list")
    INTERNAL_ATTRIBUTE = ()
    PREFIX = "atom_type"

    def __init__(self, symbol=None, 
                 analytical_mass=None, description=None, number_in_cell=None, oxidation_number=None, radius_bond=None, radius_contact=None, 
                 cromer_mann_a1=None, cromer_mann_a2=None, cromer_mann_a3=None, cromer_mann_a4=None,
                 cromer_mann_b1=None, cromer_mann_b2=None, cromer_mann_b3=None, cromer_mann_b4=None, cromer_mann_c=None,
                 scat_dispersion_real=None, scat_dispersion_imag=None, dispersion_source=None, scat_length_neutron=None,
                 scat_source=None, scat_versus_stol_list=None):
        super(AtomType, self).__init__(mandatory_attribute=self.MANDATORY_ATTRIBUTE, 
                                       optional_attribute=self.OPTIONAL_ATTRIBUTE, 
                                       internal_attribute=self.INTERNAL_ATTRIBUTE,
                                       prefix=self.PREFIX)

        self.symbol = symbol
        self.analytical_mass = analytical_mass
        self.description = description
        self.number_in_cell = number_in_cell
        self.oxidation_number = oxidation_number
        self.radius_bond = radius_bond
        self.radius_contact = radius_contact
        self.cromer_mann_a1 = cromer_mann_a1
        self.cromer_mann_a2 = cromer_mann_a2
        self.cromer_mann_a3 = cromer_mann_a3
        self.cromer_mann_a4 = cromer_mann_a4
        self.cromer_mann_b1 = cromer_mann_b1
        self.cromer_mann_b2 = cromer_mann_b2
        self.cromer_mann_b3 = cromer_mann_b3
        self.cromer_mann_b4 = cromer_mann_b4
        self.cromer_mann_c = cromer_mann_c
        self.scat_dispersion_real = scat_dispersion_real
        self.scat_dispersion_imag = scat_dispersion_imag
        self.dispersion_source = dispersion_source
        self.scat_length_neutron = scat_length_neutron
        self.scat_source = scat_source
        self.scat_versus_stol_list = scat_versus_stol_list

        if self.is_defined:
            self.form_object

    def __repr__(self):
        ls_out = ["AtomType:"]
        ls_out.append(str(self))
        return "\n".join(ls_out)     


    @property
    def symbol(self):
        """
The code used to identify the atom species (singular or plural)
representing this atom type. Normally this code is the element
symbol. The code may be composed of any character except an
underscore with the additional proviso that digits designate an
oxidation state and must be followed by a + or - character.

Examples:
-------------
C, Cu2+, H(SDS), dummy, FeNi 

Type: char

Reference:
-------------
`iucr.org <https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Iatom_type_symbol.html>`_
        """
        return getattr(self, "__symbol")
    @symbol.setter
    def symbol(self, x):
        if x is None:
            x_in = None
        else:
            x_in = str(x)
        setattr(self, "__symbol", x_in)


    @property
    def analytical_mass(self):
        """
Mass percentage of this atom type derived from chemical analysis.
Appears in list containing _atom_type_symbol.

The permitted range is 0.0 → 100.0. 
        """
        return getattr(self, "__analytical_mass")
    @analytical_mass.setter
    def analytical_mass(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__analytical_mass", x_in)


    @property
    def description(self):
        """
A description of the atom(s) designated by this atom type. In
most cases, this will be the element name and oxidation state of
a single atom  species. For disordered or nonstoichiometric
structures it will describe a combination of atom species.


Examples:

deuterium 

0.34Fe+0.66Ni 

Appears in list containing _atom_type_symbol 
Type: char

Reference:
-------------
`iucr.org <https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Iatom_type_description.html>`_
        """
        return getattr(self, "__description")
    @description.setter
    def description(self, x):
        if x is None:
            x_in = None
        else:
            x_in = str(x)
        setattr(self, "__description", x_in)



    @property
    def number_in_cell(self):
        """
Total number of atoms of this atom type in the unit cell.


Appears in list containing _atom_type_symbol 
The permitted range is 0 -> infinity 

Type: numb
        """
        return getattr(self, "__number_in_cell")
    @number_in_cell.setter
    def number_in_cell(self, x):
        if x is None:
            x_in = None
        else:
            x_in = int(x)
        setattr(self, "__number_in_cell", x_in)


    @property
    def oxidation_number(self):
        """
Formal oxidation state of this atom type in the structure.

Appears in list containing _atom_type_symbol 
The permitted range is -8 -> 8 

Enumeration default: 0 

Type: numb
        """
        return getattr(self, "__oxidation_number")
    @oxidation_number.setter
    def oxidation_number(self, x):
        if x is None:
            x_in = None
        else:
            x_in = int(x)
        setattr(self, "__oxidation_number", x_in)


    @property
    def radius_bond(self):
        """
The effective intra- and intermolecular bonding radii in
   angstroms of this atom type.

Appears in list containing _atom_type_symbol 

The permitted range is 0.0 -> 5.0 

Type: numb
        """
        return getattr(self, "__radius_bond")
    @radius_bond.setter
    def radius_bond(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__radius_bond", x_in)


    @property
    def radius_contact(self):
        """
The effective intra- and intermolecular bonding radii in
   angstroms of this atom type.

Appears in list containing _atom_type_symbol 
The permitted range is 0.0 -> 5.0 

Type: numb

Category: atom_type
        """
        return getattr(self, "__radius_contact")
    @radius_contact.setter
    def radius_contact(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__radius_contact", x_in)


    @property
    def cromer_mann_a1(self):
        """
The Cromer-Mann scattering-factor coefficients used to calculate
the scattering factors for this atom type.

Ref: International Tables for X-ray Crystallography (1974). 
     Vol. IV, Table 2.2B
or   International Tables for Crystallography (2004). Vol. C,
     Tables 6.1.1.4 and 6.1.1.5

Appears in list containing _atom_type_symbol 
Type: numb
        """
        return getattr(self, "__cromer_mann_a1")
    @cromer_mann_a1.setter
    def cromer_mann_a1(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__cromer_mann_a1", x_in)


    @property
    def cromer_mann_a2(self):
        """
see documentation for cromer_mann_a1
        """
        return getattr(self, "__cromer_mann_a2")
    @cromer_mann_a2.setter
    def cromer_mann_a2(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__cromer_mann_a2", x_in)

    @property
    def cromer_mann_a3(self):
        """
see documentation for cromer_mann_a1
        """
        return getattr(self, "__cromer_mann_a3")
    @cromer_mann_a3.setter
    def cromer_mann_a3(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__cromer_mann_a3", x_in)


    @property
    def cromer_mann_a4(self):
        """
see documentation for cromer_mann_a1
        """
        return getattr(self, "__cromer_mann_a4")
    @cromer_mann_a4.setter
    def cromer_mann_a4(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__cromer_mann_a4", x_in)


    @property
    def cromer_mann_b1(self):
        """
see documentation for cromer_mann_a1
        """
        return getattr(self, "__cromer_mann_b1")
    @cromer_mann_b1.setter
    def cromer_mann_b1(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__cromer_mann_b1", x_in)

    @property
    def cromer_mann_b2(self):
        """
see documentation for cromer_mann_a1
        """
        return getattr(self, "__cromer_mann_b2")
    @cromer_mann_b2.setter
    def cromer_mann_b2(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__cromer_mann_b2", x_in)

    @property
    def cromer_mann_b3(self):
        """
see documentation for cromer_mann_a1
        """
        return getattr(self, "__cromer_mann_b3")
    @cromer_mann_b3.setter
    def cromer_mann_b3(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__cromer_mann_b3", x_in)


    @property
    def cromer_mann_b4(self):
        """
see documentation for cromer_mann_a1
        """
        return getattr(self, "__cromer_mann_b4")
    @cromer_mann_b4.setter
    def cromer_mann_b4(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__cromer_mann_b4", x_in)


    @property
    def cromer_mann_c(self):
        """
see documentation for cromer_mann_a1
        """
        return getattr(self, "__cromer_mann_c")
    @cromer_mann_c.setter
    def cromer_mann_c(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__cromer_mann_c", x_in)


    @property
    def scat_dispersion_real(self):
        """
The imaginary and real components of the anomalous-dispersion
scattering factor, f'' and f', in electrons for this atom type
and the radiation given in _diffrn_radiation_wavelength.


Appears in list containing _atom_type_symbol 
Enumeration default: 0.0 
Type: numb
        """
        return getattr(self, "__scat_dispersion_real")
    @scat_dispersion_real.setter
    def scat_dispersion_real(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__scat_dispersion_real", x_in)


    @property
    def scat_dispersion_imag(self):
        """
The imaginary and real components of the anomalous-dispersion
scattering factor, f'' and f', in electrons for this atom type
and the radiation given in _diffrn_radiation_wavelength.


Appears in list containing _atom_type_symbol 
Enumeration default: 0.0 
Type: numb
        """
        return getattr(self, "__scat_dispersion_imag")
    @scat_dispersion_imag.setter
    def scat_dispersion_imag(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__scat_dispersion_imag", x_in)


    @property
    def dispersion_source(self):
        """
Reference to source of real and imaginary dispersion
corrections for scattering factors used for this atom type.

Example:
'International Tables Vol. IV Table 2.3.1' 

Appears in list containing _atom_type_symbol 
Type: char
        """
        return getattr(self, "__dispersion_source")
    @dispersion_source.setter
    def dispersion_source(self, x):
        if x is None:
            x_in = None
        else:
            x_in = str(x)
        setattr(self, "__dispersion_source", x_in)


    @property
    def scat_length_neutron(self):
        """
The bound coherent scattering length in femtometres for the
atom type at the isotopic composition used for the diffraction
experiment.


Appears in list containing _atom_type_symbol 
Enumeration default: 0.0 
Type: numb

FIXME: 
now in 10**-12cm
        """
        return getattr(self, "__scat_length_neutron")
    @scat_length_neutron.setter
    def scat_length_neutron(self, x):
        if x is None:
            x_in = None
        else:
            x_in = complex(x)
        setattr(self, "__scat_length_neutron", x_in)


    @property
    def scat_source(self):
        """
Reference to source of scattering factors or scattering lengths
used for this atom type.

Example:
'International Tables Vol. IV Table 2.4.6B' 

Appears in list containing _atom_type_symbol 
Type: char
        """
        return getattr(self, "__scat_source")
    @scat_source.setter
    def scat_source(self, x):
        if x is None:
            x_in = None
        else:
            x_in = str(x)
        setattr(self, "__scat_source", x_in)


    @property
    def scat_versus_stol_list(self):
        """
A table of scattering factors as a function of sin theta over
lambda. This table should be well commented to indicate the
items present. Regularly formatted lists are strongly
recommended.

Appears in list containing _atom_type_symbol 

Type: char
        """
        return getattr(self, "__scat_versus_stol_list")
    @scat_versus_stol_list.setter
    def scat_versus_stol_list(self, x):
        if x is None:
            x_in = None
        else:
            x_in = str(x)
        setattr(self, "__scat_versus_stol_list", x_in)


    def _show_message(self, s_out: str):
        warnings.warn("***  Error ***\n"+s_out, UserWarning, stacklevel=2)

    @property
    def is_variable(self):
        res = False
        return res

    def get_variables(self):
        l_variable = []
        return l_variable


class AtomTypeL(LoopConstr):
    """
AtomTypeL
============
Data items in the ATOM_TYPE category record details about
properties of the atoms that occupy the atom sites, such as the
atomic scattering factors.

Description in cif file:
:: 

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

Reference: 
--------------
`iucr.org <https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Catom_type.html>`_
FIXME:
replace attribute 'analytical_mass' by 'analytical_mass_%' as in CoreCIF
    """
    CATEGORY_KEY = ("symbol", )
    ITEM_CLASS = AtomType
    def __init__(self, item=[], loop_name=""):
        super(AtomTypeL, self).__init__(category_key=self.CATEGORY_KEY, item_class=self.ITEM_CLASS, loop_name=loop_name)
        self.item = item

    def __repr__(self) -> str:
        ls_out = []
        ls_out.append("AtomTypeL: ")
        ls_out.append(f"{str(self):}")
        return "\n".join(ls_out)
