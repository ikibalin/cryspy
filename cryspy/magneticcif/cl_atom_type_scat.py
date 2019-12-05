"""
AtomTypeScat, AtomTypeScatL classes are defined. 
"""
__author__ = 'ikibalin'
__version__ = "2019_12_03"

import os
import numpy

from pycifstar import Global


from typing import List, Tuple
from cryspy.common.cl_item_constr import ItemConstr
from cryspy.common.cl_loop_constr import LoopConstr
from cryspy.common.cl_fitable import Fitable


class AtomTypeScat(ItemConstr):
    """
AtomTypeScat
====================

Data items in the ATOM_TYPE_SCAT category describe atomic
scattering information used in crystallographic structure studies.
This category is fully defined in the core CIF dictionary.


Description in cif:
-----------------------
::

 _atom_type_scat_symbol                       O2
 _atom_type_scat_lande                        2.0
 _atom_type_scat_neutron_magnetic_j0_A1   0.99895
 _atom_type_scat_neutron_magnetic_j0_a2  12.09652
 _atom_type_scat_neutron_magnetic_j0_B1   0.28854
 _atom_type_scat_neutron_magnetic_j0_b2   0.12914
 _atom_type_scat_neutron_magnetic_j0_C1   0.11425
 _atom_type_scat_neutron_magnetic_j0_c2  -0.22968
 _atom_type_scat_neutron_magnetic_j0_D   -0.40685

Attributes:
-----------
- symbol
- lande
- neutron_magnetic_j0_a1
- neutron_magnetic_j0_a2
- neutron_magnetic_j0_b1
- neutron_magnetic_j0_b2
- neutron_magnetic_j0_c1
- neutron_magnetic_j0_c2
- neutron_magnetic_j0_d
- neutron_magnetic_j0_e
- neutron_magnetic_j2_a1
- neutron_magnetic_j2_a2
- neutron_magnetic_j2_b1
- neutron_magnetic_j2_b2 
- neutron_magnetic_j2_c1
- neutron_magnetic_j2_c2
- neutron_magnetic_j2_d
- neutron_magnetic_j2_e
- neutron_magnetic_j4_a1
- neutron_magnetic_j4_a2
- neutron_magnetic_j4_b1
- neutron_magnetic_j4_b2
- neutron_magnetic_j4_c1
- neutron_magnetic_j4_c2
- neutron_magnetic_j4_d
- neutron_magnetic_j4_e
- neutron_magnetic_j6_a1
- neutron_magnetic_j6_a2
- neutron_magnetic_j6_b1
- neutron_magnetic_j6_b2
- neutron_magnetic_j6_c1
- neutron_magnetic_j6_c2
- neutron_magnetic_j6_d
- neutron_magnetic_j6_e
- neutron_magnetic_source

Methods:
---------
- 

Reference:
---------------
`iucr.org <ftp://ftp.iucr.org/cifdics/cif_mag_0.9.7.dic.pdf>`_
FIXME:
attribute 'symbol' is not defined in Magneti CIF dictionary
attribute 'lande' is not defined in Magneti CIF dictionary

    """
    MANDATORY_ATTRIBUTE = ("symbol", )
    OPTIONAL_ATTRIBUTE = ("neutron_magnetic_j0_a1", "neutron_magnetic_j0_a2",
                          "neutron_magnetic_j0_b1", "neutron_magnetic_j0_b2",
                          "neutron_magnetic_j0_c1", "neutron_magnetic_j0_c2",
                          "neutron_magnetic_j0_d",  "neutron_magnetic_j0_e",
                          "neutron_magnetic_j2_a1", "neutron_magnetic_j2_a2",
                          "neutron_magnetic_j2_b1", "neutron_magnetic_j2_b2",
                          "neutron_magnetic_j2_c1", "neutron_magnetic_j2_c2",
                          "neutron_magnetic_j2_d",  "neutron_magnetic_j2_e",
                          "neutron_magnetic_j4_a1", "neutron_magnetic_j4_a2",
                          "neutron_magnetic_j4_b1", "neutron_magnetic_j4_b2",
                          "neutron_magnetic_j4_c1", "neutron_magnetic_j4_c2",
                          "neutron_magnetic_j4_d",  "neutron_magnetic_j4_e",
                          "neutron_magnetic_j6_a1", "neutron_magnetic_j6_a2",
                          "neutron_magnetic_j6_b1", "neutron_magnetic_j6_b2",
                          "neutron_magnetic_j6_c1", "neutron_magnetic_j6_c2",
                          "neutron_magnetic_j6_d",  "neutron_magnetic_j6_e", "neutron_magnetic_source", "lande")
    INTERNAL_ATTRIBUTE = ()
    PREFIX = "atom_type_scat"
    def __init__(self, symbol=None, neutron_magnetic_j0_a1=None, neutron_magnetic_j0_a2=None,
                 neutron_magnetic_j0_b1=None, neutron_magnetic_j0_b2=None, neutron_magnetic_j0_c1=None, neutron_magnetic_j0_c2=None,
                 neutron_magnetic_j0_d=None, neutron_magnetic_j0_e=None, neutron_magnetic_j2_a1=None, neutron_magnetic_j2_a2=None,
                 neutron_magnetic_j2_b1=None, neutron_magnetic_j2_b2=None, neutron_magnetic_j2_c1=None, neutron_magnetic_j2_c2=None,
                 neutron_magnetic_j2_d=None, neutron_magnetic_j2_e=None, neutron_magnetic_j4_a1=None, neutron_magnetic_j4_a2=None,
                 neutron_magnetic_j4_b1=None, neutron_magnetic_j4_b2=None, neutron_magnetic_j4_c1=None, neutron_magnetic_j4_c2=None,
                 neutron_magnetic_j4_d=None, neutron_magnetic_j4_e=None, neutron_magnetic_j6_a1=None, neutron_magnetic_j6_a2=None,
                 neutron_magnetic_j6_b1=None, neutron_magnetic_j6_b2=None, neutron_magnetic_j6_c1=None, neutron_magnetic_j6_c2=None,
                 neutron_magnetic_j6_d=None, neutron_magnetic_j6_e=None, neutron_magnetic_source=None, lande=None):
        super(AtomTypeScat, self).__init__(mandatory_attribute=self.MANDATORY_ATTRIBUTE, 
                                                optional_attribute=self.OPTIONAL_ATTRIBUTE, 
                                                internal_attribute=self.INTERNAL_ATTRIBUTE,
                                                prefix=self.PREFIX)

        self.symbol, self.neutron_magnetic_source = symbol, neutron_magnetic_source
        self.neutron_magnetic_j0_a1, self.neutron_magnetic_j0_a2 = neutron_magnetic_j0_a1, neutron_magnetic_j0_a2
        self.neutron_magnetic_j0_b1, self.neutron_magnetic_j0_b2 = neutron_magnetic_j0_b1, neutron_magnetic_j0_b2
        self.neutron_magnetic_j0_c1, self.neutron_magnetic_j0_c2 = neutron_magnetic_j0_c1, neutron_magnetic_j0_c2
        self.neutron_magnetic_j0_d, self.neutron_magnetic_j0_e = neutron_magnetic_j0_d, neutron_magnetic_j0_e
        self.neutron_magnetic_j2_a1, self.neutron_magnetic_j2_a2 = neutron_magnetic_j2_a1, neutron_magnetic_j2_a2
        self.neutron_magnetic_j2_b1, self.neutron_magnetic_j2_b2 = neutron_magnetic_j2_b1, neutron_magnetic_j2_b2
        self.neutron_magnetic_j2_c1, self.neutron_magnetic_j2_c2 = neutron_magnetic_j2_c1, neutron_magnetic_j2_c2
        self.neutron_magnetic_j2_d, self.neutron_magnetic_j2_e = neutron_magnetic_j2_d, neutron_magnetic_j2_e
        self.neutron_magnetic_j4_a1, self.neutron_magnetic_j4_a2 = neutron_magnetic_j4_a1, neutron_magnetic_j4_a2
        self.neutron_magnetic_j4_b1, self.neutron_magnetic_j4_b2 = neutron_magnetic_j4_b1, neutron_magnetic_j4_b2
        self.neutron_magnetic_j4_c1, self.neutron_magnetic_j4_c2 = neutron_magnetic_j4_c1, neutron_magnetic_j4_c2
        self.neutron_magnetic_j4_d, self.neutron_magnetic_j4_e = neutron_magnetic_j4_d, neutron_magnetic_j4_e
        self.neutron_magnetic_j6_a1, self.neutron_magnetic_j6_a2 = neutron_magnetic_j6_a1, neutron_magnetic_j6_a2
        self.neutron_magnetic_j6_b1, self.neutron_magnetic_j6_b2 = neutron_magnetic_j6_b1, neutron_magnetic_j6_b2
        self.neutron_magnetic_j6_c1, self.neutron_magnetic_j6_c2 = neutron_magnetic_j6_c1, neutron_magnetic_j6_c2
        self.neutron_magnetic_j6_d, self.neutron_magnetic_j6_e = neutron_magnetic_j6_d, neutron_magnetic_j6_e
        self.lande = lande

        if self.is_defined:
            self.form_object

    @property
    def symbol(self) -> str:
        """
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
    def neutron_magnetic_j0_a1(self) -> str:
        """
First, the parameters are used directly to approximate spatial averages of 
spherical Bessel functions over the electronic wave functions of unpaired 
electrons of the given atom type as a function of s = sin(θ)/λ. 

.. math::
 <jn(s)> = [A1 exp(−a2s2) + B1 exp(−b2s2) + C1 exp(−c2s2) + D] ×
 [1 if n = 0, s^{2} ifn = 2,4,6] 

The <jn(s)> are then combined to determine the spin and orbital contributions 
to the magnetic form factor of the atom. The ‘e’ parameter is a measure of error 
in the approximation.         

Reference:
---------------
`iucr.org <ftp://ftp.iucr.org/cifdics/cif_mag_0.9.7.dic.pdf>`_
        """
        return getattr(self, "__neutron_magnetic_j0_a1")
    @neutron_magnetic_j0_a1.setter
    def neutron_magnetic_j0_a1(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__neutron_magnetic_j0_a1", x_in)


    @property
    def neutron_magnetic_j0_a2(self) -> str:
        """
see documentation for neutron_magnetic_j0_a1
        """
        return getattr(self, "__neutron_magnetic_j0_a2")
    @neutron_magnetic_j0_a2.setter
    def neutron_magnetic_j0_a2(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__neutron_magnetic_j0_a2", x_in)

    @property
    def neutron_magnetic_j0_b1(self) -> str:
        """
see documentation for neutron_magnetic_j0_a1
        """
        return getattr(self, "__neutron_magnetic_j0_b1")
    @neutron_magnetic_j0_b1.setter
    def neutron_magnetic_j0_b1(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__neutron_magnetic_j0_b1", x_in)

    @property
    def neutron_magnetic_j0_b2(self) -> str:
        """
see documentation for neutron_magnetic_j0_a1
        """
        return getattr(self, "__neutron_magnetic_j0_b2")
    @neutron_magnetic_j0_b2.setter
    def neutron_magnetic_j0_b2(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__neutron_magnetic_j0_b2", x_in)


    @property
    def neutron_magnetic_j0_c1(self) -> str:
        """
see documentation for neutron_magnetic_j0_a1
        """
        return getattr(self, "__neutron_magnetic_j0_c1")
    @neutron_magnetic_j0_c1.setter
    def neutron_magnetic_j0_c1(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__neutron_magnetic_j0_c1", x_in)

    @property
    def neutron_magnetic_j0_c2(self) -> str:
        """
see documentation for neutron_magnetic_j0_a1
        """
        return getattr(self, "__neutron_magnetic_j0_c2")
    @neutron_magnetic_j0_c2.setter
    def neutron_magnetic_j0_c2(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__neutron_magnetic_j0_c2", x_in)

    @property
    def neutron_magnetic_j0_d(self) -> str:
        """
see documentation for neutron_magnetic_j0_a1
        """
        return getattr(self, "__neutron_magnetic_j0_d")
    @neutron_magnetic_j0_d.setter
    def neutron_magnetic_j0_d(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__neutron_magnetic_j0_d", x_in)

    @property
    def neutron_magnetic_j0_e(self) -> str:
        """
see documentation for neutron_magnetic_j0_a1
        """
        return getattr(self, "__neutron_magnetic_j0_e")
    @neutron_magnetic_j0_e.setter
    def neutron_magnetic_j0_e(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__neutron_magnetic_j0_e", x_in)


    @property
    def neutron_magnetic_j2_a1(self) -> str:
        """
see documentation for neutron_magnetic_j0_a1
        """
        return getattr(self, "__neutron_magnetic_j2_a1")
    @neutron_magnetic_j2_a1.setter
    def neutron_magnetic_j2_a1(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__neutron_magnetic_j2_a1", x_in)


    @property
    def neutron_magnetic_j2_a2(self) -> str:
        """
see documentation for neutron_magnetic_j0_a1
        """
        return getattr(self, "__neutron_magnetic_j2_a2")
    @neutron_magnetic_j2_a2.setter
    def neutron_magnetic_j2_a2(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__neutron_magnetic_j2_a2", x_in)

    @property
    def neutron_magnetic_j2_b1(self) -> str:
        """
see documentation for neutron_magnetic_j0_a1
        """
        return getattr(self, "__neutron_magnetic_j2_b1")
    @neutron_magnetic_j2_b1.setter
    def neutron_magnetic_j2_b1(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__neutron_magnetic_j2_b1", x_in)

    @property
    def neutron_magnetic_j2_b2(self) -> str:
        """
see documentation for neutron_magnetic_j0_a1
        """
        return getattr(self, "__neutron_magnetic_j2_b2")
    @neutron_magnetic_j2_b2.setter
    def neutron_magnetic_j2_b2(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__neutron_magnetic_j2_b2", x_in)


    @property
    def neutron_magnetic_j2_c1(self) -> str:
        """
see documentation for neutron_magnetic_j0_a1
        """
        return getattr(self, "__neutron_magnetic_j2_c1")
    @neutron_magnetic_j2_c1.setter
    def neutron_magnetic_j2_c1(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__neutron_magnetic_j2_c1", x_in)

    @property
    def neutron_magnetic_j2_c2(self) -> str:
        """
see documentation for neutron_magnetic_j0_a1
        """
        return getattr(self, "__neutron_magnetic_j2_c2")
    @neutron_magnetic_j2_c2.setter
    def neutron_magnetic_j2_c2(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__neutron_magnetic_j2_c2", x_in)

    @property
    def neutron_magnetic_j2_d(self) -> str:
        """
see documentation for neutron_magnetic_j0_a1
        """
        return getattr(self, "__neutron_magnetic_j2_d")
    @neutron_magnetic_j2_d.setter
    def neutron_magnetic_j2_d(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__neutron_magnetic_j2_d", x_in)

    @property
    def neutron_magnetic_j2_e(self) -> str:
        """
see documentation for neutron_magnetic_j0_a1
        """
        return getattr(self, "__neutron_magnetic_j2_e")
    @neutron_magnetic_j2_e.setter
    def neutron_magnetic_j2_e(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__neutron_magnetic_j2_e", x_in)



    @property
    def neutron_magnetic_j4_a1(self) -> str:
        """
see documentation for neutron_magnetic_j0_a1
        """
        return getattr(self, "__neutron_magnetic_j4_a1")
    @neutron_magnetic_j4_a1.setter
    def neutron_magnetic_j4_a1(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__neutron_magnetic_j4_a1", x_in)


    @property
    def neutron_magnetic_j4_a2(self) -> str:
        """
see documentation for neutron_magnetic_j0_a1
        """
        return getattr(self, "__neutron_magnetic_j4_a2")
    @neutron_magnetic_j4_a2.setter
    def neutron_magnetic_j4_a2(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__neutron_magnetic_j4_a2", x_in)

    @property
    def neutron_magnetic_j4_b1(self) -> str:
        """
see documentation for neutron_magnetic_j0_a1
        """
        return getattr(self, "__neutron_magnetic_j4_b1")
    @neutron_magnetic_j4_b1.setter
    def neutron_magnetic_j4_b1(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__neutron_magnetic_j4_b1", x_in)

    @property
    def neutron_magnetic_j4_b2(self) -> str:
        """
see documentation for neutron_magnetic_j0_a1
        """
        return getattr(self, "__neutron_magnetic_j4_b2")
    @neutron_magnetic_j4_b2.setter
    def neutron_magnetic_j4_b2(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__neutron_magnetic_j4_b2", x_in)


    @property
    def neutron_magnetic_j4_c1(self) -> str:
        """
see documentation for neutron_magnetic_j0_a1
        """
        return getattr(self, "__neutron_magnetic_j4_c1")
    @neutron_magnetic_j4_c1.setter
    def neutron_magnetic_j4_c1(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__neutron_magnetic_j4_c1", x_in)

    @property
    def neutron_magnetic_j4_c2(self) -> str:
        """
see documentation for neutron_magnetic_j0_a1
        """
        return getattr(self, "__neutron_magnetic_j4_c2")
    @neutron_magnetic_j4_c2.setter
    def neutron_magnetic_j4_c2(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__neutron_magnetic_j4_c2", x_in)

    @property
    def neutron_magnetic_j4_d(self) -> str:
        """
see documentation for neutron_magnetic_j0_a1
        """
        return getattr(self, "__neutron_magnetic_j4_d")
    @neutron_magnetic_j4_d.setter
    def neutron_magnetic_j4_d(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__neutron_magnetic_j4_d", x_in)

    @property
    def neutron_magnetic_j4_e(self) -> str:
        """
see documentation for neutron_magnetic_j0_a1
        """
        return getattr(self, "__neutron_magnetic_j4_e")
    @neutron_magnetic_j4_e.setter
    def neutron_magnetic_j4_e(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__neutron_magnetic_j4_e", x_in)




    @property
    def neutron_magnetic_j6_a1(self) -> str:
        """
see documentation for neutron_magnetic_j0_a1
        """
        return getattr(self, "__neutron_magnetic_j6_a1")
    @neutron_magnetic_j6_a1.setter
    def neutron_magnetic_j6_a1(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__neutron_magnetic_j6_a1", x_in)


    @property
    def neutron_magnetic_j6_a2(self) -> str:
        """
see documentation for neutron_magnetic_j0_a1
        """
        return getattr(self, "__neutron_magnetic_j6_a2")
    @neutron_magnetic_j6_a2.setter
    def neutron_magnetic_j6_a2(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__neutron_magnetic_j6_a2", x_in)

    @property
    def neutron_magnetic_j6_b1(self) -> str:
        """
see documentation for neutron_magnetic_j0_a1
        """
        return getattr(self, "__neutron_magnetic_j6_b1")
    @neutron_magnetic_j6_b1.setter
    def neutron_magnetic_j6_b1(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__neutron_magnetic_j6_b1", x_in)

    @property
    def neutron_magnetic_j6_b2(self) -> str:
        """
see documentation for neutron_magnetic_j0_a1
        """
        return getattr(self, "__neutron_magnetic_j6_b2")
    @neutron_magnetic_j6_b2.setter
    def neutron_magnetic_j6_b2(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__neutron_magnetic_j6_b2", x_in)


    @property
    def neutron_magnetic_j6_c1(self) -> str:
        """
see documentation for neutron_magnetic_j0_a1
        """
        return getattr(self, "__neutron_magnetic_j6_c1")
    @neutron_magnetic_j6_c1.setter
    def neutron_magnetic_j6_c1(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__neutron_magnetic_j6_c1", x_in)

    @property
    def neutron_magnetic_j6_c2(self) -> str:
        """
see documentation for neutron_magnetic_j0_a1
        """
        return getattr(self, "__neutron_magnetic_j6_c2")
    @neutron_magnetic_j6_c2.setter
    def neutron_magnetic_j6_c2(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__neutron_magnetic_j6_c2", x_in)

    @property
    def neutron_magnetic_j6_d(self) -> str:
        """
see documentation for neutron_magnetic_j0_a1
        """
        return getattr(self, "__neutron_magnetic_j6_d")
    @neutron_magnetic_j6_d.setter
    def neutron_magnetic_j6_d(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__neutron_magnetic_j6_d", x_in)

    @property
    def neutron_magnetic_j6_e(self) -> str:
        """
see documentation for neutron_magnetic_j0_a1
        """
        return getattr(self, "__neutron_magnetic_j6_e")
    @neutron_magnetic_j6_e.setter
    def neutron_magnetic_j6_e(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__neutron_magnetic_j6_e", x_in)


    @property
    def neutron_magnetic_source(self) -> str:
        """
Reference to the source of magnetic neutron scattering factors for 
a given atom type. 

Analogous tags: coreCIF: _atom_site.scat_source 
Example: 
;
International Tables for Crystallography (2006). Vol. C, Section 4.4.5. 
;
        """
        return getattr(self, "__neutron_magnetic_source")
    @neutron_magnetic_source.setter
    def neutron_magnetic_source(self, x):
        if x is None:
            x_in = None
        else:
            x_in = str(x)
        setattr(self, "__neutron_magnetic_source", x_in)

    @property
    def lande(self):
        """
Lande factor L:
form factor = <j0> + (2/L - 1) * <j2>

        """
        return getattr(self, "__lande")
    @lande.setter
    def lande(self, x):
        if x is None:
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__lande", x_in)

    def _show_message(self, s_out: str):
        warnings.warn("***  Error ***\n"+s_out, UserWarning, stacklevel=2)


    @property
    def is_variable(self):
        res = any([self.lande.refinement])
        return res

    def get_variables(self):
        l_variable = []
        if self.lande.refinement: l_variable.append(self.lande)
        return l_variable



class AtomTypeScatL(LoopConstr):
    """
AtomTypeScatL
====================


Data items in the ATOM_TYPE_SCAT category describe atomic
scattering information used in crystallographic structure studies.
This category is fully defined in the core CIF dictionary.

Description in cif:
-----------------------
::

 loop_
 _atom_type_scat_symbol
 _atom_type_scat_neutron_magnetic_j0_A1
 _atom_type_scat_neutron_magnetic_j0_a2
 _atom_type_scat_neutron_magnetic_j0_B1
 _atom_type_scat_neutron_magnetic_j0_b2
 _atom_type_scat_neutron_magnetic_j0_C1
 _atom_type_scat_neutron_magnetic_j0_c2
 _atom_type_scat_neutron_magnetic_j0_D
  O2   0.99895  12.09652  0.28854  0.12914  0.11425 -0.22968 -0.40685 
  N2   1.00581  13.37218 -0.05868  0.07792 -0.00444  0.01678  0.05146

Reference:
---------------
`iucr.org <ftp://ftp.iucr.org/cifdics/cif_mag_0.9.7.dic.pdf>`_
    """
    CATEGORY_KEY = ("symbol", )
    ITEM_CLASS = AtomTypeScat
    def __init__(self, item = [], loop_name=""):
        super(AtomTypeScatL, self).__init__(category_key=self.CATEGORY_KEY, item_class=self.ITEM_CLASS, loop_name=loop_name)
        self.item = item

    def __repr__(self) -> str:
        ls_out = []
        ls_out.append("AtomTypeScatL: ")
        ls_out.append(f"{str(self):}")
        return "\n".join(ls_out)
