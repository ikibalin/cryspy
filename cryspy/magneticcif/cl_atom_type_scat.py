__author__ = 'ikibalin'
__version__ = "2019_12_03"

import os
import numpy

from pycifstar import Global


from typing import List, Tuple, NoReturn
from cryspy.common.cl_item_constr import ItemConstr
from cryspy.common.cl_loop_constr import LoopConstr
from cryspy.common.cl_fitable import Fitable

import cryspy.magneticcif.CONSTANTS_AND_FUNCTIONS as CONSTANTS_AND_FUNCTIONS


class AtomTypeScat(ItemConstr):
    """
Data items in the ATOM_TYPE_SCAT category describe atomic
scattering information used in crystallographic structure studies.
This category is fully defined in the core CIF dictionary.


Description in cif::

 _atom_type_scat_symbol                        O2
 _atom_type_scat_neutron_magnetic_j0_A1   0.99895
 _atom_type_scat_neutron_magnetic_j0_a2  12.09652
 _atom_type_scat_neutron_magnetic_j0_B1   0.28854
 _atom_type_scat_neutron_magnetic_j0_b2   0.12914
 _atom_type_scat_neutron_magnetic_j0_C1   0.11425
 _atom_type_scat_neutron_magnetic_j0_c2  -0.22968
 _atom_type_scat_neutron_magnetic_j0_D   -0.40685

`<ftp://ftp.iucr.org/cifdics/cif_mag_0.9.7.dic.pdf>`_

:FIXME: attribute 'symbol' is not defined in Magnetic CIF dictionary
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
                          "neutron_magnetic_j6_d",  "neutron_magnetic_j6_e", "neutron_magnetic_source")
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
                 neutron_magnetic_j6_d=None, neutron_magnetic_j6_e=None, neutron_magnetic_source=None):
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

        if self.is_defined:
            self.form_object

    @property
    def symbol(self) -> str:
        """
        """
        return getattr(self, "__symbol")
    @symbol.setter
    def symbol(self, x):
        if ((x is None) | (x == ".")):
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

   \\left\\langle j_{n}(s)\\right\\rangle = \\left( A_{1} e^{−a_{2}s^{2}} + B_{1} e^{−b_{2}s^{2}} + C_{1} e^{−c_{2}s^{2}}  + D \\right) \\times
   \\left[\\left. 1 \\right|_{n = 0} \\vee \\left.s^{2}\\right|_{n = 2, 4, 6}\\right] 

The :math:`\\left\\langle j_{n}(s)\\right\\rangle` are then combined to determine the spin and orbital contributions 
to the magnetic form factor of the atom. The ‘e’ parameter is a measure of error 
in the approximation.         

`<ftp://ftp.iucr.org/cifdics/cif_mag_0.9.7.dic.pdf>`_
        """
        return getattr(self, "__neutron_magnetic_j0_a1")
    @neutron_magnetic_j0_a1.setter
    def neutron_magnetic_j0_a1(self, x):
        if ((x is None) | (x == ".")):
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
        if ((x is None) | (x == ".")):
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
        if ((x is None) | (x == ".")):
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
        if ((x is None) | (x == ".")):
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
        if ((x is None) | (x == ".")):
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
        if ((x is None) | (x == ".")):
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
        if ((x is None) | (x == ".")):
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
        if ((x is None) | (x == ".")):
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
        if ((x is None) | (x == ".")):
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
        if ((x is None) | (x == ".")):
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
        if ((x is None) | (x == ".")):
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
        if ((x is None) | (x == ".")):
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
        if ((x is None) | (x == ".")):
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
        if ((x is None) | (x == ".")):
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
        if ((x is None) | (x == ".")):
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
        if ((x is None) | (x == ".")):
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
        if ((x is None) | (x == ".")):
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
        if ((x is None) | (x == ".")):
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
        if ((x is None) | (x == ".")):
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
        if ((x is None) | (x == ".")):
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
        if ((x is None) | (x == ".")):
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
        if ((x is None) | (x == ".")):
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
        if ((x is None) | (x == ".")):
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
        if ((x is None) | (x == ".")):
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
        if ((x is None) | (x == ".")):
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
        if ((x is None) | (x == ".")):
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
        if ((x is None) | (x == ".")):
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
        if ((x is None) | (x == ".")):
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
        if ((x is None) | (x == ".")):
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
        if ((x is None) | (x == ".")):
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
        if ((x is None) | (x == ".")):
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
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__neutron_magnetic_j6_e", x_in)


    @property
    def neutron_magnetic_source(self) -> str:
        """
Reference to the source of magnetic neutron scattering factors for 
a given atom type. 

Example: 
;
International Tables for Crystallography (2006). Vol. C, Section 4.4.5. 
;
        """
        return getattr(self, "__neutron_magnetic_source")
    @neutron_magnetic_source.setter
    def neutron_magnetic_source(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = str(x)
        setattr(self, "__neutron_magnetic_source", x_in)

    def __repr__(self) -> str:
        ls_out = []
        ls_out.append("AtomTypeScat: ")
        ls_out.append(f"{str(self):}")
        return "\n".join(ls_out)

    def _show_message(self, s_out: str):
        warnings.warn("***  Error ***\n"+s_out, UserWarning, stacklevel=2)

    @property
    def is_variable(self):
        res = any([])
        return res

    def get_variables(self):
        l_variable = []
        return l_variable

    def calc_form_factor(self, sthovl, lande=2., kappa=1., flag_only_orbital=False):
        """
Calculate magnetic form factor in frame of Spherical model (Int.Tabl.C.p.592)

:LFactor:  is Lande factor
:lsthovl:  is list :math:`sin(\\theta)/\\lambda` in :math:`\\A^{-1}`

Calculation of magnetic form factor <j0>, <j2>, <j4>, <j6>

:coeff: is a list [A, a, B, b, C, c, D] at n = 0, 2, 4, 6

according `<https://journals.aps.org/prb/pdf/10.1103/PhysRevB.79.140405>`
mismatch with international tables where (1.0-2.0/np_factor_lande)
        """
        #not sure about kappa, it is here just for test, by default it is 1.0
        j0_av = self.calc_j0(sthovl, kappa=kappa)    
        j2_av = self.calc_j2(sthovl, kappa=kappa)    

        if flag_only_orbital:
            form_factor = (2.0/float(lande)-1.0)*j2_av
        else:
            form_factor = j0_av+(2.0/float(lande)-1.0)*j2_av
        return form_factor


    def calc_j0(self, sthovl, kappa:float=1):
        j0_A = self.neutron_magnetic_j0_a1
        j0_a = self.neutron_magnetic_j0_a2
        j0_B = self.neutron_magnetic_j0_b1
        j0_b = self.neutron_magnetic_j0_b2
        j0_C = self.neutron_magnetic_j0_c1
        j0_c = self.neutron_magnetic_j0_c2
        j0_D = self.neutron_magnetic_j0_d
        _h = (sthovl/float(kappa))**2
        j0_av = (j0_A*numpy.exp(-j0_a*_h)+
                 j0_B*numpy.exp(-j0_b*_h)+
                 j0_C*numpy.exp(-j0_c*_h)+j0_D)
        return j0_av

    def calc_j2(self, sthovl, kappa:float=1):
        j2_A = self.neutron_magnetic_j2_a1
        j2_a = self.neutron_magnetic_j2_a2
        j2_B = self.neutron_magnetic_j2_b1
        j2_b = self.neutron_magnetic_j2_b2
        j2_C = self.neutron_magnetic_j2_c1
        j2_c = self.neutron_magnetic_j2_c2
        j2_D = self.neutron_magnetic_j2_d     
        _h = (sthovl/float(kappa))**2
        j2_av = (j2_A*numpy.exp(-j2_a*_h)+
                 j2_B*numpy.exp(-j2_b*_h)+
                 j2_C*numpy.exp(-j2_c*_h)+j2_D)*_h
        return j2_av

    def load_by_symbol(self, symbol:str=None) -> NoReturn:
        if symbol is None:
            symbol = self.symbol
        else:
            self.symbol = symbol
        j0_A1, j0_a2, j0_B1, j0_b2, j0_C1, j0_c2, j0_D, j2_A1, j2_a2, j2_B1, j2_b2, j2_C1, j2_c2, j2_D = CONSTANTS_AND_FUNCTIONS.get_j0_j2_by_symbol(symbol)
        self.neutron_magnetic_j0_a1 = j0_A1 
        self.neutron_magnetic_j0_a2 = j0_a2
        self.neutron_magnetic_j0_b1 = j0_B1 
        self.neutron_magnetic_j0_b2 = j0_b2
        self.neutron_magnetic_j0_c1 = j0_C1 
        self.neutron_magnetic_j0_c2 = j0_c2
        self.neutron_magnetic_j0_d = j0_D
        self.neutron_magnetic_j2_a1 = j2_A1 
        self.neutron_magnetic_j2_a2 = j2_a2
        self.neutron_magnetic_j2_b1 = j2_B1 
        self.neutron_magnetic_j2_b2 = j2_b2
        self.neutron_magnetic_j2_c1 = j2_C1 
        self.neutron_magnetic_j2_c2 = j2_c2
        self.neutron_magnetic_j2_d = j2_D

    @classmethod
    def form_by_symbol(cls, symbol:str):
        j0_A1, j0_a2, j0_B1, j0_b2, j0_C1, j0_c2, j0_D, j2_A1, j2_a2, j2_B1, j2_b2, j2_C1, j2_c2, j2_D = CONSTANTS_AND_FUNCTIONS.get_j0_j2_by_symbol(symbol)
        item = cls(symbol=symbol, 
            neutron_magnetic_j0_a1=j0_A1, neutron_magnetic_j0_a2=j0_a2,
            neutron_magnetic_j0_b1=j0_B1, neutron_magnetic_j0_b2=j0_b2,
            neutron_magnetic_j0_c1=j0_C1, neutron_magnetic_j0_c2=j0_c2,
            neutron_magnetic_j0_d=j0_D,
            neutron_magnetic_j2_a1=j2_A1, neutron_magnetic_j2_a2=j2_a2,
            neutron_magnetic_j2_b1=j2_B1, neutron_magnetic_j2_b2=j2_b2,
            neutron_magnetic_j2_c1=j2_C1, neutron_magnetic_j2_c2=j2_c2,
            neutron_magnetic_j2_d=j2_D)
        return item

class AtomTypeScatL(LoopConstr):
    """
Data items in the ATOM_TYPE_SCAT category describe atomic
scattering information used in crystallographic structure studies.
This category is fully defined in the core CIF dictionary.

Description in cif::

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

`<ftp://ftp.iucr.org/cifdics/cif_mag_0.9.7.dic.pdf>`_
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

    def calc_form_factor(self, sthovl, l_lande: List, l_kappa: List, flag_only_orbital=False):
        form_factor = [_item.calc_form_factor(sthovl, float(lande), float(kappa), flag_only_orbital=flag_only_orbital) 
                       for _item, lande, kappa in zip(self.item, l_lande, l_kappa)]
        return numpy.array(form_factor, dtype=float)

    @classmethod
    def form_by_symbols(cls, symbols:List):
        item = []
        item_class = self.ITEM_CLASS
        for symbol in symbols:
            item.append(item_class.form_by_symbol(symbol))
        obj = cls(item=item)
        return obj
