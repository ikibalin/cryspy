__author__ = 'ikibalin'
__version__ = "2020_04_27"
import os
import math
import numpy
import re
from pycifstar import Global

import warnings
from typing import List, Tuple
from cryspy.common.cl_item_constr import ItemConstr
from cryspy.common.cl_loop_constr import LoopConstr
from cryspy.common.cl_fitable import Fitable
import cryspy.corecif.CONSTANTS_AND_FUNCTIONS as CONSTANTS_AND_FUNCTIONS 



class AtomElectronConfiguration(ItemConstr):
    """
Define the electron configuration for atoms.
The elcetron configuration is separated on the core part which is spherical
and non-spherical valence part. 

The orientation and population of valence part is described in the object
AtomRhoOrbitalValence of cryspy library.

Description in cif file::

_atom_electron_configuration_label Co1 
_atom_electron_configuration_core "1s2 2s2 2p6 3s2 3p6"
_atom_electron_configuration_valence "4s 3d"
    """    
    MANDATORY_ATTRIBUTE = ("label", )
    OPTIONAL_ATTRIBUTE = ("core", "valence")
    INTERNAL_ATTRIBUTE = ()
    PREFIX = "atom_electron_configuration"
    def __init__(self, label=None, core=None, valence=None):
        super(AtomElectronConfiguration, self).__init__(mandatory_attribute=self.MANDATORY_ATTRIBUTE, 
                                            optional_attribute=self.OPTIONAL_ATTRIBUTE, 
                                            internal_attribute=self.INTERNAL_ATTRIBUTE,
                                            prefix=self.PREFIX)

        self.label = label
        self.core = core
        self.valence = valence

        if self.is_defined:
            self.form_object

    @property
    def label(self)->str:
        return getattr(self, "__label")
    @label.setter
    def label(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = str(x)
        setattr(self, "__label", x_in)

    @property
    def core(self)->str:
        return getattr(self, "__core")
    @core.setter
    def core(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = str(x)
        setattr(self, "__core", x_in)

    @property
    def valence(self)->str:
        return getattr(self, "__valence")
    @valence.setter
    def valence(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = str(x)
        setattr(self, "__valence", x_in)

    def get_core_shells_populations(self):
        """
Give list of shels and their population for core part
        """
        l_res = []
        s_core = self.core
        if s_core is not None:
            for _s_shell in s_core.split():
                m = re.match(r"(\d*)([s,p,d,f])(\d*)", _s_shell)
                if m is not None:
                    res = (f"{m.group(1):}{m.group(2):}", int(m.group(3)))
                    l_res.append(res)
        return l_res

    def get_valence_shells(self):
        """
Give list of shels for valence part
        """
        l_res = []
        s_valence = self.valence
        if s_valence is not None:
            for _s_shell in s_valence.split():
                m = re.match(r"(\d*)([s,p,d,f])", _s_shell)
                if m is not None:
                    res = f"{m.group(1):}{m.group(2):}"
                    l_res.append(res)
        return l_res


class AtomElectronConfigurationL(LoopConstr):
    """
Define the electron configuration for atoms.
The elcetron configuration is separated on the core part which is spherical
and non-spherical valence part. 

The orientation and population of valence part is described in the object
AtomRhoOrbitalValence of cryspy library.

Description in cif file::

loop_
_atom_electron_configuration_label
_atom_electron_configuration_core
_atom_electron_configuration_valence
Co1 "1s2 2s2 2p6 3s2 3p6" .
Co2 "1s2 2s2 2p6 3s2 3p6" .
 V1 "1s2 2s2 2p6 3s2 3p6" .
 O1 "1s2" "2s 2p"
 O2 "1s2" "2s 2p"
 O3 "1s2" "2s 2p"

    """
    CATEGORY_KEY = ("label", )
    ITEM_CLASS = AtomElectronConfiguration
    def __init__(self, item=[], loop_name=""):
        super(AtomElectronConfigurationL, self).__init__(category_key=self.CATEGORY_KEY, item_class=self.ITEM_CLASS, loop_name=loop_name)
        self.item = item

    def get_core_shells_populations(self):
        """
Give list of shels and their population for core part of all atoms 
        """
        ll_res = [_item.get_core_shells_populations() for _item in self.item]
        return ll_res

