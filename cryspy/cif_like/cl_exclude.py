__author__ = 'ikibalin'
__version__ = "2019_12_10"
import os
import numpy

import warnings
from typing import List, Tuple

from cryspy.common.cl_item_constr import ItemConstr
from cryspy.common.cl_loop_constr import LoopConstr
from cryspy.common.cl_fitable import Fitable

class Exclude(ItemConstr):
    """
Description in cif file::

 _exclude_id            1
 _exclude_ttheta_min  4.0  
 _exclude_ttheta_max 12.0
 _exclude_phi_min  4.0  
 _exclude_phi_max 12.0
    """
    MANDATORY_ATTRIBUTE = ("ttheta_min", "ttheta_max")
    OPTIONAL_ATTRIBUTE = ("id", "phi_min", "phi_max")
    RELATED_CIF_MANDATORY_ATTRIBUTE = ("2theta_min", "2theta_max")
    RELATED_CIF_OPTIONAL_ATTRIBUTE = ("id", "phi_min", "phi_max")
    INTERNAL_ATTRIBUTE = ()
    PREFIX = "exclude"
    def __init__(self, id=None, ttheta_min=None, ttheta_max=None, phi_min=None, phi_max=None):
        super(Exclude, self).__init__(mandatory_attribute=self.MANDATORY_ATTRIBUTE, 
                                        optional_attribute=self.OPTIONAL_ATTRIBUTE, 
                                        internal_attribute=self.INTERNAL_ATTRIBUTE,
                                        prefix=self.PREFIX)
        self.id = id
        self.ttheta_min = ttheta_min
        self.ttheta_max = ttheta_max

        if self.is_defined:
            self.form_object

    @property
    def id(self):
        return getattr(self, "__id")
    @id.setter
    def id(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = str(x)
        setattr(self, "__id", x_in)

    @property
    def ttheta_min(self):
        return getattr(self, "__ttheta_min")
    @ttheta_min.setter
    def ttheta_min(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__ttheta_min", x_in)

    @property
    def ttheta_max(self):
        return getattr(self, "__ttheta_max")
    @ttheta_max.setter
    def ttheta_max(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__ttheta_max", x_in)

    @property
    def phi_min(self):
        return getattr(self, "__phi_min")
    @phi_min.setter
    def phi_min(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__phi_min", x_in)

    @property
    def phi_max(self):
        return getattr(self, "__phi_max")
    @phi_max.setter
    def phi_max(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__phi_max", x_in)

    def __repr__(self) -> str:
        ls_out = []
        ls_out.append(f"Exclude:\n{str(self):}")
        return "\n".join(ls_out)

    def _show_message(self, s_out: str):
        warnings.warn("***  Error ***\n"+s_out, UserWarning, stacklevel=2)

    @property
    def is_variable(self) -> bool:
        return False
    
    def get_variables(self) -> List:
        return []

class ExcludeL(LoopConstr):
    """
Description in cif file::

 loop_
 _exclude_id
 _exclude_ttheta_min
 _exclude_ttheta_max
  1   4.0  12.0 
  2  30.0  45.0 
  3  58.0  63.0 
    """
    CATEGORY_KEY = ("id", )
    ITEM_CLASS = Exclude
    def __init__(self, item=[], loop_name=""):
        super(ExcludeL, self).__init__(category_key=self.CATEGORY_KEY, item_class=self.ITEM_CLASS, loop_name=loop_name)
        self.item = item

    def __repr__(self) -> str:
        ls_out = []
        ls_out.append("ExcludeL: ")
        ls_out.append(f"{str(self):}")
        return "\n".join(ls_out)

