__author__ = 'ikibalin'
__version__ = "2019_12_10"
import os
import numpy

import warnings
from typing import List, Tuple

from cryspy.common.cl_item_constr import ItemConstr
from cryspy.common.cl_loop_constr import LoopConstr
from cryspy.common.cl_fitable import Fitable

class PdExclude(ItemConstr):
    """
Description in cif file::

 _pd_exclude_id            1
 _pd_exclude_ttheta_min  4.0  
 _pd_exclude_ttheta_max 12.0
    """
    MANDATORY_ATTRIBUTE = ("id", "ttheta_min", "ttheta_max")
    OPTIONAL_ATTRIBUTE = ()
    INTERNAL_ATTRIBUTE = ()
    PREFIX = "pd_exclude"
    def __init__(self, id=None, ttheta_min=None, ttheta_max=None):
        super(PdExclude, self).__init__(mandatory_attribute=self.MANDATORY_ATTRIBUTE, 
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
        if x is None:
            x_in = None
        else:
            x_in = str(x)
        setattr(self, "__id", x_in)

    @property
    def ttheta_min(self):
        return getattr(self, "__ttheta_min")
    @ttheta_min.setter
    def ttheta_min(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__ttheta_min", x_in)

    @property
    def ttheta_max(self):
        return getattr(self, "__ttheta_max")
    @ttheta_max.setter
    def ttheta_max(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__ttheta_max", x_in)

    def __repr__(self) -> str:
        ls_out = []
        ls_out.append(f"PdExclude:\n{str(self):}")
        return "\n".join(ls_out)

    def _show_message(self, s_out: str):
        warnings.warn("***  Error ***\n"+s_out, UserWarning, stacklevel=2)

    @property
    def is_variable(self) -> bool:
        return False
    
    def get_variables(self) -> List:
        return []

class PdExcludeL(LoopConstr):
    """
Description in cif file::

 loop_
 _pd_exclude_id
 _pd_exclude_ttheta_min
 _pd_exclude_ttheta_max
  1   4.0  12.0 
  2  30.0  45.0 
  3  58.0  63.0 
    """
    CATEGORY_KEY = ("id", )
    ITEM_CLASS = PdExclude
    def __init__(self, item=[], loop_name=""):
        super(PdExcludeL, self).__init__(category_key=self.CATEGORY_KEY, item_class=self.ITEM_CLASS, loop_name=loop_name)
        self.item = item

    def __repr__(self) -> str:
        ls_out = []
        ls_out.append("PdExcludeL: ")
        ls_out.append(f"{str(self):}")
        return "\n".join(ls_out)

