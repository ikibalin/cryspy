__author__ = 'ikibalin'
__version__ = "2019_12_11"
import os
import numpy

import warnings
from typing import List, Tuple

from cryspy.common.cl_item_constr import ItemConstr
from cryspy.common.cl_loop_constr import LoopConstr
from cryspy.common.cl_fitable import Fitable

class Range(ItemConstr):
    """
Describe the range.

Description in cif file::

 _range_ttheta_min     4.000
 _range_ttheta_max    80.000
    """
    MANDATORY_ATTRIBUTE = ("ttheta_min", "ttheta_max")
    OPTIONAL_ATTRIBUTE = ("phi_min", "phi_max")
    INTERNAL_ATTRIBUTE = ()
    PREFIX = "range"
    def __init__(self, ttheta_min=None, ttheta_max=None, phi_min=None, phi_max=None):
        super(Range, self).__init__(mandatory_attribute=self.MANDATORY_ATTRIBUTE, 
                                    optional_attribute=self.OPTIONAL_ATTRIBUTE, 
                                    internal_attribute=self.INTERNAL_ATTRIBUTE,
                                    prefix=self.PREFIX)

        self.ttheta_min = ttheta_min
        self.ttheta_max = ttheta_max
        self.phi_min = phi_min
        self.phi_max = phi_max

        if self.is_defined:
            self.form_object

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

    @property
    def phi_min(self):
        return getattr(self, "__phi_min")
    @phi_min.setter
    def phi_min(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__phi_min", x_in)

    @property
    def phi_max(self):
        return getattr(self, "__phi_max")
    @phi_max.setter
    def phi_max(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__phi_max", x_in)

    def __repr__(self):
        ls_out = []
        ls_out.append(f"Range:\n{str(self):}")
        return "\n".join(ls_out)

