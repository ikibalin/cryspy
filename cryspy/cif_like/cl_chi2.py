__author__ = 'ikibalin'
__version__ = "2019_12_11"
import os
import numpy

import warnings
from typing import List, Tuple

from cryspy.common.cl_item_constr import ItemConstr
from cryspy.common.cl_loop_constr import LoopConstr
from cryspy.common.cl_fitable import Fitable

class Chi2(ItemConstr):
    """
Describe the chi2.

Description in cif file::

 _chi2_sum  True
 _chi2_diff False
 _chi2_up   False
 _chi2_down False
    """
    MANDATORY_ATTRIBUTE = ("sum", "diff", "up", "down")
    OPTIONAL_ATTRIBUTE = ()
    INTERNAL_ATTRIBUTE = ()
    PREFIX = "chi2"
    def __init__(self, sum=None, diff=None, up=None, down=None):
        super(Chi2, self).__init__(mandatory_attribute=self.MANDATORY_ATTRIBUTE, 
                                   optional_attribute=self.OPTIONAL_ATTRIBUTE, 
                                   internal_attribute=self.INTERNAL_ATTRIBUTE,
                                   prefix=self.PREFIX)

        self.sum = sum
        self.diff = diff
        self.up = up
        self.down = down

        if self.is_defined:
            self.form_object

    @property
    def sum(self):
        return getattr(self, "__sum")
    @sum.setter
    def sum(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = (str(x).strip().lower() == "true")
        setattr(self, "__sum", x_in)

    @property
    def diff(self):
        return getattr(self, "__diff")
    @diff.setter
    def diff(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = (str(x).strip().lower() == "true")
        setattr(self, "__diff", x_in)

    @property
    def up(self):
        return getattr(self, "__up")
    @up.setter
    def up(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = (str(x).strip().lower() == "true")
        setattr(self, "__up", x_in)

    @property
    def down(self):
        return getattr(self, "__down")
    @down.setter
    def down(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = (str(x).strip().lower() == "true")
        setattr(self, "__down", x_in)

    def __repr__(self):
        ls_out = []
        ls_out.append(f"Chi2:\n{str(self):}")
        return "\n".join(ls_out)

    def _show_message(self, s_out: str):
        warnings.warn("***  Error ***\n"+s_out, UserWarning, stacklevel=2)
        
    @property
    def is_variable(self) -> bool:
        """
Output: True if there is any refined parameter
        """
        return False

    def get_variables(self) -> List:
        """
Output: the list of the refined parameters
        """
        l_variable = []
        return l_variable


