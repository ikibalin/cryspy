__author__ = 'ikibalin'
__version__ = "2019_12_10"
import os
import numpy

import warnings
from typing import List, Tuple

from cryspy.common.cl_item_constr import ItemConstr
from cryspy.common.cl_loop_constr import LoopConstr
from cryspy.common.cl_fitable import Fitable

class DiffrnRadiation(ItemConstr):
    """
Describe the beam polarization.

Description in cif file::

 _diffrn_radiation_efficiency    1.0
 _diffrn_radiation_polarization -0.87
    """
    MANDATORY_ATTRIBUTE = ("polarization", "efficiency")
    OPTIONAL_ATTRIBUTE = ()
    INTERNAL_ATTRIBUTE = ()
    PREFIX = "diffrn_radiation"
    def __init__(self, polarization=None, efficiency=None):
        super(DiffrnRadiation, self).__init__(mandatory_attribute=self.MANDATORY_ATTRIBUTE, 
                                              optional_attribute=self.OPTIONAL_ATTRIBUTE, 
                                              internal_attribute=self.INTERNAL_ATTRIBUTE,
                                              prefix=self.PREFIX)

        self.polarization = polarization
        self.efficiency = efficiency

        if self.is_defined:
            self.form_object
        
    @property
    def polarization(self):
        """
The polarization of the incident beam. 

The permitted range is -1.0 -> 1.0
        """
        return getattr(self, "__polarization")
    @polarization.setter
    def polarization(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__polarization", x_in)

    @property
    def efficiency(self):
        """
The efficiency of the efficiency. 

The permitted range is -1.0 -> 1.0
        """
        return getattr(self, "__efficiency")
    @efficiency.setter
    def efficiency(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__efficiency", x_in)


    def __repr__(self):
        ls_out = []
        ls_out.append(f"DiffrnRadiation:\n{str(self):}")
        return "\n".join(ls_out)

    def _show_message(self, s_out: str):
        warnings.warn("***  Error ***\n"+s_out, UserWarning, stacklevel=2)

    @property
    def is_variable(self) -> bool:
        """
Output: True if there is any refined parameter
        """
        res = any([self.polarization.refinement,
                   self.efficiency.refinement])
        return res

    def get_variables(self) -> List:
        """
Output: the list of the refined parameters
        """
        l_variable = []
        if self.polarization.refinement: l_variable.append(self.polarization)
        if self.efficiency.refinement: l_variable.append(self.efficiency)
        return l_variable


