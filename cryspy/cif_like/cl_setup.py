__author__ = 'ikibalin'
__version__ = "2019_12_11"
import os
import numpy

import warnings
from typing import List, Tuple

from cryspy.common.cl_item_constr import ItemConstr
from cryspy.common.cl_loop_constr import LoopConstr
from cryspy.common.cl_fitable import Fitable

class Setup(ItemConstr):
    """
Describe the setup.

Description in cif file::

 _setup_wavelength   0.84
 _setup_field        1.00
 _setup_offset_2theta -0.385
 _setup_offset_phi -0.385
    """
    MANDATORY_ATTRIBUTE = ("wavelength", )
    OPTIONAL_ATTRIBUTE = ("field", "offset_ttheta", "offset_phi")
    RELATED_CIF_MANDATORY_ATTRIBUTE = ("wavelength", )
    RELATED_CIF_OPTIONAL_ATTRIBUTE = ("field", "offset_2theta", "offset_phi")
    INTERNAL_ATTRIBUTE = ()
    PREFIX = "setup"
    def __init__(self, wavelength=None, field=None, offset_ttheta=None, offset_phi=None):
        super(Setup, self).__init__(mandatory_attribute=self.MANDATORY_ATTRIBUTE, 
                                    optional_attribute=self.OPTIONAL_ATTRIBUTE, 
                                    internal_attribute=self.INTERNAL_ATTRIBUTE,
                                    prefix=self.PREFIX)

        self.wavelength = wavelength
        self.field = field
        self.offset_ttheta = offset_ttheta
        self.offset_phi = offset_phi

        if self.is_defined:
            self.form_object

    @property
    def wavelength(self):
        return getattr(self, "__wavelength")
    @wavelength.setter
    def wavelength(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__wavelength", x_in)

    @property
    def field(self):
        return getattr(self, "__field")
    @field.setter
    def field(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__field", x_in)

    @property
    def offset_ttheta(self):
        return getattr(self, "__offset_ttheta")
    @offset_ttheta.setter
    def offset_ttheta(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__offset_ttheta", x_in)

    @property
    def offset_phi(self):
        return getattr(self, "__offset_phi")
    @offset_phi.setter
    def offset_phi(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__offset_phi", x_in)


    def __repr__(self):
        ls_out = []
        ls_out.append(f"Setup:\n{str(self):}")
        return "\n".join(ls_out)

    def _show_message(self, s_out: str):
        warnings.warn("***  Error ***\n"+s_out, UserWarning, stacklevel=2)
        
    @property
    def is_variable(self) -> bool:
        """
Output: True if there is any refined parameter
        """
        _l = [self.wavelength.refinement,
              self.offset_ttheta.refinement]
        if self.offset_phi is not None: _l.append(self.offset_phi.refinement)

        res = any(_l)
        return res

    def get_variables(self) -> List:
        """
Output: the list of the refined parameters
        """
        l_variable = []
        if self.wavelength.refinement: l_variable.append(self.wavelength)
        if self.offset_ttheta is not None:
            if self.offset_ttheta.refinement: l_variable.append(self.offset_ttheta)
        if self.offset_phi is not None:
            if self.offset_phi.refinement: l_variable.append(self.offset_phi)
        return l_variable
