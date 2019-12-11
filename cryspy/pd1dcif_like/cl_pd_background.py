__author__ = 'ikibalin'
__version__ = "2019_12_10"
import os
import numpy

import warnings
from typing import List, Tuple

from cryspy.common.cl_item_constr import ItemConstr
from cryspy.common.cl_loop_constr import LoopConstr
from cryspy.common.cl_fitable import Fitable

class PdBackground(ItemConstr):
    """
Description in cif file::

 _pd_background_ttheta       4.5
 _pd_background_intensity  256.0
    """
    MANDATORY_ATTRIBUTE = ("ttheta", "intensity")
    OPTIONAL_ATTRIBUTE = ()
    INTERNAL_ATTRIBUTE = ()
    PREFIX = "pd_background"
    def __init__(self, ttheta=None, intensity=None):
        super(PdBackground, self).__init__(mandatory_attribute=self.MANDATORY_ATTRIBUTE, 
                                           optional_attribute=self.OPTIONAL_ATTRIBUTE, 
                                           internal_attribute=self.INTERNAL_ATTRIBUTE,
                                           prefix=self.PREFIX)

        self.ttheta = ttheta
        self.intensity = intensity

        if self.is_defined:
            self.form_object

    @property
    def ttheta(self):
        return getattr(self, "__ttheta")
    @ttheta.setter
    def ttheta(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__ttheta", x_in)

    @property
    def intensity(self):
        return getattr(self, "__intensity")
    @intensity.setter
    def intensity(self, x):
        if x is None:
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__intensity", x_in)

    def __repr__(self):
        ls_out = []
        ls_out.append(f"PdBackground:\n{str(self):}")
        return "\n".join(ls_out)

    def _show_message(self, s_out: str):
        warnings.warn("***  Error ***\n"+s_out, UserWarning, stacklevel=2)

    @property
    def is_variable(self) -> bool:
        """
Output: True if there is any refined parameter
        """
        res = any([self.scale.refinement])
        return res
    
    def get_variables(self) -> List:
        """
Output: the list of the refined parameters
        """
        l_variable = []
        if self.intensity.refinement: l_variable.append(self.intensity)
        return l_variable

    def interpolate_by_points(self, tth):
        l_ttheta = self.ttheta
        l_intensity = self.intensity
        tth_b = numpy.array(l_ttheta, dtype=float)
        int_b = numpy.array(l_intensity, dtype=float)
        if len(l_ttheta) == 0:
            int_1d = numpy.zeros(tth.size, dtype=float)
        else:
            int_1d = numpy.interp(tth, tth_b, int_b)
        return int_1d

class PdBackgroundL(LoopConstr):
    """
Description in cif file::

 loop_
 _pd_background_ttheta
 _pd_background_intensity
  4.5  256.0
  40.0  158.0
  80.0  65.0
    """
    CATEGORY_KEY = ("ttheta", )
    ITEM_CLASS = PdBackground
    def __init__(self, item=[], loop_name=""):
        super(PdBackgroundL, self).__init__(category_key=self.CATEGORY_KEY, item_class=self.ITEM_CLASS, loop_name=loop_name)
        self.item = item

    def __repr__(self) -> str:
        ls_out = []
        ls_out.append("PdBackgroundL: ")
        ls_out.append(f"{str(self):}")
        return "\n".join(ls_out)
