__author__ = 'ikibalin'
__version__ = "2019_12_11"
import os
import numpy

import warnings
from typing import List, Tuple

from cryspy.common.cl_item_constr import ItemConstr
from cryspy.common.cl_loop_constr import LoopConstr
from cryspy.common.cl_fitable import Fitable

class PdPeak(ItemConstr):
    """
This section contains peak information extracted from the
measured or, if present, the processed diffractogram. Each
peak in this table will have a unique label (see _pd_peak_id).
The reflections and phases associated with each peak will be
specified in other sections (see the _pd_refln_ and
_pd_phase_ sections).

Note that peak positions are customarily determined from the
processed diffractogram and thus corrections for position
and intensity will have been previously applied.

Description in cif file::

 _pd_peak_index_h            2  
 _pd_peak_index_k            2  
 _pd_peak_index_l            0  
 _pd_peak_index_mult         4  
 _pd_peak_ttheta            17.2
 _pd_peak_intensity_up     100.0
 _pd_peak_intensity_down    90.0
 _pd_peak_width_ttheta       2.3

`Reference. <https://www.iucr.org/__data/iucr/cifdic_html/1/cif_pd.dic/Cpd_peak.html>`_
    """
    MANDATORY_ATTRIBUTE = ("index_h", "index_k", "index_l")
    OPTIONAL_ATTRIBUTE = ("index_mult", "ttheta",
                          "intensity_up", "intensity_down", "width_ttheta")
    INTERNAL_ATTRIBUTE = ()
    PREFIX = "pd_peak"
    def __init__(self, index_h=None, index_k=None, index_l=None, index_mult=None, 
                 ttheta=None, intensity_up=None, intensity_down=None,
                 width_ttheta=None):
        super(PdPeak, self).__init__(mandatory_attribute=self.MANDATORY_ATTRIBUTE, 
                                     optional_attribute=self.OPTIONAL_ATTRIBUTE, 
                                     internal_attribute=self.INTERNAL_ATTRIBUTE,
                                     prefix=self.PREFIX)

        self.index_h = index_h
        self.index_k = index_k
        self.index_l = index_l
        self.index_mult = index_mult
        self.ttheta = ttheta
        self.intensity_up = intensity_up
        self.intensity_down = intensity_down
        self.width_ttheta = width_ttheta

    @property
    def index_h(self):
        return getattr(self, "__index_h")
    @index_h.setter
    def index_h(self, x):
        if x is None:
            x_in = None
        else:
            x_in = int(x)
        setattr(self, "__index_h", x_in)

    @property
    def index_k(self):
        return getattr(self, "__index_k")
    @index_k.setter
    def index_k(self, x):
        if x is None:
            x_in = None
        else:
            x_in = int(x)
        setattr(self, "__index_k", x_in)

    @property
    def index_l(self):
        return getattr(self, "__index_l")
    @index_l.setter
    def index_l(self, x):
        if x is None:
            x_in = None
        else:
            x_in = int(x)
        setattr(self, "__index_l", x_in)

    @property
    def index_mult(self):
        return getattr(self, "__index_mult")
    @index_mult.setter
    def index_mult(self, x):
        if x is None:
            x_in = None
        else:
            x_in = int(x)
        setattr(self, "__index_mult", x_in)

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
    def intensity_up(self):
        return getattr(self, "__intensity_up")
    @intensity_up.setter
    def intensity_up(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__intensity_up", x_in)

    @property
    def intensity_down(self):
        return getattr(self, "__intensity_down")
    @intensity_down.setter
    def intensity_down(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__intensity_down", x_in)

    @property
    def width_ttheta(self):
        return getattr(self, "__width_ttheta")
    @width_ttheta.setter
    def width_ttheta(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__width_ttheta", x_in)

    def __repr__(self):
        ls_out = []
        ls_out.append(f"PdPeak:\n{str(self):}")
        return "\n".join(ls_out)

    def _show_message(self, s_out: str):
        warnings.warn("***  Error ***\n"+s_out, UserWarning, stacklevel=2)

    @property
    def is_variable(self) -> bool:
        return False
    
    def get_variables(self) -> List:
        return []





class PdPeakL(LoopConstr):
    """
This section contains peak information extracted from the
measured or, if present, the processed diffractogram. Each
peak in this table will have a unique label (see _pd_peak_id).
The reflections and phases associated with each peak will be
specified in other sections (see the _pd_refln_ and
_pd_phase_ sections).

Note that peak positions are customarily determined from the
processed diffractogram and thus corrections for position
and intensity will have been previously applied.

Description in cif file::

 loop_
 _pd_peak_index_h
 _pd_peak_index_k
 _pd_peak_index_l
 _pd_peak_index_mult
 _pd_peak_2theta
 _pd_peak_intensity_up
 _pd_peak_intensity_down
 _pd_peak_width_2theta
  2  2  0  4 17.2 100.0  90.0  2.3

`Reference. <https://www.iucr.org/__data/iucr/cifdic_html/1/cif_pd.dic/Cpd_peak.html>`_
    """
    CATEGORY_KEY = ("index_h", "index_k", "index_l")
    ITEM_CLASS = PdPeak
    def __init__(self, item=[], loop_name=""):
        super(PdPeakL, self).__init__(category_key=self.CATEGORY_KEY, item_class=self.ITEM_CLASS, loop_name=loop_name)
        self.item = item

    def __repr__(self) -> str:
        ls_out = []
        ls_out.append("PdPeakL: ")
        ls_out.append(f"{str(self):}")
        return "\n".join(ls_out)
