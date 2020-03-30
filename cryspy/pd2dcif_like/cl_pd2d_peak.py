__author__ = 'ikibalin'
__version__ = "2020_01_04"
import os
import numpy

import warnings
from typing import List, Tuple

from cryspy.common.cl_item_constr import ItemConstr
from cryspy.common.cl_loop_constr import LoopConstr
from cryspy.common.cl_fitable import Fitable

class Pd2dPeak(ItemConstr):
    """
This section contains peak information extracted from the
measured or, if present, the processed diffractogram. 

Description in cif file::

 _pd2d_peak_index_h           2
 _pd2d_peak_index_k           2
 _pd2d_peak_index_l           0
 _pd2d_peak_index_mult        4
 _pd2d_peak_2theta           17.2
 _pd2d_peak_f_nucl_sq       100.0
 _pd2d_peak_f_m_p_sin_sq    101.2
 _pd2d_peak_f_m_p_cos_sq     90.0
 _pd2d_peak_cross_sin        87.4
 _pd2d_peak_width_2theta      2.3

    """
    MANDATORY_ATTRIBUTE = ("index_h", "index_k", "index_l")
    OPTIONAL_ATTRIBUTE = ("index_mult", "ttheta",
                          "f_nucl_sq", "f_m_p_sin_sq", "f_m_p_cos_sq", "cross_sin", "width_ttheta")
    RELATED_CIF_MANDATORY_ATTRIBUTE = ("index_h", "index_k", "index_l")
    RELATED_CIF_OPTIONAL_ATTRIBUTE = ("index_mult", "2theta",
                          "f_nucl_sq", "f_m_p_sin_sq", "f_m_p_cos_sq", "cross_sin", "width_2theta")
    INTERNAL_ATTRIBUTE = ()
    PREFIX = "pd2d_peak"
    def __init__(self, index_h=None, index_k=None, index_l=None, index_mult=None, 
                 ttheta=None, f_nucl_sq=None, f_m_p_sin_sq=None, f_m_p_cos_sq=None,
                 cross_sin=None,
                 width_ttheta=None):
        super(Pd2dPeak, self).__init__(mandatory_attribute=self.MANDATORY_ATTRIBUTE, 
                                       optional_attribute=self.OPTIONAL_ATTRIBUTE, 
                                       internal_attribute=self.INTERNAL_ATTRIBUTE,
                                       prefix=self.PREFIX)

        self.index_h = index_h
        self.index_k = index_k
        self.index_l = index_l
        self.index_mult = index_mult
        self.ttheta = ttheta
        self.f_nucl_sq = f_nucl_sq
        self.f_m_p_sin_sq = f_m_p_sin_sq
        self.f_m_p_cos_sq = f_m_p_cos_sq
        self.cross_sin = cross_sin
        self.width_ttheta = width_ttheta

        if self.is_defined:
            self.form_object


    @property
    def index_h(self):
        return getattr(self, "__index_h")
    @index_h.setter
    def index_h(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = int(x)
        setattr(self, "__index_h", x_in)

    @property
    def index_k(self):
        return getattr(self, "__index_k")
    @index_k.setter
    def index_k(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = int(x)
        setattr(self, "__index_k", x_in)

    @property
    def index_l(self):
        return getattr(self, "__index_l")
    @index_l.setter
    def index_l(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = int(x)
        setattr(self, "__index_l", x_in)

    @property
    def index_mult(self):
        return getattr(self, "__index_mult")
    @index_mult.setter
    def index_mult(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = int(x)
        setattr(self, "__index_mult", x_in)

    @property
    def ttheta(self):
        return getattr(self, "__ttheta")
    @ttheta.setter
    def ttheta(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__ttheta", x_in)

    @property
    def f_nucl_sq(self):
        return getattr(self, "__f_nucl_sq")
    @f_nucl_sq.setter
    def f_nucl_sq(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__f_nucl_sq", x_in)

    @property
    def f_m_p_sin_sq(self):
        return getattr(self, "__f_m_p_sin_sq")
    @f_m_p_sin_sq.setter
    def f_m_p_sin_sq(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__f_m_p_sin_sq", x_in)


    @property
    def f_m_p_cos_sq(self):
        return getattr(self, "__f_m_p_cos_sq")
    @f_m_p_cos_sq.setter
    def f_m_p_cos_sq(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__f_m_p_cos_sq", x_in)

    @property
    def cross_sin(self):
        return getattr(self, "__cross_sin")
    @cross_sin.setter
    def cross_sin(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__cross_sin", x_in)


    @property
    def width_ttheta(self):
        return getattr(self, "__width_ttheta")
    @width_ttheta.setter
    def width_ttheta(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__width_ttheta", x_in)

    def __repr__(self):
        ls_out = []
        ls_out.append(f"Pd2dPeak:\n{str(self):}")
        return "\n".join(ls_out)


class Pd2dPeakL(LoopConstr):
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
 _pd2d_peak_index_h
 _pd2d_peak_index_k
 _pd2d_peak_index_l
 _pd2d_peak_index_mult
 _pd2d_peak_ttheta
 _pd2d_peak_f_nucl_sq
 _pd2d_peak_f_m_p_sin_sq
 _pd2d_peak_f_m_p_cos_sq
 _pd2d_peak_cross_sin
 _pd2d_peak_width_ttheta
  2  2  0  4 17.2 100.0 101.2   90.0  87.4 2.3
  2  1  3  2 19.1  78.6 101.0   92.0  82.7 5.7

    """
    CATEGORY_KEY = ("index_h", "index_k", "index_l")
    ITEM_CLASS = Pd2dPeak
    def __init__(self, item=[], loop_name=""):
        super(Pd2dPeakL, self).__init__(category_key=self.CATEGORY_KEY, item_class=self.ITEM_CLASS, loop_name=loop_name)
        self.item = item

    def __repr__(self) -> str:
        ls_out = []
        ls_out.append("Pd2dPeakL: ")
        ls_out.append(f"{str(self):}")
        return "\n".join(ls_out)