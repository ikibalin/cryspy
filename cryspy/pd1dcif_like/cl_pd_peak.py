__author__ = 'ikibalin'
__version__ = "2020_01_04"
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
 _pd_peak_2theta            17.2
 _pd_peak_intensity_up     100.0
 _pd_peak_intensity_down    90.0
 _pd_peak_width_2theta       2.3

`<https://www.iucr.org/__data/iucr/cifdic_html/1/cif_pd.dic/Cpd_peak.html>`_
    """
    MANDATORY_ATTRIBUTE = ("index_h", "index_k", "index_l")
    OPTIONAL_ATTRIBUTE = ("index_mult", "ttheta",
                          "intensity_up", "intensity_down", "width_ttheta")
    RELATED_CIF_MANDATORY_ATTRIBUTE = ("index_h", "index_k", "index_l")
    RELATED_CIF_OPTIONAL_ATTRIBUTE = ("index_mult", "2theta",
                          "intensity_up", "intensity_down", "width_2theta")
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
    def intensity_up(self):
        return getattr(self, "__intensity_up")
    @intensity_up.setter
    def intensity_up(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__intensity_up", x_in)

    @property
    def intensity_down(self):
        return getattr(self, "__intensity_down")
    @intensity_down.setter
    def intensity_down(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__intensity_down", x_in)

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
        ls_out.append(f"PdPeak:\n{str(self):}")
        return "\n".join(ls_out)



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
    INTERNAL_ATTRIBUTE = ("numpy_index_h", "numpy_index_k", 
                          "numpy_index_l", "numpy_index_mult", 
                          "numpy_ttheta", "numpy_intensity_up", 
                          "numpy_intensity_down", 
                          "numpy_width_ttheta")
    def __init__(self, item=[], loop_name=""):
        super(PdPeakL, self).__init__(category_key=self.CATEGORY_KEY, item_class=self.ITEM_CLASS, 
                                      loop_name=loop_name, internal_attribute=self.INTERNAL_ATTRIBUTE)
        self.item = item

    def __repr__(self) -> str:
        ls_out = []
        ls_out.append("PdPeakL: ")
        ls_out.append(f"{str(self):}")
        return "\n".join(ls_out)

    def get_numpy_index_h(self):
        return getattr(self, "__numpy_index_h")

    def set_numpy_index_h(self, x):
        setattr(self, "__numpy_index_h", x)

    def get_numpy_index_k(self):
        return getattr(self, "__numpy_index_k")

    def set_numpy_index_k(self, x):
        setattr(self, "__numpy_index_k", x)

    def get_numpy_index_l(self):
        return getattr(self, "__numpy_index_l")

    def set_numpy_index_l(self, x):
        setattr(self, "__numpy_index_l", x)

    def get_numpy_index_mult(self):
        return getattr(self, "__numpy_index_mult")

    def set_numpy_index_mult(self, x):
        setattr(self, "__numpy_index_mult", x)

    def get_numpy_ttheta(self):
        return getattr(self, "__numpy_ttheta")

    def set_numpy_ttheta(self, x):
        setattr(self, "__numpy_ttheta", x)

    def get_numpy_intensity_up(self):
        return getattr(self, "__numpy_intensity_up")

    def set_numpy_intensity_up(self, x):
        setattr(self, "__numpy_intensity_up", x)

    def get_numpy_intensity_down(self):
        return getattr(self, "__numpy_intensity_down")

    def set_numpy_intensity_down(self, x):
        setattr(self, "__numpy_intensity_down", x)

    def get_numpy_width_ttheta(self):
        return getattr(self, "__numpy_width_ttheta")

    def set_numpy_width_ttheta(self, x):
        setattr(self, "__numpy_width_ttheta", x)

    def transform_items_to_numpy_arrays(self):
        """
Transform items to numpy arrays (to speed up the calculations):

    numpy_index_h: 1D numpy array of ttheta, dtype=int
    numpy_index_k: 1D numpy array of ttheta, dtype=int
    numpy_index_l: 1D numpy array of ttheta, dtype=int
    numpy_index_mult: 1D numpy array of ttheta, dtype=int
    numpy_ttheta: 1D numpy array of ttheta, dtype=float
    numpy_intensity_up: 1D numpy array of ttheta, dtype=float
    numpy_intensity_down: 1D numpy array of ttheta, dtype=float
    numpy_width_ttheta: 1D numpy array of ttheta, dtype=float
        """

        setattr(self, "__numpy_index_h", numpy.array(self.index_h, dtype=int))
        setattr(self, "__numpy_index_k", numpy.array(self.index_k, dtype=int))
        setattr(self, "__numpy_index_l", numpy.array(self.index_l, dtype=int))
        setattr(self, "__numpy_index_mult", numpy.array(self.index_mult, dtype=int))
        setattr(self, "__numpy_ttheta", numpy.array(self.ttheta, dtype=float))
        setattr(self, "__numpy_intensity_up", numpy.array(self.intensity_up, dtype=float))
        setattr(self, "__numpy_intensity_down", numpy.array(self.intensity_down, dtype=float))
        setattr(self, "__numpy_intensity_widht_ttheta", numpy.array(self.widht_ttheta, dtype=float))


    def transform_numpy_arrays_to_items(self):
        """
Transform data from numpy arrays to items:

    numpy_index_h: 1D numpy array of ttheta, dtype=int
    numpy_index_k: 1D numpy array of ttheta, dtype=int
    numpy_index_l: 1D numpy array of ttheta, dtype=int
    numpy_index_mult: 1D numpy array of ttheta, dtype=int
    numpy_ttheta: 1D numpy array of ttheta, dtype=float
    numpy_intensity_up: 1D numpy array of ttheta, dtype=float
    numpy_intensity_down: 1D numpy array of ttheta, dtype=float
    numpy_width_ttheta: 1D numpy array of ttheta, dtype=float
        """
        numpy_index_h = getattr(self, "__numpy_index_h")
        if numpy_index_h is None: return
        l_item = [PdPeak(index_h=_val) for _val in numpy_index_h]

        np_val = getattr(self, "__numpy_index_k")
        if np_val is not None: 
            for _item, val in zip(l_item, np_val):
                _item.index_k = val

        np_val = getattr(self, "__numpy_index_l")
        if np_val is not None: 
            for _item, val in zip(l_item, np_val):
                _item.index_l = val

        np_val = getattr(self, "__numpy_index_mult")
        if np_val is not None: 
            for _item, val in zip(l_item, np_val):
                _item.index_mult = val

        np_val = getattr(self, "__numpy_ttheta")
        if np_val is not None: 
            for _item, val in zip(l_item, np_val):
                _item.ttheta = val

        np_val = getattr(self, "__numpy_intensity_up")
        if np_val is not None: 
            for _item, val in zip(l_item, np_val):
                _item.intensity_up = val

        np_val = getattr(self, "__numpy_intensity_down")
        if np_val is not None: 
            for _item, val in zip(l_item, np_val):
                _item.intensity_down = val

        np_val = getattr(self, "__numpy_width_ttheta")
        if np_val is not None: 
            for _item, val in zip(l_item, np_val):
                _item.width_ttheta = val

        self.item = l_item
