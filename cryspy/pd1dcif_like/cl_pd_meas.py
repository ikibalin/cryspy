__author__ = 'ikibalin'
__version__ = "2019_12_10"
import os
import numpy

import warnings
from typing import List, Tuple

from cryspy.common.cl_item_constr import ItemConstr
from cryspy.common.cl_loop_constr import LoopConstr
from cryspy.common.cl_fitable import Fitable

class PdMeas(ItemConstr):
    """
This section contains the measured diffractogram and information
about the conditions used for the measurement of the diffraction 
data set, prior to processing and application of correction
terms. While additional information may be added to the CIF
as data are processed and transported between laboratories
(possibly with the addition of a new _pd_block_id entry), the
information in this section of the CIF will rarely be changed
once data collection is complete.

Where possible, measurements in this section should have no
post-collection processing applied (normalizations, corrections,
smoothing, zero-offset corrections etc.). Such corrected
measurements should be recorded in the _pd_proc_ section.

Data sets that are measured as counts, where a standard
uncertainty can be considered equivalent to the standard
deviation and where the standard deviation can be estimated
as the square root of the number of counts recorded, should
use the _pd_meas_counts_ fields. All other intensity values
should be recorded using _pd_meas_intensity_.

Description in cif file::

 _pd_meas_2theta                 4.00
 _pd_meas_intensity_up         465.80
 _pd_meas_intensity_up_sigma   128.97
 _pd_meas_intensity_down       301.88
 _pd_meas_intensity_down_sigma 129.30

`Reference. <https://www.iucr.org/__data/iucr/cifdic_html/1/cif_pd.dic/Cpd_meas.html>`_
    """
    MANDATORY_ATTRIBUTE = ("ttheta", )
    OPTIONAL_ATTRIBUTE = ("intensity_up", "intensity_up_sigma", "intensity_down", "intensity_down_sigma", 
                          "intensity", "intensity_sigma")
    RELATED_CIF_MANDATORY_ATTRIBUTE = ("2theta", )
    RELATED_CIF_OPTIONAL_ATTRIBUTE = ("intensity_up", "intensity_up_sigma", "intensity_down", "intensity_down_sigma", 
                          "intensity", "intensity_sigma")
    INTERNAL_ATTRIBUTE = ()
    PREFIX = "pd_meas"
    def __init__(self, ttheta=None, intensity_up=None, intensity_up_sigma=None, 
                 intensity_down=None, intensity_down_sigma=None, intensity=None, intensity_sigma=None):
        super(PdMeas, self).__init__(mandatory_attribute=self.MANDATORY_ATTRIBUTE, 
                                     optional_attribute=self.OPTIONAL_ATTRIBUTE, 
                                     internal_attribute=self.INTERNAL_ATTRIBUTE,
                                     prefix=self.PREFIX)

        self.ttheta = ttheta
        self.intensity_up = intensity_up
        self.intensity_up_sigma = intensity_up_sigma
        self.intensity_down = intensity_down
        self.intensity_down_sigma = intensity_down_sigma
        self.intensity = intensity
        self.intensity_sigma = intensity_sigma

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
    def intensity_up_sigma(self):
        return getattr(self, "__intensity_up_sigma")
    @intensity_up_sigma.setter
    def intensity_up_sigma(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__intensity_up_sigma", x_in)

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
    def intensity_down_sigma(self):
        return getattr(self, "__intensity_down_sigma")
    @intensity_down_sigma.setter
    def intensity_down_sigma(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__intensity_down_sigma", x_in)

    @property
    def intensity(self):
        return getattr(self, "__intensity")
    @intensity.setter
    def intensity(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__intensity", x_in)

    @property
    def intensity_sigma(self):
        return getattr(self, "__intensity_sigma")
    @intensity_sigma.setter
    def intensity_sigma(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__intensity_sigma", x_in)

    def __repr__(self):
        ls_out = []
        ls_out.append(f"PdMeas:\n{str(self):}")
        return "\n".join(ls_out)

    def _show_message(self, s_out: str):
        warnings.warn("***  Error ***\n"+s_out, UserWarning, stacklevel=2)

    @property
    def is_variable(self) -> bool:
        return False
    
    def get_variables(self) -> List:
        return []

class PdMeasL(LoopConstr):
    """
This section contains the measured diffractogram and information
about the conditions used for the measurement of the diffraction 
data set, prior to processing and application of correction
terms. While additional information may be added to the CIF
as data are processed and transported between laboratories
(possibly with the addition of a new _pd_block_id entry), the
information in this section of the CIF will rarely be changed
once data collection is complete.

Where possible, measurements in this section should have no
post-collection processing applied (normalizations, corrections,
smoothing, zero-offset corrections etc.). Such corrected
measurements should be recorded in the _pd_proc_ section.

Data sets that are measured as counts, where a standard
uncertainty can be considered equivalent to the standard
deviation and where the standard deviation can be estimated
as the square root of the number of counts recorded, should
use the _pd_meas_counts_ fields. All other intensity values
should be recorded using _pd_meas_intensity_.

Description in cif file::

 loop_
 _pd_meas_ttheta
 _pd_meas_intensity_up
 _pd_meas_intensity_up_sigma
 _pd_meas_intensity_down
 _pd_meas_intensity_down_sigma
  4.00   465.80000   128.97000   301.88000   129.30000
  4.20   323.78000   118.22000   206.06000   120.00000

`Reference. <https://www.iucr.org/__data/iucr/cifdic_html/1/cif_pd.dic/Cpd_meas.html>`_
    """
    CATEGORY_KEY = ("ttheta", )
    ITEM_CLASS = PdMeas
    def __init__(self, item=[], loop_name=""):
        super(PdMeasL, self).__init__(category_key=self.CATEGORY_KEY, item_class=self.ITEM_CLASS, loop_name=loop_name)
        self.item = item

    def __repr__(self) -> str:
        ls_out = []
        ls_out.append("PdMeasL: ")
        ls_out.append(f"{str(self):}")
        return "\n".join(ls_out)

    @property
    def is_polarized(self):
        """
True if polaraized data are defined

and

False for unpolaraized data
        """
        flag = False
        if len(self.item) != 0:
            if ((self.item[0].intensity_up is not None) &
                (self.item[0].intensity_up is not None)):
                flag = True
        return flag