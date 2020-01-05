__author__ = 'ikibalin'
__version__ = "2019_12_10"
import os
import numpy

import warnings
from typing import List, Tuple

from cryspy.common.cl_item_constr import ItemConstr
from cryspy.common.cl_loop_constr import LoopConstr
from cryspy.common.cl_fitable import Fitable

class PdProc(ItemConstr):
    """
This section contains the diffraction data set after processing
and application of correction terms. If the data set is
reprocessed, this section may be replaced (with the addition of
a new _pd_block_id entry).

Description in cif file::
 
 _pd_proc_2theta                   4.00
 _pd_proc_2theta_corrected         3.90
 _pd_proc_d_spacing               11.20              
 _pd_proc_intensity_up_net       400.00
 _pd_proc_intensity_down_net     317.00
 _pd_proc_intensity_up_total     460.00
 _pd_proc_intensity_down_total   377.00
 _pd_proc_intensity_bkg_calc      60.00
 _pd_proc_intensity_up           465.80
 _pd_proc_intensity_up_sigma     128.97
 _pd_proc_intensity_down         301.88
 _pd_proc_intensity_down_sigma   129.30

`<https://www.iucr.org/__data/iucr/cifdic_html/1/cif_pd.dic/Cpd_proc.html>`_
    """
    MANDATORY_ATTRIBUTE = ("ttheta", )
    OPTIONAL_ATTRIBUTE = ("ttheta_corrected", "d_spacing",
                          "intensity_up_net", "intensity_down_net", "intensity_up_total", "intensity_down_total", 
                          "intensity_bkg_calc", "intensity_net", "intensity_total", 
                          "intensity_up", "intensity_up_sigma", "intensity_down", "intensity_down_sigma", 
                          "intensity", "intensity_sigma")
    INTERNAL_ATTRIBUTE = ()
    PREFIX = "pd_proc"
    def __init__(self, ttheta=None, ttheta_corrected=None, d_spacing=None, intensity_up_net=None, 
                 intensity_up_total=None,
                 intensity_down_net=None, intensity_down_total=None, intensity_bkg_calc=None, 
                 intensity_net=None, intensity_total=None,
                 intensity_up=None, intensity_up_sigma=None, 
                 intensity_down=None, intensity_down_sigma=None, intensity=None, intensity_sigma=None):
        super(PdProc, self).__init__(mandatory_attribute=self.MANDATORY_ATTRIBUTE, 
                                     optional_attribute=self.OPTIONAL_ATTRIBUTE, 
                                     internal_attribute=self.INTERNAL_ATTRIBUTE,
                                     prefix=self.PREFIX)


        self.ttheta = ttheta
        self.ttheta_corrected = ttheta_corrected
        self.d_spacing = d_spacing
        self.intensity_up_net = intensity_up_net
        self.intensity_down_net = intensity_down_net
        self.intensity_net = intensity_net
        self.intensity_up_total = intensity_up_total
        self.intensity_down_total = intensity_down_total
        self.intensity_total = intensity_total
        self.intensity_bkg_calc = intensity_bkg_calc
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
    def ttheta_corrected(self):
        return getattr(self, "__ttheta_corrected")
    @ttheta_corrected.setter
    def ttheta_corrected(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__ttheta_corrected", x_in)

    @property
    def d_spacing(self):
        return getattr(self, "__d_spacing")
    @d_spacing.setter
    def d_spacing(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__d_spacing", x_in)

    @property
    def intensity_up_net(self):
        return getattr(self, "__intensity_up_net")
    @intensity_up_net.setter
    def intensity_up_net(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__intensity_up_net", x_in)

    @property
    def intensity_up_total(self):
        return getattr(self, "__intensity_up_total")
    @intensity_up_total.setter
    def intensity_up_total(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__intensity_up_total", x_in)

    @property
    def intensity_down_net(self):
        return getattr(self, "__intensity_down_net")
    @intensity_down_net.setter
    def intensity_down_net(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__intensity_down_net", x_in)

    @property
    def intensity_down_total(self):
        return getattr(self, "__intensity_down_total")
    @intensity_down_total.setter
    def intensity_down_total(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__intensity_down_total", x_in)

    @property
    def intensity_net(self):
        return getattr(self, "__intensity_net")
    @intensity_net.setter
    def intensity_net(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__intensity_net", x_in)

    @property
    def intensity_total(self):
        return getattr(self, "__intensity_total")
    @intensity_total.setter
    def intensity_total(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__intensity_total", x_in)

    @property
    def intensity_bkg_calc(self):
        return getattr(self, "__intensity_bkg_calc")
    @intensity_bkg_calc.setter
    def intensity_bkg_calc(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__intensity_bkg_calc", x_in)

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
        ls_out.append(f"PdProc:\n{str(self):}")
        return "\n".join(ls_out)

    def _show_message(self, s_out: str):
        warnings.warn("***  Error ***\n"+s_out, UserWarning, stacklevel=2)

    @property
    def is_variable(self) -> bool:
        return False
    
    def get_variables(self) -> List:
        return []


class PdProcL(LoopConstr):
    """
This section contains the diffraction data set after processing
and application of correction terms. If the data set is
reprocessed, this section may be replaced (with the addition of
a new _pd_block_id entry).

Description in cif file::

 loop_
 _pd_proc_ttheta
 _pd_proc_ttheta_corrected
 _pd_proc_d_spacing
 _pd_proc_intensity_up_net
 _pd_proc_intensity_down_net
 _pd_proc_intensity_up_total
 _pd_proc_intensity_down_total
 _pd_proc_intensity_bkg_calc
 _pd_proc_intensity_up
 _pd_proc_intensity_up_sigma
 _pd_proc_intensity_down
 _pd_proc_intensity_down_sigma
  4.00  3.9  11.2   400.00   317.00   460.00   377.00   60.00    465.80000   128.97000   301.88000   129.30000

`Reference. <https://www.iucr.org/__data/iucr/cifdic_html/1/cif_pd.dic/Cpd_proc.html>`_
    """
    CATEGORY_KEY = ("ttheta", )
    ITEM_CLASS = PdProc
    def __init__(self, item=[], loop_name=""):
        super(PdProcL, self).__init__(category_key=self.CATEGORY_KEY, item_class=self.ITEM_CLASS, loop_name=loop_name)
        self.item = item

    def __repr__(self) -> str:
        ls_out = []
        ls_out.append("PdProcL: ")
        ls_out.append(f"{str(self):}")
        return "\n".join(ls_out)
