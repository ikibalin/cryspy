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
                          "intensity_bkg_calc", "intensity_net", "intensity_total", "intensity_diff_total", 
                          "intensity_up", "intensity_up_sigma", "intensity_down", "intensity_down_sigma", 
                          "intensity", "intensity_sigma")
    RELATED_CIF_MANDATORY_ATTRIBUTE = ("2theta", )
    RELATED_CIF_OPTIONAL_ATTRIBUTE = ("2theta_corrected", "d_spacing",
                          "intensity_up_net", "intensity_down_net", "intensity_up_total", "intensity_down_total", 
                          "intensity_bkg_calc", "intensity_net", "intensity_total", "intensity_diff_total", 
                          "intensity_up", "intensity_up_sigma", "intensity_down", "intensity_down_sigma", 
                          "intensity", "intensity_sigma")
    INTERNAL_ATTRIBUTE = ()
    PREFIX = "pd_proc"
    def __init__(self, ttheta=None, ttheta_corrected=None, d_spacing=None, intensity_up_net=None, 
                 intensity_up_total=None,
                 intensity_down_net=None, intensity_down_total=None, intensity_bkg_calc=None, 
                 intensity_net=None, intensity_total=None, intensity_diff_total=None,
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
        self.intensity_diff_total = intensity_diff_total
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
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__ttheta", x_in)

    @property
    def ttheta_corrected(self):
        return getattr(self, "__ttheta_corrected")
    @ttheta_corrected.setter
    def ttheta_corrected(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__ttheta_corrected", x_in)

    @property
    def d_spacing(self):
        return getattr(self, "__d_spacing")
    @d_spacing.setter
    def d_spacing(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__d_spacing", x_in)

    @property
    def intensity_up_net(self):
        return getattr(self, "__intensity_up_net")
    @intensity_up_net.setter
    def intensity_up_net(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__intensity_up_net", x_in)

    @property
    def intensity_up_total(self):
        return getattr(self, "__intensity_up_total")
    @intensity_up_total.setter
    def intensity_up_total(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__intensity_up_total", x_in)

    @property
    def intensity_down_net(self):
        return getattr(self, "__intensity_down_net")
    @intensity_down_net.setter
    def intensity_down_net(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__intensity_down_net", x_in)

    @property
    def intensity_down_total(self):
        return getattr(self, "__intensity_down_total")
    @intensity_down_total.setter
    def intensity_down_total(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__intensity_down_total", x_in)

    @property
    def intensity_net(self):
        return getattr(self, "__intensity_net")
    @intensity_net.setter
    def intensity_net(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__intensity_net", x_in)

    @property
    def intensity_total(self):
        return getattr(self, "__intensity_total")
    @intensity_total.setter
    def intensity_total(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__intensity_total", x_in)

    @property
    def intensity_diff_total(self):
        return getattr(self, "__intensity_diff_total")
    @intensity_diff_total.setter
    def intensity_diff_total(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__intensity_diff_total", x_in)

    @property
    def intensity_bkg_calc(self):
        return getattr(self, "__intensity_bkg_calc")
    @intensity_bkg_calc.setter
    def intensity_bkg_calc(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__intensity_bkg_calc", x_in)

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
    def intensity_up_sigma(self):
        return getattr(self, "__intensity_up_sigma")
    @intensity_up_sigma.setter
    def intensity_up_sigma(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__intensity_up_sigma", x_in)

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
    def intensity_down_sigma(self):
        return getattr(self, "__intensity_down_sigma")
    @intensity_down_sigma.setter
    def intensity_down_sigma(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__intensity_down_sigma", x_in)

    @property
    def intensity(self):
        return getattr(self, "__intensity")
    @intensity.setter
    def intensity(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__intensity", x_in)

    @property
    def intensity_sigma(self):
        return getattr(self, "__intensity_sigma")
    @intensity_sigma.setter
    def intensity_sigma(self, x):
        if ((x is None) | (x == ".")):
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
    INTERNAL_ATTRIBUTE = ("numpy_ttheta", "numpy_ttheta_corrected", 
                          "numpy_intensity_up_net", "numpy_intensity_down_net", 
                          "numpy_intensity_up_total", "numpy_intensity_down_total", 
                          "numpy_intensity_bkg_calc", 
                          "numpy_intensity_net", "numpy_intensity_total", "numpy_intensity_diff_total", 
                          "numpy_intensity_up", "numpy_intensity_up_sigma", 
                          "numpy_intensity_down", "numpy_intensity_down_sigma", 
                          "numpy_intensity", "numpy_intensity_sigma")
    def __init__(self, item=[], loop_name=""):
        super(PdProcL, self).__init__(category_key=self.CATEGORY_KEY, item_class=self.ITEM_CLASS, loop_name=loop_name,
                                      internal_attribute=self.INTERNAL_ATTRIBUTE)
        self.item = item

    def __repr__(self) -> str:
        ls_out = []
        ls_out.append("PdProcL: ")
        ls_out.append(f"{str(self):}")
        return "\n".join(ls_out)

    def get_numpy_ttheta(self):
        return getattr(self, "__numpy_ttheta")

    def set_numpy_ttheta(self, x):
        setattr(self, "__numpy_ttheta", x)

    def get_numpy_ttheta_corrected(self):
        return getattr(self, "__numpy_ttheta_corrected")

    def set_numpy_ttheta_corrected(self, x):
        setattr(self, "__numpy_ttheta_corrected", x)

    def get_numpy_intensity_up_net(self):
        return getattr(self, "__numpy_intensity_up_net")

    def set_numpy_intensity_up_net(self, x):
        setattr(self, "__numpy_intensity_up_net", x)

    def get_numpy_intensity_down_net(self):
        return getattr(self, "__numpy_intensity_down_net")

    def set_numpy_intensity_down_net(self, x):
        setattr(self, "__numpy_intensity_down_net", x)

    def get_numpy_intensity_up_total(self):
        return getattr(self, "__numpy_intensity_up_total")

    def set_numpy_intensity_up_total(self, x):
        setattr(self, "__numpy_intensity_up_total", x)

    def get_numpy_intensity_down_total(self):
        return getattr(self, "__numpy_intensity_down_total")

    def set_numpy_intensity_down_total(self, x):
        setattr(self, "__numpy_intensity_down_total", x)

    def get_numpy_intensity_bkg_calc(self):
        return getattr(self, "__numpy_intensity_bkg_calc")

    def set_numpy_intensity_bkg_calc(self, x):
        setattr(self, "__numpy_intensity_bkg_calc", x)

    def get_numpy_intensity_net(self):
        return getattr(self, "__numpy_intensity_net")

    def set_numpy_intensity_net(self, x):
        setattr(self, "__numpy_intensity_net", x)

    def get_numpy_intensity_total(self):
        return getattr(self, "__numpy_intensity_total")

    def set_numpy_intensity_total(self, x):
        setattr(self, "__numpy_intensity_total", x)

    def get_numpy_intensity_diff_total(self):
        return getattr(self, "__numpy_intensity_diff_total")

    def set_numpy_intensity_diff_total(self, x):
        setattr(self, "__numpy_intensity_diff_total", x)

    def get_numpy_intensity_up(self):
        return getattr(self, "__numpy_intensity_up")

    def set_numpy_intensity_up(self, x):
        setattr(self, "__numpy_intensity_up", x)

    def get_numpy_intensity_up_sigma(self):
        return getattr(self, "__numpy_intensity_up_sigma")

    def set_numpy_intensity_up_sigma(self, x):
        setattr(self, "__numpy_intensity_up_sigma", x)

    def get_numpy_intensity_down(self):
        return getattr(self, "__numpy_intensity_down")

    def set_numpy_intensity_down(self, x):
        setattr(self, "__numpy_intensity_down", x)

    def get_numpy_intensity_down_sigma(self):
        return getattr(self, "__numpy_intensity_down_sigma")

    def set_numpy_intensity_down_sigma(self, x):
        setattr(self, "__numpy_intensity_down_sigma", x)

    def get_numpy_intensity(self):
        return getattr(self, "__numpy_intensity")

    def set_numpy_intensity(self, x):
        setattr(self, "__numpy_intensity", x)

    def get_numpy_intensity_sigma(self):
        return getattr(self, "__numpy_intensity_sigma")

    def set_numpy_intensity_sigma(self, x):
        setattr(self, "__numpy_intensity_sigma", x)

    def transform_items_to_numpy_arrays(self):
        """
Transform items to numpy arrays (to speed up the calculations):

    numpy_ttheta: 1D numpy array of ttheta, dtype=float
    numpy_ttheta_corrected: 1D numpy array of ttheta_corrected, dtype=float 
    numpy_intensity_up_net: 1D numpy array of intensity_up_net, dtype=float
    numpy_intensity_down_net: 1D numpy array of intensity_down_net, dtype=float
    numpy_intensity_up_total: 1D numpy array of intensity_up_total, dtype=float
    numpy_intensity_down_total 1D numpy array of intensity_down_total, dtype=floatl: 1D numpy array of , dtype=float
    numpy_intensity_total 1D numpy array of intensity_total, dtype=floatl: 1D numpy array of , dtype=float
    numpy_intensity_diff_total 1D numpy array of intensity_diff_total, dtype=floatl: 1D numpy array of , dtype=float
    numpy_intensity_bkg_calc 1D numpy array of intensity_bkg_calc, dtype=float_calc: 1D numpy array of , dtype=float
    numpy_intensity_up: 1D numpy array of intensity_up, dtype=float
    numpy_intensity_up_sigma: 1D numpy array of intensity_up_sigma, dtype=float
    numpy_intensity_down: 1D numpy array of intensity_down, dtype=float
    numpy_intensity_down_sigma: 1D numpy array of intensity_down_sigma, dtype=float
    numpy_intensity: 1D numpy array of intensity, dtype=float
    numpy_intensity_sigma: 1D numpy array of intensity_sigma, dtype=float
        """

        numpy_ttheta = numpy.array(self.ttheta, dtype=float)
        setattr(self, "__numpy_ttheta", numpy_ttheta)
        numpy_ttheta_corrected = numpy.array(self.ttheta_corrected, dtype=float)
        setattr(self, "__numpy_ttheta_corrected", numpy_ttheta_corrected)
        numpy_intensity_up_net = numpy.array(self.intensity_up_net, dtype=float)
        setattr(self, "__numpy_intensity_up_net", numpy_intensity_up_net)
        numpy_intensity_down_net = numpy.array(self.intensity_down_net, dtype=float)
        setattr(self, "__numpy_intensity_down_net", numpy_intensity_down_net)
        numpy_intensity_up_total = numpy.array(self.intensity_up_total, dtype=float)
        setattr(self, "__numpy_intensity_up_total", numpy_intensity_up_total)
        numpy_intensity_down_total = numpy.array(self.intensity_down_total, dtype=float)
        setattr(self, "__numpy_intensity_down_total", numpy_intensity_down_total)
        numpy_intensity_bkg_calc = numpy.array(self.intensity_bkg_calc, dtype=float)
        setattr(self, "__numpy_intensity_bkg_calc", numpy_intensity_bkg_calc)
        numpy_intensity_net = numpy.array(self.intensity_net, dtype=float)
        setattr(self, "__numpy_intensity_net", numpy_intensity_net)
        numpy_intensity_total = numpy.array(self.intensity_total, dtype=float)
        setattr(self, "__numpy_intensity_total", numpy_intensity_total)
        numpy_intensity_diff_total = numpy.array(self.intensity_diff_total, dtype=float)
        setattr(self, "__numpy_intensity_diff_total", numpy_intensity_diff_total)
        numpy_intensity_up = numpy.array(self.intensity_up, dtype=float)
        setattr(self, "__numpy_intensity_up", numpy_intensity_up)
        numpy_intensity_up_sigma = numpy.array(self.intensity_up_sigma, dtype=float)
        setattr(self, "__numpy_intensity_up_sigma", numpy_intensity_up_sigma)
        numpy_intensity_down = numpy.array(self.intensity_down, dtype=float)
        setattr(self, "__numpy_intensity_down", numpy_intensity_down)
        numpy_intensity_down_sigma = numpy.array(self.intensity_down_sigma, dtype=float)
        setattr(self, "__numpy_intensity_down_sigma", numpy_intensity_down_sigma)
        numpy_intensity = numpy.array(self.intensity, dtype=float)
        setattr(self, "__numpy_intensity", numpy_intensity)
        numpy_intensity_sigma = numpy.array(self.intensity_sigma, dtype=float)
        setattr(self, "__numpy_intensity_sigma", numpy_intensity_sigma)

    def transform_numpy_arrays_to_items(self):
        """
Transform data from numpy arrays to items:

    numpy_ttheta: 1D numpy array of ttheta, dtype=float
    numpy_ttheta_corrected: 1D numpy array of ttheta_corrected, dtype=float 
    numpy_intensity_up_net: 1D numpy array of intensity_up_net, dtype=float
    numpy_intensity_down_net: 1D numpy array of intensity_down_net, dtype=float
    numpy_intensity_up_total: 1D numpy array of intensity_up_total, dtype=float
    numpy_intensity_down_total 1D numpy array of intensity_down_total, dtype=floatl: 1D numpy array of , dtype=float
    numpy_intensity_total 1D numpy array of intensity_total, dtype=floatl: 1D numpy array of , dtype=float
    numpy_intensity_diff_total 1D numpy array of intensity_diff_total, dtype=floatl: 1D numpy array of , dtype=float
    numpy_intensity_bkg_calc 1D numpy array of intensity_bkg_calc, dtype=float_calc: 1D numpy array of , dtype=float
    numpy_intensity_up: 1D numpy array of intensity_up, dtype=float
    numpy_intensity_up_sigma: 1D numpy array of intensity_up_sigma, dtype=float
    numpy_intensity_down: 1D numpy array of intensity_down, dtype=float
    numpy_intensity_down_sigma: 1D numpy array of intensity_down_sigma, dtype=float
    numpy_intensity: 1D numpy array of intensity, dtype=float
    numpy_intensity_sigma: 1D numpy array of intensity_sigma, dtype=float
        """
        numpy_ttheta = getattr(self, "__numpy_ttheta")
        if numpy_ttheta is None: return
        l_item = [PdProc(ttheta=_val) for _val in numpy_ttheta]
        numpy_ttheta_corrected = getattr(self, "__numpy_ttheta_corrected")
        if numpy_ttheta_corrected is not None: 
            for _item, val in zip(l_item, numpy_ttheta_corrected):
                _item.ttheta_corrected = val
        numpy_intensity_up_net = getattr(self, "__numpy_intensity_up_net")
        if numpy_intensity_up_net is not None: 
            for _item, val in zip(l_item, numpy_intensity_up_net):
                _item.intensity_up_net = val
        numpy_intensity_down_net = getattr(self, "__numpy_intensity_down_net")
        if numpy_intensity_down_net is not None: 
            for _item, val in zip(l_item, numpy_intensity_down_net):
                _item.intensity_down_net = val
        numpy_intensity_up_total = getattr(self, "__numpy_intensity_up_total")
        if numpy_intensity_up_total is not None: 
            for _item, val in zip(l_item, numpy_intensity_up_total):
                _item.intensity_up_total = val
        numpy_intensity_down_total = getattr(self, "__numpy_intensity_down_total")
        if numpy_intensity_down_total is not None: 
            for _item, val in zip(l_item, numpy_intensity_down_total):
                _item.intensity_down_total = val
        numpy_intensity_bkg_calc = getattr(self, "__numpy_intensity_bkg_calc")
        if numpy_intensity_bkg_calc is not None: 
            for _item, val in zip(l_item, numpy_intensity_bkg_calc):
                _item.intensity_bkg_calc = val
        numpy_intensity_net = getattr(self, "__numpy_intensity_net")
        if numpy_intensity_net is not None: 
            for _item, val in zip(l_item, numpy_intensity_net):
                _item.intensity_net = val
        numpy_intensity_total = getattr(self, "__numpy_intensity_total")
        if numpy_intensity_total is not None: 
            for _item, val in zip(l_item, numpy_intensity_total):
                _item.intensity_total = val
        numpy_intensity_diff_total = getattr(self, "__numpy_intensity_diff_total")
        if numpy_intensity_diff_total is not None: 
            for _item, val in zip(l_item, numpy_intensity_diff_total):
                _item.intensity_diff_total = val
        numpy_intensity_up = getattr(self, "__numpy_intensity_up")
        if numpy_intensity_up is not None: 
            for _item, val in zip(l_item, numpy_intensity_up):
                _item.intensity_up = val
        numpy_intensity_up_sigma = getattr(self, "__numpy_intensity_up_sigma")
        if numpy_intensity_up_sigma is not None: 
            for _item, val in zip(l_item, numpy_intensity_up_sigma):
                _item.intensity_up_sigma = val
        numpy_intensity_down = getattr(self, "__numpy_intensity_down")
        if numpy_intensity_down is not None: 
            for _item, val in zip(l_item, numpy_intensity_down):
                _item.intensity_down = val
        numpy_intensity_down_sigma = getattr(self, "__numpy_intensity_down_sigma")
        if numpy_intensity_down_sigma is not None: 
            for _item, val in zip(l_item, numpy_intensity_down_sigma):
                _item.intensity_down_sigma = val
        numpy_intensity = getattr(self, "__numpy_intensity")
        if numpy_intensity is not None: 
            for _item, val in zip(l_item, numpy_intensity):
                _item.intensity = val
        numpy_intensity_sigma = getattr(self, "__numpy_intensity_sigma")
        if numpy_intensity_sigma is not None: 
            for _item, val in zip(l_item, numpy_intensity_sigma):
                _item.intensity_sigma = val
        self.item = l_item
