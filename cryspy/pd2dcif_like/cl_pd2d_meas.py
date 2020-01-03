__author__ = 'ikibalin'
__version__ = "2020_01_03"
import os
import numpy

import warnings
from typing import List, Tuple

from cryspy.common.cl_item_constr import ItemConstr
from cryspy.common.cl_loop_constr import LoopConstr
from cryspy.common.cl_fitable import Fitable

class Pd2dMeas(ItemConstr):
    """
This section contains the measured diffractogram and information
about the conditions used for the measurement of the diffraction 
data set, prior to processing and application of correction
terms.

Description in cif file::

 _pd2d_meas_ttheta_phi_intensity_up
 ;
 ;

 _pd2d_meas_ttheta_phi_intensity_up_sigma
 ;
 ;

 _pd2d_meas_ttheta_phi_intensity_down
 ;
 ;

 _pd2d_meas_ttheta_phi_intensity_down_sigma
 ;
 ;
    """
    MANDATORY_ATTRIBUTE = ("ttheta_phi_intensity_up", "ttheta_phi_intensity_up_sigma", 
                           "ttheta_phi_intensity_down", "ttheta_phi_intensity_down_sigma")
    OPTIONAL_ATTRIBUTE = ()
    INTERNAL_ATTRIBUTE = ("ttheta", "phi", "intensity_up", "intensity_up_sigma", "intensity_down", "intensity_down_sigma")
    PREFIX = "pd2d_meas"
    def __init__(self, ttheta_phi_intensity_up=None, ttheta_phi_intensity_up_sigma=None, 
                       ttheta_phi_intensity_down=None, ttheta_phi_intensity_down_sigma=None):
        super(Pd2dMeas, self).__init__(mandatory_attribute=self.MANDATORY_ATTRIBUTE, 
                                       optional_attribute=self.OPTIONAL_ATTRIBUTE, 
                                       internal_attribute=self.INTERNAL_ATTRIBUTE,
                                       prefix=self.PREFIX)

        self.ttheta_phi_intensity_up = ttheta_phi_intensity_up
        self.ttheta_phi_intensity_up_sigma = ttheta_phi_intensity_up_sigma
        self.ttheta_phi_intensity_down = ttheta_phi_intensity_down
        self.ttheta_phi_intensity_down_sigma = ttheta_phi_intensity_down_sigma

        if self.is_defined:
            self.form_object

    @property
    def ttheta_phi_intensity_up(self):
        self.form_ttheta_phi_intensity_up
        return getattr(self, "__ttheta_phi_intensity_up")
    @ttheta_phi_intensity_up.setter
    def ttheta_phi_intensity_up(self, x):
        if x is None:
            x_in = None
        else:
            x_in = str(x)
        setattr(self, "__ttheta_phi_intensity_up", x_in)
        if x_in is not None:
            self.form_object

    @property
    def ttheta_phi_intensity_up_sigma(self):
        self.form_ttheta_phi_intensity_up_sigma
        return getattr(self, "__ttheta_phi_intensity_up_sigma")
    @ttheta_phi_intensity_up_sigma.setter
    def ttheta_phi_intensity_up_sigma(self, x):
        if x is None:
            x_in = None
        else:
            x_in = str(x)
        setattr(self, "__ttheta_phi_intensity_up_sigma", x_in)
        if x_in is not None:
            self.form_object

    @property
    def ttheta_phi_intensity_down(self):
        self.form_ttheta_phi_intensity_down
        return getattr(self, "__ttheta_phi_intensity_down")
    @ttheta_phi_intensity_down.setter
    def ttheta_phi_intensity_down(self, x):
        if x is None:
            x_in = None
        else:
            x_in = str(x)
        setattr(self, "__ttheta_phi_intensity_down", x_in)
        if x_in is not None:
            self.form_object

    @property
    def ttheta_phi_intensity_down_sigma(self):
        self.form_ttheta_phi_intensity_down_sigma
        return getattr(self, "__ttheta_phi_intensity_down_sigma")
    @ttheta_phi_intensity_down_sigma.setter
    def ttheta_phi_intensity_down_sigma(self, x):
        if x is None:
            x_in = None
        else:
            x_in = str(x)
        setattr(self, "__ttheta_phi_intensity_down_sigma", x_in)
        if x_in is not None:
            self.form_object

    @property
    def ttheta(self):
        return getattr(self, "__ttheta")

    @property
    def phi(self):
        return getattr(self, "__phi")

    @property
    def intensity_up(self):
        return getattr(self, "__intensity_up")
    @property
    def intensity_up_sigma(self):
        return getattr(self, "__intensity_up_sigma")
    @property
    def intensity_down(self):
        return getattr(self, "__intensity_down")
    @property
    def intensity_down_sigma(self):
        return getattr(self, "__intensity_down_sigma")

    def __repr__(self):
        ls_out = []
        ls_out.append(f"Pd2dMeas:\n{str(self):}")
        return "\n".join(ls_out)
            
    def _show_message(self, s_out: str):
        warnings.warn("***  Error ***\n"+s_out, UserWarning, stacklevel=2)

    @property
    def form_object(self) -> bool:
        flag = True
        l_1 = (self.ttheta_phi_intensity_up).strip().split("\n")
        l_2 = (self.ttheta_phi_intensity_up_sigma).strip().split("\n")
        l_3 = (self.ttheta_phi_intensity_down).strip().split("\n")
        l_4 = (self.ttheta_phi_intensity_down_sigma).strip().split("\n")

        l_ttheta = numpy.array([_ for _ in l_1[0].strip().split()[1:]], dtype=float)
        l_phi, ll_intensity_up, ll_intensity_up_sigma = [], [], []
        ll_intensity_down, ll_intensity_down_sigma = [], []
        for line_1, line_2, line_3, line_4 in zip(l_1[1:], l_2[1:], l_3[1:], l_4[1:]):
            _l_1 = line_1.strip().split()
            _l_2 = line_2.strip().split()
            _l_3 = line_3.strip().split()
            _l_4 = line_4.strip().split()
            l_phi.append(float(_l_1[0]))
            ll_intensity_up.append(_l_1[1:])
            ll_intensity_up_sigma.append(_l_2[1:])
            ll_intensity_down.append(_l_3[1:])
            ll_intensity_down_sigma.append(_l_4[1:])

        ll_intensity_up = numpy.array(ll_intensity_up, dtype=float).transpose()
        ll_intensity_up_sigma = numpy.array(ll_intensity_up_sigma, dtype=float).transpose()
        ll_intensity_down = numpy.array(ll_intensity_down, dtype=float).transpose()
        ll_intensity_down_sigma = numpy.array(ll_intensity_down_sigma, dtype=float).transpose()

        setattr(self,"__ttheta", l_ttheta)
        setattr(self,"__phi", numpy.array(l_phi, dtype=float))
        setattr(self,"__intensity_up", ll_intensity_up)
        setattr(self,"__intensity_up_sigma", ll_intensity_up_sigma)
        setattr(self,"__intensity_down", ll_intensity_down)
        setattr(self,"__intensity_down_sigma", ll_intensity_down_sigma)
        return flag

    @property
    def form_ttheta_phi_intensity_up(self) -> bool:
        if ((self.phi is not None) & (self.ttheta is not None) & (self.intensity_up is not None)):
            ls_out = []
            ls_out.append("{:12} ".format(len(self.phi)) + " ".join(["{:6.2f}      ".format(_) for _ in self.ttheta]))
            ll_intensity = self.intensity_up
            ll_intensity = [[ll_intensity[_2][_1] for _2 in range(len(ll_intensity))] for _1 in range(len(ll_intensity[0]))]
            for phi, l_intensity in zip(self.phi, ll_intensity):
                ls_out.append("{:12.2f} ".format(phi) + " ".join(["{:12}".format(_.print_with_sigma) for _ in l_intensity]))
            setattr(self, "__ttheta_phi_intensity_up", "\n".join(ls_out))

    @property
    def form_ttheta_phi_intensity_up_sigma(self) -> bool:
        if ((self.phi is not None) & (self.ttheta is not None) & (self.intensity_up is not None)):
            ls_out = []
            ls_out.append("{:12} ".format(len(self.phi)) + " ".join(["{:6.2f}      ".format(_) for _ in self.ttheta]))
            ll_intensity = self.intensity_up_sigma
            ll_intensity = [[ll_intensity[_2][_1] for _2 in range(len(ll_intensity))] for _1 in range(len(ll_intensity[0]))]
            for phi, l_intensity in zip(self.phi, ll_intensity):
                ls_out.append("{:12.2f} ".format(phi) + " ".join(["{:12}".format(_.print_with_sigma) for _ in l_intensity]))
            setattr(self, "__ttheta_phi_intensity_up_sigma", "\n".join(ls_out))

    @property
    def form_ttheta_phi_intensity_down(self) -> bool:
        if ((self.phi is not None) & (self.ttheta is not None) & (self.intensity_up is not None)):
            ls_out = []
            ls_out.append("{:12} ".format(len(self.phi)) + " ".join(["{:6.2f}      ".format(_) for _ in self.ttheta]))
            ll_intensity = self.intensity_down
            ll_intensity = [[ll_intensity[_2][_1] for _2 in range(len(ll_intensity))] for _1 in range(len(ll_intensity[0]))]
            for phi, l_intensity in zip(self.phi, ll_intensity):
                ls_out.append("{:12.2f} ".format(phi) + " ".join(["{:12}".format(_.print_with_sigma) for _ in l_intensity]))
            setattr(self, "__ttheta_phi_intensity_down", "\n".join(ls_out))

    @property
    def form_ttheta_phi_intensity_down_sigma(self) -> bool:
        if ((self.phi is not None) & (self.ttheta is not None) & (self.intensity_up is not None)):
            ls_out = []
            ls_out.append("{:12} ".format(len(self.phi)) + " ".join(["{:6.2f}      ".format(_) for _ in self.ttheta]))
            ll_intensity = self.intensity_down_sigma
            ll_intensity = [[ll_intensity[_2][_1] for _2 in range(len(ll_intensity))] for _1 in range(len(ll_intensity[0]))]
            for phi, l_intensity in zip(self.phi, ll_intensity):
                ls_out.append("{:12.2f} ".format(phi) + " ".join(["{:12}".format(_.print_with_sigma) for _ in l_intensity]))
            setattr(self, "__ttheta_phi_intensity_down_sigma", "\n".join(ls_out))
