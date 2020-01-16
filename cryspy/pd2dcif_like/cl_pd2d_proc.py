__author__ = 'ikibalin'
__version__ = "2019_12_10"
import os
import numpy

import warnings
from typing import List, Tuple

from cryspy.common.cl_item_constr import ItemConstr
from cryspy.common.cl_loop_constr import LoopConstr
from cryspy.common.cl_fitable import Fitable


class Pd2dProc(ItemConstr):
    """
This section contains the diffraction data set after processing
and application of correction terms. 

Description in cif file::

 _pd2d_proc_2theta_phi_intensity_up_net
 ;
      2    4.5     40.0     80.0
 -3.000 -356.0   -350.0   -400.0
 41.000 -357.0   -350.0   -400.0
 ;
 _pd2d_proc_2theta_phi_intensity_down_net
 ;
      2    4.5     40.0     80.0
 -3.000 -356.0   -350.0   -400.0
 41.000 -357.0   -350.0   -400.0
 ;
 _pd2d_proc_2theta_phi_intensity_up_total
 ;
      2    4.5     40.0     80.0
 -3.000 -356.0   -350.0   -400.0
 41.000 -357.0   -350.0   -400.0
 ;
 _pd2d_proc_2theta_phi_intensity_down_total
 ;
      2    4.5     40.0     80.0
 -3.000 -356.0   -350.0   -400.0
 41.000 -357.0   -350.0   -400.0
 ;
 _pd2d_proc_2theta_phi_intensity_bkg_calc
 ;
      2    4.5     40.0     80.0
 -3.000 -356.0   -350.0   -400.0
 41.000 -357.0   -350.0   -400.0
 ;
 _pd2d_proc_2theta_phi_intensity_up
 ;
      2    4.5     40.0     80.0
 -3.000 -356.0   -350.0   -400.0
 41.000 -357.0   -350.0   -400.0
 ;
 _pd2d_proc_2theta_phi_intensity_up_sigma
 ;
      2    4.5     40.0     80.0
 -3.000 -356.0   -350.0   -400.0
 41.000 -357.0   -350.0   -400.0
 ;
 _pd2d_proc_2theta_phi_intensity_down
 ;
      2    4.5     40.0     80.0
 -3.000 -356.0   -350.0   -400.0
 41.000 -357.0   -350.0   -400.0
 ;
 _pd2d_proc_2theta_phi_intensity_down_sigma
 ;
      2    4.5     40.0     80.0
 -3.000 -356.0   -350.0   -400.0
 41.000 -357.0   -350.0   -400.0
 ;
    """
    MANDATORY_ATTRIBUTE = ("ttheta_phi_intensity_up_net", "ttheta_phi_intensity_down_net",
        "ttheta_phi_intensity_up_total", "ttheta_phi_intensity_down_total", "ttheta_phi_intensity_bkg_calc",
        "ttheta_phi_intensity_up", "ttheta_phi_intensity_up_sigma", "ttheta_phi_intensity_down",
        "ttheta_phi_intensity_down_sigma")
    OPTIONAL_ATTRIBUTE = ()
    RELATED_CIF_MANDATORY_ATTRIBUTE = ("2theta_phi_intensity_up_net", "2theta_phi_intensity_down_net",
        "2theta_phi_intensity_up_total", "2theta_phi_intensity_down_total", "2theta_phi_intensity_bkg_calc",
        "2theta_phi_intensity_up", "2theta_phi_intensity_up_sigma", "2theta_phi_intensity_down",
        "2theta_phi_intensity_down_sigma")
    RELATED_CIF_OPTIONAL_ATTRIBUTE = ()
    INTERNAL_ATTRIBUTE = ("ttheta", "phi", "ttheta_corrected", "d_spacing",
                          "intensity_up_net", "intensity_down_net", "intensity_up_total", "intensity_down_total", 
                          "intensity_bkg_calc", "intensity_net", "intensity_total", 
                          "intensity_up", "intensity_up_sigma", "intensity_down", "intensity_down_sigma", 
                          "intensity", "intensity_sigma")
    PREFIX = "pd2d_proc"
    def __init__(self, ttheta_phi_intensity_up_net=None, ttheta_phi_intensity_down_net=None, 
                 ttheta_phi_intensity_up_total=None, ttheta_phi_intensity_down_total=None, 
                 ttheta_phi_intensity_bkg_calc=None, ttheta_phi_intensity_up=None, 
                 ttheta_phi_intensity_up_sigma=None, ttheta_phi_intensity_down=None, 
                 ttheta_phi_intensity_down_sigma=None):

        super(Pd2dProc, self).__init__(mandatory_attribute=self.MANDATORY_ATTRIBUTE, 
                                       optional_attribute=self.OPTIONAL_ATTRIBUTE, 
                                       internal_attribute=self.INTERNAL_ATTRIBUTE,
                                       prefix=self.PREFIX)

        self.ttheta_phi_intensity_up_net = ttheta_phi_intensity_up_net
        self.ttheta_phi_intensity_down_net = ttheta_phi_intensity_down_net
        self.ttheta_phi_intensity_up_total = ttheta_phi_intensity_up_total
        self.ttheta_phi_intensity_down_total = ttheta_phi_intensity_down_total
        self.ttheta_phi_intensity_bkg_calc = ttheta_phi_intensity_bkg_calc
        self.ttheta_phi_intensity_up = ttheta_phi_intensity_up
        self.ttheta_phi_intensity_up_sigma = ttheta_phi_intensity_up_sigma
        self.ttheta_phi_intensity_down = ttheta_phi_intensity_down
        self.ttheta_phi_intensity_down_sigma = ttheta_phi_intensity_down_sigma

        if self.is_defined:
            self.form_object

    @property
    def ttheta_phi_intensity_up_net(self):
        self.form_ttheta_phi_intensity_up_net
        return getattr(self, "__ttheta_phi_intensity_up_net")
    @ttheta_phi_intensity_up_net.setter
    def ttheta_phi_intensity_up_net(self, x):
        if x is None:
            x_in = None
        else:
            x_in = str(x)
        setattr(self, "__ttheta_phi_intensity_up_net", x_in)
        if x_in is not None:
            self.form_object

    @property
    def ttheta_phi_intensity_down_net(self):
        self.form_ttheta_phi_intensity_down_net
        return getattr(self, "__ttheta_phi_intensity_down_net")
    @ttheta_phi_intensity_down_net.setter
    def ttheta_phi_intensity_down_net(self, x):
        if x is None:
            x_in = None
        else:
            x_in = str(x)
        setattr(self, "__ttheta_phi_intensity_down_net", x_in)
        if x_in is not None:
            self.form_object

    @property
    def ttheta_phi_intensity_up_total(self):
        self.form_ttheta_phi_intensity_up_total
        return getattr(self, "__ttheta_phi_intensity_up_total")
    @ttheta_phi_intensity_up_total.setter
    def ttheta_phi_intensity_up_total(self, x):
        if x is None:
            x_in = None
        else:
            x_in = str(x)
        setattr(self, "__ttheta_phi_intensity_up_total", x_in)
        if x_in is not None:
            self.form_object

    @property
    def ttheta_phi_intensity_down_total(self):
        self.form_ttheta_phi_intensity_down_total
        return getattr(self, "__ttheta_phi_intensity_down_total")
    @ttheta_phi_intensity_down_total.setter
    def ttheta_phi_intensity_down_total(self, x):
        if x is None:
            x_in = None
        else:
            x_in = str(x)
        setattr(self, "__ttheta_phi_intensity_down_total", x_in)
        if x_in is not None:
            self.form_object

    @property
    def ttheta_phi_intensity_bkg_calc(self):
        self.form_ttheta_phi_intensity_bkg_calc
        return getattr(self, "__ttheta_phi_intensity_bkg_calc")
    @ttheta_phi_intensity_bkg_calc.setter
    def ttheta_phi_intensity_bkg_calc(self, x):
        if x is None:
            x_in = None
        else:
            x_in = str(x)
        setattr(self, "__ttheta_phi_intensity_bkg_calc", x_in)
        if x_in is not None:
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
    def intensity_up_net(self):
        return getattr(self, "__intensity_up_net")
    @property
    def intensity_down_net(self):
        return getattr(self, "__intensity_down_net")
    @property
    def intensity_up_total(self):
        return getattr(self, "__intensity_up_total")
    @property
    def intensity_down_total(self):
        return getattr(self, "__intensity_down_total")
    @property
    def intensity_bkg_calc(self):
        return getattr(self, "__intensity_bkg_calc")
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
        ls_out.append(f"Pd2dProc:\n{str(self):}")
        return "\n".join(ls_out)


    @property
    def form_object(self) -> bool:
        flag = True
        if any([self.ttheta_phi_intensity_up_net is None, 
                self.ttheta_phi_intensity_down_net is None, 
                self.ttheta_phi_intensity_up_total is None, 
                self.ttheta_phi_intensity_down_total is None, 
                self.ttheta_phi_intensity_bkg_calc is None, 
                self.ttheta_phi_intensity_up is None, 
                self.ttheta_phi_intensity_up_sigma is None,
                self.ttheta_phi_intensity_down is None, 
                self.ttheta_phi_intensity_down_sigma is None]):
            return False
        l_1 = (self.ttheta_phi_intensity_up_net).strip().split("\n")
        l_2 = (self.ttheta_phi_intensity_down_net).strip().split("\n")
        l_3 = (self.ttheta_phi_intensity_up_total).strip().split("\n")
        l_4 = (self.ttheta_phi_intensity_down_total).strip().split("\n")
        l_5 = (self.ttheta_phi_intensity_bkg_calc).strip().split("\n")
        l_6 = (self.ttheta_phi_intensity_up).strip().split("\n")
        l_7 = (self.ttheta_phi_intensity_up_sigma).strip().split("\n")
        l_8 = (self.ttheta_phi_intensity_down).strip().split("\n")
        l_9 = (self.ttheta_phi_intensity_down_sigma).strip().split("\n")

        l_ttheta = numpy.array([_ for _ in l_1[0].strip().split()[1:]], dtype=float)
        l_phi, ll_intensity_up, ll_intensity_up_sigma = [], [], []
        ll_intensity_down, ll_intensity_down_sigma = [], []
        ll_intensity_up_net, ll_intensity_down_net = [], []
        ll_intensity_up_total, ll_intensity_down_total = [], []
        ll_intensity_bkg_calc = []
        for line_1, line_2, line_3, line_4, line_5, line_6, line_7, line_8, line_9 in zip(
            l_1[1:], l_2[1:], l_3[1:], l_4[1:], l_5[1:], l_6[1:], l_7[1:], l_8[1:], l_9[1:]):
            _l_1 = line_1.strip().split()
            _l_2 = line_2.strip().split()
            _l_3 = line_3.strip().split()
            _l_4 = line_4.strip().split()
            _l_5 = line_5.strip().split()
            _l_6 = line_6.strip().split()
            _l_7 = line_7.strip().split()
            _l_8 = line_8.strip().split()
            _l_9 = line_9.strip().split()
            l_phi.append(float(_l_1[0]))
            ll_intensity_up_net.append(_l_1[1:])
            ll_intensity_down_net.append(_l_2[1:])
            ll_intensity_up_total.append(_l_3[1:])
            ll_intensity_down_total.append(_l_4[1:])
            ll_intensity_bkg_calc.append(_l_5[1:])
            ll_intensity_up.append(_l_6[1:])
            ll_intensity_up_sigma.append(_l_7[1:])
            ll_intensity_down.append(_l_8[1:])
            ll_intensity_down_sigma.append(_l_9[1:])

        ll_intensity_up_net = numpy.array(ll_intensity_up_net, dtype=float).transpose()
        ll_intensity_down_net = numpy.array(ll_intensity_down_net, dtype=float).transpose()
        ll_intensity_up_total = numpy.array(ll_intensity_up_total, dtype=float).transpose()
        ll_intensity_down_total = numpy.array(ll_intensity_down_total, dtype=float).transpose()
        ll_intensity_bkg_calc = numpy.array(ll_intensity_bkg_calc, dtype=float).transpose()
        ll_intensity_up = numpy.array(ll_intensity_up, dtype=float).transpose()
        ll_intensity_up_sigma = numpy.array(ll_intensity_up_sigma, dtype=float).transpose()
        ll_intensity_down = numpy.array(ll_intensity_down, dtype=float).transpose()
        ll_intensity_down_sigma = numpy.array(ll_intensity_down_sigma, dtype=float).transpose()

        setattr(self,"__ttheta", l_ttheta)
        setattr(self,"__phi", numpy.array(l_phi, dtype=float))
        setattr(self,"__intensity_up_net", ll_intensity_up_net)
        setattr(self,"__intensity_down_net", ll_intensity_down_net)
        setattr(self,"__intensity_up_total", ll_intensity_up_total)
        setattr(self,"__intensity_down_total", ll_intensity_down_total)
        setattr(self,"__intensity_bkg_calc", ll_intensity_bkg_calc)
        setattr(self,"__intensity_up", ll_intensity_up)
        setattr(self,"__intensity_up_sigma", ll_intensity_up_sigma)
        setattr(self,"__intensity_down", ll_intensity_down)
        setattr(self,"__intensity_down_sigma", ll_intensity_down_sigma)
        return flag

    @property
    def form_ttheta_phi_intensity_up_net(self) -> bool:
        if ((self.phi is not None) & (self.ttheta is not None) & (self.intensity_up_net is not None)):
            ls_out = []
            ls_out.append("{:12} ".format(len(self.phi)) + " ".join(["{:6.2f}      ".format(_) for _ in self.ttheta]))
            ll_intensity = self.intensity_up_net
            ll_intensity = [[ll_intensity[_2][_1] for _2 in range(len(ll_intensity))] for _1 in range(len(ll_intensity[0]))]
            for phi, l_intensity in zip(self.phi, ll_intensity):
                ls_out.append("{:12.2f} ".format(phi) + " ".join(["{:12}".format(_) for _ in l_intensity]))
            setattr(self, "__ttheta_phi_intensity_up_net", "\n".join(ls_out))

    @property
    def form_ttheta_phi_intensity_down_net(self) -> bool:
        if ((self.phi is not None) & (self.ttheta is not None) & (self.intensity_down_net is not None)):
            ls_out = []
            ls_out.append("{:12} ".format(len(self.phi)) + " ".join(["{:6.2f}      ".format(_) for _ in self.ttheta]))
            ll_intensity = self.intensity_down_net
            ll_intensity = [[ll_intensity[_2][_1] for _2 in range(len(ll_intensity))] for _1 in range(len(ll_intensity[0]))]
            for phi, l_intensity in zip(self.phi, ll_intensity):
                ls_out.append("{:12.2f} ".format(phi) + " ".join(["{:12}".format(_) for _ in l_intensity]))
            setattr(self, "__ttheta_phi_intensity_down_net", "\n".join(ls_out))

    @property
    def form_ttheta_phi_intensity_up_total(self) -> bool:
        if ((self.phi is not None) & (self.ttheta is not None) & (self.intensity_up_total is not None)):
            ls_out = []
            ls_out.append("{:12} ".format(len(self.phi)) + " ".join(["{:6.2f}      ".format(_) for _ in self.ttheta]))
            ll_intensity = self.intensity_up_total
            ll_intensity = [[ll_intensity[_2][_1] for _2 in range(len(ll_intensity))] for _1 in range(len(ll_intensity[0]))]
            for phi, l_intensity in zip(self.phi, ll_intensity):
                ls_out.append("{:12.2f} ".format(phi) + " ".join(["{:12}".format(_) for _ in l_intensity]))
            setattr(self, "__ttheta_phi_intensity_up_total", "\n".join(ls_out))

    @property
    def form_ttheta_phi_intensity_down_total(self) -> bool:
        if ((self.phi is not None) & (self.ttheta is not None) & (self.intensity_down_total is not None)):
            ls_out = []
            ls_out.append("{:12} ".format(len(self.phi)) + " ".join(["{:6.2f}      ".format(_) for _ in self.ttheta]))
            ll_intensity = self.intensity_down_total
            ll_intensity = [[ll_intensity[_2][_1] for _2 in range(len(ll_intensity))] for _1 in range(len(ll_intensity[0]))]
            for phi, l_intensity in zip(self.phi, ll_intensity):
                ls_out.append("{:12.2f} ".format(phi) + " ".join(["{:12}".format(_) for _ in l_intensity]))
            setattr(self, "__ttheta_phi_intensity_down_total", "\n".join(ls_out))

    @property
    def form_ttheta_phi_intensity_bkg_calc(self) -> bool:
        if ((self.phi is not None) & (self.ttheta is not None) & (self.intensity_bkg_calc is not None)):
            ls_out = []
            ls_out.append("{:12} ".format(len(self.phi)) + " ".join(["{:6.2f}      ".format(_) for _ in self.ttheta]))
            ll_intensity = self.intensity_bkg_calc
            ll_intensity = [[ll_intensity[_2][_1] for _2 in range(len(ll_intensity))] for _1 in range(len(ll_intensity[0]))]
            for phi, l_intensity in zip(self.phi, ll_intensity):
                ls_out.append("{:12.2f} ".format(phi) + " ".join(["{:12}".format(_) for _ in l_intensity]))
            setattr(self, "__ttheta_phi_intensity_bkg_calc", "\n".join(ls_out))


    @property
    def form_ttheta_phi_intensity_up(self) -> bool:
        if ((self.phi is not None) & (self.ttheta is not None) & (self.intensity_up is not None)):
            ls_out = []
            ls_out.append("{:12} ".format(len(self.phi)) + " ".join(["{:6.2f}      ".format(_) for _ in self.ttheta]))
            ll_intensity = self.intensity_up
            ll_intensity = [[ll_intensity[_2][_1] for _2 in range(len(ll_intensity))] for _1 in range(len(ll_intensity[0]))]
            for phi, l_intensity in zip(self.phi, ll_intensity):
                ls_out.append("{:12.2f} ".format(phi) + " ".join(["{:12}".format(_) for _ in l_intensity]))
            setattr(self, "__ttheta_phi_intensity_up", "\n".join(ls_out))

    @property
    def form_ttheta_phi_intensity_up_sigma(self) -> bool:
        if ((self.phi is not None) & (self.ttheta is not None) & (self.intensity_up_sigma is not None)):
            ls_out = []
            ls_out.append("{:12} ".format(len(self.phi)) + " ".join(["{:6.2f}      ".format(_) for _ in self.ttheta]))
            ll_intensity = self.intensity_up_sigma
            ll_intensity = [[ll_intensity[_2][_1] for _2 in range(len(ll_intensity))] for _1 in range(len(ll_intensity[0]))]
            for phi, l_intensity in zip(self.phi, ll_intensity):
                ls_out.append("{:12.2f} ".format(phi) + " ".join(["{:12}".format(_) for _ in l_intensity]))
            setattr(self, "__ttheta_phi_intensity_up_sigma", "\n".join(ls_out))

    @property
    def form_ttheta_phi_intensity_down(self) -> bool:
        if ((self.phi is not None) & (self.ttheta is not None) & (self.intensity_down is not None)):
            ls_out = []
            ls_out.append("{:12} ".format(len(self.phi)) + " ".join(["{:6.2f}      ".format(_) for _ in self.ttheta]))
            ll_intensity = self.intensity_down
            ll_intensity = [[ll_intensity[_2][_1] for _2 in range(len(ll_intensity))] for _1 in range(len(ll_intensity[0]))]
            for phi, l_intensity in zip(self.phi, ll_intensity):
                ls_out.append("{:12.2f} ".format(phi) + " ".join(["{:12}".format(_) for _ in l_intensity]))
            setattr(self, "__ttheta_phi_intensity_down", "\n".join(ls_out))

    @property
    def form_ttheta_phi_intensity_down_sigma(self) -> bool:
        if ((self.phi is not None) & (self.ttheta is not None) & (self.intensity_down_sigma is not None)):
            ls_out = []
            ls_out.append("{:12} ".format(len(self.phi)) + " ".join(["{:6.2f}      ".format(_) for _ in self.ttheta]))
            ll_intensity = self.intensity_down_sigma
            ll_intensity = [[ll_intensity[_2][_1] for _2 in range(len(ll_intensity))] for _1 in range(len(ll_intensity[0]))]
            for phi, l_intensity in zip(self.phi, ll_intensity):
                ls_out.append("{:12.2f} ".format(phi) + " ".join(["{:12}".format(_) for _ in l_intensity]))
            setattr(self, "__ttheta_phi_intensity_down_sigma", "\n".join(ls_out))
