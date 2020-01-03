__author__ = 'ikibalin'
__version__ = "2020_01_03"
import os
import numpy

import warnings
from typing import List, Tuple

from cryspy.common.cl_item_constr import ItemConstr
from cryspy.common.cl_loop_constr import LoopConstr
from cryspy.common.cl_fitable import Fitable


class Pd2dBackground(ItemConstr):
    """
Description in cif file::

 _pd2d_background_ttheta_phi_intensity
 ;
      2    4.5     40.0     80.0
 -3.000 -350.0   -350.0   -400.0
 41.000 -350.0   -350.0   -400.0
 ;
    """
    MANDATORY_ATTRIBUTE = ("ttheta_phi_intensity", )
    OPTIONAL_ATTRIBUTE = ()
    INTERNAL_ATTRIBUTE = ("ttheta", "phi", "intensity")
    PREFIX = "pd2d_background"
    def __init__(self, ttheta_phi_intensity=None):
        super(Pd2dBackground, self).__init__(mandatory_attribute=self.MANDATORY_ATTRIBUTE, 
                                           optional_attribute=self.OPTIONAL_ATTRIBUTE, 
                                           internal_attribute=self.INTERNAL_ATTRIBUTE,
                                           prefix=self.PREFIX)

        self.ttheta_phi_intensity = ttheta_phi_intensity

        if self.is_defined:
            self.form_object

    @property
    def ttheta_phi_intensity(self):
        self.form_ttheta_phi_intensity
        return getattr(self, "__ttheta_phi_intensity")
    @ttheta_phi_intensity.setter
    def ttheta_phi_intensity(self, x):
        if x is None:
            x_in = None
        else:
            x_in = str(x)
        setattr(self, "__ttheta_phi_intensity", x_in)
        if x_in is not None:
            self.form_object

    @property
    def ttheta(self):
        return getattr(self, "__ttheta")

    @property
    def phi(self):
        return getattr(self, "__phi")

    @property
    def intensity(self):
        return getattr(self, "__intensity")

    def __repr__(self):
        ls_out = []
        ls_out.append(f"Pd2dBackground:\n{str(self):}")
        return "\n".join(ls_out)

    def _show_message(self, s_out: str):
        warnings.warn("***  Error ***\n"+s_out, UserWarning, stacklevel=2)

    @property
    def form_object(self) -> bool:
        flag = True
        l_1 = (self.ttheta_phi_intensity).strip().split("\n")

        l_ttheta = [float(_) for _ in l_1[0].strip().split()[1:]]
        l_phi, ll_intensity = [], []
        for line in l_1[1:]:
            l_1 = line.strip().split()
            l_phi.append(float(l_1[0]))
            ll_intensity.append(l_1[1:])
        ll_intensity = [[Fitable.from_object(ll_intensity[_2][_1]) for _2 in range(len(ll_intensity))] for _1 in range(len(ll_intensity[0]))]

        setattr(self,"__ttheta", l_ttheta)
        setattr(self,"__phi", l_phi)
        setattr(self,"__intensity", ll_intensity)
        return flag

    @property
    def form_ttheta_phi_intensity(self) -> bool:
        if ((self.phi is not None) & (self.ttheta is not None) & (self.intensity is not None)):
            ls_out = []
            ls_out.append("{:12} ".format(len(self.phi)) + " ".join(["{:6.2f}      ".format(_) for _ in self.ttheta]))
            ll_intensity = self.intensity
            ll_intensity = [[ll_intensity[_2][_1] for _2 in range(len(ll_intensity))] for _1 in range(len(ll_intensity[0]))]
            for phi, l_intensity in zip(self.phi, ll_intensity):
                ls_out.append("{:12.2f} ".format(phi) + " ".join(["{:12}".format(_.print_with_sigma) for _ in l_intensity]))
            setattr(self, "__ttheta_phi_intensity", "\n".join(ls_out))

    @property
    def is_variable(self):
        res = any([[_2.refinement for _2 in _1] for _1 in self.intensity])
        return res
    
    def get_variables(self):
        l_variable = []
        for _1 in self.intensity:
            for _2 in _1:
                if _2.refinement:
                    l_variable.append(_2)
        return l_variable


    def interpolate_by_points(self, tth, phi):
        l_phi_b = self.phi
        l_tth_b = self.ttheta
        ll_int_b = self.intensity
        ll_int_b = [[float(ll_int_b[_2][_1]) for _2 in range(len(ll_int_b))] for _1 in range(len(ll_int_b[0]))]
        if len(l_tth_b) == 0:
            int_2d = numpy.zeros((tth.size, phi.size), dtype=float)
        else:
            phi_b = numpy.array(l_phi_b, dtype=float)
            tth_b = numpy.array(l_tth_b, dtype=float)
            int_b = numpy.array(ll_int_b, dtype=float)
            func = scipy.interpolate.interp2d(tth_b, phi_b, int_b)
            #tth_2d, phi_2d = numpy.meshgrid(tth, phi, indexing="ij")
            int_2d = func(tth, phi)
            int_2d = int_2d.transpose()
        return int_2d

