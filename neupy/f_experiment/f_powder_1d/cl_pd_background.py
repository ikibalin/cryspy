"""
define class PdBackground which describes the background in 1d powder diffraction experiment
"""
__author__ = 'ikibalin'
__version__ = "2019_09_06"
import os
import numpy


from pystar import CIFglobal
from neupy.f_common.cl_fitable import Fitable


class PdBackground(object):
    """
    Example:

    loop_
    _pd_background_2theta
    _pd_background_intensity
     4.5  256.0
     40.0  158.0
     80.0  65.0
    """
    def __init__(self, ttheta=[], intensity=[]):
        super(PdBackground, self).__init__()
        self.__pd_background_2theta = []
        self.__pd_background_intensity = []
        self.ttheta = ttheta
        self.intensity = intensity

    @property
    def ttheta(self):
        return self.__pd_background_2theta
    @ttheta.setter
    def ttheta(self, l_x):
        l_x_in = []
        for x in l_x:
            x_in = float(x)
            l_x_in.append(x_in)
        self.__pd_background_2theta = l_x_in
        len_x = len(l_x_in)

        len_1 = len(self.intensity)
        if len_1 > len_x:
            self.__pd_background_intensity = self.__pd_background_intensity[:len_x]
        elif len_1 < len_x:
            l_fitable = [Fitable(value=0., name="_pd_background_intensity") for hh in range(len_x-len_1)] #default
            self.__pd_background_intensity.extend(l_fitable)

    @property
    def intensity(self):
        return tuple(self.__pd_background_intensity)
    @intensity.setter
    def intensity(self, l_x):
        l_fitable = []
        for x in l_x:
            if isinstance(x, Fitable):
                x_in = x
            else:
                x_in = Fitable()
                flag = x_in.take_it(x)
            l_fitable.append(x_in)
        len_x = len(l_fitable)
        len_1 = len(self.ttheta)
        if len_1 < len_x:
            l_fitable = l_fitable[:len_1]
        elif len_1 > len_x:
            l_fitable.extend([Fitable(value=0., name="_pd_background_intensity") for hh in range(len_1-len_x)])
        self.__pd_background_intensity = l_fitable



            
    def __repr__(self):
        ls_out = ["PdBackground:"]
        ls_out.append(" ttheta intensity")
        for _1, _2 in zip(self.ttheta, self.intensity):
            ls_out.append("{:7.2f} {:}".format(_1, _2.print_with_sigma))
        return "\n".join(ls_out)

    @property
    def to_cif(self):
        ls_out = []
        if self.is_defined:
            ls_out.append("loop_")
            ls_out.append("_pd_background_2theta")
            ls_out.append("_pd_background_intensity")
            for _1, _2 in zip(self.ttheta, self.intensity):
                ls_out.append("{:} {:}".format(_1, _2))
        return "\n".join(ls_out)

    def from_cif(self, string: str):
        cif_global = CIFglobal()
        flag = cif_global.take_from_string(string)
        if not flag:
            return False
        flag = False
        flag = cif_global.is_prefix("_pd_background")
        if flag:
            cif_loop = cif_global["_pd_background"]
            l_name = cif_loop.names
            if "_pd_background_2theta" in l_name:
                self.ttheta = [float(_1) for _1 in cif_loop["_pd_background_2theta"]]
            if "_pd_background_intensity" in l_name:
                self.intensity = [str(_1) for _1 in cif_loop["_pd_background_intensity"]]
        else:
            self.ttheta, self.intensity = [], []
        return True

    @property
    def is_defined(self):
        cond = all([self.ttheta is not None, self.intensity is not None])
        return cond

    @property
    def is_variable(self):
        res = any([_1.refinement for _1 in self.intensity])
        return res
    
    def get_variables(self):
        l_variable = []
        l_variable.extend([_1 for _1 in self.intensity if _1.refinement])
        return l_variable

    def _show_message(self, s_out: str):
        print("***  Error ***")
        print(s_out)

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