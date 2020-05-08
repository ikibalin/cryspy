"""
define class Pd2dBackground which describes the background in 2d powder diffraction experiment
"""
__author__ = 'ikibalin'
__version__ = "2019_09_10"
import os
import numpy
import scipy
import scipy.interpolate


from pycifstar import Global
from cryspy.f_common.cl_fitable import Fitable


class Pd2dBackground(object):
    """
    Example:

    _pd2d_background_2theta_phi_intensity
    ;
         2    4.5     40.0     80.0
    -3.000 -350.0   -350.0   -400.0
    41.000 -350.0   -350.0   -400.0
    ;
    """
    def __init__(self, ttheta=[], phi=[], intensity=[[]]):
        super(Pd2dBackground, self).__init__()
        self.__pd2d_background_phi = []
        self.__pd2d_background_2theta = []
        self.__pd2d_background_2theta_phi_intensity = [[]]
        self.phi = phi
        self.ttheta = ttheta
        self.intensity = intensity

    @property
    def phi(self):
        return tuple(self.__pd2d_background_phi)
    @phi.setter
    def phi(self, l_x):
        l_x_in = []
        for x in l_x:
            x_in = float(x)
            l_x_in.append(x_in)
        self.__pd2d_background_phi = l_x_in
        len_x = len(l_x_in)

        for ind, l_intensity in enumerate(self.intensity):
            len_1 = len(l_intensity)
            if len_1 > len_x:
                self.__pd2d_background_2theta_phi_intensity[ind] = self.__pd2d_background_2theta_phi_intensity[ind][:len_x]
            elif len_1 < len_x:
                l_fitable = [Fitable(value=0., name="_pd2d_background_2theta_phi_intensity") for hh in range(len_x-len_1)] #default
                self.__pd2d_background_2theta_phi_intensity[ind].extend(l_fitable)




    @property
    def ttheta(self):
        return tuple(self.__pd2d_background_2theta)
    @ttheta.setter
    def ttheta(self, l_x):
        l_x_in = []
        for x in l_x:
            x_in = float(x)
            l_x_in.append(x_in)
        self.__pd2d_background_2theta = l_x_in

        len_x = len(l_x_in)
        len_1 = len(self.intensity)
        if len_1 > len_x:
            self.__pd2d_background_2theta_phi_intensity = self.__pd2d_background_2theta_phi_intensity[:len_x]
        elif len_1 < len_x:
            l_fitable = [[] for hh in range(len_x-len_1)] #default
            self.__pd2d_background_2theta_phi_intensity.extend(l_fitable)        

    @property
    def intensity(self):
        return tuple(self.__pd2d_background_2theta_phi_intensity)
    @intensity.setter
    def intensity(self, ll_x):
        ll_fitable = []
        for l_x in ll_x:
            l_fitable = []
            for x in l_x:
                if isinstance(x, Fitable):
                    x_in = x
                else:
                    x_in = Fitable()
                    x_in.name = "_pd2d_background_2theta_phi_intensity"
                    flag = x_in.take_it(x)
                l_fitable.append(x_in)
            ll_fitable.append(l_fitable)
        len_x = len(ll_fitable)

        len_1 = len(self.ttheta)
        if len_1 < len_x:
            ll_fitable = ll_fitable[:len_1]
        elif len_1 > len_x:
            ll_fitable.extend([[] for hh in range(len_1-len_x)])

        len_1 = len(self.phi)
        for l_fitable in ll_fitable:
            len_x = len(l_fitable)
            if len_1 < len_x:
                l_fitable = l_fitable[:len_1]
            elif len_1 > len_x:
                l_fitable.extend([Fitable(value=0., name="_pd2d_background_2theta_phi_intensity") for hh in range(len_1-len_x)])

        self.__pd2d_background_2theta_phi_intensity = ll_fitable



            
    def __repr__(self):
        ls_out = ["Pd2dBackground:"]
        ls_out.append("{:12} ".format(len(self.phi)) + " ".join(["{:6.2f}      ".format(_) for _ in self.ttheta]))
        ll_intensity = self.intensity
        ll_intensity = [[ll_intensity[_2][_1] for _2 in range(len(ll_intensity))] for _1 in range(len(ll_intensity[0]))]
        for phi, l_intensity in zip(self.phi, ll_intensity):
            ls_out.append("{:12.2f} ".format(phi) + " ".join(["{:12}".format(_.print_with_sigma) for _ in l_intensity]))
        return "\n".join(ls_out)

    @property
    def to_cif(self):
        ls_out = []
        if self.is_defined:
            ls_out.append("_pd2d_background_2theta_phi_intensity")
            ls_out.append(";")
            ls_out.append("{:12} ".format(len(self.phi)) + " ".join(["{:6.2f}      ".format(_) for _ in self.ttheta]))
            ll_intensity = self.intensity
            ll_intensity = [[ll_intensity[_2][_1] for _2 in range(len(ll_intensity))] for _1 in range(len(ll_intensity[0]))]
            for phi, l_intensity in zip(self.phi, ll_intensity):
                ls_out.append("{:12.2f} ".format(phi) + " ".join(["{:12}".format(_.print_with_sigma) for _ in l_intensity]))
            ls_out.append(";")
        return "\n".join(ls_out)

    def from_cif(self, string: str):
        cif_global = Global()
        flag = cif_global.take_from_string(string)
        if not flag:
            return False
        flag = False
        flag = cif_global.is_prefix("_pd2d_background_2theta_phi_intensity")
        if flag:
            cif_value = cif_global["_pd2d_background_2theta_phi_intensity"]
            string = cif_value.value
            l_1 = string.strip().split("\n")
            l_ttheta = [float(_) for _ in l_1[0].strip().split()[1:]]
            l_phi, ll_intensity = [], []
            for line in l_1[1:]:
                l_1 = line.strip().split()
                l_phi.append(float(l_1[0]))
                ll_intensity.append(l_1[1:])
            ll_intensity = [[ll_intensity[_2][_1] for _2 in range(len(ll_intensity))] for _1 in range(len(ll_intensity[0]))]
            self.phi = l_phi
            self.ttheta = l_ttheta
            self.intensity = ll_intensity
        else:
            self.phi, self.ttheta, self.intensity = [], [], [[]]
        return True

    @property
    def is_defined(self):
        cond = all([self.phi is not None, self.ttheta is not None, self.intensity is not None])
        return cond

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

    def _show_message(self, s_out: str):
        print("***  Error ***")
        print(s_out)

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

