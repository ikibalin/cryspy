"""
define classe Pd2dExclude2Theta which describes the excluded regions in 2d powder diffraction experiment
"""
__author__ = 'ikibalin'
__version__ = "2019_09_10"
import os
import numpy


from pystar import CIFglobal


class Pd2dExclude2Theta(object):
    """
    Example:

    loop_
    _pd2d_exclude_phi_min
    _pd2d_exclude_phi_max
    _pd2d_exclude_tth_min
    _pd2d_exclude_tth_max
     -10.0 60.0 4.0 10.0 
    """
    def __init__(self, ttheta_min=[], ttheta_max=[], phi_min=[], phi_max=[]
                 ):
        super(Pd2dExclude2Theta, self).__init__()
        self.__pd2d_exclude_2theta_min = None
        self.__pd2d_exclude_2theta_max = None
        self.__pd2d_exclude_phi_min = None
        self.__pd2d_exclude_phi_max = None

        self.ttheta_min = ttheta_min
        self.ttheta_max = ttheta_max
        self.phi_min = phi_min
        self.phi_max = phi_max

    @property
    def ttheta_min(self):
        return self.__pd2d_exclude_2theta_min
    @ttheta_min.setter
    def ttheta_min(self, l_x):
        l_x_in = []
        for x in l_x:
            if isinstance(x, float):
                x_in = x
            else:
                x_in = float(x)
            l_x_in.append(x_in)
        np_x_in = numpy.array(l_x_in, dtype=float)
        self.__pd2d_exclude_2theta_min = np_x_in

    @property
    def ttheta_max(self):
        return self.__pd2d_exclude_2theta_max
    @ttheta_max.setter
    def ttheta_max(self, l_x):
        l_x_in = []
        for x in l_x:
            if isinstance(x, float):
                x_in = x
            else:
                x_in = float(x)
            l_x_in.append(x_in)
        np_x_in = numpy.array(l_x_in, dtype=float)
        self.__pd2d_exclude_2theta_max = np_x_in


    @property
    def phi_min(self):
        return self.__pd2d_exclude_phi_min
    @phi_min.setter
    def phi_min(self, l_x):
        l_x_in = []
        for x in l_x:
            if isinstance(x, float):
                x_in = x
            else:
                x_in = float(x)
            l_x_in.append(x_in)
        np_x_in = numpy.array(l_x_in, dtype=float)
        self.__pd2d_exclude_phi_min = np_x_in

    @property
    def phi_max(self):
        return self.__pd2d_exclude_phi_max
    @phi_max.setter
    def phi_max(self, l_x):
        l_x_in = []
        for x in l_x:
            if isinstance(x, float):
                x_in = x
            else:
                x_in = float(x)
            l_x_in.append(x_in)
        np_x_in = numpy.array(l_x_in, dtype=float)
        self.__pd2d_exclude_phi_max = np_x_in



    def __repr__(self):
        ls_out = ["Pd2dExclude2Theta:"]
        ls_out.append("ttheta_min ttheta_max    phi_min    phi_max")
        for _1, _2, _3, _4 in zip(self.ttheta_min, self.ttheta_max, self.phi_min, self.phi_max):
            ls_out.append("{:10.2f} {:10.2f} {:10.2f} {:10.2f}".format(_1, _2, _3, _4))
        return "\n".join(ls_out)

    @property
    def to_cif(self):
        ls_out = []
        if self.is_defined:
            ls_out.append("loop_")
            ls_out.append("_pd2d_exclude_ttheta_min")
            ls_out.append("_pd2d_exclude_ttheta_max")
            ls_out.append("_pd2d_exclude_phi_min")
            ls_out.append("_pd2d_exclude_phi_max")
            for _1, _2, _3, _4 in zip(self.ttheta_min, self.ttheta_max, self.phi_min, self.phi_max):
                ls_out.append("{:} {:} {:} {:}".format(_1, _2, _3, _4))
        return "\n".join(ls_out)

    def from_cif(self, string: str):
        cif_global = CIFglobal()
        flag = cif_global.take_from_string(string)
        if not flag:
            return False
        flag = False
        flag = cif_global.is_prefix("_pd2d_exclude")
        if flag:
            cif_loop = cif_global["_pd2d_exclude"]
            l_name = cif_loop.names
            if "_pd2d_exclude_2theta_min" in l_name:
                self.ttheta_min = [float(_1) for _1 in cif_loop["_pd2d_exclude_2theta_min"]]
            if "_pd2d_exclude_2theta_max" in l_name:
                self.ttheta_max = [float(_1) for _1 in cif_loop["_pd2d_exclude_2theta_max"]]
            if "_pd2d_exclude_phi_min" in l_name:
                self.phi_min = [float(_1) for _1 in cif_loop["_pd2d_exclude_phi_min"]]
            if "_pd2d_exclude_phi_max" in l_name:
                self.phi_max = [float(_1) for _1 in cif_loop["_pd2d_exclude_phi_max"]]
        else:
            self.ttheta_min, self.ttheta_max, self.phi_min, self.phi_max = [], [], [], []
        return True

    @property
    def is_defined(self):
        cond = all([self.ttheta_min is not None, self.ttheta_max is not None, self.phi_min is not None, self.phi_max is not None])
        return cond

    @property
    def is_variable(self):
        return False
    
    def get_variables(self):
        return []

    def _show_message(self, s_out: str):
        print("***  Error ***")
        print(s_out)
