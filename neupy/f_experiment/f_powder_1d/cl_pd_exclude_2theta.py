"""
define classe PdExclude2Theta which describes the excluded regions in 1d powder diffraction experiment
"""
__author__ = 'ikibalin'
__version__ = "2019_09_06"
import os
import numpy


from pystar import CIFglobal


class PdExclude2Theta(object):
    """
    Example:

    loop_
    _pd_exclude_2theta_min
    _pd_exclude_2theta_max
     4.0  12.0 
     30.0 45.0 
     58.0 63.0 
    """
    def __init__(self, min=[], max=[]
                 ):
        super(PdExclude2Theta, self).__init__()
        self.__pd_exclude_2theta_min = None
        self.__pd_exclude_2theta_max = None

        self.min = min
        self.max = max

    @property
    def min(self):
        return self.__pd_exclude_2theta_min
    @min.setter
    def min(self, l_x):
        l_x_in = []
        for x in l_x:
            if isinstance(x, float):
                x_in = x
            else:
                x_in = float(x)
            l_x_in.append(x_in)
        np_x_in = numpy.array(l_x_in, dtype=float)
        self.__pd_exclude_2theta_min = np_x_in


    @property
    def max(self):
        return self.__pd_exclude_2theta_max
    @max.setter
    def max(self, l_x):
        l_x_in = []
        for x in l_x:
            if isinstance(x, float):
                x_in = x
            else:
                x_in = float(x)
            l_x_in.append(x_in)
        np_x_in = numpy.array(l_x_in, dtype=float)
        self.__pd_exclude_2theta_max = np_x_in


    def __repr__(self):
        ls_out = ["PdExclude2Theta:"]
        ls_out.append("    min     max")
        for _1, _2 in zip(self.min, self.max):
            ls_out.append("{:7.2f} {:7.2f}".format(_1, _2))
        return "\n".join(ls_out)

    @property
    def to_cif(self):
        ls_out = []
        if self.is_defined:
            ls_out.append("loop_")
            ls_out.append("_pd_exclude_2theta_min")
            ls_out.append("_pd_exclude_2theta_max")
            for _1, _2 in zip(self.min, self.max):
                ls_out.append("{:} {:}".format(_1, _2))
        return "\n".join(ls_out)

    def from_cif(self, string: str):
        cif_global = CIFglobal()
        flag = cif_global.take_from_string(string)
        if not flag:
            return False
        flag = False
        flag = cif_global.is_prefix("_pd_exclude_2theta")
        if flag:
            cif_loop = cif_global["_pd_exclude_2theta"]
            l_name = cif_loop.names
            if "_pd_exclude_2theta_min" in l_name:
                self.min = [float(_1) for _1 in cif_loop["_pd_exclude_2theta_min"]]
            if "_pd_exclude_2theta_max" in l_name:
                self.max = [float(_1) for _1 in cif_loop["_pd_exclude_2theta_max"]]
        else:
            self.min, self.max = [], []
        return True

    @property
    def is_defined(self):
        cond = all([self.min is not None, self.max is not None])
        return cond

    @property
    def is_variable(self):
        return False
    
    def get_variables(self):
        return []

    def _show_message(self, s_out: str):
        print("***  Error ***")
        print(s_out)
