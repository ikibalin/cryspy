"""
define class PdPhase which describes the phase contributions in 1d powder diffraction experiment
"""
__author__ = 'ikibalin'
__version__ = "2019_09_06"
import os
import numpy


from pycifstar import Global
from cryspy.f_common.cl_fitable import Fitable


class PdPhase(object):
    """
    Example:

    loop_
    _pd_phase_label
    _pd_phase_scale
    _pd_phase_igsize
     Fe3O4 0.02381 0.0
    """
    def __init__(self, label=[], scale=[], igsize=[]):
        super(PdPhase, self).__init__()
        self.__pd_phase_label = []
        self.__pd_phase_scale = []
        self.__pd_phase_igsize = []

        self.label = label
        self.scale = scale
        self.igsize = igsize

    @property
    def label(self):
        return self.__pd_phase_label
    @label.setter
    def label(self, l_x):
        l_x_in = []
        for x in l_x:
            x_in = str(x).strip()
            l_x_in.append(x_in)
        self.__pd_phase_label = l_x_in
        len_x = len(l_x_in)

        len_1 = len(self.scale)
        if len_1 > len_x:
            self.__pd_phase_scale = self.__pd_phase_scale[:len_x]
        elif len_1 < len_x:
            l_fitable = [Fitable(value=1., name="_pd_phase_scale") for hh in range(len_x-len_1)] #default
            self.__pd_phase_scale.extend(l_fitable)

        len_1 = len(self.igsize)
        if len_1 > len_x:
            self.__pd_phase_igsize = self.__pd_phase_igsize[:len_x]
        elif len_1 < len_x:
            l_fitable = [Fitable(value=0., name="_pd_phase_igsize") for hh in range(len_x-len_1)] #default
            self.__pd_phase_igsize.extend(l_fitable)

    @property
    def scale(self):
        return tuple(self.__pd_phase_scale)
    @scale.setter
    def scale(self, l_x):
        l_fitable = []
        for x in l_x:
            if isinstance(x, Fitable):
                x_in = x
            else:
                x_in = Fitable()
                x_in.name = "_pd_phase_scale"
                flag = x_in.take_it(x)
            l_fitable.append(x_in)
        len_x = len(l_fitable)
        len_1 = len(self.label)
        if len_1 < len_x:
            l_fitable = l_fitable[:len_1]
        elif len_1 > len_x:
            l_fitable.extend([Fitable(value=1., name="_pd_phase_scale") for hh in range(len_1-len_x)])
        self.__pd_phase_scale = l_fitable

    @property
    def igsize(self):
        return tuple(self.__pd_phase_igsize)
    @igsize.setter
    def igsize(self, l_x):
        l_fitable = []
        for x in l_x:
            if isinstance(x, Fitable):
                x_in = x
            else:
                x_in = Fitable()
                x_in.name = "_pd_phase_igsize"
                flag = x_in.take_it(x)
            l_fitable.append(x_in)
        len_x = len(l_fitable)
        len_1 = len(self.label)
        if len_1 < len_x:
            l_fitable = l_fitable[:len_1]
        elif len_1 > len_x:
            l_fitable.extend([Fitable(value=0., name="_pd_phase_igsize") for hh in range(len_1-len_x)])
        self.__pd_phase_igsize = l_fitable


            
    def __repr__(self):
        ls_out = ["PdPhase:"]
        ls_out.append(" label   scale        igsize      ")
        for _1, _2, _3 in zip(self.label, self.scale, self.igsize):
            ls_out.append(" {:7} {:12} {:12}".format(_1, _2.print_with_sigma, _3.print_with_sigma))
        return "\n".join(ls_out)

    @property
    def to_cif(self):
        ls_out = []
        if self.is_defined:
            ls_out.append("loop_")
            ls_out.append("_pd_phase_label")
            ls_out.append("_pd_phase_scale")
            ls_out.append("_pd_phase_igsize")
            for _1, _2, _3 in zip(self.label, self.scale, self.igsize):
                ls_out.append("{:} {:} {:}".format(_1, _2.print_with_sigma, _3.print_with_sigma))
        return "\n".join(ls_out)

    def from_cif(self, string: str):
        cif_global = Global()
        flag = cif_global.take_from_string(string)
        if not flag:
            return False
        flag = False
        flag = cif_global.is_prefix("_pd_phase")
        if flag:
            cif_loop = cif_global["_pd_phase"]
            l_name = cif_loop.names
            if "_pd_phase_label" in l_name:
                self.label = [_1.strip() for _1 in cif_loop["_pd_phase_label"]]
            if "_pd_phase_scale" in l_name:
                self.scale = [_1.strip() for _1 in cif_loop["_pd_phase_scale"]]
            if "_pd_phase_igsize" in l_name:
                self.igsize = [_1.strip() for _1 in cif_loop["_pd_phase_igsize"]]
        else:
            self.label = []
        return True

    @property
    def is_defined(self):
        cond = all([self.label is not None, self.scale is not None, self.igsize is not None])
        return cond

    @property
    def is_variable(self):
        res = (any([_1.refinement for _1 in self.scale]) |
               any([_1.refinement for _1 in self.igsize]))
        return res
    
    def get_variables(self):
        l_variable = []
        l_variable.extend([_1 for _1 in self.scale if _1.refinement])
        l_variable.extend([_1 for _1 in self.igsize if _1.refinement])
        return l_variable

    def _show_message(self, s_out: str):
        print("***  Error ***")
        print(s_out)
