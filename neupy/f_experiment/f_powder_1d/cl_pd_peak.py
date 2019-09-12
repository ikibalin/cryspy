"""
define class PdPeak
"""
__author__ = 'ikibalin'
__version__ = "2019_09_06"
import os
import numpy


from pystar import Global


class PdPeak(object):
    """
    This section contains peak information extracted from the
    measured or, if present, the processed diffractogram. Each
    peak in this table will have a unique label (see _pd_peak_id).
    The reflections and phases associated with each peak will be
    specified in other sections (see the _pd_refln_ and
    _pd_phase_ sections).

    Note that peak positions are customarily determined from the
    processed diffractogram and thus corrections for position
    and intensity will have been previously applied.

    Example:

    loop_
    _pd_peak_index_h
    _pd_peak_index_k
    _pd_peak_index_l
    _pd_peak_index_mult
    _pd_peak_2theta
    _pd_peak_intensity_up
    _pd_peak_intensity_down
    _pd_peak_width_2theta

    reference: https://www.iucr.org/__data/iucr/cifdic_html/1/cif_pd.dic/Cpd_peak.html
    """
    def __init__(self, h=[], k=[], l=[], mult=[], ttheta=[], up=[], down=[], width_2theta=[]):
        super(PdPeak, self).__init__()
        self.__pd_peak_index_h = None
        self.__pd_peak_index_k = None
        self.__pd_peak_index_l = None
        self.__pd_peak_index_mult = None
        self.__pd_peak_2theta = None
        self.__pd_peak_intensity_up = None
        self.__pd_peak_intensity_down = None
        self.__pd_peak_width_2theta = None

        self.h = h
        self.k = k
        self.l = l
        self.mult = mult
        self.ttheta = ttheta
        self.up = up
        self.down = down
        self.width_2theta = width_2theta

    @property
    def h(self):
        return self.__pd_peak_index_h
    @h.setter
    def h(self, l_x):
        l_x_in = []
        for x in l_x:
            if isinstance(x, float):
                x_in = x
            else:
                x_in = int(round(x))
            l_x_in.append(x_in)
        np_x_in = numpy.array(l_x_in, dtype=int)
        self.__pd_peak_index_h = np_x_in

    @property
    def k(self):
        return self.__pd_peak_index_k
    @k.setter
    def k(self, l_x):
        l_x_in = []
        for x in l_x:
            if isinstance(x, float):
                x_in = x
            else:
                x_in = int(round(x))
            l_x_in.append(x_in)
        np_x_in = numpy.array(l_x_in, dtype=int)
        self.__pd_peak_index_k = np_x_in

    @property
    def l(self):
        return self.__pd_peak_index_l
    @l.setter
    def l(self, l_x):
        l_x_in = []
        for x in l_x:
            if isinstance(x, float):
                x_in = x
            else:
                x_in = int(round(x))
            l_x_in.append(x_in)
        np_x_in = numpy.array(l_x_in, dtype=int)
        self.__pd_peak_index_l = np_x_in

    @property
    def mult(self):
        return self.__pd_peak_index_mult
    @mult.setter
    def mult(self, l_x):
        l_x_in = []
        for x in l_x:
            if isinstance(x, float):
                x_in = x
            else:
                x_in = int(round(x))
            l_x_in.append(x_in)
        np_x_in = numpy.array(l_x_in, dtype=int)
        self.__pd_peak_index_mult = np_x_in

    @property
    def ttheta(self):
        return self.__pd_peak_2theta
    @ttheta.setter
    def ttheta(self, l_x):
        l_x_in = []
        for x in l_x:
            if isinstance(x, float):
                x_in = x
            else:
                x_in = float(x)
            l_x_in.append(x_in)
        np_x_in = numpy.array(l_x_in, dtype=float)
        self.__pd_peak_2theta = np_x_in


    @property
    def up(self):
        return self.__pd_peak_intensity_up
    @up.setter
    def up(self, l_x):
        l_x_in = []
        for x in l_x:
            if isinstance(x, float):
                x_in = x
            else:
                x_in = float(x)
            l_x_in.append(x_in)
        np_x_in = numpy.array(l_x_in, dtype=float)
        self.__pd_peak_intensity_up = np_x_in

    @property
    def down(self):
        return self.__pd_peak_intensity_down
    @down.setter
    def down(self, l_x):
        l_x_in = []
        for x in l_x:
            if isinstance(x, float):
                x_in = x
            else:
                x_in = float(x)
            l_x_in.append(x_in)
        np_x_in = numpy.array(l_x_in, dtype=float)
        self.__pd_peak_intensity_down = np_x_in

    @property
    def width_2theta(self):
        return self.__pd_peak_width_2theta
    @width_2theta.setter
    def width_2theta(self, l_x):
        l_x_in = []
        for x in l_x:
            if isinstance(x, float):
                x_in = x
            else:
                x_in = float(x)
            l_x_in.append(x_in)
        np_x_in = numpy.array(l_x_in, dtype=float)
        self.__pd_peak_width_2theta = np_x_in

    def __repr__(self):
        ls_out = ["PdPeak:"]
        ls_out.append("    h     k     l  mult    ttheta           up         down width_2theta")
        for _1, _2, _3, _4, _5, _6, _7, _8 in zip(self.h, self.k, self.l, self.mult, self.ttheta, self.up, self.down, self.width_2theta):
            ls_out.append("{:5} {:5} {:5} {:5} {:9.2f} {:12.5f} {:12.5f} {:12.5f}".format(_1, _2, _3, _4, _5, _6, _7, _8))
        return "\n".join(ls_out)

    @property
    def to_cif(self):
        ls_out = []
        if self.is_defined:
            ls_out.append("loop_")
            ls_out.append("_pd_peak_index_h")
            ls_out.append("_pd_peak_index_k")
            ls_out.append("_pd_peak_index_l")
            ls_out.append("_pd_peak_mult")
            ls_out.append("_pd_peak_2theta")
            ls_out.append("_pd_peak_intensity_up")
            ls_out.append("_pd_peak_intensity_down")
            ls_out.append("_pd_peak_width_2theta")
            for _1, _2, _3, _4, _5, _6, _7, _8 in zip(self.h, self.k, self.l, self.mult, self.ttheta, self.up, self.down, self.width_2theta):
                ls_out.append("{:} {:} {:} {:} {:.3f} {:.5f} {:.5f} {:.3f}".format(_1, _2, _3, _4, _5, _6, _7, _8))
        return "\n".join(ls_out)

    def from_cif(self, string: str):
        cif_global = Global()
        flag = cif_global.take_from_string(string)
        if not flag:
            return False
        flag = False
        flag = cif_global.is_prefix("_pd_peak")
        if flag:
            cif_loop = cif_global["_pd_peak"]
            l_name = cif_loop.names
            if "_pd_peak_index_h" in l_name:
                self.h = [int(_1) for _1 in cif_loop["_pd_peak_index_h"]]
            if "_pd_peak_index_k" in l_name:
                self.k = [int(_1) for _1 in cif_loop["_pd_peak_index_k"]]
            if "_pd_peak_index_l" in l_name:
                self.l = [int(_1) for _1 in cif_loop["_pd_peak_index_l"]]
            if "_pd_peak_index_mult" in l_name:
                self.mult = [int(_1) for _1 in cif_loop["_pd_peak_index_mult"]]
            if "_pd_peak_2theta" in l_name:
                self.ttheta = [float(_1) for _1 in cif_loop["_pd_peak_2theta"]]
            if "_pd_peak_intensity_up" in l_name:
                self.up = [float(_1) for _1 in cif_loop["_pd_peak_intensity_up"]]
            if "_pd_peak_intensity_down" in l_name:
                self.down = [float(_1) for _1 in cif_loop["_pd_peak_intensity_down"]]
            if "_pd_peak_width_2theta" in l_name:
                self.width_2theta = [float(_1) for _1 in cif_loop["_pd_peak_width_2theta"]]
        else:
            self.h, self.k, self.l, self.mult, self.ttheta, self.up, self.down, self.width_2theta = [], [], [], [], [], [], [], []
        return True

    @property
    def is_defined(self):
        cond = all([self.h is not None, self.k is not None, self.l is not None, self.mult is not None, 
                    self.ttheta is not None, self.up is not None, self.down is not None,
                    self.width_2theta is not None])
        return cond

    @property
    def is_variable(self):
        return False
    
    def get_variables(self):
        return []

    def _show_message(self, s_out: str):
        print("***  Error ***")
        print(s_out)
