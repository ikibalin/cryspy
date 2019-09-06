"""
define class Refln
"""
__author__ = 'ikibalin'
__version__ = "2019_09_06"
import os
import numpy


from pystar import CIFglobal


class Refln(object):
    """
    Data items in the REFLN category record details about the
    reflections used to determine the ATOM_SITE data items.

    The REFLN data items refer to individual reflections and must
    be included in looped lists.

    The REFLNS data items specify the parameters that apply to all
    reflections. The REFLNS data items are not looped.


    Example:

    loop_
     _refln_index_h
     _refln_index_k
     _refln_index_l
     _refln_d_spacing
     _refln_A_calc
     _refln_B_calc
     _refln_F_squared_calc
     _refln_chi_11_A_calc
     _refln_chi_22_A_calc
     _refln_chi_33_A_calc
     _refln_chi_12_A_calc
     _refln_chi_13_A_calc
     _refln_chi_23_A_calc
     _refln_chi_11_B_calc
     _refln_chi_22_B_calc
     _refln_chi_33_B_calc
     _refln_chi_12_B_calc
     _refln_chi_13_B_calc
     _refln_chi_23_B_calc
    
    reference: https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Crefln.html
    """
    NOT FINISHED
    def __init__(self, ttheta=[], up=[], down=[], up_sigma=[], down_sigma=[]
                 ):
        super(Refln, self).__init__()
        self._refln_index_h = None
        self._refln_index_k = None
        self._refln_index_l = None
        self.__pd_meas_intensity_down = None
        self.__pd_meas_intensity_down_sigma = None

        self.ttheta = ttheta
        self.up = up
        self.up_sigma = up_sigma
        self.down = down
        self.down_sigma = down_sigma

    @property
    def ttheta(self):
        return self.__pd_meas_angle_2theta
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
        self.__pd_meas_angle_2theta = np_x_in


    @property
    def up(self):
        return self.__pd_meas_intensity_up
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
        self.__pd_meas_intensity_up = np_x_in

    @property
    def up_sigma(self):
        return self.__pd_meas_intensity_up_sigma
    @up_sigma.setter
    def up_sigma(self, l_x):
        l_x_in = []
        for x in l_x:
            if isinstance(x, float):
                x_in = x
            else:
                x_in = float(x)
            l_x_in.append(x_in)
        np_x_in = numpy.array(l_x_in, dtype=float)
        self.__pd_meas_intensity_up_sigma = np_x_in

    @property
    def down(self):
        return self.__pd_meas_intensity_down
    @down.setter
    def updown(self, l_x):
        l_x_in = []
        for x in l_x:
            if isinstance(x, float):
                x_in = x
            else:
                x_in = float(x)
            l_x_in.append(x_in)
        np_x_in = numpy.array(l_x_in, dtype=float)
        self.__pd_meas_intensity_down = np_x_in

    @property
    def down_sigma(self):
        return self.__pd_meas_intensity_down_sigma
    @down_sigma.setter
    def down_sigma(self, l_x):
        l_x_in = []
        for x in l_x:
            if isinstance(x, float):
                x_in = x
            else:
                x_in = float(x)
            l_x_in.append(x_in)
        np_x_in = numpy.array(l_x_in, dtype=float)
        self.__pd_meas_intensity_down_sigma = np_x_in


            
    def __repr__(self):
        ls_out = ["PdMeas:"]
        ls_out.append(" ttheta           up     up_sigma         down   down_sigma")
        for _1, _2, _3, _4, _5 in zip(self.ttheta, self.up, self.up_sigma, self.down, self.down_sigma):
            ls_out.append("{:7.2f} {:12.5f} {:12.5f} {:12.5f} {:12.5f}".format(_1, _2, _3, _4, _5))
        return "\n".join(ls_out)

    @property
    def to_cif(self):
        ls_out = []
        flag_calc = self.up_calc is not None
        if self.is_defined:
            ls_out.append("loop_")
            ls_out.append("_pd_meas_angle_2theta")
            ls_out.append("_pd_meas_intensity_up")
            ls_out.append("_pd_meas_intensity_up_sigma")
            ls_out.append("_pd_meas_intensity_down")
            ls_out.append("_pd_meas_intensity_down_sigma")
            for _1, _2, _3, _4, _5 in zip(self.ttheta, self.up, self.up_sigma, self.down, self.down_sigma):
                ls_out.append("{:} {:} {:} {:} {:}".format(_1, _2, _3, _4, _5))
        return "\n".join(ls_out)

    def from_cif(self, string: str):
        cif_global = CIFglobal()
        flag = cif_global.take_from_string(string)
        if not flag:
            return False
        flag = False
        flag = cif_global.is_prefix("_pd_meas")
        if flag:
            cif_loop = cif_global["_pd_meas"]
            l_name = cif_loop.names
            if "_pd_meas_angle_2theta" in l_name:
                self.ttheta = [float(_1) for _1 in cif_loop["_pd_meas_angle_2theta"]]
            if "_pd_meas_intensity_up" in l_name:
                self.up = [float(_1) for _1 in cif_loop["_pd_meas_intensity_up"]]
            if "_pd_meas_intensity_up_sigma" in l_name:
                self.up_sigma = [float(_1) for _1 in cif_loop["_pd_meas_intensity_up_sigma"]]
            if "_pd_meas_intensity_down" in l_name:
                self.down = [float(_1) for _1 in cif_loop["_pd_meas_intensity_down"]]
            if "_pd_meas_intensity_down_sigma" in l_name:
                self.down_sigma = [float(_1) for _1 in cif_loop["_pd_meas_intensity_down_sigma"]]
        else:
            self.ttheta, self.up, self.up_sigma, self.down, self.down_sigma = [], [], [], [], []
        return True

    @property
    def is_defined(self):
        cond = all([self.ttheta is not None, self.up is not None, self.up_sigma is not None, self.down is not None, self.down_sigma is not None])
        return cond

    @property
    def is_variable(self):
        return False
    
    def get_variables(self):
        return []

    def _show_message(self, s_out: str):
        print("***  Error ***")
        print(s_out)
