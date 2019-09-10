"""
define class Refln
"""
__author__ = 'ikibalin'
__version__ = "2019_09_09"
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
     _refln_chi_11_A_calc
     _refln_chi_12_A_calc
     _refln_chi_13_A_calc
     _refln_chi_21_A_calc
     _refln_chi_22_A_calc
     _refln_chi_23_A_calc
     _refln_chi_31_A_calc
     _refln_chi_32_A_calc
     _refln_chi_33_A_calc
     _refln_chi_11_B_calc
     _refln_chi_12_B_calc
     _refln_chi_13_B_calc
     _refln_chi_21_B_calc
     _refln_chi_22_B_calc
     _refln_chi_23_B_calc
     _refln_chi_31_B_calc
     _refln_chi_32_B_calc
     _refln_chi_33_B_calc
    
    reference: https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Crefln.html
    """
    def __init__(self, h=[], k=[], l=[], f_nucl=[], sft_11=[], sft_12=[], sft_13=[], sft_21=[], sft_22=[], sft_23=[], sft_31=[], sft_32=[], sft_33=[]):
        super(Refln, self).__init__()
        self.__refln_index_h = None
        self.__refln_index_k = None
        self.__refln_index_l = None
        self.__refln_d_spacing = None
        self.__refln_f_nucl_calc = None
        self.__refln_sft_11_calc = None
        self.__refln_sft_12_calc = None
        self.__refln_sft_13_calc = None
        self.__refln_sft_21_calc = None
        self.__refln_sft_22_calc = None
        self.__refln_sft_23_calc = None
        self.__refln_sft_31_calc = None
        self.__refln_sft_32_calc = None
        self.__refln_sft_33_calc = None
        self.h = h
        self.k = k
        self.l = l
        self.f_nucl = f_nucl
        self.sft_11 = sft_11
        self.sft_12 = sft_12
        self.sft_13 = sft_13
        self.sft_21 = sft_21
        self.sft_22 = sft_22
        self.sft_23 = sft_23
        self.sft_31 = sft_31
        self.sft_32 = sft_32
        self.sft_33 = sft_33

    @property
    def h(self):
        return self.__refln_index_h
    @h.setter
    def h(self, l_x):
        l_x_in = []
        for x in l_x:
            if isinstance(x, int):
                x_in = x
            else:
                x_in = int(x)
            l_x_in.append(x_in)
        np_x_in = numpy.array(l_x_in, dtype=int)
        self.__refln_index_h = np_x_in

    @property
    def k(self):
        return self.__refln_index_k
    @k.setter
    def k(self, l_x):
        l_x_in = []
        for x in l_x:
            if isinstance(x, int):
                x_in = x
            else:
                x_in = int(x)
            l_x_in.append(x_in)
        np_x_in = numpy.array(l_x_in, dtype=int)
        self.__refln_index_k = np_x_in

    @property
    def l(self):
        return self.__refln_index_l
    @l.setter
    def l(self, l_x):
        l_x_in = []
        for x in l_x:
            if isinstance(x, int):
                x_in = x
            else:
                x_in = int(x)
            l_x_in.append(x_in)
        np_x_in = numpy.array(l_x_in, dtype=int)
        self.__refln_index_l = np_x_in

    @property
    def f_nucl(self):
        return self.__refln_f_nucl_calc
    @f_nucl.setter
    def f_nucl(self, l_x):
        l_x_in = []
        for x in l_x:
            if isinstance(x, complex):
                x_in = x
            else:
                x_in = complex(x)
            l_x_in.append(x_in)
        np_x_in = numpy.array(l_x_in, dtype=complex)
        self.__refln_f_nucl_calc = np_x_in

    @property
    def sft_11(self):
        return self.__refln_sft_11_calc
    @sft_11.setter
    def sft_11(self, l_x):
        l_x_in = []
        for x in l_x:
            if isinstance(x, complex):
                x_in = x
            else:
                x_in = complex(x)
            l_x_in.append(x_in)
        np_x_in = numpy.array(l_x_in, dtype=complex)
        self.__refln_sft_11_calc = np_x_in

    @property
    def sft_12(self):
        return self.__refln_sft_12_calc
    @sft_12.setter
    def sft_12(self, l_x):
        l_x_in = []
        for x in l_x:
            if isinstance(x, complex):
                x_in = x
            else:
                x_in = complex(x)
            l_x_in.append(x_in)
        np_x_in = numpy.array(l_x_in, dtype=complex)
        self.__refln_sft_12_calc = np_x_in

    @property
    def sft_13(self):
        return self.__refln_sft_13_calc
    @sft_13.setter
    def sft_13(self, l_x):
        l_x_in = []
        for x in l_x:
            if isinstance(x, complex):
                x_in = x
            else:
                x_in = complex(x)
            l_x_in.append(x_in)
        np_x_in = numpy.array(l_x_in, dtype=complex)
        self.__refln_sft_13_calc = np_x_in

    @property
    def sft_21(self):
        return self.__refln_sft_21_calc
    @sft_21.setter
    def sft_21(self, l_x):
        l_x_in = []
        for x in l_x:
            if isinstance(x, complex):
                x_in = x
            else:
                x_in = complex(x)
            l_x_in.append(x_in)
        np_x_in = numpy.array(l_x_in, dtype=complex)
        self.__refln_sft_21_calc = np_x_in

    @property
    def sft_22(self):
        return self.__refln_sft_22_calc
    @sft_22.setter
    def sft_22(self, l_x):
        l_x_in = []
        for x in l_x:
            if isinstance(x, complex):
                x_in = x
            else:
                x_in = complex(x)
            l_x_in.append(x_in)
        np_x_in = numpy.array(l_x_in, dtype=complex)
        self.__refln_sft_22_calc = np_x_in

    @property
    def sft_23(self):
        return self.__refln_sft_23_calc
    @sft_23.setter
    def sft_23(self, l_x):
        l_x_in = []
        for x in l_x:
            if isinstance(x, complex):
                x_in = x
            else:
                x_in = complex(x)
            l_x_in.append(x_in)
        np_x_in = numpy.array(l_x_in, dtype=complex)
        self.__refln_sft_23_calc = np_x_in

    @property
    def sft_31(self):
        return self.__refln_sft_31_calc
    @sft_31.setter
    def sft_31(self, l_x):
        l_x_in = []
        for x in l_x:
            if isinstance(x, complex):
                x_in = x
            else:
                x_in = complex(x)
            l_x_in.append(x_in)
        np_x_in = numpy.array(l_x_in, dtype=complex)
        self.__refln_sft_31_calc = np_x_in

    @property
    def sft_32(self):
        return self.__refln_sft_32_calc
    @sft_32.setter
    def sft_32(self, l_x):
        l_x_in = []
        for x in l_x:
            if isinstance(x, complex):
                x_in = x
            else:
                x_in = complex(x)
            l_x_in.append(x_in)
        np_x_in = numpy.array(l_x_in, dtype=complex)
        self.__refln_sft_32_calc = np_x_in

    @property
    def sft_33(self):
        return self.__refln_sft_33_calc
    @sft_33.setter
    def sft_33(self, l_x):
        l_x_in = []
        for x in l_x:
            if isinstance(x, complex):
                x_in = x
            else:
                x_in = complex(x)
            l_x_in.append(x_in)
        np_x_in = numpy.array(l_x_in, dtype=complex)
        self.__refln_sft_33_calc = np_x_in

            
    def __repr__(self):
        ls_out = ["Refln:"]
        ls_out.append("   h   k   l           A           B    sft_11_A    sft_12_A    sft_13_A      sft_11_B    sft_12_B    sft_13_B")
        ls_out.append("                                        sft_21_A    sft_22_A    sft_23_A      sft_21_B    sft_22_B    sft_23_B")
        ls_out.append("                                        sft_31_A    sft_32_A    sft_33_A      sft_31_B    sft_32_B    sft_33_B")
        for _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13 in zip(self.h, self.k, self.l, self.f_nucl, 
                                                                          self.sft_11, self.sft_12, self.sft_13, 
                                                                          self.sft_21, self.sft_22, self.sft_23, 
                                                                          self.sft_31, self.sft_32, self.sft_33):
            ls_out.append("\n{:4}{:4}{:4}{:12.5f}{:12.5f}{:12.5f}{:12.5f}{:12.5f}  {:12.5f}{:12.5f}{:12.5f}".format(
                _1, _2, _3, _4.real, _4.imag, _5.real, _6.real, _7.real, _5.imag, _6.imag, _7.imag))
            ls_out.append("                                    {:12.5f}{:12.5f}{:12.5f}  {:12.5f}{:12.5f}{:12.5f}".format(
                _8.real, _9.real, _10.real, _8.imag, _9.imag, _10.imag))
            ls_out.append("                                    {:12.5f}{:12.5f}{:12.5f}  {:12.5f}{:12.5f}{:12.5f}".format(
                _11.real, _12.real, _13.real, _11.imag, _12.imag, _13.imag))
        
        ls_out.append("*  all structure factors are in 10**-12cm")
        ls_out.append("** structure factor tensor given in crystallographic Cartesian coordinate system\n")
        return "\n".join(ls_out)

    @property
    def to_cif(self):
        ls_out = []
        if self.is_defined:
            ls_out.append("loop_")
            ls_out.append("_refln_index_h")
            ls_out.append("_refln_index_k")
            ls_out.append("_refln_index_l")
            #ls_out.append("_refln_d_spacing")
            ls_out.append("_refln_A_calc")
            ls_out.append("_refln_B_calc")
            ls_out.append("_refln_chi_11_A_calc")
            ls_out.append("_refln_chi_12_A_calc")
            ls_out.append("_refln_chi_13_A_calc")
            ls_out.append("_refln_chi_21_A_calc")
            ls_out.append("_refln_chi_22_A_calc")
            ls_out.append("_refln_chi_23_A_calc")
            ls_out.append("_refln_chi_31_A_calc")
            ls_out.append("_refln_chi_32_A_calc")
            ls_out.append("_refln_chi_33_A_calc")
            ls_out.append("_refln_chi_11_B_calc")
            ls_out.append("_refln_chi_12_B_calc")
            ls_out.append("_refln_chi_13_B_calc")
            ls_out.append("_refln_chi_21_B_calc")
            ls_out.append("_refln_chi_22_B_calc")
            ls_out.append("_refln_chi_23_B_calc")
            ls_out.append("_refln_chi_31_B_calc")
            ls_out.append("_refln_chi_32_B_calc")
            ls_out.append("_refln_chi_33_B_calc")
            for _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13 in zip(self.h, self.k, self.l, self.f_nucl, 
                                                                          self.sft_11, self.sft_12, self.sft_13, 
                                                                          self.sft_21, self.sft_22, self.sft_23, 
                                                                          self.sft_31, self.sft_32, self.sft_33):
                ls_out.append("{:} {:} {:} {:.5f} {:.5f} {:.5f} {:.5f} {:.5f} {:.5f} {:.5f} {:.5f} {:.5f} {:.5f} {:.5f} {:.5f} {:.5f} {:.5f} {:.5f} {:.5f} {:.5f} {:.5f} {:.5f} {:.5f}".format(
                _1, _2, _3, _4.real, _4.imag, _5.real, _6.real, _7.real, _8.real, _9.real, _10.real, _11.real, _12.real, _13.real,
                _5.imag, _6.imag, _7.imag, _8.imag, _9.imag, _10.imag, _11.imag, _12.imag, _13.imag))
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

            l_a, l_b, l_chi_11_a, l_chi_12_a, l_chi_13_a, l_chi_21_a, l_chi_22_a, l_chi_23_a = [], [], [], [], [], [], [], []
            l_chi_31_a, l_chi_32_a, l_chi_33_a = [], [], []
            l_chi_11_b, l_chi_12_b, l_chi_13_b, l_chi_21_b, l_chi_22_b, l_chi_23_b = [], [], [], [], [], []
            l_chi_31_b, l_chi_32_b, l_chi_33_b = [], [], []
            if "_refln_index_h" in l_name:
                self.h = [int(_1) for _1 in cif_loop["_refln_index_h"]]
            if "_refln_index_k" in l_name:
                self.k = [int(_1) for _1 in cif_loop["_refln_index_k"]]
            if "_refln_index_l" in l_name:
                self.l = [int(_1) for _1 in cif_loop["_refln_index_l"]]
            if "_refln_A_calc" in l_name:
                l_a = [float(_1) for _1 in cif_loop["_refln_A_calc"]]
            if "_refln_B_calc" in l_name:
                l_b = [float(_1) for _1 in cif_loop["_refln_B_calc"]]
            if "_refln_chi_11_A_calc" in l_name:
                l_chi_11_a = [float(_1) for _1 in cif_loop["_refln_chi_11_A_calc"]]
            if "_refln_chi_12_A_calc" in l_name:
                l_chi_12_a = [float(_1) for _1 in cif_loop["_refln_chi_12_A_calc"]]
            if "_refln_chi_13_A_calc" in l_name:
                l_chi_13_a = [float(_1) for _1 in cif_loop["_refln_chi_13_A_calc"]]
            if "_refln_chi_21_A_calc" in l_name:
                l_chi_21_a = [float(_1) for _1 in cif_loop["_refln_chi_21_A_calc"]]
            if "_refln_chi_22_A_calc" in l_name:
                l_chi_22_a = [float(_1) for _1 in cif_loop["_refln_chi_22_A_calc"]]
            if "_refln_chi_23_A_calc" in l_name:
                l_chi_23_a = [float(_1) for _1 in cif_loop["_refln_chi_23_A_calc"]]
            if "_refln_chi_31_A_calc" in l_name:
                l_chi_31_a = [float(_1) for _1 in cif_loop["_refln_chi_31_A_calc"]]
            if "_refln_chi_32_A_calc" in l_name:
                l_chi_32_a = [float(_1) for _1 in cif_loop["_refln_chi_32_A_calc"]]
            if "_refln_chi_33_A_calc" in l_name:
                l_chi_33_a = [float(_1) for _1 in cif_loop["_refln_chi_33_A_calc"]]
            if "_refln_chi_11_B_calc" in l_name:
                l_chi_11_b = [float(_1) for _1 in cif_loop["_refln_chi_11_B_calc"]]
            if "_refln_chi_12_B_calc" in l_name:
                l_chi_12_b = [float(_1) for _1 in cif_loop["_refln_chi_12_B_calc"]]
            if "_refln_chi_13_B_calc" in l_name:
                l_chi_13_b = [float(_1) for _1 in cif_loop["_refln_chi_13_B_calc"]]
            if "_refln_chi_21_B_calc" in l_name:
                l_chi_21_b = [float(_1) for _1 in cif_loop["_refln_chi_21_B_calc"]]
            if "_refln_chi_22_B_calc" in l_name:
                l_chi_22_b = [float(_1) for _1 in cif_loop["_refln_chi_22_B_calc"]]
            if "_refln_chi_23_B_calc" in l_name:
                l_chi_23_b = [float(_1) for _1 in cif_loop["_refln_chi_23_B_calc"]]
            if "_refln_chi_31_B_calc" in l_name:
                l_chi_31_b = [float(_1) for _1 in cif_loop["_refln_chi_31_B_calc"]]
            if "_refln_chi_32_B_calc" in l_name:
                l_chi_32_b = [float(_1) for _1 in cif_loop["_refln_chi_32_B_calc"]]
            if "_refln_chi_33_B_calc" in l_name:
                l_chi_33_b = [float(_1) for _1 in cif_loop["_refln_chi_33_B_calc"]]
            self.f_nucl = [complex(_1, _2) for _1, _2 in zip(l_a, l_b)]
            self.chi_11 = [complex(_1, _2) for _1, _2 in zip(l_chi_11_a, l_chi_11_b)]
            self.chi_12 = [complex(_1, _2) for _1, _2 in zip(l_chi_12_a, l_chi_12_b)]
            self.chi_13 = [complex(_1, _2) for _1, _2 in zip(l_chi_13_a, l_chi_13_b)]
            self.chi_21 = [complex(_1, _2) for _1, _2 in zip(l_chi_21_a, l_chi_21_b)]
            self.chi_22 = [complex(_1, _2) for _1, _2 in zip(l_chi_22_a, l_chi_22_b)]
            self.chi_23 = [complex(_1, _2) for _1, _2 in zip(l_chi_23_a, l_chi_23_b)]
            self.chi_31 = [complex(_1, _2) for _1, _2 in zip(l_chi_31_a, l_chi_31_b)]
            self.chi_32 = [complex(_1, _2) for _1, _2 in zip(l_chi_32_a, l_chi_32_b)]
            self.chi_33 = [complex(_1, _2) for _1, _2 in zip(l_chi_33_a, l_chi_33_b)]

        else:
            self.h, self.k, self.l, self.f_nucl = [], [], [], []
            self.chi_11, self.chi_12, self.chi_13 = [], [], []
            self.chi_21, self.chi_22, self.chi_23 = [], [], []
            self.chi_31, self.chi_32, self.chi_33 = [], [], []
        return True

    @property
    def is_defined(self):
        cond = all([self.h is not None, self.k is not None, self.l is not None, self.f_nucl is not None, 
                    self.sft_11 is not None, self.sft_12 is not None, self.sft_13 is not None,
                    self.sft_21 is not None, self.sft_22 is not None, self.sft_23 is not None,
                    self.sft_31 is not None, self.sft_32 is not None, self.sft_33 is not None])
        return cond

    @property
    def is_variable(self):
        return False
    
    def get_variables(self):
        return []

    def _show_message(self, s_out: str):
        print("***  Error ***")
        print(s_out)
