"""
:Constants:

:Functions:

"""
import os
import numpy
import warnings
from fractions import Fraction
from typing import List, Tuple

F_FORMMAG = os.path.join(os.path.dirname(__file__), "formmag.tab")


def _load_handbook_m(self):
    f_name = os.path.join(self.__f_dir_prog, "tables", "formmag.tab")
    fid = open(f_name, 'r')
    lcont = fid.readlines()
    fid.close()
    lcont = [line for line in lcont if line.startswith("F")]
    ldcard = []
    for line in lcont:
        lhelp = line.strip().split()
        dcard = {"type_m": lhelp[1], "order": int(lhelp[2]),
                 "A": float(lhelp[3]),"a": float(lhelp[4]),
                 "B": float(lhelp[5]),"b": float(lhelp[6]),
                 "C": float(lhelp[7]),"c": float(lhelp[8]),
                 "D": float(lhelp[9])}
        ldcard.append(dcard)
    self.__handbook_mag = ldcard




def j0j2(self):
    flag = any([self.__j0_A is None, self.__j0_a is None, self.__j0_B is None, self.__j0_b is None, 
                self.__j0_C is None, self.__j0_c is None, self.__j0_D is None, 
                self.__j2_A is None, self.__j2_a is None, self.__j2_B is None, self.__j2_b is None, 
                self.__j2_C is None, self.__j2_c is None, self.__j2_D is None])
    if flag:
        l_j0_A, l_j0_a, l_j0_B, l_j0_b, l_j0_C, l_j0_c, l_j0_D = [], [], [], [], [], [], []
        l_j2_A, l_j2_a, l_j2_B, l_j2_b, l_j2_C, l_j2_c, l_j2_D = [], [], [], [], [], [], []
        for type_symbol in self.type_symbol:
            j0_A, j0_a, j0_B, j0_b, j0_C, j0_c, j0_D, j2_A, j2_a, j2_B, j2_b, j2_C, j2_c, j2_D = self._get_j0j2(type_symbol)
            l_j0_A.append(j0_A)
            l_j0_a.append(j0_a)
            l_j0_B.append(j0_B)
            l_j0_b.append(j0_b)
            l_j0_C.append(j0_C)
            l_j0_c.append(j0_c)
            l_j0_D.append(j0_D)
            l_j2_A.append(j2_A)
            l_j2_a.append(j2_a)
            l_j2_B.append(j2_B)
            l_j2_b.append(j2_b)
            l_j2_C.append(j2_C)
            l_j2_c.append(j2_c)
            l_j2_D.append(j2_D)
        self.__j0_A, self.__j0_a, self.__j0_B, self.__j0_b = tuple(l_j0_A), tuple(l_j0_a), tuple(l_j0_B), tuple(l_j0_b)
        self.__j0_C, self.__j0_c, self.__j0_D = tuple(l_j0_C), tuple(l_j0_c), tuple(l_j0_D)
        self.__j2_A, self.__j2_a, self.__j2_B, self.__j2_b = tuple(l_j2_A), tuple(l_j2_a), tuple(l_j2_B), tuple(l_j2_b)
        self.__j2_C, self.__j2_c, self.__j2_D = tuple(l_j2_C), tuple(l_j2_c), tuple(l_j2_D)
        l_res = tuple([l_j0_A, l_j0_a, l_j0_B, l_j0_b, l_j0_C, l_j0_c, l_j0_D, 
                       l_j2_A, l_j2_a, l_j2_B, l_j2_b, l_j2_C, l_j2_c, l_j2_D])
    else:
        j0_A, j0_a, j0_B, j0_b = self.__j0_A, self.__j0_a, self.__j0_B, self.__j0_b
        j0_C, j0_c, j0_D = self.__j0_C, self.__j0_c, self.__j0_D
        j2_A, j2_a, j2_B, j2_b = self.__j2_A, self.__j2_a, self.__j2_B, self.__j2_b
        j2_C, j2_c, j2_D = self.__j2_C, self.__j2_c, self.__j2_D
        l_res = tuple([j0_A, j0_a, j0_B, j0_b, j0_C, j0_c, j0_D, 
                       j2_A, j2_a, j2_B, j2_b, j2_C, j2_c, j2_D])
    return l_res
