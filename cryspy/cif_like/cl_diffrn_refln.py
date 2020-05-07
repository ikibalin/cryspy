__author__ = 'ikibalin'
__version__ = "2020_01_02"

import os
import math
import numpy
from pycifstar import Global

import warnings
from typing import List, Tuple
from cryspy.common.cl_item_constr import ItemConstr
from cryspy.common.cl_loop_constr import LoopConstr


class DiffrnRefln(ItemConstr):
    """
Data items in the DIFFRN_REFLN category record details about
the intensities measured in the diffraction experiment.

The DIFFRN_REFLN data items refer to individual intensity
measurements and must be included in looped lists.

Description in cif file::

 _diffrn_refln_index_h     0 
 _diffrn_refln_index_k     0
 _diffrn_refln_index_l     8
 _diffrn_refln_fr          0.64545
 _diffrn_refln_fr_sigma    0.01329

`<https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Cdiffrn_refln.html>`_
    """
    MANDATORY_ATTRIBUTE = ("index_h", "index_k", "index_l")
    OPTIONAL_ATTRIBUTE = ("fr", "fr_sigma", "fr_calc", "intensity_up_calc", "intensity_down_calc")
    INTERNAL_ATTRIBUTE = ()
    PREFIX = "diffrn_refln"

    def __init__(self, index_h=None, index_k=None, index_l=None,
                 fr=None, fr_sigma=None, intensity_up_calc=None, intensity_down_calc=None):
        super(DiffrnRefln, self).__init__(mandatory_attribute=self.MANDATORY_ATTRIBUTE,
                                          optional_attribute=self.OPTIONAL_ATTRIBUTE,
                                          internal_attribute=self.INTERNAL_ATTRIBUTE,
                                          prefix=self.PREFIX)
        self.index_h = index_h
        self.index_k = index_k
        self.index_l = index_l
        self.fr = fr
        self.intensity_up_calc = intensity_up_calc
        self.intensity_down_calc = intensity_down_calc

        if self.is_defined:
            self.form_object

    @property
    def index_h(self):
        """
Miller indices of the reflection. The values of the Miller
indices in the DIFFN_REFLN category must correspond to the cell
defined by the cell lengths and cell angles in the CELL category.        
        """
        return getattr(self, "__index_h")

    @index_h.setter
    def index_h(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = int(x)
        setattr(self, "__index_h", x_in)

    @property
    def index_k(self):
        """
Miller indices of the reflection. The values of the Miller
indices in the DIFFN_REFLN category must correspond to the cell
defined by the cell lengths and cell angles in the CELL category.        
        """
        return getattr(self, "__index_k")

    @index_k.setter
    def index_k(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = int(x)
        setattr(self, "__index_k", x_in)

    @property
    def index_l(self):
        """
Miller indices of the reflection. The values of the Miller
indices in the DIFFN_REFLN category must correspond to the cell
defined by the cell lengths and cell angles in the CELL category.        
        """
        return getattr(self, "__index_l")

    @index_l.setter
    def index_l(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = int(x)
        setattr(self, "__index_l", x_in)

    @property
    def fr(self):
        """
Measured flip ratio
        """
        return getattr(self, "__fr")

    @fr.setter
    def fr(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__fr", x_in)

    @property
    def fr_sigma(self):
        """
Error bar of measured flip ratio
        """
        return getattr(self, "__fr_sigma")

    @fr_sigma.setter
    def fr_sigma(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__fr_sigma", x_in)

    @property
    def fr_calc(self):
        """
Calculated flip ratio
        """
        return getattr(self, "__fr_calc")

    @fr_calc.setter
    def fr_calc(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__fr_calc", x_in)

    @property
    def intensity_up_calc(self):
        """
Calculated intensity up
        """
        return getattr(self, "__intensity_up_calc")

    @intensity_up_calc.setter
    def intensity_up_calc(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__intensity_up_calc", x_in)

    @property
    def intensity_down_calc(self):
        """
Calculated intensity down
        """
        return getattr(self, "__intensity_down_calc")

    @intensity_down_calc.setter
    def intensity_down_calc(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__intensity_down_calc", x_in)


class DiffrnReflnL(LoopConstr):
    """
Data items in the DIFFRN_REFLN category record details about
the intensities measured in the diffraction experiment.

The DIFFRN_REFLN data items refer to individual intensity
measurements and must be included in looped lists.

Description in cif file::

 loop_
 _diffrn_refln_index_h
 _diffrn_refln_index_k
 _diffrn_refln_index_l
 _diffrn_refln_fr
 _diffrn_refln_fr_sigma
     0    0    8   0.64545   0.01329 
     2    0    6   1.75682   0.0454  

`<https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Cdiffrn_refln.html>`_
    """
    CATEGORY_KEY = ("index_h", "index_k", "index_l")
    ITEM_CLASS = DiffrnRefln

    def __init__(self, item=[], loop_name=""):
        super(DiffrnReflnL, self).__init__(category_key=self.CATEGORY_KEY, item_class=self.ITEM_CLASS,
                                           loop_name=loop_name)
        self.item = item

    def report_agreement_factor_exp(self):
        """
Make a report about experimental agreement factor in string format
        """

        l_chi_sq_exp, l_ag_f_exp = [], []
        l_diff = []
        n_3s, n_2s, n_1s, n_0s = 0, 0, 0, 0
        l_hkl = [(int(_1), int(_2), int(_3)) for _1, _2, _3 in zip(self.index_h, self.index_k, self.index_l)]
        for _hkl, _fr_1, _fr_sigma_1 in zip(l_hkl, self.fr, self.fr_sigma):
            _diff = abs(_fr_1-1.)/float(_fr_sigma_1)
            if _diff >= 3.:
                n_3s +=1
            elif _diff >= 2.:
                n_2s += 1
            elif _diff >= 1.:
                n_1s += 1
            else: 
                n_0s += 1
            if _fr_1 <= 1.:
                l_diff.append(1./_fr_1-1.)
            else:
                l_diff.append(_fr_1-1.)

            _mhkl = (-1 * _hkl[0], -1 * _hkl[1], -1 * _hkl[2])
            if _mhkl in l_hkl:
                ind_mhkl = l_hkl.index(_mhkl)
                _fr_2, _fr_sigma_2 = self.fr[ind_mhkl], self.fr_sigma[ind_mhkl]
                _fr_sigma = 1. / (_fr_sigma_1 ** (-2) + _fr_sigma_2 ** (-2)) ** 0.5
                _fr_average = (_fr_1 * _fr_sigma_1 ** (-2) + _fr_2 * _fr_sigma_2 ** (-2)) * _fr_sigma ** 2
                delta_fr = abs(_fr_1 - _fr_average)
                chi_sq_exp = (delta_fr / _fr_sigma_1) ** 2
                l_chi_sq_exp.append(chi_sq_exp)
                if math.isclose(_fr_1-1., 0):
                    ag_f_exp = abs((_fr_1 - _fr_average) / (_fr_1 - 1.+_fr_sigma_1))
                else:
                    ag_f_exp = abs((_fr_1 - _fr_average) / (_fr_1 - 1.))
                l_ag_f_exp.append(ag_f_exp)
                # print("hkl: {:4} {:4} {:4}".format(_hkl[0], _hkl[1], _hkl[2]))
                # print("chi_sq_exp: {:.3f} ".format(chi_sq_exp))
                # print("ag_f_exp: {:.3f} ".format(ag_f_exp))
        ls_out = []
        n = n_0s + n_1s + n_2s + n_3s
        ls_out.append(f"\n Number of measured reflections: {n:4}")
        ls_out.append(f"      range of h is {min([_[0] for _ in l_hkl]):3},{max([_[0] for _ in l_hkl]):3} ")
        ls_out.append(f"               k is {min([_[1] for _ in l_hkl]):3},{max([_[1] for _ in l_hkl]):3} ")
        ls_out.append(f"               l is {min([_[2] for _ in l_hkl]):3},{max([_[2] for _ in l_hkl]):3} ")
        ls_out.append(f" max(FR_exp - 1) is {max(l_diff):5.3f} ")
        ls_out.append(f" N+1 > |FR_exp - 1|/FR_sigma > N: ")
        ls_out.append(f" |FR_exp - 1|/FR_sigma < 1: {n_0s:4}, {100*float(n_0s)/float(n):5.1f}%  ")
        ls_out.append(f"                     N = 1: {n_1s:4}, {100*float(n_1s)/float(n):5.1f}%  ")
        ls_out.append(f"                     N = 2: {n_2s:4}, {100 * float(n_2s) / float(n):5.1f}% ")
        ls_out.append(f" |FR_exp - 1|/FR_sigma > 3: {n_3s:4}, {100 * float(n_3s) / float(n):5.1f}% ")

        n_friedel = len(l_chi_sq_exp)
        ls_out.append(f"Total number of Friedel reflections is {n_friedel:}.")
        if n_friedel != 0:
            ls_out.append(f"  (|FR_exp-FR_av.|/|FR_sigma|)^2  is {sum(l_chi_sq_exp) / n_friedel:.2f}")
            ls_out.append(f"   |FR_exp-FR_av.|/|FR_exp-1| per reflection is {(100*sum(l_ag_f_exp)/n_friedel):.2f}% ")

        return "\n".join(ls_out)

    def report_chi_sq_exp(self):
        """
Make a report about experimental chi_sq in string format
        """
        ls_out = []
        l_hkl = [(_1, _2, _3) for _1, _2, _3 in zip(self.index_h, self.index_k, self.index_l)]
        l_fr, l_fr_sigma, l_fr_calc = self.fr, self.fr_sigma, self.fr_calc
        if (l_fr is None) & (l_fr_sigma is None) & (l_fr_calc is None):
            return "\n".join(ls_out)
        n = len(l_fr)
        n_1s, n_2s, n_3s = 0, 0, 0
        l_chi_sq, l_worsest, l_af_f, l_af_r = [], [], [], []
        for _hkl, _fr, _fr_sigma, _fr_calc in zip(l_hkl, l_fr, l_fr_sigma, l_fr_calc):
            _diff = abs(float(_fr-_fr_calc)/float(_fr_sigma))
            l_chi_sq.append(_diff**2)
            l_af_f.append(abs(float(_fr - _fr_calc) / float(_fr)))
            if _fr != 1.:
                l_af_r.append(abs(float(_fr-_fr_calc)/float(_fr-1)))
            if _diff <= 1.:
                n_1s +=1
            elif _diff <= 2.:
                n_2s += 1
            elif _diff <= 3.:
                n_3s += 1
            else:
                l_worsest.append((_hkl, _fr, _fr_sigma, _fr_calc, _diff))
        ls_out.append(f"Total number of reflections is {n:}.")
        ls_out.append(f"  (|FR_exp-FR_mod|/|FR_sigma|)^2 per reflection is {sum(l_chi_sq)/float(n):.2f}")
        ls_out.append(f"   |FR_exp-FR_mod|/|FR_exp|      per reflection is {100*sum(l_af_f) / float(n):.2f}%")
        ls_out.append(f"   |FR_exp-FR_mod|/|FR_exp-1|    per reflection is {100*sum(l_af_r) / float(n):.2f}%")
        ls_out.append(f"           (reflections with FR_exp = 1 are excluded)")
        n_worsest = len(l_worsest)
        ls_out.append(f"\nReflections in range  ")
        ls_out.append(f" (N-1)*FR_sigma < |FR_exp - FR_mod| < N*FR_sigma: ")
        ls_out.append(f"      N = 1: {n_1s:}/{n:} ={100*float(n_1s)/float(n):5.1f}% ({2*34.1:4.1f}%, three sigma rule) ")
        ls_out.append(f"      N = 2: {n_2s:}/{n:} ={100 * float(n_2s) / float(n):5.1f}% ({2*(13.6):4.1f}%, three sigma rule)")
        ls_out.append(f"      N = 3: {n_3s:}/{n:} ={100 * float(n_3s) / float(n):5.1f}% ({2*(2.1):4.1f}%, three sigma rule)")
        ls_out.append(f"      N > 3: {n_worsest:}/{n:} ={100 * float(n_worsest) / float(n):5.1f}% ({2*(0.1):4.1f}%, three sigma rule)")
        l_worsest.sort(key=lambda x: x[4], reverse=True)
        if len(l_worsest) > 1:
            if n_worsest > 10: n_worsest = 10
            ls_out.append("\nThe ten worsest reflections:")
            ls_out.append("  h  k  l       FR FR_sigma  FR_calc  diff")
            for (_hkl, _fr, _fr_sigma, _fr_calc, _diff) in l_worsest[:n_worsest]:
                ls_out.append(f"{_hkl[0]:3}{_hkl[1]:3}{_hkl[2]:3}{_fr:9.5f}{_fr_sigma:9.5f}{_fr_calc:9.5f}{_diff:6.1f}")
        return "\n".join(ls_out)
