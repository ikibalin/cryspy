__author__ = 'ikibalin'
__version__ = "2020_01_02"
import os
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
    OPTIONAL_ATTRIBUTE = ("fr", "fr_sigma", "intensity_up_calc", "intensity_down_calc") 
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
        if x is None:
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
        if x is None:
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
        if x is None:
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
        if x is None:
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
        if x is None:
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
        if x is None:
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
        if x is None:
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
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__intensity_down_calc", x_in)

    def __repr__(self):
        ls_out = []
        ls_out.append("DiffrnRefln:")
        ls_out.append(f"{str(self):}")
        return "\n".join(ls_out)



    def _show_message(self, s_out: str):
        warnings.warn("***  Error ***\n"+s_out, UserWarning, stacklevel=2)

    @property
    def is_variable(self):
        return False
    
    def get_variables(self):
        return []


    def print_agreement_factor_exp(self):
        l_chi_sq_exp, l_ag_f_exp = [], []
        l_hkl = [(int(_1), int(_2), int(_3))for _1, _2, _3 in zip(self.h, self.k, self.l)]
        for _hkl, _fr_1, _fr_sigma_1 in zip(l_hkl, self.fr, self.fr_sigma):
            _mhkl = (-1*_hkl[0], -1*_hkl[1], -1*_hkl[2])
            if _mhkl in l_hkl:
                ind_mhkl = l_hkl.index(_mhkl)
                _fr_2, _fr_sigma_2 = self.fr[ind_mhkl], self.fr_sigma[ind_mhkl]
                _fr_sigma = 1./(_fr_sigma_1**(-2)+_fr_sigma_2**(-2))**0.5
                _fr_average = (_fr_1*_fr_sigma_1**(-2)+_fr_2*_fr_sigma_2**(-2))*_fr_sigma**2
                delta_fr = abs(_fr_1-_fr_average)
                chi_sq_exp = (delta_fr/_fr_sigma_1)**2
                l_chi_sq_exp.append(chi_sq_exp)
                ag_f_exp = abs((_fr_1-_fr_average)/(_fr_1-1.))
                l_ag_f_exp.append(ag_f_exp)
                #print("hkl: {:4} {:4} {:4}".format(_hkl[0], _hkl[1], _hkl[2]))
                #print("chi_sq_exp: {:.3f} ".format(chi_sq_exp))
                #print("ag_f_exp: {:.3f} ".format(ag_f_exp))
        ls_out = []
        n_friedel = len(l_chi_sq_exp)
        ls_out.append("number of Friedel reflections is {:}".format(n_friedel))
        if n_friedel != 0:
            ls_out.append("agreement factor_exp/n is {:.3f}".format(sum(l_ag_f_exp)/n_friedel))
            ls_out.append("chi_sq_exp/n is {:.3f}".format(sum(l_chi_sq_exp)/n_friedel))
        return "\n".join(ls_out) 


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
        super(DiffrnReflnL, self).__init__(category_key=self.CATEGORY_KEY, item_class=self.ITEM_CLASS, loop_name=loop_name)
        self.item = item

    def __repr__(self) -> str:
        ls_out = []
        ls_out.append("DiffrnReflnL: ")
        ls_out.append(f"{str(self):}")
        return "\n".join(ls_out)