__author__ = 'ikibalin'
__version__ = "2019_12_06"
import os
import numpy
from pycifstar import Global

import warnings
from typing import List, Tuple
from cryspy.common.cl_item_constr import ItemConstr
from cryspy.common.cl_loop_constr import LoopConstr


class ReflnSusceptibility(ItemConstr):
    """
ReflnSusceptibility
===================
Data items in the REFLN_SUSCEPTIBILITY category record details about the
susceptibility structure factor tensor and magnetization structure factor 
tensor calculated for different reflections.


Description in cif file:
-------------------------
::

  _refln_susceptibility.index_h         2
  _refln_susceptibility.index_k         0
  _refln_susceptibility.index_l         0
  _refln_susceptibility.d_spacing       4.5540
  _refln_susceptibility.sintlambda      0.1235
  _refln_susceptibility.chi_11_calc      0.+0j
  _refln_susceptibility.chi_12_calc      0.+0j
  _refln_susceptibility.chi_13_calc      0.+0j
  _refln_susceptibility.chi_21_calc      0.+0j
  _refln_susceptibility.chi_22_calc      0.+0j
  _refln_susceptibility.chi_23_calc      0.+0j
  _refln_susceptibility.chi_31_calc      0.+0j
  _refln_susceptibility.chi_32_calc      0.+0j
  _refln_susceptibility.chi_33_calc      0.+0j
  _refln_susceptibility.moment_11_calc   0.+0j
  _refln_susceptibility.moment_12_calc   0.+0j
  _refln_susceptibility.moment_13_calc   0.+0j
  _refln_susceptibility.moment_21_calc   0.+0j
  _refln_susceptibility.moment_22_calc   0.+0j
  _refln_susceptibility.moment_23_calc   0.+0j
  _refln_susceptibility.moment_31_calc   0.+0j
  _refln_susceptibility.moment_32_calc   0.+0j
  _refln_susceptibility.moment_33_calc   0.+0j
 

FIXME:
the attribute 'sint/lambda'is replaced by 'sintlambda'
    """
    MANDATORY_ATTRIBUTE = ("index_h", "index_k", "index_l")
    OPTIONAL_ATTRIBUTE = ("sintlambda", "d_spacing", 
                          "chi_11_calc", "chi_12_calc", "chi_13_calc", 
                          "chi_21_calc", "chi_22_calc", "chi_23_calc", 
                          "chi_31_calc", "chi_32_calc", "chi_33_calc", 
                          "moment_11_calc", "moment_12_calc", "moment_13_calc", 
                          "moment_21_calc", "moment_22_calc", "moment_23_calc", 
                          "moment_31_calc", "moment_32_calc", "moment_33_calc")
    INTERNAL_ATTRIBUTE = ()
    PREFIX = "refln_susceptibility"
    def __init__(self, index_h=None, index_k=None, index_l=None, sintlambda=None, d_spacing=None, 
                 chi_11_calc=None, chi_12_calc=None, chi_13_calc=None, 
                 chi_21_calc=None, chi_22_calc=None, chi_23_calc=None, 
                 chi_31_calc=None, chi_32_calc=None, chi_33_calc=None, 
                 moment_11_calc=None, moment_12_calc=None, moment_13_calc=None, 
                 moment_21_calc=None, moment_22_calc=None, moment_23_calc=None, 
                 moment_31_calc=None, moment_32_calc=None, moment_33_calc=None):
        super(ReflnSusceptibility, self).__init__(mandatory_attribute=self.MANDATORY_ATTRIBUTE,
                                                  optional_attribute=self.OPTIONAL_ATTRIBUTE, 
                                                  internal_attribute=self.INTERNAL_ATTRIBUTE,
                                                  prefix=self.PREFIX)
        self.index_h = index_h
        self.index_k = index_k
        self.index_l = index_l
        self.sintlambda = sintlambda
        self.d_spacing = d_spacing
        self.chi_11_calc, self.chi_12_calc, self.chi_13_calc = chi_11_calc, chi_12_calc, chi_13_calc
        self.chi_21_calc, self.chi_22_calc, self.chi_23_calc = chi_21_calc, chi_22_calc, chi_23_calc
        self.chi_31_calc, self.chi_32_calc, self.chi_33_calc = chi_31_calc, chi_32_calc, chi_33_calc
        self.moment_11_calc, self.moment_12_calc, self.moment_13_calc = moment_11_calc, moment_12_calc, moment_13_calc
        self.moment_21_calc, self.moment_22_calc, self.moment_23_calc = moment_21_calc, moment_22_calc, moment_23_calc
        self.moment_31_calc, self.moment_32_calc, self.moment_33_calc = moment_31_calc, moment_32_calc, moment_33_calc
        if self.is_defined:
            self.form_object

    @property
    def index_h(self):
        """
Miller indices of the reflection. The values of the Miller
indices in the REFLN category must correspond to the cell
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
indices in the REFLN category must correspond to the cell
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
indices in the REFLN category must correspond to the cell
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
    def d_spacing(self):
        """
The d spacing in angstroms for this reflection. This is related
to the (sin theta)/lambda value by the expression
_refln_susceptibility_d_spacing = 2/(_refln_sintlambda)

Appears in list containing _refln_index_ 
The permitted range is 0.0 -> infinity

Type: numb        
        """
        return getattr(self, "__d_spacing")
    @d_spacing.setter
    def d_spacing(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__d_spacing", x_in)
    @property
    def sintlambda(self):
        """
The (sin theta)/lambda value in reciprocal angstroms for this
   reflection.

Appears in list containing _refln_index_ 
The permitted range is 0.0 -> infinity

Type: numb        
        """
        return getattr(self, "__sintlambda")
    @sintlambda.setter
    def sintlambda(self, x):
        if x is None:
            x_in = None
        else:
            x_in = float(x)
        setattr(self, "__sintlambda", x_in)


    @property
    def chi_11_calc(self):
        return getattr(self, "__chi_11_calc")
    @chi_11_calc.setter
    def chi_11_calc(self, x):
        if x is None:
            x_in = None
        else:
            x_in = complex(x)
        setattr(self, "__chi_11_calc", x_in)

    @property
    def chi_12_calc(self):
        return getattr(self, "__chi_12_calc")
    @chi_12_calc.setter
    def chi_12_calc(self, x):
        if x is None:
            x_in = None
        else:
            x_in = complex(x)
        setattr(self, "__chi_12_calc", x_in)

    @property
    def chi_13_calc(self):
        return getattr(self, "__chi_13_calc")
    @chi_13_calc.setter
    def chi_13_calc(self, x):
        if x is None:
            x_in = None
        else:
            x_in = complex(x)
        setattr(self, "__chi_13_calc", x_in)

    @property
    def chi_21_calc(self):
        return getattr(self, "__chi_21_calc")
    @chi_21_calc.setter
    def chi_21_calc(self, x):
        if x is None:
            x_in = None
        else:
            x_in = complex(x)
        setattr(self, "__chi_21_calc", x_in)

    @property
    def chi_22_calc(self):
        return getattr(self, "__chi_22_calc")
    @chi_22_calc.setter
    def chi_22_calc(self, x):
        if x is None:
            x_in = None
        else:
            x_in = complex(x)
        setattr(self, "__chi_22_calc", x_in)

    @property
    def chi_23_calc(self):
        return getattr(self, "__chi_23_calc")
    @chi_23_calc.setter
    def chi_23_calc(self, x):
        if x is None:
            x_in = None
        else:
            x_in = complex(x)
        setattr(self, "__chi_23_calc", x_in)

    @property
    def chi_31_calc(self):
        return getattr(self, "__chi_31_calc")
    @chi_31_calc.setter
    def chi_31_calc(self, x):
        if x is None:
            x_in = None
        else:
            x_in = complex(x)
        setattr(self, "__chi_31_calc", x_in)

    @property
    def chi_32_calc(self):
        return getattr(self, "__chi_32_calc")
    @chi_32_calc.setter
    def chi_32_calc(self, x):
        if x is None:
            x_in = None
        else:
            x_in = complex(x)
        setattr(self, "__chi_32_calc", x_in)

    @property
    def chi_33_calc(self):
        return getattr(self, "__chi_33_calc")
    @chi_33_calc.setter
    def chi_33_calc(self, x):
        if x is None:
            x_in = None
        else:
            x_in = complex(x)
        setattr(self, "__chi_33_calc", x_in)

    @property
    def moment_11_calc(self):
        return getattr(self, "__moment_11_calc")
    @moment_11_calc.setter
    def moment_11_calc(self, x):
        if x is None:
            x_in = None
        else:
            x_in = complex(x)
        setattr(self, "__moment_11_calc", x_in)

    @property
    def moment_12_calc(self):
        return getattr(self, "__moment_12_calc")
    @moment_12_calc.setter
    def moment_12_calc(self, x):
        if x is None:
            x_in = None
        else:
            x_in = complex(x)
        setattr(self, "__moment_12_calc", x_in)

    @property
    def moment_13_calc(self):
        return getattr(self, "__moment_13_calc")
    @moment_13_calc.setter
    def moment_13_calc(self, x):
        if x is None:
            x_in = None
        else:
            x_in = complex(x)
        setattr(self, "__moment_13_calc", x_in)

    @property
    def moment_21_calc(self):
        return getattr(self, "__moment_21_calc")
    @moment_21_calc.setter
    def moment_21_calc(self, x):
        if x is None:
            x_in = None
        else:
            x_in = complex(x)
        setattr(self, "__moment_21_calc", x_in)

    @property
    def moment_22_calc(self):
        return getattr(self, "__moment_22_calc")
    @moment_22_calc.setter
    def moment_22_calc(self, x):
        if x is None:
            x_in = None
        else:
            x_in = complex(x)
        setattr(self, "__moment_22_calc", x_in)

    @property
    def moment_23_calc(self):
        return getattr(self, "__moment_23_calc")
    @moment_23_calc.setter
    def moment_23_calc(self, x):
        if x is None:
            x_in = None
        else:
            x_in = complex(x)
        setattr(self, "__moment_23_calc", x_in)

    @property
    def moment_31_calc(self):
        return getattr(self, "__moment_31_calc")
    @moment_31_calc.setter
    def moment_31_calc(self, x):
        if x is None:
            x_in = None
        else:
            x_in = complex(x)
        setattr(self, "__moment_31_calc", x_in)

    @property
    def moment_32_calc(self):
        return getattr(self, "__moment_32_calc")
    @moment_32_calc.setter
    def moment_32_calc(self, x):
        if x is None:
            x_in = None
        else:
            x_in = complex(x)
        setattr(self, "__moment_32_calc", x_in)

    @property
    def moment_33_calc(self):
        return getattr(self, "__moment_33_calc")
    @moment_33_calc.setter
    def moment_33_calc(self, x):
        if x is None:
            x_in = None
        else:
            x_in = complex(x)
        setattr(self, "__moment_33_calc", x_in)


    def __repr__(self):
        ls_out = []
        ls_out.append("ReflnSusceptibility:")
        ls_out.append(f"{str(self):}")
        return "\n".join(ls_out)


    def _show_message(self, s_out: str):
        warnings.warn("***  Error ***\n"+s_out, UserWarning, stacklevel=2)

    @property
    def is_variable(self):
        return False
    
    def get_variables(self):
        return []


class ReflnSusceptibilityL(LoopConstr):
    """
ReflnSusceptibilityL
=====================
Data items in the REFLN_SUSCEPTIBILITY category record details about the
susceptibility structure factor tensor and magnetization structure factor 
tensor calculated for different reflections.


Description in cif file:
-------------------------
::

 loop_
 _refln_susceptibility.index_h    
 _refln_susceptibility.index_k    
 _refln_susceptibility.index_l    
 _refln_susceptibility.d_spacing       
 _refln_susceptibility.sintlambda      
 _refln_susceptibility.chi_11_calc     
 _refln_susceptibility.chi_12_calc     
 _refln_susceptibility.chi_13_calc     
 _refln_susceptibility.chi_21_calc     
 _refln_susceptibility.chi_22_calc     
 _refln_susceptibility.chi_23_calc     
 _refln_susceptibility.chi_31_calc     
 _refln_susceptibility.chi_32_calc     
 _refln_susceptibility.chi_33_calc     
 _refln_susceptibility.moment_11_calc  
 _refln_susceptibility.moment_12_calc  
 _refln_susceptibility.moment_13_calc  
 _refln_susceptibility.moment_21_calc  
 _refln_susceptibility.moment_22_calc  
 _refln_susceptibility.moment_23_calc  
 _refln_susceptibility.moment_31_calc  
 _refln_susceptibility.moment_32_calc  
 _refln_susceptibility.moment_33_calc  
 2 0 0 4.52 0.123 0+0j 0+0j 0+0j 0+0j 0+0j 0+0j 0+0j 0+0j 0+0j  0+0j 0+0j 0+0j 0+0j 0+0j 0+0j 0+0j 0+0j 0+0j
 0 2 0 4.52 0.123 0+0j 0+0j 0+0j 0+0j 0+0j 0+0j 0+0j 0+0j 0+0j  0+0j 0+0j 0+0j 0+0j 0+0j 0+0j 0+0j 0+0j 0+0j

Reference:
`iucr.org <https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Crefln.html>`_
    """
    CATEGORY_KEY = ("index_h", "index_k", "index_l")
    ITEM_CLASS = ReflnSusceptibility
    def __init__(self, item=[], loop_name=""):
        super(ReflnSusceptibilityL, self).__init__(category_key=self.CATEGORY_KEY, item_class=self.ITEM_CLASS, loop_name=loop_name)
        self.item = item

    def __repr__(self) -> str:
        ls_out = []
        ls_out.append("ReflnSusceptibilityL: ")
        ls_out.append(f"{str(self):}")
        return "\n".join(ls_out)
