__author__ = 'ikibalin'
__version__ = "2019_12_11"
import os
import numpy

import warnings
from typing import List, Tuple

from cryspy.common.cl_item_constr import ItemConstr
from cryspy.common.cl_loop_constr import LoopConstr
from cryspy.common.cl_fitable import Fitable

class Texture(ItemConstr):
    """
Describe the setup.

Description in cif file::

 _texture_g_1 0.1239
 _texture_g_2 0.94211
 _texture_h_ax -0.66119
 _texture_k_ax -0.0541
 _texture_l_ax 3.0613
    """
    MANDATORY_ATTRIBUTE = ("g_1", "g_2", "h_ax", "k_ax", "l_ax")
    OPTIONAL_ATTRIBUTE = ()
    INTERNAL_ATTRIBUTE = ()
    PREFIX = "texture"
    def __init__(self, g_1=None, g_2=None, h_ax=None, k_ax=None, l_ax=None):
        super(Texture, self).__init__(mandatory_attribute=self.MANDATORY_ATTRIBUTE, 
                                    optional_attribute=self.OPTIONAL_ATTRIBUTE, 
                                    internal_attribute=self.INTERNAL_ATTRIBUTE,
                                    prefix=self.PREFIX)

        self.g_1 = g_1
        self.g_2 = g_2
        self.h_ax = h_ax
        self.k_ax = k_ax
        self.l_ax = l_ax

        if self.is_defined:
            self.form_object

    @property
    def g_1(self):
        return getattr(self, "__g_1")
    @g_1.setter
    def g_1(self, x):
        if x is None:
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__g_1", x_in)


    @property
    def g_2(self):
        return getattr(self, "__g_2")
    @g_2.setter
    def g_2(self, x):
        if x is None:
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__g_2", x_in)

    @property
    def h_ax(self):
        return getattr(self, "__h_ax")
    @h_ax.setter
    def h_ax(self, x):
        if x is None:
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__h_ax", x_in)

    @property
    def k_ax(self):
        return getattr(self, "__k_ax")
    @k_ax.setter
    def k_ax(self, x):
        if x is None:
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__k_ax", x_in)

    @property
    def l_ax(self):
        return getattr(self, "__l_ax")
    @l_ax.setter
    def l_ax(self, x):
        if x is None:
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__l_ax", x_in)

    def __repr__(self):
        ls_out = []
        ls_out.append(f"Texture:\n{str(self):}")
        return "\n".join(ls_out)

        
    @property
    def is_variable(self) -> bool:
        """
Output: True if there is any refined parameter
        """
        _l = [self.g_1.refinement,
              self.g_2.refinement, 
              self.h_ax.refinement,
              self.k_ax.refinement, 
              self.l_ax.refinement]
        res = any(_l)
        return res

    def get_variables(self) -> List:
        """
Output: the list of the refined parameters
        """
        l_variable = []
        if self.g_1.refinement: l_variable.append(self.g_1)
        if self.g_2.refinement: l_variable.append(self.g_2)
        if self.h_ax.refinement: l_variable.append(self.h_ax)
        if self.k_ax.refinement: l_variable.append(self.k_ax)
        if self.l_ax.refinement: l_variable.append(self.l_ax)
        return l_variable
