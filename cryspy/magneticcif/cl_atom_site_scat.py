__author__ = 'ikibalin'
__version__ = "2019_12_03"

import os
import numpy

from pycifstar import Global


from typing import List, Tuple
from cryspy.common.cl_item_constr import ItemConstr
from cryspy.common.cl_loop_constr import LoopConstr
from cryspy.common.cl_fitable import Fitable

from .cl_atom_type_scat import AtomTypeScat, AtomTypeScatL

class AtomSiteScat(ItemConstr):
    """

Description in cif::

 _atom_site_scat_label                         O2
 _atom_site_scat_lande                        2.0
 _atom_site_scat_kappa                        2.0
    """
    MANDATORY_ATTRIBUTE = ("label", )
    OPTIONAL_ATTRIBUTE = ("lande", "kappa")
    INTERNAL_ATTRIBUTE = ("atom_type_scat")
    PREFIX = "atom_site_scat"
    def __init__(self, label=None, lande=2.0, kappa=1.0):
        super(AtomSiteScat, self).__init__(mandatory_attribute=self.MANDATORY_ATTRIBUTE, 
                                           optional_attribute=self.OPTIONAL_ATTRIBUTE, 
                                           internal_attribute=self.INTERNAL_ATTRIBUTE,
                                           prefix=self.PREFIX)

        self.label = label 
        self.lande = lande
        self.kappa = kappa

        if self.is_defined:
            self.form_object

    @property
    def label(self) -> str:
        return getattr(self, "__label")
    @label.setter
    def label(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = str(x)
        setattr(self, "__label", x_in)


    @property
    def lande(self):
        """
Lande factor L:
form factor = <j0> + (2/L - 1) * <j2>

        """
        return getattr(self, "__lande")
    @lande.setter
    def lande(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__lande", x_in)


    @property
    def kappa(self):
        """
kappa factor
        """
        return getattr(self, "__kappa")
    @kappa.setter
    def kappa(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__kappa", x_in)


    @property
    def atom_type_scat(self):
        """
atom_type_scat
        """
        return getattr(self, "__atom_type_scat")

    def __repr__(self) -> str:
        ls_out = []
        ls_out.append("AtomSiteScat: ")
        ls_out.append(f"{str(self):}")
        return "\n".join(ls_out)

    def _show_message(self, s_out: str):
        warnings.warn("***  Error ***\n"+s_out, UserWarning, stacklevel=2)


    @property
    def is_variable(self):
        res = any([self.lande.refinement,
                   self.kappa.refinement])
        return res

    def get_variables(self):
        l_variable = []
        if self.lande is not None:
            if self.lande.refinement: l_variable.append(self.lande)
        if self.kappa is not None:
            if self.kappa.refinement: l_variable.append(self.kappa)
        return l_variable

    def calc_form_factor(self, sthovl):
        atom_type_scat = self.atom_type_scat
        form_factor = atom_type_scat.calc_form_factor(sthovl, lande=self.lande, kappa=self.kappa)
        return form_factor

    def load_atom_type_scat_by_symbol(self, symbol:str):
        flag = True
        _a_t_s = AtomTypeScat.form_by_symbol(symbol)
        setattr(self, "__atom_type_scat", _a_t_s)
        return flag

class AtomSiteScatL(LoopConstr):
    """
Description in cif::

 loop_
 _atom_site_scat_label
 _atom_site_scat_lande
 _atom_site_scat_kappa
  Fe3A    2.00 1.00
  Fe3B    2.00 1.00
    """
    CATEGORY_KEY = ("label", )
    ITEM_CLASS = AtomSiteScat
    def __init__(self, item = [], loop_name=""):
        super(AtomSiteScatL, self).__init__(category_key=self.CATEGORY_KEY, item_class=self.ITEM_CLASS, loop_name=loop_name)
        self.item = item

    def __repr__(self) -> str:
        ls_out = []
        ls_out.append("AtomSiteScatL: ")
        ls_out.append(f"{str(self):}")
        return "\n".join(ls_out)

    def calc_form_factor(self, sthovl):
        form_factor = [_item.calc_form_factor(sthovl) for _item in self.item]
        return numpy.array(list(zip(*form_factor)), dtype=float)
    
    def load_atom_type_scat_by_atom_site(self, atom_site):
        flag = all([_item.load_atom_type_scat_by_symbol(atom_site[_item.label].type_symbol) for _item in self.item])
        return flag
