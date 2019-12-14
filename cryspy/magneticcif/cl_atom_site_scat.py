__author__ = 'ikibalin'
__version__ = "2019_12_03"

import os
import numpy

from pycifstar import Global


from typing import List, Tuple
from cryspy.common.cl_item_constr import ItemConstr
from cryspy.common.cl_loop_constr import LoopConstr
from cryspy.common.cl_fitable import Fitable


class AtomSiteScat(ItemConstr):
    """

Description in cif::

 _atom_site_scat_label                         O2
 _atom_site_scat_lande                        2.0
 _atom_site_scat_kappa                        2.0
    """
    MANDATORY_ATTRIBUTE = ("label", )
    OPTIONAL_ATTRIBUTE = ("lande", "kappa")
    INTERNAL_ATTRIBUTE = ()
    PREFIX = "atom_site_scat"
    def __init__(self, label=None, lande=None, kappa=None):
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
        if x is None:
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
        if x is None:
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
        if x is None:
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__kappa", x_in)

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
        if self.lande.refinement: l_variable.append(self.lande)
        if self.kappa.refinement: l_variable.append(self.kappa)
        return l_variable



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
