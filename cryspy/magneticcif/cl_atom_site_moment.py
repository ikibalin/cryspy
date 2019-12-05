"""
AtomSiteMoment, AtomSiteMomentL classes are defined. 
"""
__author__ = 'ikibalin'
__version__ = "2019_12_03"

import os
import numpy

from pycifstar import Global


from typing import List, Tuple
from cryspy.common.cl_item_constr import ItemConstr
from cryspy.common.cl_loop_constr import LoopConstr
from cryspy.common.cl_fitable import Fitable

class AtomSiteMoment(ItemConstr):
    """
AtomSiteMomentL
===================
This category provides a loop for presenting the magnetic moments 
of atoms in one of several coordinate systems. This is a child category 
of the AtomSite category, so that the magnetic moments can either be 
listed alongside the non-magnetic atom properties in the main AtomSite loop
(not realized) or be listed in a separate loop (realized)

Category key: 
--------------------
atom_site_moment_label

Description in cif file:
---------------------------
::

 _atom_site_moment_label              Fe3A
 _atom_site_moment_crystalaxis_x       4.8
 _atom_site_moment_crystalaxis_y       0.0 
 _atom_site_moment_crystalaxis_z       0.0
 
 Reference:
---------------
`iucr.org <ftp://ftp.iucr.org/cifdics/cif_mag_0.9.7.dic.pdf>`_
    """    
    MANDATORY_ATTRIBUTE = ("label", )
    OPTIONAL_ATTRIBUTE = ("crystalaxis_x", "crystalaxis_y", "crystalaxis_z")
    INTERNAL_ATTRIBUTE = ()
    PREFIX = "atom_site_moment"
    def __init__(self, label=None, crystalaxis_x=None, crystalaxis_y=None, crystalaxis_z=None):
        super(AtomSiteMoment, self).__init__(mandatory_attribute=self.MANDATORY_ATTRIBUTE, 
                                             optional_attribute=self.OPTIONAL_ATTRIBUTE, 
                                             internal_attribute=self.INTERNAL_ATTRIBUTE,
                                             prefix=self.PREFIX)
        self.label = label
        self.crystalaxis_x = crystalaxis_x
        self.crystalaxis_y = crystalaxis_y
        self.crystalaxis_z = crystalaxis_z

        if self.is_defined:
            self.form_object
        
    def __repr__(self):
        ls_out = ["AtomSiteMoment:"]
        ls_out.append(str(self))
        return "\n".join(ls_out)

    @property
    def label(self):
        """
This label in a unique identifier for a particular site in the asymmetric unit of the crystal unit cell
Type: char
        """
        return getattr(self, "__label")
    @label.setter
    def label(self, x):
        if x is None:
            x_in = None
        else:
            x_in = str(x)
        setattr(self, "__label", x_in)

    @property
    def crystalaxis_x(self):
        """
The atom-site magnetic moment vector specified as a 
projection onto the axes of the unit cell.
in mu_B

Default: 0.

Type: float
        """
        return getattr(self, "__crystalaxis_x")
    @crystalaxis_x.setter
    def crystalaxis_x(self, x):
        if x is None:
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__crystalaxis_x", x_in)


    @property
    def crystalaxis_y(self):
        """
The atom-site magnetic moment vector specified as a 
projection onto the axes of the unit cell.
in mu_B

Default: 0.

Type: float
        """
        return getattr(self, "__crystalaxis_y")
    @crystalaxis_y.setter
    def crystalaxis_y(self, x):
        if x is None:
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__crystalaxis_y", x_in)

    @property
    def crystalaxis_z(self):
        """
The atom-site magnetic moment vector specified as a 
projection onto the axes of the unit cell.
in mu_B

Default: 0.

Type: float
        """
        return getattr(self, "__crystalaxis_z")
    @crystalaxis_z.setter
    def crystalaxis_z(self, x):
        if x is None:
            x_in = None
        elif isinstance(x, Fitable):
            x_in = x
        else:
            x_in = Fitable()
            flag = x_in.take_it(x)
        setattr(self, "__crystalaxis_z", x_in)

    def _show_message(self, s_out: str):
        warnings.warn("***  Error ***\n"+s_out, UserWarning, stacklevel=2)


    @property
    def is_variable(self):
        res = any([self.crystalaxis_x.refinement])
        res = any([self.crystalaxis_y.refinement])
        res = any([self.crystalaxis_z.refinement])
        return res

    def get_variables(self):
        l_variable = []
        if self.crystalaxis_x.refinement: l_variable.append(self.crystalaxis_x)
        if self.crystalaxis_y.refinement: l_variable.append(self.crystalaxis_y)
        if self.crystalaxis_z.refinement: l_variable.append(self.crystalaxis_z)
        return l_variable


    @property
    def moment(self):
        np_x = numpy.array(self.crystalaxis_x, dtype=float)
        np_y = numpy.array(self.crystalaxis_y, dtype=float)
        np_z = numpy.array(self.crystalaxis_z, dtype=float)
        np_moment = numpy.sqrt(numpy.square(np_x)+numpy.square(np_y)+numpy.square(np_z))
        return np_moment



    def calc_zeeman(self, field_cryst):
        """
        !!! It is valid only for cubic crystals
        moment*sin(moment^field) at angle decreesing 
        """
        h_a, h_b, h_c = field_cryst[0], field_cryst[1], field_cryst[2]
        if abs(h_a)+abs(h_b)+abs(h_c) ==0.:
            np_val = numpy.zeros(len(self.label), dtype=float)
            return np_val
        else:
            mod_h = (abs(h_a)**2+abs(h_b)**2+abs(h_c)**2)**0.5
            h_a_n, h_b_n, h_c_n = h_a/mod_h, h_b/mod_h, h_c/mod_h
            np_moment_sq = numpy.square(self.moment)
            np_x = numpy.array(self.crystalaxis_x, dtype=float)
            np_y = numpy.array(self.crystalaxis_y, dtype=float)
            np_z = numpy.array(self.crystalaxis_z, dtype=float)
            np_val = h_a_n*np_x+h_b_n*np_y+h_c_n*np_z
            np_res = numpy.sqrt(np_moment_sq-numpy.square(np_val))
        return np_res

class AtomSiteMomentL(LoopConstr):
    """
AtomSiteMomentL
===================
This category provides a loop for presenting the magnetic moments 
of atoms in one of several coordinate systems. This is a child category 
of the AtomSite category, so that the magnetic moments can either be 
listed alongside the non-magnetic atom properties in the main AtomSite loop
(not realized) or be listed in a separate loop (realized)

Category key: 
--------------------
atom_site_moment_label

Description in cif file:
---------------------------
::

 loop_                                     
 _atom_site_moment_label
 _atom_site_moment_crystalaxis_x
 _atom_site_moment_crystalaxis_y
 _atom_site_moment_crystalaxis_z
 Fe3A 4.8  0.0  0.0
 Fe3B 0.0 -4.5  0.0
 
 Reference:
---------------
`iucr.org <ftp://ftp.iucr.org/cifdics/cif_mag_0.9.7.dic.pdf>`_
    """    
    CATEGORY_KEY = ("label", )
    ITEM_CLASS = AtomSiteMoment
    def __init__(self, item = [], loop_name=""):
        super(AtomSiteMomentL, self).__init__(category_key=self.CATEGORY_KEY, item_class=self.ITEM_CLASS, loop_name=loop_name)
        self.item = item

    def __repr__(self) -> str:
        ls_out = []
        ls_out.append("AtomSiteMomentL: ")
        ls_out.append(f"{str(self):}")
        return "\n".join(ls_out)
