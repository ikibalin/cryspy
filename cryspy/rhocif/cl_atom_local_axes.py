__author__ = 'ikibalin'
__version__ = "2019_12_05"
import os
import numpy
from pycifstar import Global

import warnings
from typing import List, Tuple
from cryspy.common.cl_item_constr import ItemConstr
from cryspy.common.cl_loop_constr import LoopConstr
from cryspy.common.cl_fitable import Fitable
import cryspy.corecif.CONSTANTS_AND_FUNCTIONS as CONSTANTS_AND_FUNCTIONS 


class AtomLocalAxes(ItemConstr):
    """
This category allows the definition of local axes around each
atom in terms of vectors between neighbouring atoms.
High-resolution X-ray diffraction methods enable the
determination of the electron density distribution in crystal
lattices and molecules, which in turn allows for a
characterization of chemical interactions (Coppens, 1997;
Koritsanszky & Coppens, 2001). This is accomplished by the
construction of a mathematical model of the charge density
in a crystal and then by fitting the parameters of such a
model to the experimental pattern of diffracted X-rays. The
model on which this dictionary is based is the so-called
multipole formalism proposed by Hansen & Coppens (1978). In
this model, the electron density in a crystal is described
by a sum of aspherical "pseudoatoms" where the pseudoatom
density has the form defined in the _atom_rho_multipole_* items.
Each pseudoatom density consists of terms representing the
core density, the spherical part of the valence density and
the deviation of the valence density from sphericity. The
continuous electron density in the crystal is then modelled
as a sum of atom-centred charge distributions. Once the
experimental electron density has been established, the
"atoms in molecules" theory of Bader (1990) provides tools for
the interpretation of the density distribution in terms of its
topological properties.


Description in cif file::


  _atom_local_axes_atom_label     Ni2+(1)
  _atom_local_axes_atom0          DUM0
  _atom_local_axes_ax1            Z
  _atom_local_axes_atom1          Ni2+(1)
  _atom_local_axes_atom2          N(1)
  _atom_local_axes_ax2            X
      
`<https://www.iucr.org/__data/iucr/cifdic_html/1/cif_rho.dic/Catom_local_axes.html>`_
    """    
    MANDATORY_ATTRIBUTE = ("atom_label", "atom0", "ax1", "atom1", "atom2", "ax2")
    OPTIONAL_ATTRIBUTE = ()
    INTERNAL_ATTRIBUTE = ()
    ACCESIBLE_AX = ("X", "Y", "Z", "-X", "-Y", "-Z", "+X", "+Y", "+Z", "x", "y", "z", "-x", "-y", "-z", "+x", "+y", "+z")
    REPLACED_AX = ("X", "Y", "Z", "-X", "-Y", "-Z", "X", "Y", "Z", "X", "Y", "Z", "-X", "-Y", "-Z", "X", "Y", "Z")
    PREFIX = "atom_local_axes"

    def __init__(self, atom_label=None, atom0=None, 
                 ax1=None, atom1=None, atom2=None, ax2=None):
        super(AtomLocalAxes, self).__init__(mandatory_attribute=self.MANDATORY_ATTRIBUTE, 
                                            optional_attribute=self.OPTIONAL_ATTRIBUTE, 
                                            internal_attribute=self.INTERNAL_ATTRIBUTE,
                                            prefix=self.PREFIX)

        self.atom_label = atom_label
        self.atom0 = atom0
        self.ax1 = ax1
        self.atom1 = atom1
        self.atom2 = atom2
        self.ax2 = ax2

        if self.is_defined:
            self.form_object

    @property
    def atom_label(self)->str:
        """
This item is used to identify an atom for which a local axis
system is to be defined.  Its value must be identical to one
of the values given in the _atom_site_label list.

:Type: char

`<https://www.iucr.org/__data/iucr/cifdic_html/1/cif_rho.dic/Iatom_local_axes_atom_label.html>`_
        """
        return getattr(self, "__atom_label")
    @atom_label.setter
    def atom_label(self, x):
        if x is None:
            x_in = None
        else:
            x_in = str(x)
        setattr(self, "__atom_label", x_in)

    @property
    def atom0(self):
        """
Specifies 'atom0' in the definition of a local axis frame.
The definition employs three atom-site labels, 'atom0', 'atom1'
and 'atom2', and two axis labels, 'ax1' and 'ax2', having values
'+/-X', '+/-Y' or '+/-Z'. For the atom defined by
'_atom_local_axes_atom_label', whose nuclear position is taken
as the origin, local axis 'ax1' is the vector from the origin to
atom0, axis 'ax2' is perpendicular to 'ax1' and lies in the
plane of 'ax1' and a vector
passing through the origin parallel to the vector atom1 -> atom2
(its positive direction making an acute angle with the vector
parallel to atom1 -> atom2), and a right-handed orthonormal
vector triplet is formed from the vector product of these two
vectors. In most cases, atom1 will be the same as the atom
specified by _atom_local_axes_atom_label. One or more 'dummy'
atoms (with arbitrary labels) may be used in the vector
definitions, specified with zero occupancy in the _atom_site_
description.  The values of *_atom0, *_atom1 and *_atom2 must
be identical to values given in the _atom_site_label list.

:Type: char

`<https://www.iucr.org/__data/iucr/cifdic_html/1/cif_rho.dic/Iatom_local_axes_atom0.html>`_
        """
        return getattr(self, "__atom0")
    @atom0.setter
    def atom0(self, x):
        if x is None:
            x_in = None
        else:
            x_in = str(x)
        setattr(self, "__atom0", x_in)


    @property
    def atom1(self):
        """
Specifies 'atom1' in the definition of a local axis frame.
The definition employs three atom-site labels, 'atom0', 'atom1'
and 'atom2', and two axis labels, 'ax1' and 'ax2', having values
'+/-X', '+/-Y' or '+/-Z'. For the atom defined by
'_atom_local_axes_atom_label', whose nuclear position is taken
as the origin, local axis 'ax1' is the vector from the origin to
atom0, axis 'ax2' is perpendicular to 'ax1' and lies in the
plane of 'ax1' and a vector
passing through the origin parallel to the vector
atom1 -> atom2 (its positive direction making an acute angle
with the vector parallel to atom1 -> atom2), and a right-handed
orthonormal vector triplet is formed from the vector product
of these two vectors. In most cases, atom1 will be the same
as the atom specified by _atom_local_axes_atom_label. One or
more 'dummy' atoms (with arbitrary labels) may be used in the
vector definitions, specified with zero occupancy in the
_atom_site_ description.  The values of *_atom0, *_atom1 and
*_atom2 must be identical to values given in the
_atom_site_label list.


:Type: char

`<https://www.iucr.org/__data/iucr/cifdic_html/1/cif_rho.dic/Iatom_local_axes_atom1.html>`_
        """
        return getattr(self, "__atom1")
    @atom1.setter
    def atom1(self, x):
        if x is None:
            x_in = None
        else:
            x_in = str(x)
        setattr(self, "__atom1", x_in)


    @property
    def atom2(self):
        """
Specifies 'atom2' in the definition of a local axis frame.
The definition employs three atom-site labels, 'atom0', 'atom1'
and 'atom2', and two axis labels, 'ax1' and 'ax2', having values
'+/-X', '+/-Y' or '+/-Z'. For the atom defined by
'_atom_local_axes_atom_label', whose nuclear position is taken
as the origin, local axis 'ax1' is the vector from the origin to
atom0, axis 'ax2' is perpendicular to 'ax1' and lies in the
plane of 'ax1' and a vector
passing through the origin parallel to the vector atom1 -> atom2
(its positive direction making an acute angle with the vector
parallel to atom1 -> atom2), and a right-handed orthonormal
vector triplet is formed from the vector product of these
two vectors. In most cases, atom1 will be the same as the
atom specified by _atom_local_axes_atom_label. One or more
'dummy' atoms (with arbitrary labels) may be used in the
vector definitions, specified with zero occupancy in the
_atom_site_ description.  The values of *_atom0, *_atom1 and
*_atom2 must be identical to values given in the
_atom_site_label list.

:Type: char

`<https://www.iucr.org/__data/iucr/cifdic_html/1/cif_rho.dic/Iatom_local_axes_atom2.html>`_
        """
        return getattr(self, "__atom2")
    @atom2.setter
    def atom2(self, x):
        if x is None:
            x_in = None
        else:
            x_in = str(x)
        setattr(self, "__atom2", x_in)


    @property
    def ax1(self):
        """
Specifies 'ax1' in the definition of a local axis frame.
The definition employs three atom-site labels, 'atom0', 'atom1'
and 'atom2', and two axis labels, 'ax1' and 'ax2', having values
'+/-X', '+/-Y' or '+/-Z'. For the atom defined by
'_atom_local_axes_atom_label', whose nuclear position is taken
as the origin, local axis 'ax1' is the vector from the origin to
atom0, axis 'ax2' is perpendicular to 'ax1' and lies in the
plane of 'ax1' and a vector
passing through the origin parallel to the vector atom1 -> atom2
(its positive direction making an acute angle with the vector
parallel to atom1 -> atom2), and a right-handed orthonormal
vector triplet is formed from the vector product of these two
vectors. In most cases, atom1 will be the same as the atom
specified by _atom_local_axes_atom_label. One or more 'dummy'
atoms (with arbitrary labels) may be used in the vector
definitions, specified with zero occupancy in the _atom_site_
description.  The values of *_atom0, *_atom1 and *_atom2 must
be identical to values given in the _atom_site_label list.

Appears in list containing _atom_local_axes_atom_label 

The data value must be one of the following:

x, X, y, Y, z, Z, 
+x, +X, +y, +Y, +z, +Z, 
-x, -X, -y, -Y, -z, -Z 

`<https://www.iucr.org/__data/iucr/cifdic_html/1/cif_rho.dic/Iatom_local_axes_ax1.html>`_
        """
        return getattr(self, "__ax1")
    @ax1.setter
    def ax1(self, x):
        if x is None:
            x_in = None
        else:
            if not(x in self.ACCESIBLE_AX):
                warnings.warn(f"ax1 '{x:}' is not supported", UserWarning, stacklevel=2)
                x_in = None            
            else:
                x_in = self.REPLACED_AX[self.ACCESIBLE_AX.index(x)]
        setattr(self, "__ax1", x_in)


    @property
    def ax2(self):
        """
Specifies 'ax2' in the definition of a local axis frame.
The definition employs three atom-site labels, 'atom0', 'atom1'
and 'atom2', and two axis labels, 'ax1' and 'ax2', having values
'+/-X', '+/-Y' or '+/-Z'. For the atom defined by
'_atom_local_axes_atom_label', whose nuclear position is taken
as the origin, local axis 'ax1' is the vector from the origin to
atom0, axis 'ax2' is perpendicular to 'ax1' and lies in the
plane of 'ax1' and a vector
passing through the origin parallel to the vector atom1 -> atom2
(its positive direction making an acute angle with the vector
parallel to atom1 -> atom2), and a right-handed orthonormal
vector triplet is formed from the vector product of these two
vectors. In most cases, atom1 will be the same as the atom
specified by _atom_local_axes_atom_label. One or more 'dummy'
atoms (with arbitrary labels) may be used in the vector
definitions, specified with zero occupancy in the _atom_site_
description.  The values of *_atom0, *_atom1 and *_atom2 must
be identical to values given in the _atom_site_label list.


Appears in list containing _atom_local_axes_atom_label 

The data value must be one of the following:

x, X, y, Y, z, Z, 
+x, +X, +y, +Y, +z, +Z, 
-x, -X, -y, -Y, -z, -Z 

`<https://www.iucr.org/__data/iucr/cifdic_html/1/cif_rho.dic/Iatom_local_axes_ax2.html>`_
        """
        return getattr(self, "__ax2")
    @ax2.setter
    def ax2(self, x):
        if x is None:
            x_in = None
        else:
            if not(x in self.ACCESIBLE_AX):
                warnings.warn(f"ax2 '{x:}' is not supported", UserWarning, stacklevel=2)
                x_in = None            
            else:
                x_in = self.REPLACED_AX[self.ACCESIBLE_AX.index(x)]
        setattr(self, "__ax2", x_in)


    def __repr__(self) -> str:
        ls_out = []
        ls_out.append("AtomLocalAxes: ")
        ls_out.append(f"{str(self):}")
        return "\n".join(ls_out)


class AtomLocalAxesL(LoopConstr):
    """
This category allows the definition of local axes around each
atom in terms of vectors between neighbouring atoms.
High-resolution X-ray diffraction methods enable the
determination of the electron density distribution in crystal
lattices and molecules, which in turn allows for a
characterization of chemical interactions (Coppens, 1997;
Koritsanszky & Coppens, 2001). This is accomplished by the
construction of a mathematical model of the charge density
in a crystal and then by fitting the parameters of such a
model to the experimental pattern of diffracted X-rays. The
model on which this dictionary is based is the so-called
multipole formalism proposed by Hansen & Coppens (1978). In
this model, the electron density in a crystal is described
by a sum of aspherical "pseudoatoms" where the pseudoatom
density has the form defined in the _atom_rho_multipole_* items.
Each pseudoatom density consists of terms representing the
core density, the spherical part of the valence density and
the deviation of the valence density from sphericity. The
continuous electron density in the crystal is then modelled
as a sum of atom-centred charge distributions. Once the
experimental electron density has been established, the
"atoms in molecules" theory of Bader (1990) provides tools for
the interpretation of the density distribution in terms of its
topological properties.


Description in cif file::

  loop_
  _atom_local_axes_atom_label
  _atom_local_axes_atom0
  _atom_local_axes_ax1
  _atom_local_axes_atom1
  _atom_local_axes_atom2
  _atom_local_axes_ax2
      Ni2+(1)  DUM0      Z    Ni2+(1)  N(1)      X

`<https://www.iucr.org/__data/iucr/cifdic_html/1/cif_rho.dic/Catom_local_axes.html>`_
    """
    CATEGORY_KEY = ("atom_label", )
    ITEM_CLASS = AtomLocalAxes
    def __init__(self, item=[], loop_name=""):
        super(AtomLocalAxesL, self).__init__(category_key=self.CATEGORY_KEY, item_class=self.ITEM_CLASS, loop_name=loop_name)
        self.item = item

    def __repr__(self) -> str:
        ls_out = []
        ls_out.append("AtomLocalAxesL: ")
        ls_out.append(f"{str(self):}")
        return "\n".join(ls_out)
