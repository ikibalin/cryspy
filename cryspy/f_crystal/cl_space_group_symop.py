"""
Defines SpaceGoupSymopEl and SpaceGoupSymop classes. 

"""

__author__ = 'ikibalin'
__version__ = "2019_11_26"

import warnings

from typing import List, Tuple
from cryspy.f_common.cl_item_constr import ItemConstr
from cryspy.f_common.cl_loop_constr import LoopConstr

class SpaceGroupSymopEl(ItemConstr):
    """
SpaceGroupSymopEl
=================

Contains information about the symmetry operations of the
space group.

Description in cif file:
-------------------------
_space_group_symop.id  1    
_space_group_symop.operation_xyz   x,y,z              
_space_group_symop.operation_description 'identity mapping'

Attributes:
-----------
- id
- operation_xyz
- operation_description
- generator_xyz
- sg_id

Methods:
---------
-  

    """
    MANDATORY_ATTRIBUTE = ("id", "operation_xyz")
    OPTIONAL_ATTRIBUTE = ("operation_description", "generator_xyz", "sg_id")
    INTERNAL_ATTRIBUTE = ()
    PREFIX = "space_group_symop"
    def __init__(self, id=None, operation_xyz=None, operation_description=None, generator_xyz=None, sg_id=None):
        super(SpaceGroupSymopEl, self).__init__(mandatory_attribute=self.MANDATORY_ATTRIBUTE, 
                                                optional_attribute=self.OPTIONAL_ATTRIBUTE, 
                                                internal_attribute=self.INTERNAL_ATTRIBUTE,
                                                prefix=self.PREFIX)
        self.id = id
        self.operation_xyz = operation_xyz
        self.operation_description = operation_description
        self.generator_xyz = generator_xyz
        self.sg_id = sg_id
        
    @property
    def id(self) -> str:
        """
An arbitrary identifier that uniquely labels each symmetry opera-
tion in the list.
        """
        return getattr(self, "__id")
    @id.setter
    def id(self, x):
        if x is None:
            x_in = None
        else:
            x_in = str(x)
        setattr(self, "__id", x_in)

    @property
    def operation_xyz(self) -> str:
        """
A parsable string giving one of the symmetry operations of the
space group in algebraic form. If W is a matrix representation of
the rotational part of the symmetry operation defined by the posi-
tions and signs of x, y and z, and w is a column of translations
defined by the fractions, an equivalent positionx is generated from
a given position x by


..math::
x^{'} = Wx + w

When a list of symmetry operations is given, it is assumed that the
list contains all the operations of the space group (including the
identity operation) as given by the representatives of the general
position in International Tables for Crystallography Volume A.

Reference:
------------
International Tables for Crystallography (2002).
Volume A, Space-group symmetry, edited by Th. Hahn, 5th ed.
Dordrecht: Kluwer Academic Publishers.

Example: 
-----------
'x,1/2-y,1/2+z'
        """
        return getattr(self, "__operation_xyz")
    @operation_xyz.setter
    def operation_xyz(self, x):
        if x is None:
            x_in = None
        else:
            x_in = "".join(str(x).split())
        setattr(self, "__operation_xyz", x_in)

    @property
    def operation_description(self) -> str:
        """
An optional text description of a particular symmetry operation of
the space group
        """
        return getattr(self, "__operation_description")
    @operation_description.setter
    def operation_description(self, x):
        if x is None:
            x_in = None
        else:
            x_in = str(x)
        setattr(self, "__operation_description", x_in)

    @property
    def generator_xyz(self):
        """
A parsable string giving one of the symmetry generators of the
space group in algebraic form. If W is a matrix representation of
the rotational part of the generator defined by the positions and
signs of x, y and z, and w is a column of translations defined by
the fractions, an equivalent position x ? is generated from a given
position x by

..math::
x^{'} = Wx + w

When a list of symmetry generators is given, it is assumed that the
completelist of symmetryoperationsofthespacegroup(including
the identity operation) can be generated through repeated multipli-
cation of the generators, that is, (W 3 ,w 3 ) is an operation of the
space group if (W 2 ,w 2 ) and (W 1 ,w 1 ) [where (W 1 ,w 1 ) is applied
first] are either operations or generators and

..math::
W_{3} = W_{2} \times  W_{1},

..math::
w_{3} = W_{2} \times w_{1} + w_{2}

Reference:
-------------
International Tables for Crystallography (2002).
Volume A, Space-group symmetry, edited by Th. Hahn, 5th ed.
Dordrecht: Kluwer Academic Publishers.

Example: 
---------
'x,1/2-y,1/2+z'
        """
        return getattr(self, "__generator_xyz")
    @generator_xyz.setter
    def generator_xyz(self, x) -> str:
        if x is None:
            x_in = None
        else:
            x_in = "".join(str(x).split())
        setattr(self, "__generator_xyz", x_in)

    @property
    def sg_id(self):
        """
A child of _space_group.id allowing the symmetry operation to
be identified with a particular space group
        """
        return getattr(self, "__sg_id")

    @sg_id.setter
    def sg_id(self, x):
        if x is None:
            x_in = None
        else:
            x_in = int(x)
        setattr(self, "__sg_id", x_in)
        




class SpaceGroupSymop(LoopConstr):
    """
SpaceGroupSymop
=================

Contains information about the symmetry operations of the
space group.

Description in cif file:
-------------------------

 loop_
_space_group_symop.id
_space_group_symop.operation_xyz
_space_group_symop.operation_description
  1    x,y,z              'identity mapping'
  2    -x,-y,-z           'inversion'
  3    -x,1/2+y,1/2-z
              '2-fold screw rotation with axis in (0,y,1/4)'
  4    x,1/2-y,1/2+z
            'c glide reflection through the plane (x,1/4,y)'

Attributes:
-----------
- id
- operation_xyz
- operation_description
- generator_xyz


Mandatory attribute: 
---------------------
 - id (category key, 1st)
 - operation_xyz


Optional attribute: 
---------------------
- operation_description
- generator_xyz


Class methods:
---------
- sg_id


reference: https://www.iucr.org/__data/iucr/cifdic_html/2/cif_sym.dic/Cspace_group_symop.html
    """
    CATEGORY_KEY = ("id", )
    ITEM_CLASS = SpaceGroupSymopEl
    def __init__(self, item = [], label=""):
        super(SpaceGroupSymop, self).__init__(category_key=self.CATEGORY_KEY, item_class=self.ITEM_CLASS, label=label)
        self.item = item



