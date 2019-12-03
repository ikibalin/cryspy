"""
Defines SpaceGroupSymop and SpaceGroupSymopL classes. 

"""

__author__ = 'ikibalin'
__version__ = "2019_11_26"

import numpy
from fractions import Fraction
import warnings
import copy

from typing import List, Tuple
from cryspy.f_common.cl_item_constr import ItemConstr
from cryspy.f_common.cl_loop_constr import LoopConstr
import cryspy.f_crystal.symcif.CONSTANTS_AND_FUNCTIONS as CONSTANTS_AND_FUNCTIONS

class SpaceGroupSymop(ItemConstr):
    """
SpaceGroupSymop
=================

Contains information about the symmetry operations of the
space group.

Description in cif file:
-------------------------
_space_group_symop.id  1    
_space_group_symop.operation_xyz   x,y,z              
_space_group_symop.operation_description 'identity mapping'

Mandatory attributes:
-----------
- operation_xyz

Optional attributes:
-----------
- id
- operation_description
- generator_xyz
- sg_id

Internal attributes (only for reading):
 - r, r_11, r_12, r_13, r_21, r_22, r_23, r_31, r_32, r_33 -- rotational matrix
 - b, b_1, b_2, b_3 -- translation vector
 
Methods:
---------
- get_symop_inversed(pcentr=numpy.zeros(3, dtype=float)) 
- get_symops_by_centring_type(centring_type)
- get_symops_by_generator_xyz(generator_xyz)
- get_coords_xyz_by_coord_xyz(coord_xyz)


    """
    MANDATORY_ATTRIBUTE = ("operation_xyz", )
    OPTIONAL_ATTRIBUTE = ("id", "operation_description", "generator_xyz", "sg_id")
    INTERNAL_ATTRIBUTE = ("r", "b", "r_11", "r_12", "r_13", "r_21", "r_22", "r_23", "r_31", "r_32", "r_33", "b_1", "b_2", "b_3")
    PREFIX = "space_group_symop"
    def __init__(self, id=None, operation_xyz=None, operation_description=None, generator_xyz=None, sg_id=None):
        super(SpaceGroupSymop, self).__init__(mandatory_attribute=self.MANDATORY_ATTRIBUTE, 
                                                optional_attribute=self.OPTIONAL_ATTRIBUTE, 
                                                internal_attribute=self.INTERNAL_ATTRIBUTE,
                                                prefix=self.PREFIX)
        self.id = id
        self.operation_xyz = operation_xyz
        self.operation_description = operation_description
        self.generator_xyz = generator_xyz
        self.sg_id = sg_id
        if self.is_defined:
            self.form_object

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

    @property
    def r_11(self)->Fraction:
        return getattr(self, "__r_11")
    @property
    def r_12(self)->Fraction:
        return getattr(self, "__r_12")
    @property
    def r_13(self)->Fraction:
        return getattr(self, "__r_13")
    @property
    def r_21(self)->Fraction:
        return getattr(self, "__r_21")
    @property
    def r_22(self)->Fraction:
        return getattr(self, "__r_22")
    @property
    def r_23(self)->Fraction:
        return getattr(self, "__r_23")
    @property
    def r_31(self)->Fraction:
        return getattr(self, "__r_31")
    @property
    def r_32(self)->Fraction:
        return getattr(self, "__r_32")
    @property
    def r_33(self)->Fraction:
        return getattr(self, "__r_33")
    @property
    def b_1(self)->Fraction:
        return getattr(self, "__b_1")
    @property
    def b_2(self):
        return getattr(self, "__b_2")
    @property
    def b_3(self)->Fraction:
        return getattr(self, "__b_3")
    @property
    def r(self):
        return getattr(self, "__r")
    @property
    def b(self):
        return getattr(self, "__b")

    @property
    def form_object(self)->bool:
        flag = True
        operation_xyz = self.operation_xyz
        if operation_xyz is None:
            return False
        r, b = CONSTANTS_AND_FUNCTIONS.transform_string_to_r_b(operation_xyz, labels=("x", "y", "z"))
        setattr(self, "__r", r)
        setattr(self, "__b", b)
        setattr(self, "__r_11", r[0, 0])
        setattr(self, "__r_12", r[0, 1])
        setattr(self, "__r_13", r[0, 2])
        setattr(self, "__r_21", r[1, 0])
        setattr(self, "__r_22", r[1, 1])
        setattr(self, "__r_23", r[1, 2])
        setattr(self, "__r_31", r[2, 0])
        setattr(self, "__r_32", r[2, 1])
        setattr(self, "__r_33", r[2, 2])
        setattr(self, "__b_1", b[0])
        setattr(self, "__b_2", b[1])
        setattr(self, "__b_3", b[2])
        return flag


    def get_symop_inversed(self, pcentr=numpy.zeros(3, float)):
        r, b = self.r, self.b
        r_new = -1*r
        b_new = -1*b+2*pcentr
        _symop = CONSTANTS_AND_FUNCTIONS.transform_r_b_to_string(r_new, b_new, labels=("x", "y", "z"))
        _item = SpaceGroupSymop(operation_xyz=_symop)
        return _item

    def get_symops_by_centring_type(self, centring_type:str)->Tuple[str]:
        r, b = self.r, self.b
        r_new = r
        shift = CONSTANTS_AND_FUNCTIONS.get_shift_by_centring_type(centring_type)
        symops = []
        for _shift in shift:
            b_new =  b + numpy.array(_shift, dtype=Fraction)
            _symop = CONSTANTS_AND_FUNCTIONS.transform_r_b_to_string(r_new, b_new, labels=("x", "y", "z"))
            symop = SpaceGroupSymop(operation_xyz=_symop)
            symops.append(symop)
        return symops

    def get_symops_by_generator_xyz(self, generator_xyz:str)->Tuple[str]:
        r, b = self.r, self.b
        symop = []
        W_g, w_g = CONSTANTS_AND_FUNCTIONS.transform_string_to_r_b(generator_xyz, labels=("x", "y", "z"))
        
        flag = True
        W_i = copy.deepcopy(r)
        w_i = copy.deepcopy(b)%1
        _symop = CONSTANTS_AND_FUNCTIONS.transform_r_b_to_string(W_i, w_i, labels=("x", "y", "z"))
        symop.append(_symop)
        i_stop = 1
        while flag:
            i_stop += 1
            w_i = (CONSTANTS_AND_FUNCTIONS.mult_matrix_vector(W_g, w_i) + w_g)%1 # w3 = W2 x w1 + w2
            W_i = CONSTANTS_AND_FUNCTIONS.mult_matrixes(W_g, W_i)          # W3 = W2 x W1
            _symop = CONSTANTS_AND_FUNCTIONS.transform_r_b_to_string(W_i, w_i, labels=("x", "y", "z"))
            if ((_symop in symop) | (i_stop > 100)):
                flag = False
            else:
                symop.append(_symop)
        res = []
        gen_orig = ""
        if not(self.generator_xyz is None): gen_orig = f"{self.generator_xyz:}"
        res =[SpaceGroupSymop(id=f"{_i+1:}", operation_xyz=_symop, generator_xyz=f"{gen_orig:}_{generator_xyz:}") for _i, _symop in enumerate(symop)]
        return res

    def get_coords_xyz_by_coord_xyz(self, coord_xyz:str)->Tuple:
        _item = SpaceGroupSymop(operation_xyz=coord_xyz)
        generator_xyz = self.operation_xyz
        res = [_.operation_xyz for _ in _item.get_symops_by_generator_xyz(generator_xyz)]
        return tuple(res)

class SpaceGroupSymopL(LoopConstr):
    """
SpaceGroupSymopL
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
- sg_id


Class methods:
---------
- create_by_generators_xyz
- get_coords_xyz_by_coord_xyz


reference: https://www.iucr.org/__data/iucr/cifdic_html/2/cif_sym.dic/Cspace_group_symop.html
    """
    CATEGORY_KEY = ("id", )
    ITEM_CLASS = SpaceGroupSymop
    def __init__(self, item = [], loop_name=""):
        super(SpaceGroupSymopL, self).__init__(category_key=self.CATEGORY_KEY, item_class=self.ITEM_CLASS, loop_name=loop_name)
        self.item = item

    @classmethod
    def create_by_generators_xyz(cls, generators_xyz:List):
        _item = SpaceGroupSymop(operation_xyz="x,y,z")
        items_1 = [_item]
        for generator_xyz in generators_xyz:
            items_2 = []
            for _item in items_1:
                items_2.extend(_item.get_symops_by_generator_xyz(generator_xyz))
            items_1 = items_2
        symop = []
        item_uniq=[]
        _i = 0
        for _item in items_1:
            if not(_item.operation_xyz in symop):
                _i += 1
                symop.append(_item.operation_xyz)
                _item.id=f"{_i:}"
                item_uniq.append(_item)
        _obj = cls(item=item_uniq)
        return _obj


    def get_coords_xyz_by_coord_xyz(self, coord_xyz:str)->Tuple:
        item = self.item
        res = []
        for _item in item:
            res.extend(_item.get_coords_xyz_by_coord_xyz(coord_xyz))
        return tuple(frozenset(res))

