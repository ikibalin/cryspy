__author__ = 'ikibalin'
__version__ = "2019_11_26"

import os
import numpy
from pycifstar import Global


import warnings

from typing import List, Tuple
from cryspy.common.cl_item_constr import ItemConstr

import cryspy.symcif.CONSTANTS_AND_FUNCTIONS as CONSTANTS_AND_FUNCTIONS
from cryspy.symcif.cl_space_group_symop import SpaceGroupSymop, SpaceGroupSymopL
from cryspy.symcif.cl_space_group_wyckoff import SpaceGroupWyckoff, SpaceGroupWyckoffL


class SpaceGroup(ItemConstr):
    """
Contains all the data items that refer to the space group as a
whole, such as its name, Laue group etc. 

Space-group types are identified by their number as listed in
International Tables for Crystallography Volume A, or by their
Schoenflies symbol. Specific settings of the space groups can
be identified by their Hall symbol, by specifying their
symmetry operations or generators, or by giving the
transformation that relates the specific setting to the
reference setting based on International Tables Volume A and
stored in this dictionary.

The commonly used Hermann-Mauguin symbol determines the
space-group type uniquely but several different Hermann-Mauguin
symbols may refer to the same space-group type. A
Hermann-Mauguin symbol contains information on the choice of
the basis, but not on the choice of origin.

:Ref: International Tables for Crystallography (2002). Volume A,
      Space-group symmetry, edited by Th. Hahn, 5th ed.
      Dordrecht: Kluwer Academic Publishers.

Description in cif file::

    _space_group.id                    1
    _space_group.name_H-M_ref         'C 2/c'
    _space_group.name_Schoenflies      C2h.6
    _space_group.IT_number             15
    _space_group.name_Hall           '-C 2yc'
    _space_group.Bravais_type          mS
    _space_group.Laue_class            2/m
    _space_group.crystal_system        monoclinic
    _space_group.centring_type         C
    _space_group.Patterson_name_H-M  'C 2/m'

:Mandatory attribute: no

:Optional attribute: - name_hm_alt
                     - it_number
                     - it_coordinate_system_code
                     - id
                     - bravais_type
                     - laue_class
                     - patterson_name_hm
                     - centring_type
                     - crystal_system
                     - name_hm_alt_description
                     - name_hm_full
                     - name_hm_ref
                     - name_hall
                     - name_hall_schoenflies
                     - point_group_hm
                     - reference_setting
                     - transform_pp_abc
                     - transform_qq_xyz

:Methods: get_symop

`<https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Cspace_group.html>`_
    """
    MANDATORY_ATTRIBUTE = ()
    OPTIONAL_ATTRIBUTE = ("id", "name_hm_alt", "name_hm_alt_description", "it_number", "name_hm_ref", "it_coordinate_system_code", "name_hm_full", "name_hall", "name_schoenflies",
     "point_group_hm", "laue_class", "patterson_name_hm", 
    "centring_type", "bravais_type", "crystal_system", "reference_setting", "transform_pp_abc", "transform_qq_xyz")
    RELATED_CIF_MANDATORY_ATTRIBUTE = ()
    RELATED_CIF_OPTIONAL_ATTRIBUTE = ("id", "name_H-M_alt", "name_H-M_alt_description", "IT_number", "name_H-M_ref", "IT_coordinate_system_code", "name_H-M_full", "name_Hall", "name_Schoenflies",
     "point_group_H-M", "Laue_class", "Patterson_name_H-M", 
    "centring_type", "Bravais_type", "crystal_system", "reference_setting", "transform_pp_abc", "transform_qq_xyz")
    INTERNAL_ATTRIBUTE = ("centrosymmetry", "pcentr", "reduced_space_group_symop", "full_space_group_symop", "space_group_wyckoff",
                          "shift")
    PREFIX = "space_group"
    ACCESIBLE_IT_NUMBER = CONSTANTS_AND_FUNCTIONS.ACCESIBLE_IT_NUMBER
    ACCESIBLE_BRAVAIS_TYPE = CONSTANTS_AND_FUNCTIONS.ACCESIBLE_BRAVAIS_TYPE
    ACCESIBLE_IT_COORDINATE_SYSTEM_CODE = CONSTANTS_AND_FUNCTIONS.ACCESIBLE_IT_COORDINATE_SYSTEM_CODE
    ACCESIBLE_LAUE_CLASS = CONSTANTS_AND_FUNCTIONS.ACCESIBLE_LAUE_CLASS
    ACCESIBLE_CENTRING_TYPE = CONSTANTS_AND_FUNCTIONS.ACCESIBLE_CENTRING_TYPE
    ACCESIBLE_CRYSTAL_SYSTEM = CONSTANTS_AND_FUNCTIONS.ACCESIBLE_CRYSTAL_SYSTEM
    ACCESIBLE_NAME_HM_REF = CONSTANTS_AND_FUNCTIONS.ACCESIBLE_NAME_HM_SHORT
    ACCESIBLE_NAME_HM_ALT = frozenset(set(CONSTANTS_AND_FUNCTIONS.ACCESIBLE_NAME_HM_FULL).union(set(CONSTANTS_AND_FUNCTIONS.ACCESIBLE_NAME_HM_EXTENDED)))
    ACCESIBLE_NAME_SCHOENFLIES = CONSTANTS_AND_FUNCTIONS.ACCESIBLE_NAME_SCHOENFLIES
    ACCESIBLE_NAME_HALL = CONSTANTS_AND_FUNCTIONS.ACCESIBLE_NAME_HALL_SHORT
    ACCESIBLE_REFERENCE_SETTING = CONSTANTS_AND_FUNCTIONS.ACCESIBLE_REFERENCE_SETTING
    
    def __init__(self, name_hm_alt=None, it_number=None, it_coordinate_system_code=None, id=None, bravais_type=None,
    laue_class=None, patterson_name_hm=None, centring_type=None, crystal_system=None, name_hm_alt_description=None, 
    name_hm_full=None, name_hm_ref=None, name_hall=None, name_schoenflies=None, point_group_hm=None,
    reference_setting=None, transform_pp_abc=None, transform_qq_xyz=None):
        
        super(SpaceGroup, self).__init__(mandatory_attribute=self.MANDATORY_ATTRIBUTE, 
                                         optional_attribute=self.OPTIONAL_ATTRIBUTE, 
                                         internal_attribute=self.INTERNAL_ATTRIBUTE,
                                         prefix=self.PREFIX)

        
        self.name_hm_alt = name_hm_alt
        self.it_number = it_number
        self.it_coordinate_system_code = it_coordinate_system_code
        self.id = id
        self.bravais_type = bravais_type
        self.laue_class = laue_class
        self.patterson_name_hm = patterson_name_hm
        self.centring_type = centring_type
        self.crystal_system = crystal_system
        self.name_hm_alt_description = name_hm_alt_description
        self.name_hm_full = name_hm_full
        self.name_hm_ref = name_hm_ref
        self.name_hall = name_hall
        self.name_schoenflies = name_schoenflies
        self.point_group_hm = point_group_hm
        self.reference_setting = reference_setting
        self.transform_pp_abc = transform_pp_abc
        self.transform_qq_xyz = transform_qq_xyz
        if self.is_defined:
            self.form_object

        
    @property
    def name_hm_alt(self) -> str:
        """
_space_group.name_H-M_alt allows for an alternative Hermann–
Mauguin symbol to be given. The way in which this item is used
is determined by the user and should be described in the item
_space_group.name_H-M_alt_description. It may, for exam-
ple, be used to give one of the extended Hermann–Mauguin
symbols given in Table 4.3.2.1 of International Tables for
Crystallography Volume A (2002) or a full Hermann–Mauguin
symbol for an unconventional setting. Each component of the
space-group name is separated by a space or an underscore
character. The use of a space is strongly recommended. The
underscore is only retained because it was used in older CIFs.
It should not be used in new CIFs. Subscripts should appear
without special symbols. Bars should be given as negative
signs before the numbers to which they apply. The commonly
used Hermann–Mauguin symbol determines the space-group
type uniquely, but a given space-group type may be described
by more than one Hermann–Mauguin symbol. The space-
group type is best described using _space_group.IT_number or
_space_group.name_Schoenflies . The Hermann–Mauguin sym-
bol may contain information on the choice of basis but does not
contain information on the choice of origin. To define the setting
uniquely, use _space_group.name_Hall , list the symmetry oper-
ations or generators, or give the transformation that relates the
setting to the reference setting defined in this dictionary under
_space_group.reference_setting.

:Reference: International Tables for Crystallography (2002).
            Volume A, Space-group symmetry, edited by Th. Hahn, 5th ed.
            Dordrecht: Kluwer Academic Publishers.

(full tuple stored in SpaceGroup.ACCESIBLE_NAME_HM_ALT)
        """
        return getattr(self, "__name_hm_alt")
    @name_hm_alt.setter
    def name_hm_alt(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = str(x)
            if not(x_in in self.ACCESIBLE_NAME_HM_ALT):
                warnings.warn(f"name_hm_alt type '{x_in:}' is not supported", UserWarning, stacklevel=2)
                x_in = None
        setattr(self, "__name_hm_alt", x_in)

    @property
    def it_number(self) -> int:
        """
The number as assigned in International Tables for Crystal-
lography Volume A, specifying the proper affine class (i.e. the
orientation-preserving affine class) of space groups (crystallo-
graphic space-group type) to which the space group belongs. This
number defines the space-group type but not the coordinate system
in which it is expressed.

:Reference: International Tables for Crystallography (2002).
            Volume A, Space-group symmetry, edited by Th. Hahn, 5th ed.
            Dordrecht: Kluwer Academic Publishers.

The permitted range is (1,230)
        """
        return getattr(self, "__it_number")
    @it_number.setter
    def it_number(self, x):

        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = int(x)
            if not((x_in >=1)&(x_in <=230)):
                warnings.warn(f"it_number '{x_in:}' is not supported", UserWarning, stacklevel=2)
                x_in = None              
        setattr(self, "__it_number", x_in)

    @property
    def it_coordinate_system_code(self) -> str:
        """
A qualifier taken from the enumeration list identifying which set-
ting in International Tables for Crystallography Volume A (2002)
(IT) is used. See IT Table 4.3.2.1, Section 2.2.16, Table 2.2.16.1,
Section 2.2.16.1 and Fig. 2.2.6.4. This item is not computer-
interpretable and cannot be used to define the coordinate system.

Use _space_group.transform_ * instead.

:Reference: International Tables for Crystallography (2002).
            Volume A, Space-group symmetry, edited by Th. Hahn, 5th ed.
            Dordrecht: Kluwer Academic Publishers.

The data value must be one of the following:

   b1, b2, b3, -b1, -b2, -b3, c1, c2, c3, 
   -c1, -c2, -c3, a1, a2, a3, -a1, -a2, -a3, 
   abc, ba-c, cab, -cba, bca, a-cb, 1abc, 1ba-c, 1cab, 1-cba, 1bca, 1a-cb, 
   2abc, 2ba-c, 2cab, 2-cba, 2bca, 2a-cb, 1, 2, h, r, 
        """
        return getattr(self, "__it_coordinate_system_code")
    @it_coordinate_system_code.setter
    def it_coordinate_system_code(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            if not(x in self.ACCESIBLE_IT_COORDINATE_SYSTEM_CODE):
                warnings.warn(f"it_coordinate_system_code '{x_in:}' is not supported", UserWarning, stacklevel=2)
                x_in = None            
            x_in = str(x)
        setattr(self, "__it_coordinate_system_code", x_in)

    @property
    def id(self) -> str:
        """
This is an identifier needed if _space_group .* items are looped.

:TODO: looped _space_group is not introduced
        """
        return getattr(self, "__id")
    @id.setter
    def id(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = str(x)
        setattr(self, "__id", x_in)

    @property
    def bravais_type(self) -> str:
        """
The symbol denoting the lattice type (Bravais type) to which the
translational subgroup (vector lattice) of the space group belongs.
It consists of a lower-case letter indicating the crystal system
followed by an upper-case letter indicating the lattice centring.

The setting-independent symbol mS replaces the setting-dependent
symbols mB and mC, and the setting-independent symbol oS
replaces the setting-dependent symbols oA, oB and oC.

:Reference: International Tables for Crystallography (2002).
            Volume A, Space-group symmetry, edited by Th. Hahn, 5th ed.,
            p. 15. Dordrecht: Kluwer Academic Publishers.

The data value must be one of the following:

    aP
    mP mS
    oP oS oI oF
    tP tI
    hP hR
    cP cI cF
        """
        return getattr(self, "__bravais_type")
    @bravais_type.setter
    def bravais_type(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = "".join(str(x).strip().split())
            if not(x_in in self.ACCESIBLE_BRAVAIS_TYPE):
                warnings.warn(f"bravais_type '{x_in:}' is not supported", UserWarning, stacklevel=2)
                x_in = None
        setattr(self, "__bravais_type", x_in)


    @property
    def laue_class(self) -> str:
        """
The Hermann–Mauguin symbol of the geometric crystal class of
the point group of the space group where a centre of inversion is
added if not already present.

The data value must be one of the following:

    -1, 2/m, mmm, 4/m, 4/mmm, 
    -3, -3m, 6/m, 6/mmm, m-3, m-3m
"""
        return getattr(self, "__laue_class")
    @laue_class.setter
    def laue_class(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = str(x)
            if not(x_in in self.ACCESIBLE_LAUE_CLASS):
                warnings.warn(f"laue_class '{x_in:}' is not supported", UserWarning, stacklevel=2)
                x_in = None            
        setattr(self, "__laue_class", x_in)

    @property
    def patterson_name_hm(self) -> str:
        """
The Hermann–Mauguin symbol of the type of that centrosym-
metric symmorphic space group to which the Patterson function
belongs; see Table 2.2.5.1 in International Tables for Crystallog-
raphy Volume A (2002). A space separates each symbol referring
to different axes. Underscores may replace the spaces, but this use
is discouraged. Subscripts should appear without special symbols.
Bars should be given as negative signs before the number to which
they apply.

:Reference: International Tables for Crystallography (2002).
            Volume A, Space-group symmetry, edited by Th. Hahn, 5th ed.,
            Table 2.2.5.1. Dordrecht: Kluwer Academic Publishers.

:Examples: ‘P -1’, ‘P 2/m’, ‘C 2/m’, ‘P m m m’, ‘C m m m’, ‘I m m m’, ‘F m m m’,
           ‘P 4/m’, ‘I 4/m’, ‘P 4/m m m’, ‘I 4/m m m’, ‘P -3’, ‘R -3’, ‘P -3 m 1’,
           ‘R -3 m’, ‘P -3 1 m’, ‘P 6/m’, ‘P 6/m m m’, ‘P m -3’, ‘I m -3’, ‘F m -3’,
           ‘P m -3 m’, ‘I m -3 m’, ‘F m -3 m’.
        """
        return getattr(self, "__patterson_name_hm")
    @patterson_name_hm.setter
    def patterson_name_hm(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = str(x).strip()
        setattr(self, "__patterson_name_hm", x_in)

    @property
    def centring_type(self) -> str:
        """
Symbol for the lattice centring. This symbol may be dependent on
the coordinate system chosen.

The data value must be one of the following:

    :"P": primitive no centring
    :"A": A-face centred (0, 1/2, 1/2)
    :"B": B-face centred (1/2, 0, 1/2)
    :"C": C-face centred (1/2, 1/2, 0)
    :"F": all faces centred (0, 1/2, 1/2), (1/2, 0, 1/2), (1/2, 1/2, 0)
    :"I": body centred (1/2, 1/2, 1/2)
    :"R": rhombohedral obverse centred (2/3, 1/3, 1/3), (1/3, 2/3, 2/3)
    :"Rrev": rhombohedral reverse centred (1/3, 2/3, 1/3), (2/3, 1/3, 2/3)
    :"H": hexagonal centred (2/3, 1/3, 0), (1/3, 2/3, 0)
        """
        return getattr(self, "__centring_type")
    @centring_type.setter
    def centring_type(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = str(x)
            if not(x_in in self.ACCESIBLE_CENTRING_TYPE):
                warnings.warn(f"centring_type '{x_in:}' is not supported", UserWarning, stacklevel=2)
                x_in = None            
        setattr(self, "__centring_type", x_in)

    @property
    def crystal_system(self) -> str:
        """
The name of the system of geometric crystal classes of space
groups (crystal system) to which the space group belongs. Note
that crystals with the hR lattice type belong to the trigonal system.

The data value must be one of the following:

    triclinic, monoclinic, orthorhombic, 
    tetragonal, trigonal, hexagonal, cubic   
    """
        return getattr(self, "__crystal_system")
    @crystal_system.setter
    def crystal_system(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = str(x)
            if not(x_in in self.ACCESIBLE_CRYSTAL_SYSTEM):
                warnings.warn(f"crystal_system '{x_in:}' is not supported", UserWarning, stacklevel=2)
                x_in = None            
        setattr(self, "__crystal_system", x_in)

    @property
    def name_hm_alt_description(self) -> str:
        """
A free-text description of the code appearing in _space_
group.name_H-M_alt
        """
        return getattr(self, "__name_hm_alt_description")
    @name_hm_alt_description.setter
    def name_hm_alt_description(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = str(x)
        setattr(self, "__name_hm_alt_description", x_in)

    @property
    def name_hm_full(self) -> str:
        """
The full international Hermann–Mauguin space-group symbol as
defined in Section 2.2.3 and given as the second item of the sec-
ondlineof eachof the space-grouptablesof Part 7of International
Tables for Crystallography Volume A (2002). Each component of
the space-group name is separated by a space or an underscore
character. The use of a space is strongly recommended. The under-
score is only retained because it was used in old CIFs. It should
not be used in new CIFs. Subscripts should appear without special
symbols. Bars should be given as negative signs before the num-
bers to which they apply. The commonly used Hermann–Mauguin
symbol determines the space-group type uniquely but a given
space-group type may be described by more than one Hermann–
Mauguin symbol. The space-group type is best described using
_space_group.IT_number or _space_group.name_Schoenflies .
The full international Hermann–Mauguin symbol contains infor-
mation about the choice of basis for monoclinic and orthorhombic
space groups but does not give information about the choice of ori-
gin. To define the setting uniquely use _space_group.name_Hall ,
list the symmetry operations or generators, or give the transforma-
tion relating the setting used to the reference setting defined in this
dictionary under _space_group.reference_setting .
Reference: International Tables for Crystallography (2002).
Volume A, Space-group symmetry, edited by Th. Hahn, 5th ed.
Dordrecht: Kluwer Academic Publishers.

:Example: ‘P 21/n 21/m 21/a’ (full symbol for Pnma)
        """
        return getattr(self, "__name_hm_full")
    @name_hm_full.setter
    def name_hm_full(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = str(x)
        setattr(self, "__name_hm_full", x_in)

    @property
    def name_hm_ref(self) -> str:
        """
The short international Hermann–Mauguin space-group sym-
bol as defined in Section 2.2.3 and given as the first item of
each space-group table in Part 7 of International Tables for
Crystallography Volume A (2002). Each component of the space-
group name is separated by a space or an underscore charac-
ter. The use of a space is strongly recommended. The under-
score is only retained because it was used in old CIFs. It
should not be used in new CIFs. Subscripts should appear with-
out special symbols. Bars should be given as negative signs
before the numbers to which they apply. The short international
Hermann–Mauguin symbol determines the space-group type
uniquely. However, the space-group type is better described using
_space_group.IT_number or _space_group.name_Schoenflies .
The short international Hermann–Mauguin symbol contains
no information on the choice of basis or origin. To define
the setting uniquely use _space_group.name_Hall , list the
symmetry operations or generators, or give the transforma-
tion that relates the setting to the reference setting defined
in this dictionary under _space_group.reference_setting .
_space_group.name_H-M_alt may be used to give the Hermann–
Mauguin symbol corresponding to the setting used. In the enumer-
ation list below, each possible value is identified by space-group
number and Schoenflies symbol.

:Reference: International Tables for Crystallography (2002).
            Volume A, Space-group symmetry, edited by Th. Hahn, 5th ed.
            Dordrecht: Kluwer Academic Publishers.

The data value must be one of the following:

    ’P 1’, ’P 2’, ...
    (full tuple stored in SpaceGroup.ACCESIBLE_NAME_HM_REF)

        """
        return getattr(self, "__name_hm_ref")
    @name_hm_ref.setter
    def name_hm_ref(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = str(x)
            if not(x_in in self.ACCESIBLE_NAME_HM_REF):
                warnings.warn(f"name_hm_ref '{x_in:}' is not supported", UserWarning, stacklevel=2)
                x_in = None            
        setattr(self, "__name_hm_ref", x_in)

    @property
    def name_hall(self) -> str:
        """
Space-group symbol defined by Hall. _space_group.name_Hall
uniquely defines the space group and its reference to a particu-
lar coordinate system. Each component of the space-group name
is separated by a space or an underscore character. The use of a
space is strongly recommended. The underscore is only retained
because it was used in old CIFs. It should not be used in new CIFs.

:References: Hall, S. R. (1981). Acta Cryst. A37, 517–525; erra-
             tum (1981), A37, 921. International Tables for Crystallography
             (2001). Volume B, Reciprocal space, edited by U. Shmueli, 2nd
             ed., Appendix 1.4.2. Dordrecht: Kluwer Academic Publishers.

:Examples: ‘P 2c -2ac’, ‘-I 4bd 2ab 3’, ...
            (full tuple stored in SpaceGroup.ACCESIBLE_NAME_HALL)

        """
        return getattr(self, "__name_hall")
    @name_hall.setter
    def name_hall(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = str(x)
            if not(x_in in self.ACCESIBLE_NAME_HALL):
                warnings.warn(f"name_hall '{x_in:}' is not supported", UserWarning, stacklevel=2)
                x_in = None            
        setattr(self, "__name_hall", x_in)

    @property
    def name_schoenflies(self) -> str:
        """
The Schoenflies symbol as listed in International Tables for
Crystallography Volume A denoting the proper affine class (i.e.
orientation-preserving affine class) of space groups (space-group
type) to which the space group belongs. This symbol defines the
space-group type independently of the coordinate system in which
the space group is expressed. The symbol is given with a period,
‘.’, separating the Schoenflies point group and the superscript.

:Reference: International Tables for Crystallography (2002).
            Volume A, Space-group symmetry, edited by Th. Hahn, 5th ed.
            Dordrecht: Kluwer Academic Publishers.

The data value must be one of the following:

    C1.1, Ci.1, ...
    (full tuple stored in SpaceGroup.ACCESIBLE_NAME_SCHOENFLIES)
        """
        return getattr(self, "__name_schoenflies")
    @name_schoenflies.setter
    def name_schoenflies(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = str(x)
            if not(x_in in self.ACCESIBLE_NAME_SCHOENFLIES):
                warnings.warn(f"name_schoenflies '{x_in:}' is not supported", UserWarning, stacklevel=2)
                x_in = None            
        setattr(self, "__name_schoenflies", x_in)

    @property
    def point_group_hm(self) -> str:
        """
The Hermann–Mauguin symbol denoting the geometric crystal
class of space groups to which the space group belongs, and the
geometric crystal class of point groups to which the point group of
the space group belongs.

:Examples: ‘-4’, ‘4/m’. 
        """
        return getattr(self, "__point_group_hm")
    @point_group_hm.setter
    def point_group_hm(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = str(x)
        setattr(self, "__point_group_hm", x_in)

    @property
    def reference_setting(self) -> str:
        """
The reference setting of a given space group is the setting chosen
by the International Union of Crystallography as a unique setting
to which other settings can be referred using the transformation
matrix column pair given in _space_group.transform_Pp_abc
and _space_group.transform_Qq_xyz . The settings are given
in the enumeration list in the form ‘ _space_group.IT_number:
_space_group.name_Hall ’. The space-group number defines the
space-group type and the Hall symbol specifies the symmetry gen-
erators referred to the reference coordinate system. The 230 refer-
ence settings chosen are identical to the settings listed in Inter-
national Tables for Crystallography Volume A (2002). For the
space groups where more than one setting is given in Interna-
tional Tables, the following choices have been made. For mono-
clinic space groups: unique axis b and cell choice 1. For space
groups with two origins: origin choice 2 (origin at inversion cen-
tre, indicated by adding :2 to the Hermann–Mauguin symbol in
the enumeration list). For rhombohedral space groups: hexago-
nal axes (indicated by adding :h to the Hermann–Mauguin sym-
bol in the enumeration list). Based on the symmetry table of
R. W. Grosse-Kunstleve, ETH, Zurich. The enumeration list may
be extracted from the dictionary and stored as a separate CIF that
can be referred to as required.

In the enumeration list below, each reference setting is identi-
fied by Schoenfliessymbol and by the Hermann–Mauguinsymbol,
augmented by :2 or :h suffixes as described above.

:References: International Tables for Crystallography (2002).
             Volume A, Space-group symmetry, edited by Th. Hahn, 5th
             ed. Dordrecht: Kluwer Academic Publishers. Grosse-Kunstleve,
             R. W. (2001). Xtal System of Crystallographic Programs, Sys-
             tem Documentation. http://xtal.crystal.uwa.edu.au/man/xtal3.7-
             228.html (or follow links to D

The data value must be one of the following:

    ’001:P 1’, ’002:-P 1’, ... 
    (full tuple stored in SpaceGroup.ACCESIBLE_REFERENCE_SETTING)
        """
        return getattr(self, "__reference_setting")
    @reference_setting.setter
    def reference_setting(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = str(x)
            if not(x_in in self.ACCESIBLE_REFERENCE_SETTING):
                warnings.warn(f"reference_setting '{x_in:}' is not supported", UserWarning, stacklevel=2)
                x_in = None            
        setattr(self, "__reference_setting", x_in)

    @property
    def transform_pp_abc(self) -> str:
        """
This item specifies the transformation (P,p) of the basis vectors
from the setting used in the CIF (a,b,c) to the reference setting
given in _space_group.reference_setting (a' ,b' ,c'). The value
is given in Jones–Faithful notation corresponding to the rotational
matrix P combined with the origin shift vector p in the expression

.. math::

   (a', b', c') = (a, b, c) P + p.

P is a post-multiplication matrix of a row (a,b,c) of column vec-
tors. It is related to the inverse transformation (Q,q) by

.. math::

   P = Q^{−1},

.. math::

   p = Pq = − (Q^{−1}) q.

These transformations are applied as follows: atomic coordi-
nates :math:`(x', y', z') = Q (x,y,z) + q`, Miller indices 
:math:`(h', k', l') = (h, k, l) P`, 
symmetry operations
:math:`W' = (Q,q) W (P,p)`, basis vectors

.. math::

   (a', b', c') = (a,b,c) P + p.

This item is given as a character string involving the characters
a, b and c with commas separating the expressions for the a', b'
and c' vectors. The numeric values may be given as integers, frac-
tions or real numbers. Multiplication is implicit, division must be
explicit. White space within the string is optional.

:Examples: - ‘-b+c, a+c, -a+b+c’ (R3:r to R3:h), 
           - ‘a-1/4, b-1/4, c-1/4’ (Pnnn:1 to Pnnn:2), 
           - ‘b-1/2, c-1/2, a-1/2’ (Bbab:1 toCcca:2). 
        """
        return getattr(self, "__transform_pp_abc")
    @transform_pp_abc.setter
    def transform_pp_abc(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = str(x)
        setattr(self, "__transform_pp_abc", x_in)

    @property
    def transform_qq_xyz(self) -> str:
        """
This item specifies the transformation (Q,q) of the atomic
coordinates from the setting used in the CIF [(x, y, z) referred
to the basis vectors (a, b, c)] to the reference setting given
in _space_group.reference_setting [(x', y', z') referred to the
basis vectors (a', b', c')]. The value given in Jones–Faithful nota-
tion corresponds to the rotational matrix Q combined with the ori-
gin shift vector q in the expression

.. math::
  
   (x', y', z') = Q(x,y,z) + q.

Q is a pre-multiplication matrix of the column vector (x,y,z). It is
related to the inverse transformation (P,p) by

.. math::
   
   P = Q^{−1} ,

.. math::
   
   p = Pq = − (Q^{−1}) q,

where the P and Q transformations are applied as follows: atomic
coordinates :math:`(x', y', z') = Q(x,y,z)+q`, Miller indices 
:math:`(h', k', l') = (h,k,l)P`, 
symmetry operations 
:math:`W' = (Q, q) W (P, p)`, 
basis vectors

.. math::

   (a', b', c') = (a, b, c) P + p.

This item is given as a character string involving the characters
math:`x`, :math:`y` and :math:`z` with commas separating the expressions 
for the :math:`x'`, :math:`y'`
and :math:`z'` components. The numeric values may be given as integers,
fractions or real numbers. Multiplication is implicit, division must
be explicit. White space within the string is optional.

:Examples: - ‘-x/3+2y/3-z/3, -2x/3+y/3+z/3, x/3+y/3+z/3’ (R3:r to R3:h),
           - ‘x+1/4,y+1/4,z+1/4’ (Pnnn:1 to Pnnn:2), 
           - ‘z+1/2,x+1/2,y+1/2’ (Bbab:1 to Ccca:2).
        """
        return getattr(self, "__transform_qq_xyz")
    @transform_qq_xyz.setter
    def transform_qq_xyz(self, x):
        if ((x is None) | (x == ".")):
            x_in = None
        else:
            x_in = str(x)
        setattr(self, "__transform_qq_xyz", x_in)

    @property
    def centrosymmetry(self)->bool:
        return getattr(self, "__centrosymmetry")

    @property
    def pcentr(self):
        return getattr(self, "__pcentr")

    @property
    def reduced_space_group_symop(self):
        return getattr(self, "__reduced_space_group_symop")

    @property
    def full_space_group_symop(self):
        return getattr(self, "__full_space_group_symop")

    @property
    def space_group_wyckoff(self):
        return getattr(self, "__space_group_wyckoff")

    @property
    def shift(self):
        return getattr(self, "__shift")

    @property
    def is_defined(self)->bool:
        flag =  self.define_it_number_and_it_coordinate_system_code()
        return flag 
    @property
    def form_object(self)->bool:
        flag = self.form_object_by_it_number_it_coordinate_system_code()
        setattr(self, "__shift", CONSTANTS_AND_FUNCTIONS.get_shift_by_centring_type(self.centring_type))
        return flag

    def define_it_number_and_it_coordinate_system_code(self):
        flag = True
        it_coordinate_system_code = self.it_coordinate_system_code
        if (self.is_defined_attribute("it_number") & self.is_defined_attribute("it_coordinate_system_code")):
            it_number = self.it_number
            it_coordinate_system_code = self.it_coordinate_system_code
        elif self.is_defined_attribute("name_hm_alt"):
            _res = CONSTANTS_AND_FUNCTIONS.get_type_hm(self.name_hm_alt)
            if "extended" in _res:
                it_number, it_coordinate_system_codes = CONSTANTS_AND_FUNCTIONS.get_it_number_it_coordinate_system_codes_by_name_hm_extended(self.name_hm_alt)
                if  (not(it_coordinate_system_code in it_coordinate_system_codes)):
                    it_coordinate_system_code = CONSTANTS_AND_FUNCTIONS.auto_choose_it_coordinate_system_code(it_number, it_coordinate_system_codes)
            elif "full" in _res:
                it_number = CONSTANTS_AND_FUNCTIONS.get_it_number_by_name_hm_full(self.name_hm_alt)
        elif self.is_defined_attribute("name_hm_ref"):
            it_number = CONSTANTS_AND_FUNCTIONS.get_it_number_by_name_hm_short(self.name_hm_ref)
        elif self.is_defined_attribute("name_hall"):
            it_number = CONSTANTS_AND_FUNCTIONS.get_it_number_by_name_hall(self.name_hall)
        elif self.is_defined_attribute("it_number"):
            it_number = self.it_number
            it_coordinate_system_code = self.it_coordinate_system_code
        elif self.is_defined_attribute("name_schoenflies"):
            it_number = CONSTANTS_AND_FUNCTIONS.get_it_number_by_name_schoenflies(self.name_schoenflies)
        else:
            flag = False

        if flag:
            it_coordinate_system_codes = CONSTANTS_AND_FUNCTIONS.get_it_coordinate_system_codes_by_it_number(it_number)
            setattr(self, "it_number", it_number)
            if  (not(it_coordinate_system_code in it_coordinate_system_codes)):
                it_coordinate_system_code = CONSTANTS_AND_FUNCTIONS.auto_choose_it_coordinate_system_code(it_number, it_coordinate_system_codes)
            
            setattr(self, "it_coordinate_system_code", it_coordinate_system_code)
        return flag

    def form_object_by_it_number_it_coordinate_system_code(self):
        """
        TODO: Solve the problem with centring_type for hexagonal systems (Ex.: 166 spcace group)
        """
        bravais_type, laue_class, patterson_name_hm, centring_type, crystal_system = None, None, None, None, None
        name_hm_extended, name_hm_full, name_hm_short = None, None, None
        name_hall, name_schoenflies, point_group_hm = None, None, None
        pcentr = None
        reduced_space_group_symop, full_space_group_symop = None, None
        space_group_wyckoff = None

        lattice_type = None
        generators = ()        
        flag = self.define_it_number_and_it_coordinate_system_code()
        it_number = self.it_number
        it_c_s_c = self.it_coordinate_system_code
        crystal_system = CONSTANTS_AND_FUNCTIONS.get_crystal_system_by_it_number(it_number)
        name_hm_extended = CONSTANTS_AND_FUNCTIONS.get_name_hm_extended_by_it_number_it_coordinate_system_code(it_number, it_c_s_c)
        if (name_hm_extended is not None): 
            centring_type = CONSTANTS_AND_FUNCTIONS.get_centring_type_by_name_hm_extended(name_hm_extended)
        if ((centring_type is not None) & (crystal_system is not None)): 
            bravais_type = CONSTANTS_AND_FUNCTIONS.get_bravais_type_by_centring_type_crystal_system(centring_type, crystal_system)              
        name_hm_short = CONSTANTS_AND_FUNCTIONS.get_name_hm_short_by_it_number(it_number)
        if (name_hm_short is not None):
            lattice_type = CONSTANTS_AND_FUNCTIONS.get_lattice_type_by_name_hm_short(name_hm_short)
        name_hm_full = CONSTANTS_AND_FUNCTIONS.get_name_hm_full_by_it_number(it_number)
        name_hall = CONSTANTS_AND_FUNCTIONS.get_name_hall_by_it_number(it_number)
        if name_hall is not None:
            centrosymmetry = CONSTANTS_AND_FUNCTIONS.get_centrosymmetry_by_name_hall(name_hall)
        name_schoenflies = CONSTANTS_AND_FUNCTIONS.get_name_schoenflies_by_it_number(it_number)
        if name_schoenflies is not None:
            laue_class = CONSTANTS_AND_FUNCTIONS.get_laue_class_by_name_schoenflies(name_schoenflies)
            point_group_hm = CONSTANTS_AND_FUNCTIONS.get_point_group_hm_short_by_name_schoenflies(name_schoenflies)
            if point_group_hm is not None:
                generators = CONSTANTS_AND_FUNCTIONS.get_generators_by_point_group_hm(point_group_hm)
            if ((lattice_type is not None) & (laue_class is not None)):
                patterson_name_hm = CONSTANTS_AND_FUNCTIONS.get_patterson_name_hm_by_lattice_type_laue_class(lattice_type, laue_class)

        pcentr, reduced_space_group_symop, full_space_group_symop, space_group_wyckoff = self.form_reduced_full_space_group_symop_space_group_wyckoff(it_number, it_c_s_c, centrosymmetry, centring_type)

        if pcentr is not None:
            setattr(self, "__pcentr", pcentr)
        if reduced_space_group_symop is not None:
            setattr(self, "__reduced_space_group_symop", reduced_space_group_symop)
        if full_space_group_symop is not None:
            setattr(self, "__full_space_group_symop", full_space_group_symop)
        if space_group_wyckoff is not None:
            setattr(self, "__space_group_wyckoff", space_group_wyckoff)

        if centrosymmetry is not None:
            setattr(self, "__centrosymmetry", centrosymmetry)
        if name_hm_extended is not None: 
            setattr(self, "name_hm_alt", name_hm_extended)
            setattr(self, "name_hm_alt_description", "The extended international Hermann-Mauguin space-group symbol.")
        if name_hm_full is not None: 
            setattr(self, "name_hm_full", name_hm_full)
        if name_hm_short is not None: setattr(self, "name_hm_ref", name_hm_short)
        if name_hall is not None: setattr(self, "name_hall", name_hall)
        if name_schoenflies is not None: setattr(self, "name_schoenflies", name_schoenflies)

        if point_group_hm is not None: setattr(self, "point_group_hm", point_group_hm)
        if laue_class is not None: setattr(self, "laue_class", laue_class)
        if patterson_name_hm is not None: setattr(self, "patterson_name_hm", patterson_name_hm)
        if centring_type is not None: setattr(self, "centring_type", centring_type)
        if bravais_type is not None: setattr(self, "bravais_type", bravais_type)
        if crystal_system is not None: setattr(self, "crystal_system", crystal_system)
        return flag

    def form_reduced_full_space_group_symop_space_group_wyckoff(self, it_number:int, it_coordinate_system_code:str, centrosymmetry:bool, centring_type:str):
        pcentr, reduced_space_group_symop, full_space_group_symop = None, None, None

        symop, pcentr, _multiplicity, _letter, _site_symmetry, _l_coords_xyz_2 = CONSTANTS_AND_FUNCTIONS.get_symop_pcentr_multiplicity_letter_site_symmetry_coords_xyz_2(it_number, it_coordinate_system_code)
        item_reduced = []
        for _i_symop, _symop in enumerate(symop):
            _item = SpaceGroupSymop(id=f"{_i_symop+1:}", operation_xyz=_symop, operation_description=f"reduced",sg_id=it_number)
            item_reduced.append(_item)
        reduced_space_group_symop = SpaceGroupSymopL(item=item_reduced)

        item_full = []
        for _item in item_reduced:
            _symop = _item.operation_xyz
            _item_full = SpaceGroupSymop(operation_xyz=_symop)
            item_full.append(_item_full)
        if centrosymmetry:
            item_add = []
            for _item in item_full:
                _item_add = _item.get_symop_inversed(pcentr)
                item_add.append(_item_add)
            item_full.extend(item_add)
        item_new = []
        for _item in item_full:
            _items_new = _item.get_symops_by_centring_type(centring_type)
            item_new.extend(_items_new)
        for _i, _item in enumerate(item_new):
            _item.id = f"{_i+1:}"
            _item.sg_id = f"{it_number:}"
            _item.operation_description = "full"
        full_space_group_symop =  SpaceGroupSymopL(item=item_new)

        _i_numb = 0
        item = []
        for _1, _2, _3, _4 in zip(_multiplicity, _letter, _site_symmetry, _l_coords_xyz_2):
            _i_numb += 1
            _item = SpaceGroupWyckoff(id=f"{_i_numb:}", multiplicity=_1, coord_xyz=_4[0], letter=_2, site_symmetry=_3)
            _item.set_it_coord_xyz(_4)
            _item.set_centring_type(centring_type)
            _item.form_object
            item.append(_item)
        space_group_wyckoff = SpaceGroupWyckoffL(item=item)
        return pcentr, reduced_space_group_symop, full_space_group_symop, space_group_wyckoff

    @staticmethod
    def get_it_coordinate_system_codes_by_it_number(it_number:int)->str:
        return CONSTANTS_AND_FUNCTIONS.get_it_coordinate_system_codes_by_it_number(it_number)

    @staticmethod
    def get_default_it_coordinate_system_code_by_it_number(it_number:int)->str:
        return CONSTANTS_AND_FUNCTIONS.get_default_it_coordinate_system_code_by_it_number(it_number)
        
    def _show_message(self, s_out: str):
        warnings.warn("***  Error ***\n"+s_out, UserWarning, stacklevel=2)
        
        
    def calc_hkl_equiv(self, h, k, l):
        """
Give equivalent reflections of hkl and its multiplicity
        """
        r_11 = self.reduced_space_group_symop.r_11
        r_12 = self.reduced_space_group_symop.r_12
        r_13 = self.reduced_space_group_symop.r_13
        r_21 = self.reduced_space_group_symop.r_21
        r_22 = self.reduced_space_group_symop.r_22
        r_23 = self.reduced_space_group_symop.r_23
        r_31 = self.reduced_space_group_symop.r_31
        r_32 = self.reduced_space_group_symop.r_32
        r_33 = self.reduced_space_group_symop.r_33

        h_s = r_11*h + r_21*k + r_31*l 
        k_s = r_12*h + r_22*k + r_32*l 
        l_s = r_13*h + r_23*k + r_33*l 
        
        hkl_s = numpy.vstack([h_s, k_s, l_s])
        hkl_s = numpy.hstack([hkl_s,-1*hkl_s])
        hkl_s = hkl_s.astype(int)
        hkl_s_un = numpy.unique(hkl_s, axis=1)
        multiplicity = int(round(hkl_s.shape[1]*1./hkl_s_un.shape[1]))
        h_s, k_s, l_s = hkl_s_un[0, :], hkl_s_un[1, :], hkl_s_un[2, :]
        return h_s, k_s, l_s, multiplicity


    def calc_xyz_mult(self, x, y, z):
        """
give unique x,y,z elements and calculate multiplicity for given x,y,z fract
        """
        wyckoff = self.space_group_wyckoff.get_wyckoff_for_fract(x, y, z)
        np_r = numpy.array(wyckoff.full_r, dtype=float)
        np_b = numpy.array(wyckoff.full_b, dtype=float)
        x_s = (np_r[:, 0, 0]*x + np_r[:, 0, 1]*y + np_r[:, 0, 2]*z + np_b[:, 0])%1
        y_s = (np_r[:, 1, 0]*x + np_r[:, 1, 1]*y + np_r[:, 1, 2]*z + np_b[:, 1])%1
        z_s = (np_r[:, 2, 0]*x + np_r[:, 2, 1]*y + np_r[:, 2, 2]*z + np_b[:, 2])%1

        #l_shift = self.shift
        #l_x, l_y, l_z = [], [], []
        #for _shift in l_shift:
        #    l_x.extend(numpy.mod(x_s+_shift[0],1))
        #    l_y.extend(numpy.mod(y_s+_shift[1],1))
        #    l_z.extend(numpy.mod(z_s+_shift[2],1))
        multiplicity = wyckoff.multiplicity
        #x_out = numpy.array(l_x, dtype=float)
        #y_out = numpy.array(l_y, dtype=float)
        #z_out = numpy.array(l_z, dtype=float)

        return x_s, y_s, z_s, multiplicity
    

    def calc_symop_for_xyz(self, x_in, y_in, z_in):
        x, y, z = x_in%1., y_in%1., z_in%1.

        symop = self.full_space_group_symop
        e_11, e_12, e_13, e_1 = numpy.array(symop.r_11, dtype=float), numpy.array(symop.r_12, dtype=float), numpy.array(symop.r_13, dtype=float), numpy.array(symop.b_1, dtype=float) 
        e_21, e_22, e_23, e_2 = numpy.array(symop.r_21, dtype=float), numpy.array(symop.r_22, dtype=float), numpy.array(symop.r_23, dtype=float), numpy.array(symop.b_2, dtype=float) 
        e_31, e_32, e_33, e_3 = numpy.array(symop.r_31, dtype=float), numpy.array(symop.r_32, dtype=float), numpy.array(symop.r_33, dtype=float), numpy.array(symop.b_3, dtype=float) 
        
        x_s = numpy.round(numpy.mod(e_11*x + e_12*y + e_13*z + e_1, 1), 5)
        y_s = numpy.round(numpy.mod(e_21*x + e_22*y + e_23*z + e_2, 1), 5)
        z_s = numpy.round(numpy.mod(e_31*x + e_32*y + e_33*z + e_3, 1), 5)
            
        xyz_s = numpy.vstack([x_s, y_s, z_s])
        
        xyz_s_un, unique_inverse = numpy.unique(xyz_s, return_inverse=True, axis=1)
        x_s, y_s, z_s = xyz_s_un[0, :], xyz_s_un[1, :], xyz_s_un[2, :]
        ind = (numpy.where((x-x_s)**2+(y-y_s)**2+(z-z_s)**2 < 0.00001))[0][0]

        flag = unique_inverse == ind

        item_out = [_item for _item, _flag in zip(symop.item, flag) if _flag]
        symop_out = SpaceGroupSymopL(item=item_out)
        return symop_out


    def calc_el_symm_for_xyz(self, x_in, y_in, z_in):
        """
FIXME: should be deleted
        """
        x, y, z = x_in%1., y_in%1., z_in%1.

        symop = self.full_space_group_symop
        e_11, e_12, e_13, e_1 = numpy.array(symop.r_11, dtype=float), numpy.array(symop.r_12, dtype=float), numpy.array(symop.r_13, dtype=float), numpy.array(symop.b_1, dtype=float) 
        e_21, e_22, e_23, e_2 = numpy.array(symop.r_21, dtype=float), numpy.array(symop.r_22, dtype=float), numpy.array(symop.r_23, dtype=float), numpy.array(symop.b_2, dtype=float) 
        e_31, e_32, e_33, e_3 = numpy.array(symop.r_31, dtype=float), numpy.array(symop.r_32, dtype=float), numpy.array(symop.r_33, dtype=float), numpy.array(symop.b_3, dtype=float) 
        
        x_s = numpy.round(numpy.mod(e_11*x + e_12*y + e_13*z + e_1, 1), 5)
        y_s = numpy.round(numpy.mod(e_21*x + e_22*y + e_23*z + e_2, 1), 5)
        z_s = numpy.round(numpy.mod(e_31*x + e_32*y + e_33*z + e_3, 1), 5)
            
        xyz_s = numpy.vstack([x_s, y_s, z_s])
        
        xyz_s_un, unique_inverse = numpy.unique(xyz_s, return_inverse=True, axis=1)
        x_s, y_s, z_s = xyz_s_un[0, :], xyz_s_un[1, :], xyz_s_un[2, :]
        ind = (numpy.where((x-x_s)**2+(y-y_s)**2+(z-z_s)**2 < 0.00001))[0][0]

        flag = unique_inverse == ind
        o_11, o_12, o_13 = e_11[flag], e_12[flag], e_13[flag]
        o_21, o_22, o_23 = e_21[flag], e_22[flag], e_23[flag]
        o_31, o_32, o_33 = e_31[flag], e_32[flag], e_33[flag]
        o_1, o_2, o_3 = e_1[flag], e_2[flag], e_3[flag]
        return o_11, o_12, o_13, o_21, o_22, o_23, o_31, o_32, o_33, o_1, o_2, o_3


    def calc_f_hkl_by_f_hkl_as(self, h, k, l, f_hkl_as):
        """
Calculate the structure factor giving the structure factor 
in asymmetric unit cell.
        """
        np_h, np_k, np_l = numpy.array(h, dtype=float), numpy.array(k, dtype=float), numpy.array(l, dtype=float)
        shift = numpy.array(self.shift, dtype=float)
        centr = self.centrosymmetry

        orig_x, orig_y, orig_z = shift[:, 0], shift[:, 1], shift[:, 2]

        #orig_x = [hh[0] for hh in shift]
        #orig_y = [hh[1] for hh in shift]
        #orig_z = [hh[2] for hh in shift]
        
        np_h_2d, np_orig_x_2d = numpy.meshgrid(h, orig_x, indexing = "ij")
        np_k_2d, np_orig_y_2d = numpy.meshgrid(k, orig_y, indexing = "ij")
        np_l_2d, np_orig_z_2d = numpy.meshgrid(l, orig_z, indexing = "ij")

        hh = (2*numpy.pi*1j*(np_h_2d*np_orig_x_2d+np_k_2d*np_orig_y_2d+np_l_2d*np_orig_z_2d)).astype(complex)
        np_orig_as = numpy.exp(hh)
        _hh = np_orig_as.sum(axis=1)
        if len(f_hkl_as.shape)==2:
           _hh = _hh[:, numpy.newaxis] 
        f_hkl_1 = f_hkl_as*_hh*1./len(shift)

        if (centr):
            p_centr = numpy.array(self.pcentr, dtype=float)
            hh = (2.*2.*numpy.pi*1j* (np_h*p_centr[0]+np_k*p_centr[1]+np_l*p_centr[2])).astype(complex)
            if len(f_hkl_as.shape)==2:
                hh = hh[:, numpy.newaxis]
            f_hkl = 0.5*(f_hkl_1+f_hkl_1.conjugate()*numpy.exp(hh))
            
        else:
            f_hkl = f_hkl_1
        return f_hkl


    def calc_asymmetric_cell(self, n_a, n_b, n_c):
        """
give the numbers in asymmetric cell

 :n_a: the number of points along a axis
 :n_b: the numper of points along b axis
 :n_c: the numper of points along c axis

na, n_b, nc should be divided on 24: 8 and 3

:output: - l_coord is a list of coordinates in asymmetric cell (frac_x = n_x/n_a and so on)
         - l_symm contains a list of symmetry given as (n_symm, centr, n_orig)
        """

        n_a_new = int(round(n_a/24))*24
        n_b_new = int(round(n_b/24))*24
        n_c_new = int(round(n_c/24))*24

        #print("na: {:}, n_b: {:}, n_c: {:}".format(n_a_new, n_b_new, n_c_new))

        l_el_symm = self.el_symm
        f_centr = self.centr
        p_centr = self.p_centr
        l_orig = self.orig
        l_coord = []
        
        spgr_choice = self.spgr_choice 
        spgr_name = self.spgr_name 
        spgr_number = self.spgr_number
        
        if (spgr_number==227) & (spgr_choice=="2"):
            n_a_new = int(round(n_a/8))*8
            for n_x in range(-n_a_new//8, 3*n_a_new//8+1):
                for n_y in range(-n_a_new//8, 0+1):
                    for n_z in range(-n_a_new//4, 0+1):
                        cond_1 = (n_y < min([n_a_new//4-n_x, n_x]))
                        cond_2 = (n_z >= -n_y-n_a_new//4)
                        cond_3 = (n_z <= n_y)
                        if (cond_1 & cond_2 & cond_3):
                            coord_x, coord_y = float(n_x)/float(n_a_new), float(n_y)/float(n_a_new)
                            coord_z = float(n_z)/float(n_a_new)
                            #print(" {:3} {:3} {:3}".format(n_x, n_y, n_z), " ", " {:6.3f} {:6.3f} {:6.3f}".format(coord_x, coord_y, coord_z))
                            l_coord.append((coord_x, coord_y, coord_z))
        return l_coord
    
    def calc_rotated_matrix_for_position(self, m_chi, x, y, z):

        symop = self.full_space_group_symop
        e_11, e_12, e_13, e_1 = numpy.array(symop.r_11, dtype=float), numpy.array(symop.r_12, dtype=float), numpy.array(symop.r_13, dtype=float), numpy.array(symop.b_1, dtype=float) 
        e_21, e_22, e_23, e_2 = numpy.array(symop.r_21, dtype=float), numpy.array(symop.r_22, dtype=float), numpy.array(symop.r_23, dtype=float), numpy.array(symop.b_2, dtype=float) 
        e_31, e_32, e_33, e_3 = numpy.array(symop.r_31, dtype=float), numpy.array(symop.r_32, dtype=float), numpy.array(symop.r_33, dtype=float), numpy.array(symop.b_3, dtype=float) 
        
        x_s = numpy.round(numpy.mod(e_11*x + e_12*y + e_13*z + e_1, 1), 5)
        y_s = numpy.round(numpy.mod(e_21*x + e_22*y + e_23*z + e_2, 1), 5)
        z_s = numpy.round(numpy.mod(e_31*x + e_32*y + e_33*z + e_3, 1), 5)
        #o_11, o_12, o_13, o_21, o_22, o_23, o_31, o_32, o_33, o_1, o_2, o_3 = self.calc_el_symm_for_xyz(x, y, z)
        #np_x, np_y, np_z, mult = self.calc_xyz_mult(x, y, z)
        
        l_ind, l_xyz = [], []
        _ind = 0
        for _x, _y, _z in zip(x_s, y_s, z_s):
            if (_x, _y, _z) not in l_xyz:
                l_ind.append(_ind)
                l_xyz.append((_x, _y, _z))
            _ind += 1
        l_res = []
        for _ind, _xyz in zip(l_ind, l_xyz):
            _11, _12, _13 = e_11[_ind], e_12[_ind], e_13[_ind]
            _21, _22, _23 = e_21[_ind], e_22[_ind], e_23[_ind]
            _31, _32, _33 = e_31[_ind], e_32[_ind], e_33[_ind]
            _1, _2, _3 = e_1[_ind], e_2[_ind], e_3[_ind]
            matrix_r = numpy.array([[_11, _12, _13], [_21, _22, _23], 
                               [_31, _32, _33]], dtype=float)
            matrix_rt = matrix_r.transpose()
            r_chi = numpy.matmul(matrix_r, m_chi)
            matrix_chi_rot = numpy.matmul(r_chi, matrix_rt)
            l_res.append((_xyz, matrix_chi_rot))
        return l_res
