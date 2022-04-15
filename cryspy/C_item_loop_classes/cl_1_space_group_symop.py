from typing import NoReturn, Tuple, List
import copy
import numpy
from fractions import Fraction

from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary
from cryspy.A_functions_base.function_1_strings import \
    transform_string_to_r_b, transform_r_b_to_string
from cryspy.A_functions_base.function_2_space_group import \
    get_shift_by_centring_type, mult_matrix_vector, mult_matrixes

from cryspy.A_functions_base.function_2_sym_elems import \
    form_symm_elems_by_b_i_r_ij

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN


class SpaceGroupSymop(ItemN):
    """Symmetry space group operation.
    
    Contains information about the symmetry operations of the
    space group.

    Attributes
    ----------
        - operation_xyz (mandatory)
        - id, operation_description, generator_xyz, sg_id (optional)
        - r, b, r_11, r_12, r_13, r_21, r_22, r_23, r_31, r_32, r_33, b_1, b_2,
          b_3 (internal attributes)

    Methods
    -------------
        - define_by_el_symm() (class method)
        - get_symop_inversed()
        - get_symops_by_centring_type()
        - get_symops_by_generator_xyz()
        - get_coords_xyz_by_coord_xyz()
        
    """
    ATTR_MANDATORY_NAMES = ("operation_xyz", )
    ATTR_MANDATORY_TYPES = (str, )
    ATTR_MANDATORY_CIF = ("operation_xyz", )

    ATTR_OPTIONAL_NAMES = ("id", "operation_description", "generator_xyz",
                           "sg_id")
    ATTR_OPTIONAL_TYPES = (str, str, str, int)
    ATTR_OPTIONAL_CIF = ("id", "operation_description", "generator_xyz",
                         "sg_id")

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ("r", "b", "r_11", "r_12", "r_13", "r_21", "r_22",
                      "r_23", "r_31", "r_32", "r_33", "b_1", "b_2", "b_3")
    ATTR_INT_PROTECTED_NAMES = ()

    # parameters considered are refined parameters
    ATTR_REF = ()
    ATTR_SIGMA = tuple([f"{_h:}_sigma" for _h in ATTR_REF])
    ATTR_CONSTR_FLAG = tuple([f"{_h:}_constraint" for _h in ATTR_REF])
    ATTR_REF_FLAG = tuple([f"{_h:}_refinement" for _h in ATTR_REF])
    ATTR_CONSTR_MARK = tuple([f"{_h:}_mark" for _h in ATTR_REF])

    # constraints on the parameters
    D_CONSTRAINTS = {}

    # default values for the parameters
    D_DEFAULT = {}
    for key in ATTR_SIGMA:
        D_DEFAULT[key] = 0.
    for key in (ATTR_CONSTR_FLAG + ATTR_REF_FLAG):
        D_DEFAULT[key] = False
    for key in ATTR_CONSTR_MARK:
        D_DEFAULT[key] = ""

    PREFIX = "space_group_symop"

    def __init__(self, **kwargs) -> NoReturn:
        super(SpaceGroupSymop, self).__init__()

        # defined for any integer and float parameters
        D_MIN = {}

        # defined for ani integer and float parameters
        D_MAX = {}

        self.__dict__["D_MIN"] = D_MIN
        self.__dict__["D_MAX"] = D_MAX
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)

    def form_object(self) -> NoReturn:
        operation_xyz = self.operation_xyz
        if operation_xyz is None:
            return False
        r, b = transform_string_to_r_b(operation_xyz, labels=("x", "y", "z"))

        self.__dict__["r"] = r
        self.__dict__["b"] = b
        self.__dict__["r_11"] = r[0, 0]
        self.__dict__["r_12"] = r[0, 1]
        self.__dict__["r_13"] = r[0, 2]
        self.__dict__["r_21"] = r[1, 0]
        self.__dict__["r_22"] = r[1, 1]
        self.__dict__["r_23"] = r[1, 2]
        self.__dict__["r_31"] = r[2, 0]
        self.__dict__["r_32"] = r[2, 1]
        self.__dict__["r_33"] = r[2, 2]
        self.__dict__["b_1"] = b[0]
        self.__dict__["b_2"] = b[1]
        self.__dict__["b_3"] = b[2]

    @classmethod
    def define_by_el_symm(cls, el_symm):
        """Define object by element of symmetry (class method).
        
        Arguments
        ---------
            - el_symm is [b_1, r_11, r_12, r_13,
                          b_2, r_21, r_22, r_23,
                          b_3, r_31, r_32, r_33]
        """
        [b_1, r_11, r_12, r_13, b_2, r_21, r_22, r_23, b_3, r_31, r_32, r_33] \
            = el_symm
        r = numpy.array([[r_11, r_12, r_13], [r_21, r_22, r_23],
                         [r_31, r_32, r_33]], dtype=int)
        b = numpy.array([b_1, b_2, b_3], dtype=float)
        symop = transform_r_b_to_string(r, b, labels=("x", "y", "z"))
        obj = cls(operation_xyz=symop)
        return obj

    def get_symop_inversed(self, pcentr=numpy.zeros(3, float)):
        r, b = self.r, self.b
        r_new = -1*r
        b_new = -1*b+2*pcentr
        _symop = transform_r_b_to_string(r_new, b_new, labels=("x", "y", "z"))
        _item = SpaceGroupSymop(operation_xyz=_symop)
        return _item

    def get_symops_by_centring_type(self, centring_type:str)->Tuple[str]:
        r, b = self.r, self.b
        r_new = r
        shift = get_shift_by_centring_type(centring_type)
        symops = []
        for _shift in shift:
            b_new =  b + numpy.array(_shift, dtype=Fraction)
            _symop = transform_r_b_to_string(r_new, b_new, labels=("x", "y", "z"))
            symop = SpaceGroupSymop(operation_xyz=_symop)
            symops.append(symop)
        return symops

    def get_symops_by_generator_xyz(self, generator_xyz:str)->Tuple[str]:
        r, b = self.r, self.b
        symop = []
        W_g, w_g = transform_string_to_r_b(generator_xyz, labels=("x", "y", "z"))
        
        flag = True
        W_i = copy.deepcopy(r)
        w_i = copy.deepcopy(b)%1
        _symop = transform_r_b_to_string(W_i, w_i, labels=("x", "y", "z"))
        symop.append(_symop)
        i_stop = 1
        while flag:
            i_stop += 1
            w_i = (mult_matrix_vector(W_g, w_i) + w_g)%1 # w3 = W2 x w1 + w2
            W_i = mult_matrixes(W_g, W_i)          # W3 = W2 x W1
            _symop = transform_r_b_to_string(W_i, w_i, labels=("x", "y", "z"))
            if ((_symop in symop) | (i_stop > 100)):
                flag = False
            else:
                symop.append(_symop)
        res = []
        gen_orig = ""
        if not(self.generator_xyz is None): gen_orig = f"{self.generator_xyz:}"
        res =[SpaceGroupSymop(id=f"{_i+1:}", operation_xyz=_symop, generator_xyz=f"{gen_orig:}_{generator_xyz:}") for _i, _symop in enumerate(symop)]
        return res

    def get_coords_xyz_by_coord_xyz(self, coord_xyz: str) -> Tuple:
        _item = SpaceGroupSymop(operation_xyz=coord_xyz)
        generator_xyz = self.operation_xyz
        res = [_.operation_xyz for _ in _item.get_symops_by_generator_xyz(generator_xyz)]
        return tuple(res)



class SpaceGroupSymopL(LoopN):
    """Symmetry space group operations.

    Contains information about the symmetry operations of the space group.
    """
    ITEM_CLASS = SpaceGroupSymop
    ATTR_INDEX = "id"
    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(SpaceGroupSymopL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name

    @classmethod
    def create_by_generators_xyz(cls, generators_xyz: List):
        """
        Create list of symmetry operators by generators.

        Parameters
        ----------
        cls : TYPE
            DESCRIPTION.
        generators_xyz : List
            DESCRIPTION.

        Returns
        -------
        _obj : TYPE
            DESCRIPTION.

        """
        _item = SpaceGroupSymop(operation_xyz="x,y,z")
        items_1 = [_item]
        for generator_xyz in generators_xyz:
            items_2 = []
            for _item in items_1:
                items_2.extend(_item.get_symops_by_generator_xyz(
                    generator_xyz))
            items_1 = items_2
        symop = []
        item_uniq = []
        _i = 0
        for _item in items_1:
            if not(_item.operation_xyz in symop):
                _i += 1
                symop.append(_item.operation_xyz)
                _item.id = f"{_i:}"
                item_uniq.append(_item)
        _obj = cls()
        _obj.items = item_uniq
        return _obj

    def get_coords_xyz_by_coord_xyz(self, coord_xyz: str) -> Tuple:
        """
        Calculate new fraction for given fraction.

        Parameters
        ----------
        coord_xyz : str
            DESCRIPTION.

        Returns
        -------
        Tuple
            DESCRIPTION.

        """
        items = self.items
        res = []
        for _item in items:
            res.extend(_item.get_coords_xyz_by_coord_xyz(coord_xyz))
        return tuple(frozenset(res))

    def get_symop_for_x1_x2(self, fract_xyz_1: numpy.ndarray,
                            fract_xyz_2: numpy.ndarray):
        """
        Calculate symmetry operators transforming one fraction to another.

        Parameters
        ----------
        fract_xyz_1 : numpy.ndarray
            DESCRIPTION.
        fract_xyz_2 : numpy.ndarray
            DESCRIPTION.

        Returns
        -------
        res : TYPE
            DESCRIPTION.

        """
        np_fract_xyz_1 = numpy.array(fract_xyz_1, dtype=float)
        np_fract_xyz_2 = numpy.array(fract_xyz_2, dtype=float)
        res = SpaceGroupSymopL()
        l_item = []
        for item_s_g_s in self.items:
            m_r = numpy.array(item_s_g_s.r, dtype=float)
            v_b = numpy.array(item_s_g_s.b, dtype=float)
            if numpy.allclose(numpy.mod(numpy.dot(m_r, np_fract_xyz_1)+v_b, 1),
                              numpy.mod(np_fract_xyz_2, 1)):
                l_item.append(item_s_g_s)
        res.items = l_item
        return res

    def get_symm_elems(self):
        r_ij = (self.r_11, self.r_12, self.r_13, self.r_21, self.r_22, self.r_23, self.r_31, self.r_32, self.r_33)
        b_i = (self.b_1, self.b_2, self.b_3)
        symm_elems = form_symm_elems_by_b_i_r_ij(b_i, r_ij)
        return symm_elems
    
# s_cont = """
# loop_
# _space_group_symop.id
# _space_group_symop.operation_xyz
# _space_group_symop.operation_description
#   1    x,y,z              'identity mapping'
#   2    -x,-y,-z           'inversion'
#   3    -x,1/2+y,1/2-z '2-fold screw rotation with axis in (0,y,1/4)'
#   4    x,1/2-y,1/2+z  'c glide reflection through the plane (x,1/4,y)'
# """

# obj = SpaceGroupSymopL.from_cif(s_cont)
# print(obj, end="\n\n")
# print(obj["3"], end="\n\n")
# print(obj["3"].r_11, end="\n\n")
# print(obj.numpy_operation_xyz, end="\n\n")
# print(obj.numpy_r_11.astype(int), end="\n\n")
