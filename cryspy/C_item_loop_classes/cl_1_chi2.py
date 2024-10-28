from typing import NoReturn

from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN


class Chi2(ItemN):
    """
    Choise of the experimental data for the refinement procedure.

    Attributes
    ----------
        - sum  True
        - diff True
        - up   False
        - down False
        - asymmetry False # For flip ratios
    """

    ATTR_MANDATORY_NAMES = ()
    ATTR_MANDATORY_TYPES = ()
    ATTR_MANDATORY_CIF = ()

    ATTR_OPTIONAL_NAMES = ("sum", "diff", "up", "down", "asymmetry")
    ATTR_OPTIONAL_TYPES = (bool, bool, bool, bool, bool)
    ATTR_OPTIONAL_CIF = ("sum", "diff", "up", "down", "asymmetry")

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ()
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
    D_DEFAULT = {}# 'sum': True, 'diff': True, 'up': False, 'down': False
    for key in ATTR_SIGMA:
        D_DEFAULT[key] = 0.
    for key in (ATTR_CONSTR_FLAG + ATTR_REF_FLAG):
        D_DEFAULT[key] = False
    for key in ATTR_CONSTR_MARK:
        D_DEFAULT[key] = ""

    PREFIX = "chi2"

    def __init__(self, **kwargs) -> NoReturn:
        super(Chi2, self).__init__()

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

    def get_dictionary(self):
        res = {}
        if self.is_attribute("sum"):
            res["flag_chi_sq_sum"] = self.sum
        if self.is_attribute("diff"):
            res["flag_chi_sq_difference"] = self.diff
        if self.is_attribute("asymmetry"):
            res["flag_asymmetry"] = self.asymmetry
        return res


class Chi2L(LoopN):
    """
    Description of chi2 in loop.

    """
    ITEM_CLASS = Chi2
    ATTR_INDEX = None
    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(Chi2L, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name

# s_cont = """
# loop_
#   _chi2_sum
#   _chi2_diff
#   _chi2_up
#   _chi2_down
#   False True  False False
#   False False True True
# """

# """


# val_2 = Cell()
# val_2.length_a = 3.
# val_2.angle_alpha = 750

# """
# obj = Chi2L.from_cif(s_cont)
# print(obj, end="\n\n")
# print(obj.get_variable_names(), end="\n\n")

