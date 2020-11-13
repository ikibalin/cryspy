"""MEMParameters, MEMParametersL classes."""
from typing import NoReturn
from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN


class MEMParameters(ItemN):
    """Description of MEMParameters.

    Attributes
    ----------
        - points_a, points_b, points_c, chi_ferro (mandatory)
        - chi_antiferro, prior_density, method (mandatory)
        - gof_desired (mandatory)

    Accessible Values
    -----------------
        method is "2channel" (default) or "tensorMEM"
        prior_density is "uniform" (default) or "core"
    """

    ATTR_MANDATORY_NAMES = ("points_a", "points_b", "points_c", "chi_ferro",
                            "chi_antiferro", "prior_density", "method",
                            "gof_desired")
    ATTR_MANDATORY_TYPES = (int, int, int, float, float, str, str, float)
    ATTR_MANDATORY_CIF = ("points_a", "points_b", "points_c", "chi_ferro",
                          "chi_antiferro", "prior_density", "method",
                          "GoF_desired")

    ATTR_OPTIONAL_NAMES = ()
    ATTR_OPTIONAL_TYPES = ()
    ATTR_OPTIONAL_CIF = ()

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ()
    ATTR_INT_PROTECTED_NAMES = ()

    # parameters considered are refined parameters
    ATTR_REF = ("chi_ferro", "chi_antiferro")
    ATTR_SIGMA = tuple([f"{_h:}_sigma" for _h in ATTR_REF])
    ATTR_CONSTR_FLAG = tuple([f"{_h:}_constraint" for _h in ATTR_REF])
    ATTR_REF_FLAG = tuple([f"{_h:}_refinement" for _h in ATTR_REF])

    # formats if cif format
    D_FORMATS = {'chi_ferro': "{:.5f}", 'chi_antiferro': "{:.5f}",
                 'gof_desired': "{:.2f}"}

    # constraints on the parameters
    D_CONSTRAINTS = {"prior_density": ["core", "uniform"],
                     "method": ["tensorMEM", "2channel"]}

    # default values for the parameters
    D_DEFAULT = {"points_a": 48, "points_b": 48, "points_c": 48,
                 "chi_ferro": 0., "chi_antiferro": 0.,
                 "prior_density": "uniform", "method": "2channel",
                 "gof_desired": 1.}
    for key in ATTR_SIGMA:
        D_DEFAULT[key] = 0.
    for key in (ATTR_CONSTR_FLAG + ATTR_REF_FLAG):
        D_DEFAULT[key] = False

    PREFIX = "mem_parameters"

    def __init__(self, **kwargs) -> NoReturn:
        super(MEMParameters, self).__init__()

        # defined for any integer and float parameters
        D_MIN = {"points_a": 0, "points_b": 0, "points_c": 0, "chi_ferro": 0.,
                 "gof_desired": 0.}

        # defined for ani integer and float parameters
        D_MAX = {"chi_antiferro": 0}

        self.__dict__["D_MIN"] = D_MIN
        self.__dict__["D_MAX"] = D_MAX
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)


class MEMParametersL(LoopN):
    """Description of chi2 in loop."""

    ITEM_CLASS = MEMParameters
    ATTR_INDEX = None

    def __init__(self, loop_name=None) -> NoReturn:
        super(MEMParametersL, self).__init__()
        self.__dict__["items"] = []
        self.__dict__["loop_name"] = loop_name


# s_cont = """
# loop_
#   _mem_parameters_points_a
#   _mem_parameters_points_b
#   _mem_parameters_points_c
#   _mem_parameters_chi_ferro
#   _mem_parameters_chi_antiferro
#   _mem_parameters_prior_density
#   48 48 48 0.0   0.0  core
#   64 32 48 1.0  -3.0() uniform
# """

# obj = MEMParametersL.from_cif(s_cont)
# print(obj, end="\n\n")
# print(obj.get_variable_names(), end="\n\n")
