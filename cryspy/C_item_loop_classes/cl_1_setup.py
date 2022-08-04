"""Setup and SetupL classes."""
from typing import NoReturn

from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN


class Setup(ItemN):
    """Experimental diffraction setup (constant wavelength).

    Attributes
    ----------
        - wavelength (mandatory) (in Angstrems)
        - field (optional) (in Tesla)
        - radiation (optional) (neutrons by default, or X-rays)
        - offset_ttheta (optional for powder 1d and 2d) (in degrees)
        - offset_phi (optional for powder 2d) (in degrees)
        - ratio_lambdaover2 (optional, for single diffraction)
        - k (0. for neutrons, 0.5 for characteristic X-ray, 0.1 for synchrotron radiation)
        - cthm (cos**2 (2 theta_M)) (for calculation of Lorentrz polarization factor)
    """

    ATTR_MANDATORY_NAMES = ()
    ATTR_MANDATORY_TYPES = ()
    ATTR_MANDATORY_CIF = ()

    ATTR_OPTIONAL_NAMES = ("wavelength", "field", "offset_ttheta", "offset_phi", "offset_gamma", "offset_nu",
                           "ratio_lambdaover2", "radiation", "k", "cthm", "temperature")
    ATTR_OPTIONAL_TYPES = (float, float, float, float, float, float, float, str, float, float, float)
    ATTR_OPTIONAL_CIF = ("wavelength", "field", "offset_2theta", "offset_phi", "offset_gamma", "offset_nu",
                         "ratio_lambda/2", "radiation", "K", "cthm", "temperature")

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ()
    ATTR_INT_PROTECTED_NAMES = ()

    # parameters considered are refined parameters
    ATTR_REF = ("wavelength", "offset_ttheta", "offset_phi", "offset_gamma", "offset_nu",
                "ratio_lambdaover2")
    ATTR_SIGMA = tuple([f"{_h:}_sigma" for _h in ATTR_REF])
    ATTR_CONSTR_FLAG = tuple([f"{_h:}_constraint" for _h in ATTR_REF])
    ATTR_REF_FLAG = tuple([f"{_h:}_refinement" for _h in ATTR_REF])
    ATTR_CONSTR_MARK = tuple([f"{_h:}_mark" for _h in ATTR_REF])

    # formats if cif format
    D_FORMATS = {'wavelength': "{:.4f}", 'field': "{:.2f}",
                 'offset_ttheta': "{:.3f}", 'offset_phi': "{:.3f}",
                 'offset_gamma': "{:.3f}", 'offset_nu': "{:.3f}",
                 "ratio_lambdaover2": "{:.3f}", "k": "{:.1f}", "cthm": "{:.5f}",
                 'temperature': "{:.2f}",}

    # constraints on the parameters
    D_CONSTRAINTS = {"radiation": ["neutrons", "X-rays"]}

    # default values for the parameters
    D_DEFAULT = {"offset_2theta": 0., "radiation": "neutrons", "k":0., "cthm": 0.91}
    for key in ATTR_SIGMA:
        D_DEFAULT[key] = 0.
    for key in (ATTR_CONSTR_FLAG + ATTR_REF_FLAG):
        D_DEFAULT[key] = False
    for key in ATTR_CONSTR_MARK:
        D_DEFAULT[key] = ""

    PREFIX = "setup"

    def __init__(self, **kwargs) -> NoReturn:
        super(Setup, self).__init__()

        # defined for any integer and float parameters
        D_MIN = {"wavelength": 0., "ratio_lambdaover2": 0.}

        # defined for ani integer and float parameters
        D_MAX = {"ratio_lambdaover2": 1.}

        self.__dict__["D_MIN"] = D_MIN
        self.__dict__["D_MAX"] = D_MAX
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)


class SetupL(LoopN):
    """Experimental diffraction setup (constant wavelength).
    """

    ITEM_CLASS = Setup
    ATTR_INDEX = None

    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(SetupL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name
