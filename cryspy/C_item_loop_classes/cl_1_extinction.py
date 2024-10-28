from typing import NoReturn
import numpy

from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary
from cryspy.A_functions_base.function_3_extinction import \
    calc_extinction_2

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN


class Extinction(ItemN):
    """
    Extinction class.

    Attributes
    ----------
        - mosaicity 100.0 # in minutes
        - radius    50.0 # in micrometers
        - model     "gauss" or "lorentz"



    """

    ATTR_MANDATORY_NAMES = ("model", "mosaicity", "radius")
    ATTR_MANDATORY_TYPES = (str, float, float)
    ATTR_MANDATORY_CIF = ("model", "mosaicity", "radius")

    ATTR_OPTIONAL_NAMES = ()
    ATTR_OPTIONAL_TYPES = ()
    ATTR_OPTIONAL_CIF = ()

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ()
    ATTR_INT_PROTECTED_NAMES = ()

    # parameters considered are refined parameters
    ATTR_REF = ("mosaicity", "radius")
    ATTR_SIGMA = tuple([f"{_h:}_sigma" for _h in ATTR_REF])
    ATTR_CONSTR_FLAG = tuple([f"{_h:}_constraint" for _h in ATTR_REF])
    ATTR_REF_FLAG = tuple([f"{_h:}_refinement" for _h in ATTR_REF])
    ATTR_CONSTR_MARK = tuple([f"{_h:}_mark" for _h in ATTR_REF])

    # constraints on the parameters
    D_CONSTRAINTS = {"model": ["gauss", "lorentz"]}

    # default values for the parameters
    D_DEFAULT = {}
    for key in ATTR_SIGMA:
        D_DEFAULT[key] = 0.
    for key in (ATTR_CONSTR_FLAG + ATTR_REF_FLAG):
        D_DEFAULT[key] = False
    for key in ATTR_CONSTR_MARK:
        D_DEFAULT[key] = ""

    PREFIX = "extinction"

    def __init__(self, **kwargs) -> NoReturn:
        super(Extinction, self).__init__()

        # defined for any integer and float parameters
        D_MIN = {"mosaicity": 0., "radius": 0.}

        # defined for ani integer and float parameters
        D_MAX = {}

        self.__dict__["D_MIN"] = D_MIN
        self.__dict__["D_MAX"] = D_MAX
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)

    def calc_extinction(self, cell, h, k, l, f_sq, wavelength:float,
                        flag_derivatives:bool=False):
        """
        f_sq in 10-12cm
        extinction for spherical model

        Parameters
        ----------
        cell : TYPE
            DESCRIPTION.
        h : TYPE
            DESCRIPTION.
        k : TYPE
            DESCRIPTION.
        l : TYPE
            DESCRIPTION.
        f_sq : TYPE
            DESCRIPTION.
        wavelength : TYPE
            DESCRIPTION.
        flag_derivatives : TYPE, optional
            DESCRIPTION. The default is False.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        radius, mosaicity, model = self.radius, self.mosaicity, self.model
        volume_unit_cell = cell.volume
        sthovl = cell.calc_sthovl(h, k, l)
        y_ext, dder = calc_extinction_2(radius, mosaicity, model, f_sq,
                                        volume_unit_cell, sthovl, wavelength)
        if flag_derivatives:
            return y_ext, dder
        else:
            return y_ext

    def get_dictionary(self):
        res = {}
        model_extinction = self.model
        radius = self.radius
        mosaicity = self.mosaicity
        res["extinction_model"] = model_extinction
        res["extinction_radius"] = numpy.array([radius], dtype=float)
        res["extinction_mosaicity"] = numpy.array([mosaicity], dtype=float)
        res["flags_extinction_radius"] = numpy.array([self.radius_refinement], dtype=bool)
        res["flags_extinction_mosaicity"] = numpy.array([self.mosaicity_refinement], dtype=bool)
        return res

class ExtinctionL(LoopN):
    """
    Description of extinction in loop.

    """
    ITEM_CLASS = Extinction
    ATTR_INDEX = None
    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(ExtinctionL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name

# s_cont = """
#   loop_
#  _extinction_mosaicity
#  _extinction_radius
#  _extinction_model
#  100() 50 gauss
#  3 7() lorentz
# """

# obj = ExtinctionL.from_cif(s_cont)
# print(obj, end="\n\n")
