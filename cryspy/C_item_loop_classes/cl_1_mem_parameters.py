"""MEMParameters, MEMParametersL classes."""
from typing import NoReturn

import numpy

from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary

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

    ATTR_MANDATORY_NAMES = ("points_a", "points_b", "points_c", )
    ATTR_MANDATORY_TYPES = (int, int, int, )
    ATTR_MANDATORY_CIF = ("points_a", "points_b", "points_c", )

    ATTR_OPTIONAL_NAMES = ("channel_plus_minus",
                          "channel_chi", "magnetization_plus", "magnetization_minus",
                          "prior_density", "use_asymmetry", "only_magnetic_basins",
                          "gof_desired", "spin_density_to_den_file", "magnetization_density_to_den_file")
    ATTR_OPTIONAL_TYPES = (bool, bool, float, float, str, bool, bool, float, str, str)
    ATTR_OPTIONAL_CIF = ("channel_plus_minus",
                          "channel_chi", "magnetization_plus", "magnetization_minus",
                          "prior_density", "use_asymmetry", "only_magnetic_basins",
                          "GoF_desired", "spin_density_to_den_file", "magnetization_density_to_den_file")

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ()
    ATTR_INT_PROTECTED_NAMES = ()

    # parameters considered are refined parameters
    ATTR_REF = ("magnetization_plus", "magnetization_minus")
    ATTR_SIGMA = tuple([f"{_h:}_sigma" for _h in ATTR_REF])
    ATTR_CONSTR_FLAG = tuple([f"{_h:}_constraint" for _h in ATTR_REF])
    ATTR_REF_FLAG = tuple([f"{_h:}_refinement" for _h in ATTR_REF])
    ATTR_CONSTR_MARK = tuple([f"{_h:}_mark" for _h in ATTR_REF])

    # formats if cif format
    D_FORMATS = {'magnetization_plus': "{:.5f}", 'magnetization_minus': "{:.5f}",
                 'gof_desired': "{:.2f}"}

    # constraints on the parameters
    D_CONSTRAINTS = {"prior_density": ["core", "uniform"],
                     }

    # default values for the parameters
    D_DEFAULT = {"points_a": 48, "points_b": 48, "points_c": 48,
                 "channel_plus_minus": True,
                 "magnetization_plus": 4., "magnetization_minus": -1.,
                 "channel_chi": False,
                 "only_magnetic_basins": False,
                 "prior_density": "uniform", "use_asymmetry": True,
                 "gof_desired": 1.,
                 "spin_density_to_den_file": "mem_spin_density.den",
                 "magnetization_density_to_den_file": "mem_magnetization_density.den",}
    for key in ATTR_SIGMA:
        D_DEFAULT[key] = 0.
    for key in (ATTR_CONSTR_FLAG + ATTR_REF_FLAG):
        D_DEFAULT[key] = False
    for key in ATTR_CONSTR_MARK:
        D_DEFAULT[key] = ""

    PREFIX = "mem_parameters"

    def __init__(self, **kwargs) -> NoReturn:
        super(MEMParameters, self).__init__()

        # defined for any integer and float parameters
        D_MIN = {"points_a": 1, "points_b": 1, "points_c": 1, "magnetization_plus": 0.,
                 "gof_desired": 0.}

        # defined for ani integer and float parameters
        D_MAX = {"magnetization_minus": 0}

        self.__dict__["D_MIN"] = D_MIN
        self.__dict__["D_MAX"] = D_MAX
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)

    def get_dictionary(self):
        numb_a = self.points_a
        numb_b = self.points_b
        numb_c = self.points_c
        points_abc = numpy.array([numb_a, numb_b, numb_c], dtype=int)
        
        channel_plus_minus = self.channel_plus_minus
        channel_chi = self.channel_chi
        magnetization_plus = self.magnetization_plus
        magnetization_minus = self.magnetization_minus
        flag_uniform_prior_density = self.prior_density == "uniform"
        flag_asymmetry = self.use_asymmetry
        flag_only_magnetic_basins = self.only_magnetic_basins
        gof_desired = self.gof_desired

        file_spin_density = self.spin_density_to_den_file
        file_magnetization_density = self.magnetization_density_to_den_file
        dict_mem = {
            'points_abc': points_abc,
            "channel_plus_minus": channel_plus_minus,
            "channel_chi": channel_chi,
            "magnetization_plus": magnetization_plus,
            "magnetization_minus": magnetization_minus,
            "flag_uniform_prior_density": flag_uniform_prior_density,
            "flag_asymmetry": flag_asymmetry,
            "flag_only_magnetic_basins": flag_only_magnetic_basins,
            "gof_desired": gof_desired,
            "file_spin_density": file_spin_density,
            "file_magnetization_density": file_magnetization_density,
            }
        return dict_mem

class MEMParametersL(LoopN):
    """Description of chi2 in loop."""

    ITEM_CLASS = MEMParameters
    ATTR_INDEX = None

    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(MEMParametersL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name

