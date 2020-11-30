import numpy
from typing import NoReturn
from cryspy.A_functions_base.function_1_matrices import\
    calc_product_matrices, calc_product_matrix_vector
from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN
from cryspy.C_item_loop_classes.cl_1_cell import Cell


class TOFProfile(ItemN):
    """
    Profile parameters-I, T.O.F

    Attributes
    ----------
        - sigma0, sigma1, sigma2, gamma0, gamma1, gamma2 (mandatory)
          size_g, size_l, strain_g, strain_l (optional)

    The variance of the Gaussian component of the peak shape in TOF neutrons is
    given by:
 
        H_G**2 = (sigma_2+size_g)*d**4 + (sigma_1+strain_g)*d**2 + sigma_0
        H_L**2 = (gamma_2+size_l)*d**4 + (sigma_1+strain_l)*d**2 + sigma_0
        
    where d is the d-spacing in angstroms.

    Units
    -----
        sigma2, gamma2 are (microsecs/Ang.**2)**2;
        sigma1, gamma1 are (microsecs/Ang.**2);
        sigma0, gamma0 are (microsecs**2);

    """
    ATTR_MANDATORY_NAMES = ("sigma2", "sigma1", "sigma0",
                            "gamma2", "gamma1", "gamma0")
    ATTR_MANDATORY_TYPES = (float, float, float, float, float, float)
    ATTR_MANDATORY_CIF = ("sigma2", "sigma1", "sigma0",
                          "gamma2", "gamma1", "gamma0")

    ATTR_OPTIONAL_NAMES = ("size_g", "size_l", "strain_g", "strain_l")
    ATTR_OPTIONAL_TYPES = (float, float, float, float)
    ATTR_OPTIONAL_CIF = ("size_g", "size_l", "strain_g", "strain_l")

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ()
    ATTR_INT_PROTECTED_NAMES = ()

    # parameters considered are refined parameters
    ATTR_REF = ("sigma2", "sigma1", "sigma0", "gamma2", "gamma1",
                "gamma0", "size_g", "size_l", "strain_g", "strain_l")
    ATTR_SIGMA = tuple([f"{_h:}_sigma" for _h in ATTR_REF])
    ATTR_CONSTR_FLAG = tuple([f"{_h:}_constraint" for _h in ATTR_REF])
    ATTR_REF_FLAG = tuple([f"{_h:}_refinement" for _h in ATTR_REF])

    # formats if cif format
    D_FORMATS = {"sigma0": "{:.5f}", "sigma1": "{:.5f}", "sigma2": "{:.5f}",
                 "gamma0": "{:.5f}", "gamma1": "{:.5f}", "gamma2": "{:.5f}",
                 "size_g": "{:.5f}", "size_l": "{:.5f}", "strain_g": "{:.5f}",
                 "strain_l": "{:.5f}"}

    # constraints on the parameters
    D_CONSTRAINTS = {}

    # default values for the parameters
    D_DEFAULT = {"sigma0": 0., "sigma1": 0., "sigma2": 0.,
                 "gamma0": 0., "gamma1": 0., "gamma2": 0.}

    for key in ATTR_SIGMA:
        D_DEFAULT[key] = 0.
    for key in (ATTR_CONSTR_FLAG + ATTR_REF_FLAG):
        D_DEFAULT[key] = False

    PREFIX = "tof_profile"

    def __init__(self, **kwargs) -> NoReturn:
        super(TOFProfile, self).__init__()

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



class TOFProfileL(LoopN):
    """
    Description of AtomSite in loop.

    """
    ITEM_CLASS = TOFProfile
    ATTR_INDEX = None
    def __init__(self, loop_name = None) -> NoReturn:
        super(TOFProfileL, self).__init__()
        self.__dict__["items"] = []
        self.__dict__["loop_name"] = loop_name
   

# s_cont = """
#   loop_
#   _tof_profile_sigma2 
#   _tof_profile_sigma1
#   _tof_profile_sigma0 
#   _tof_profile_gamma2 
#   _tof_profile_gamma1 
#   _tof_profile_gamma0
#   0.000     61.677      6.195  0.000     0.604     0.000 
# """

# obj = TOFProfileL.from_cif(s_cont)
# print(obj, end="\n\n")
