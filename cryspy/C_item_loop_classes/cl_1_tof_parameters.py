import numpy
import scipy
import scipy.special

from typing import NoReturn

from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN


# FIXME: estimate d by time for epithermal neutrons
class TOFParameters(ItemN):
    """Parameters of the reflexion positions in time-of-flight experiments.

    Attributes
    ----------
        - zero, dtt1, ttheta_bank (mandatory)
        - neutrons, dtt2, zerot, dtt1t, dtt2t, width, x_cross, field,
          exinction (optional)

    for thermal neutrons
    
        time = zero + dtt1 * d + dtt2 * d**2

    or for epithermal neutrons
    
        time_e = zero + dtt1 * d 
        
        time_t = zerot + dtt1t * d - dtt2t / d 
        
        n_cross = 0.5*erfc(Width * (x_cross - 1/d))
        
        time = n_cross * time_e + (1-n_cross) time_t

    time is given in microseconds.
    """
    ATTR_MANDATORY_NAMES = ("zero", "dtt1", "ttheta_bank")
    ATTR_MANDATORY_TYPES = (float, float, float)
    ATTR_MANDATORY_CIF = ("Zero", "Dtt1", "2theta_bank")

    ATTR_OPTIONAL_NAMES = ('neutrons', "dtt2", "zerot", "dtt1t", "dtt2t",
                           "width", "x_cross", "field", "extinction")
    ATTR_OPTIONAL_TYPES = (str, float, float, float, float, float, float,
                           float, float)
    ATTR_OPTIONAL_CIF = ('neutrons', "dtt2", "zerot", "dtt1t", "dtt2t",
                         "width", "x_cross", "field", "extinction")

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ()
    ATTR_INT_PROTECTED_NAMES = ()

    # parameters considered are refined parameters
    ATTR_REF = ("zero", "dtt1", "dtt2", "zerot", "dtt1t", "dtt2t")
    ATTR_SIGMA = tuple([f"{_h:}_sigma" for _h in ATTR_REF])
    ATTR_CONSTR_FLAG = tuple([f"{_h:}_constraint" for _h in ATTR_REF])
    ATTR_REF_FLAG = tuple([f"{_h:}_refinement" for _h in ATTR_REF])
    ATTR_CONSTR_MARK = tuple([f"{_h:}_mark" for _h in ATTR_REF])

    # formats if cif format
    D_FORMATS = {
        "zero": "{:.5f}", "dtt1": "{:.5f}", "dtt2": "{:.5f}",
        "zerot": "{:.5f}", "dtt1t": "{:.5f}", "dtt2t": "{:.5f}",
        "ttheta_bank": "{:.2f}", "width": "{:.2f}", "x_cross": "{:.2f}",
        "field": "{:.3f}", "extinction": "{:.5f}"}

    # constraints on the parameters
    D_CONSTRAINTS = {"neutrons": ["thermal", "epithermal"]}

    # default values for the parameters
    D_DEFAULT = {"zero": 0., "neutrons": "thermal"}

    for key in ATTR_SIGMA:
        D_DEFAULT[key] = 0.
    for key in (ATTR_CONSTR_FLAG + ATTR_REF_FLAG):
        D_DEFAULT[key] = False
    for key in ATTR_CONSTR_MARK:
        D_DEFAULT[key] = ""

    PREFIX = "tof_parameters"

    def __init__(self, **kwargs) -> NoReturn:
        super(TOFParameters, self).__init__()

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

    def calc_time_by_d(self, d):
        """Calculate time by given d
        
        time = zero + dtt1 * d + dtt2 * d**2
        (if dtt2 is defined)

        time = zero + dtt1 * d
        (if dtt2 is not defined)

        Parameters
        ----------
        d : TYPE
            DESCRIPTION.

        Returns
        -------
        time : TYPE
            DESCRIPTION.

        """
        if self.neutrons == "epithermal":
            time_e = self.zero + self.dtt1 * d
            time_t = self.zerot + self.dtt1t * d - self.dtt2t / d 
            n_cross = 0.5*scipy.special.erfc(self.width * (self.x_cross - 1/d))
            time = n_cross * time_e + (1-n_cross) * time_t
        else:  # self.neutrons == "thermal"
            time = self.zero + self.dtt1 * d + self.dtt2 * d**2
        return time

    def calc_d_min_max(self, time):
        """Calculate d_min, d_max 
        

        Parameters
        ----------
        time : TYPE
            DESCRIPTION.

        Returns
        -------
        d_min : TYPE
            DESCRIPTION.
        d_max : TYPE
            DESCRIPTION.

        """
        time_min = numpy.min(time)
        time_max = numpy.max(time)
        if self.neutrons == "epithermal":
            raise AttributeError("function calc_d_min_max is not realized for \
epithermal neutrons")
            d_min = (time_min-self.zero)/self.dtt1
            d_max = (time_max-self.zero)/self.dtt1
        else:  # self.neutrons == "thermal"
            det_sq_min = self.dtt1**2 - 4.*self.dtt2*(self.zero - time_min)
            det_sq_max = self.dtt1**2 - 4.*self.dtt2*(self.zero - time_max)
            d_max = (-self.dtt1+det_sq_max**0.5)/(2.*self.dtt2)
            d_min = (-self.dtt1+det_sq_min**0.5)/(2.*self.dtt2)
        return d_min, d_max
    
    def calc_d_by_time(self, time):
        """Calculate d by given time 
        
        Relation between d and time is
        
        time = zero + dtt1 * d + dtt2 * d**2
        (if dtt2 is defined)

        time = zero + dtt1 * d
        (if dtt2 is not defined)

        Parameters
        ----------
        time : TYPE
            DESCRIPTION.

        Returns
        -------
        d : TYPE
            DESCRIPTION.

        """
        if self.neutrons == "epithermal":
            raise AttributeError("function calc_d_min_max is not realized for \
epithermal neutrons")
            d = (time - self.zero)/self.dtt1
        else:  # self.neutrons == "thermal"
            det = numpy.sqrt(self.dtt1**2 - 4.*(self.zero-time)*self.dtt2)
            if self.dtt2 < 0.:
                d = (det-self.dtt1)/(2.*self.dtt2)
            elif self.dtt2 > 0.:
                d = (-det-self.dtt1)/(2.*self.dtt2)
            else:
                d = (time - self.zero)/self.dtt1
        return d

    def calc_time_by_sthovl(self, sthovl):
        """Calculate time by given sin(theta)/lambda
        
        Returns
        -------
        time : TYPE
            DESCRIPTION.

        """
        d = 0.5/sthovl
        time = self.calc_time_by_d(d)
        return time

    def get_dictionary(self):
        res = {}
        res["ttheta_bank"] = getattr(self, "ttheta_bank") * numpy.pi/180
        res["neutron_type"] = getattr(self, "neutrons")
        res["zero"] = numpy.array([getattr(self, "zero"), ], dtype=float)
        res["flags_zero"] = numpy.array([getattr(self, "zero_refinement"), ], dtype=bool)
        res["dtt1"] = numpy.array([getattr(self, "dtt1"), ], dtype=float)
        res["flags_dtt1"] = numpy.array([getattr(self, "dtt1_refinement"), ], dtype=bool)
        if self.is_attribute("dtt2"):
            res["dtt2"] = numpy.array([getattr(self, "dtt2"), ], dtype=float)
            res["flags_dtt2"] = numpy.array([getattr(self, "dtt2_refinement"), ], dtype=bool)
        if self.is_attribute("zerot"):
            res["zerot"] = numpy.array([getattr(self, "zerot"), ], dtype=float)
            res["flags_zerot"] = numpy.array([getattr(self, "zerot_refinement"), ], dtype=bool)
        if self.is_attribute("dtt1t"):
            res["dtt1t"] = numpy.array([getattr(self, "dtt1t"), ], dtype=float)
            res["flags_dtt1t"] = numpy.array([getattr(self, "dtt1t_refinement"), ], dtype=bool)
        if self.is_attribute("dtt2t"):
            res["dtt2t"] = numpy.array([getattr(self, "dtt2t"), ], dtype=float)
            res["flags_dtt2t"] = numpy.array([getattr(self, "dtt2t_refinement"), ], dtype=bool)
        return res

class TOFParametersL(LoopN):
    """Parameters of the reflexion positions in time-of-flight experiments.

    """
    ITEM_CLASS = TOFParameters
    ATTR_INDEX = None
    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(TOFParametersL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name
   

# s_cont = """
# _tof_parameters_zero 2.921
# _tof_parameters_dtt1 6167.247
# _tof_parameters_dtt2 -2.280
# _tof_parameters_2theta_bank 145.
# """

# obj = TOFParameters.from_cif(s_cont)
# print(obj, end="\n\n")
