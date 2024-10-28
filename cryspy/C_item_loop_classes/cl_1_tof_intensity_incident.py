import numpy
import matplotlib.pyplot as plt
from typing import NoReturn

from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN


class TOFIntensityIncident(ItemN):
    """Correction of incident intensity for time-of-flight experiment.

    Attributes
    ----------
        - a1, a2, a3, a4, a5, a6, a7, a8 (mandatory)
        - spectrum (optional)
        

    spectrum is "Maxwell" (default), "Empirical-Exponents"
        """
    ATTR_MANDATORY_NAMES = ("a0", "a1", "a2", "a3", "a4", "a5", "a6", "a7",
                            "a8")
    ATTR_MANDATORY_TYPES = (float, float, float, float, float, float, float,
                            float, float)
    ATTR_MANDATORY_CIF = ("a0", "a1", "a2", "a3", "a4", "a5", "a6", "a7", "a8")

    ATTR_OPTIONAL_NAMES = ("spectrum", )
    ATTR_OPTIONAL_TYPES = (str, )
    ATTR_OPTIONAL_CIF = ("spectrum", )

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ()
    ATTR_INT_PROTECTED_NAMES = ()

    # parameters considered are refined parameters
    ATTR_REF = ("a0", "a1", "a2", "a3", "a4", "a5", "a6", "a7", "a8")
    ATTR_SIGMA = tuple([f"{_h:}_sigma" for _h in ATTR_REF])
    ATTR_CONSTR_FLAG = tuple([f"{_h:}_constraint" for _h in ATTR_REF])
    ATTR_REF_FLAG = tuple([f"{_h:}_refinement" for _h in ATTR_REF])
    ATTR_CONSTR_MARK = tuple([f"{_h:}_mark" for _h in ATTR_REF])

    # formats if cif format
    D_FORMATS = {}

    # constraints on the parameters
    D_CONSTRAINTS = {"spectrum": ["Maxwell", "Empirical-Exponents"]}

    # default values for the parameters
    D_DEFAULT = {"spectrum": "Maxwell", "a0": 0., "a1": 2339., "a2": 35.59e6,
                 "a3": 1., "a4": 0.2015e-6, "a5": 0.592, "a6": 0.1421e-6,
                 "a7": 0., "a8": 0.}  # from J. Appl. Cryst. (1982). 15, 581-589

    for key in ATTR_SIGMA:
        D_DEFAULT[key] = 0.
    for key in (ATTR_CONSTR_FLAG + ATTR_REF_FLAG):
        D_DEFAULT[key] = False
    for key in ATTR_CONSTR_MARK:
        D_DEFAULT[key] = ""

    PREFIX = "tof_intensity_incident"

    def __init__(self, **kwargs) -> NoReturn:
        super(TOFIntensityIncident, self).__init__()

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

    def calc_spectrum(self, time):
        """Calculate spectrum for time in microseconds
        
        Two spectra can be calculated:
            (i) Maxwell (by default);
            
            spectrum = A0 + A1*exp[-A2 time**2]/time**5 + A3*exp[-A4 time**2] + 
                A5*exp[-A6 time**3] + A7*exp[-A8 time**4]
            
            (ii) Empirical-Exponents.

            spectrum = A0 + A1*exp[-A2 time] + A3*exp[-A4 time**2] + 
                A5*exp[-A6 time**3] + A7*exp[-A8 time**4]
            
        Parameters
        ----------
        time : TYPE
            DESCRIPTION.

        Returns
        -------
        res : TYPE
            DESCRIPTION.

        """
        exp = numpy.exp
        time_sq = numpy.square(time)
        time_4 = numpy.square(time_sq)
        if self.spectrum == "Empirical-Exponents":
            res = self.a0 + self.a1 * exp(-self.a2 * time) + \
                self.a3 * exp(-self.a4 * time_sq) + \
                self.a5 * exp(-self.a6 * time*time_sq) + \
                self.a7 * exp(-self.a8 * time_4) 
        else:  # self.spectrum == "Maxwell"
            res = self.a0 + self.a1 * exp(-self.a2 * time_sq)/(time*time_4) + \
                self.a3 * exp(-self.a4 * time_sq) + \
                self.a5 * exp(-self.a6 * time*time_sq) + \
                self.a7 * exp(-self.a8 * numpy.square(time_sq)) 
        return res


    def get_coefficients(self):
        coefficients = numpy.array([
            self.a0, self.a1, self.a2, self.a3, 
            self.a4, self.a5, self.a6, self.a7, 
            self.a8], dtype=float)
        return coefficients

    def get_flags_coefficients(self):
        flags_coefficients = numpy.array([
            self.a0_refinement, self.a1_refinement, self.a2_refinement, self.a3_refinement, 
            self.a4_refinement, self.a5_refinement, self.a6_refinement, self.a7_refinement, 
            self.a8_refinement], dtype=bool)
        return flags_coefficients

    def plots(self):
        return [self.plot_spectrum()]
    
    def plot_spectrum(self):
        time_min, time_max = 4000, 18000
        time = numpy.linspace(time_min, time_max, 100)
        spectrum = self.calc_spectrum(time)
        
        fig, ax = plt.subplots()
        ax.set_title("Incident spectrum")
        ax.set_xlabel("time (mc. sec.)")
        ax.set_xlim(time_min, time_max)
        ax.set_ylabel('Incident spectrup')
        ax.plot(time, spectrum, "k-", label=self.spectrum)
        ax.legend(loc='upper right')
        fig.tight_layout()
        return (fig, ax)
    
    def get_dictionary(self):
        res = {}
        res["spectrum_incident_type"] = tof_intensity_incident.spectrum
        res["spectrum_incident_coefficients"] = tof_intensity_incident.get_coefficients()
        res["flags_spectrum_incident_coefficients"] = tof_intensity_incident.get_flags_coefficients()
        return res
        

class TOFIntensityIncidentL(LoopN):
    """Correction of incident intensity for time-of-flight experiment.

    """
    ITEM_CLASS = TOFIntensityIncident
    ATTR_INDEX = None
    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(TOFIntensityIncidentL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name
   

# s_cont = """
# _tof_intensity_incident_a1 0
# _tof_intensity_incident_a2 0
# _tof_intensity_incident_a3 0
# _tof_intensity_incident_a4 0
# _tof_intensity_incident_a5 0
# _tof_intensity_incident_a6 0
# _tof_intensity_incident_a7 0
# _tof_intensity_incident_a8 0
# """

# obj = TOFCorrectionIntensityIncident.from_cif(s_cont)
# print(obj, end="\n\n")
