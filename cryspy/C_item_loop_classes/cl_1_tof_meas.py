from typing import NoReturn
import numpy
import matplotlib.pyplot as plt

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN


class TOFMeas(ItemN):
    """Measured data for time-of-flight experiments.

    Attributes
    ----------
        - time (mandatory)
        - intensity_up, intensity_up_sigma, intensity_down,
          intensity_down_sigma, intensity, intensity_sigma (optional)

    time is given in microseconds.
    """
    ATTR_MANDATORY_NAMES = ("time", )
    ATTR_MANDATORY_TYPES = (float, )
    ATTR_MANDATORY_CIF = ("time", )

    ATTR_OPTIONAL_NAMES = ("intensity_up", "intensity_up_sigma",
                           "intensity_down", "intensity_down_sigma",
                           "intensity", "intensity_sigma")
    ATTR_OPTIONAL_TYPES = (float, float, float, float, float, float)
    ATTR_OPTIONAL_CIF = ("intensity_up", "intensity_up_sigma",
                         "intensity_down", "intensity_down_sigma",
                         "intensity", "intensity_sigma")

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

    # formats if cif format
    D_FORMATS = {'intensity_up': "{:.2f}", 'intensity_up_sigma': "{:.2f}",
                 'intensity_down': "{:.2f}", 'intensity_down_sigma': "{:.2f}",
                 'intensity': "{:.2f}", 'intensity_sigma': "{:.2f}",
                 'time': "{:.3f}"}

    # constraints on the parameters
    D_CONSTRAINTS = {}

    # default values for the parameters
    D_DEFAULT = {}
    for key in ATTR_SIGMA:
        D_DEFAULT[key] = 0.
    for key in (ATTR_CONSTR_FLAG + ATTR_REF_FLAG):
        D_DEFAULT[key] = False

    PREFIX = "tof_meas"

    def __init__(self, **kwargs) -> NoReturn:
        super(TOFMeas, self).__init__()

        # defined for any integer and float parameters
        D_MIN = {"time": 0.}

        # defined for ani integer and float parameters
        D_MAX = {}

        self.__dict__["D_MIN"] = D_MIN
        self.__dict__["D_MAX"] = D_MAX
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)

    def apply_constraints(self):
        """Apply constraints."""
        keys = self.__dict__.keys()
        if (("intensity_up" in keys) & ("intensity_up_sigma" not in keys)):
            self.intensity_up_sigma = (self.intensity_up)**0.5
        if (("intensity_down" in keys) & ("intensity_down_sigma" not in keys)):
            self.intensity_down_sigma = (self.intensity_down)**0.5
        if (("intensity" in keys) & ("intensity_sigma" not in keys)):
            self.intensity_sigma = (self.intensity)**0.5


class TOFMeasL(LoopN):
    """Measured data for time-of-flight experiments.

    This section contains the measured diffractogram and information
    about the conditions used for the measurement of the diffraction 
    data set, prior to processing and application of correction
    terms. While additional information may be added to the CIF
    as data are processed and transported between laboratories
    (possibly with the addition of a new _pd_block_id entry), the
    information in this section of the CIF will rarely be changed
    once data collection is complete.

    Where possible, measurements in this section should have no
    post-collection processing applied (normalizations, corrections,
    smoothing, zero-offset corrections etc.). Such corrected
    measurements should be recorded in the _pd_proc_ section.

    Data sets that are measured as counts, where a standard
    uncertainty can be considered equivalent to the standard
    deviation and where the standard deviation can be estimated
    as the square root of the number of counts recorded, should
    use the _pd_meas_counts_ fields. All other intensity values
    should be recorded using _pd_meas_intensity_.
    """

    ITEM_CLASS = TOFMeas
    ATTR_INDEX = "time"

    def __init__(self, loop_name: str = None) -> NoReturn:
        super(TOFMeasL, self).__init__()
        self.__dict__["items"] = []
        self.__dict__["loop_name"] = loop_name

    def is_polarized(self):
        """Give True if polaraized data are defined."""
        flag = False
        if len(self.items) != 0:
            items_0 = self.items[0]
            flag_up = items_0.is_attribute("intensity_up")
            flag_down = items_0.is_attribute("intensity_down")
            flag = flag_up & flag_down
        return flag

    def apply_constraints(self):
        """Redefined applyied constraints."""
        for item in self.items:
            item.apply_constraints()

    def plots(self):
        return [self.plot_up_down(), self.plot_sum(), self.plot_diff()]

    def plot_up_down(self):
        """Plot experimental intensity up and down vs. 2 theta (degrees)
        """
        fig, ax = plt.subplots()
        ax.set_title("I_up and I_down")
        ax.set_xlabel("Time (microseconds)")
        ax.set_ylabel('Intensity')

        if (self.is_attribute("time") & self.is_attribute("intensity_up") & 
            self.is_attribute("intensity_up_sigma") &
            self.is_attribute("intensity_down") & 
            self.is_attribute("intensity_down_sigma")):
            np_time = numpy.array(self.time, dtype=float)
            np_up = numpy.array(self.intensity_up, dtype=float)
            np_sup = numpy.array(self.intensity_up_sigma, dtype=float)
            np_down = numpy.array(self.intensity_down, dtype=float)
            np_sdown = numpy.array(self.intensity_down_sigma, dtype=float)

            ax.plot(np_time, np_up, "r-", alpha=0.2)
            ax.errorbar(np_time, np_up, yerr=np_sup, fmt="ro", alpha=0.2,
                        label="up")
            ax.plot(np_time, np_down, "b-", alpha=0.2)
            ax.errorbar(np_time, np_down, yerr=np_sdown, fmt="bs", alpha=0.2,
                        label="down")
        else:
            return
        ax.legend(loc='upper right')
        fig.tight_layout()
        return (fig, ax)
    
    def plot_sum(self):
        """Plot experimental unpolarized intensity vs. time
        """
        fig, ax = plt.subplots()
        ax.set_title("Unpolarized intensity: I_up + I_down")
        ax.set_xlabel("Time (microseconds)")
        ax.set_ylabel('Intensity')

        if (self.is_attribute("time") & self.is_attribute("intensity_up") & 
            self.is_attribute("intensity_up_sigma") &
            self.is_attribute("intensity_down") & 
            self.is_attribute("intensity_down_sigma")):
            np_time = numpy.array(self.time, dtype=float)
            np_up = numpy.array(self.intensity_up, dtype=float)
            np_sup = numpy.array(self.intensity_up_sigma, dtype=float)
            np_down = numpy.array(self.intensity_down, dtype=float)
            np_sdown = numpy.array(self.intensity_down_sigma, dtype=float)
            np_sum = np_up + np_down
            np_ssum = numpy.sqrt(numpy.square(np_sup)+numpy.square(np_sdown))
            ax.plot(np_time, np_sum, "k-", alpha=0.2)
            ax.errorbar(np_time, np_sum, yerr=np_ssum, fmt="ko", alpha=0.2,
                        label="experiment")
        elif (self.is_attribute("time") & self.is_attribute("intensity") & 
              self.is_attribute("intensity_sigma")):
            np_time = numpy.array(self.time, dtype=float)
            np_sum = numpy.array(self.intensity, dtype=float)
            np_ssum = numpy.array(self.intensity_sigma, dtype=float)
            ax.plot(np_time, np_sum, "k-", alpha=0.2)
            ax.errorbar(np_time, np_sum, yerr=np_ssum, fmt="ko", alpha=0.2,
                        label="experiment")
        ax.legend(loc='upper right')
        fig.tight_layout()
        return (fig, ax)

    def plot_diff(self):
        """Plot experimental polarized intensity vs. time
        """
        if not(self.is_attribute("time") & self.is_attribute("intensity_up") & 
               self.is_attribute("intensity_up_sigma") &
               self.is_attribute("intensity_down") & 
               self.is_attribute("intensity_down_sigma")):
            return
        fig, ax = plt.subplots()
        ax.set_title("Polarized intensity: I_up - I_down")
        ax.set_xlabel("Time (microseconds)")
        ax.set_ylabel('Intensity')
            
        np_time = numpy.array(self.time, dtype=float)
        np_up = numpy.array(self.intensity_up, dtype=float)
        np_sup = numpy.array(self.intensity_up_sigma, dtype=float)
        np_down = numpy.array(self.intensity_down, dtype=float)
        np_sdown = numpy.array(self.intensity_down_sigma, dtype=float)
        np_diff = np_up - np_down
        np_sdiff = numpy.sqrt(numpy.square(np_sup)+numpy.square(np_sdown))

        ax.plot([np_time.min(), np_time.max()], [0., 0.], "k:")
        ax.plot(np_time, np_diff, "k-", alpha=0.2)
        ax.errorbar(np_time, np_diff, yerr=np_sdiff, fmt="ko", alpha=0.2,
                    label="experiment")
        ax.legend(loc='upper right')
        fig.tight_layout()
        return (fig, ax)

# s_cont = """
#   loop_
#   _tof_meas_time
#   _tof_meas_intensity_up
#   _tof_meas_intensity_up_sigma
#   _tof_meas_intensity_down
#   _tof_meas_intensity_down_sigma
# 3001.589    40409.0    462.0    40409.0    462.0    
# 3003.09     40171.0    460.0    40409.0    462.0    
# 3004.591    39733.0    458.0   40409.0    462.0    
# """

# obj = TOFMeasL.from_cif(s_cont)
# print(obj, end="\n\n")
