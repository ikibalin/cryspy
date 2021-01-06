from typing import NoReturn
import numpy
import matplotlib.pyplot as plt

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN


class PdMeas(ItemN):
    """Measured data point.

    Attributes
    ----------
        - ttheta (mandatory)
        - intensity_up, intensity_up_sigma, intensity_down,
          intensity_down_sigma, intensity, intensity_sigma (optional)
    """
    ATTR_MANDATORY_NAMES = ("ttheta", )
    ATTR_MANDATORY_TYPES = (float, )
    ATTR_MANDATORY_CIF = ("2theta", )

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

    # constraints on the parameters
    D_CONSTRAINTS = {}

    # default values for the parameters
    D_DEFAULT = {}
    for key in ATTR_SIGMA:
        D_DEFAULT[key] = 0.
    for key in (ATTR_CONSTR_FLAG + ATTR_REF_FLAG):
        D_DEFAULT[key] = False

    PREFIX = "pd_meas"

    def __init__(self, **kwargs) -> NoReturn:
        super(PdMeas, self).__init__()

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

    def apply_constraints(self):
        """Apply constraints."""
        keys = self.__dict__.keys()
        if (("intensity_up" in keys) & ("intensity_up_sigma" not in keys)):
            self.intensity_up_sigma = (self.intensity_up)**0.5
        if (("intensity_down" in keys) & ("intensity_down_sigma" not in keys)):
            self.intensity_down_sigma = (self.intensity_down)**0.5
        if (("intensity" in keys) & ("intensity_sigma" not in keys)):
            self.intensity_sigma = (self.intensity)**0.5


class PdMeasL(LoopN):
    """Measured data points.

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

    ITEM_CLASS = PdMeas
    ATTR_INDEX = "ttheta"

    def __init__(self, loop_name: str = None) -> NoReturn:
        super(PdMeasL, self).__init__()
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
        return [self.plot_up_down(), self.plot_sum_diff()]

    def plot_up_down(self):
        """Plot experimental intensity up and down vs. 2 theta (degrees)
        """
        fig, ax = plt.subplots()
        ax.set_title("I_up and I_down")
        ax.set_xlabel("2 theta (degrees)")
        ax.set_ylabel('Intensity')

        if (self.is_attribute("ttheta") & self.is_attribute("intensity_up") & 
            self.is_attribute("intensity_up_sigma") &
            self.is_attribute("intensity_down") & 
            self.is_attribute("intensity_down_sigma")):
            np_tth = numpy.array(self.ttheta, dtype=float)
            np_up = numpy.array(self.intensity_up, dtype=float)
            np_sup = numpy.array(self.intensity_up_sigma, dtype=float)
            np_down = numpy.array(self.intensity_down, dtype=float)
            np_sdown = numpy.array(self.intensity_down_sigma, dtype=float)

            ax.plot(np_tth, np_up, "r-", alpha=0.2)
            ax.errorbar(np_tth, np_up, yerr=np_sup, fmt="ro", alpha=0.2,
                        label="up")
            ax.plot(np_tth, np_down, "b-", alpha=0.2)
            ax.errorbar(np_tth, np_down, yerr=np_sdown, fmt="bs", alpha=0.2,
                        label="down")
        else:
            return
        ax.legend(loc='upper right')
        fig.tight_layout()
        return (fig, ax)

    def plot_sum_diff(self):

        if (self.is_attribute("ttheta") & self.is_attribute("intensity_up") & 
            self.is_attribute("intensity_up_sigma") &
            self.is_attribute("intensity_down") & 
            self.is_attribute("intensity_down_sigma")):
            fig, axs = plt.subplots(2)
            ax_1, ax_2 = axs
            ax_1.set_title("Unpolarized intensity: I_up + I_down")
            ax_1.set_xlabel("2 theta (degrees)")
            ax_1.set_ylabel('Intensity')            

            ax_2.set_title("Polarized intensity: I_up - I_down")
            ax_2.set_xlabel("2 theta (degrees)")
            ax_2.set_ylabel('Intensity')            
            
            np_tth = numpy.array(self.ttheta, dtype=float)
            np_up = numpy.array(self.intensity_up, dtype=float)
            np_sup = numpy.array(self.intensity_up_sigma, dtype=float)
            np_down = numpy.array(self.intensity_down, dtype=float)
            np_sdown = numpy.array(self.intensity_down_sigma, dtype=float)
            np_sum = np_up + np_down
            np_ssum = numpy.sqrt(numpy.square(np_sup)+numpy.square(np_sdown))
            ax_1.plot(np_tth, np_sum, "k-", alpha=0.2)
            ax_1.errorbar(np_tth, np_sum, yerr=np_ssum, fmt="ko", alpha=0.2,
                        label="experiment")
            ax_1.legend(loc='upper right')

            np_diff = np_up - np_down
            np_sdiff = numpy.sqrt(numpy.square(np_sup)+numpy.square(np_sdown))
            ax_2.plot([np_tth.min(), np_tth.max()], [0., 0.], "k:")
            ax_2.plot(np_tth, np_diff, "k-", alpha=0.2)
            ax_2.errorbar(np_tth, np_diff, yerr=np_sdiff, fmt="ko", alpha=0.2,
                            label="experiment")
            ax_2.legend(loc='upper right')
            fig.tight_layout()
            return (fig, ax_1)
        else:
            return self.plot_sum()

    
    def plot_sum(self):
        """Plot experimental unpolarized intensity vs. 2 theta (degrees)
        """
        fig, ax = plt.subplots()
        ax.set_title("Unpolarized intensity: I_up + I_down")
        ax.set_xlabel("2 theta (degrees)")
        ax.set_ylabel('Intensity')

        if (self.is_attribute("ttheta") & self.is_attribute("intensity_up") & 
            self.is_attribute("intensity_up_sigma") &
            self.is_attribute("intensity_down") & 
            self.is_attribute("intensity_down_sigma")):
            np_tth = numpy.array(self.ttheta, dtype=float)
            np_up = numpy.array(self.intensity_up, dtype=float)
            np_sup = numpy.array(self.intensity_up_sigma, dtype=float)
            np_down = numpy.array(self.intensity_down, dtype=float)
            np_sdown = numpy.array(self.intensity_down_sigma, dtype=float)
            np_sum = np_up + np_down
            np_ssum = numpy.sqrt(numpy.square(np_sup)+numpy.square(np_sdown))
            ax.plot(np_tth, np_sum, "k-", alpha=0.2)
            ax.errorbar(np_tth, np_sum, yerr=np_ssum, fmt="ko", alpha=0.2,
                        label="experiment")
        elif (self.is_attribute("ttheta") & self.is_attribute("intensity") & 
              self.is_attribute("intensity_sigma")):
            np_tth = numpy.array(self.ttheta, dtype=float)
            np_sum = numpy.array(self.intensity, dtype=float)
            np_ssum = numpy.array(self.intensity_sigma, dtype=float)
            ax.plot(np_tth, np_sum, "k-", alpha=0.2)
            ax.errorbar(np_tth, np_sum, yerr=np_ssum, fmt="ko", alpha=0.2,
                        label="experiment")
        ax.legend(loc='upper right')
        fig.tight_layout()
        return (fig, ax)

    def plot_diff(self):
        """Plot experimental polarized intensity vs. 2 theta (degrees)
        """
        if not(self.is_attribute("ttheta") & self.is_attribute("intensity_up") & 
               self.is_attribute("intensity_up_sigma") &
               self.is_attribute("intensity_down") & 
               self.is_attribute("intensity_down_sigma")):
            return
        fig, ax = plt.subplots()
        ax.set_title("Polarized intensity: I_up - I_down")
        ax.set_xlabel("2 theta (degrees)")
        ax.set_ylabel('Intensity')
            
        np_tth = numpy.array(self.ttheta, dtype=float)
        np_up = numpy.array(self.intensity_up, dtype=float)
        np_sup = numpy.array(self.intensity_up_sigma, dtype=float)
        np_down = numpy.array(self.intensity_down, dtype=float)
        np_sdown = numpy.array(self.intensity_down_sigma, dtype=float)
        np_diff = np_up - np_down
        np_sdiff = numpy.sqrt(numpy.square(np_sup)+numpy.square(np_sdown))
        ax.plot([np_tth.min(), np_tth.max()], [0., 0.], "k:")
        ax.plot(np_tth, np_diff, "k-", alpha=0.2)
        ax.errorbar(np_tth, np_diff, yerr=np_sdiff, fmt="ko", alpha=0.2,
                        label="experiment")
        ax.legend(loc='upper right')
        fig.tight_layout()
        return (fig, ax)


# s_cont = """
#  loop_
#  _pd_meas_ttheta
#  _pd_meas_intensity_up
#  _pd_meas_intensity_up_sigma
#  _pd_meas_intensity_down
#  _pd_meas_intensity_down_sigma
#   4.00   465.80000   128.97000   301.88000   129.30000
#   4.20   323.78000   118.22000   206.06000   120.00000
# """

# obj = PdMeasL.from_cif(s_cont)
# print(obj, end="\n\n")
