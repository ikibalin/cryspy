from typing import NoReturn
import numpy
import matplotlib.pyplot as plt

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN


class TOFProc(ItemN):
    """Calculated intensity for time-of-flight experiment.

    Attributes
    ----------
        - time_corrected, d_spacing, intensity_up_net,
          intensity_down_net, intensity_up_total, intensity_down_total, 
          intensity_bkg_calc, intensity_net, intensity_total,
          intensity_diff_total, intensity_up, intensity_up_sigma,
          intensity_down, intensity_down_sigma, intensity,
          intensity_sigma (optional)

    time is given in microseconds.
    """
    ATTR_MANDATORY_NAMES = ("time", )
    ATTR_MANDATORY_TYPES = (float, )
    ATTR_MANDATORY_CIF = ("time", )

    ATTR_OPTIONAL_NAMES = (
        "time_corrected", "d_spacing", "intensity_up_net",
        "intensity_down_net", "intensity_up_total", "intensity_down_total", 
        "intensity_bkg_calc", "intensity_net", "intensity_total",
        "intensity_diff_total", "intensity_up", "intensity_up_sigma",
        "intensity_down", "intensity_down_sigma", "intensity",
        "intensity_sigma", "excluded")
    ATTR_OPTIONAL_TYPES = (float, float, float, float, float, float, float,
                           float, float, float, float, float, float, float, 
                           float, float, bool)
    ATTR_OPTIONAL_CIF = (
        "time_corrected", "d_spacing", "intensity_up_net",
        "intensity_down_net", "intensity_up_total", "intensity_down_total", 
        "intensity_bkg_calc", "intensity_net", "intensity_total",
        "intensity_diff_total", "intensity_up", "intensity_up_sigma",
        "intensity_down", "intensity_down_sigma", "intensity",
        "intensity_sigma", "excluded")

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
    D_FORMATS = {
        "time_corrected": "{:.2f}", "d_spacing": "{:.5f}",
        "intensity_up_net": "{:.2f}", "intensity_down_net": "{:.2f}",
        "intensity_up_total": "{:.2f}", "intensity_down_total": "{:.2f}",
        "intensity_bkg_calc": "{:.2f}", "intensity_net": "{:.2f}",
        "intensity_total": "{:.2f}", "intensity_diff_total": "{:.2f}",
        "intensity_down": "{:.2f}", "intensity_down_sigma": "{:.2f}",
        "intensity": "{:.2f}", "intensity_sigma": "{:.2f}"}

    # constraints on the parameters
    D_CONSTRAINTS = {}

    # default values for the parameters
    D_DEFAULT = {"excluded": False}
    for key in ATTR_SIGMA:
        D_DEFAULT[key] = 0.
    for key in (ATTR_CONSTR_FLAG + ATTR_REF_FLAG):
        D_DEFAULT[key] = False

    PREFIX = "tof_proc"

    def __init__(self, **kwargs) -> NoReturn:
        super(TOFProc, self).__init__()

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


class TOFProcL(LoopN):
    """Calculated intensity for time-of-flight experiment.
    """
    ITEM_CLASS = TOFProc
    ATTR_INDEX = "time"
    def __init__(self, loop_name = None) -> NoReturn:
        super(TOFProcL, self).__init__()
        self.__dict__["items"] = []
        self.__dict__["loop_name"] = loop_name

    def plots(self):
        return [self.plot_sum(), self.plot_diff()]
    
    def plot_sum(self):
        """Plot unpolarized intensity vs. time
        """
        fig, ax = plt.subplots()
        ax.set_title("Unpolarized intensity: I_up + I_down")
        ax.set_xlabel("Time (microseconds)")
        ax.set_ylabel('Intensity')

        if (self.is_attribute("time") & self.is_attribute("intensity_up") & 
                self.is_attribute("intensity_up_sigma") &
                self.is_attribute("intensity_down") & 
                self.is_attribute("intensity_down_sigma") &
                self.is_attribute("intensity_up_total") &
                self.is_attribute("intensity_down_total")):
            np_excl = numpy.array(self.excluded, dtype=bool)
            np_notexcl = numpy.logical_not(np_excl)
            np_time = numpy.array(self.time, dtype=float)
            np_up = numpy.array(self.intensity_up, dtype=float)
            np_sup = numpy.array(self.intensity_up_sigma, dtype=float)
            np_up_mod = numpy.array(self.intensity_up_total, dtype=float)
            np_down = numpy.array(self.intensity_down, dtype=float)
            np_sdown = numpy.array(self.intensity_down_sigma, dtype=float)
            np_down_mod = numpy.array(self.intensity_down_total, dtype=float)
            np_sum = np_up + np_down
            np_sum_mod = np_up_mod + np_down_mod
            np_ssum = numpy.sqrt(numpy.square(np_sup)+numpy.square(np_sdown))
            ax.plot(np_time, np_sum_mod, "k-", label="model")
            ax.errorbar(np_time[np_notexcl], np_sum[np_notexcl], yerr=np_ssum[np_notexcl], fmt="ko", alpha=0.2, label="experiment")
            ax.errorbar(np_time[np_excl], np_sum[np_excl], yerr=np_ssum[np_excl], fmt="rs", alpha=0.2, label="excluded")

            y_min_d, y_max_d = ax.get_ylim()
            param = y_min_d-(np_sum - np_sum_mod).max()
            coeff = np_notexcl.astype(int)

            ax.plot([np_time.min(), np_time.max()], [param, param], "k:")
            ax.plot(np_time, coeff*(np_sum - np_sum_mod)+param, "r-", alpha=0.7,
                    label="difference")
        elif (self.is_attribute("time") & self.is_attribute("intensity") & 
              self.is_attribute("intensity_total") &
              self.is_attribute("intensity_sigma")):
            np_excl = numpy.array(self.excluded, dtype=bool)
            np_notexcl = numpy.logical_not(np_excl)
            np_time = numpy.array(self.time, dtype=float)
            np_sum = numpy.array(self.intensity, dtype=float)
            np_sum_mod = numpy.array(self.intensity_total, dtype=float)
            np_ssum = numpy.array(self.intensity_sigma, dtype=float)
            ax.plot(np_time, np_sum_mod, "k-", label="model")
            ax.errorbar(np_time[np_notexcl], np_sum[np_notexcl], yerr=np_ssum[np_notexcl], fmt="ko", alpha=0.2, label="experiment")
            ax.errorbar(np_time[np_excl], np_sum[np_excl], yerr=np_ssum[np_excl], fmt="rs", alpha=0.2, label="excluded")

            y_min_d, y_max_d = ax.get_ylim()
            param = y_min_d-(np_sum - np_sum_mod).max()
            coeff = np_notexcl.astype(int)

            ax.plot([np_time.min(), np_time.max()], [param, param], "k:")
            ax.plot(np_time, coeff*(np_sum - np_sum_mod)+param, "r-", alpha=0.7,
                    label="difference")
        ax.legend(loc='upper right')
        fig.tight_layout()
        return (fig, ax)

    def plot_diff(self):
        """Plot polarized intensity vs. time
        """
        if not(self.is_attribute("time") & self.is_attribute("intensity_up") & 
               self.is_attribute("intensity_up_sigma") &
               self.is_attribute("intensity_down") & 
               self.is_attribute("intensity_down_sigma") &
               self.is_attribute("intensity_up_total") &
               self.is_attribute("intensity_down_total")):
            return
        fig, ax = plt.subplots()
        ax.set_title("Polarized intensity: I_up - I_down")
        ax.set_xlabel("Time (microseconds)")
        ax.set_ylabel('Intensity')
            
        np_time = numpy.array(self.time, dtype=float)
        np_up = numpy.array(self.intensity_up, dtype=float)
        np_sup = numpy.array(self.intensity_up_sigma, dtype=float)
        np_up_mod = numpy.array(self.intensity_up_total, dtype=float)
        np_down = numpy.array(self.intensity_down, dtype=float)
        np_sdown = numpy.array(self.intensity_down_sigma, dtype=float)
        np_down_mod = numpy.array(self.intensity_down_total, dtype=float)
        np_diff = np_up - np_down
        np_diff_mod = np_up_mod - np_down_mod
        np_sdiff = numpy.sqrt(numpy.square(np_sup)+numpy.square(np_sdown))

        ax.plot([np_time.min(), np_time.max()], [0., 0.], "b:")
        ax.plot(np_time, np_diff_mod, "k-",
                    label="model")
        ax.errorbar(np_time, np_diff, yerr=np_sdiff, fmt="ko", alpha=0.2,
                    label="experiment")

        y_min_d, y_max_d = ax.get_ylim()
        param = y_min_d-(np_diff-np_diff_mod).max()

        ax.plot([np_time.min(), np_time.max()], [param, param], "k:")
        ax.plot(np_time, np_diff-np_diff_mod+param, "r-", alpha=0.7,
                    label="difference")
        ax.legend(loc='upper right')
        fig.tight_layout()
        return (fig, ax)


# s_cont = """
#   loop_
#   _tof_proc_time
#   _tof_proc_time_corrected
#   _tof_proc_d_spacing
#   _tof_proc_intensity_up_net
#   _tof_proc_intensity_down_net
#   _tof_proc_intensity_up_total
#   _tof_proc_intensity_down_total
#   _tof_proc_intensity_bkg_calc
#   _tof_proc_intensity_up
#   _tof_proc_intensity_up_sigma
#   _tof_proc_intensity_down
#   _tof_proc_intensity_down_sigma
#   4.00  3.9  11.2   400.00   317.00   460.00   377.00   60.00    465.80000   128.97000   301.88000   129.30000"""

# obj = TOFProcL.from_cif(s_cont)
# print(obj, end="\n\n")
