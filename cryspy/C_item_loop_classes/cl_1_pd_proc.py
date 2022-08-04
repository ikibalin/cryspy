from typing import NoReturn
import numpy
import matplotlib.pyplot as plt
import matplotlib.widgets as wdg
import scipy
import scipy.optimize

from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN

from .cl_1_pd_background import PdBackgroundL

class PdProc(ItemN):
    """Calculated one point for diffraction pattern.

    This section contains the diffraction data set after processing
    and application of correction terms. If the data set is
    reprocessed, this section may be replaced (with the addition of
    a new _pd_block_id entry).
    """
    ATTR_MANDATORY_NAMES = ("ttheta", )
    ATTR_MANDATORY_TYPES = (float, )
    ATTR_MANDATORY_CIF = ("2theta", )

    ATTR_OPTIONAL_NAMES = (
        "ttheta_corrected", "d_spacing", "intensity_plus_net",
        "intensity_minus_net", "intensity_plus_total", "intensity_minus_total", 
        "intensity_bkg_calc", "intensity_net", "intensity_total",
        "intensity_diff_total", "intensity_plus", "intensity_plus_sigma",
        "intensity_minus", "intensity_minus_sigma", "intensity",
        "intensity_sigma", "excluded")
    ATTR_OPTIONAL_TYPES = (float, float, float, float, float, float, float,
                           float, float, float, float, float, float, float, 
                           float, float, bool)
    ATTR_OPTIONAL_CIF = (
        "2theta_corrected", "d_spacing", "intensity_plus_net",
        "intensity_minus_net", "intensity_plus_total", "intensity_minus_total", 
        "intensity_bkg_calc", "intensity_net", "intensity_total",
        "intensity_diff_total", "intensity_plus", "intensity_plus_sigma",
        "intensity_minus", "intensity_minus_sigma", "intensity",
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
    ATTR_CONSTR_MARK = tuple([f"{_h:}_mark" for _h in ATTR_REF])

    # formats if cif format
    D_FORMATS = {"ttheta": "{:.2f}",
        "ttheta_corrected": "{:.2f}", "d_spacing": "{:.5f}",
        "intensity_plus_net": "{:.2f}", "intensity_minus_net": "{:.2f}",
        "intensity_plus_total": "{:.2f}", "intensity_minus_total": "{:.2f}",
        "intensity_bkg_calc": "{:.2f}", "intensity_net": "{:.2f}",
        "intensity_total": "{:.2f}", "intensity_diff_total": "{:.2f}",
        "intensity_minus": "{:.2f}", "intensity_minus_sigma": "{:.2f}",
        "intensity": "{:.2f}", "intensity_sigma": "{:.2f}"}

    # constraints on the parameters
    D_CONSTRAINTS = {}

    # default values for the parameters
    D_DEFAULT = {"excluded": False}
    for key in ATTR_SIGMA:
        D_DEFAULT[key] = 0.
    for key in (ATTR_CONSTR_FLAG + ATTR_REF_FLAG):
        D_DEFAULT[key] = False
    for key in ATTR_CONSTR_MARK:
        D_DEFAULT[key] = ""

    PREFIX = "pd_proc"

    def __init__(self, **kwargs) -> NoReturn:
        super(PdProc, self).__init__()

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


class PdProcL(LoopN):
    """Calculated diffraction pattern.

    This section contains the diffraction data set after processing
    and application of correction terms. If the data set is
    reprocessed, this section may be replaced (with the addition of
    a new _pd_block_id entry).

    """
    ITEM_CLASS = PdProc
    ATTR_INDEX = "ttheta"
    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(PdProcL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name
    
    def calc_chi_sq(self):
        pass

    def plots(self):
        return [self.plot_sum_diff()]
    
    def plot_sum_diff(self):
        """Plot unpolarized intensity vs. 2 theta
        """

        if (self.is_attribute("ttheta") and self.is_attribute("intensity_plus") and 
                self.is_attribute("intensity_plus_sigma") and
                self.is_attribute("intensity_minus") and 
                self.is_attribute("intensity_minus_sigma") and
                self.is_attribute("intensity_plus_net") and
                self.is_attribute("intensity_minus_net") and
                self.is_attribute("intensity_bkg_calc")):

            fig = plt.figure(constrained_layout=False)
            gs = fig.add_gridspec(nrows=3, ncols=1, hspace=0, height_ratios=[7,1,7])

            ax_1 = fig.add_subplot(gs[0, 0])
            ax_3 = fig.add_subplot(gs[1, 0], sharex=ax_1)
            ax_3.set_yticks([])
            ax_2 = fig.add_subplot(gs[2, 0], sharex=ax_1)

            # fig, axs = plt.subplots(nrows = 2, sharex=True)
            # ax_1, ax_2 = axs
            ax_2.set_xlabel(r"$2\theta$ (degrees)")
            ax_2.set_ylabel('Intensity (arb.u.)')
            ax_1.set_ylabel('Intensity (arb.u.)')

            np_excl = numpy.array(self.excluded, dtype=bool)
            np_notexcl = numpy.logical_not(np_excl)
            np_tth = numpy.array(self.ttheta, dtype=float)
            np_up = numpy.array(self.intensity_plus, dtype=float)
            np_sup = numpy.array(self.intensity_plus_sigma, dtype=float)
            np_up_mod = numpy.array(self.intensity_plus_net, dtype=float)
            np_down = numpy.array(self.intensity_minus, dtype=float)
            np_sdown = numpy.array(self.intensity_minus_sigma, dtype=float)
            np_down_mod = numpy.array(self.intensity_minus_net, dtype=float)
            np_bkg = numpy.array(self.intensity_bkg_calc, dtype=float)
            np_sum = np_up + np_down 
            np_sum_mod = np_up_mod + np_down_mod + np_bkg
            np_diff = np_up - np_down 
            np_diff_mod = np_up_mod - np_down_mod
            np_ssum = numpy.sqrt(numpy.square(np_sup)+numpy.square(np_sdown))
            chi_sq_points = numpy.nansum(numpy.square((np_sum - np_sum_mod)/np_ssum)[np_notexcl])/numpy.sum(np_notexcl) 
            chi_sq_diff_points = numpy.nansum(numpy.square((np_diff - np_diff_mod)/np_ssum))/np_diff_mod.size

            ax_1.set_title(f"Unpolarized ($\chi^2/n = ${chi_sq_points:.2f}) and polarized ($\chi^2/n = ${chi_sq_diff_points:.2f}) signals")
            ax_1.fill_between(np_tth, (np_sum-np_ssum), (np_sum+np_ssum), where=np_notexcl, color="k", alpha=0.4, label="experiment")
            if numpy.any(np_excl):
                ax_1.fill_between(np_tth, (np_sum-np_ssum), (np_sum+np_ssum), where=np_excl, color="r", alpha=0.5, label="excluded")
            # ax_button = plt.axes([0.25, 0.1, 0.08, 0.05])
            # grid_button = wdg.Button(ax_button, 'Grid', color='white', hovercolor='grey')
            # ax_1.errorbar(np_tth[np_notexcl], np_sum[np_notexcl], yerr=np_ssum[np_notexcl], color="k", alpha=0.2, label="experiment")
            # ax_1.errorbar(np_tth[np_excl], np_sum[np_excl], yerr=np_ssum[np_excl], color="r", alpha=0.2, label="excluded")

            y_min_d, y_max_d = ax_1.get_ylim()
            param = y_min_d-((np_sum - np_sum_mod)[np_notexcl]).max()
            coeff = np_notexcl.astype(int)
            ax_1.plot(np_tth, coeff*(np_sum - np_sum_mod)+param, "r-", alpha=0.5, label="difference")
            ax_1.plot(np_tth, np_sum_mod, "-", color="blue", label="model", linewidth=1)
            ax_1.plot([np_tth.min(), np_tth.max()], [param, param], "k:")
            ax_1.plot(np_tth, np_bkg, "b:", label="background")
            ax_1.legend(loc='upper right')

            ax_2.plot([np_tth.min(), np_tth.max()], [0., 0.], "b:")

            ax_2.fill_between(np_tth, np_diff-np_ssum, np_diff+np_ssum, color="k", alpha=0.4, label="experiment")

            # ax_2.errorbar(np_tth, np_diff, yerr=np_ssum, color="k", alpha=0.2, label="experiment") # fmt="k."

            y_min_d, y_max_d = ax_2.get_ylim()
            param = y_min_d-(np_diff-np_diff_mod).max()
            ax_2.plot([np_tth.min(), np_tth.max()], [param, param], "k:")
            ax_2.plot(np_tth, np_diff-np_diff_mod+param, "r-", alpha=0.5, label="difference")
            ax_2.plot(np_tth, np_diff_mod, "b-", label="model", linewidth=1)
            ax_2.legend(loc='upper right')

        elif (self.is_attribute("ttheta") and 
                self.is_attribute("intensity") and
                self.is_attribute("intensity_plus_net") and
                self.is_attribute("intensity_minus_net") and
                self.is_attribute("intensity_bkg_calc")):

            fig = plt.figure(constrained_layout=False)
            gs = fig.add_gridspec(nrows=2, ncols=1, hspace=0, height_ratios=[14,1])

            ax_1 = fig.add_subplot(gs[0, 0])
            ax_3 = fig.add_subplot(gs[1, 0], sharex=ax_1)
            ax_3.set_yticks([])

            ax_1.set_xlabel(r"$2\theta$ (degrees)")
            ax_1.set_ylabel('Intensity (arb.u.)')

            np_excl = numpy.array(self.excluded, dtype=bool)
            np_notexcl = numpy.logical_not(np_excl)
            np_tth = numpy.array(self.ttheta, dtype=float)
            np_sum = numpy.array(self.intensity, dtype=float)
            np_ssum = numpy.array(self.intensity_sigma, dtype=float)

            np_up_mod = numpy.array(self.intensity_plus_net, dtype=float)
            np_down_mod = numpy.array(self.intensity_minus_net, dtype=float)
            np_bkg = numpy.array(self.intensity_bkg_calc, dtype=float)
            np_sum_mod = np_up_mod + np_down_mod + np_bkg

            chi_sq_points = numpy.nansum(numpy.square((np_sum - np_sum_mod)/np_ssum)[np_notexcl])/numpy.sum(np_notexcl) 
            ax_1.set_title(f"Unpolarized signal $\chi^2/n = ${chi_sq_points:.2f}")

            ax_1.fill_between(np_tth, (np_sum-np_ssum), (np_sum+np_ssum), where=np_notexcl, color="k", alpha=0.4, label="experiment")
            ax_1.fill_between(np_tth, (np_sum-np_ssum), (np_sum+np_ssum), where=np_excl, color="r", alpha=0.5, label="excluded")

            # ax_1.errorbar(np_tth[np_notexcl], np_sum[np_notexcl], yerr=np_ssum[np_notexcl], fmt="k.", alpha=0.2, label="experiment")
            # ax_1.errorbar(np_tth[np_excl], np_sum[np_excl], yerr=np_ssum[np_excl], fmt="rs", alpha=0.2, label="excluded")

            y_min_d, y_max_d = ax_1.get_ylim()
            param = y_min_d-(np_sum - np_sum_mod).max()
            coeff = np_notexcl.astype(int)

            ax_1.plot([np_tth.min(), np_tth.max()], [param, param], "k:")
            ax_1.plot(np_tth, coeff*(np_sum - np_sum_mod)+param, "r-", alpha=0.5, label="difference")
            ax_1.plot(np_tth, np_sum_mod, "b-", label="model", linewidth=1)

            if (self.is_attribute("ttheta") and
                    self.is_attribute("intensity_bkg_calc")):
                np_tth = numpy.array(self.ttheta, dtype=float)
                np_bkg = numpy.array(self.intensity_bkg_calc, dtype=float)
                ax_1.plot(np_tth, np_bkg, "b:", label="background")
            ax_1.legend(loc='upper right')
            
        fig.tight_layout()
        return (fig, ax_1)

    def plot_diff(self):
        """Plot polarized intensity vs. 2 theta
        """

        if not(self.is_attribute("ttheta") &
               self.is_attribute("intensity_plus") & 
               self.is_attribute("intensity_plus_sigma") &
               self.is_attribute("intensity_minus") & 
               self.is_attribute("intensity_minus_sigma") &
               self.is_attribute("intensity_plus_net") &
               self.is_attribute("intensity_minus_net")):
            return 
        fig, ax = plt.subplots()
        ax.set_title("Polarized intensity: I_up - I_down")
        ax.set_xlabel("2 theta (degrees)")
        ax.set_ylabel('Intensity')
            
        np_tth = numpy.array(self.ttheta, dtype=float)
        np_up = numpy.array(self.intensity_plus, dtype=float)
        np_sup = numpy.array(self.intensity_plus_sigma, dtype=float)
        np_up_mod = numpy.array(self.intensity_plus_net, dtype=float)
        np_down = numpy.array(self.intensity_minus, dtype=float)
        np_sdown = numpy.array(self.intensity_minus_sigma, dtype=float)
        np_down_mod = numpy.array(self.intensity_minus_net, dtype=float)
        np_diff = np_up - np_down
        np_diff_mod = np_up_mod - np_down_mod
        np_sdiff = numpy.sqrt(numpy.square(np_sup)+numpy.square(np_sdown))

        chi_sq_points = numpy.nansum(numpy.square((np_diff - np_diff_mod)/np_sdiff))/(np_sdiff.size) 
        ax.set_title(f"Polarized signal $\chi^2/n = ${chi_sq_points:.2f}")

        ax.plot([np_tth.min(), np_tth.max()], [0., 0.], "b:")
        ax.fill_between(np_tth, (np_diff-np_sdiff), (np_diff+np_sdiff), color="k", alpha=0.4, label="experiment")        
        # ax.errorbar(np_tth, np_diff, yerr=np_sdiff, fmt="ko", alpha=0.2, label="experiment")
        
        y_min_d, y_max_d = ax.get_ylim()
        param = y_min_d-(np_diff-np_diff_mod).max()
        
        ax.plot([np_tth.min(), np_tth.max()], [param, param], "k:")
        ax.plot(np_tth, np_diff-np_diff_mod+param, "r-", alpha=0.5,
                    label="difference")
        ax.plot(np_tth, np_diff_mod, "b-",
                    label="model", linewidth=2)
        ax.legend(loc='upper right')
        fig.tight_layout()
        return (fig, ax)

    def estimate_background(self,
            pd_background: PdBackgroundL, offset_ttheta: float = 0):
        ttheta = numpy.array(self.ttheta, dtype=float)
        if self.is_attribute("intensity_plus"):
            signal_exp_sum = numpy.array(self.intensity_plus, dtype=float) + \
                numpy.array(self.intensity_minus, dtype=float)
            signal_exp_sum_sigma = numpy.sqrt(
                numpy.square(numpy.array(self.intensity_plus_sigma, dtype=float)) +
                numpy.square(numpy.array(self.intensity_minus_sigma, dtype=float)))
        else:
            signal_exp_sum = numpy.array(self.intensity, dtype=float)
            signal_exp_sum_sigma = numpy.array(self.intensity_sigma, dtype=float)

        signal_sum_net = numpy.array(self.intensity_plus_net, dtype=float) + \
            numpy.array(self.intensity_minus_net, dtype=float) 
        included_points = numpy.logical_not(numpy.array(self.excluded, dtype=bool))
        model_bkgr = numpy.array(self.intensity_bkg_calc, dtype=float)

        variable_names = pd_background.get_variable_names()
        param_0 = [pd_background.get_variable_by_name(name) for name in variable_names]
        def temp_func(param):
            for name, value in zip(variable_names, param):
                pd_background.set_variable_by_name(name, value)
            model_bkgr = pd_background.interpolate_by_points(ttheta+offset_ttheta)
            chi_sq_sum = numpy.nansum(numpy.square((signal_exp_sum - signal_sum_net - model_bkgr)/signal_exp_sum_sigma)[included_points])
            return chi_sq_sum
        
        res = scipy.optimize.minimize(temp_func, param_0, method="Nelder-Mead")
        
        self.numpy_intensity_bkg_calc = pd_background.interpolate_by_points(ttheta+offset_ttheta)
        self.numpy_to_items()
        return res