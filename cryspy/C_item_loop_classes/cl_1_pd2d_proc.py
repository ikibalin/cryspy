from typing import NoReturn
import numpy
import scipy
from scipy.interpolate import interp2d
from scipy.optimize import minimize

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.widgets as wdg

from cryspy.A_functions_base.powder_diffraction_const_wavelength import \
    calc_ttheta_phi_by_gamma_nu

from cryspy.A_functions_base.function_1_gamma_nu import \
    recal_int_to_gammanu_grid

from cryspy.B_parent_classes.cl_1_item import ItemN

from .cl_1_pd2d_background import Pd2dBackground

na = numpy.newaxis

class Pd2dProc(ItemN):
    """
    Calculated 2D diffraction pattern (polarized neutrons).

    Attributes
    ----------
        - gamma_nu_intensity_plus_net, gamma_nu_intensity_minus_net,
          gamma_nu_intensity_plus_total, gamma_nu_intensity_minus_total,
          gamma_nu_intensity_bkg_calc, gamma_nu_intensity_plus,
          gamma_nu_intensity_plus_sigma, gamma_nu_intensity_minus,
          gamma_nu_intensity_minus_sigma (mandatory)

    Internal
    --------
        - gamma, nu, intensity_plus_net, intensity_minus_net,
          intensity_plus_total, intensity_minus_total, intensity_bkg_calc,
          intensity_plus, intensity_plus_sigma, intensity_minus,
          intensity_minus_sigma

    """
    ATTR_MANDATORY_NAMES = (
        "gamma_nu_intensity_plus_net", "gamma_nu_intensity_minus_net",
        "gamma_nu_intensity_bkg_calc", "gamma_nu_excluded_points", )
    ATTR_MANDATORY_TYPES = (str, str, str, str, )

    ATTR_MANDATORY_CIF = (
        "gamma_nu_intensity_plus_net", "gamma_nu_intensity_minus_net",
        "gamma_nu_intensity_bkg_calc", "gamma_nu_excluded_points", )

    ATTR_OPTIONAL_NAMES = (
        "gamma_nu_intensity_plus", "gamma_nu_intensity_plus_sigma",
        "gamma_nu_intensity_minus", "gamma_nu_intensity_minus_sigma",
        "gamma_nu_intensity", "gamma_nu_intensity_sigma", )
    ATTR_OPTIONAL_TYPES = (str, str, str, str, str, str,  )
    ATTR_OPTIONAL_CIF = (
        "gamma_nu_intensity_plus", "gamma_nu_intensity_plus_sigma",
        "gamma_nu_intensity_minus", "gamma_nu_intensity_minus_sigma",
        "gamma_nu_intensity", "gamma_nu_intensity_sigma", )

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = (
        "gamma", "nu", "intensity_plus_net", "intensity_minus_net",
        "intensity_bkg_calc",
        "intensity_plus", "intensity_plus_sigma",
        "intensity_minus", "intensity_minus_sigma",
        "intensity", "intensity_sigma",
        "excluded_points")
    ATTR_INT_PROTECTED_NAMES = ()

    # parameters considered are refined parameters
    ATTR_REF = ()
    ATTR_SIGMA = tuple([f"{_h:}_sigma" for _h in ATTR_REF])
    ATTR_CONSTR_FLAG = tuple([f"{_h:}_constraint" for _h in ATTR_REF])
    ATTR_REF_FLAG = tuple([f"{_h:}_refinement" for _h in ATTR_REF])
    ATTR_CONSTR_MARK = tuple([f"{_h:}_mark" for _h in ATTR_REF])

    # constraints on the parameters
    D_CONSTRAINTS = {}

    # default values for the parameters
    D_DEFAULT = {}
    for key in ATTR_SIGMA:
        D_DEFAULT[key] = 0.
    for key in (ATTR_CONSTR_FLAG + ATTR_REF_FLAG):
        D_DEFAULT[key] = False
    for key in ATTR_CONSTR_MARK:
        D_DEFAULT[key] = ""

    PREFIX = "pd2d_proc"

    def __init__(self, **kwargs) -> NoReturn:
        super(Pd2dProc, self).__init__()

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

    def form_object(self) -> NoReturn:
        flag = True
        if not(self.is_attribute("gamma_nu_intensity_plus_net")):
            return False
        l_1 = (self.gamma_nu_intensity_plus_net).strip().split("\n")
        l_2 = (self.gamma_nu_intensity_minus_net).strip().split("\n")
        l_3 = (self.gamma_nu_intensity_bkg_calc).strip().split("\n")
        l_8 = (self.gamma_nu_excluded_points).strip().split("\n")

        flag_polarized, flag_unpolarized = False, False
        if self.is_attribute("gamma_nu_intensity_plus"):
            l_4 = (self.gamma_nu_intensity_plus).strip().split("\n")
            l_5 = (self.gamma_nu_intensity_plus_sigma).strip().split("\n")
            l_6 = (self.gamma_nu_intensity_minus).strip().split("\n")
            l_7 = (self.gamma_nu_intensity_minus_sigma).strip().split("\n")
            flag_polarized = True
        elif self.is_attribute("gamma_nu_intensity"):
            l_4 = (self.gamma_nu_intensity).strip().split("\n")
            l_5 = (self.gamma_nu_intensity_sigma).strip().split("\n")
            flag_unpolarized = True
        else:
            return False
        l_gamma = numpy.array([_ for _ in l_1[0].strip().split()[1:]], dtype=float)
        l_nu = []
        ll_intensity_plus_net, ll_intensity_minus_net = [], []
        ll_excluded_points = []
        ll_intensity_bkg_calc = []
        for line_1, line_2, line_3, line_8 in zip(l_1[1:], l_2[1:], l_3[1:], l_8[1:]):
            _l_1 = line_1.strip().split()
            _l_2 = line_2.strip().split()
            _l_3 = line_3.strip().split()
            _l_8 = line_8.strip().split()
            l_nu.append(float(_l_1[0]))
            ll_intensity_plus_net.append(_l_1[1:])
            ll_intensity_minus_net.append(_l_2[1:])
            ll_intensity_bkg_calc.append(_l_3[1:])
            ll_excluded_points.append(_l_8[1:])

        if flag_polarized:
            ll_intensity_plus, ll_intensity_plus_sigma = [], []
            ll_intensity_minus, ll_intensity_minus_sigma = [], []
            for line_4, line_5, line_6, line_7 in zip(l_4[1:], l_5[1:], l_6[1:], l_7[1:]):
                _l_4 = line_4.strip().split()
                _l_5 = line_5.strip().split()
                _l_6 = line_6.strip().split()
                _l_7 = line_7.strip().split()
                ll_intensity_plus.append(_l_4[1:])
                ll_intensity_plus_sigma.append(_l_5[1:])
                ll_intensity_minus.append(_l_6[1:])
                ll_intensity_minus_sigma.append(_l_7[1:])
        elif flag_unpolarized:
            ll_intensity, ll_intensity_sigma = [], []
            for line_4, line_5 in zip(l_4[1:], l_5[1:]):
                _l_4 = line_4.strip().split()
                _l_5 = line_5.strip().split()
                ll_intensity.append(_l_4[1:])
                ll_intensity_sigma.append(_l_5[1:])

        ll_intensity_plus_net = numpy.array(ll_intensity_plus_net, dtype=float).transpose()
        ll_intensity_minus_net = numpy.array(ll_intensity_minus_net, dtype=float).transpose()
        ll_intensity_bkg_calc = numpy.array(ll_intensity_bkg_calc, dtype=float).transpose()
        ll_excluded_points = numpy.array([[hhh == "True" for hhh in hh ]for hh in ll_excluded_points], dtype=bool).transpose()

        if flag_polarized:
            ll_intensity_plus = numpy.array(ll_intensity_plus, dtype=float).transpose()
            ll_intensity_plus_sigma = numpy.array(ll_intensity_plus_sigma, dtype=float).transpose()
            ll_intensity_minus = numpy.array(ll_intensity_minus, dtype=float).transpose()
            ll_intensity_minus_sigma = numpy.array(ll_intensity_minus_sigma, dtype=float).transpose()
        elif flag_unpolarized:
            ll_intensity = numpy.array(ll_intensity, dtype=float).transpose()
            ll_intensity_sigma = numpy.array(ll_intensity_sigma, dtype=float).transpose()

        self.__dict__["gamma"] = l_gamma
        self.__dict__["nu"] = numpy.array(l_nu, dtype=float)
        self.__dict__["intensity_plus_net"] = ll_intensity_plus_net
        self.__dict__["intensity_minus_net"] = ll_intensity_minus_net
        self.__dict__["intensity_bkg_calc"] = ll_intensity_bkg_calc
        self.__dict__["excluded_points"] = ll_excluded_points

        if flag_polarized:
            self.__dict__["intensity_plus"] = ll_intensity_plus
            self.__dict__["intensity_plus_sigma"] = ll_intensity_plus_sigma
            self.__dict__["intensity_minus"] = ll_intensity_minus
            self.__dict__["intensity_minus_sigma"] = ll_intensity_minus_sigma
        elif flag_unpolarized:
            self.__dict__["intensity"] = ll_intensity
            self.__dict__["intensity_sigma"] = ll_intensity_sigma
        return flag

    def form_gamma_nu_intensity_plus_net(self) -> bool:
        if ((self.nu is not None) & (self.gamma is not None) & (self.intensity_plus_net is not None)):
            ls_out = []
            ls_out.append("{:12} ".format(len(self.nu)) + " ".join(["{:6.2f}      ".format(_) for _ in self.gamma]))
            ll_intensity = self.intensity_plus_net
            ll_intensity = [[ll_intensity[_2][_1] for _2 in range(len(ll_intensity))] for _1 in range(len(ll_intensity[0]))]
            for nu, l_intensity in zip(self.nu, ll_intensity):
                ls_out.append("{:12.2f} ".format(nu) + " ".join(["{:12}".format(_) for _ in l_intensity]))
            self.__dict__["gamma_nu_intensity_plus_net"] = "\n".join(ls_out)

    def form_gamma_nu_intensity_minus_net(self) -> bool:
        if ((self.nu is not None) & (self.gamma is not None) & (self.intensity_minus_net is not None)):
            ls_out = []
            ls_out.append("{:12} ".format(len(self.nu)) + " ".join(["{:6.2f}      ".format(_) for _ in self.gamma]))
            ll_intensity = self.intensity_minus_net
            ll_intensity = [[ll_intensity[_2][_1] for _2 in range(len(ll_intensity))] for _1 in range(len(ll_intensity[0]))]
            for nu, l_intensity in zip(self.nu, ll_intensity):
                ls_out.append("{:12.2f} ".format(nu) + " ".join(["{:12}".format(_) for _ in l_intensity]))
            self.__dict__["gamma_nu_intensity_minus_net"] = "\n".join(ls_out)

    def form_gamma_nu_excluded_points(self) -> bool:
        if ((self.nu is not None) & (self.gamma is not None) & (self.excluded_points is not None)):
            ls_out = []
            ls_out.append("{:12} ".format(len(self.nu)) + " ".join(["{:6.2f}      ".format(_) for _ in self.gamma]))
            ll_intensity = self.excluded_points
            ll_intensity = [["True" if ll_intensity[_2][_1] else "False" for _2 in range(len(ll_intensity))] for _1 in range(len(ll_intensity[0]))]
            for nu, l_intensity in zip(self.nu, ll_intensity):
                ls_out.append("{:12.2f} ".format(nu) + " ".join(["{:12}".format(_) for _ in l_intensity]))
            self.__dict__["gamma_nu_excluded_points"] = "\n".join(ls_out)

    def form_gamma_nu_intensity_bkg_calc(self) -> bool:
        if ((self.nu is not None) & (self.gamma is not None) & (self.intensity_bkg_calc is not None)):
            ls_out = []
            ls_out.append("{:12} ".format(len(self.nu)) + " ".join(["{:6.2f}      ".format(_) for _ in self.gamma]))
            ll_intensity = self.intensity_bkg_calc
            ll_intensity = [[ll_intensity[_2][_1] for _2 in range(len(ll_intensity))] for _1 in range(len(ll_intensity[0]))]
            for nu, l_intensity in zip(self.nu, ll_intensity):
                ls_out.append("{:12.2f} ".format(nu) + " ".join(["{:12}".format(_) for _ in l_intensity]))
            self.__dict__["gamma_nu_intensity_bkg_calc"] = "\n".join(ls_out)

    def form_gamma_nu_intensity_plus(self) -> bool:
        if ((self.nu is not None) & (self.gamma is not None) & (self.intensity_plus is not None)):
            ls_out = []
            ls_out.append("{:12} ".format(len(self.nu)) + " ".join(["{:6.2f}      ".format(_) for _ in self.gamma]))
            ll_intensity = self.intensity_plus
            ll_intensity = [[ll_intensity[_2][_1] for _2 in range(len(ll_intensity))] for _1 in range(len(ll_intensity[0]))]
            for nu, l_intensity in zip(self.nu, ll_intensity):
                ls_out.append("{:12.2f} ".format(nu) + " ".join(["{:12}".format(_) for _ in l_intensity]))
            self.__dict__["gamma_nu_intensity_plus"] = "\n".join(ls_out)

    def form_gamma_nu_intensity_plus_sigma(self) -> bool:
        if ((self.nu is not None) & (self.gamma is not None) & (self.intensity_plus_sigma is not None)):
            ls_out = []
            ls_out.append("{:12} ".format(len(self.nu)) + " ".join(["{:6.2f}      ".format(_) for _ in self.gamma]))
            ll_intensity = self.intensity_plus_sigma
            ll_intensity = [[ll_intensity[_2][_1] for _2 in range(len(ll_intensity))] for _1 in range(len(ll_intensity[0]))]
            for nu, l_intensity in zip(self.nu, ll_intensity):
                ls_out.append("{:12.2f} ".format(nu) + " ".join(["{:12}".format(_) for _ in l_intensity]))
            self.__dict__["gamma_nu_intensity_plus_sigma"] = "\n".join(ls_out)

    def form_gamma_nu_intensity_minus(self) -> bool:
        if ((self.nu is not None) & (self.gamma is not None) & (self.intensity_minus is not None)):
            ls_out = []
            ls_out.append("{:12} ".format(len(self.nu)) + " ".join(["{:6.2f}      ".format(_) for _ in self.gamma]))
            ll_intensity = self.intensity_minus
            ll_intensity = [[ll_intensity[_2][_1] for _2 in range(len(ll_intensity))] for _1 in range(len(ll_intensity[0]))]
            for nu, l_intensity in zip(self.nu, ll_intensity):
                ls_out.append("{:12.2f} ".format(nu) + " ".join(["{:12}".format(_) for _ in l_intensity]))
            self.__dict__["gamma_nu_intensity_minus"] = "\n".join(ls_out)

    def form_gamma_nu_intensity_minus_sigma(self) -> bool:
        if ((self.nu is not None) & (self.gamma is not None) & (self.intensity_minus_sigma is not None)):
            ls_out = []
            ls_out.append("{:12} ".format(len(self.nu)) + " ".join(["{:6.2f}      ".format(_) for _ in self.gamma]))
            ll_intensity = self.intensity_minus_sigma
            ll_intensity = [[ll_intensity[_2][_1] for _2 in range(len(ll_intensity))] for _1 in range(len(ll_intensity[0]))]
            for nu, l_intensity in zip(self.nu, ll_intensity):
                ls_out.append("{:12.2f} ".format(nu) + " ".join(["{:12}".format(_) for _ in l_intensity]))
            self.__dict__["gamma_nu_intensity_minus_sigma"] = "\n".join(ls_out)

    def form_gamma_nu_intensity(self) -> bool:
        if ((self.nu is not None) & (self.gamma is not None) & (self.intensity is not None)):
            ls_out = []
            ls_out.append("{:12} ".format(len(self.nu)) + " ".join(["{:6.2f}      ".format(_) for _ in self.gamma]))
            ll_intensity = self.intensity
            ll_intensity = [[ll_intensity[_2][_1] for _2 in range(len(ll_intensity))] for _1 in range(len(ll_intensity[0]))]
            for nu, l_intensity in zip(self.nu, ll_intensity):
                ls_out.append("{:12.2f} ".format(nu) + " ".join(["{:12}".format(_) for _ in l_intensity]))
            self.__dict__["gamma_nu_intensity"] = "\n".join(ls_out)

    def form_gamma_nu_intensity_sigma(self) -> bool:
        if ((self.nu is not None) & (self.gamma is not None) & (self.intensity_sigma is not None)):
            ls_out = []
            ls_out.append("{:12} ".format(len(self.nu)) + " ".join(["{:6.2f}      ".format(_) for _ in self.gamma]))
            ll_intensity = self.intensity_sigma
            ll_intensity = [[ll_intensity[_2][_1] for _2 in range(len(ll_intensity))] for _1 in range(len(ll_intensity[0]))]
            for nu, l_intensity in zip(self.nu, ll_intensity):
                ls_out.append("{:12.2f} ".format(nu) + " ".join(["{:12}".format(_) for _ in l_intensity]))
            self.__dict__["gamma_nu_intensity_sigma"] = "\n".join(ls_out)


    def plots(self):
        return [self.plot_gamma_nu()]
    
    def plot_gamma_nu(self):
        if not(self.is_attribute("gamma") & self.is_attribute("nu") &
               self.is_attribute("intensity_plus_net") & 
               self.is_attribute("intensity_minus_net")):
            return

        gamma = self.gamma 
        nu = self.nu 
        signal_plus = self.intensity_plus_net
        signal_minus = self.intensity_minus_net
        bkg_calc = self.intensity_bkg_calc
        signal_sum = signal_plus + signal_minus + bkg_calc
        signal_difference = signal_plus - signal_minus

        flag_polarized, flag_unpolarized = False, False
        if self.is_attribute("intensity_plus"):
            flag_polarized = True
        if self.is_attribute("intensity"):
            flag_unpolarized = True

        if flag_polarized:
            signal_exp_plus = self.intensity_plus
            signal_exp_minus = self.intensity_minus
            signal_exp_plus_sigma = self.intensity_plus_sigma
            signal_exp_minus_sigma = self.intensity_minus_sigma
            signal_exp_sum = signal_exp_plus + signal_exp_minus
            signal_exp_difference = signal_exp_plus - signal_exp_minus
            signal_exp_sum_sigma_sq = numpy.square(signal_exp_plus_sigma)+numpy.square(signal_exp_minus_sigma)
            signal_exp_sum_sigma = numpy.sqrt(signal_exp_sum_sigma_sq)
        elif flag_unpolarized:
            signal_exp_sum = self.intensity
            signal_exp_sum_sigma = self.intensity_sigma
            signal_exp_sum_sigma_sq = numpy.square(signal_exp_sum_sigma)

        excluded_points = numpy.logical_or(
            self.excluded_points,
            numpy.isnan(signal_exp_sum))
        included_points = numpy.logical_not(excluded_points)

        signal_em_sum = numpy.concatenate([numpy.flip(signal_sum, axis=1), signal_exp_sum], axis=1)
        if flag_polarized:
            signal_em_difference = numpy.concatenate([numpy.flip(signal_difference, axis=1), signal_exp_difference], axis=1)

        max_val = max([numpy.nanmax(signal_exp_sum[included_points]), numpy.nanmax(signal_sum[included_points])])
        min_val = numpy.nanmin(signal_em_sum)
        max_val = max_val - (max_val-min_val)*0.25
        n_sigma_threshold = 20.
        
        alpha_excl_points = numpy.zeros(signal_exp_sum.shape, dtype=int)
        alpha_excl_points[excluded_points]=200
        alpha_em_sum = numpy.concatenate([numpy.zeros_like(alpha_excl_points), alpha_excl_points], axis=1).transpose()

        zz_sum = numpy.stack([
            numpy.concatenate([numpy.zeros_like(alpha_excl_points), alpha_excl_points], axis=1).transpose(),
            numpy.concatenate([numpy.zeros_like(alpha_excl_points), alpha_excl_points], axis=1).transpose(),
            numpy.concatenate([numpy.zeros_like(alpha_excl_points), alpha_excl_points], axis=1).transpose(),
            alpha_em_sum],axis=2)

        extent = [gamma.min(), gamma.max(), nu.min(), nu.max()]
        cmap_sum = plt.get_cmap("turbo") # BuPu

        if flag_polarized:
            fig, axs = plt.subplots(2,2, sharex=True)
            ax1, ax2, ax3, ax4 = axs[0,0], axs[0,1], axs[1,0], axs[1,1]
        elif flag_unpolarized:
            fig, axs = plt.subplots(2, 1, sharex=True)
            ax1, ax3= axs[0], axs[1]

        norm_1 = matplotlib.colors.Normalize(vmax=max_val, vmin=min_val)
        img_unpol = ax1.imshow(signal_em_sum.transpose(), origin="lower", aspect="auto", cmap=cmap_sum, norm= norm_1, extent=extent)
        ax1.imshow(zz_sum, origin="lower", aspect="auto", alpha=0.8, extent=extent)

        ax1.set_xticks([])
        ax1.set_xlabel(r"$\gamma$ (deg.)")
        ax1.set_yticks([])
        ax1.set_ylabel("  Model       Experiment")

        if flag_polarized:
            cmap_difference = plt.get_cmap("turbo") # BrBG
            max_val = numpy.nanmax(numpy.abs(signal_em_difference))*0.50
            norm_2 = matplotlib.colors.Normalize(vmax=max_val, vmin=-max_val)
            ax2.imshow(signal_em_difference.transpose(), origin="lower", aspect="auto",
                      norm=norm_2,
                      cmap=cmap_difference,
                      extent=extent)

            ax2.set_xticks([])
            ax2.set_xlabel(r"$\gamma$ (deg.)")
            ax2.set_yticks([])

        ax1.set_ylabel("  Model       Experiment")

        chi_sq_per_n_sum = numpy.nansum(numpy.square(
                (signal_exp_sum-signal_sum)*included_points/signal_exp_sum_sigma))/numpy.sum(included_points)

        ax1.set_title(r"Unpolarized signal $\chi^2/n=$"+f"{chi_sq_per_n_sum:.2f}")


        ttheta, signal_projection_sum, signal_projection_exp_sum, signal_projection_exp_sum_sigma, \
            signal_projection_difference, signal_projection_exp_difference, signal_projection_exp_difference_sigma = \
            self.calc_projections_sum_difference()

        ax3.fill_between(
            ttheta[:gamma.size],(signal_projection_exp_sum-signal_projection_exp_sum_sigma)[:gamma.size],
            (signal_projection_exp_sum+signal_projection_exp_sum_sigma)[:gamma.size], color="k", alpha=0.4, label="experiment")

        # ax3.errorbar(
        #     ttheta[:gamma.size], signal_projection_exp_sum[:gamma.size],
        #     yerr=signal_projection_exp_sum_sigma[:gamma.size], fmt="ko", alpha=0.2, label="experiment")
        ax3.plot(ttheta[:gamma.size], signal_projection_sum[:gamma.size], "b-", label="model", linewidth=1)

        y_min_d, y_max_d = ax3.get_ylim()
        param = y_min_d-numpy.nanmax((signal_projection_exp_sum-signal_projection_sum)[:gamma.size])
        ax3.plot([ttheta[:gamma.size].min(), ttheta[:gamma.size].max()], [param, param], "k:")
        ax3.plot(ttheta[:gamma.size], (signal_projection_exp_sum-signal_projection_sum)[:gamma.size]+param,
            "r-", alpha=0.5) # , label="difference"

        # ax3.plot(ttheta[:gamma.size], signal_projection_exp_sum[:gamma.size], "b.", label="experiment")
        # ax3.plot(ttheta[:gamma.size], (signal_projection_exp_sum-signal_projection_sum)[:gamma.size], "r-")
        ax3.legend()
        ax3.set_xlabel(r"$2\theta$ (deg.)")
        # ax3.get_yaxis().set_visible(False)

        if flag_polarized:
            chi_sq_per_n_difference = numpy.nansum(numpy.square(
                    (signal_exp_difference-signal_difference)/signal_exp_sum_sigma))/numpy.product(signal_exp_difference.shape)
            ax2.set_title(r"Polarized signal $\chi^2/n=$"+f"{chi_sq_per_n_difference:.2f}")

            ax4.fill_between(
                ttheta[:gamma.size], (signal_projection_exp_difference-signal_projection_exp_difference_sigma)[:gamma.size],
                (signal_projection_exp_difference+signal_projection_exp_difference_sigma)[:gamma.size], color="k", alpha=0.4, label="experiment")
            # ax4.errorbar(
            #     ttheta[:gamma.size], signal_projection_exp_difference[:gamma.size],
            #     yerr=signal_projection_exp_difference_sigma[:gamma.size], fmt="ko", alpha=0.2, label="experiment")
            ax4.plot(ttheta[:gamma.size], signal_projection_difference[:gamma.size], "b-", label="model", linewidth=1)

            y_min_d, y_max_d = ax4.get_ylim()
            param = y_min_d-numpy.nanmax((signal_projection_exp_difference-signal_projection_difference)[:gamma.size])
            ax4.plot([ttheta[:gamma.size].min(), ttheta[:gamma.size].max()], [param, param], "k:")
            ax4.plot(ttheta[:gamma.size], (signal_projection_exp_difference-signal_projection_difference)[:gamma.size]+param,
                "r-", alpha=0.5) #, label="difference"
            ax4.legend()
            ax4.set_xlabel(r"$2\theta$ (deg.)")
            # ax4.get_yaxis().set_visible(False)
        
        return (fig, ax1)

    def calc_projections_sum_difference(self):
        gamma = self.gamma 
        nu = self.nu 

        gamma_rad, nu_rad = gamma*numpy.pi/180., nu*numpy.pi/180.
        ttheta_rad, phi_rad = calc_ttheta_phi_by_gamma_nu(
            gamma_rad[:, na], nu_rad[na, :], flag_gamma=False, flag_nu=False)[:2]
        ttheta_2d = ttheta_rad * 180./numpy.pi

        ttheta_step = numpy.round(gamma[1]-gamma[0], decimals=5)
        ttheta_min = gamma.min()
        ttheta_index = numpy.round((ttheta_2d-ttheta_min)/ttheta_step, decimals=0).astype(int)

        ttheta = numpy.linspace(
            ttheta_min,
            ttheta_step*ttheta_index.max()+ttheta_min,
            ttheta_index.max(), endpoint=False)

        signal_plus = self.intensity_plus_net
        signal_minus = self.intensity_minus_net
        bkg_calc = self.intensity_bkg_calc
        signal_sum = signal_plus + signal_minus + bkg_calc
        signal_difference = signal_plus - signal_minus

        flag_polarized, flag_unpolarized = False, False
        if self.is_attribute("intensity_plus"):
            flag_polarized = True
        elif self.is_attribute("intensity"):
            flag_unpolarized = True
        
        if flag_polarized:
            signal_exp_plus = self.intensity_plus
            signal_exp_minus = self.intensity_minus
            signal_exp_plus_sigma = self.intensity_plus_sigma
            signal_exp_minus_sigma = self.intensity_minus_sigma
            signal_exp_sum = signal_exp_plus + signal_exp_minus
            signal_exp_difference = signal_exp_plus - signal_exp_minus
            signal_exp_sum_sigma_sq = numpy.square(signal_exp_plus_sigma)+numpy.square(signal_exp_minus_sigma)
        elif flag_unpolarized:
            signal_exp_sum = self.intensity
            signal_exp_sum_sigma_sq = numpy.square(self.intensity_sigma)

        excluded_points = numpy.logical_or(
            self.excluded_points,
            numpy.isnan(signal_exp_sum))
        included_points = numpy.logical_not(excluded_points)

        signal_projection_sum = numpy.zeros_like(ttheta)
        signal_projection_exp_sum = numpy.zeros_like(ttheta)
        signal_projection_difference = numpy.zeros_like(ttheta)
        signal_projection_exp_sum_sigma = numpy.zeros_like(ttheta)
        signal_projection_exp_difference = numpy.zeros_like(ttheta)
        signal_projection_exp_difference_sigma = numpy.zeros_like(ttheta)

        if flag_polarized:
            flag_nonan = numpy.logical_not(numpy.isnan(signal_exp_difference))
        elif flag_polarized:
            flag_nonan = numpy.logical_not(numpy.isnan(signal_exp_sum))

        for ind in range(ttheta.size):
            flag_1 = ttheta_index == ind
            flag = numpy.logical_and(flag_1, included_points)
            points = numpy.sum(flag)
            if points != 0:
                signal_projection_sum[ind] += numpy.sum(signal_sum[flag])/points
                signal_projection_exp_sum[ind] += numpy.sum(signal_exp_sum[flag])/points
                signal_projection_exp_sum_sigma[ind] += numpy.sqrt(numpy.sum(signal_exp_sum_sigma_sq[flag]))/points
            else:
                signal_projection_sum[ind] = numpy.nan
                signal_projection_exp_sum[ind] = numpy.nan
                signal_projection_exp_sum_sigma[ind] = numpy.nan
            if flag_polarized:
                flag = numpy.logical_and(flag_1, flag_nonan)
                points = numpy.sum(flag)
                if points != 0:
                    signal_projection_difference[ind] += numpy.sum(signal_difference[flag])/points
                    signal_projection_exp_difference[ind] += numpy.sum(signal_exp_difference[flag])/points
                    signal_projection_exp_difference_sigma[ind] += numpy.sqrt(numpy.sum(signal_exp_sum_sigma_sq[flag]))/points
                else:
                    signal_projection_difference[ind] = numpy.nan
                    signal_projection_exp_difference[ind] = numpy.nan
                    signal_projection_exp_difference_sigma[ind] = numpy.nan
        
        
        return ttheta, signal_projection_sum, signal_projection_exp_sum, signal_projection_exp_sum_sigma, \
            signal_projection_difference, signal_projection_exp_difference, signal_projection_exp_difference_sigma

    def plot_diff_total(self):
        if not(self.is_attribute("gamma") & self.is_attribute("nu") &
               self.is_attribute("intensity_plus_total") & 
               self.is_attribute("intensity_minus_total")):
            return
        fig, ax = plt.subplots()
        ax.set_title("Polarized intensity")
        ax.set_xlabel("2 theta (degrees)")
        ax.set_ylabel('nu (degrees)')

        np_tth = numpy.array(self.gamma, dtype=float)
        np_nu = numpy.array(self.nu, dtype=float)
        np_up = self.intensity_plus_net
        np_down = self.intensity_minus_net
        np_diff = np_up - np_down
        
        norm = matplotlib.colors.Normalize(vmax=abs(np_diff).max(),
                                           vmin=-abs(np_diff).max())
        # cmap = matplotlib.PRGn
        
        extent =(np_tth.min(), np_tth.max(), np_nu.min(), np_nu.max())
        plt.set_cmap('rainbow')
        ax.imshow(np_diff.transpose(), interpolation='bilinear', extent=extent,
                  alpha=0.9, origin="lower", norm=norm)        
        
        # blk = '#000000'
        # ax.contour(np_tth, np_nu, np_sum.transpose(),
        #             #levels=[0.1, 0.5, 1., 5., 10., 50.],
        #             levels=5,
        #             colors=[blk, blk, blk, blk, blk, blk],
        #             linewidths=0.5)
        return (fig, ax)


    def estimate_gamma_nu_offsets(self):
        coeff_difference = 0.5
        gamma = self.gamma
        nu = self.nu
        flag_polarized, flag_unpolarized = False, False
        if self.is_attribute("intensity_plus"):
            flag_polarized = True
        elif self.is_attribute("intensity"):
            flag_unpolarized = True

        signal_sum = self.intensity_plus_net + self.intensity_minus_net + self.intensity_bkg_calc
        signal_difference = self.intensity_plus_net - self.intensity_minus_net
        
        if flag_polarized:
            signal_exp_sum = self.intensity_plus + self.intensity_minus 
            signal_exp_difference = self.intensity_plus - self.intensity_minus
            signal_exp_sum_sigma = numpy.sqrt(
                numpy.square(self.intensity_plus_sigma) +
                numpy.square(self.intensity_minus_sigma))
        elif flag_unpolarized:
            signal_exp_sum = self.intensity
            signal_exp_sum_sigma = self.intensity_sigma

        included_points = numpy.logical_not(self.excluded_points)

        def calc_chi_sq(param):
            offset_gamma = param[0]
            offset_nu = param[1]
            f_sum = interp2d(gamma, nu, signal_sum.transpose(), kind="linear")
            model_sum = f_sum(gamma-offset_gamma, nu-offset_nu).transpose()
            chi_sq_sum = numpy.nansum(numpy.square((model_sum-signal_exp_sum)/signal_exp_sum_sigma)[included_points])
            
            if flag_polarized:
                f_difference = interp2d(gamma, nu, signal_difference.transpose(), kind="linear")
                model_difference = f_difference(gamma-offset_gamma, nu-offset_nu).transpose()
                chi_sq_difference = numpy.nansum(numpy.square((model_difference-signal_exp_difference)/signal_exp_sum_sigma))
                chi_sq = (1. - coeff_difference) * chi_sq_sum + coeff_difference* chi_sq_difference
            elif flag_unpolarized:
                chi_sq = chi_sq_sum
            return chi_sq

        param = numpy.array([0., 0.], dtype=float)
        res = scipy.optimize.minimize(calc_chi_sq, param, method="Nelder-Mead")
        return res
    
    def estimate_background(self,
            pd2d_background: Pd2dBackground, offset_gamma: float = 0, offset_nu: float = 0):
        gamma = self.gamma
        nu = self.nu

        if self.is_attribute("intensity_plus"):
            signal_exp_sum = self.intensity_plus + self.intensity_minus 
            signal_exp_sum_sigma = numpy.sqrt(
                numpy.square(self.intensity_plus_sigma) +
                numpy.square(self.intensity_minus_sigma))
        elif self.is_attribute("intensity"):
            signal_exp_sum = self.intensity 
            signal_exp_sum_sigma = self.intensity_sigma

        signal_sum_net = self.intensity_plus_net + self.intensity_minus_net 
        included_points = numpy.logical_not(self.excluded_points)

        variable_names = pd2d_background.get_variable_names()
        param_0 = [pd2d_background.get_variable_by_name(name) for name in variable_names]
        def temp_func(param):
            for name, value in zip(variable_names, param):
                pd2d_background.set_variable_by_name(name, value)
            model_bkgr = pd2d_background.interpolate_by_points(
                gamma+offset_gamma, nu+offset_nu)
            chi_sq_sum = numpy.nansum(numpy.square((signal_exp_sum - signal_sum_net - model_bkgr)/signal_exp_sum_sigma)[included_points])
            return chi_sq_sum

        res = scipy.optimize.minimize(temp_func, param_0, method="Nelder-Mead")
        pd2d_background.form_gamma_nu_intensity()
        self.intensity_bkg_calc = pd2d_background.interpolate_by_points(gamma+offset_gamma, nu+offset_nu)
        self.form_gamma_nu_intensity_bkg_calc()
        return res