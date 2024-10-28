from typing import NoReturn
import numpy
import matplotlib
import matplotlib.pyplot as plt


from cryspy.A_functions_base.powder_diffraction_const_wavelength import \
    calc_ttheta_phi_by_gamma_nu

from cryspy.A_functions_base.function_1_gamma_nu import \
    recal_int_to_gammanu_grid
    
from cryspy.B_parent_classes.cl_1_item import ItemN

na = numpy.newaxis

class Pd2dMeas(ItemN):
    """
    Measured 2D diffraction pattern (polarized neutrons).
    
    This section contains the measured diffractogram and information
    about the conditions used for the measurement of the diffraction 
    data set, prior to processing and application of correction
    terms.

    """
    ATTR_MANDATORY_NAMES = ()
    ATTR_MANDATORY_TYPES = ()
    # ("matrix", "matrix", "matrix", "matrix")
    ATTR_MANDATORY_CIF = ()

    ATTR_OPTIONAL_NAMES = (
        "ttheta_phi_intensity_plus", "ttheta_phi_intensity_plus_sigma",
        "ttheta_phi_intensity_minus", "ttheta_phi_intensity_minus_sigma",
        "gamma_nu_intensity_plus", "gamma_nu_intensity_plus_sigma",
        "gamma_nu_intensity_minus", "gamma_nu_intensity_minus_sigma",
        "gamma_nu_intensity", "gamma_nu_intensity_sigma")
    ATTR_OPTIONAL_TYPES = (
        str, str, str, str, str, str, str, str, str, str)
    ATTR_OPTIONAL_CIF = (
        "2theta_phi_intensity_plus", "2theta_phi_intensity_plus_sigma",
        "2theta_phi_intensity_minus", "2theta_phi_intensity_minus_sigma",
        "gamma_nu_intensity_plus", "gamma_nu_intensity_plus_sigma",
        "gamma_nu_intensity_minus", "gamma_nu_intensity_minus_sigma",
        "gamma_nu_intensity", "gamma_nu_intensity_sigma")

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ("ttheta", "phi", "intensity_plus", "intensity_plus_sigma",
                      "intensity_minus", "intensity_minus_sigma", "gamma", "nu",
                      "gn_intensity_plus", "gn_intensity_plus_sigma",
                      "gn_intensity_minus", "gn_intensity_minus_sigma",
                      "gn_intensity", "gn_intensity_sigma")
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

    PREFIX = "pd2d_meas"

    def __init__(self, **kwargs) -> NoReturn:
        super(Pd2dMeas, self).__init__()

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
        """Redefine form_object."""
        flag_tth_phi = (self.is_attribute("ttheta_phi_intensity_plus") &
                self.is_attribute("ttheta_phi_intensity_plus_sigma") &
                self.is_attribute("ttheta_phi_intensity_minus") &
                self.is_attribute("ttheta_phi_intensity_minus_sigma"))
        flag_gamma_nu = (
                self.is_attribute("gamma_nu_intensity_plus") &
                self.is_attribute("gamma_nu_intensity_plus_sigma") &
                self.is_attribute("gamma_nu_intensity_minus") &
                self.is_attribute("gamma_nu_intensity_minus_sigma"))
        flag_gamma_nu_unpolarized = (
                self.is_attribute("gamma_nu_intensity") &
                self.is_attribute("gamma_nu_intensity_sigma"))

        if flag_tth_phi:
            l_1 = (self.ttheta_phi_intensity_plus).strip().split("\n")
            l_2 = (self.ttheta_phi_intensity_plus_sigma).strip().split("\n")
            l_3 = (self.ttheta_phi_intensity_minus).strip().split("\n")
            l_4 = (self.ttheta_phi_intensity_minus_sigma).strip().split("\n")

            l_ttheta = numpy.array([_ for _ in l_1[0].strip().split()[1:]],
                                   dtype=float)
            l_phi, ll_intensity_plus, ll_intensity_plus_sigma = [], [], []
            ll_intensity_minus, ll_intensity_minus_sigma = [], []
            for line_1, line_2, line_3, line_4 in zip(l_1[1:], l_2[1:], l_3[1:],
                                                      l_4[1:]):
                _l_1 = line_1.strip().split()
                _l_2 = line_2.strip().split()
                _l_3 = line_3.strip().split()
                _l_4 = line_4.strip().split()
                l_phi.append(float(_l_1[0]))
                ll_intensity_plus.append(_l_1[1:])
                ll_intensity_plus_sigma.append(_l_2[1:])
                ll_intensity_minus.append(_l_3[1:])
                ll_intensity_minus_sigma.append(_l_4[1:])

            ll_intensity_plus = numpy.array(ll_intensity_plus, dtype=float).transpose()
            ll_intensity_plus_sigma = numpy.array(ll_intensity_plus_sigma, dtype=float).transpose()
            ll_intensity_minus = numpy.array(ll_intensity_minus, dtype=float).transpose()
            ll_intensity_minus_sigma = numpy.array(ll_intensity_minus_sigma, dtype=float).transpose()

            self.__dict__["ttheta"] = l_ttheta
            self.__dict__["phi"] = numpy.array(l_phi, dtype=float)
            self.__dict__["intensity_plus"] = ll_intensity_plus
            self.__dict__["intensity_plus_sigma"] = ll_intensity_plus_sigma
            self.__dict__["intensity_minus"] = ll_intensity_minus
            self.__dict__["intensity_minus_sigma"] = ll_intensity_minus_sigma

        if flag_gamma_nu:
            l_1 = (self.gamma_nu_intensity_plus).strip().split("\n")
            l_2 = (self.gamma_nu_intensity_plus_sigma).strip().split("\n")
            l_3 = (self.gamma_nu_intensity_minus).strip().split("\n")
            l_4 = (self.gamma_nu_intensity_minus_sigma).strip().split("\n")

            l_gamma = numpy.array([_ for _ in l_1[0].strip().split()[1:]],
                                   dtype=float)
            l_nu, ll_intensity_plus, ll_intensity_plus_sigma = [], [], []
            ll_intensity_minus, ll_intensity_minus_sigma = [], []
            for line_1, line_2, line_3, line_4 in zip(l_1[1:], l_2[1:], l_3[1:],
                                                      l_4[1:]):
                _l_1 = line_1.strip().split()
                _l_2 = line_2.strip().split()
                _l_3 = line_3.strip().split()
                _l_4 = line_4.strip().split()
                l_nu.append(float(_l_1[0]))
                ll_intensity_plus.append(_l_1[1:])
                ll_intensity_plus_sigma.append(_l_2[1:])
                ll_intensity_minus.append(_l_3[1:])
                ll_intensity_minus_sigma.append(_l_4[1:])

            ll_intensity_plus = numpy.array(ll_intensity_plus, dtype=float).transpose()
            ll_intensity_plus_sigma = numpy.array(ll_intensity_plus_sigma, dtype=float).transpose()
            ll_intensity_minus = numpy.array(ll_intensity_minus, dtype=float).transpose()
            ll_intensity_minus_sigma = numpy.array(ll_intensity_minus_sigma, dtype=float).transpose()

            self.__dict__["gamma"] = l_gamma
            self.__dict__["nu"] = numpy.array(l_nu, dtype=float)
            self.__dict__["gn_intensity_plus"] = ll_intensity_plus
            self.__dict__["gn_intensity_plus_sigma"] = ll_intensity_plus_sigma
            self.__dict__["gn_intensity_minus"] = ll_intensity_minus
            self.__dict__["gn_intensity_minus_sigma"] = ll_intensity_minus_sigma

        if flag_gamma_nu_unpolarized:
            l_1 = (self.gamma_nu_intensity).strip().split("\n")
            l_2 = (self.gamma_nu_intensity_sigma).strip().split("\n")

            l_gamma = numpy.array([_ for _ in l_1[0].strip().split()[1:]],
                                   dtype=float)
            l_nu, ll_intensity, ll_intensity_sigma = [], [], []
            for line_1, line_2 in zip(l_1[1:], l_2[1:]):
                _l_1 = line_1.strip().split()
                _l_2 = line_2.strip().split()
                l_nu.append(float(_l_1[0]))
                ll_intensity.append(_l_1[1:])
                ll_intensity_sigma.append(_l_2[1:])

            ll_intensity = numpy.array(ll_intensity, dtype=float).transpose()
            ll_intensity_sigma = numpy.array(ll_intensity_sigma, dtype=float).transpose()

            self.__dict__["gamma"] = l_gamma
            self.__dict__["nu"] = numpy.array(l_nu, dtype=float)
            self.__dict__["gn_intensity"] = ll_intensity
            self.__dict__["gn_intensity_sigma"] = ll_intensity_sigma

    def form_ttheta_phi_intensity_plus(self) -> bool:
        if (self.is_attribute("phi") & self.is_attribute("ttheta") &
            self.is_attribute("intensity_plus")):
            ls_out = []
            ls_out.append("{:12} ".format(len(self.phi)) + " ".join(["{:6.2f}      ".format(_) for _ in self.ttheta]))
            ll_intensity = self.intensity_plus
            ll_intensity = [[ll_intensity[_2][_1] for _2 in range(len(ll_intensity))] for _1 in range(len(ll_intensity[0]))]
            for phi, l_intensity in zip(self.phi, ll_intensity):
                ls_out.append("{:12.2f} ".format(phi) + " ".join(["{:12}".format(_) for _ in l_intensity]))
            self.__dict__["ttheta_phi_intensity_plus"] = "\n".join(ls_out)

    def form_ttheta_phi_intensity_plus_sigma(self) -> bool:
        if (self.is_attribute("phi") & self.is_attribute("ttheta") &
            self.is_attribute("intensity_plus_sigma")):
            ls_out = []
            ls_out.append("{:12} ".format(len(self.phi)) + " ".join(["{:6.2f}      ".format(_) for _ in self.ttheta]))
            ll_intensity = self.intensity_plus_sigma
            ll_intensity = [[ll_intensity[_2][_1] for _2 in range(len(ll_intensity))] for _1 in range(len(ll_intensity[0]))]
            for phi, l_intensity in zip(self.phi, ll_intensity):
                ls_out.append("{:12.2f} ".format(phi) + " ".join(["{:12}".format(_) for _ in l_intensity]))
            self.__dict__["ttheta_phi_intensity_plus_sigma"] = "\n".join(ls_out)

    def form_ttheta_phi_intensity_minus(self) -> bool:
        if (self.is_attribute("phi") & self.is_attribute("ttheta") &
            self.is_attribute("intensity_minus")):
            ls_out = []
            ls_out.append("{:12} ".format(len(self.phi)) + " ".join(["{:6.2f}      ".format(_) for _ in self.ttheta]))
            ll_intensity = self.intensity_minus
            ll_intensity = [[ll_intensity[_2][_1] for _2 in range(len(ll_intensity))] for _1 in range(len(ll_intensity[0]))]
            for phi, l_intensity in zip(self.phi, ll_intensity):
                ls_out.append("{:12.2f} ".format(phi) + " ".join(["{:12}".format(_) for _ in l_intensity]))
            self.__dict__["ttheta_phi_intensity_minus"] = "\n".join(ls_out)

    def form_ttheta_phi_intensity_minus_sigma(self) -> bool:
        if (self.is_attribute("phi") & self.is_attribute("ttheta") &
            self.is_attribute("intensity_minus_sigma")):
            ls_out = []
            ls_out.append("{:12} ".format(len(self.phi)) + " ".join(["{:6.2f}      ".format(_) for _ in self.ttheta]))
            ll_intensity = self.intensity_minus_sigma
            ll_intensity = [[ll_intensity[_2][_1] for _2 in range(len(ll_intensity))] for _1 in range(len(ll_intensity[0]))]
            for phi, l_intensity in zip(self.phi, ll_intensity):
                ls_out.append("{:12.2f} ".format(phi) + " ".join(["{:12}".format(_) for _ in l_intensity]))
            self.__dict__["ttheta_phi_intensity_minus_sigma"] = "\n".join(ls_out)

    def form_gamma_nu_intensity_plus(self) -> bool:
        if (self.is_attribute("nu") & self.is_attribute("gamma") &
            self.is_attribute("gn_intensity_plus")):
            ls_out = []
            ls_out.append("{:12} ".format(len(self.nu)) + " ".join(["{:6.2f}      ".format(_) for _ in self.gamma]))
            ll_intensity = self.gn_intensity_plus
            ll_intensity = [[ll_intensity[_2][_1] for _2 in range(len(ll_intensity))] for _1 in range(len(ll_intensity[0]))]
            for nu, l_intensity in zip(self.nu, ll_intensity):
                ls_out.append("{:12.2f} ".format(nu) + " ".join(["{:12}".format(_) for _ in l_intensity]))
            self.__dict__["gamma_nu_intensity_plus"] = "\n".join(ls_out)

    def form_gamma_nu_intensity_plus_sigma(self) -> bool:
        if (self.is_attribute("nu") & self.is_attribute("gamma") &
            self.is_attribute("gn_intensity_plus_sigma")):
            ls_out = []
            ls_out.append("{:12} ".format(len(self.nu)) + " ".join(["{:6.2f}      ".format(_) for _ in self.gamma]))
            ll_intensity = self.gn_intensity_plus_sigma
            ll_intensity = [[ll_intensity[_2][_1] for _2 in range(len(ll_intensity))] for _1 in range(len(ll_intensity[0]))]
            for nu, l_intensity in zip(self.nu, ll_intensity):
                ls_out.append("{:12.2f} ".format(nu) + " ".join(["{:12}".format(_) for _ in l_intensity]))
            self.__dict__["gamma_nu_intensity_plus_sigma"] = "\n".join(ls_out)

    def form_gamma_nu_intensity_minus(self) -> bool:
        if (self.is_attribute("nu") & self.is_attribute("gamma") &
            self.is_attribute("gn_intensity_minus")):
            ls_out = []
            ls_out.append("{:12} ".format(len(self.nu)) + " ".join(["{:6.2f}      ".format(_) for _ in self.gamma]))
            ll_intensity = self.gn_intensity_minus
            ll_intensity = [[ll_intensity[_2][_1] for _2 in range(len(ll_intensity))] for _1 in range(len(ll_intensity[0]))]
            for nu, l_intensity in zip(self.nu, ll_intensity):
                ls_out.append("{:12.2f} ".format(nu) + " ".join(["{:12}".format(_) for _ in l_intensity]))
            self.__dict__["gamma_nu_intensity_minus"] = "\n".join(ls_out)

    def form_gamma_nu_intensity_minus_sigma(self) -> bool:
        if (self.is_attribute("nu") & self.is_attribute("gamma") &
            self.is_attribute("gn_intensity_minus_sigma")):
            ls_out = []
            ls_out.append("{:12} ".format(len(self.nu)) + " ".join(["{:6.2f}      ".format(_) for _ in self.gamma]))
            ll_intensity = self.gn_intensity_minus_sigma
            ll_intensity = [[ll_intensity[_2][_1] for _2 in range(len(ll_intensity))] for _1 in range(len(ll_intensity[0]))]
            for nu, l_intensity in zip(self.nu, ll_intensity):
                ls_out.append("{:12.2f} ".format(nu) + " ".join(["{:12}".format(_) for _ in l_intensity]))
            self.__dict__["gamma_nu_intensity_minus_sigma"] = "\n".join(ls_out)


    def form_gamma_nu_intensity(self) -> bool:
        if (self.is_attribute("nu") & self.is_attribute("gamma") &
            self.is_attribute("gn_intensity")):
            ls_out = []
            ls_out.append("{:12} ".format(len(self.nu)) + " ".join(["{:6.2f}      ".format(_) for _ in self.gamma]))
            ll_intensity = self.gn_intensity
            ll_intensity = [[ll_intensity[_2][_1] for _2 in range(len(ll_intensity))] for _1 in range(len(ll_intensity[0]))]
            for nu, l_intensity in zip(self.nu, ll_intensity):
                ls_out.append("{:12.2f} ".format(nu) + " ".join(["{:12}".format(_) for _ in l_intensity]))
            self.__dict__["gamma_nu_intensity"] = "\n".join(ls_out)

    def form_gamma_nu_intensity_sigma(self) -> bool:
        if (self.is_attribute("nu") & self.is_attribute("gamma") &
            self.is_attribute("gn_intensity_sigma")):
            ls_out = []
            ls_out.append("{:12} ".format(len(self.nu)) + " ".join(["{:6.2f}      ".format(_) for _ in self.gamma]))
            ll_intensity = self.gn_intensity_sigma
            ll_intensity = [[ll_intensity[_2][_1] for _2 in range(len(ll_intensity))] for _1 in range(len(ll_intensity[0]))]
            for nu, l_intensity in zip(self.nu, ll_intensity):
                ls_out.append("{:12.2f} ".format(nu) + " ".join(["{:12}".format(_) for _ in l_intensity]))
            self.__dict__["gamma_nu_intensity_sigma"] = "\n".join(ls_out)


    def recalc_to_gamma_nu_grid(self):
        l_tth_grid = numpy.array(self.ttheta)*numpy.pi/180.
        l_phi_grid = numpy.array(self.phi)*numpy.pi/180.
        int_u = numpy.array(self.intensity_plus, dtype=float).transpose()
        int_d = numpy.array(self.intensity_minus, dtype=float).transpose()
        int_sum = int_u + int_d
        int_diff = int_u - int_d

        min_tth, max_tth = min(l_tth_grid), max(l_tth_grid)
        min_phi, max_phi = min(l_phi_grid), max(l_phi_grid)

        min_gamma, max_gamma = min_tth, max_tth
        num_gamma = len(l_tth_grid)

        min_nu, max_nu = -10.*numpy.pi/180., 15.*numpy.pi/180.
        num_nu = len(l_phi_grid)

        l_gamma_grid = numpy.linspace(min_gamma, max_gamma, num=num_gamma)
        l_nu_grid = numpy.linspace(min_nu, max_nu, num=num_nu)

        int_u_out = numpy.array(recal_int_to_gammanu_grid(l_tth_grid, l_phi_grid, int_u, l_gamma_grid, l_nu_grid), dtype=float)
        int_d_out = numpy.array(recal_int_to_gammanu_grid(l_tth_grid, l_phi_grid, int_d, l_gamma_grid, l_nu_grid), dtype=float)
        int_sum_out = numpy.array(recal_int_to_gammanu_grid(l_tth_grid, l_phi_grid, int_sum, l_gamma_grid, l_nu_grid), dtype=float)
        int_diff_out = numpy.array(recal_int_to_gammanu_grid(l_tth_grid, l_phi_grid, int_diff, l_gamma_grid, l_nu_grid), dtype=float)

        return l_gamma_grid*180./numpy.pi, l_nu_grid*180./numpy.pi, [int_u_out, int_d_out, int_sum_out, int_diff_out]
        
    def plots(self):
        return [self.plot_gamma_nu(), ]

    def plot_gamma_nu(self):
        flag_polarized = (self.is_attribute("gamma") & self.is_attribute("nu") &
               self.is_attribute("gn_intensity_plus") & 
               self.is_attribute("gn_intensity_minus"))
        flag_unpolarized = (self.is_attribute("gamma") & self.is_attribute("nu") &
               self.is_attribute("gn_intensity"))

        if not(self.is_attribute("gamma") & self.is_attribute("nu")):
            return

        gamma = self.gamma 
        nu = self.nu 
        if flag_polarized:
            signal_exp_plus = self.gn_intensity_plus
            signal_exp_minus = self.gn_intensity_minus
            signal_exp_plus_sigma = self.gn_intensity_plus_sigma
            signal_exp_minus_sigma = self.gn_intensity_minus_sigma
            signal_exp_sum = signal_exp_plus + signal_exp_minus
            signal_exp_difference = signal_exp_plus - signal_exp_minus
            signal_exp_sum_sigma_sq = numpy.square(signal_exp_plus_sigma)+numpy.square(signal_exp_minus_sigma)
        elif flag_unpolarized:
            signal_exp_sum = self.gn_intensity
            signal_exp_sum_sigma_sq = numpy.square(self.gn_intensity_sigma)
        else:
            return

        max_val = numpy.nanmax(signal_exp_sum)
        min_val = numpy.nanmin(signal_exp_sum)
        max_val = max_val - (max_val-min_val)*0.25

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

        signal_projection_exp_sum = numpy.zeros_like(ttheta)
        signal_projection_exp_difference = numpy.zeros_like(ttheta)
        signal_projection_exp_sum_sigma = numpy.zeros_like(ttheta)
        signal_projection_exp_difference_sigma = numpy.zeros_like(ttheta)

        if flag_polarized:
            flag_nonan = numpy.logical_not(numpy.isnan(signal_exp_difference))
        elif flag_unpolarized:
            flag_nonan = numpy.logical_not(numpy.isnan(signal_exp_sum))
        for ind in range(ttheta.size):
            flag_1 = ttheta_index == ind
            flag = numpy.logical_and(flag_1, flag_nonan)
            points = numpy.sum(flag)
            if points != 0:
                signal_projection_exp_sum[ind] += numpy.sum(signal_exp_sum[flag])/points
                signal_projection_exp_sum_sigma[ind] += numpy.sqrt(numpy.sum(signal_exp_sum_sigma_sq[flag]))/points
            else:
                signal_projection_exp_sum[ind] = numpy.nan
                signal_projection_exp_sum_sigma[ind] = numpy.nan

            if flag_polarized:
                flag = numpy.logical_and(flag_1, flag_nonan)
                points = numpy.sum(flag)
                if points != 0:
                    signal_projection_exp_difference[ind] += numpy.sum(signal_exp_difference[flag])/points
                    signal_projection_exp_difference_sigma[ind] += numpy.sqrt(numpy.sum(signal_exp_sum_sigma_sq[flag]))/points
                else:
                    signal_projection_exp_difference[ind] = numpy.nan
                    signal_projection_exp_difference_sigma[ind] = numpy.nan


        extent = [gamma.min(), gamma.max(), nu.min(), nu.max()]

        cmap_sum = plt.get_cmap("turbo") # BuPu
        if flag_polarized:
            fig, axs = plt.subplots(2, 2, sharex=True)
            ax1, ax2, ax3, ax4 = axs[0,0], axs[0,1], axs[1,0], axs[1,1]
        elif flag_unpolarized:
            fig, axs = plt.subplots(2, 1, sharex=True)
            ax1, ax3= axs[0], axs[1]

        norm_1 = matplotlib.colors.Normalize(vmax=max_val, vmin=min_val)
        hh = numpy.copy(signal_exp_sum)
        hh[numpy.isnan(hh)]=0
        ax1.imshow(hh.transpose(), origin="lower", aspect="auto", cmap=cmap_sum, norm= norm_1, extent=extent)

        ax1.set_xticks([])
        ax1.set_xlabel(r"$\gamma$ (deg.)")
        ax1.set_yticks([])
        ax1.set_ylabel(r"$\nu$ (deg.)")

        if flag_polarized:
            cmap_difference = plt.get_cmap("turbo") # BrBG
            max_val = numpy.nanmax(numpy.abs(signal_exp_difference))*0.50
            norm_2 = matplotlib.colors.Normalize(vmax=max_val, vmin=-max_val)
            hh = numpy.copy(signal_exp_difference)
            hh[numpy.isnan(hh)]=0
            ax2.imshow(hh.transpose(), origin="lower", aspect="auto",
                      norm=norm_2,
                      cmap=cmap_difference,
                      extent=extent)
            # ax2.imshow(zz_difference, origin="lower", aspect="auto", alpha=0.8, extent=extent)
            # ax2.contour(zz_difference[:,:,3], levels=[50, 70])

            # ax2.get_xaxis().set_visible(False)
            # ax2.get_yaxis().set_visible(False)
            ax2.set_xticks([])
            ax2.set_xlabel(r"$\gamma$ (deg.)")
            ax2.set_yticks([])
            ax2.set_ylabel(r"$\nu$ (deg.)")

        ax1.set_title(r"Unpolarized signal")


        ax3.fill_between(
            ttheta[:gamma.size],(signal_projection_exp_sum-signal_projection_exp_sum_sigma)[:gamma.size],
            (signal_projection_exp_sum+signal_projection_exp_sum_sigma)[:gamma.size], color="k", alpha=0.4)

        # ax3.errorbar(
        #     ttheta[:gamma.size], signal_projection_exp_sum[:gamma.size],
        #     yerr=signal_projection_exp_sum_sigma[:gamma.size], fmt="ko", alpha=0.2)
        # ax3.legend()
        ax3.set_xlabel(r"$2\theta$ (deg.)")

        if flag_polarized:
            ax2.set_title(r"Polarized signal")

            ax4.fill_between(
                ttheta[:gamma.size],(signal_projection_exp_difference-signal_projection_exp_difference_sigma)[:gamma.size],
                (signal_projection_exp_difference+signal_projection_exp_difference_sigma)[:gamma.size], color="k", alpha=0.4)

            # ax4.errorbar(
            #     ttheta[:gamma.size], signal_projection_exp_difference[:gamma.size],
            #     yerr=signal_projection_exp_difference_sigma[:gamma.size], fmt="ko", alpha=0.2)

            # ax4.legend()
            ax4.set_xlabel(r"$2\theta$ (deg.)")
        
        return (fig, ax1)

