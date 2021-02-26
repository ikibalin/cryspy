from typing import NoReturn
import numpy
import matplotlib
import matplotlib.pyplot as plt

from cryspy.A_functions_base.function_1_gamma_nu import \
    recal_int_to_gammanu_grid

from cryspy.B_parent_classes.cl_1_item import ItemN

class Pd2dProc(ItemN):
    """
    Calculated 2D diffraction pattern (polarized neutrons).

    Attributes
    ----------
        - ttheta_phi_intensity_up_net, ttheta_phi_intensity_down_net,
          ttheta_phi_intensity_up_total, ttheta_phi_intensity_down_total,
          ttheta_phi_intensity_bkg_calc, ttheta_phi_intensity_up,
          ttheta_phi_intensity_up_sigma, ttheta_phi_intensity_down,
          ttheta_phi_intensity_down_sigma (mandatory)

    Internal
    --------
        - ttheta, phi, intensity_up_net, intensity_down_net,
          intensity_up_total, intensity_down_total, intensity_bkg_calc,
          intensity_up, intensity_up_sigma, intensity_down,
          intensity_down_sigma

    """
    ATTR_MANDATORY_NAMES = (
        "ttheta_phi_intensity_up_net", "ttheta_phi_intensity_down_net",
        "ttheta_phi_intensity_up_total", "ttheta_phi_intensity_down_total",
        "ttheta_phi_intensity_bkg_calc", "ttheta_phi_intensity_up",
        "ttheta_phi_intensity_up_sigma", "ttheta_phi_intensity_down",
        "ttheta_phi_intensity_down_sigma")
    ATTR_MANDATORY_TYPES = (str, str, str, str, str, str, str, str, str)
    # ("matrix", "matrix", "matrix", "matrix")
    ATTR_MANDATORY_CIF = (
        "2theta_phi_intensity_up_net", "2theta_phi_intensity_down_net",
        "2theta_phi_intensity_up_total", "2theta_phi_intensity_down_total",
        "2theta_phi_intensity_bkg_calc", "2theta_phi_intensity_up",
        "2theta_phi_intensity_up_sigma", "2theta_phi_intensity_down",
        "2theta_phi_intensity_down_sigma")

    ATTR_OPTIONAL_NAMES = ()
    ATTR_OPTIONAL_TYPES = ()
    ATTR_OPTIONAL_CIF = ()

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = (
        "ttheta", "phi", "intensity_up_net", "intensity_down_net",
        "intensity_up_total", "intensity_down_total", "intensity_bkg_calc",
        "intensity_up", "intensity_up_sigma", "intensity_down",
        "intensity_down_sigma")
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
        if any([self.ttheta_phi_intensity_up_net is None,
                self.ttheta_phi_intensity_down_net is None,
                self.ttheta_phi_intensity_up_total is None,
                self.ttheta_phi_intensity_down_total is None,
                self.ttheta_phi_intensity_bkg_calc is None,
                self.ttheta_phi_intensity_up is None,
                self.ttheta_phi_intensity_up_sigma is None,
                self.ttheta_phi_intensity_down is None,
                self.ttheta_phi_intensity_down_sigma is None]):
            return False
        l_1 = (self.ttheta_phi_intensity_up_net).strip().split("\n")
        l_2 = (self.ttheta_phi_intensity_down_net).strip().split("\n")
        l_3 = (self.ttheta_phi_intensity_up_total).strip().split("\n")
        l_4 = (self.ttheta_phi_intensity_down_total).strip().split("\n")
        l_5 = (self.ttheta_phi_intensity_bkg_calc).strip().split("\n")
        l_6 = (self.ttheta_phi_intensity_up).strip().split("\n")
        l_7 = (self.ttheta_phi_intensity_up_sigma).strip().split("\n")
        l_8 = (self.ttheta_phi_intensity_down).strip().split("\n")
        l_9 = (self.ttheta_phi_intensity_down_sigma).strip().split("\n")

        l_ttheta = numpy.array([_ for _ in l_1[0].strip().split()[1:]], dtype=float)
        l_phi, ll_intensity_up, ll_intensity_up_sigma = [], [], []
        ll_intensity_down, ll_intensity_down_sigma = [], []
        ll_intensity_up_net, ll_intensity_down_net = [], []
        ll_intensity_up_total, ll_intensity_down_total = [], []
        ll_intensity_bkg_calc = []
        for line_1, line_2, line_3, line_4, line_5, line_6, line_7, line_8, line_9 in zip(
            l_1[1:], l_2[1:], l_3[1:], l_4[1:], l_5[1:], l_6[1:], l_7[1:], l_8[1:], l_9[1:]):
            _l_1 = line_1.strip().split()
            _l_2 = line_2.strip().split()
            _l_3 = line_3.strip().split()
            _l_4 = line_4.strip().split()
            _l_5 = line_5.strip().split()
            _l_6 = line_6.strip().split()
            _l_7 = line_7.strip().split()
            _l_8 = line_8.strip().split()
            _l_9 = line_9.strip().split()
            l_phi.append(float(_l_1[0]))
            ll_intensity_up_net.append(_l_1[1:])
            ll_intensity_down_net.append(_l_2[1:])
            ll_intensity_up_total.append(_l_3[1:])
            ll_intensity_down_total.append(_l_4[1:])
            ll_intensity_bkg_calc.append(_l_5[1:])
            ll_intensity_up.append(_l_6[1:])
            ll_intensity_up_sigma.append(_l_7[1:])
            ll_intensity_down.append(_l_8[1:])
            ll_intensity_down_sigma.append(_l_9[1:])

        ll_intensity_up_net = numpy.array(ll_intensity_up_net, dtype=float).transpose()
        ll_intensity_down_net = numpy.array(ll_intensity_down_net, dtype=float).transpose()
        ll_intensity_up_total = numpy.array(ll_intensity_up_total, dtype=float).transpose()
        ll_intensity_down_total = numpy.array(ll_intensity_down_total, dtype=float).transpose()
        ll_intensity_bkg_calc = numpy.array(ll_intensity_bkg_calc, dtype=float).transpose()
        ll_intensity_up = numpy.array(ll_intensity_up, dtype=float).transpose()
        ll_intensity_up_sigma = numpy.array(ll_intensity_up_sigma, dtype=float).transpose()
        ll_intensity_down = numpy.array(ll_intensity_down, dtype=float).transpose()
        ll_intensity_down_sigma = numpy.array(ll_intensity_down_sigma, dtype=float).transpose()

        self.__dict__["ttheta"] = l_ttheta
        self.__dict__["phi"] = numpy.array(l_phi, dtype=float)
        self.__dict__["intensity_up_net"] = ll_intensity_up_net
        self.__dict__["intensity_down_net"] = ll_intensity_down_net
        self.__dict__["intensity_up_total"] = ll_intensity_up_total
        self.__dict__["intensity_down_total"] = ll_intensity_down_total
        self.__dict__["intensity_bkg_calc"] = ll_intensity_bkg_calc
        self.__dict__["intensity_up"] = ll_intensity_up
        self.__dict__["intensity_up_sigma"] = ll_intensity_up_sigma
        self.__dict__["intensity_down"] = ll_intensity_down
        self.__dict__["intensity_down_sigma"] = ll_intensity_down_sigma
        return flag

    def form_ttheta_phi_intensity_up_net(self) -> bool:
        if ((self.phi is not None) & (self.ttheta is not None) & (self.intensity_up_net is not None)):
            ls_out = []
            ls_out.append("{:12} ".format(len(self.phi)) + " ".join(["{:6.2f}      ".format(_) for _ in self.ttheta]))
            ll_intensity = self.intensity_up_net
            ll_intensity = [[ll_intensity[_2][_1] for _2 in range(len(ll_intensity))] for _1 in range(len(ll_intensity[0]))]
            for phi, l_intensity in zip(self.phi, ll_intensity):
                ls_out.append("{:12.2f} ".format(phi) + " ".join(["{:12}".format(_) for _ in l_intensity]))
            self.__dict__["ttheta_phi_intensity_up_net"] = "\n".join(ls_out)

    def form_ttheta_phi_intensity_down_net(self) -> bool:
        if ((self.phi is not None) & (self.ttheta is not None) & (self.intensity_down_net is not None)):
            ls_out = []
            ls_out.append("{:12} ".format(len(self.phi)) + " ".join(["{:6.2f}      ".format(_) for _ in self.ttheta]))
            ll_intensity = self.intensity_down_net
            ll_intensity = [[ll_intensity[_2][_1] for _2 in range(len(ll_intensity))] for _1 in range(len(ll_intensity[0]))]
            for phi, l_intensity in zip(self.phi, ll_intensity):
                ls_out.append("{:12.2f} ".format(phi) + " ".join(["{:12}".format(_) for _ in l_intensity]))
            self.__dict__["ttheta_phi_intensity_down_net"] = "\n".join(ls_out)

    def form_ttheta_phi_intensity_up_total(self) -> bool:
        if ((self.phi is not None) & (self.ttheta is not None) & (self.intensity_up_total is not None)):
            ls_out = []
            ls_out.append("{:12} ".format(len(self.phi)) + " ".join(["{:6.2f}      ".format(_) for _ in self.ttheta]))
            ll_intensity = self.intensity_up_total
            ll_intensity = [[ll_intensity[_2][_1] for _2 in range(len(ll_intensity))] for _1 in range(len(ll_intensity[0]))]
            for phi, l_intensity in zip(self.phi, ll_intensity):
                ls_out.append("{:12.2f} ".format(phi) + " ".join(["{:12}".format(_) for _ in l_intensity]))
            self.__dict__["ttheta_phi_intensity_up_total"] = "\n".join(ls_out)

    def form_ttheta_phi_intensity_down_total(self) -> bool:
        if ((self.phi is not None) & (self.ttheta is not None) & (self.intensity_down_total is not None)):
            ls_out = []
            ls_out.append("{:12} ".format(len(self.phi)) + " ".join(["{:6.2f}      ".format(_) for _ in self.ttheta]))
            ll_intensity = self.intensity_down_total
            ll_intensity = [[ll_intensity[_2][_1] for _2 in range(len(ll_intensity))] for _1 in range(len(ll_intensity[0]))]
            for phi, l_intensity in zip(self.phi, ll_intensity):
                ls_out.append("{:12.2f} ".format(phi) + " ".join(["{:12}".format(_) for _ in l_intensity]))
            self.__dict__["ttheta_phi_intensity_down_total"] = "\n".join(ls_out)

    def form_ttheta_phi_intensity_bkg_calc(self) -> bool:
        if ((self.phi is not None) & (self.ttheta is not None) & (self.intensity_bkg_calc is not None)):
            ls_out = []
            ls_out.append("{:12} ".format(len(self.phi)) + " ".join(["{:6.2f}      ".format(_) for _ in self.ttheta]))
            ll_intensity = self.intensity_bkg_calc
            ll_intensity = [[ll_intensity[_2][_1] for _2 in range(len(ll_intensity))] for _1 in range(len(ll_intensity[0]))]
            for phi, l_intensity in zip(self.phi, ll_intensity):
                ls_out.append("{:12.2f} ".format(phi) + " ".join(["{:12}".format(_) for _ in l_intensity]))
            self.__dict__["ttheta_phi_intensity_bkg_calc"] = "\n".join(ls_out)


    def form_ttheta_phi_intensity_up(self) -> bool:
        if ((self.phi is not None) & (self.ttheta is not None) & (self.intensity_up is not None)):
            ls_out = []
            ls_out.append("{:12} ".format(len(self.phi)) + " ".join(["{:6.2f}      ".format(_) for _ in self.ttheta]))
            ll_intensity = self.intensity_up
            ll_intensity = [[ll_intensity[_2][_1] for _2 in range(len(ll_intensity))] for _1 in range(len(ll_intensity[0]))]
            for phi, l_intensity in zip(self.phi, ll_intensity):
                ls_out.append("{:12.2f} ".format(phi) + " ".join(["{:12}".format(_) for _ in l_intensity]))
            self.__dict__["ttheta_phi_intensity_up"] = "\n".join(ls_out)

    def form_ttheta_phi_intensity_up_sigma(self) -> bool:
        if ((self.phi is not None) & (self.ttheta is not None) & (self.intensity_up_sigma is not None)):
            ls_out = []
            ls_out.append("{:12} ".format(len(self.phi)) + " ".join(["{:6.2f}      ".format(_) for _ in self.ttheta]))
            ll_intensity = self.intensity_up_sigma
            ll_intensity = [[ll_intensity[_2][_1] for _2 in range(len(ll_intensity))] for _1 in range(len(ll_intensity[0]))]
            for phi, l_intensity in zip(self.phi, ll_intensity):
                ls_out.append("{:12.2f} ".format(phi) + " ".join(["{:12}".format(_) for _ in l_intensity]))
            self.__dict__["ttheta_phi_intensity_up_sigma"] = "\n".join(ls_out)

    def form_ttheta_phi_intensity_down(self) -> bool:
        if ((self.phi is not None) & (self.ttheta is not None) & (self.intensity_down is not None)):
            ls_out = []
            ls_out.append("{:12} ".format(len(self.phi)) + " ".join(["{:6.2f}      ".format(_) for _ in self.ttheta]))
            ll_intensity = self.intensity_down
            ll_intensity = [[ll_intensity[_2][_1] for _2 in range(len(ll_intensity))] for _1 in range(len(ll_intensity[0]))]
            for phi, l_intensity in zip(self.phi, ll_intensity):
                ls_out.append("{:12.2f} ".format(phi) + " ".join(["{:12}".format(_) for _ in l_intensity]))
            self.__dict__["ttheta_phi_intensity_down"] = "\n".join(ls_out)

    def form_ttheta_phi_intensity_down_sigma(self) -> bool:
        if ((self.phi is not None) & (self.ttheta is not None) & (self.intensity_down_sigma is not None)):
            ls_out = []
            ls_out.append("{:12} ".format(len(self.phi)) + " ".join(["{:6.2f}      ".format(_) for _ in self.ttheta]))
            ll_intensity = self.intensity_down_sigma
            ll_intensity = [[ll_intensity[_2][_1] for _2 in range(len(ll_intensity))] for _1 in range(len(ll_intensity[0]))]
            for phi, l_intensity in zip(self.phi, ll_intensity):
                ls_out.append("{:12.2f} ".format(phi) + " ".join(["{:12}".format(_) for _ in l_intensity]))
            self.__dict__["ttheta_phi_intensity_down_sigma"] = "\n".join(ls_out)

    def recalc_to_gamma_nu_grid(self):
        l_tth_grid = numpy.array(self.ttheta)*numpy.pi/180.
        l_phi_grid = numpy.array(self.phi)*numpy.pi/180.
        int_u = numpy.array(self.intensity_up, dtype=float).transpose()
        int_d = numpy.array(self.intensity_down, dtype=float).transpose()
        int_sum = int_u + int_d
        int_diff = int_u - int_d

        int_u_m = numpy.array(self.intensity_up_total, dtype=float).transpose()
        int_d_m = numpy.array(self.intensity_down_total, dtype=float).transpose()
        int_sum_m = int_u_m + int_d_m
        int_diff_m = int_u_m - int_d_m

        min_tth, max_tth = min(l_tth_grid), max(l_tth_grid)
        min_phi, max_phi = min(l_phi_grid), max(l_phi_grid)
        
        min_gamma, max_gamma = min_tth, max_tth
        num_gamma = len(l_tth_grid)

        min_nu, max_nu = -3.*numpy.pi/180., 15.*numpy.pi/180.
        num_nu = len(l_phi_grid)

        l_gamma_grid = numpy.linspace(min_gamma, max_gamma, num=num_gamma)
        l_nu_grid = numpy.linspace(min_nu, max_nu, num=num_nu)

        int_u_out = numpy.array(recal_int_to_gammanu_grid(l_tth_grid, l_phi_grid, int_u, l_gamma_grid, l_nu_grid), dtype=float)
        int_d_out = numpy.array(recal_int_to_gammanu_grid(l_tth_grid, l_phi_grid, int_d, l_gamma_grid, l_nu_grid), dtype=float)
        int_sum_out = numpy.array(recal_int_to_gammanu_grid(l_tth_grid, l_phi_grid, int_sum, l_gamma_grid, l_nu_grid), dtype=float)
        int_diff_out = numpy.array(recal_int_to_gammanu_grid(l_tth_grid, l_phi_grid, int_diff, l_gamma_grid, l_nu_grid), dtype=float)

        int_u_m_out = numpy.array(recal_int_to_gammanu_grid(l_tth_grid, l_phi_grid, int_u_m, l_gamma_grid, l_nu_grid), dtype=float)
        int_d_m_out = numpy.array(recal_int_to_gammanu_grid(l_tth_grid, l_phi_grid, int_d_m, l_gamma_grid, l_nu_grid), dtype=float)
        int_sum_m_out = numpy.array(recal_int_to_gammanu_grid(l_tth_grid, l_phi_grid, int_sum_m, l_gamma_grid, l_nu_grid), dtype=float)
        int_diff_m_out = numpy.array(recal_int_to_gammanu_grid(l_tth_grid, l_phi_grid, int_diff_m, l_gamma_grid, l_nu_grid), dtype=float)

        return l_gamma_grid*180./numpy.pi, l_nu_grid*180./numpy.pi, [int_u_out, int_d_out, int_sum_out, int_diff_out, int_u_m_out, int_d_m_out, int_sum_m_out, int_diff_m_out]

        
    def plots(self):
        return [self.plot_projection_sum(), self.plot_projection_diff()]
    
    def plot_projection_sum(self):
        """Plot experimental unpolarized intensity vs. 2 theta (degrees)
        """
        if not(self.is_attribute("ttheta") & self.is_attribute("intensity_up") & 
               self.is_attribute("intensity_up_sigma") &
               self.is_attribute("intensity_down") & 
               self.is_attribute("intensity_down_sigma") &
               self.is_attribute("intensity_up_total") &
               self.is_attribute("intensity_down_total")):
            return
        fig, ax = plt.subplots()
        ax.set_title("Projection of unpolarized intensity on 2 theta axis")
        ax.set_xlabel("2 theta (degrees)")
        ax.set_ylabel('Intensity')

        np_tth = numpy.array(self.ttheta, dtype=float)
        np_up = self.intensity_up
        np_sup = self.intensity_up_sigma
        np_down = self.intensity_down
        np_sdown = self.intensity_down_sigma
        np_sum = np_up + np_down
        np_ssum = numpy.sqrt(numpy.square(np_sup)+numpy.square(np_sdown))

        np_up_m = self.intensity_up_total
        np_down_m = self.intensity_down_total
        np_sum_m = np_up_m + np_down_m

        np_sum_1d = numpy.where(numpy.isnan(np_sum), 0., np_sum).sum(axis=1)
        np_ssum_1d = ((numpy.where(numpy.isnan(np_sum), 0., np_ssum)**2).sum(
            axis=1))**0.5

        np_sum_m_1d = numpy.where(numpy.isnan(np_sum), 0., np_sum_m).sum(
            axis=1)

        ax.plot(np_tth, np_sum_m_1d, "k-", label="model")
        ax.errorbar(np_tth, np_sum_1d, yerr=np_ssum_1d, fmt="ko", alpha=0.2,
                    label="experiment")

        y_min_d, y_max_d = ax.get_ylim()
        param = y_min_d-(np_sum_1d-np_sum_m_1d).max()

        ax.plot([np_tth.min(), np_tth.max()], [param, param], "k:")
        ax.plot(np_tth, np_sum_1d - np_sum_m_1d+param, "r-", alpha=0.7,
                label="difference")

        if self.is_attribute("intensity_bkg_calc"):
            np_bkg = self.intensity_bkg_calc
            np_bkg_1d = numpy.where(numpy.isnan(np_sum), 0., np_bkg).sum(
                axis=1)
            ax.plot(np_tth, 2*np_bkg_1d, "b:", label="background")

        ax.legend(loc='upper right')
        fig.tight_layout()
        return (fig, ax)

    def plot_projection_diff(self):
        """Plot experimental polarized intensity vs. 2 theta (degrees)
        """
        if not(self.is_attribute("ttheta") & self.is_attribute("intensity_up") & 
               self.is_attribute("intensity_up_sigma") &
               self.is_attribute("intensity_down") & 
               self.is_attribute("intensity_down_sigma")):
            return
        fig, ax = plt.subplots()
        ax.set_title("Projection of polarized intensity on 2 theta axis")
        ax.set_xlabel("2 theta (degrees)")
        ax.set_ylabel('Intensity')


        np_tth = numpy.array(self.ttheta, dtype=float)
        np_up = self.intensity_up
        np_sup = self.intensity_up_sigma
        np_down = self.intensity_down
        np_sdown = self.intensity_down_sigma
        np_diff = np_up - np_down
        np_sdiff = numpy.sqrt(numpy.square(np_sup)+numpy.square(np_sdown))

        np_diff_1d = numpy.where(numpy.isnan(np_diff), 0., np_diff).sum(axis=1)
        np_sdiff_1d = ((numpy.where(numpy.isnan(np_diff), 0., np_sdiff)**2
                        ).sum(axis=1))**0.5

        np_up_m = self.intensity_up_total
        np_down_m = self.intensity_down_total
        np_diff_m = np_up_m - np_down_m

        np_diff_m_1d = numpy.where(numpy.isnan(np_diff), 0., np_diff_m).sum(
            axis=1)

        ax.plot([np_tth.min(), np_tth.max()], [0., 0.], "b:")
        ax.plot(np_tth, np_diff_m_1d, "k-", label="model")
        ax.errorbar(np_tth, np_diff_1d, yerr=np_sdiff_1d, fmt="ko", alpha=0.2,
                        label="experiment")
        y_min_d, y_max_d = ax.get_ylim()
        param = y_min_d-(np_diff_1d-np_diff_m_1d).max()

        ax.plot([np_tth.min(), np_tth.max()], [param, param], "k:")
        ax.plot(np_tth, np_diff_1d-np_diff_m_1d+param, "r-", alpha=0.7,
                    label="difference")

        ax.legend(loc='upper right')
        fig.tight_layout()
        return (fig, ax)

    def plot_ttheta_phi(self):
        if not(self.is_attribute("ttheta") & self.is_attribute("phi") &
               self.is_attribute("intensity_up_total") & 
               self.is_attribute("intensity_down_total")):
            return

        fig, axs = plt.subplots(2, 2)
        ax_1, ax_2 = axs[:, 1]
        ax_3, ax_4 = axs[:, 0]
        ax_1.set_title("Unpolarized intensity (model)")
        ax_1.set_xlabel("2 theta (degrees)")
        ax_1.set_ylabel('phi (degrees)')

        np_tth = numpy.array(self.ttheta, dtype=float)
        np_phi = numpy.array(self.phi, dtype=float)
        np_up_m = self.intensity_up_total
        np_down_m = self.intensity_down_total
        np_sum_m = np_up_m + np_down_m
        np_diff_m = np_up_m - np_down_m

        np_up = self.intensity_up
        np_down = self.intensity_down
        np_sum = np_up + np_down
        np_diff = np_up - np_down
        
        max_val = numpy.nanmax(np_sum_m)
        min_val = numpy.nanmin(np_sum_m)
        norm_1 = matplotlib.colors.Normalize(
            vmax=max_val, vmin=min_val)
        # cmap = matplotlib.PRGn
        
        max_val = numpy.nanmax(numpy.abs(np_sum_m))
        norm_1 = matplotlib.colors.Normalize(
            vmax=max_val, vmin=-max_val)
        # cmap = matplotlib.PRGn
        
        extent =(np_tth.min(), np_tth.max(), np_phi.min(), np_phi.max())
        plt.set_cmap('rainbow')
        ax_1.imshow(np_sum_m.transpose(), interpolation='bilinear',
                    extent=extent, alpha=0.9, origin="lower", norm=norm_1)

        ax_2.set_title("Polarized intensity (model)")
        ax_2.set_xlabel("2 theta (degrees)")
        ax_2.set_ylabel('phi (degrees)')

        max_val = numpy.nanmax(numpy.abs(np_diff_m))
        norm_2 = matplotlib.colors.Normalize(
            vmax=max_val, vmin=-max_val)
        ax_2.imshow(np_diff_m.transpose(), interpolation='bilinear',
                    extent=extent, alpha=0.9, origin="lower", norm=norm_2)

        ax_3.set_title("Unpolarized intensity (experiment)")
        ax_3.set_xlabel("2 theta (degrees)")
        ax_3.set_ylabel('phi (degrees)')
        
        plt.set_cmap('rainbow')
        ax_3.imshow(np_sum.transpose(), interpolation='bilinear',
                    extent=extent, alpha=0.9, origin="lower", norm=norm_1)

        ax_4.set_title("Polarized intensity (experiment)")
        ax_4.set_xlabel("2 theta (degrees)")
        ax_4.set_ylabel('phi (degrees)')

        ax_4.imshow(np_diff.transpose(), interpolation='bilinear',
                    extent=extent, alpha=0.9, origin="lower", norm=norm_2)

        return (fig, ax_1)

    def plot_diff_total(self):
        if not(self.is_attribute("ttheta") & self.is_attribute("phi") &
               self.is_attribute("intensity_up_total") & 
               self.is_attribute("intensity_down_total")):
            return
        fig, ax = plt.subplots()
        ax.set_title("Polarized intensity")
        ax.set_xlabel("2 theta (degrees)")
        ax.set_ylabel('phi (degrees)')

        np_tth = numpy.array(self.ttheta, dtype=float)
        np_phi = numpy.array(self.phi, dtype=float)
        np_up = self.intensity_up_total
        np_down = self.intensity_down_total
        np_diff = np_up - np_down
        
        norm = matplotlib.colors.Normalize(vmax=abs(np_diff).max(),
                                           vmin=-abs(np_diff).max())
        # cmap = matplotlib.PRGn
        
        extent =(np_tth.min(), np_tth.max(), np_phi.min(), np_phi.max())
        plt.set_cmap('rainbow')
        ax.imshow(np_diff.transpose(), interpolation='bilinear', extent=extent,
                  alpha=0.9, origin="lower", norm=norm)        
        
        # blk = '#000000'
        # ax.contour(np_tth, np_phi, np_sum.transpose(),
        #             #levels=[0.1, 0.5, 1., 5., 10., 50.],
        #             levels=5,
        #             colors=[blk, blk, blk, blk, blk, blk],
        #             linewidths=0.5)
        return (fig, ax)

    def plot_gamma_nu(self):
        if not(self.is_attribute("ttheta") & self.is_attribute("phi") &
               self.is_attribute("intensity_up_total") & 
               self.is_attribute("intensity_down_total")):
            return
        np_gamma, np_nu, ints = self.recalc_to_gamma_nu_grid()
        [np_u, np_d, np_sum, np_diff,
         np_u_m, np_d_m, np_sum_m, np_diff_m]= ints
        fig, axs = plt.subplots(2, 2)
        ax_1, ax_2 = axs[:, 1]
        ax_3, ax_4 = axs[:, 0]
        ax_1.set_title("Unpolarized intensity (model)")
        ax_1.set_xlabel("gamma (degrees)")
        ax_1.set_ylabel('nu (degrees)')

        max_val = numpy.nanmax(np_sum_m)
        min_val = numpy.nanmin(np_sum_m)
        norm_1 = matplotlib.colors.Normalize(
            vmax=max_val, vmin=min_val)
        # cmap = matplotlib.PRGn
        
        extent =(np_gamma.min(), np_gamma.max(),
                 np_nu.min(), np_nu.max())
        plt.set_cmap('rainbow')
        ax_1.imshow(np_sum_m, interpolation='bilinear', extent=extent,
                  alpha=0.9, origin="lower", norm=norm_1)

        ax_2.set_title("Polarized intensity (model)")
        ax_2.set_xlabel("gamma (degrees)")
        ax_2.set_ylabel('nu (degrees)')

        max_val = numpy.nanmax(numpy.abs(np_diff_m))
        norm_2 = matplotlib.colors.Normalize(
            vmax=max_val, vmin=-max_val)
        ax_2.imshow(np_diff_m, interpolation='bilinear', extent=extent,
                  alpha=0.9, origin="lower", norm=norm_2)


        ax_3.set_title("Unpolarized intensity (experiment)")
        ax_3.set_xlabel("gamma (degrees)")
        ax_3.set_ylabel('nu (degrees)')
        
        plt.set_cmap('rainbow')
        ax_3.imshow(np_sum, interpolation='bilinear', extent=extent,
                  alpha=0.9, origin="lower", norm=norm_1)

        ax_4.set_title("Polarized intensity (experiment)")
        ax_4.set_xlabel("gamma (degrees)")
        ax_4.set_ylabel('nu (degrees)')

        ax_4.imshow(np_diff, interpolation='bilinear', extent=extent,
                  alpha=0.9, origin="lower", norm=norm_2)

        
        # blk = '#000000'
        # ax.contour(np_tth, np_phi, np_sum.transpose(),
        #             #levels=[0.1, 0.5, 1., 5., 10., 50.],
        #             levels=5,
        #             colors=[blk, blk, blk, blk, blk, blk],
        #             linewidths=0.5)
        return (fig, ax_1)

# s_cont = """
#  _pd2d_proc_2theta_phi_intensity_up_net
#  ;
#       2    4.5     40.0     80.0
#  -3.000 -356.0   -350.0   -400.0
#  41.000 -357.0   -350.0   -400.0
#  ;
#  _pd2d_proc_2theta_phi_intensity_down_net
#  ;
#       2    4.5     40.0     80.0
#  -3.000 -356.0   -350.0   -400.0
#  41.000 -357.0   -350.0   -400.0
#  ;
#  _pd2d_proc_2theta_phi_intensity_up_total
#  ;
#       2    4.5     40.0     80.0
#  -3.000 -356.0   -350.0   -400.0
#  41.000 -357.0   -350.0   -400.0
#  ;
#  _pd2d_proc_2theta_phi_intensity_down_total
#  ;
#       2    4.5     40.0     80.0
#  -3.000 -356.0   -350.0   -400.0
#  41.000 -357.0   -350.0   -400.0
#  ;
#  _pd2d_proc_2theta_phi_intensity_bkg_calc
#  ;
#       2    4.5     40.0     80.0
#  -3.000 -356.0   -350.0   -400.0
#  41.000 -357.0   -350.0   -400.0
#  ;
#  _pd2d_proc_2theta_phi_intensity_up
#  ;
#       2    4.5     40.0     80.0
#  -3.000 -356.0   -350.0   -400.0
#  41.000 -357.0   -350.0   -400.0
#  ;
#  _pd2d_proc_2theta_phi_intensity_up_sigma
#  ;
#       2    4.5     40.0     80.0
#  -3.000 -356.0   -350.0   -400.0
#  41.000 -357.0   -350.0   -400.0
#  ;
#  _pd2d_proc_2theta_phi_intensity_down
#  ;
#       2    4.5     40.0     80.0
#  -3.000 -356.0   -350.0   -400.0
#  41.000 -357.0   -350.0   -400.0
#  ;
#  _pd2d_proc_2theta_phi_intensity_down_sigma
#  ;
#       2    4.5     40.0     80.0
#  -3.000 -356.0   -350.0   -400.0
#  41.000 -357.0   -350.0   -400.0
#  ;

# """

# obj = Pd2dProc.from_cif(s_cont)
# print(obj, end="\n\n")
# print(obj.ttheta, end="\n\n")
# print(obj.phi, end="\n\n")
# print(obj.intensity_bkg_calc, end="\n\n")
