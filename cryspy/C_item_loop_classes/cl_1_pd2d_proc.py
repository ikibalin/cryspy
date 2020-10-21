import numpy
from typing import NoReturn
from cryspy.A_functions_base.function_1_gamma_nu import \
    recal_int_to_gammanu_grid
from cryspy.B_parent_classes.cl_1_item import ItemN

class Pd2dProc(ItemN):
    """
    Pd2dProc class

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

        min_nu, max_nu = -10.*numpy.pi/180., 15.*numpy.pi/180.
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
