import numpy
import scipy
import scipy.special
from typing import NoReturn

from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary
from cryspy.A_functions_base.function_1_tof import tof_Jorgensen, \
    tof_Jorgensen_VonDreele, calc_hpv_eta, calc_sigma_gamma

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN


class TOFProfile(ItemN):
    """Profile parameters for time-of-flight experiment.

    Attributes
    ----------
        - sigma0, sigma1, sigma2, gamma0, gamma1, gamma2
          alpha0, alpha1, beta0, beta1 (mandatory)
        - size_g, size_l, strain_g, strain_l, peak_shape (optional)

    The variance of the Gaussian component of the peak shape in TOF neutrons is
    given by:

        H_G**2 = (sigma_2+size_g)*d**4 + (sigma_1+strain_g)*d**2 + sigma_0
        H_L = (gamma_2+size_l)*d**2 + (gamma_1+strain_l)*d + gamma_0

    where d is the d-spacing in angstroms.

    Units
    -----
        sigma2, gamma2 are (microsecs/Ang.**2)**2;
        sigma1, gamma1 are (microsecs/Ang.**2);
        sigma0, gamma0 are (microsecs**2);

    """
    ATTR_MANDATORY_NAMES = ("peak_shape", )
    ATTR_MANDATORY_TYPES = (str, )
    ATTR_MANDATORY_CIF = ("peak_shape", )

    ATTR_OPTIONAL_NAMES = (
        "sigma0", "sigma1", "sigma2", "size_g", "strain_g",
        "gamma0", "gamma1", "gamma2", "size_l", "strain_l",
        "alpha0", "alpha1", "alpha2", "beta0", "beta1", "beta00", "beta01", "beta10",
        "r01", "r02", "r03")
    ATTR_OPTIONAL_TYPES = (
        float, float, float, float, float,
        float, float, float, float, float,
        float, float, float, float, float, float, float, float,
        float, float, float)
    ATTR_OPTIONAL_CIF = (
        "sigma0", "sigma1", "sigma2", "size_g", "strain_g",
        "gamma0", "gamma1", "gamma2", "size_l", "strain_l",
        "alpha0", "alpha1", "alpha2", "beta0", "beta1", "beta00", "beta01", "beta10",
        "r01", "r02", "r03")

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ()
    ATTR_INT_PROTECTED_NAMES = ()

    # parameters considered are refined parameters
    ATTR_REF = ("sigma2", "sigma1", "sigma0", "gamma2", "gamma1", "gamma0",
                "alpha0", "alpha1", "beta0", "beta1", "size_g", "size_l",
                "strain_g", "strain_l", "beta00", "beta01", "beta10",
                "r01", "r02", "r03", "alpha2")
    ATTR_SIGMA = tuple([f"{_h:}_sigma" for _h in ATTR_REF])
    ATTR_CONSTR_FLAG = tuple([f"{_h:}_constraint" for _h in ATTR_REF])
    ATTR_REF_FLAG = tuple([f"{_h:}_refinement" for _h in ATTR_REF])
    ATTR_CONSTR_MARK = tuple([f"{_h:}_mark" for _h in ATTR_REF])

    # formats if cif format
    D_FORMATS = {"sigma0": "{:.5f}", "sigma1": "{:.5f}", "sigma2": "{:.5f}",
                 "gamma0": "{:.5f}", "gamma1": "{:.5f}", "gamma2": "{:.5f}",
                 "alpha0": "{:.5f}", "alpha1": "{:.5f}", "beta0": "{:.5f}",
                 "beta1": "{:.5f}", "size_g": "{:.5f}", "size_l": "{:.5f}",
                 "strain_g": "{:.5f}", "strain_l": "{:.5f}",
                 "beta00": "{:.5f}", "beta01": "{:.5f}", "beta10": "{:.5f}",
                 "r01": "{:.5f}", "r02": "{:.5f}", "r03": "{:.5f}",
                 "alpha2": "{:.5f}"}

    # constraints on the parameters
    D_CONSTRAINTS = {"peak_shape": ["pseudo-Voigt", "Gauss", "type0m"]}

    # default values for the parameters
    D_DEFAULT = {"peak_shape": "pseudo-Voigt"}

    for key in ATTR_SIGMA:
        D_DEFAULT[key] = 0.
    for key in (ATTR_CONSTR_FLAG + ATTR_REF_FLAG):
        D_DEFAULT[key] = False
    for key in ATTR_CONSTR_MARK:
        D_DEFAULT[key] = ""

    PREFIX = "tof_profile"

    def __init__(self, **kwargs) -> NoReturn:
        super(TOFProfile, self).__init__()

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

    def to_cif(self, separator: str = "_", flag_all_attributes:bool=False) -> str:
        """
        Redefined function to_cif
        """
        d_attributes_peak_shape =  {
            "pseudo-Voigt": [
                "peak_shape",
                "sigma0", "sigma1", "sigma2", "size_g", "strain_g",
                "gamma0", "gamma1", "gamma2", "strain_l", "size_l",
                "alpha0", "alpha1",
                "beta0", "beta1"],
            "Gauss": [
                "peak_shape",
                "sigma0", "sigma1", "sigma2", "size_g", "strain_g",
                "alpha0", "alpha1",
                "beta0", "beta1"],
            "type0m": [
                "peak_shape",
                "sigma0", "sigma1", "sigma2",
                "gamma0", "gamma1", "gamma2",
                "alpha1", "alpha2",
                "beta00", "beta01", "beta10",
                "r01", "r02", "r03"]}
        peak_shape = getattr(self, "peak_shape")
        l_attr_print = d_attributes_peak_shape[peak_shape]
        # REDO: the lines below are fully copied from cl_1_item.py, but the original code should be re-used instead of copying!!!
        ls_out = []
        prefix = self.PREFIX
        for name, name_cif in zip(self.ATTR_NAMES, self.ATTR_CIF):
            if name in l_attr_print:
                if name_cif == "":
                    name_cif = name
                flag_value = self.is_attribute(name)
                if (flag_value | flag_all_attributes):
                    if flag_value:
                        s_val = str(getattr(self, name))
                    else:
                        s_val = "."
                    l_s_val = s_val.split("\n")
                    if len(l_s_val) > 1:
                        ls_out.append(f"_{prefix:}{separator:}{name_cif:}")
                        ls_out.append(";")
                        ls_out.extend(l_s_val)
                        ls_out.append(";")
                    else:
                        if len(s_val.split(" ")) > 1:
                            ls_out.append(
                                f"_{prefix:}{separator:}{name_cif:} \"{s_val:}\"")
                        else:
                            s_val = getattr(self, f"{name:}_as_string")
                            ls_out.append(
                                f"_{prefix:}{separator:}{name_cif:} {s_val:}")
        return "\n".join(ls_out)
    # def calc_hpveta(self, h_g, h_l):
    #     """
    #     ttheta in radians, could be array
    #     pseudo-Voight function
    #     """
    #     h_pv, eta = calc_hpv_eta(h_g, h_l)
    #     return h_pv, eta
    #
    # def calc_agbg(self, h_pv):
    #     a_g = (2./h_pv)*(numpy.log(2.)/numpy.pi)**0.5
    #     b_g = 4*numpy.log(2)/(h_pv**2)
    #     return a_g, b_g
    #
    # def calc_albl(self, h_pv):
    #     a_l = 2./(numpy.pi*h_pv)
    #     b_l = 4./(h_pv**2)
    #     return a_l, b_l
    #
    #
    # def calc_peak_shape_function(
    #         self, d, time, time_hkl, size_g: float = 0., strain_g: float = 0.,
    #         size_l: float = 0., strain_l: float = 0.):
    #     """Calculate peak-shape function F(DELTA)
    #
    #     For Gauss peak-shape:
    #         F(DELTA) = alpha * beta / (alpha+beta) * [exp(u) erfc(y) + exp(v) erfc(z)]
    #
    #     For pseudo-Voigt peak-shape:
    #         F(DELTA) = ...
    #
    #     alpha = alpha0 + alpha1 / d
    #     beta = beta0 + beta1 / d**4
    #     sigma = sigma0 + sigma1 * d
    #
    #     u = alpha/2 * (alpha * sigma**2 + 2 DELTA)
    #     v = beta/2 * (beta * sigma**2 - 2 DELTA)
    #     y = (alpha * sigma**2 + DELTA)/(sigma * 2**0.5)
    #     z = (beta * sigma**2 - DELTA)/(sigma * 2**0.5)
    #
    #     DELTA = time - time_hkl
    #
    #     """
    #
    #     alpha = self.alpha0 + self.alpha1 / d
    #     beta = self.beta0 + self.beta1 / d**4
    #
    #     if self.peak_shape == "Gauss":
    #         sigma = self.sigma0 + self.sigma1 * d
    #         res_2d = tof_Jorgensen(alpha, beta, sigma, time, time_hkl)
    #     else:  # self.peak_shape == "pseudo-Voigt"
    #         if self.is_attribute("size_g"):
    #             size_g_s = self.size_g+size_g
    #         else:
    #             size_g_s = size_g
    #         if self.is_attribute("size_l"):
    #             size_l_s = self.size_l+size_l
    #         else:
    #             size_l_s = size_l
    #         if self.is_attribute("strain_g"):
    #             strain_g_s = self.strain_g+strain_g
    #         else:
    #             strain_g_s = strain_g
    #         if self.is_attribute("strain_l"):
    #             strain_l_s = self.strain_l+strain_l
    #         else:
    #             strain_l_s = strain_l
    #
    #         sigma, gamma = calc_sigma_gamma(
    #             d, self.sigma0, self.sigma1, self.sigma2, self.gamma0,
    #             self.gamma1, self.gamma2, size_g=size_g_s, size_l=size_l_s,
    #             strain_g=strain_g_s, strain_l=strain_l_s)
    #
    #         res_2d = tof_Jorgensen_VonDreele(
    #             alpha, beta, sigma, gamma, time, time_hkl)
    #     return res_2d

class TOFProfileL(LoopN):
    """Profile parameters for time-of-flight experiment.

    """
    ITEM_CLASS = TOFProfile
    ATTR_INDEX = None
    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(TOFProfileL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name

