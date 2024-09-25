import numpy
import scipy
import scipy.special
from typing import NoReturn

from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary

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

    def get_dictionary(self):
        res = {}
        if self.is_attribute("peak_shape"):
            peak_shape = self.peak_shape
        else:
            peak_shape = "pseudo-Voigt"
        res["profile_peak_shape"] = peak_shape
        l_sigma = []
        l_sigma_refinement = []
        for numb in range(3):
            if self.is_attribute(f"sigma{numb:}"):
                l_sigma.append(getattr(self, f"sigma{numb:}"))
                l_sigma_refinement.append(getattr(self, f"sigma{numb:}_refinement"))
            else:
                break
        res["profile_sigmas"] = numpy.array(l_sigma, dtype=float)
        res["flags_profile_sigmas"] = numpy.array(l_sigma_refinement, dtype=float)
        if peak_shape in ["pseudo-Voigt", "type0m"]:
            l_gamma = []
            l_gamma_refinement = []
            for numb in range(3):
                if self.is_attribute(f"gamma{numb:}"):
                    l_gamma.append(getattr(self, f"gamma{numb:}"))
                    l_gamma_refinement.append(getattr(self, f"gamma{numb:}_refinement"))
                else:
                    break
            res["profile_gammas"] = numpy.array(l_gamma, dtype=float)
            res["flags_profile_gammas"] = numpy.array(l_gamma_refinement, dtype=float)
        if peak_shape in ["pseudo-Voigt", "Gauss", ]:
            l_alpha = []
            l_alpha_refinement = []
            for numb in range(2):
                if self.is_attribute(f"alpha{numb:}"):
                    l_alpha.append(getattr(self, f"alpha{numb:}"))
                    l_alpha_refinement.append(getattr(self, f"alpha{numb:}_refinement"))
                else:
                    break
            l_beta = []
            l_beta_refinement = []
            for numb in range(2):
                if self.is_attribute(f"beta{numb:}"):
                    l_beta.append(getattr(self, f"beta{numb:}"))
                    l_beta_refinement.append(getattr(self, f"beta{numb:}_refinement"))
                else:
                    break
            res["profile_alphas"] = numpy.array(l_alpha, dtype=float)
            res["flags_profile_alphas"] = numpy.array(l_alpha_refinement, dtype=float)
            res["profile_betas"] = numpy.array(l_beta, dtype=float)
            res["flags_profile_betas"] = numpy.array(l_beta_refinement, dtype=float)
        elif peak_shape in ["type0m", ]:
            profile_alphas = numpy.array([
                getattr(self, "alpha1"), getattr(self, "alpha2")
                ], dtype=float)
            flags_profile_alphas = numpy.array([
                getattr(self, "alpha1_refinement"), getattr(self, "alpha2_refinement")
                ], dtype=bool)
            profile_betas = numpy.array([
                getattr(self, "beta00"), getattr(self, "beta01"), getattr(self, "beta10")
                ], dtype=float)
            flags_profile_betas = numpy.array([
                getattr(self, "beta00_refinement"),
                getattr(self, "beta01_refinement"),
                getattr(self, "beta10_refinement")
                ], dtype=bool)
            profile_rs = numpy.array([
                getattr(self, "r01"), getattr(self, "r02"), getattr(self, "r03")
                ], dtype=float)
            flags_profile_rs = numpy.array([
                getattr(self, "r01_refinement"),
                getattr(self, "r02_refinement"),
                getattr(self, "r03_refinement")
                ], dtype=bool)
            res["profile_alphas"] = profile_alphas
            res["flags_profile_alphas"] = flags_profile_alphas
            res["profile_betas"] = profile_betas
            res["flags_profile_betas"] = flags_profile_betas
            res["profile_rs"] = profile_rs
            res["flags_profile_rs"] = flags_profile_rs
        return res

class TOFProfileL(LoopN):
    """Profile parameters for time-of-flight experiment.

    """
    ITEM_CLASS = TOFProfile
    ATTR_INDEX = None
    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(TOFProfileL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name

