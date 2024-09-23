from typing import NoReturn
import math
import numpy
import matplotlib.pyplot as plt

from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN


class DiffrnRefln(ItemN):
    """
    The flip ratios measured in the single diffraction experiment.

    Data items in the DIFFRN_REFLN category record details about
    the intensities measured in the diffraction experiment.

    The DIFFRN_REFLN data items refer to individual intensity
    measurements and must be included in looped lists.

    Attributes:
        - 

    """
    ATTR_MANDATORY_NAMES = ("index_h", "index_k", "index_l")
    ATTR_MANDATORY_TYPES = (int, int, int)
    ATTR_MANDATORY_CIF = ("index_h", "index_k", "index_l")

    ATTR_OPTIONAL_NAMES = ("fr", "fr_sigma", "fr_calc", "intensity", "intensity_sigma", "intensity_calc",
                           "excluded", "sintlambda")
    ATTR_OPTIONAL_TYPES = (float, float, float, float, float, float, bool, float)
    ATTR_OPTIONAL_CIF = ("fr", "fr_sigma", "fr_calc", "intensity", "intensity_sigma", "intensity_calc",
                         "excluded", "sint/lambda")

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
    D_FORMATS = {"fr_calc": "{:.5f}", "intensity_calc": "{:.2f}",
                 "sintlambda": "{:.5f}",
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

    PREFIX = "diffrn_refln"

    def __init__(self, **kwargs) -> NoReturn:
        super(DiffrnRefln, self).__init__()

        # defined for any integer and float parameters
        D_MIN = {"fr": 0., "fr_sigma": 0., "fr_calc": 0.,
                 "intensity_calc": 0.,
                 "sintlambda": 0.,
                 "intensity": 0., "intensity_sigma":0.}

        # defined for ani integer and float parameters
        D_MAX = {}

        self.__dict__["D_MIN"] = D_MIN
        self.__dict__["D_MAX"] = D_MAX
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)


class DiffrnReflnL(LoopN):
    """
    Flip ratios measured in the single diffraction experiment.

    """
    ITEM_CLASS = DiffrnRefln
    ATTR_INDEX = None
    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(DiffrnReflnL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name

    def report(self):
        return self.report_agreement_factor_exp() + "\n" + self.report_chi_sq_exp()

    def report_agreement_factor_exp(self):
        """
        Make a report about experimental agreement factor in string format.
        """

        l_chi_sq_exp, l_ag_f_exp = [], []
        l_diff = []
        n_3s, n_2s, n_1s, n_0s = 0, 0, 0, 0
        l_hkl = [(int(_1), int(_2), int(_3)) for _1, _2, _3 in
                 zip(self.index_h, self.index_k, self.index_l)]
        for _hkl, _fr_1, _fr_sigma_1 in zip(l_hkl, self.fr, self.fr_sigma):
            _diff = abs(_fr_1-1.)/float(_fr_sigma_1)
            if _diff >= 3.:
                n_3s += 1
            elif _diff >= 2.:
                n_2s += 1
            elif _diff >= 1.:
                n_1s += 1
            else:
                n_0s += 1
            if _fr_1 <= 1.:
                if math.isclose(_fr_1, 0.):
                    l_diff.append(1./(_fr_1+_fr_sigma_1)-1.)
                else:
                    l_diff.append(1./_fr_1-1.)
            else:
                l_diff.append(_fr_1-1.)

            _mhkl = (-1 * _hkl[0], -1 * _hkl[1], -1 * _hkl[2])
            if _mhkl in l_hkl:
                ind_mhkl = l_hkl.index(_mhkl)
                _fr_2, _fr_sigma_2 = self.fr[ind_mhkl], self.fr_sigma[ind_mhkl]
                _fr_sigma = 1. / (_fr_sigma_1 ** (-2) + _fr_sigma_2 ** (-2)) ** 0.5
                _fr_average = (_fr_1 * _fr_sigma_1 ** (-2) + _fr_2 * _fr_sigma_2 ** (-2)) * _fr_sigma ** 2
                delta_fr = abs(_fr_1 - _fr_average)
                chi_sq_exp = (delta_fr / _fr_sigma_1) ** 2
                l_chi_sq_exp.append(chi_sq_exp)
                if math.isclose(_fr_1-1., 0):
                    ag_f_exp = abs((_fr_1 - _fr_average) / (_fr_1 - 1.+_fr_sigma_1))
                else:
                    ag_f_exp = abs((_fr_1 - _fr_average) / (_fr_1 - 1.))
                l_ag_f_exp.append(ag_f_exp)
                # print("hkl: {:4} {:4} {:4}".format(_hkl[0], _hkl[1], _hkl[2]))
                # print("chi_sq_exp: {:.3f} ".format(chi_sq_exp))
                # print("ag_f_exp: {:.3f} ".format(ag_f_exp))
        ls_out = ["# Experimental agreement factor\n"]
        n = n_0s + n_1s + n_2s + n_3s
        ls_out.append(f"\n| Number of measured reflections:| {n:4}|")
        ls_out.append(f"  |    range of h is |{min([_[0] for _ in l_hkl]):3},{max([_[0] for _ in l_hkl]):3} |")
        ls_out.append(f"  |             k is |{min([_[1] for _ in l_hkl]):3},{max([_[1] for _ in l_hkl]):3} |")
        ls_out.append(f"  |             l is |{min([_[2] for _ in l_hkl]):3},{max([_[2] for _ in l_hkl]):3} |")
        ls_out.append(f" |max(FR_exp - 1) is |{max(l_diff):5.3f} |")
        ls_out.append("\n N+1 > abs(FR_exp - 1)/FR_sigma > N: ")
        ls_out.append(f" |abs(FR_exp - 1)/FR_sigma < 1: |{n_0s:4}| {100*float(n_0s)/float(n):5.1f}% | ")
        ls_out.append(f" |                    N = 1: |{n_1s:4}| {100*float(n_1s)/float(n):5.1f}%  |")
        ls_out.append(f" |                    N = 2: |{n_2s:4}| {100 * float(n_2s) / float(n):5.1f}% |")
        ls_out.append(f" |abs(FR_exp - 1)/FR_sigma > 3: |{n_3s:4}| {100 * float(n_3s) / float(n):5.1f}% |")

        n_friedel = len(l_chi_sq_exp)
        ls_out.append(f"\n|Total number of Friedel reflections is |{n_friedel:}.|")
        if n_friedel != 0:
            ls_out.append(f"|  (abs(FR_exp-FR_av.)/FR_sigma)^2  is| {sum(l_chi_sq_exp) / n_friedel:.2f}|")
            ls_out.append(f"|   abs(FR_exp-FR_av.)/abs(FR_exp-1) per reflection is |{(100*sum(l_ag_f_exp)/n_friedel):.2f}% |")

        return "\n".join(ls_out)

    def calc_chi_sq_points(self) -> (float, int):
        """Calculate chi_sq and points if fr, fr_sigma, fr_calc are given."""
        if (self.is_attribute("fr") & self.is_attribute("fr_sigma") &
                self.is_attribute("fr_calc")):

            np_excl = numpy.array(self.excluded, dtype=bool)
            if numpy.all(np_excl):
                return (0., 0)

            flag_in = numpy.logical_not(np_excl)
            np_fr = numpy.array(self.fr, dtype=float)[flag_in]
            np_fr_calc = numpy.array(self.fr_calc, dtype=float)[flag_in]
            np_fr_sigma = numpy.array(self.fr_sigma, dtype=float)[flag_in]

            chi_sq = numpy.sum(numpy.square((np_fr-np_fr_calc)/np_fr_sigma))
            points = np_fr.size
            return chi_sq, points

        return (0., 0)


    def report_chi_sq_exp(self, cell=None) -> str:
        """
        Make a report about experimental chi_sq in string format.

        cell is unit cell object.
        """
        ls_out = []
        l_hkl = [(_1, _2, _3) for _1, _2, _3 in zip(self.index_h, self.index_k, self.index_l)]
        l_fr, l_fr_sigma, l_fr_calc = self.fr, self.fr_sigma, self.fr_calc
        if (l_fr is None) & (l_fr_sigma is None) & (l_fr_calc is None):
            return "\n".join(ls_out)
        n = len(l_fr)
        n_1s, n_2s, n_3s = 0, 0, 0
        l_chi_sq, l_worsest, l_af_f, l_af_r = [], [], [], []
        for _hkl, _fr, _fr_sigma, _fr_calc in zip(l_hkl, l_fr, l_fr_sigma, l_fr_calc):
            _diff = abs(float(_fr-_fr_calc)/float(_fr_sigma))
            l_chi_sq.append(_diff**2)
            if math.isclose(_fr, 0.):
                l_af_f.append(abs(float(_fr - _fr_calc) / float(_fr_sigma)))
            else:
                l_af_f.append(abs(float(_fr - _fr_calc) / float(_fr)))
            if not(math.isclose(_fr, 1.)):
                l_af_r.append(abs(float(_fr-_fr_calc)/float(_fr-1.)))
            if _diff <= 1.:
                n_1s +=1
            elif _diff <= 2.:
                n_2s += 1
            elif _diff <= 3.:
                n_3s += 1
            else:
                l_worsest.append((_hkl, _fr, _fr_sigma, _fr_calc, _diff))
        ls_out.append("# Chi_sq experimental")
        ls_out.append(f"Total number of reflections is {n:}")
        ls_out.append(f"|  (abs(FR_exp-FR_mod)/FR_sigma)^2 per reflection is |{sum(l_chi_sq)/float(n):.2f}|")
        ls_out.append(f"|   abs(FR_exp-FR_mod)/FR_exp      per reflection is |{100*sum(l_af_f) / float(n):.2f}%|")
        ls_out.append(f"|   abs(FR_exp-FR_mod)/abs(FR_exp-1)    per reflection is |{100*sum(l_af_r) / float(n):.2f}%|")
        ls_out.append("           (reflections with FR_exp = 1 are excluded)")
        n_worsest = len(l_worsest)
        ls_out.append("\n## Reflections in range  ")
        ls_out.append(" (N-1)*FR_sigma < abs(FR_exp - FR_mod) < N*FR_sigma: ")
        ls_out.append(f"|      N = 1:| {n_1s:}/{n:} ={100*float(n_1s)/float(n):5.1f}% ({2*34.1:4.1f}%, three sigma rule) |")
        ls_out.append(f"|      N = 2:| {n_2s:}/{n:} ={100 * float(n_2s) / float(n):5.1f}% ({2*(13.6):4.1f}%, three sigma rule)|")
        ls_out.append(f"|      N = 3:| {n_3s:}/{n:} ={100 * float(n_3s) / float(n):5.1f}% ({2*(2.1):4.1f}%, three sigma rule)|")
        ls_out.append(f"|      N > 3:| {n_worsest:}/{n:} ={100 * float(n_worsest) / float(n):5.1f}% ({2*(0.1):4.1f}%, three sigma rule)|")
        l_worsest.sort(key=lambda x: x[4], reverse=True)
        if len(l_worsest) > 1:
            if n_worsest > 10: n_worsest = 10
            ls_out.append("\n## The ten worsest reflections:")
            ls_out.append("|  h | k | l  |     FR |FR_sigma | FR_calc | diff|")
            for (_hkl, _fr, _fr_sigma, _fr_calc, _diff) in l_worsest[:n_worsest]:
                ls_out.append(f"|{_hkl[0]:3}|{_hkl[1]:3}|{_hkl[2]:3}|{_fr:9.5f}|{_fr_sigma:9.5f}|{_fr_calc:9.5f}|{_diff:6.1f}|")
        if cell is not None:
            np_h = numpy.array(self.index_h, dtype=int)
            np_k = numpy.array(self.index_k, dtype=int)
            np_l = numpy.array(self.index_l, dtype=int)
            np_sthovl = cell.calc_sthovl(index_h=np_h, index_k=np_k,
                                         index_l=np_l)
            np_fr, np_fr_sigma = numpy.array(l_fr, dtype=float), numpy.array(l_fr_sigma, dtype=float) 
            np_fr_calc = numpy.array(l_fr_calc, dtype=float)
            np_diff = numpy.abs((np_fr-np_fr_calc)/np_fr_sigma)
            n_bins = 10
            np_val, np_sthovl_bins = numpy.histogram(np_sthovl, bins=n_bins)
            ls_out.append("\n## Distribution of reflection in sin(theta)/lambda range")
            ls_out.append("|sthovl_1 |sthovl_2 | (Exp-Mod)/Sigma |n_points|")
            for _1, _2, _3 in zip(np_sthovl_bins[:-1], np_sthovl_bins[1:], np_val):
                np_flag = numpy.logical_and(np_sthovl >= _1, np_sthovl < _2) 
                res = np_diff[np_flag]
                if res.size == 0:
                    ls_out.append(f"|{_1:8.3f} |{_2:8.3f}|       None    |      0|")
                else:
                    ls_out.append(f"|{_1:8.3f} |{_2:8.3f}|    {np_diff[np_flag].mean():7.3f} |{_3:10}|")
            numpy.histogram(np_sthovl, bins=n_bins)
        return "\n".join(ls_out)

    def plots(self):
        return [self.plot_fr_vs_fr_calc(),
                self.plot_asymmetry_vs_asymmetry_calc(),
                self.plot_intensity_vs_intensity_calc(),]
    
    def plot_fr_vs_fr_calc(self):
        """Plot experimental fr vs. fr_calc
        """
        if not(self.is_attribute("fr") & self.is_attribute("fr_sigma") &
               self.is_attribute("fr_calc")):
            return 

        fig, ax = plt.subplots()
        np_excl = numpy.array(self.excluded, dtype=bool)
        flag_in = numpy.logical_not(np_excl)
        if numpy.all(np_excl):
            ax.set_title("Flip Ratio: I_plus / I_minus")
            fr_min = 0.
            fr_max = 10.
        else:
            np_fr = numpy.array(self.fr, dtype=float)[flag_in]
            np_fr_calc = numpy.array(self.fr_calc, dtype=float)[flag_in]
            np_fr_sigma = numpy.array(self.fr_sigma, dtype=float)[flag_in]

            chi_sq_per_n = numpy.square((np_fr-np_fr_calc)/np_fr_sigma).sum()/np_fr.shape[0]
            ax.set_title(r"Flip Ratio: $I_{plus} / I_{minus}$, $\chi^2/n=$" + f"{chi_sq_per_n:.2f}.")
    
            np_fr_1 = np_fr - np_fr_sigma
            np_fr_2 = np_fr + np_fr_sigma
    
            fr_min = min([min(np_fr_1), min(self.fr_calc)])
            fr_max = max([max(np_fr_2), max(self.fr_calc)])
            ax.plot([fr_min, fr_max], [fr_min, fr_max], "k:")
            ax.errorbar(np_fr_calc, np_fr, yerr=np_fr_sigma, fmt="ko", alpha=0.2)

        flag_excl = numpy.logical_not(flag_in)
        np_fr_excl = numpy.array(self.fr, dtype=float)[flag_excl]
        np_fr_calc_excl = numpy.array(self.fr_calc, dtype=float)[flag_excl]
        np_fr_sigma_excl = numpy.array(self.fr_sigma, dtype=float)[flag_excl]
        ax.errorbar(np_fr_calc_excl, np_fr_excl, yerr=np_fr_sigma_excl,
                    fmt="rs", alpha=0.2, label="excluded")

        ax.set_xlim(fr_min, fr_max)
        ax.set_ylim(fr_min, fr_max)
        
        ax.set_xlabel("Flip ratio (model)")
        ax.set_ylabel('Flip ratio (experiment)')
        ax.set_aspect(1)
        fig.tight_layout()
        return (fig, ax)

    def plot_asymmetry_vs_asymmetry_calc(self):
        """Plot experimental fr vs. fr_calc
        """
        if not(self.is_attribute("fr") & self.is_attribute("fr_sigma") &
               self.is_attribute("fr_calc")):
            return 

        fig, ax = plt.subplots()
        np_excl = numpy.array(self.excluded, dtype=bool)
        flag_in = numpy.logical_not(np_excl)
        np_fr = numpy.array(self.fr, dtype=float)[flag_in]
        np_fr_calc = numpy.array(self.fr_calc, dtype=float)[flag_in]
        np_fr_sigma = numpy.array(self.fr_sigma, dtype=float)[flag_in]

        asymmetry = (np_fr - 1.)/(np_fr + 1.)
        asymmetry_sigma = np_fr_sigma*numpy.sqrt(2.)*numpy.sqrt(numpy.square(np_fr)+1.) / \
            numpy.square(np_fr+1.)
        asymmetry_calc = (np_fr_calc - 1.)/(np_fr_calc + 1.)
        chi_sq_per_n = numpy.square((asymmetry-asymmetry_calc)/asymmetry_sigma).sum()/asymmetry.size

        ax.set_title(r"Asymmetry parameter: $\frac{I_{plus}-I_{minus}}{I_{plus}+I_{minus}}$, $\chi^2/n=$" + f"{chi_sq_per_n:.2f}.")

        ax.plot([-1, 1], [-1, 1], "k:")
        ax.errorbar(asymmetry_calc, asymmetry, yerr=asymmetry_sigma, fmt="ko",
                    alpha=0.2)
        ax.set_xlim(-1, 1)
        ax.set_ylim(-1, 1)

        ax.set_xlabel("Asymmetry (model)")
        ax.set_ylabel('Asymmetry (experiment)')
        ax.set_aspect(1)
        fig.tight_layout()
        return (fig, ax)

    def plot_intensity_vs_intensity_calc(self):
        """Plot experimental intensity vs. intensity_calc
        """
        if not(self.is_attribute("intensity") & self.is_attribute("intensity_sigma") &
               self.is_attribute("intensity_calc")):
            return 

        fig, ax = plt.subplots()
        np_excl = numpy.array(self.excluded, dtype=bool)
        flag_in = numpy.logical_not(np_excl)
        if numpy.all(np_excl):
            ax.set_title("Intensity")
            intensity_min = 0.
            intensity_max = 10.
        else:
            np_intensity = numpy.array(self.intensity, dtype=float)[flag_in]
            np_intensity_calc = numpy.array(self.intensity_calc, dtype=float)[flag_in]
            np_intensity_sigma = numpy.array(self.intensity_sigma, dtype=float)[flag_in]

            chi_sq_per_n = numpy.square((np_intensity-np_intensity_calc)/np_intensity_sigma).sum()/np_intensity.shape[0]
            agreement_factor = numpy.abs(np_intensity-np_intensity_calc).sum()/numpy.abs(np_intensity).sum() * 100.
            ax.set_title(r"Intensity, $\chi^2/n=$" + f"{chi_sq_per_n:.2f} "+ r" $RF^2=$"+f"{agreement_factor:.1f}%")
    
            np_intensity_1 = np_intensity - np_intensity_sigma
            np_intensity_2 = np_intensity + np_intensity_sigma
    
            intensity_min = min([min(np_intensity_1), min(self.intensity_calc)])
            intensity_max = max([max(np_intensity_2), max(self.intensity_calc)])
            ax.plot([intensity_min, intensity_max], [intensity_min, intensity_max], "k:")
            ax.errorbar(np_intensity_calc, np_intensity, yerr=np_intensity_sigma, fmt="ko", alpha=0.2)

        flag_excl = numpy.logical_not(flag_in)
        np_intensity_excl = numpy.array(self.intensity, dtype=float)[flag_excl]
        np_intensity_calc_excl = numpy.array(self.intensity_calc, dtype=float)[flag_excl]
        np_intensity_sigma_excl = numpy.array(self.intensity_sigma, dtype=float)[flag_excl]
        ax.errorbar(np_intensity_calc_excl, np_intensity_excl, yerr=np_intensity_sigma_excl,
                    fmt="rs", alpha=0.2, label="excluded")

        ax.set_xlim(intensity_min, intensity_max)
        ax.set_ylim(intensity_min, intensity_max)
        
        ax.set_xlabel("Intensity (model)")
        ax.set_ylabel('Intensity (experiment)')
        ax.set_aspect(1)
        fig.tight_layout()
        return (fig, ax)


    def include_all_points(self):
        for item in self.items:
            item.excluded = False
        name = "numpy_excluded"
        if name in self.__dict__.keys():
            del self.__dict__[name]

    def exclude_by_fr_max(self, fr_max: float = 3.):
        """Exclude points by rules
        1./fr_max <= FR <= fr_max
        """
        if fr_max >= 1.:
            fr_1 = 1./fr_max
            fr_2 = fr_max
        else:
            fr_1 = fr_max
            fr_2 = 1./fr_max

        for item in self.items:
            if ((item.fr < fr_1) | (item.fr > fr_2)):
                item.excluded = True

        name = "numpy_excluded"
        if name in self.__dict__.keys():
            del self.__dict__[name]

    def get_dictionary(self):
        res = {}
        index_hkl = numpy.array([
            self.index_h, self.index_k,
            self.index_l], dtype=float)
        refln_fr_excl = numpy.array(self.excluded, dtype=bool)
        res["index_hkl"] = index_hkl
        res["flip_ratio_excluded"] = refln_fr_excl
        try:
            refln_fr_es = numpy.array([
                self.fr, self.fr_sigma], dtype=float)
            res["flip_ratio_es"] = refln_fr_es
        except AttributeError:
            refln_intensity_es = numpy.array([
                self.intensity, self.intensity_sigma], dtype=float)
            res["intensity_es"] = refln_intensity_es
        return res

# s_cont = """
#   loop_
#   _diffrn_refln_index_h
#   _diffrn_refln_index_k
#   _diffrn_refln_index_l
#   _diffrn_refln_fr
#   _diffrn_refln_fr_sigma
#       0    0    8   0.64545   0.01329 
#       2    0    6   1.75682   0.0454
# """

# obj = DiffrnReflnL.from_cif(s_cont)
# print(obj, end="\n\n")
# print(obj.report_agreement_factor_exp(), end="\n\n")
# print(obj.numpy_fr_sigma, end="\n\n")
