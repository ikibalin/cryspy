import numpy
from typing import NoReturn

from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN


class TOFBackground(ItemN):
    """Background description for time-of-flight experiment.

    Attributes
    ----------
        - time_max, coeff1 (mandatory)
        - coeff2, coeff3, coeff4, coeff5, coeff6, coeff7, coeff8,
          coeff9, coeff10, coeff11, coeff12, coeff13, coeff14, coeff15,
          coeff16, coeff17, coeff18, id (optional)
    """
    ATTR_MANDATORY_NAMES = ("time_max", "coeff1")
    ATTR_MANDATORY_TYPES = (float, float)
    ATTR_MANDATORY_CIF = ("time_max", "coeff1")

    ATTR_OPTIONAL_NAMES = (
        "coeff2", "coeff3", "coeff4", "coeff5", "coeff6", "coeff7", "coeff8",
        "coeff9", "coeff10", "coeff11", "coeff12", "coeff13", "coeff14",
        "coeff15", "coeff16", "coeff17", "coeff18", "id")

    ATTR_OPTIONAL_TYPES = (
        float, float, float, float, float, float, float, float, float, float,
        float, float, float, float, float, float, float, str)
    ATTR_OPTIONAL_CIF = (
        "coeff2", "coeff3", "coeff4", "coeff5", "coeff6", "coeff7","coeff8",
        "coeff9", "coeff10", "coeff11", "coeff12", "coeff13", "coeff14",
        "coeff15", "coeff16", "coeff17", "coeff18", "id")


    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ()
    ATTR_INT_PROTECTED_NAMES = ()

    # parameters considered are refined parameters
    ATTR_REF = (
        "coeff1", "coeff2", "coeff3", "coeff4", "coeff5", "coeff6", "coeff7",
        "coeff8", "coeff9", "coeff10", "coeff11", "coeff12", "coeff13",
        "coeff14", "coeff15", "coeff16", "coeff17", "coeff18")
    ATTR_SIGMA = tuple([f"{_h:}_sigma" for _h in ATTR_REF])
    ATTR_CONSTR_FLAG = tuple([f"{_h:}_constraint" for _h in ATTR_REF])
    ATTR_REF_FLAG = tuple([f"{_h:}_refinement" for _h in ATTR_REF])
    ATTR_CONSTR_MARK = tuple([f"{_h:}_mark" for _h in ATTR_REF])

    # formats if cif format
    D_FORMATS = {'time_max': "{:.2f}",
                 'coeff1': "{:.5f}", 'coeff2': "{:.5f}", 'coeff3': "{:.5f}",
                 'coeff4': "{:.5f}", 'coeff5': "{:.5f}", 'coeff6': "{:.5f}",
                 'coeff7': "{:.5f}", 'coeff8': "{:.5f}", 'coeff9': "{:.5f}",
                 'coeff10': "{:.5f}", 'coeff11': "{:.5f}", 'coeff12': "{:.5f}",
                 'coeff13': "{:.5f}", 'coeff14': "{:.5f}", 'coeff15': "{:.5f}",
                 'coeff16': "{:.5f}", 'coeff17': "{:.5f}", 'coeff18': "{:.5f}"}

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

    PREFIX = "tof_background"

    def __init__(self, **kwargs) -> NoReturn:
        super(TOFBackground, self).__init__()

        # defined for any integer and float parameters
        D_MIN = {"time_max": 0.01}

        # defined for ani integer and float parameters
        D_MAX = {}

        self.__dict__["D_MIN"] = D_MIN
        self.__dict__["D_MAX"] = D_MAX
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)

    def calc_background(self, time):
        """
        y = \sum_{i=1}^{18} coeff_i cos (pi i * (time/time_max))

        Parameters
        ----------
        time : TYPE
            DESCRIPTION.

        Returns
        -------
        y: background points

        """
        time_rel = (time/self.time_max)*numpy.pi
        np_cos = numpy.cos(0.*time_rel)
        res = self.coeff1 * np_cos
        if self.is_attribute("coeff2"):
            np_cos = numpy.cos(1.*time_rel)
            res += self.coeff2 * np_cos
        if self.is_attribute("coeff3"):
            np_cos = numpy.cos(2.*time_rel)
            res += self.coeff3 * np_cos
        if self.is_attribute("coeff4"):
            np_cos = numpy.cos(3.*time_rel)
            res += self.coeff4 * np_cos
        if self.is_attribute("coeff5"):
            np_cos = numpy.cos(4.*time_rel)
            res += self.coeff5 * np_cos
        if self.is_attribute("coeff6"):
            np_cos = numpy.cos(5.*time_rel)
            res += self.coeff6 * np_cos
        if self.is_attribute("coeff7"):
            np_cos = numpy.cos(6.*time_rel)
            res += self.coeff7 * np_cos
        if self.is_attribute("coeff8"):
            np_cos = numpy.cos(7.*time_rel)
            res += self.coeff8 * np_cos
        if self.is_attribute("coeff9"):
            np_cos = numpy.cos(8.*time_rel)
            res += self.coeff9 * np_cos
        if self.is_attribute("coeff10"):
            np_cos = numpy.cos(9.*time_rel)
            res += self.coeff10 * np_cos
        if self.is_attribute("coeff11"):
            np_cos = numpy.cos(10.*time_rel)
            res += self.coeff11 * np_cos
        if self.is_attribute("coeff12"):
            np_cos = numpy.cos(11.*time_rel)
            res += self.coeff12 * np_cos
        if self.is_attribute("coeff13"):
            np_cos = numpy.cos(12.*time_rel)
            res += self.coeff13 * np_cos
        if self.is_attribute("coeff14"):
            np_cos = numpy.cos(13.*time_rel)
            res += self.coeff14 * np_cos
        if self.is_attribute("coeff15"):
            np_cos = numpy.cos(14.*time_rel)
            res += self.coeff15 * np_cos
        if self.is_attribute("coeff16"):
            np_cos = numpy.cos(15.*time_rel)
            res += self.coeff16 * np_cos
        if self.is_attribute("coeff17"):
            np_cos = numpy.cos(16.*time_rel)
            res += self.coeff17 * np_cos
        if self.is_attribute("coeff18"):
            np_cos = numpy.cos(17.*time_rel)
            res += self.coeff18 * np_cos
        return res

    def get_coefficients(self):
        l_coeff = []
        last_number = 1
        for numb in range(1, 19):
            if self.is_attribute(f"coeff{numb:}"):
                coeff = getattr(self, f"coeff{numb:}")
                try:
                    coeff = float(coeff)
                    last_number = numb
                except TypeError:
                    coeff = 0
            else:
                coeff = 0.
            l_coeff.append(coeff)
        coefficients = numpy.array(l_coeff[:last_number], dtype=float)
        return coefficients

    def get_flags_coefficients(self):
        l_flag_coeff = []
        last_number = 1
        for numb in range(1, 19):
            if self.is_attribute(f"coeff{numb:}"):
                flag_coeff = getattr(self, f"coeff{numb:}_refinement")
                last_number = numb
            else:
                flag_coeff = False
            l_flag_coeff.append(flag_coeff)
        flags_coefficients = numpy.array(l_flag_coeff[:last_number], dtype=bool)
        return flags_coefficients

    def get_dictionary(self):
        res = {}
        res["background_time_max"] = getattr(self, "time_max")
        res["background_coefficients"] = self.get_coefficients()
        res["flags_background_coefficients"] = self.get_flags_coefficients()
        return res

class TOFBackgroundL(LoopN):
    """Bacground description for time-of-flight experiment.

    """
    ITEM_CLASS = TOFBackground
    ATTR_INDEX = "id"
    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(TOFBackgroundL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name


# s_cont = """
# _tof_background_time_max 30000
# _tof_background_coeff1 24832.850
# _tof_background_coeff2 6139.244
# _tof_background_coeff3 8063.472
# _tof_background_coeff4 3125.050
# _tof_background_coeff5 2566.956
# _tof_background_coeff6 311.077
# _tof_background_coeff7 837.348
# _tof_background_coeff8 -103.742
# _tof_background_coeff9 -11.806
# """

# obj = TOFBackground.from_cif(s_cont)
# print(obj, end="\n\n")
# time = numpy.linspace(1, 30000, 10)
# y = obj.calc_background(time)
# print("time bbackground")
# for h1, h2 in zip(time, y):
#     print(f"{h1:10.2f} {h2:10.2f}")
