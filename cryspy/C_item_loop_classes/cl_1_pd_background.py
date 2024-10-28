import numpy
from typing import NoReturn
import scipy

from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN

from cryspy.C_item_loop_classes.cl_1_pd_meas import PdMeasL

na = numpy.newaxis

class PdBackground(ItemN):
    """Background point in powder 1d experiment.
    """
    ATTR_MANDATORY_NAMES = ("ttheta", "intensity")
    ATTR_MANDATORY_TYPES = (float, float)
    ATTR_MANDATORY_CIF = ("2theta", "intensity")

    ATTR_OPTIONAL_NAMES = ()
    ATTR_OPTIONAL_TYPES = ()
    ATTR_OPTIONAL_CIF = ()

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ()
    ATTR_INT_PROTECTED_NAMES = ()

    # parameters considered are refined parameters
    ATTR_REF = ("intensity", )
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

    PREFIX = "pd_background"

    def __init__(self, **kwargs) -> NoReturn:
        super(PdBackground, self).__init__()

        # defined for any integer and float parameters
        D_MIN = {"ttheta": 0}

        # defined for ani integer and float parameters
        D_MAX = {}

        self.__dict__["D_MIN"] = D_MIN
        self.__dict__["D_MAX"] = D_MAX
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)


class PdBackgroundL(LoopN):
    """Background points in powder 1d experiment.
    """
    ITEM_CLASS = PdBackground
    ATTR_INDEX = "ttheta"
    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(PdBackgroundL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name

    def interpolate_by_points(self, tth):
        l_ttheta = self.ttheta
        l_intensity = self.intensity
        tth_b = numpy.array(l_ttheta, dtype=float)
        int_b = numpy.array(l_intensity, dtype=float)
        if len(l_ttheta) == 0:
            int_1d = numpy.zeros(tth.size, dtype=float)
        else:
            int_1d = numpy.interp(tth, tth_b, int_b)
        return int_1d

    def define_points(self, pd_meas: PdMeasL, step_ttheta: float = 10.):
        ttheta = numpy.array(pd_meas.ttheta, dtype=float)
        if pd_meas.items[0].is_attribute("intensity_plus"):
            intensity = (
                numpy.array(pd_meas.intensity_plus, dtype=float) +
                numpy.array(pd_meas.intensity_minus, dtype=float))
            intensity_sigma = numpy.sqrt(
                numpy.square(numpy.array(pd_meas.intensity_plus_sigma, dtype=float)) +
                numpy.square(numpy.array(pd_meas.intensity_minus_sigma, dtype=float)))
        else:
            intensity = numpy.array(pd_meas.intensity, dtype=float)
            intensity_sigma = numpy.array(pd_meas.intensity_sigma, dtype=float)

        step_n = int(step_ttheta/(ttheta[1]-ttheta[0]))
        if step_n == 0:
            step_n = 1
        intensity_bkg = estimate_background(intensity, step_n=step_n)

        delta_0 = 0.
        res = scipy.optimize.minimize(
            lambda p: calc_chi_sq(intensity_bkg, intensity, intensity_sigma, p),
            delta_0)
        delta_0 = res["x"]
        intensity_bkg += delta_0

        # n_points = int((ttheta.max()-ttheta.min())/step_ttheta + 2)
        # ttheta_bkgr = numpy.linspace(ttheta.min(), ttheta.max(), n_points, endpoint=True)
        #
        # flags = numpy.abs(ttheta[na, :]-ttheta_bkgr[:, na])<step_ttheta
        #
        # int_bkrg = numpy.array([numpy.nanmin(intensity_bkg[flag], axis=0) for flag in flags], dtype=float)
        # int_bkrg = numpy.round(int_bkrg, decimals=5)
        self.numpy_ttheta = ttheta
        self.numpy_intensity = intensity_bkg
        self.items = []
        self.numpy_to_items()
        return

    def get_dictionary(self):
        res = {}
        res["background_ttheta"] = numpy.array(self.ttheta, dtype=float)*numpy.pi/180.
        res["background_intensity"] = numpy.array(self.intensity, dtype=float)
        res["flags_background_intensity"] = numpy.array(self.intensity_refinement, dtype=bool)
        return res

def estimate_background(y_exp, step_n: int = 10):
    n_points = numpy.arange(1, step_n+1)
    y_aver = y_exp.sum()/y_exp.size
    y_min = y_exp.min()
    flag = y_exp < y_aver + 2*(y_aver-y_min)

    # y_bkg_es = numpy.zeros_like(y_exp)
    # y_bkg_es[flag] = y_exp[flag]
    # y_bkg_es[numpy.logical_not(flag)] = y_aver+2*(y_aver-y_min)
    y_bkg_es = numpy.copy(y_exp)

    n_total = y_bkg_es.size
    for iii in range(100):
        y_bkg_new = numpy.zeros_like(y_bkg_es)
        for i_point in range(n_total):
            ind_left = i_point - n_points
            ind_right = i_point + n_points

            flag_left = ind_left < 0
            if numpy.any(flag_left):
                ind_left[flag_left] = 0
                ind_right[flag_left] = 0

            flag_right = ind_right >= n_total
            if numpy.any(flag_right):
                ind_right[flag_right] = n_total-1
                ind_left[flag_right] = n_total-1

            y_left = y_bkg_es[ind_left]
            y_right = y_bkg_es[ind_right]

            y_bkg_new[i_point] = (y_left+y_right).sum()/(2*step_n)


        if iii != 99:
            flag = y_bkg_new > y_bkg_es
            y_bkg_new[flag] = y_bkg_es[flag]
        y_bkg_es = y_bkg_new
    y_bkg = y_bkg_es
    return y_bkg


def calc_chi_sq(y_bkg, y_exp, sy_exp, delta):
    y_diff = (y_exp-(y_bkg+delta))/sy_exp
    flag = y_diff > 1
    y_diff[flag] = 1.
    wss = numpy.sum(numpy.square(y_diff))
    return wss