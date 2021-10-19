import numpy
from typing import NoReturn

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
            intensity = 0.5*(
                numpy.array(pd_meas.intensity_plus, dtype=float) + 
                numpy.array(pd_meas.intensity_minus, dtype=float))
        else:
            intensity = numpy.array(pd_meas.intensity, dtype=float)

        n_points = int((ttheta.max()-ttheta.min())/step_ttheta + 2)
        ttheta_bkgr = numpy.linspace(ttheta.min(), ttheta.max(), n_points, endpoint=True)
        
        flags = numpy.abs(ttheta[na, :]-ttheta_bkgr[:, na])<step_ttheta
        
        int_bkrg = numpy.array([numpy.nanmin(intensity[flag], axis=0) for flag in flags], dtype=float)
        int_bkrg = numpy.round(int_bkrg,5)
        self.numpy_ttheta = ttheta_bkgr
        self.numpy_intensity = int_bkrg
        self.items = []
        self.numpy_to_items()
        return 