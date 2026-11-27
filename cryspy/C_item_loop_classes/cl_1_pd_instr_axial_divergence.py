"""Asymmetry of the powder peaks from the axial divergence (Finger-Cox-Jephcoat).

_pd_instr_axial_divergence_s_l     0.020(1)
_pd_instr_axial_divergence_d_l     0.030(2)
_pd_instr_axial_divergence_ttheta_max 60.

See cryspy/A_functions_base/powder_diffraction_axial_divergence.py.
"""
from typing import NoReturn
import numpy

from cryspy.B_parent_classes.cl_1_item import ItemN


class PdInstrAxialDivergence(ItemN):
    """
    Asymmetry of the peaks of a 1d powder diffractometer from the axial divergence.

    The physical model of Finger, Cox and Jephcoat (J. Appl. Cryst. 27
    (1994) 892) of the asymmetry at low and high angles caused by the
    axial divergence of the beam on the curved Debye-Scherrer cone, with
    two parameters of the geometry instead of the empirical p1-p4 of
    pd_instr_reflex_asymmetry: below 90 degrees the peaks get a tail to
    low angles, above 90 degrees to high angles. The same parameters as
    S_L and D_L of FullProf and S/L and H/L of GSAS-II.

    Attributes
    ----------
        - s_l: half height of the sample over the distance sample-detector
          (refinable, dimensionless)
        - d_l: half height of the detector or of the receiving slit over the
          same distance (refinable, dimensionless)
        - ttheta_max (optional): reflections above this 2theta, degrees,
          keep the symmetric profile (FullProf AsymLim); 180 by default
    """
    ATTR_MANDATORY_NAMES = ("s_l", "d_l")
    ATTR_MANDATORY_TYPES = (float, float)
    ATTR_MANDATORY_CIF = ("s_l", "d_l")

    ATTR_OPTIONAL_NAMES = ("ttheta_max", )
    ATTR_OPTIONAL_TYPES = (float, )
    ATTR_OPTIONAL_CIF = ("ttheta_max", )

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ()
    ATTR_INT_PROTECTED_NAMES = ()

    ATTR_REF = ("s_l", "d_l")
    ATTR_SIGMA = tuple([f"{_h:}_sigma" for _h in ATTR_REF])
    ATTR_CONSTR_FLAG = tuple([f"{_h:}_constraint" for _h in ATTR_REF])
    ATTR_REF_FLAG = tuple([f"{_h:}_refinement" for _h in ATTR_REF])
    ATTR_CONSTR_MARK = tuple([f"{_h:}_mark" for _h in ATTR_REF])

    D_FORMATS = {"s_l": "{:.5f}", "d_l": "{:.5f}", "ttheta_max": "{:.2f}"}

    D_CONSTRAINTS = {}

    D_DEFAULT = {}
    for _h in ATTR_SIGMA:
        D_DEFAULT[_h] = 0.
    for _h in (ATTR_CONSTR_FLAG + ATTR_REF_FLAG):
        D_DEFAULT[_h] = False
    for _h in ATTR_CONSTR_MARK:
        D_DEFAULT[_h] = ""

    PREFIX = "pd_instr_axial_divergence"

    def __init__(self, **kwargs) -> NoReturn:
        super(PdInstrAxialDivergence, self).__init__()
        self.__dict__["D_MIN"] = {"s_l": 0., "d_l": 0., "ttheta_max": 0.}
        self.__dict__["D_MAX"] = {"ttheta_max": 180.}
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)

    def get_dictionary(self):
        ttheta_max = float(self.ttheta_max) if self.is_attribute("ttheta_max") else 180.
        return {
            "axial_divergence_parameters": numpy.array(
                [self.s_l, self.d_l], dtype=float),
            "flags_axial_divergence_parameters": numpy.array(
                [self.s_l_refinement, self.d_l_refinement], dtype=bool),
            "axial_divergence_ttheta_max": numpy.radians(ttheta_max),
        }

    def take_parameters_from_dictionary(self, ddict_diffrn, l_parameter_name: list = None,
                                        l_sigma: list = None):
        if "axial_divergence_parameters" in ddict_diffrn.keys():
            hh = ddict_diffrn["axial_divergence_parameters"]
            self.s_l = float(hh[0])
            self.d_l = float(hh[1])
        if l_parameter_name is not None:
            for name, sigma in zip(l_parameter_name, l_sigma):
                if name[0] == "axial_divergence_parameters":
                    attr = ("s_l_sigma", "d_l_sigma")[name[1][0]]
                    setattr(self, attr, float(sigma))
