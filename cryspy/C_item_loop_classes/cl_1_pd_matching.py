"""Matching mode of a powder pattern (item pd_matching).

In the matching mode (profile matching, Le Bail/Pawley fit) the integrated
intensities of the reflections are not calculated from the crystal
structure:

- ``refine_intensities``: the integrated intensities of all reflections of
  every phase are fitted to the pattern (linear least squares with
  intensities >= 0) at every calculation, for the current peak positions
  and shapes. The cell, the zero shift, the resolution, the asymmetry and
  the background can be refined at the same time. The fitted intensities
  are written to the loop pd_peak of every phase;
- ``fixed_intensities``: the integrated intensities are taken from the loop
  pd_peak of every phase and kept fixed; the peak positions and shapes (the
  cell, the zero shift, the resolution, ...) are refined;
- ``off`` (the same as no item): the Rietveld calculation.

_pd_matching_mode refine_intensities
"""
from typing import NoReturn

from cryspy.B_parent_classes.cl_1_item import ItemN

MATCHING_MODES = ("off", "refine_intensities", "fixed_intensities")


class PdMatching(ItemN):
    """Matching mode of a powder pattern.

    Attributes
    ----------
        - mode (mandatory): "off", "refine_intensities" (fit the integrated
          intensities of the reflections, >= 0, and write them to pd_peak)
          or "fixed_intensities" (keep the intensities of pd_peak and refine
          the peak positions and shapes)
    """
    ATTR_MANDATORY_NAMES = ("mode", )
    ATTR_MANDATORY_TYPES = (str, )
    ATTR_MANDATORY_CIF = ("mode", )

    ATTR_OPTIONAL_NAMES = ()
    ATTR_OPTIONAL_TYPES = ()
    ATTR_OPTIONAL_CIF = ()

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ()
    ATTR_INT_PROTECTED_NAMES = ()

    ATTR_REF = ()
    ATTR_SIGMA = ()
    ATTR_CONSTR_FLAG = ()
    ATTR_REF_FLAG = ()
    ATTR_CONSTR_MARK = ()

    D_CONSTRAINTS = {"mode": list(MATCHING_MODES)}
    D_DEFAULT = {}

    PREFIX = "pd_matching"

    def __init__(self, **kwargs) -> NoReturn:
        super(PdMatching, self).__init__()

        self.__dict__["D_MIN"] = {}
        self.__dict__["D_MAX"] = {}
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)

    def get_dictionary(self):
        res = {}
        if self.is_attribute("mode"):
            res["matching_mode"] = self.mode
        return res

    def report(self) -> str:
        """What the matching mode does (markdown); '' when it is off."""
        if not self.is_attribute("mode") or self.mode == "off":
            return ""
        if self.mode == "refine_intensities":
            s_mode = ("the integrated intensities of the reflections are fitted "
                      "to the pattern (>= 0) and written to pd_peak")
        else:
            s_mode = ("the integrated intensities of the reflections are taken "
                      "from pd_peak and kept fixed")
        return (f"# Matching mode\n\nMode `{self.mode:}`: {s_mode:}; the atoms "
                "and the scale factor of the phases are not used.")
