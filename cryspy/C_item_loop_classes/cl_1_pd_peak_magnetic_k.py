"""Magnetic satellites of a powder pattern (loop pd_peak_magnetic_k).

Written after a calculation for every phase with propagation vectors
(atom_site_moment_Fourier); the loop name is the phase label. Equivalent
satellites are merged into one row (index_mult). The indices are the
fractional h + k (reciprocal lattice units); wave_vector_seq_id is the
number of the propagation vector (atom_site_Fourier_wave_vector).

loop_Ho2BaNiO5
_pd_peak_magnetic_k_index_h
_pd_peak_magnetic_k_index_k
_pd_peak_magnetic_k_index_l
_pd_peak_magnetic_k_wave_vector_seq_id
_pd_peak_magnetic_k_index_mult
_pd_peak_magnetic_k_2theta
_pd_peak_magnetic_k_intensity
 0.5 0.0 0.5 1 4 12.74 101.37
"""
from typing import NoReturn

from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN


class PdPeakMagneticK(ItemN):
    """One magnetic satellite of a powder pattern.

    Attributes
    ----------
        - index_h, index_k, index_l (mandatory): fractional indices h + k
        - wave_vector_seq_id: number of the propagation vector
        - index_mult: number of merged equivalent satellites
        - ttheta: position 2theta (degrees)
        - intensity: summed intensity with the scale and Lorentz factors
    """
    ATTR_MANDATORY_NAMES = ("index_h", "index_k", "index_l")
    ATTR_MANDATORY_TYPES = (float, float, float)
    ATTR_MANDATORY_CIF = ("index_h", "index_k", "index_l")

    ATTR_OPTIONAL_NAMES = (
        "wave_vector_seq_id", "index_mult", "ttheta", "intensity")
    ATTR_OPTIONAL_TYPES = (int, int, float, float)
    ATTR_OPTIONAL_CIF = (
        "wave_vector_seq_id", "index_mult", "2theta", "intensity")

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ()
    ATTR_INT_PROTECTED_NAMES = ()

    ATTR_REF = ()
    ATTR_SIGMA = tuple([f"{_h:}_sigma" for _h in ATTR_REF])
    ATTR_CONSTR_FLAG = tuple([f"{_h:}_constraint" for _h in ATTR_REF])
    ATTR_REF_FLAG = tuple([f"{_h:}_refinement" for _h in ATTR_REF])
    ATTR_CONSTR_MARK = tuple([f"{_h:}_mark" for _h in ATTR_REF])

    D_FORMATS = {
        "index_h": "{:.4f}", "index_k": "{:.4f}", "index_l": "{:.4f}",
        "ttheta": "{:.3f}", "intensity": "{:.2f}"}

    D_CONSTRAINTS = {}
    D_DEFAULT = {}

    PREFIX = "pd_peak_magnetic_k"

    def __init__(self, **kwargs) -> NoReturn:
        super(PdPeakMagneticK, self).__init__()

        self.__dict__["D_MIN"] = {}
        self.__dict__["D_MAX"] = {}
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)


class PdPeakMagneticKL(LoopN):
    """Magnetic satellites of a powder pattern, one loop per phase."""
    ITEM_CLASS = PdPeakMagneticK
    ATTR_INDEX = None

    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(PdPeakMagneticKL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name
