"""Magnetic correlation length of the phases of a powder experiment.

MagneticCorrelation (loop magnetic_correlation): one row for each phase
whose magnetic reflections are broader than the nuclear ones. The magnetic
intensity (|F_M|^2, the paramagnetic terms) and the nuclear-magnetic
interference of the phase are convolved with an extra Lorentzian of the
half width 1/xi in Q; the nuclear intensity |F_N|^2 keeps the profile of
the phase. See cryspy/A_functions_base/magnetic_correlation.py.

loop_
_magnetic_correlation_phase_label
_magnetic_correlation_length
Fe3O4  150.(3)

or, for an ellipsoid of correlation lengths (angstrom^2, U_ij convention;
aniso_11, aniso_22, aniso_33 are the squared lengths along a*, b*, c*):

loop_
_magnetic_correlation_phase_label
_magnetic_correlation_aniso_11
_magnetic_correlation_aniso_22
_magnetic_correlation_aniso_33
_magnetic_correlation_aniso_12
_magnetic_correlation_aniso_13
_magnetic_correlation_aniso_23
TbTe3  40000.(500) 40000.(500) 900.(40) 0. 0. 0.
"""
from typing import NoReturn
import numpy

from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN

ANISO_NAMES = ("aniso_11", "aniso_22", "aniso_33", "aniso_12", "aniso_13", "aniso_23")


class MagneticCorrelation(ItemN):
    """
    Magnetic correlation length of one phase of a powder experiment.

    Attributes
    ----------
        - phase_label (mandatory): label of the phase of the experiment
        - length: isotropic correlation length xi, angstrom (the half width
          of the magnetic Lorentzian in Q is 1/xi; the Scherrer size of
          FullProf is 2 xi); 0 means infinite (refinable)
        - aniso_11, aniso_22, aniso_33, aniso_12, aniso_13, aniso_23:
          ellipsoid of correlation lengths, angstrom^2, in the convention
          of U_ij (refinable); when given, the length is not used
    """
    ATTR_MANDATORY_NAMES = ("phase_label", )
    ATTR_MANDATORY_TYPES = (str, )
    ATTR_MANDATORY_CIF = ("phase_label", )

    ATTR_OPTIONAL_NAMES = ("length", ) + ANISO_NAMES
    ATTR_OPTIONAL_TYPES = (float, ) * 7
    ATTR_OPTIONAL_CIF = ("length", ) + ANISO_NAMES

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ()
    ATTR_INT_PROTECTED_NAMES = ()

    ATTR_REF = ("length", ) + ANISO_NAMES
    ATTR_SIGMA = tuple([f"{_h:}_sigma" for _h in ATTR_REF])
    ATTR_CONSTR_FLAG = tuple([f"{_h:}_constraint" for _h in ATTR_REF])
    ATTR_REF_FLAG = tuple([f"{_h:}_refinement" for _h in ATTR_REF])
    ATTR_CONSTR_MARK = tuple([f"{_h:}_mark" for _h in ATTR_REF])

    D_FORMATS = {"length": "{:.2f}"}
    for _h in ANISO_NAMES:
        D_FORMATS[_h] = "{:.2f}"

    D_CONSTRAINTS = {}

    D_DEFAULT = {}
    for _h in ATTR_SIGMA:
        D_DEFAULT[_h] = 0.
    for _h in (ATTR_CONSTR_FLAG + ATTR_REF_FLAG):
        D_DEFAULT[_h] = False
    for _h in ATTR_CONSTR_MARK:
        D_DEFAULT[_h] = ""

    PREFIX = "magnetic_correlation"

    def __init__(self, **kwargs) -> NoReturn:
        super(MagneticCorrelation, self).__init__()
        self.__dict__["D_MIN"] = {"length": 0.}
        self.__dict__["D_MAX"] = {}
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)

    def is_anisotropic(self) -> bool:
        """True when the ellipsoid of correlation lengths is given."""
        return any(self.is_attribute(name) for name in ANISO_NAMES)

    def get_aniso(self):
        """Components Xi_11, Xi_22, Xi_33, Xi_12, Xi_13, Xi_23 (angstrom^2)."""
        return [float(getattr(self, name)) if self.is_attribute(name) else 0.
                for name in ANISO_NAMES]


class MagneticCorrelationL(LoopN):
    """Magnetic correlation lengths of the phases of a powder experiment."""
    ITEM_CLASS = MagneticCorrelation
    ATTR_INDEX = "phase_label"

    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(MagneticCorrelationL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name

    def get_dictionary(self):
        n_item = len(self.items)
        length = numpy.array(
            [float(item.length) if item.is_attribute("length") else 0.
             for item in self.items], dtype=float)
        flags_length = numpy.array(
            [bool(item.length_refinement) for item in self.items], dtype=bool)
        aniso = numpy.array([item.get_aniso() for item in self.items],
                            dtype=float).reshape(n_item, 6).transpose()
        flags_aniso = numpy.array(
            [[bool(getattr(item, f"{name:}_refinement")) for name in ANISO_NAMES]
             for item in self.items], dtype=bool).reshape(n_item, 6).transpose()
        res = {
            "magnetic_correlation_phase_name": numpy.array(
                [item.phase_label.lower() for item in self.items], dtype=str),
            "magnetic_correlation_is_aniso": numpy.array(
                [item.is_anisotropic() for item in self.items], dtype=bool),
            "magnetic_correlation_length": length,
            "flags_magnetic_correlation_length": flags_length,
            "magnetic_correlation_aniso": aniso,
            "flags_magnetic_correlation_aniso": flags_aniso,
        }
        return res

    def take_parameters_from_dictionary(self, ddict_diffrn, l_parameter_name: list = None,
                                        l_sigma: list = None):
        keys = ddict_diffrn.keys()
        if "magnetic_correlation_length" in keys:
            for item, value in zip(self.items, ddict_diffrn["magnetic_correlation_length"]):
                if item.is_attribute("length"):
                    item.length = float(value)
        # the symmetrized tensor (dict_in_out of the calculation) is written
        # over the refined one: its dependent components follow the others
        for key in ("magnetic_correlation_aniso", "magnetic_correlation_aniso_symmetrized"):
            if key in keys:
                aniso = ddict_diffrn[key]
                for i_item, item in enumerate(self.items):
                    if not item.is_anisotropic():
                        continue
                    for i_comp, name in enumerate(ANISO_NAMES):
                        setattr(item, name, float(aniso[i_comp, i_item]))
        if l_parameter_name is not None:
            for name, sigma in zip(l_parameter_name, l_sigma):
                if name[0] == "magnetic_correlation_length":
                    self.items[name[1][0]].length_sigma = float(sigma)
                elif name[0] == "magnetic_correlation_aniso":
                    i_comp, i_item = name[1]
                    setattr(self.items[i_item], f"{ANISO_NAMES[i_comp]:}_sigma",
                            float(sigma))

    def take_row(self, phase_label: str):
        """Row of the phase or None."""
        for item in self.items:
            if item.phase_label.lower() == phase_label.lower():
                return item
        return None
