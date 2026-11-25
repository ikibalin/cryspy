"""Random stacking faults of a powder phase (loop stacking_fault).

StackingFault: one row for each phase whose layers normal to
(normal_h normal_k normal_l) are shifted at random by the fault vector
(shift_x shift_y shift_z) with the probability `probability` (refinable)
between two neighbour layers. The reflections H with H.t not an integer
are spread along the normal of the layers (asymmetric or broad lines),
the others stay sharp. See cryspy/A_functions_base/stacking_faults.py.

loop_
_stacking_fault_phase_label
_stacking_fault_normal_h
_stacking_fault_normal_k
_stacking_fault_normal_l
_stacking_fault_shift_x
_stacking_fault_shift_y
_stacking_fault_shift_z
_stacking_fault_probability
Pbcn  0 0 2  0.0 0.5 0.0  0.19(1)
"""
from typing import NoReturn
import numpy

from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN


class StackingFault(ItemN):
    """
    Random stacking faults of the layers of a powder phase.

    Attributes
    ----------
        - phase_label (mandatory): label of the phase of the experiment
        - normal_h, normal_k, normal_l (mandatory): Miller indices of the
          layer planes; (0 0 2) are layers at c/2, two per cell
        - shift_x, shift_y, shift_z: fault vector, fractional coordinates
          of the cell (the extra shift of the layers above a fault)
        - probability: probability of a fault between two neighbour
          layers, 0 ... 1 (refinable)
    """
    ATTR_MANDATORY_NAMES = ("phase_label", "normal_h", "normal_k", "normal_l")
    ATTR_MANDATORY_TYPES = (str, int, int, int)
    ATTR_MANDATORY_CIF = ("phase_label", "normal_h", "normal_k", "normal_l")

    ATTR_OPTIONAL_NAMES = ("shift_x", "shift_y", "shift_z", "probability")
    ATTR_OPTIONAL_TYPES = (float, float, float, float)
    ATTR_OPTIONAL_CIF = ("shift_x", "shift_y", "shift_z", "probability")

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ()
    ATTR_INT_PROTECTED_NAMES = ()

    ATTR_REF = ("probability", )
    ATTR_SIGMA = tuple([f"{_h:}_sigma" for _h in ATTR_REF])
    ATTR_CONSTR_FLAG = tuple([f"{_h:}_constraint" for _h in ATTR_REF])
    ATTR_REF_FLAG = tuple([f"{_h:}_refinement" for _h in ATTR_REF])
    ATTR_CONSTR_MARK = tuple([f"{_h:}_mark" for _h in ATTR_REF])

    D_FORMATS = {"shift_x": "{:.5f}", "shift_y": "{:.5f}", "shift_z": "{:.5f}",
                 "probability": "{:.5f}"}

    D_CONSTRAINTS = {}

    D_DEFAULT = {"shift_x": 0., "shift_y": 0., "shift_z": 0., "probability": 0.1}
    for _h in ATTR_SIGMA:
        D_DEFAULT[_h] = 0.
    for _h in (ATTR_CONSTR_FLAG + ATTR_REF_FLAG):
        D_DEFAULT[_h] = False
    for _h in ATTR_CONSTR_MARK:
        D_DEFAULT[_h] = ""

    PREFIX = "stacking_fault"

    def __init__(self, **kwargs) -> NoReturn:
        super(StackingFault, self).__init__()
        self.__dict__["D_MIN"] = {"probability": 0.}
        self.__dict__["D_MAX"] = {"probability": 1.}
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)

    def get_normal_index(self):
        """Miller indices (h, k, l) of the layer planes."""
        return numpy.array([self.normal_h, self.normal_k, self.normal_l], dtype=int)

    def get_shift(self):
        """Fault vector t, fractional coordinates of the cell."""
        return numpy.array([self.shift_x, self.shift_y, self.shift_z], dtype=float)


class StackingFaultL(LoopN):
    """Random stacking faults of the layers of the powder phases."""
    ITEM_CLASS = StackingFault
    ATTR_INDEX = "phase_label"

    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(StackingFaultL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name

    def get_dictionary(self):
        n_item = len(self.items)
        return {
            "stacking_fault_phase_name": numpy.array(
                [item.phase_label.lower() for item in self.items], dtype=str),
            "stacking_fault_normal_index": numpy.array(
                [item.get_normal_index() for item in self.items],
                dtype=int).reshape(n_item, 3).transpose(),
            "stacking_fault_shift": numpy.array(
                [item.get_shift() for item in self.items],
                dtype=float).reshape(n_item, 3).transpose(),
            "stacking_fault_probability": numpy.array(
                [float(item.probability) for item in self.items], dtype=float),
            "flags_stacking_fault_probability": numpy.array(
                [bool(item.probability_refinement) for item in self.items], dtype=bool),
        }

    def take_parameters_from_dictionary(self, ddict_diffrn, l_parameter_name: list = None,
                                        l_sigma: list = None):
        key = "stacking_fault_probability"
        if key in ddict_diffrn.keys():
            for item, value in zip(self.items, ddict_diffrn[key]):
                item.probability = float(value)
        if l_parameter_name is not None:
            for name, sigma in zip(l_parameter_name, l_sigma):
                if name[0] == key:
                    self.items[name[1][0]].probability_sigma = float(sigma)

    def take_row(self, phase_label: str):
        """Row of the phase or None."""
        for item in self.items:
            if item.phase_label.lower() == phase_label.lower():
                return item
        return None
