"""Layers without correlations between them in a powder experiment (Warren line shape).

LayerCorrelation (loop layer_correlation): one row for each phase whose
order (magnetic, nuclear or both) exists only inside the layers normal to
(normal_h normal_k normal_l): the layers are not correlated along the
normal and scatter along rods (hk) instead of Bragg reflections, which
gives the asymmetric saw-tooth line of Warren (sharp front, long tail to
high angles). The order of the phase given by `contribution` is taken away
from its Bragg reflections and calculated as rods; the in-plane
correlation length `length` (angstrom) and the `scale` are refinable,
the moments come from the crystal (refinable there).
See cryspy/A_functions_base/layer_correlations.py.

loop_
_layer_correlation_phase_label
_layer_correlation_normal_h
_layer_correlation_normal_k
_layer_correlation_normal_l
_layer_correlation_contribution
_layer_correlation_length
_layer_correlation_scale
TbTe3  0 1 0  magnetic  85.(4)  0.83(2)
"""
from typing import NoReturn
import numpy

from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN

CONTRIBUTIONS = ("magnetic", "nuclear", "nuclear and magnetic")


class LayerCorrelation(ItemN):
    """
    Layers of a phase without correlations between them (Warren line shape).

    Attributes
    ----------
        - phase_label (mandatory): label of the phase of the experiment
        - normal_h, normal_k, normal_l (mandatory): Miller indices of the
          layer planes (the normal of the layers is the reciprocal vector
          h a* + k b* + l c*)
        - contribution: "magnetic" (default), "nuclear" or
          "nuclear and magnetic": the order that exists only in the layers
        - length: correlation length in the plane of the layers, angstrom
          (Warren's mean layer dimension, refinable)
        - scale: scale factor of the rods, the meaning of the scale of the
          phase (refinable)
    """
    ATTR_MANDATORY_NAMES = ("phase_label", "normal_h", "normal_k", "normal_l")
    ATTR_MANDATORY_TYPES = (str, int, int, int)
    ATTR_MANDATORY_CIF = ("phase_label", "normal_h", "normal_k", "normal_l")

    ATTR_OPTIONAL_NAMES = ("contribution", "length", "scale")
    ATTR_OPTIONAL_TYPES = (str, float, float)
    ATTR_OPTIONAL_CIF = ("contribution", "length", "scale")

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ()
    ATTR_INT_PROTECTED_NAMES = ()

    ATTR_REF = ("length", "scale")
    ATTR_SIGMA = tuple([f"{_h:}_sigma" for _h in ATTR_REF])
    ATTR_CONSTR_FLAG = tuple([f"{_h:}_constraint" for _h in ATTR_REF])
    ATTR_REF_FLAG = tuple([f"{_h:}_refinement" for _h in ATTR_REF])
    ATTR_CONSTR_MARK = tuple([f"{_h:}_mark" for _h in ATTR_REF])

    D_FORMATS = {"length": "{:.2f}", "scale": "{:.5f}"}

    D_CONSTRAINTS = {"contribution": CONTRIBUTIONS}

    D_DEFAULT = {"contribution": "magnetic", "length": 100., "scale": 1.}
    for _h in ATTR_SIGMA:
        D_DEFAULT[_h] = 0.
    for _h in (ATTR_CONSTR_FLAG + ATTR_REF_FLAG):
        D_DEFAULT[_h] = False
    for _h in ATTR_CONSTR_MARK:
        D_DEFAULT[_h] = ""

    PREFIX = "layer_correlation"

    def __init__(self, **kwargs) -> NoReturn:
        super(LayerCorrelation, self).__init__()
        self.__dict__["D_MIN"] = {"length": 1.}
        self.__dict__["D_MAX"] = {}
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)

    def get_normal_index(self):
        """Miller indices (h, k, l) of the layer planes."""
        return numpy.array([self.normal_h, self.normal_k, self.normal_l], dtype=int)


class LayerCorrelationL(LoopN):
    """Layers of the phases without correlations between them (Warren line shape)."""
    ITEM_CLASS = LayerCorrelation
    ATTR_INDEX = "phase_label"

    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(LayerCorrelationL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name

    def get_dictionary(self):
        n_item = len(self.items)
        res = {
            "layer_correlation_phase_name": numpy.array(
                [item.phase_label.lower() for item in self.items], dtype=str),
            "layer_correlation_normal_index": numpy.array(
                [item.get_normal_index() for item in self.items],
                dtype=int).reshape(n_item, 3).transpose(),
            "layer_correlation_contribution": numpy.array(
                [item.contribution for item in self.items], dtype=str),
            "layer_correlation_length": numpy.array(
                [float(item.length) for item in self.items], dtype=float),
            "flags_layer_correlation_length": numpy.array(
                [bool(item.length_refinement) for item in self.items], dtype=bool),
            "layer_correlation_scale": numpy.array(
                [float(item.scale) for item in self.items], dtype=float),
            "flags_layer_correlation_scale": numpy.array(
                [bool(item.scale_refinement) for item in self.items], dtype=bool),
        }
        return res

    def take_parameters_from_dictionary(self, ddict_diffrn, l_parameter_name: list = None,
                                        l_sigma: list = None):
        keys = ddict_diffrn.keys()
        for name in ("length", "scale"):
            key = f"layer_correlation_{name:}"
            if key in keys:
                for item, value in zip(self.items, ddict_diffrn[key]):
                    setattr(item, name, float(value))
        if l_parameter_name is not None:
            for name, sigma in zip(l_parameter_name, l_sigma):
                if name[0] in ("layer_correlation_length", "layer_correlation_scale"):
                    setattr(self.items[name[1][0]], f"{name[0][18:]:}_sigma", float(sigma))

    def take_row(self, phase_label: str):
        """Row of the phase or None."""
        for item in self.items:
            if item.phase_label.lower() == phase_label.lower():
                return item
        return None
