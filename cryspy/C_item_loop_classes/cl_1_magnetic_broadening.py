"""Extra broadening of the magnetic reflections of the phases of a powder experiment.

MagneticBroadening (loop magnetic_broadening): one row for each phase
whose magnetic reflections are broader (or narrower) than the nuclear
ones. The magnetic intensity (|F_M|^2, the paramagnetic terms) and the
nuclear-magnetic interference of the phase are convolved with an extra
Lorentzian, whose half width in Q is `width` (inverse angstrom,
1/xi_mag - 1/xi_nucl); the nuclear intensity |F_N|^2 keeps the profile of
the phase. width = 0 gives the magnetic reflections the nuclear profile.
See cryspy/A_functions_base/magnetic_broadening.py.

loop_
_magnetic_broadening_phase_label
_magnetic_broadening_width
Fe3O4  0.0125(4)

or, when the width depends on the direction of hkl (inverse angstrom, U_ij
convention; width_11, width_22, width_33 are the widths along a*, b*, c*):

loop_
_magnetic_broadening_phase_label
_magnetic_broadening_width_11
_magnetic_broadening_width_22
_magnetic_broadening_width_33
_magnetic_broadening_width_12
_magnetic_broadening_width_13
_magnetic_broadening_width_23
TbTe3  0.004(1) 0.004(1) 0.031(2) 0. 0. 0.
"""
from typing import NoReturn
import numpy

from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN

ANISO_NAMES = ("width_11", "width_22", "width_33", "width_12", "width_13", "width_23")


class MagneticBroadening(ItemN):
    """
    Extra broadening of the magnetic reflections of one phase of a powder experiment.

    Attributes
    ----------
        - phase_label (mandatory): label of the phase of the experiment
        - width: extra half width of the magnetic reflections in Q, inverse
          angstrom, 1/xi_mag - 1/xi_nucl; 0 (default) for the profile of
          the nuclear reflections, > 0 for shorter magnetic correlations,
          < 0 for longer ones (refinable)
        - width_11, width_22, width_33, width_12, width_13, width_23:
          width depending on the direction of hkl, inverse angstrom, in the
          convention of U_ij (refinable); when given, width is not used
    """
    ATTR_MANDATORY_NAMES = ("phase_label", )
    ATTR_MANDATORY_TYPES = (str, )
    ATTR_MANDATORY_CIF = ("phase_label", )

    ATTR_OPTIONAL_NAMES = ("width", ) + ANISO_NAMES
    ATTR_OPTIONAL_TYPES = (float, ) * 7
    ATTR_OPTIONAL_CIF = ("width", ) + ANISO_NAMES

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ()
    ATTR_INT_PROTECTED_NAMES = ()

    ATTR_REF = ("width", ) + ANISO_NAMES
    ATTR_SIGMA = tuple([f"{_h:}_sigma" for _h in ATTR_REF])
    ATTR_CONSTR_FLAG = tuple([f"{_h:}_constraint" for _h in ATTR_REF])
    ATTR_REF_FLAG = tuple([f"{_h:}_refinement" for _h in ATTR_REF])
    ATTR_CONSTR_MARK = tuple([f"{_h:}_mark" for _h in ATTR_REF])

    D_FORMATS = {}
    for _h in ATTR_REF:
        D_FORMATS[_h] = "{:.5f}"

    D_CONSTRAINTS = {}

    D_DEFAULT = {}
    for _h in ATTR_SIGMA:
        D_DEFAULT[_h] = 0.
    for _h in (ATTR_CONSTR_FLAG + ATTR_REF_FLAG):
        D_DEFAULT[_h] = False
    for _h in ATTR_CONSTR_MARK:
        D_DEFAULT[_h] = ""

    PREFIX = "magnetic_broadening"

    def __init__(self, **kwargs) -> NoReturn:
        super(MagneticBroadening, self).__init__()
        self.__dict__["D_MIN"] = {}
        self.__dict__["D_MAX"] = {}
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)

    def is_anisotropic(self) -> bool:
        """True when the width depends on the direction of hkl (width_ij given)."""
        return any(self.is_attribute(name) for name in ANISO_NAMES)

    def get_aniso(self):
        """Components W_11, W_22, W_33, W_12, W_13, W_23 (inverse angstrom)."""
        return [float(getattr(self, name)) if self.is_attribute(name) else 0.
                for name in ANISO_NAMES]


class MagneticBroadeningL(LoopN):
    """Extra broadening of the magnetic reflections of the phases of a powder experiment."""
    ITEM_CLASS = MagneticBroadening
    ATTR_INDEX = "phase_label"

    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(MagneticBroadeningL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name

    def get_dictionary(self):
        n_item = len(self.items)
        width = numpy.array(
            [float(item.width) if item.is_attribute("width") else 0.
             for item in self.items], dtype=float)
        flags_width = numpy.array(
            [bool(item.width_refinement) and not item.is_anisotropic()
             for item in self.items], dtype=bool)
        aniso = numpy.array([item.get_aniso() for item in self.items],
                            dtype=float).reshape(n_item, 6).transpose()
        flags_aniso = numpy.array(
            [[bool(getattr(item, f"{name:}_refinement")) for name in ANISO_NAMES]
             for item in self.items], dtype=bool).reshape(n_item, 6).transpose()
        res = {
            "magnetic_broadening_phase_name": numpy.array(
                [item.phase_label.lower() for item in self.items], dtype=str),
            "magnetic_broadening_is_aniso": numpy.array(
                [item.is_anisotropic() for item in self.items], dtype=bool),
            "magnetic_broadening_width": width,
            "flags_magnetic_broadening_width": flags_width,
            "magnetic_broadening_aniso": aniso,
            "flags_magnetic_broadening_aniso": flags_aniso,
        }
        return res

    def take_parameters_from_dictionary(self, ddict_diffrn, l_parameter_name: list = None,
                                        l_sigma: list = None):
        keys = ddict_diffrn.keys()
        if "magnetic_broadening_width" in keys:
            for item, value in zip(self.items, ddict_diffrn["magnetic_broadening_width"]):
                if item.is_attribute("width") and not item.is_anisotropic():
                    item.width = float(value)
        # the symmetrized tensor (dict_in_out of the calculation) is written
        # over the refined one: its dependent components follow the others
        for key in ("magnetic_broadening_aniso", "magnetic_broadening_aniso_symmetrized"):
            if key in keys:
                aniso = ddict_diffrn[key]
                for i_item, item in enumerate(self.items):
                    if not item.is_anisotropic():
                        continue
                    for i_comp, name in enumerate(ANISO_NAMES):
                        setattr(item, name, float(aniso[i_comp, i_item]))
        if l_parameter_name is not None:
            for name, sigma in zip(l_parameter_name, l_sigma):
                if name[0] == "magnetic_broadening_width":
                    self.items[name[1][0]].width_sigma = float(sigma)
                elif name[0] == "magnetic_broadening_aniso":
                    i_comp, i_item = name[1]
                    setattr(self.items[i_item], f"{ANISO_NAMES[i_comp]:}_sigma",
                            float(sigma))

    def take_row(self, phase_label: str):
        """Row of the phase or None."""
        for item in self.items:
            if item.phase_label.lower() == phase_label.lower():
                return item
        return None
