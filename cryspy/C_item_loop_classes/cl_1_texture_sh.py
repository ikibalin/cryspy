"""Texture of the phases of a powder experiment by spherical harmonics.

TextureSh (loop texture_sh): one row for each textured phase: the maximal
order of the harmonics, the sample symmetry and the zeros of the sample
orientation. TextureShCoefficient (loop texture_sh_coefficient): the
coefficients C_l^mn (l, m: sample, n: crystal). SampleGoniometer (item
sample_goniometer): the geometry and the goniometer angles of the
measurement. The model and the conventions are those of GSAS-II, see
cryspy/A_functions_base/texture_spherical_harmonics.py.
"""
from typing import NoReturn
import numpy

from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary
from cryspy.A_functions_base.texture_spherical_harmonics import SAMPLE_SYMMETRIES

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN

DEG = numpy.pi / 180.
# omega, chi, phi, azimuth (degrees) without sample_goniometer: the cylinder
# axis of the texture is vertical (normal to the horizontal scattering plane)
DEFAULT_GONIOMETER_ANGLES = (0., 90., 0., 0.)


def form_d_default(attr_ref, d_default: dict) -> dict:
    """Default sigmas, flags and marks of the refined parameters."""
    for key in attr_ref:
        d_default[f"{key:}_sigma"] = 0.
        d_default[f"{key:}_constraint"] = False
        d_default[f"{key:}_refinement"] = False
        d_default[f"{key:}_mark"] = ""
    return d_default


class TextureSh(ItemN):
    """
    Texture of one phase by spherical harmonics (GSAS-II model).

    Attributes
    ----------
        - phase_label (mandatory): label of the phase of the experiment
        - order_max (mandatory): maximal even order L of the harmonics
        - sample_symmetry: 'cylindrical' (fibre texture, default), '2/m',
          'mmm' (rolling) or '-1' (none)
        - sample_omega, sample_chi, sample_phi: zeros of the sample
          orientation, degrees (refinable)
    """
    ATTR_MANDATORY_NAMES = ("phase_label", "order_max")
    ATTR_MANDATORY_TYPES = (str, int)
    ATTR_MANDATORY_CIF = ("phase_label", "order_max")

    ATTR_OPTIONAL_NAMES = ("sample_symmetry", "sample_omega", "sample_chi", "sample_phi")
    ATTR_OPTIONAL_TYPES = (str, float, float, float)
    ATTR_OPTIONAL_CIF = ("sample_symmetry", "sample_omega", "sample_chi", "sample_phi")

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ()
    ATTR_INT_PROTECTED_NAMES = ()

    ATTR_REF = ("sample_omega", "sample_chi", "sample_phi")
    ATTR_SIGMA = tuple([f"{_h:}_sigma" for _h in ATTR_REF])
    ATTR_CONSTR_FLAG = tuple([f"{_h:}_constraint" for _h in ATTR_REF])
    ATTR_REF_FLAG = tuple([f"{_h:}_refinement" for _h in ATTR_REF])
    ATTR_CONSTR_MARK = tuple([f"{_h:}_mark" for _h in ATTR_REF])

    D_FORMATS = {"sample_omega": "{:.3f}", "sample_chi": "{:.3f}", "sample_phi": "{:.3f}"}

    D_CONSTRAINTS = {"sample_symmetry": list(SAMPLE_SYMMETRIES)}

    D_DEFAULT = form_d_default(ATTR_REF, {"sample_symmetry": "cylindrical"})

    PREFIX = "texture_sh"

    D_REFERENCES = {"phase_label": "phase.label"}

    def __init__(self, **kwargs) -> NoReturn:
        super(TextureSh, self).__init__()
        self.__dict__["D_MIN"] = {"order_max": 2}
        self.__dict__["D_MAX"] = {"order_max": 34}
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)

    def get_sample_angles(self):
        """Zeros of the sample orientation omega, chi, phi in degrees."""
        return [float(getattr(self, name)) if self.is_attribute(name) else 0.
                for name in ("sample_omega", "sample_chi", "sample_phi")]


class TextureShL(LoopN):
    """Texture of the phases by spherical harmonics."""
    ITEM_CLASS = TextureSh
    ATTR_INDEX = "phase_label"

    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(TextureShL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name

    def get_dictionary(self):
        angles = numpy.array([item.get_sample_angles() for item in self.items],
                             dtype=float).reshape(-1, 3).transpose() * DEG
        flags = numpy.array(
            [[bool(getattr(item, f"{name:}_refinement")) for name in TextureSh.ATTR_REF]
             for item in self.items], dtype=bool).reshape(-1, 3).transpose()
        res = {
            "texture_sh_phase_name": numpy.array(
                [item.phase_label.lower() for item in self.items], dtype=str),
            "texture_sh_order_max": numpy.array(
                [int(item.order_max) for item in self.items], dtype=int),
            "texture_sh_sample_symmetry": numpy.array(
                [item.sample_symmetry if item.is_attribute("sample_symmetry") else "cylindrical"
                 for item in self.items], dtype=str),
            "texture_sh_sample_angles": angles,
            "flags_texture_sh_sample_angles": flags,
        }
        return res

    def take_parameters_from_dictionary(self, ddict_diffrn, l_parameter_name: list = None,
                                        l_sigma: list = None):
        keys = ddict_diffrn.keys()
        if "texture_sh_sample_angles" in keys:
            angles = ddict_diffrn["texture_sh_sample_angles"] / DEG
            for i_item, item in enumerate(self.items):
                for i_angle, name in enumerate(TextureSh.ATTR_REF):
                    setattr(item, name, float(angles[i_angle, i_item]))
        if l_parameter_name is not None:
            for name, sigma in zip(l_parameter_name, l_sigma):
                if name[0] == "texture_sh_sample_angles":
                    i_angle, i_item = name[1]
                    setattr(self.items[i_item], f"{TextureSh.ATTR_REF[i_angle]:}_sigma",
                            float(sigma) / DEG)


class TextureShCoefficient(ItemN):
    """
    Coefficient C_l^mn of the spherical harmonics texture of a phase.

    Attributes
    ----------
        - phase_label (mandatory): label of the phase
        - l (mandatory): even order of the harmonics
        - m (mandatory): index of the sample harmonic (0 for the
          cylindrical symmetry)
        - n (mandatory): index of the crystal harmonic (number of the
          harmonic for the cubic classes)
        - value (mandatory): the coefficient (refinable)
    """
    ATTR_MANDATORY_NAMES = ("phase_label", "l", "m", "n", "value")
    ATTR_MANDATORY_TYPES = (str, int, int, int, float)
    ATTR_MANDATORY_CIF = ("phase_label", "l", "m", "n", "value")

    ATTR_OPTIONAL_NAMES = ()
    ATTR_OPTIONAL_TYPES = ()
    ATTR_OPTIONAL_CIF = ()

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ()
    ATTR_INT_PROTECTED_NAMES = ()

    ATTR_REF = ("value", )
    ATTR_SIGMA = tuple([f"{_h:}_sigma" for _h in ATTR_REF])
    ATTR_CONSTR_FLAG = tuple([f"{_h:}_constraint" for _h in ATTR_REF])
    ATTR_REF_FLAG = tuple([f"{_h:}_refinement" for _h in ATTR_REF])
    ATTR_CONSTR_MARK = tuple([f"{_h:}_mark" for _h in ATTR_REF])

    D_FORMATS = {"value": "{:.5f}"}

    D_CONSTRAINTS = {}

    D_DEFAULT = form_d_default(ATTR_REF, {})

    PREFIX = "texture_sh_coefficient"

    D_REFERENCES = {"phase_label": "phase.label"}
    REQUIRED_CATEGORIES = ("texture_sh", )

    def __init__(self, **kwargs) -> NoReturn:
        super(TextureShCoefficient, self).__init__()
        self.__dict__["D_MIN"] = {"l": 2}
        self.__dict__["D_MAX"] = {"l": 34}
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)


class TextureShCoefficientL(LoopN):
    """Coefficients C_l^mn of the spherical harmonics texture."""
    ITEM_CLASS = TextureShCoefficient
    ATTR_INDEX = None

    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(TextureShCoefficientL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name

    def get_dictionary(self):
        res = {
            "texture_sh_coefficient_phase_name": numpy.array(
                [item.phase_label.lower() for item in self.items], dtype=str),
            "texture_sh_lmn": numpy.array(
                [[int(item.l), int(item.m), int(item.n)] for item in self.items],
                dtype=int).reshape(-1, 3).transpose(),
            "texture_sh_coefficient": numpy.array(
                [float(item.value) for item in self.items], dtype=float),
            "flags_texture_sh_coefficient": numpy.array(
                [bool(item.value_refinement) for item in self.items], dtype=bool),
        }
        return res

    def take_parameters_from_dictionary(self, ddict_diffrn, l_parameter_name: list = None,
                                        l_sigma: list = None):
        keys = ddict_diffrn.keys()
        if "texture_sh_coefficient" in keys:
            for item, value in zip(self.items, ddict_diffrn["texture_sh_coefficient"]):
                item.value = float(value)
        if l_parameter_name is not None:
            for name, sigma in zip(l_parameter_name, l_sigma):
                if name[0] == "texture_sh_coefficient":
                    self.items[name[1][0]].value_sigma = float(sigma)

    def take_coefficients(self, phase_label: str):
        """Items of the coefficients of the phase."""
        return [item for item in self.items
                if item.phase_label.lower() == phase_label.lower()]


class SampleGoniometer(ItemN):
    """
    Geometry and goniometer angles of the sample in a powder experiment.

    They give the direction of the scattering vector in the sample frame
    for the spherical harmonics texture (GSAS-II sample parameters).

    Attributes
    ----------
        - geometry: 'debye_scherrer' (default, transmission) or
          'bragg_brentano' (omega follows theta)
        - omega, chi, phi: goniometer angles of the sample, degrees
          (0, 90, 0 by default: the cylinder axis of the texture is
          vertical, normal to the horizontal scattering plane, as for a
          vertical capillary, a wire or the field of a 2D experiment)
        - azimuth: azimuth of the detector around the beam, degrees (0 for
          the scattering in the horizontal plane; added to the angle phi of
          the points of a 2D detector)
    """
    ATTR_MANDATORY_NAMES = ()
    ATTR_MANDATORY_TYPES = ()
    ATTR_MANDATORY_CIF = ()

    ATTR_OPTIONAL_NAMES = ("geometry", "omega", "chi", "phi", "azimuth")
    ATTR_OPTIONAL_TYPES = (str, float, float, float, float)
    ATTR_OPTIONAL_CIF = ("geometry", "omega", "chi", "phi", "azimuth")

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

    D_FORMATS = {"omega": "{:.3f}", "chi": "{:.3f}", "phi": "{:.3f}", "azimuth": "{:.3f}"}

    D_CONSTRAINTS = {"geometry": ["debye_scherrer", "bragg_brentano"]}

    D_DEFAULT = {"geometry": "debye_scherrer"}

    PREFIX = "sample_goniometer"

    def __init__(self, **kwargs) -> NoReturn:
        super(SampleGoniometer, self).__init__()
        self.__dict__["D_MIN"] = {}
        self.__dict__["D_MAX"] = {}
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)

    def get_dictionary(self):
        angles = [float(getattr(self, name)) if self.is_attribute(name) else default
                  for name, default in zip(("omega", "chi", "phi", "azimuth"),
                                           DEFAULT_GONIOMETER_ANGLES)]
        geometry = self.geometry if self.is_attribute("geometry") else "debye_scherrer"
        return {"goniometer_geometry": str(geometry),
                "goniometer_angles": numpy.array(angles, dtype=float) * DEG}
