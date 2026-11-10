"""ExptlAbsorpt and ExptlCrystal: absorption correction of single crystal data.

The names follow the categories EXPTL_ABSORPT and EXPTL_CRYSTAL of the core
CIF dictionary (_exptl_absorpt_coefficient_mu, _exptl_absorpt_correction_type,
_exptl_absorpt_correction_T_min, _exptl_crystal_size_max, ...). CrysPy adds
the orientation of the ellipsoid (_exptl_crystal_ellipsoid_euler_*), the side
of the detector (_exptl_absorpt_diffracted_beam) and mu at half the wavelength
(_exptl_absorpt_coefficient_mu_half) for the lambda/2 contamination.
"""
from typing import NoReturn
import numpy

from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary
from cryspy.A_functions_base.absorption import (
    calc_euler_zyz_by_rotation_matrix, calc_rotation_matrix_by_euler_zyz)

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN


class ExptlAbsorpt(ItemN):
    """
    Absorption of the neutrons in the crystal and its correction.

    Attributes
    ----------
        - correction_type (mandatory): 'none' (no correction), 'sphere'
          (spherical crystal of the radius _exptl_crystal_size_rad) or
          'gaussian' (ellipsoid of the sizes _exptl_crystal_size_max, _mid,
          _min, numerical Gaussian integration over its volume)
        - coefficient_mu: linear absorption coefficient mu, 1/mm (absorption
          and scattering); calculated from the composition when it is not
          given (Diffrn.calc_absorption_coefficient)
        - coefficient_mu_half: mu at half the wavelength, 1/mm, used for the
          lambda/2 contamination (by default equal to coefficient_mu)
        - diffracted_beam: 'left' (default) or 'right', side of the detector
          of the normal beam geometry
        - correction_t_min, correction_t_max: smallest and largest
          transmission factors of the measured reflections (calculated)
        - process_details: how mu and the correction were obtained
          (calculated)
    """

    ATTR_MANDATORY_NAMES = ("correction_type", )
    ATTR_MANDATORY_TYPES = (str, )
    ATTR_MANDATORY_CIF = ("correction_type", )

    ATTR_OPTIONAL_NAMES = ("coefficient_mu", "coefficient_mu_half", "diffracted_beam",
                           "correction_t_min", "correction_t_max", "process_details")
    ATTR_OPTIONAL_TYPES = (float, float, str, float, float, str)
    ATTR_OPTIONAL_CIF = ("coefficient_mu", "coefficient_mu_half", "diffracted_beam",
                         "correction_T_min", "correction_T_max", "process_details")

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ()
    ATTR_INT_PROTECTED_NAMES = ()

    # parameters considered are refined parameters
    ATTR_REF = ("coefficient_mu", )
    ATTR_SIGMA = tuple([f"{_h:}_sigma" for _h in ATTR_REF])
    ATTR_CONSTR_FLAG = tuple([f"{_h:}_constraint" for _h in ATTR_REF])
    ATTR_REF_FLAG = tuple([f"{_h:}_refinement" for _h in ATTR_REF])
    ATTR_CONSTR_MARK = tuple([f"{_h:}_mark" for _h in ATTR_REF])

    # formats if cif format
    D_FORMATS = {"coefficient_mu": "{:.5f}", "coefficient_mu_half": "{:.5f}",
                 "correction_t_min": "{:.4f}", "correction_t_max": "{:.4f}"}

    # constraints on the parameters
    D_CONSTRAINTS = {"correction_type": ["none", "sphere", "gaussian"],
                     "diffracted_beam": ["left", "right"]}

    # default values for the parameters
    D_DEFAULT = {"correction_type": "none"}
    for key in ATTR_SIGMA:
        D_DEFAULT[key] = 0.
    for key in (ATTR_CONSTR_FLAG + ATTR_REF_FLAG):
        D_DEFAULT[key] = False
    for key in ATTR_CONSTR_MARK:
        D_DEFAULT[key] = ""

    PREFIX = "exptl_absorpt"

    def __init__(self, **kwargs) -> NoReturn:
        super(ExptlAbsorpt, self).__init__()

        # defined for any integer and float parameters
        D_MIN = {"coefficient_mu": 0., "coefficient_mu_half": 0.,
                 "correction_t_min": 0., "correction_t_max": 0.}

        # defined for ani integer and float parameters
        D_MAX = {"correction_t_min": 1., "correction_t_max": 1.}

        self.__dict__["D_MIN"] = D_MIN
        self.__dict__["D_MAX"] = D_MAX
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)

    def get_dictionary(self):
        res = {"absorption_correction_type": str(self.correction_type)}
        if self.is_attribute("coefficient_mu"):
            mu = float(self.coefficient_mu)
            res["absorption_mu"] = numpy.array([mu], dtype=float)
            res["flags_absorption_mu"] = numpy.array([self.coefficient_mu_refinement], dtype=bool)
            if self.is_attribute("coefficient_mu_half") and mu > 0.:
                res["absorption_mu_ratio_half"] = float(self.coefficient_mu_half) / mu
        if self.is_attribute("diffracted_beam"):
            res["absorption_diffracted_beam"] = str(self.diffracted_beam)
        return res

    def take_parameters_from_dictionary(self, ddict_diffrn, l_parameter_name: list = None,
                                        l_sigma: list = None):
        keys = ddict_diffrn.keys()
        if "absorption_mu" in keys and self.is_attribute("coefficient_mu"):
            self.coefficient_mu = abs(float(ddict_diffrn["absorption_mu"][0]))
        if "transmission" in keys:
            transmission = numpy.asarray(ddict_diffrn["transmission"], dtype=float)
            if transmission.size > 0:
                self.correction_t_min = float(transmission.min())
                self.correction_t_max = float(transmission.max())
        if l_parameter_name is not None:
            for name, sigma in zip(l_parameter_name, l_sigma):
                if name[0] == "absorption_mu":
                    self.coefficient_mu_sigma = float(sigma)


class ExptlAbsorptL(LoopN):
    """Absorption corrections in a loop."""

    ITEM_CLASS = ExptlAbsorpt
    ATTR_INDEX = None

    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(ExptlAbsorptL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name


class ExptlCrystal(ItemN):
    """
    Size and shape of the crystal.

    Attributes
    ----------
        - size_max, size_mid, size_min: largest, middle and smallest
          dimensions of the crystal, mm (the full lengths of the axes of the
          ellipsoid, not the semi-axes)
        - size_rad: radius of a spherical crystal, mm
        - ellipsoid_euler_alpha, ellipsoid_euler_beta, ellipsoid_euler_gamma:
          Euler angles (degrees, z-y-z) of the axes of the ellipsoid in the
          Cartesian system of the crystal (x || a*, z || c): the columns of
          R = Rz(alpha) Ry(beta) Rz(gamma) are the directions of size_max,
          size_mid and size_min
        - density_diffrn: density calculated from the crystal structure,
          g/cm^3 (calculated)
        - description: description of the crystal (shape, colour)
    """

    ATTR_MANDATORY_NAMES = ()
    ATTR_MANDATORY_TYPES = ()
    ATTR_MANDATORY_CIF = ()

    ATTR_OPTIONAL_NAMES = ("size_max", "size_mid", "size_min", "size_rad",
                           "ellipsoid_euler_alpha", "ellipsoid_euler_beta",
                           "ellipsoid_euler_gamma", "density_diffrn", "description")
    ATTR_OPTIONAL_TYPES = (float, float, float, float, float, float, float, float, str)
    ATTR_OPTIONAL_CIF = ("size_max", "size_mid", "size_min", "size_rad",
                         "ellipsoid_euler_alpha", "ellipsoid_euler_beta",
                         "ellipsoid_euler_gamma", "density_diffrn", "description")

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ()
    ATTR_INT_PROTECTED_NAMES = ()

    # parameters considered are refined parameters
    ATTR_REF = ("size_max", "size_mid", "size_min", "size_rad",
                "ellipsoid_euler_alpha", "ellipsoid_euler_beta", "ellipsoid_euler_gamma")
    ATTR_SIGMA = tuple([f"{_h:}_sigma" for _h in ATTR_REF])
    ATTR_CONSTR_FLAG = tuple([f"{_h:}_constraint" for _h in ATTR_REF])
    ATTR_REF_FLAG = tuple([f"{_h:}_refinement" for _h in ATTR_REF])
    ATTR_CONSTR_MARK = tuple([f"{_h:}_mark" for _h in ATTR_REF])

    # formats if cif format
    D_FORMATS = {"size_max": "{:.4f}", "size_mid": "{:.4f}", "size_min": "{:.4f}",
                 "size_rad": "{:.4f}", "ellipsoid_euler_alpha": "{:.3f}",
                 "ellipsoid_euler_beta": "{:.3f}", "ellipsoid_euler_gamma": "{:.3f}",
                 "density_diffrn": "{:.4f}"}

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

    PREFIX = "exptl_crystal"

    def __init__(self, **kwargs) -> NoReturn:
        super(ExptlCrystal, self).__init__()

        # defined for any integer and float parameters
        D_MIN = {"size_max": 0., "size_mid": 0., "size_min": 0., "size_rad": 0.,
                 "density_diffrn": 0.}

        # defined for ani integer and float parameters
        D_MAX = {}

        self.__dict__["D_MIN"] = D_MIN
        self.__dict__["D_MAX"] = D_MAX
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)

    def get_sizes(self):
        """Full lengths (mm) of the axes of the ellipsoid: size_max,
        size_mid, size_min; a missing one is taken from the others, a
        sphere has 2 size_rad. None when no size is given."""
        l_size = [float(getattr(self, name)) if self.is_attribute(name) else None
                  for name in ("size_max", "size_mid", "size_min")]
        l_known = [size for size in l_size if size is not None]
        if len(l_known) == 0:
            if self.is_attribute("size_rad"):
                return numpy.full(3, 2. * float(self.size_rad))
            return None
        fill = float(numpy.mean(l_known))
        return numpy.array([fill if size is None else size for size in l_size], dtype=float)

    def get_radius(self):
        """Radius (mm) of the sphere: size_rad, or the radius of the sphere
        of the same volume as the ellipsoid; None when no size is given."""
        if self.is_attribute("size_rad"):
            return float(self.size_rad)
        sizes = self.get_sizes()
        if sizes is None:
            return None
        return 0.5 * float(numpy.prod(sizes))**(1. / 3.)

    def get_euler_angles(self):
        """Euler angles alpha, beta, gamma (degrees) of the axes of the
        ellipsoid; zeros for the missing ones."""
        return numpy.array([float(getattr(self, f"ellipsoid_euler_{name:}"))
                            if self.is_attribute(f"ellipsoid_euler_{name:}") else 0.
                            for name in ("alpha", "beta", "gamma")], dtype=float)

    def set_euler_angles(self, euler_angles):
        """Set the Euler angles alpha, beta, gamma (degrees)."""
        for name, value in zip(("alpha", "beta", "gamma"), euler_angles):
            setattr(self, f"ellipsoid_euler_{name:}", round(float(value), 4) + 0.)

    def set_axes_by_directions(self, direction_max: list, direction_mid: list):
        """Orient the ellipsoid by the directions of its longest and middle
        axes.

        Parameters
        ----------
        direction_max : list
            Direction of the longest axis (size_max) as three numbers in the
            Cartesian system of the crystal (x || a*, z || c).
        direction_mid : list
            Direction of the middle axis (size_mid); its part along
            direction_max is removed.
        """
        v_max = numpy.asarray(direction_max, dtype=float)
        v_mid = numpy.asarray(direction_mid, dtype=float)
        if v_max.size != 3 or v_mid.size != 3:
            raise ValueError("The directions should be given by three numbers.")
        e_max = v_max / numpy.linalg.norm(v_max)
        v_mid = v_mid - numpy.dot(v_mid, e_max) * e_max
        norm_mid = numpy.linalg.norm(v_mid)
        if norm_mid < 1e-8:
            raise ValueError("The directions of the longest and middle axes should differ.")
        e_mid = v_mid / norm_mid
        matrix = numpy.stack([e_max, e_mid, numpy.cross(e_max, e_mid)], axis=1)
        self.set_euler_angles(calc_euler_zyz_by_rotation_matrix(matrix))

    def get_axes_matrix(self):
        """Directions of size_max, size_mid, size_min as the columns of a
        rotation matrix (3, 3), Cartesian system of the crystal."""
        return calc_rotation_matrix_by_euler_zyz(self.get_euler_angles())

    def get_dictionary(self):
        res = {}
        radius = self.get_radius()
        if radius is not None:
            res["absorption_radius"] = numpy.array([radius], dtype=float)
            res["flags_absorption_radius"] = numpy.array([self.size_rad_refinement], dtype=bool)
        sizes = self.get_sizes()
        if sizes is not None:
            res["absorption_size"] = sizes
            res["flags_absorption_size"] = numpy.array(
                [self.is_attribute(name) and getattr(self, f"{name:}_refinement")
                 for name in ("size_max", "size_mid", "size_min")], dtype=bool)
        res["absorption_euler"] = self.get_euler_angles()
        res["flags_absorption_euler"] = numpy.array(
            [self.is_attribute(f"ellipsoid_euler_{name:}")
             and getattr(self, f"ellipsoid_euler_{name:}_refinement")
             for name in ("alpha", "beta", "gamma")], dtype=bool)
        return res

    def take_parameters_from_dictionary(self, ddict_diffrn, l_parameter_name: list = None,
                                        l_sigma: list = None):
        keys = ddict_diffrn.keys()
        if "absorption_radius" in keys and self.is_attribute("size_rad"):
            self.size_rad = abs(float(ddict_diffrn["absorption_radius"][0]))
        if "absorption_size" in keys:
            for name, value in zip(("size_max", "size_mid", "size_min"), ddict_diffrn["absorption_size"]):
                if self.is_attribute(name):
                    setattr(self, name, abs(float(value)))
        if "absorption_euler" in keys:
            for name, value in zip(("alpha", "beta", "gamma"), ddict_diffrn["absorption_euler"]):
                if self.is_attribute(f"ellipsoid_euler_{name:}"):
                    setattr(self, f"ellipsoid_euler_{name:}", float(value))
        if l_parameter_name is not None:
            for name, sigma in zip(l_parameter_name, l_sigma):
                if name[0] == "absorption_radius":
                    self.size_rad_sigma = float(sigma)
                elif name[0] == "absorption_size":
                    i_comp = int(numpy.atleast_1d(name[1])[0])
                    setattr(self, f"{('size_max', 'size_mid', 'size_min')[i_comp]:}_sigma", float(sigma))
                elif name[0] == "absorption_euler":
                    i_comp = int(numpy.atleast_1d(name[1])[0])
                    setattr(self, f"ellipsoid_euler_{('alpha', 'beta', 'gamma')[i_comp]:}_sigma",
                            float(sigma))


class ExptlCrystalL(LoopN):
    """Crystals in a loop."""

    ITEM_CLASS = ExptlCrystal
    ATTR_INDEX = None

    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(ExptlCrystalL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name
