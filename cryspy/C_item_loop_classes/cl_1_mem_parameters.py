"""MEMParameters, MEMParametersL classes."""

from typing import NoReturn

import numpy

from cryspy.A_functions_base.function_1_objects import form_items_by_dictionary

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN


class MEMParameters(ItemN):
    """MEMParameters contain description of the parameters controlling MEM:

    - points_a, points_b, points_c are the number of points along axis a, b, c.
    - channel_col is True for spin density reconstruction when induced moments are along field
    - channel_ani is True for magnetisation density reconstruction when induced magnetic moment are defined
        as chi_a * H for paramagnetic compounds and as m_a for magnetically ordered compoud.
        The subindex 'a' is the index of magnetic ion.
    - channel_nucl is True for nuclear density reconstruction.
    - magnetization_plus, magnetization_minus define the total magnetization of unit cell at spin density
        reconstruction for spins along and opposit to the applied field
    - prior density can be 'uniform' for flat prior distribution in the unit cell
        or 'core' for spherical / multipole distribution (the corresponding parameters
        in the crystal should be defined)
    - use_asymmetry is True when asymmetry parameters are used for flipping ratio data during MEM reconstruction
        procedure. In the case of unpolarized single crystal data in this case the root of intensity is used
        for the reconstruction procedure.
    - gof_desired is the value of chi_sq at which the MEM reconstruction procedure should be stopped.
    - spin_density_to_den_file if the name of file with extension .den in which the reconstructed spin density should be saved
    - magnetization_density_to_den_file is the name of file with extension .den in which the reconstructed magnetisation density should be saved

        Attributes
        ----------
            - points_a, points_b, points_c, chi_ferro (mandatory)
            - chi_antiferro, prior_density, method (mandatory)
            - gof_desired (mandatory)

        Accessible Values
        -----------------
            method is "2channel" (default) or "tensorMEM"
            prior_density is "uniform" (default) or "core"
    """

    ATTR_MANDATORY_NAMES = (
        "points_a",
        "points_b",
        "points_c",
    )
    ATTR_MANDATORY_TYPES = (
        int,
        int,
        int,
    )
    ATTR_MANDATORY_CIF = (
        "points_a",
        "points_b",
        "points_c",
    )

    ATTR_OPTIONAL_NAMES = (
        "channel_col",
        "channel_ani",
        "channel_nucl",
        "magnetization_plus",
        "magnetization_minus",
        "prior_density",
        "use_asymmetry",
        "gof_desired",
        "spin_density_to_den_file",
        "magnetization_density_to_den_file",
        "nuclear_density_to_den_file",
    )
    ATTR_OPTIONAL_TYPES = (
        bool,
        bool,
        bool,
        float,
        float,
        str,
        bool,
        float,
        str,
        str,
        str,
    )
    ATTR_OPTIONAL_CIF = (
        "channel_col",
        "channel_ani",
        "channel_nucl",
        "magnetization_plus",
        "magnetization_minus",
        "prior_density",
        "use_asymmetry",
        "GoF_desired",
        "spin_density_to_den_file",
        "magnetization_density_to_den_file",
        "nuclear_density_to_den_file",
    )

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ()
    ATTR_INT_PROTECTED_NAMES = ()

    # parameters considered are refined parameters
    ATTR_REF = ("magnetization_plus", "magnetization_minus")
    ATTR_SIGMA = tuple([f"{_h:}_sigma" for _h in ATTR_REF])
    ATTR_CONSTR_FLAG = tuple([f"{_h:}_constraint" for _h in ATTR_REF])
    ATTR_REF_FLAG = tuple([f"{_h:}_refinement" for _h in ATTR_REF])
    ATTR_CONSTR_MARK = tuple([f"{_h:}_mark" for _h in ATTR_REF])

    # formats if cif format
    D_FORMATS = {
        "magnetization_plus": "{:.5f}",
        "magnetization_minus": "{:.5f}",
        "gof_desired": "{:.2f}",
    }

    # constraints on the parameters
    D_CONSTRAINTS = {
        "prior_density": ["core", "uniform"],
    }

    # default values for the parameters
    D_DEFAULT = {
        "points_a": 48,
        "points_b": 48,
        "points_c": 48,
        "channel_col": True,
        "magnetization_plus": 4.0,
        "magnetization_minus": -1.0,
        "channel_ani": False,
        "channel_nucl": False,
        "prior_density": "uniform",
        "use_asymmetry": True,
        "gof_desired": 1.0,
        # "spin_density_to_den_file": "mem_spin_density.den",
        # "magnetization_density_to_den_file": "mem_magnetization_density.den",
    }
    for key in ATTR_SIGMA:
        D_DEFAULT[key] = 0.0
    for key in ATTR_CONSTR_FLAG + ATTR_REF_FLAG:
        D_DEFAULT[key] = False
    for key in ATTR_CONSTR_MARK:
        D_DEFAULT[key] = ""

    PREFIX = "mem_parameters"

    def __init__(self, **kwargs) -> NoReturn:
        super(MEMParameters, self).__init__()

        # defined for any integer and float parameters
        D_MIN = {
            "points_a": 1,
            "points_b": 1,
            "points_c": 1,
            "magnetization_plus": 0.0,
            "gof_desired": 0.0,
        }

        # defined for ani integer and float parameters
        D_MAX = {"magnetization_minus": 0}

        self.__dict__["D_MIN"] = D_MIN
        self.__dict__["D_MAX"] = D_MAX
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)

    def get_dictionary(self):
        numb_a = self.points_a
        numb_b = self.points_b
        numb_c = self.points_c
        points_abc = numpy.array([numb_a, numb_b, numb_c], dtype=int)

        channel_col = self.channel_col
        channel_ani = self.channel_ani
        channel_nucl = self.channel_nucl
        magnetization_plus = self.magnetization_plus
        magnetization_minus = self.magnetization_minus
        flag_uniform_prior_density = self.prior_density == "uniform"
        flag_asymmetry = self.use_asymmetry
        gof_desired = self.gof_desired

        file_spin_density = getattr(self, "spin_density_to_den_file", None)
        file_magnetization_density = getattr(
            self, "magnetization_density_to_den_file", None
        )
        file_nuclear_density = getattr(
            self, "nuclear_density_to_den_file", None
        )
        dict_mem = {
            "points_abc": points_abc,
            "channel_col": channel_col,
            "channel_ani": channel_ani,
            "channel_nucl": channel_nucl,
            "magnetization_plus": magnetization_plus,
            "magnetization_minus": magnetization_minus,
            "flag_uniform_prior_density": flag_uniform_prior_density,
            "flag_asymmetry": flag_asymmetry,
            "gof_desired": gof_desired,
            "file_spin_density": file_spin_density,
            "file_magnetization_density": file_magnetization_density,
            "file_nuclear_density": file_nuclear_density,
        }
        return dict_mem


class MEMParametersL(LoopN):
    """Description of chi2 in loop."""

    ITEM_CLASS = MEMParameters
    ATTR_INDEX = None

    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(MEMParametersL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(
            self.ITEM_CLASS, kwargs
        )
        self.__dict__["loop_name"] = loop_name
