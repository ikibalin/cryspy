"""AtomSiteMomentFourier, AtomSiteMomentFourierL classes are given."""
from typing import NoReturn

from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN


class AtomSiteMomentFourier(ItemN):
    """Fourier component of the magnetic moment for a propagation vector.

The moment of the atom in the lattice cell R_L is
m = sum over k of S_k exp(-2 pi i k.R_L) (FullProf convention), where
S_k = (real + i imag) exp(-2 pi i phase). The components are given along
the unit vectors of the crystal axes a, b, c in Bohr magnetons, like
FullProf Rx, Ry, Rz, Ix, Iy, Iz and MagPh. For k equivalent to -k, for
example k = (1/2, 0, 0), the real part is the moment itself.

k_number refers to the magnetic_structure_k loop. The label must match the
label in the atom_site loop, and the magnetic form factor is taken from the
atom_site_scat loop. The symmetry of the crystal (nuclear and magnetic
together) is given by the space_group_symop_magn_k and
space_group_symop_magn_k_centering loops. Atom positions are reduced to the
first cell [0, 1) and the moment phase refers to that cell.
    """

    ATTR_MANDATORY_NAMES = ("label", )
    ATTR_MANDATORY_TYPES = (str, )
    ATTR_MANDATORY_CIF = ("label", )

    ATTR_OPTIONAL_NAMES = (
        "k_number", "real_x", "real_y", "real_z", "imag_x", "imag_y",
        "imag_z", "phase")
    ATTR_OPTIONAL_TYPES = (int, float, float, float, float, float, float,
                           float)
    ATTR_OPTIONAL_CIF = ATTR_OPTIONAL_NAMES

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ()
    ATTR_INT_PROTECTED_NAMES = ()

    # parameters considered are refined parameters
    ATTR_REF = ("real_x", "real_y", "real_z", "imag_x", "imag_y", "imag_z",
                "phase")
    ATTR_SIGMA = tuple([f"{_h:}_sigma" for _h in ATTR_REF])
    ATTR_CONSTR_FLAG = tuple([f"{_h:}_constraint" for _h in ATTR_REF])
    ATTR_REF_FLAG = tuple([f"{_h:}_refinement" for _h in ATTR_REF])
    ATTR_CONSTR_MARK = tuple([f"{_h:}_mark" for _h in ATTR_REF])

    D_FORMATS = {_h: "{:.5f}" for _h in ATTR_REF}

    D_CONSTRAINTS = {}

    D_DEFAULT = {"k_number": 1, "real_x": 0., "real_y": 0., "real_z": 0.,
                 "imag_x": 0., "imag_y": 0., "imag_z": 0., "phase": 0.}
    for key in ATTR_SIGMA:
        D_DEFAULT[key] = 0.
    for key in (ATTR_CONSTR_FLAG + ATTR_REF_FLAG):
        D_DEFAULT[key] = False
    for key in ATTR_CONSTR_MARK:
        D_DEFAULT[key] = ""

    PREFIX = "atom_site_moment_fourier"

    def __init__(self, **kwargs) -> NoReturn:
        super(AtomSiteMomentFourier, self).__init__()
        self.__dict__["D_MIN"] = {}
        self.__dict__["D_MAX"] = {}
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)


class AtomSiteMomentFourierL(LoopN):
    """Fourier components of the magnetic moments for propagation vectors.

See AtomSiteMomentFourier for the convention.
    """

    ITEM_CLASS = AtomSiteMomentFourier
    ATTR_INDEX = "label"

    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(AtomSiteMomentFourierL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name
