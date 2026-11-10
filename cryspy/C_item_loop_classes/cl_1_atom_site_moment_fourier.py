"""AtomSiteMomentFourier, AtomSiteMomentFourierL classes are given."""
from typing import NoReturn

from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN


class AtomSiteMomentFourier(ItemN):
    """Fourier component of the magnetic moment modulation (mCIF/msCIF).

The moment of the atom at the average position r (fractional, the lattice
translation included) is

    M(x4) = cos * cos(2 pi x4) + sin * sin(2 pi x4),   x4 = k.r,

per axis (x, y, z: along the unit vectors of a, b, c in Bohr magnetons, as
``_atom_site_moment.crystalaxis_x``). One row is given per atom, axis and
wave vector (``wave_vector_seq_id`` refers to atom_site_Fourier_wave_vector,
or to cell_wave_vector without it). ``cos_symmform`` and ``sin_symmform``
(e.g. ``mxs1``, ``-mxs1``, ``0``) are the symmetry constraints of the
superspace group: equal symbols are equal parameters, ``0`` is zero. These
are the tags of the Bilbao MAGNDATA database:
``_atom_site_moment_Fourier.atom_site_label``, ``.axis``,
``.wave_vector_seq_id`` and ``_atom_site_moment_Fourier_param.cos``,
``.sin``, ``.cos_symmform``, ``.sin_symmform``.

The label must match the label in the atom_site loop, and the magnetic form
factor is taken from the atom_site_scat loop.
    """

    ATTR_MANDATORY_NAMES = ("atom_site_label", "axis")
    ATTR_MANDATORY_TYPES = (str, str)
    ATTR_MANDATORY_CIF = ("atom_site_label", "axis")

    ATTR_OPTIONAL_NAMES = (
        "id", "wave_vector_seq_id", "cos", "sin", "cos_symmform", "sin_symmform")
    ATTR_OPTIONAL_TYPES = (str, int, float, float, str, str)
    ATTR_OPTIONAL_CIF = (
        "id", "wave_vector_seq_id", "param.cos", "param.sin",
        "param.cos_symmform", "param.sin_symmform")

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ()
    ATTR_INT_PROTECTED_NAMES = ()

    # parameters considered are refined parameters
    ATTR_REF = ("cos", "sin")
    ATTR_SIGMA = tuple([f"{_h:}_sigma" for _h in ATTR_REF])
    ATTR_CONSTR_FLAG = tuple([f"{_h:}_constraint" for _h in ATTR_REF])
    ATTR_REF_FLAG = tuple([f"{_h:}_refinement" for _h in ATTR_REF])
    ATTR_CONSTR_MARK = tuple([f"{_h:}_mark" for _h in ATTR_REF])

    D_FORMATS = {"cos": "{:.5f}", "sin": "{:.5f}"}

    D_CONSTRAINTS = {"axis": ["x", "y", "z", "a", "b", "c"]}

    D_DEFAULT = {"wave_vector_seq_id": 1, "cos": 0., "sin": 0.}
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
    """Fourier components of the magnetic moment modulations (mCIF/msCIF).

See AtomSiteMomentFourier for the convention.
    """

    ITEM_CLASS = AtomSiteMomentFourier
    ATTR_INDEX = None

    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(AtomSiteMomentFourierL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name

    def take_rows(self, label: str) -> dict:
        """{axis: item} of the atom with the label."""
        return {item.axis: item for item in self.items if item.atom_site_label == label}
