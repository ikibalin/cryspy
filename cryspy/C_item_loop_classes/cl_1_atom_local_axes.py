from typing import NoReturn

from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN


class AtomLocalAxes(ItemN):
    """
    This category allows the definition of local axes around each
    atom in terms of vectors between neighbouring atoms.
    High-resolution X-ray diffraction methods enable the
    determination of the electron density distribution in crystal
    lattices and molecules, which in turn allows for a
    characterization of chemical interactions (Coppens, 1997;
    Koritsanszky & Coppens, 2001). This is accomplished by the
    construction of a mathematical model of the charge density
    in a crystal and then by fitting the parameters of such a
    model to the experimental pattern of diffracted X-rays. The
    model on which this dictionary is based is the so-called
    multipole formalism proposed by Hansen & Coppens (1978). In
    this model, the electron density in a crystal is described
    by a sum of aspherical "pseudoatoms" where the pseudoatom
    density has the form defined in the _atom_rho_multipole_* items.
    Each pseudoatom density consists of terms representing the
    core density, the spherical part of the valence density and
    the deviation of the valence density from sphericity. The
    continuous electron density in the crystal is then modelled
    as a sum of atom-centred charge distributions. Once the
    experimental electron density has been established, the
    "atoms in molecules" theory of Bader (1990) provides tools for
    the interpretation of the density distribution in terms of its
    topological properties.

    Mandatory attributes:
        - wavelength

    Optional attributes:
        - field
        - offset_ttheta
        - offset_phi

    """
    ATTR_MANDATORY_NAMES = ("atom_label", "atom0", "ax1", "atom1", "atom2",
                            "ax2")
    ATTR_MANDATORY_TYPES = (str, str, str, str, str, str)
    ATTR_MANDATORY_CIF = ("atom_label", "atom0", "ax1", "atom1", "atom2",
                            "ax2")

    ATTR_OPTIONAL_NAMES = ()
    ATTR_OPTIONAL_TYPES = ()
    ATTR_OPTIONAL_CIF = ()

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ()
    ATTR_INT_PROTECTED_NAMES = ()

    # parameters considered are refined parameters
    ATTR_REF = ()
    ATTR_SIGMA = tuple([f"{_h:}_sigma" for _h in ATTR_REF])
    ATTR_CONSTR_FLAG = tuple([f"{_h:}_constraint" for _h in ATTR_REF])
    ATTR_REF_FLAG = tuple([f"{_h:}_refinement" for _h in ATTR_REF])
    ATTR_CONSTR_MARK = tuple([f"{_h:}_mark" for _h in ATTR_REF])

    # constraints on the parameters
    D_CONSTRAINTS = {"ax1": ["X", "Y", "Z", "-X", "-Y", "-Z"],
                     "ax2": ["X", "Y", "Z", "-X", "-Y", "-Z"]}

    # default values for the parameters
    D_DEFAULT = {}
    for key in ATTR_SIGMA:
        D_DEFAULT[key] = 0.
    for key in (ATTR_CONSTR_FLAG + ATTR_REF_FLAG):
        D_DEFAULT[key] = False
    for key in ATTR_CONSTR_MARK:
        D_DEFAULT[key] = ""

    PREFIX = "atom_local_axes"

    def __init__(self, **kwargs) -> NoReturn:
        super(AtomLocalAxes, self).__init__()

        # defined for any integer and float parameters
        D_MIN = {}

        # defined for ani integer and float parameters
        D_MAX = {}

        self.__dict__["D_MIN"] = D_MIN
        self.__dict__["D_MAX"] = D_MAX
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)


class AtomLocalAxesL(LoopN):
    """
    This category allows the definition of local axes around each
    atom in terms of vectors between neighbouring atoms.
    High-resolution X-ray diffraction methods enable the
    determination of the electron density distribution in crystal
    lattices and molecules, which in turn allows for a
    characterization of chemical interactions (Coppens, 1997;
    Koritsanszky & Coppens, 2001). This is accomplished by the
    construction of a mathematical model of the charge density
    in a crystal and then by fitting the parameters of such a
    model to the experimental pattern of diffracted X-rays. The
    model on which this dictionary is based is the so-called
    multipole formalism proposed by Hansen & Coppens (1978). In
    this model, the electron density in a crystal is described
    by a sum of aspherical "pseudoatoms" where the pseudoatom
    density has the form defined in the _atom_rho_multipole_* items.
    Each pseudoatom density consists of terms representing the
    core density, the spherical part of the valence density and
    the deviation of the valence density from sphericity. The
    continuous electron density in the crystal is then modelled
    as a sum of atom-centred charge distributions. Once the
    experimental electron density has been established, the
    "atoms in molecules" theory of Bader (1990) provides tools for
    the interpretation of the density distribution in terms of its
    topological properties.

    """
    ITEM_CLASS = AtomLocalAxes
    ATTR_INDEX = "atom_label"
    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(AtomLocalAxesL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name

# s_cont = """
#   loop_
#   _atom_local_axes_atom_label
#   _atom_local_axes_atom0
#   _atom_local_axes_ax1
#   _atom_local_axes_atom1
#   _atom_local_axes_atom2
#   _atom_local_axes_ax2
#       Ni2+(1)  DUM0      Z    Ni2+(1)  N(1)      X
#   """

# obj = AtomLocalAxesL.from_cif(s_cont)
# print(obj, end="\n\n")
# print(obj["Ni2+(1)"], end="\n\n")
