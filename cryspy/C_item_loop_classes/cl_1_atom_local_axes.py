from typing import NoReturn
import numpy
from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN


class AtomLocalAxes(ItemN):
    """AtomLocalAxes defines local axes around each atom using vectors between neighbouring atoms.

This is used to characterise chemical interactions (Coppens, 1997; Koritsanszky & Coppens, 2001).  
The process involves constructing a mathematical model of the charge or magnetisation density within a crystal and then fitting the model’s parameters. 
This dictionary is based on the multipole formalism proposed by Hansen and Coppens in 1978.

In this model, the electron density within a crystal is described as a sum of aspherical “pseudoatoms”.  
The pseudoatom density is defined in the atom_rho_multipole loop. 

For a more precise definition of local axes, dummy atoms (defined in the atom_site loop) can be used.
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
    D_CONSTRAINTS = {"ax1": ["X", "Y", "Z", "x", "y", "z"],
                     "ax2": ["X", "Y", "Z", "x", "y", "z"]}

    # default values for the parameters
    D_DEFAULT = {}
    for key in ATTR_SIGMA:
        D_DEFAULT[key] = 0.
    for key in (ATTR_CONSTR_FLAG + ATTR_REF_FLAG):
        D_DEFAULT[key] = False
    for key in ATTR_CONSTR_MARK:
        D_DEFAULT[key] = ""

    PREFIX = "atom_local_axes"

    D_REFERENCES = {"atom_label": "atom_site.label", "atom0": "atom_site.label",
                    "atom1": "atom_site.label", "atom2": "atom_site.label"}

    @classmethod
    def propose_values(cls, data, d_value: dict) -> dict:
        """atom1 is the atom itself; atom0 and atom2 are proposed among the
        nearest atoms (a dummy atom first, as it is put there to give an
        axis), ax1 is Z when atom0 is along c, else X; to be confirmed."""
        label = str(d_value.get("atom_label", ""))
        l_near = take_nearest_atoms(data, label)
        if len(l_near) == 0:
            return {}
        l_option = [near[1] for near in l_near]
        atom0, along_c, vector0 = l_near[0][1], l_near[0][2], l_near[0][4]
        # the second direction must not be along the first one
        l_other = [near for near in l_near[1:]
                   if abs(numpy.dot(near[4], vector0)) / (numpy.linalg.norm(near[4])
                                                          * numpy.linalg.norm(vector0)) < COS_SECOND_AXIS]
        atom2 = l_other[0][1] if l_other else ""
        ax1, ax2 = ("Z", "X") if along_c else ("X", "Z")
        reason = "the nearest atoms" + (" (a dummy atom first)" if l_near[0][3] else "")
        return {
            "atom1": {"proposal": label, "options": [label], "sure": True, "reason": "the atom itself"},
            "atom0": {"proposal": atom0, "options": l_option, "sure": False, "reason": reason},
            "ax1": {"proposal": ax1, "options": ["X", "Y", "Z"], "sure": False,
                    "reason": f"{atom0} is {'along c' if along_c else 'not along c'}"},
            "atom2": {"proposal": atom2, "options": l_option, "sure": False, "reason": reason},
            "ax2": {"proposal": ax2, "options": ["X", "Y", "Z"], "sure": False,
                    "reason": f"the axis normal to {ax1}"}}

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


# atoms nearer than this (Å) are proposed for the local axes
MAX_NEAR_DISTANCE = 3.
# cos of the angle below which a direction is taken as parallel to c (~5°)
COS_PARALLEL = 0.996
# the second axis is proposed only when it makes more than ~25° with the first
COS_SECOND_AXIS = 0.9


def take_nearest_atoms(data, label: str, max_number: int = 5) -> list:
    """(distance, label, along_c, is_dummy, Cartesian vector) of the atoms of atom_site
    nearest to the atom (lattice translations taken into account), dummy
    atoms first among those within 1 Å."""
    from cryspy.A_functions_base.unit_cell import calc_m_m_by_unit_cell_parameters
    if not (data.is_attribute("atom_site") and data.is_attribute("cell")) or \
            label not in list(data.atom_site.label):
        return []
    m_m = calc_m_m_by_unit_cell_parameters(data.cell.get_unit_cell_parameters())[0].reshape(3, 3)
    atom = data.atom_site[label]
    xyz = numpy.array([atom.fract_x, atom.fract_y, atom.fract_z], dtype=float)
    l_near = []
    for item in data.atom_site.items:
        diff = numpy.array([item.fract_x, item.fract_y, item.fract_z], dtype=float) - xyz
        diff -= numpy.round(diff)
        vector = m_m @ diff
        distance = float(numpy.linalg.norm(vector))
        if 1e-3 < distance <= MAX_NEAR_DISTANCE:
            is_dummy = str(item.type_symbol).strip().lower() == "dummy"
            l_near.append((distance, str(item.label), bool(abs(vector[2]) / distance > COS_PARALLEL),
                           is_dummy, vector))
    l_near.sort(key=lambda near: (not (near[3] and near[0] <= 1.), near[0], near[1]))
    return l_near[:max_number]


class AtomLocalAxesL(LoopN):
    """AtomLocalAxes defines local axes around each atom using vectors between neighbouring atoms.

This is used to characterise chemical interactions (Coppens, 1997; Koritsanszky & Coppens, 2001).  
The process involves constructing a mathematical model of the charge or magnetisation density within a crystal and then fitting the model’s parameters. 
This dictionary is based on the multipole formalism proposed by Hansen and Coppens in 1978.

In this model, the electron density within a crystal is described as a sum of aspherical “pseudoatoms”.  
The pseudoatom density is defined in the atom_rho_multipole loop. 

For a more precise definition of local axes, dummy atoms (defined in the atom_site loop) can be used.
    """
    ITEM_CLASS = AtomLocalAxes
    ATTR_INDEX = "atom_label"
    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(AtomLocalAxesL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name

    def get_dictionary(self):
        d_out = super(AtomLocalAxesL, self).get_dictionary()
        d_out["atom_local_axes_label"] = numpy.array( self.atom_label, dtype=str)
        return d_out

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
