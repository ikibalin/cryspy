import numpy
from typing import NoReturn

from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN

class AtomType(ItemN):
    """Details about properties of the atom.

    Data items in the ATOM_TYPE category record details about
    properties of the atoms that occupy the atom sites, such as the
    atomic scattering factors.

    Mandatory attributes:
        - label
        - type_symbol
        - fract_x
        - fract_y
        - fract_z

    Optional attributes:
        - occupancy
        - adp_type
        - u_iso_or_equiv
        - u_equiv_geom_mean
        - b_iso_or_equiv
        - multiplicity
        - wyckoff_symbol
        - cartn_x
        - cartn_y
        - cartn_z

    Internal attributes:
        - scat_length_neutron

    Internal protected attributes:
        - space_group_wyckoff
        - constr_number
    """
    ATTR_MANDATORY_NAMES = ("symbol", )
    ATTR_MANDATORY_TYPES = (str, )
    ATTR_MANDATORY_CIF = ("symbol", )

    ATTR_OPTIONAL_NAMES = (
        "analytical_mass", "description", "number_in_cell", "oxidation_number",
        "radius_bond", "radius_contact", "cromer_mann_a1", "cromer_mann_a2",
        "cromer_mann_a3", "cromer_mann_a4", "cromer_mann_b1", "cromer_mann_b2",
        "cromer_mann_b3", "cromer_mann_b4", "cromer_mann_c",
        "scat_dispersion_real", "scat_dispersion_imag", "dispersion_source",
        "scat_length_neutron", "scat_source", "scat_versus_stol_list")
    ATTR_OPTIONAL_TYPES = (float, str, int, float, float, float, float, float,
                           float, float, float, float, float, float, float,
                           float, float, str, complex, str, float)
    ATTR_OPTIONAL_CIF = (
        "analytical_mass", "description", "number_in_cell", "oxidation_number",
        "radius_bond", "radius_contact", "cromer_mann_a1", "cromer_mann_a2",
        "cromer_mann_a3", "cromer_mann_a4", "cromer_mann_b1", "cromer_mann_b2",
        "cromer_mann_b3", "cromer_mann_b4", "cromer_mann_c",
        "scat_dispersion_real", "scat_dispersion_imag", "dispersion_source",
        "scat_length_neutron", "scat_source", "scat_versus_stol_list")

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
    D_CONSTRAINTS = {}

    # default values for the parameters
    D_DEFAULT = {}
    for key in ATTR_SIGMA:
        D_DEFAULT[key] = 0.
    for key in (ATTR_CONSTR_FLAG + ATTR_REF_FLAG):
        D_DEFAULT[key] = False
    for key in ATTR_CONSTR_MARK:
        D_DEFAULT[key] = ""

    PREFIX = "atom_type"

    def __init__(self, **kwargs) -> NoReturn:
        super(AtomType, self).__init__()

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



class AtomTypeL(LoopN):
    """Details about properties of the atoms.
    """
    ITEM_CLASS = AtomType
    ATTR_INDEX = "symbol"
    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(AtomTypeL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name



# s_cont = """
#  loop_
#  _atom_type_symbol
#  _atom_type_oxidation_number
#  _atom_type_number_in_cell
#  _atom_type_scat_dispersion_real
#  _atom_type_scat_dispersion_imag
#  _atom_type_scat_source
#    C 0 72  .017  .009  International_Tables_Vol_IV_Table_2.2B
#    H 0 100  0     0    International_Tables_Vol_IV_Table_2.2B
#    O 0 12  .047  .032  International_Tables_Vol_IV_Table_2.2B
#    N 0 4   .029  .018  International_Tables_Vol_IV_Table_2.2B
#   """

# obj = AtomTypeL.from_cif(s_cont)
# print(obj, end="\n\n")
# print(obj["O"].scat_length_neutron, end="\n\n")
