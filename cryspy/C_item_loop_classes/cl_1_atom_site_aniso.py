"""Desctiption of AtomSiteAniso and AtomSiteAnisoL."""
import numpy
from typing import NoReturn

from cryspy.A_functions_base.function_1_atomic_vibrations import \
    calc_beta_by_u, vibration_constraints
from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN


class AtomSiteAniso(ItemN):
    """
    AtomSiteAniso class.

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

    ATTR_MANDATORY_NAMES = ("label", )
    ATTR_MANDATORY_TYPES = (str, )
    ATTR_MANDATORY_CIF = ("label", )

    ATTR_OPTIONAL_NAMES = ("u_11", "u_22", "u_33", "u_12", "u_13", "u_23",
                           "b_11", "b_22", "b_33", "b_12", "b_13", "b_23")
    ATTR_OPTIONAL_TYPES = (float, float, float, float, float, float,
                           float, float, float, float, float, float)
    ATTR_OPTIONAL_CIF = ("U_11", "U_22", "U_33", "U_12", "U_13", "U_23",
                         "B_11", "B_22", "B_33", "B_12", "B_13", "B_23")

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ("scat_length_neutron", )
    ATTR_INT_PROTECTED_NAMES = ("space_group_wyckoff", "constr_number")

    # parameters considered are refined parameters
    ATTR_REF = ("u_11", "u_22", "u_33", "u_12", "u_13", "u_23",
                "b_11", "b_22", "b_33", "b_12", "b_13", "b_23")
    ATTR_SIGMA = tuple([f"{_h:}_sigma" for _h in ATTR_REF])
    ATTR_CONSTR_FLAG = tuple([f"{_h:}_constraint" for _h in ATTR_REF])
    ATTR_REF_FLAG = tuple([f"{_h:}_refinement" for _h in ATTR_REF])
    ATTR_CONSTR_MARK = tuple([f"{_h:}_mark" for _h in ATTR_REF])

    # formats if cif format
    D_FORMATS = {"u_11": "{:.5f}", "u_22": "{:.5f}", "u_33": "{:.5f}",
                 "u_12": "{:.5f}", "u_13": "{:.5f}", "u_23": "{:.5f}"}

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

    PREFIX = "atom_site_aniso"

    def __init__(self, **kwargs) -> NoReturn:
        super(AtomSiteAniso, self).__init__()

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

    def calc_beta(self, cell):
        """
        Calculate $beta_{ij}$ from $U_{ij}$.

        Parameters
        ----------
        cell : TYPE
            DESCRIPTION.

        Returns
        -------
        beta_11 : TYPE
            DESCRIPTION.
        beta_22 : TYPE
            DESCRIPTION.
        beta_33 : TYPE
            DESCRIPTION.
        beta_12 : TYPE
            DESCRIPTION.
        beta_13 : TYPE
            DESCRIPTION.
        beta_23 : TYPE
            DESCRIPTION.

        """
        u_i = (self.u_11, self.u_22, self.u_33,
               self.u_12, self.u_13, self.u_23)
        beta_11, beta_22, beta_33, beta_12, beta_13, beta_23 = \
            calc_beta_by_u(u_i, cell)
        return beta_11, beta_22, beta_33, beta_12, beta_13, beta_23

    def apply_space_group_constraint(self, atom_site, space_group):
        """
        Constraints of the space group.

        According to table 1 in Peterse, Palm, Acta Cryst.(1966), 20, 147
        """

        l_numb = atom_site.calc_constr_number(space_group)
        label_aniso = self.label
        label = atom_site.label
        index = label.index(label_aniso)

        self.__dict__["u_11_constraint"] = False
        self.__dict__["u_22_constraint"] = False
        self.__dict__["u_33_constraint"] = False
        self.__dict__["u_12_constraint"] = False
        self.__dict__["u_13_constraint"] = False
        self.__dict__["u_23_constraint"] = False

        numb = l_numb[index]
        
        u_i = (self.u_11, self.u_22, self.u_33, self.u_12, self.u_13,
               self.u_23)
        u_sigma_i = (self.u_11_sigma, self.u_22_sigma, self.u_33_sigma,
                     self.u_12_sigma, self.u_13_sigma, self.u_23_sigma)
        u_ref_i = (self.u_11_refinement, self.u_22_refinement,
                   self.u_33_refinement, self.u_12_refinement,
                   self.u_13_refinement, self.u_23_refinement)
        
        u_i, u_sigma_i, u_ref_i, u_constr_i = vibration_constraints(
            numb, u_i, u_sigma_i, u_ref_i)

        self.__dict__["u_11"], self.__dict__["u_22"], self.__dict__["u_33"], \
            self.__dict__["u_12"], self.__dict__["u_13"], \
            self.__dict__["u_23"] = u_i

        self.__dict__["u_11_sigma"], self.__dict__["u_22_sigma"], \
            self.__dict__["u_33_sigma"], self.__dict__["u_12_sigma"], \
            self.__dict__["u_13_sigma"], self.__dict__["u_23_sigma"] = \
            u_sigma_i

        self.__dict__["u_11_refinement"], self.__dict__["u_22_refinement"], \
            self.__dict__["u_33_refinement"], \
            self.__dict__["u_12_refinement"], \
            self.__dict__["u_13_refinement"], \
            self.__dict__["u_23_refinement"] = u_ref_i

        self.__dict__["u_11_constraint"], self.__dict__["u_22_constraint"], \
            self.__dict__["u_33_constraint"], \
            self.__dict__["u_12_constraint"], \
            self.__dict__["u_13_constraint"], \
            self.__dict__["u_23_constraint"] = u_constr_i


class AtomSiteAnisoL(LoopN):
    """Description of AtomSite in loop."""

    ITEM_CLASS = AtomSiteAniso
    ATTR_INDEX = "label"

    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(AtomSiteAnisoL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name

    def calc_beta(self, cell):
        """
        Calculate $beta_{ij}$ from $U_{ij}$.

        Output:
            numpy array of beta (11, 22, 33, 12, 13, 23)
        """
        l_beta = []
        for item in self.items:
            beta_11, beta_22, beta_33, beta_12, beta_13, beta_23 = \
                item.calc_beta(cell)
            l_beta.append((beta_11, beta_22, beta_33, beta_12, beta_13,
                           beta_23))
        return numpy.array(l_beta, dtype=float)

    def apply_space_group_constraint(self, atom_site, space_group):
        """Space group constraint."""
        for item in self.items:
            item.apply_space_group_constraint(atom_site, space_group)

# s_cont = """
#  loop_
#  _atom_site_aniso_label
#  _atom_site_aniso_U_11
#  _atom_site_aniso_U_22
#  _atom_site_aniso_U_33
#  _atom_site_aniso_U_12
#  _atom_site_aniso_U_13
#  _atom_site_aniso_U_23
#   O1   .071(1) .076(1) .0342(9) .008(1)   .0051(9) -.0030(9)
#   C2   .060(2) .072(2) .047(1)  .002(2)   .013(1)  -.009(1)
#   """

# obj = AtomSiteAnisoL.from_cif(s_cont)
# print(obj, end="\n\n")
# print(obj["O1"], end="\n\n")
