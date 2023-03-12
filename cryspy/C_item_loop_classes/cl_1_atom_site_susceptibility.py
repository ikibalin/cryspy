"""Describe AtomSiteSusceptibility, AtomSiteSusceptibilityL."""
from typing import NoReturn
import math
import numpy

from cryspy.A_functions_base.function_1_algebra import calc_m_sigma
from cryspy.A_functions_base.function_1_atomic_vibrations import \
    vibration_constraints
from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary

import cryspy.A_functions_base.local_susceptibility as local_susceptibility

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN

na = numpy.newaxis

# FIXME: temporary solution which is enough slow
def calc_independent_lines(m):
    n_x, n_y = m.shape
    l_ind = [0, ]
    for ind_x in range(1, n_x):
        det_m = numpy.linalg.det(m[l_ind+[ind_x, ], :][:, l_ind+[ind_x, ]])
        if numpy.isclose(det_m, 0.):
            pass
        else:
            l_ind.append(ind_x)
    return l_ind


def calc_constr_matrix(matrix):
    m_x, m_y = matrix.shape
    m_rank = numpy.linalg.matrix_rank(matrix)
    unity = numpy.diag(m_x*[1.])
    m_u = matrix - unity


    for i_x in range(m_x-1, 0, -1):
        coeff_1 = m_u[i_x, i_x]
        if not(numpy.isclose(coeff_1, 0., atol=1e-10)):
            line_1 = m_u[i_x, :]/coeff_1
            m_u[i_x, :] = line_1
            for i_y in range(i_x-1, -1, -1):
                coeff_2 = m_u[i_y, i_x]
                line_2 = m_u[i_y, :] - line_1*coeff_2
                m_u[i_y, :]  = line_2

    res = numpy.round(unity-m_u, 10)
    return res

class AtomSiteSusceptibility(ItemN):
    """Magnetic properties of the atom that occupy the atom site.

    Data items in the ATOM_SITE_MAGNETISM_ANISO category record details about
    magnetic properties of the atoms that occupy the atom sites.

    Attributes
    ----------
        -
    """

    ATTR_MANDATORY_NAMES = ("label", )
    ATTR_MANDATORY_TYPES = (str, )
    ATTR_MANDATORY_CIF = ("label", )

    ATTR_OPTIONAL_NAMES = (
        "chi_type", "chi_11", "chi_22", "chi_33", "chi_12",
        "chi_13", "chi_23", )
    ATTR_OPTIONAL_TYPES = (str, float, float, float, float, float, float, )
    ATTR_OPTIONAL_CIF = (
        "chi_type",  "chi_11", "chi_22", "chi_33", "chi_12",
        "chi_13", "chi_23", )

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ()
    ATTR_INT_PROTECTED_NAMES = ()

    # parameters considered are refined parameters
    ATTR_REF = ("chi_11", "chi_22", "chi_33", "chi_12", "chi_13", "chi_23",)
    ATTR_SIGMA = tuple([f"{_h:}_sigma" for _h in ATTR_REF])
    ATTR_CONSTR_FLAG = tuple([f"{_h:}_constraint" for _h in ATTR_REF])
    ATTR_REF_FLAG = tuple([f"{_h:}_refinement" for _h in ATTR_REF])
    ATTR_CONSTR_MARK = tuple([f"{_h:}_mark" for _h in ATTR_REF])

    # formats if cif format
    D_FORMATS = {"chi_11": "{:.5f}", "chi_22": "{:.5f}", "chi_33": "{:.5f}",
                 "chi_12": "{:.5f}", "chi_13": "{:.5f}", "chi_23": "{:.5f}",}

    # constraints on the parameters
    D_CONSTRAINTS = {"chi_type": ["Ciso", "Cani"],
                     }

    # default values for the parameters
    D_DEFAULT = {}
    for key in ATTR_SIGMA:
        D_DEFAULT[key] = 0.
    for key in (ATTR_CONSTR_FLAG + ATTR_REF_FLAG):
        D_DEFAULT[key] = False
    for key in ATTR_CONSTR_MARK:
        D_DEFAULT[key] = ""

    PREFIX = "atom_site_susceptibility"

    def __init__(self, **kwargs) -> NoReturn:
        super(AtomSiteSusceptibility, self).__init__()

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

    def apply_space_group_constraint(self, atom_site, space_group, cell):
        """
        Space group constraints.

        According to table 1 in Peterse, Palm, Acta Cryst.(1966), 20, 147
        """
        label_aniso = self.label
        label = atom_site.label
        index = label.index(label_aniso)

        item_as = atom_site.items[index]
        sc_chi = item_as.calc_sc_chi(space_group, cell)
        # it should be checked
        sc_chi = calc_constr_matrix(sc_chi)

        # l_numb = atom_site.calc_constr_number(space_group)

        flag_chi = self.is_attribute("chi_type")
        if flag_chi:
            flag_chi = self.chi_type.lower().startswith("cani")

        if flag_chi:
            self.__dict__["chi_11_constraint"] = False
            self.__dict__["chi_22_constraint"] = False
            self.__dict__["chi_33_constraint"] = False
            self.__dict__["chi_12_constraint"] = False
            self.__dict__["chi_13_constraint"] = False
            self.__dict__["chi_23_constraint"] = False

        if flag_chi:
            chi_i = numpy.array([
                self.chi_11, self.chi_22, self.chi_33,
                self.chi_12, self.chi_13, self.chi_23], dtype=float)

            chi_sigma_i = numpy.array([
                self.chi_11_sigma, self.chi_22_sigma, self.chi_33_sigma,
                self.chi_12_sigma, self.chi_13_sigma, self.chi_23_sigma], dtype=float)

            chi_ref_i = numpy.array([
                self.chi_11_refinement, self.chi_22_refinement, self.chi_33_refinement,
                self.chi_12_refinement, self.chi_13_refinement, self.chi_23_refinement], dtype=bool)

            chi_i_c = sc_chi.dot(chi_i)
            chi_sigma_i_c = numpy.sqrt(numpy.square(sc_chi).dot(numpy.square(chi_sigma_i)))
            l_ind_independent = calc_independent_lines(sc_chi)

            chi_con_i_c = numpy.ones(chi_i_c.shape, dtype=bool)
            chi_con_i_c[l_ind_independent] = False
            sc_chi_bool = numpy.logical_not(numpy.isclose(numpy.round(sc_chi, decimals=5), 0.))
            # chi_con_i_c = numpy.triu((sc_chi_bool[na, :, :] == sc_chi_bool[:, na, :]).all(axis=2), k=1).any(axis=0)
            chi_ref_i_c = sc_chi_bool.dot(chi_ref_i) * numpy.logical_not(chi_con_i_c)
            # # I am not quite sure. It is to fix two parameters among three.
            # chi_ref_i_c = numpy.logical_and(chi_ref_i_c, chi_ref_i)

            self.__dict__["chi_11"], self.__dict__["chi_22"], \
                self.__dict__["chi_33"], self.__dict__["chi_12"], \
                self.__dict__["chi_13"], self.__dict__["chi_23"] = \
                    chi_i_c[0], chi_i_c[1], chi_i_c[2], chi_i_c[3], chi_i_c[4], chi_i_c[5]

            self.__dict__["chi_11_sigma"], self.__dict__["chi_22_sigma"], \
                self.__dict__["chi_33_sigma"], self.__dict__["chi_12_sigma"], \
                self.__dict__["chi_13_sigma"], self.__dict__["chi_23_sigma"] =\
                chi_sigma_i_c[0], chi_sigma_i_c[1], chi_sigma_i_c[2], chi_sigma_i_c[3], \
                chi_sigma_i_c[4], chi_sigma_i_c[5]

            self.__dict__["chi_11_refinement"], \
                self.__dict__["chi_22_refinement"], \
                self.__dict__["chi_33_refinement"], \
                self.__dict__["chi_12_refinement"], \
                self.__dict__["chi_13_refinement"], \
                self.__dict__["chi_23_refinement"] = chi_ref_i_c[0], chi_ref_i_c[1], chi_ref_i_c[2],\
                    chi_ref_i_c[3], chi_ref_i_c[4], chi_ref_i_c[5]

            self.__dict__["chi_11_constraint"], \
                self.__dict__["chi_22_constraint"], \
                self.__dict__["chi_33_constraint"], \
                self.__dict__["chi_12_constraint"], \
                self.__dict__["chi_13_constraint"], \
                self.__dict__["chi_23_constraint"] = chi_con_i_c[0], chi_con_i_c[1], chi_con_i_c[2],\
                    chi_con_i_c[3], chi_con_i_c[4], chi_con_i_c[5]

    def apply_chi_iso_constraint(self, cell):
        """Isotropic constraint on susceptibility."""
        c_a = cell.cos_a
        s_ib = cell.sin_ib
        s_ig = cell.sin_ig
        c_ib = cell.cos_ib
        c_ig = cell.cos_ig
        # not sure, it is better to check
        if not(self.is_attribute("chi_type")):
            return
        chi_type = self.chi_type
        if chi_type.lower().startswith("ciso"):
            self.__dict__["chi_22"] = self.chi_11
            self.__dict__["chi_33"] = self.chi_11
            self.__dict__["chi_12"] = 0. # self.chi_11*c_ig
            self.__dict__["chi_13"] = 0. # self.chi_11*c_ib
            self.__dict__["chi_23"] = 0. # self.chi_11*(c_ib*c_ig-s_ib*s_ig*c_a)
            self.__dict__["chi_22_sigma"] = self.chi_11_sigma
            self.__dict__["chi_33_sigma"] = self.chi_11_sigma
            self.__dict__["chi_12_sigma"] = 0. # self.chi_11_sigma * c_ig
            self.__dict__["chi_13_sigma"] = 0. # self.chi_11_sigma * c_ib
            self.__dict__["chi_23_sigma"] = 0. # self.chi_11_sigma * (c_ib*c_ig-s_ib*s_ig*c_a)
            self.__dict__["chi_22_refinement"] = False
            self.__dict__["chi_33_refinement"] = False
            self.__dict__["chi_12_refinement"] = False
            self.__dict__["chi_13_refinement"] = False
            self.__dict__["chi_23_refinement"] = False
            self.__dict__["chi_22_constraint"] = True
            self.__dict__["chi_33_constraint"] = True
            self.__dict__["chi_12_constraint"] = True
            self.__dict__["chi_13_constraint"] = True
            self.__dict__["chi_23_constraint"] = True

    def calc_main_axes_of_magnetization_ellipsoid(self, cell):
        """Susceptibility along the main axes of magnetization ellipsoid.

        Arguments
        ---------
            - cell

        Output
        ------
            - moments is main axes of ellipsoid in mu_B/T
            - moments_sigma is sigmas for main axes of ellipsoid
            - rot_matrix is directions for moments
                for moments[0] direction is rot_matrix[:, 0]
                for moments[1] direction is rot_matrix[:, 1]
                for moments[2] direction is rot_matrix[:, 2]

        The main axes are given in Cartezian coordinate system (x||a*, z||c).
        """

        ucp = numpy.array([
            cell.length_a, cell.length_b, cell.length_c,
            cell.angle_alpha*numpy.pi/180., cell.angle_beta*numpy.pi/180.,
            cell.angle_gamma*numpy.pi/180.], dtype=float)
        flag_ucp = False

        q_rn = numpy.array([
            self.chi_11, self.chi_22, self.chi_33,
            self.chi_12, self.chi_13, self.chi_23], dtype=float)
        flag_q_rn = True

        # FIXME: add errorbars calculations
        q_rn_sigma = numpy.array([
            self.chi_11_sigma, self.chi_22_sigma, self.chi_33_sigma,
            self.chi_12_sigma, self.chi_13_sigma, self.chi_23_sigma], dtype=float)

        moments, eig_fields, eig_moments = local_susceptibility.calc_magnetization_ellipsoid_axes(
            q_rn, ucp)

        moments_sigma = numpy.zeros_like(moments)
        return moments, moments_sigma, eig_fields, eig_moments

    def get_flags_susceptibility(self):
        res = numpy.array([
            self.chi_11_refinement, self.chi_22_refinement, self.chi_33_refinement,
            self.chi_12_refinement, self.chi_13_refinement, self.chi_23_refinement], dtype=bool)
        return res



class AtomSiteSusceptibilityL(LoopN):
    """Magnetic properties of the atoms that occupy the atom sites.

    Methods
    -------
        - apply_space_group_constraint
        - apply_chi_iso_constraint
    """

    ITEM_CLASS = AtomSiteSusceptibility
    ATTR_INDEX = "label"

    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(AtomSiteSusceptibilityL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name

    def apply_space_group_constraint(self, atom_site, space_group, cell):
        """Apply space group constraint."""
        for item in self.items:
            item.apply_space_group_constraint(atom_site, space_group, cell)

    def apply_chi_iso_constraint(self, cell):
        """Apply isotropic constraint on susceptibility."""
        for item in self.items:
            item.apply_chi_iso_constraint(cell)


    def calc_main_axes_of_magnetization_ellipsoid(self, cell):
        """Susceptibility along the main axes of magnetization ellipsoid.

        Arguments
        ---------
            - cell

        Output
        ------
            - l_moments is main axes of ellipsoid in mu_B/T for each atom
            - l_moments_sigma is sigmas for main axes of ellipsoid for each
              atom
            - l_rot_matrix is directions for moments
                for moments[0] direction is rot_matrix[:, 0]
                for moments[1] direction is rot_matrix[:, 1]
                for moments[2] direction is rot_matrix[:, 2]

        The main axes are given in Cartezian coordinate system (x||a*, z||c).
        """
        l_moments, l_moments_sigma, l_eig_fields, l_eig_moments = [], [], [], []
        for item in self.items:
            moments, moments_sigma, eig_fields, eig_moments = \
                item.calc_main_axes_of_magnetization_ellipsoid(cell)
            l_moments.append(moments)
            l_moments_sigma.append(moments_sigma)
            l_eig_fields.append(eig_fields)
            l_eig_moments.append(eig_moments)
        return l_moments, l_moments_sigma, l_eig_fields, l_eig_moments


    def get_flags_susceptibility(self):
        flags_susceptibility = numpy.stack([item.get_flags_susceptibility() for item in self.items], axis=0).transpose()
        return flags_susceptibility

