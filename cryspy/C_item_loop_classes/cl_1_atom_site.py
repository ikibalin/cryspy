import numpy
from typing import NoReturn

from cryspy.A_functions_base.symmetry_constraints import calc_sc_chi
from cryspy.A_functions_base.matrix_operations import calc_m1_m2_inv_m1
from cryspy.A_functions_base.unit_cell import calc_m_m_by_unit_cell_parameters

from cryspy.A_functions_base.function_1_scat_length_neutron import \
    get_scat_length_neutron
from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN

class AtomSite(ItemN):
    """Describe the atom site.

    Attributes
    ----------
        - label, type_symbol, fract_x, fract_y, fract_z (mandatory)
        - occupancy, adp_type, u_iso_or_equiv, u_equiv_geom_mean,
          b_iso_or_equiv, multiplicity, wyckoff_symbol, cartn_x, cartn_y,
          cartn_z (optional)
        - scat_length_neutron (internal)
        - space_group_wyckoff, constr_number (internal protected)
    """
    ATTR_MANDATORY_NAMES = ("label", "type_symbol", "fract_x", "fract_y",
                            "fract_z")
    ATTR_MANDATORY_TYPES = (str, str, float, float, float)
    ATTR_MANDATORY_CIF = ("label", "type_symbol", "fract_x", "fract_y",
                          "fract_z")

    ATTR_OPTIONAL_NAMES = (
        "occupancy", "adp_type",  "u_iso_or_equiv", "u_equiv_geom_mean",
        "b_iso_or_equiv", "multiplicity", "wyckoff_symbol", "cartn_x",
        "cartn_y", "cartn_z")
    ATTR_OPTIONAL_TYPES = (float, str, float, float, float, int, str, float,
                           float, float)
    ATTR_OPTIONAL_CIF = (
        "occupancy", "adp_type",  "U_iso_or_equiv", "U_equiv_geom_mean",
        "B_iso_or_equiv", "multiplicity", "Wyckoff_symbol", "Cartn_x",
        "Cartn_y", "Cartn_z")

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ("scat_length_neutron", )
    ATTR_INT_PROTECTED_NAMES = ("space_group_wyckoff", "constr_number", "sc_chi")

    # parameters considered are refined parameters
    ATTR_REF = ("fract_x", "fract_y", "fract_z", "occupancy", "u_iso_or_equiv",
                "u_equiv_geom_mean", "b_iso_or_equiv")
    ATTR_SIGMA = tuple([f"{_h:}_sigma" for _h in ATTR_REF])
    ATTR_CONSTR_FLAG = tuple([f"{_h:}_constraint" for _h in ATTR_REF])
    ATTR_REF_FLAG = tuple([f"{_h:}_refinement" for _h in ATTR_REF])
    ATTR_CONSTR_MARK = tuple([f"{_h:}_mark" for _h in ATTR_REF])

    # formats if cif format
    D_FORMATS = {'fract_x': "{:.6f}", 'fract_y': "{:.6f}",
                 'fract_z': "{:.6f}"}

    # constraints on the parameters
    D_CONSTRAINTS = {"adp_type": ["Uani", "Uiso", "Uovl", "Umpe", "Bani",
                                  "Biso", "Bovl"]}

    # default values for the parameters
    D_DEFAULT = {"occupancy": 1.}
    for key in ATTR_SIGMA:
        D_DEFAULT[key] = 0.
    for key in (ATTR_CONSTR_FLAG + ATTR_REF_FLAG):
        D_DEFAULT[key] = False
    for key in ATTR_CONSTR_MARK:
        D_DEFAULT[key] = ""

    PREFIX = "atom_site"

    def __init__(self, **kwargs) -> NoReturn:
        super(AtomSite, self).__init__()

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

    def form_object(self) -> bool:
        if self.is_attribute("space_group_wyckoff"):
            space_group_wyckoff = self.space_group_wyckoff
            self.__dict__["multiplicity"] = space_group_wyckoff.multiplicity
            self.__dict__["wyckoff_symbol"] = space_group_wyckoff.letter

            self.__dict__["fract_x_constraint"] = False
            self.__dict__["fract_y_constraint"] = False
            self.__dict__["fract_z_constraint"] = False
            xyz = numpy.array([self.fract_x, self.fract_y, self.fract_z],
                              dtype=float)

            r = space_group_wyckoff.r
            if r[0, 0] == 0:
                self.__dict__["fract_x_constraint"] = True
                self.__dict__["fract_x_refinement"] = False
            if r[1, 1] == 0:
                self.__dict__["fract_y_constraint"] = True
                self.__dict__["fract_y_refinement"] = False
            if r[2, 2] == 0:
                self.__dict__["fract_z_constraint"] = True
                self.__dict__["fract_z_refinement"] = False

            xyz_new = space_group_wyckoff.give_default_xyz(xyz)
            self.__dict__["fract_x"] = float(xyz_new[0])
            self.__dict__["fract_y"] = float(xyz_new[1])
            self.__dict__["fract_z"] = float(xyz_new[2])

        type_n = self.type_symbol
        self.__dict__["scat_length_neutron"] = get_scat_length_neutron(type_n)

    def define_space_group_wyckoff(self, space_group_wyckoff_l) -> bool:
        """Define space group by Wyckoff."""
        fract_x, fract_y, fract_z = self.fract_x, self.fract_y, self.fract_z
        x, y, z = float(fract_x), float(fract_y), float(fract_z)
        _id = space_group_wyckoff_l.get_id_for_fract(x, y, z)
        self.__dict__["space_group_wyckoff"] = space_group_wyckoff_l[_id]

    def apply_constraints(self, space_group_wyckoff_l) -> bool:
        """Apply constraints."""
        flag = not(self.is_attribute("space_group_wyckoff"))
        if flag:
            self.define_space_group_wyckoff(space_group_wyckoff_l)
        if self.is_defined():
            self.form_object()

    def calc_sc_chi(self, space_group, cell):
        """Calc sc_chi
        """
        # if self.is_attribute("sc_chi"):
        #     return self.sc_chi
        unit_cell_parameters = cell.get_unit_cell_parameters()
        x, y, z = self.fract_x, self.fract_y, self.fract_z
        o_11, o_12, o_13, o_21, o_22, o_23, o_31, o_32, o_33, o_3, o_2, o_3 =\
            space_group.calc_el_symm_for_xyz(x,y,z)
        zero_val = numpy.zeros(o_11.shape, dtype=int)
        symm_elems = numpy.array([zero_val, zero_val, zero_val, zero_val,
            o_11, o_12, o_13, o_21, o_22, o_23, o_31, o_32, o_33], dtype=int)

        sc_chi = calc_sc_chi(symm_elems, unit_cell_parameters, flag_unit_cell_parameters=False)[0]
        self.__dict__["sc_chi"] = sc_chi
        return sc_chi

    def calc_constr_number(self, space_group):
        """
        According to table 1 in Peterse, Palm, Acta Cryst.(1966), 20, 147.
        """
        if self.is_attribute("constr_number"):
            return self.constr_number
        x, y, z = self.fract_x, self.fract_y, self.fract_z
        o_11, o_12, o_13, o_21, o_22, o_23, o_31, o_32, o_33, o_3, o_2, o_3 =\
            space_group.calc_el_symm_for_xyz(x,y,z)
        b_11, b_22, b_33, b_12, b_13, b_23 = 107, 181, 41, 7, 19, 1
        i_11 = (o_11*o_11*b_11 + o_12*o_12*b_22 + o_13*o_13*b_33 +
                o_11*o_12*b_12 + o_11*o_13*b_13 + o_12*o_13*b_23 +
                o_12*o_11*b_12 + o_13*o_11*b_13 + o_13*o_12*b_23)
        i_22 = (o_21*o_21*b_11 + o_22*o_22*b_22 + o_23*o_23*b_33 +
                o_21*o_22*b_12 + o_21*o_23*b_13 + o_22*o_23*b_23 +
                o_22*o_21*b_12 + o_23*o_21*b_13 + o_23*o_22*b_23)
        i_33 = (o_31*o_31*b_11 + o_32*o_32*b_22 + o_33*o_33*b_33 +
                o_31*o_32*b_12 + o_31*o_33*b_13 + o_32*o_33*b_23 +
                o_32*o_31*b_12 + o_33*o_31*b_13 + o_33*o_32*b_23)
        i_12 = (o_11*o_21*b_11 + o_12*o_22*b_22 + o_13*o_23*b_33 +
                o_11*o_22*b_12 + o_11*o_23*b_13 + o_12*o_23*b_23 +
                o_12*o_21*b_12 + o_13*o_21*b_13 + o_13*o_22*b_23)
        i_13 = (o_11*o_31*b_11 + o_12*o_32*b_22 + o_13*o_33*b_33 +
                o_11*o_32*b_12 + o_11*o_33*b_13 + o_12*o_33*b_23 +
                o_12*o_31*b_12 + o_13*o_31*b_13 + o_13*o_32*b_23)
        i_23 = (o_21*o_31*b_11 + o_22*o_32*b_22 + o_23*o_33*b_33 +
                o_21*o_32*b_12 + o_21*o_33*b_13 + o_22*o_33*b_23 +
                o_22*o_31*b_12 + o_23*o_31*b_13 + o_23*o_32*b_23)
        r_11, r_22, r_33, r_12, r_13, r_23 = i_11.sum(), i_22.sum(),\
            i_33.sum(), i_12.sum(), i_13.sum(), i_23.sum()
        if r_13 == 0:
            if r_23 == 0:
                if r_12 == 0:
                    if r_11 == r_22:
                        if r_22 == r_33:
                            numb = 17
                        else:
                            numb = 8
                    elif r_22 == r_33:
                        numb = 12
                    else:
                        numb = 4
                elif r_11 == r_22:
                    if r_22 == 2*r_12:
                        numb = 16
                    else:
                        numb = 5
                elif r_22 == 2*r_12:
                    numb = 14
                else:
                    numb = 2
            elif r_22 == r_33:
                numb = 9
            else:
                numb = 3
        elif r_23 == 0:
            if r_22 == 2*r_12:
                numb = 13
            else:
                numb = 1
        elif r_23 == r_13:
            if r_22 == r_33:
                numb = 18
            else:
                numb = 6
        elif r_12 == r_13:
            numb = 10
        elif r_11 == 22:
            numb = 7
        elif r_22 == r_33:
            numb = 11
        elif r_22 == 2*r_12:
            numb = 15
        else:
            numb = 0  # no constraint
        self.__dict__["constr_number"] = numb
        return numb

    def report(self) -> str:
        """Report."""
        s_out = ""
        if ((self.is_attribute("scat_length_neutron")) &
                (self.is_attribute("label"))):
            s_out = \
                f'|{self.label.rjust(10):} | {self.scat_length_neutron: .3f}|'
        return s_out

    def get_flags_atom_fract_xyz(self):
        res = numpy.array([self.fract_x_refinement, self.fract_y_refinement, self.fract_z_refinement], dtype=bool)
        return res

    def get_flags_atom_occupancy(self):
        res = self.occupancy_refinement
        return res

    def get_flags_atom_b_iso(self):
        res = False
        if not(self.is_attribute("adp_type")):
            res = False
        elif self.adp_type in ["Uiso", "Uovl", "Umpe", ]:
            res = self.u_iso_or_equiv_refinement
        elif self.adp_type in ["Biso", "Bovl", ]:
            res = self.b_iso_or_equiv_refinement
        return res


class AtomSiteL(LoopN):
    """Describe the atom sites in crystal.

    Attributes
    ----------
        - label, type_symbol, fract_x, fract_y, fract_z (mandatory)
        - occupancy, adp_type, u_iso_or_equiv, u_equiv_geom_mean,
          b_iso_or_equiv, multiplicity, wyckoff_symbol, cartn_x, cartn_y,
          cartn_z (optional)
        - scat_length_neutron (internal)
        - space_group_wyckoff, constr_number (internal protected)
    """
    ITEM_CLASS = AtomSite
    ATTR_INDEX = "label"

    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(AtomSiteL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name

    def apply_constraints(self, space_group_wyckoff_l) -> bool:
        """Apply constraints."""
        for item in self.items:
            item.apply_constraints(space_group_wyckoff_l)

    def calc_constr_number(self, space_group):
        """Calc constr number."""
        l_numb = [item.calc_constr_number(space_group) for item in self.items]
        return l_numb

    def report(self) -> str:
        """Report."""
        ls_out = ["# Scattering amplitude"]
        ls_out.append("| type| b_scat|")
        ls_out.append("|-----|-------|")
        ls_out.extend([item.report() for item in self.items])
        return "\n".join(ls_out)

    def get_flags_atom_fract_xyz(self):
        flags_atom_fract_xyz = numpy.stack([item.get_flags_atom_fract_xyz() for item in self.items], axis=0).transpose()
        return flags_atom_fract_xyz

    def get_flags_atom_occupancy(self):
        flags_atom_occupancy = numpy.stack([item.get_flags_atom_occupancy() for item in self.items], axis=0)
        return flags_atom_occupancy

    def get_flags_atom_b_iso(self):
        flags_atom_b_iso = numpy.stack([item.get_flags_atom_b_iso() for item in self.items], axis=0)
        return flags_atom_b_iso


# s_cont = """
#  loop_
#  _atom_site_label
#  _atom_site_type_symbol
#  _atom_site_fract_x
#  _atom_site_fract_y
#  _atom_site_fract_z
#  _atom_site_adp_type
#  _atom_site_B_iso_or_equiv
#  _atom_site_occupancy
#   Fe3A   Fe  0.12500 0.12500 0.12500  Uani   0.0   1.0
#   Fe3B   Fe  0.50000 0.50000 0.50000  Uani   0.0   1.0
#   O1     O   0.25521 0.25521 0.25521  Uiso   0.0   1.0
#   """

# obj = AtomSiteL.from_cif(s_cont)
# print(obj, end="\n\n")
# print(obj["O1"].scat_length_neutron, end="\n\n")
