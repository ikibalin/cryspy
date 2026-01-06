"""AtomSiteMoment, AtomSiteMomentL classes are given."""
import numpy
from typing import NoReturn

from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary

from cryspy.A_functions_base.matrix_operations import calc_m_v

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN
from cryspy.C_item_loop_classes.cl_1_cell import Cell
from cryspy.A_functions_base.unit_cell import calc_m_m_norm_by_unit_cell_parameters

class AtomSiteMoment(ItemN):
    """AtomSiteMoment class.

    This category provides a loop for presenting the magnetic moments
    of atoms in one of several coordinate systems. This is a child category
    of the AtomSite category, so that the magnetic moments can either be
    listed alongside the non-magnetic atom properties in the main AtomSite loop
    (not realized) or be listed in a separate loop (realized)

    Attributes
    ----------
        - label, symmform (mandatory)
        - cartn_x, cartn_y, cartn_z, crystalaxis_x, crystalaxis_y,
          crystalaxis_z, modulation_flag, refinement_flags_magnetic,
          spherical_azimuthal, spherical_modulus, spherical_polar
    """

    ATTR_MANDATORY_NAMES = ("label",)
    ATTR_MANDATORY_TYPES = (str, )
    ATTR_MANDATORY_CIF = ("label", )

    ATTR_OPTIONAL_NAMES = (
        "cartn_x", "cartn_y", "cartn_z", "crystalaxis_x", "crystalaxis_y",
        "crystalaxis_z",  "symmform", "modulation_flag", "refinement_flags_magnetic",
        "spherical_azimuthal", "spherical_modulus", "spherical_polar")
    ATTR_OPTIONAL_TYPES = (float, float, float, float, float, float, str,
                           str, str, float, float, float)
    ATTR_OPTIONAL_CIF = (
        "cartn_x", "cartn_y", "cartn_z", "crystalaxis_x", "crystalaxis_y",
        "crystalaxis_z",  "symmform", "modulation_flag", "refinement_flags_magnetic",
        "spherical_azimuthal", "spherical_modulus", "spherical_polar")

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ()
    ATTR_INT_PROTECTED_NAMES = ()

    # parameters considered are refined parameters
    ATTR_REF = ("crystalaxis_x", "crystalaxis_y", "crystalaxis_z")
    ATTR_SIGMA = tuple([f"{_h:}_sigma" for _h in ATTR_REF])
    ATTR_CONSTR_FLAG = tuple([f"{_h:}_constraint" for _h in ATTR_REF])
    ATTR_REF_FLAG = tuple([f"{_h:}_refinement" for _h in ATTR_REF])
    ATTR_CONSTR_MARK = tuple([f"{_h:}_mark" for _h in ATTR_REF])

    # constraints on the parameters
    D_CONSTRAINTS = {"modulation_flag": ["yes", "y", "no", "n"],
                     "refinement_flags_magnetic": ["S", "M", "A", "SM", "SA",
                                                   "MA", "SMA"]}

    # default values for the parameters
    D_DEFAULT = {}
    for key in ATTR_SIGMA:
        D_DEFAULT[key] = 0.
    for key in (ATTR_CONSTR_FLAG + ATTR_REF_FLAG):
        D_DEFAULT[key] = False
    for key in ATTR_CONSTR_MARK:
        D_DEFAULT[key] = ""

    PREFIX = "atom_site_moment"

    def __init__(self, **kwargs) -> NoReturn:
        super(AtomSiteMoment, self).__init__()

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

    def calc_moment_in_xyz(self, cell:Cell):
        """Calculate moment.
        """
        unit_cell_parameters = cell.get_unit_cell_parameters()
        m_norm = calc_m_m_norm_by_unit_cell_parameters(unit_cell_parameters, flag_unit_cell_parameters=False)[0]
        m_dn = numpy.array([self.crystalaxis_x,
                            self.crystalaxis_y,
                            self.crystalaxis_z], dtype=float)
        m_xyz = calc_m_v(m_dn, m_norm, flag_m=False, flag_v=True)[0]
        return m_xyz

    def calc_moment(self, cell:Cell):
        """Calculate moment.
        """
        m_xyz = self.calc_moment_in_xyz(cell)
        np_moment = numpy.sqrt(numpy.square(m_xyz[0]) + numpy.square(m_xyz[1]) +
                               numpy.square(m_xyz[2]))
        return np_moment

    def calc_zeeman(self, field_cryst, cell:Cell):
        """Calc Zeeman energy.

        FIXME: IK I don't remember why I call it zeeman energy here as it looks like M_perp component
        field_cryst is given in XYZ crystal axes
        """
        h_a, h_b, h_c = field_cryst[0], field_cryst[1], field_cryst[2]

        if abs(h_a)+abs(h_b)+abs(h_c) == 0.:
            np_val = numpy.zeros(len(self.label), dtype=float)
            return np_val
        else:
            mod_h = (abs(h_a)**2+abs(h_b)**2+abs(h_c)**2)**0.5
            h_a_n, h_b_n, h_c_n = h_a/mod_h, h_b/mod_h, h_c/mod_h
            m_xyz = self.calc_moment_in_xyz(cell)

            np_moment_sq = self.calc_moment(cell)**2
            np_val = h_a_n*m_xyz[0]+h_b_n*m_xyz[1]+h_c_n*m_xyz[2]
            np_res = numpy.sqrt(np_moment_sq-numpy.square(np_val))
        return np_res

    def to_cif(self, separator: str = ".") -> str:
        return super(AtomSiteMoment, self).to_cif(separator=separator)

class AtomSiteMomentL(LoopN):
    """AtomSiteMomentL class.

    AtomSiteMomentL category provides a loop for presenting the magnetic
    moments  of atoms in one of several coordinate systems. This is a child
    category of the AtomSite category, so that the magnetic moments can either
    be listed alongside the non-magnetic atom properties in the main AtomSite
    loop (not realized) or be listed in a separate loop (realized)
    """

    ITEM_CLASS = AtomSiteMoment
    ATTR_INDEX = "label"

    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(AtomSiteMomentL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name

    def to_cif(self, separator: str = ".") -> str:
        return super(AtomSiteMomentL, self).to_cif(separator=separator)
    
    def calc_moment(self, cell):
        """Get moment in mu_B by unit cell parameters.
        """
        l_res = [item.calc_moment(cell) for item in self.items]
        return l_res
# s_cont = """
# loop_
# _atom_site_moment_label
# _atom_site_moment_crystalaxis_x
# _atom_site_moment_crystalaxis_y
# _atom_site_moment_crystalaxis_z
# _atom_site_moment_symmform
# Tm1_1   -6.44   0.0     0.0   mx,0,0
# Tm1_2   -6.44   0.0     0.0   mx,0,0
# Tm1_3   -6.44   -6.44   0.0   mx,my,0
#   """

# obj = AtomSiteMomentL.from_cif(s_cont)
# print(obj, end="\n\n")
# print(obj["Tm1_2"], end="\n\n")
