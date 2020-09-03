"""Description of classes Section, SectionL."""
from typing import NoReturn
import numpy
from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN

from cryspy.C_item_loop_classes.cl_1_cell import Cell
from cryspy.C_item_loop_classes.cl_1_atom_site import AtomSiteL
from cryspy.C_item_loop_classes.cl_1_space_group_symop import SpaceGroupSymop


class Section(ItemN):
    """
    Section class.

    Describe information concerning the density point.
    """

    ATTR_MANDATORY_NAMES = ("id", "size_x", "size_y", "atom_center",
                            "atom_along_axis_x", "atom_x_operation_xyz",
                            "atom_along_axis_y", "atom_y_operation_xyz",
                            "points_x", "points_y", "url_out")
    ATTR_MANDATORY_TYPES = (str, float, float, str, str, str, str, str, int,
                            int, str)
    ATTR_MANDATORY_CIF = ("id", "size_x", "size_y", "atom_center",
                          "atom_along_axis_x", "atom_x_operation_xyz",
                          "atom_along_axis_y", "atom_y_operation_xyz",
                          "points_x", "points_y", "url_out")

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

    # formats if cif format
    D_FORMATS = {"size_x": "{:.2f}", "size_y": "{:.2f}"}

    # constraints on the parameters
    D_CONSTRAINTS = {}

    # default values for the parameters
    D_DEFAULT = {}
    for key in ATTR_SIGMA:
        D_DEFAULT[key] = 0.
    for key in (ATTR_CONSTR_FLAG + ATTR_REF_FLAG):
        D_DEFAULT[key] = False

    PREFIX = "section"

    def __init__(self, **kwargs) -> NoReturn:
        super(Section, self).__init__()

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

    def calc_axes_x_y_z(self, cell: Cell, atom_site: AtomSiteL):
        """
        Calculate three vectors of axes: pos_x, pos_y, pos_z.

        Arguments
        ---------
            - cell
            - atom_site (loop)
        """
        atom_o_label = self.atom_center
        _item_a_s_o = atom_site[atom_o_label]
        fract_a_o_xyz = numpy.array([_item_a_s_o.fract_x, _item_a_s_o.fract_y,
                                     _item_a_s_o.fract_z], dtype=float)
        pos_a_o_xyz = cell.calc_position_by_coordinate(
            fract_a_o_xyz[0], fract_a_o_xyz[1], fract_a_o_xyz[2])
        pos_a_o_xyz = numpy.array(pos_a_o_xyz, dtype=float)

        atom_x_label = self.atom_along_axis_x
        atom_x_operation_xyz = self.atom_x_operation_xyz
        s_g_s_a_x = SpaceGroupSymop(operation_xyz=atom_x_operation_xyz)
        _item_a_s_x = atom_site[atom_x_label]
        fract_a_x_xyz = numpy.array([_item_a_s_x.fract_x, _item_a_s_x.fract_y,
                                     _item_a_s_x.fract_z], dtype=float)
        fract_a_x_xyz = numpy.dot(s_g_s_a_x.r, fract_a_x_xyz) + s_g_s_a_x.b
        pos_a_x_xyz = cell.calc_position_by_coordinate(
            fract_a_x_xyz[0], fract_a_x_xyz[1], fract_a_x_xyz[2])
        pos_a_x_xyz = numpy.array(pos_a_x_xyz, dtype=float)

        atom_y_label = self.atom_along_axis_y
        atom_y_operation_xyz = self.atom_y_operation_xyz
        s_g_s_a_y = SpaceGroupSymop(operation_xyz=atom_y_operation_xyz)
        _item_a_s_y = atom_site[atom_y_label]
        fract_a_y_xyz = numpy.array([_item_a_s_y.fract_x, _item_a_s_y.fract_y,
                                     _item_a_s_y.fract_z], dtype=float)
        fract_a_y_xyz = numpy.dot(s_g_s_a_y.r, fract_a_y_xyz) + s_g_s_a_y.b
        pos_a_y_xyz = cell.calc_position_by_coordinate(
            fract_a_y_xyz[0], fract_a_y_xyz[1], fract_a_y_xyz[2])
        pos_a_y_xyz = numpy.array(pos_a_y_xyz, dtype=float)

        v_pos_x = (pos_a_x_xyz-pos_a_o_xyz)/((
            numpy.square(pos_a_x_xyz-pos_a_o_xyz)).sum())**0.5
        v_pos_a_y = (pos_a_y_xyz-pos_a_o_xyz)/((
            numpy.square(pos_a_y_xyz-pos_a_o_xyz)).sum())**0.5
        v_pos_y_not_norm = v_pos_a_y - (v_pos_a_y * v_pos_x).sum() * v_pos_x
        v_pos_y = v_pos_y_not_norm / (
            numpy.square(v_pos_y_not_norm).sum())**0.5

        v_pos_z = numpy.cross(v_pos_x, v_pos_y)

        return v_pos_x, v_pos_y, v_pos_z

    def calc_fractions(self, cell: Cell, atom_site: AtomSiteL) -> \
            numpy.ndarray:
        """
        Give a numpy.nd_array of fractions: fractions_xyz.

        Arguments
        ---------
            - cell
            - atom_site
        """
        size_x, size_y = self.size_x, self.size_y
        points_x, points_y = self.points_x, self.points_y

        atom_o_label = self.atom_center
        _item_a_s_o = atom_site[atom_o_label]
        fract_a_o_xyz = numpy.array([_item_a_s_o.fract_x, _item_a_s_o.fract_y,
                                     _item_a_s_o.fract_z], dtype=float)
        pos_a_o_xyz = cell.calc_position_by_coordinate(
            fract_a_o_xyz[0], fract_a_o_xyz[1], fract_a_o_xyz[2])
        pos_a_o_xyz = numpy.array(pos_a_o_xyz, dtype=float)

        v_pos_x, v_pos_y, v_pos_z = self.calc_axes_x_y_z(cell, atom_site)

        v_delta_pos_x = v_pos_x * size_x / float(points_x)
        v_delta_pos_y = v_pos_y * size_y / float(points_y)
        np_x_2d, np_y_2d = numpy.meshgrid(range(-points_x//2, points_x//2),
                                          range(-points_y//2, points_y//2),
                                          indexing="ij")
        np_x_1d, np_y_1d = np_x_2d.flatten(), np_y_2d.flatten()

        pos_xyz = np_x_1d[numpy.newaxis, :]*v_delta_pos_x[:, numpy.newaxis] + \
            np_y_1d[numpy.newaxis, :] * v_delta_pos_y[:, numpy.newaxis] + \
            pos_a_o_xyz[:, numpy.newaxis]
        fract_x, fract_y, fract_z = cell.calc_coordinate_by_position(
            pos_xyz[0, :], pos_xyz[1, :], pos_xyz[2, :])
        return fract_x, fract_y, fract_z


class SectionL(LoopN):
    """SectionL class.

    Describe information concerning the density point.
    """

    ITEM_CLASS = Section
    ATTR_INDEX = "id"

    def __init__(self, loop_name=None) -> NoReturn:
        super(SectionL, self).__init__()
        self.__dict__["items"] = []
        self.__dict__["loop_name"] = loop_name


# s_cont = """
# loop_
# _section_id
# _section_size_x
# _section_size_y
# _section_atom_center
# _section_atom_along_axis_x
# _section_atom_along_axis_x_symop
# _section_atom_along_axis_y
# _section_atom_along_axis_y_symop
# _section_points_x
# _section_points_y
# _section_url_out
#   core 5 5 Ho1 O1 x,y,z O2 x,y,z 100 100 s_ho1.dat
#   """

# obj = SectionL.from_cif(s_cont)
# print(obj, end="\n\n")
# print(obj["core"], end="\n\n")
