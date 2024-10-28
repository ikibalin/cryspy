"""DiffrnOrientMatrix and DiffrnOrientMatrixL classes."""
from typing import NoReturn
import numpy

from cryspy.A_functions_base.function_1_matrices import\
    calc_product_matrices, calc_product_matrix_vector
from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN

from cryspy.C_item_loop_classes.cl_1_cell import Cell


class DiffrnOrientMatrix(ItemN):
    """UB-matrix desription.

    Data items in the DIFFRN_ORIENT_MATRIX category record details
    about the orientation matrix used in the measurement of the
    diffraction intensities.

    By default the UB matrix are given in ccsl notation

    Attributes
    ----------
        - ub_11, ub_12, ub_13, ub_21, ub_22, ub_23, ub_31, ub_32, ub_33
          (mandatory)
        - "id", "type" (optional)
        - u_11, u_12, u_13, u_21, u_22, u_23, u_31, u_32, u_33 (internal)
        - cell (protected)

    Constraints
    -----------
        - type: "CCSL", "6T2@LLB", "5C1@LLB", "BusingLevy", "-YXZ", "X-YZ",
                "XYZ"
    """

    ATTR_MANDATORY_NAMES = ("ub_11", "ub_12", "ub_13", "ub_21", "ub_22",
                            "ub_23", "ub_31", "ub_32", "ub_33")
    ATTR_MANDATORY_TYPES = (float, float, float, float, float, float,
                            float, float, float)
    ATTR_MANDATORY_CIF = ("ub_11", "ub_12", "ub_13", "ub_21", "ub_22", "ub_23",
                          "ub_31", "ub_32", "ub_33")

    ATTR_OPTIONAL_NAMES = ("id", "type")
    ATTR_OPTIONAL_TYPES = (str, str)
    ATTR_OPTIONAL_CIF = ("id", "type")

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ("u_11", "u_12", "u_13", "u_21", "u_22", "u_23",
                      "u_31", "u_32", "u_33", "ub", "u")
    ATTR_INT_PROTECTED_NAMES = ("cell", )

    # parameters considered are refined parameters
    ATTR_REF = ()
    ATTR_SIGMA = tuple([f"{_h:}_sigma" for _h in ATTR_REF])
    ATTR_CONSTR_FLAG = tuple([f"{_h:}_constraint" for _h in ATTR_REF])
    ATTR_REF_FLAG = tuple([f"{_h:}_refinement" for _h in ATTR_REF])
    ATTR_CONSTR_MARK = tuple([f"{_h:}_mark" for _h in ATTR_REF])

    # constraints on the parameters
    D_CONSTRAINTS = {"type": ["CCSL", "6T2@LLB", "5C1@LLB", "BusingLevy",
                              "-YXZ", "X-YZ", "XYZ"]}

    # default values for the parameters
    D_DEFAULT = {"type": "CCSL"}
    for key in ATTR_SIGMA:
        D_DEFAULT[key] = 0.
    for key in (ATTR_CONSTR_FLAG + ATTR_REF_FLAG):
        D_DEFAULT[key] = False
    for key in ATTR_CONSTR_MARK:
        D_DEFAULT[key] = ""

    PREFIX = "diffrn_orient_matrix"

    def __init__(self, **kwargs) -> NoReturn:
        super(DiffrnOrientMatrix, self).__init__()

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

    def __recalc_u_cell(self):

        ub_ccsl = self.calc_ub_ccsl()

        v_b_1 = ub_ccsl[:, 0]
        v_b_2 = ub_ccsl[:, 1]
        v_b_3 = ub_ccsl[:, 2]

        b_1 = float(numpy.sqrt(((v_b_1*v_b_1).sum())))
        b_2 = float(numpy.sqrt(((v_b_2*v_b_2).sum())))
        b_3 = float(numpy.sqrt(((v_b_3*v_b_3).sum())))

        c_i_g = (v_b_1*v_b_2).sum()/(b_1*b_2)
        c_i_b = (v_b_1*v_b_3).sum()/(b_1*b_3)
        c_i_a = (v_b_2*v_b_3).sum()/(b_2*b_3)

        i_a, i_b = numpy.arccos(c_i_a), numpy.arccos(c_i_b)
        i_g = numpy.arccos(c_i_g)

        s_i_a, s_i_b, s_i_g = numpy.sin(i_a), numpy.sin(i_b), numpy.sin(i_g)

        i_vol = b_1*b_2*b_3*(1-c_i_a**2-c_i_b**2-c_i_g**2+2.*c_i_a*c_i_b*c_i_g
                             )**0.5

        c_a = (c_i_b*c_i_g-c_i_a)/(s_i_b*s_i_g)
        c_b = (c_i_a*c_i_g-c_i_b)/(s_i_a*s_i_g)
        c_g = (c_i_a*c_i_b-c_i_g)/(s_i_a*s_i_b)

        a_1 = b_2*b_3*s_i_a/i_vol
        a_2 = b_1*b_3*s_i_b/i_vol
        a_3 = b_1*b_2*s_i_g/i_vol

        alpha, beta = numpy.arccos(c_a), numpy.arccos(c_b)
        gamma = numpy.arccos(c_g)

        if not(self.is_attribute("cell")):
            cell = Cell(length_a=a_1, length_b=a_2, length_c=a_3,
                        angle_alpha=float(alpha*180./numpy.pi),
                        angle_beta=float(beta*180./numpy.pi),
                        angle_gamma=float(gamma*180./numpy.pi))
            self.__dict__["cell"] = cell
        else:
            cell = self.cell
            cell.length_a, cell.length_b, cell.length_c = a_1, a_2, a_3
            cell.angle_alpha = float(alpha*180./numpy.pi)
            cell.angle_beta = float(beta*180./numpy.pi)
            cell.angle_gamma = float(gamma*180./numpy.pi)
        cell.form_object
        m_b = cell.m_b
        m_ib = numpy.linalg.inv(m_b)
        u_ccsl = numpy.dot(ub_ccsl,  m_ib)

        self.__dict__["u_11"] = float(u_ccsl[0, 0])
        self.__dict__["u_12"] = float(u_ccsl[0, 1])
        self.__dict__["u_13"] = float(u_ccsl[0, 2])
        self.__dict__["u_21"] = float(u_ccsl[1, 0])
        self.__dict__["u_22"] = float(u_ccsl[1, 1])
        self.__dict__["u_23"] = float(u_ccsl[1, 2])
        self.__dict__["u_31"] = float(u_ccsl[2, 0])
        self.__dict__["u_32"] = float(u_ccsl[2, 1])
        self.__dict__["u_33"] = float(u_ccsl[2, 2])

    def calc_ub_ccsl(self):
        """Calculate UB Matrix in CCSL notation.

        Based on the given UB matrix with defined type variable.
        """
        s_notation = self.type
        a1, a2, a3 = self.ub_11, self.ub_12, self.ub_13
        b1, b2, b3 = self.ub_21, self.ub_22, self.ub_23
        c1, c2, c3 = self.ub_31, self.ub_32, self.ub_33

        ub_ccsl = numpy.array([[a1, a2, a3],
                               [b1, b2, b3],
                               [c1, c2, c3]], dtype=float)
        if s_notation in ["XYZ", "6T2@LLB"]:
            ub_ccsl = numpy.array([[a1, a2, a3],
                                   [b1, b2, b3],
                                   [c1, c2, c3]], dtype=float)
        elif s_notation in ["-YXZ", "5C1@LLB"]:
            ub_ccsl = numpy.array([[-b1, -b2, -b3],
                                   [a1, a2, a3],
                                   [c1, c2, c3]], dtype=float)
            # ub_ccsl = numpy.array([[-b1,-b2,-b3],
            #                       [ a1, a2, a3],
            #                       [ c1, c2, c3]], dtype=float)
        elif s_notation in ["Y-XZ", "BusingLevy"]:
            ub_ccsl = numpy.array([[b1, b2, b3],
                                   [-a1, -a2, -a3],
                                   [c1, c2, c3]], dtype=float)
        return ub_ccsl

    def calc_matrix_from_ccsl(self, matrix_ccsl):
        """Calculate UB Matrix (or U Matrix).

        Calculations are done in the notation defined by type variable
        based on the UB Matrix (or U Matrix) in CCSL notation.
        """
        s_notation = self.get_notation()
        a1, a2, a3 = matrix_ccsl[0, 0], matrix_ccsl[0, 1], matrix_ccsl[0, 2]
        b1, b2, b3 = matrix_ccsl[1, 0], matrix_ccsl[1, 1], matrix_ccsl[1, 2]
        c1, c2, c3 = matrix_ccsl[2, 0], matrix_ccsl[2, 1], matrix_ccsl[2, 2]

        # the matrices are the same as in self.calc_ub_ccsl()
        # because Op. x Op. = 1
        ub = numpy.array([[a1, a2, a3],
                          [b1, b2, b3],
                          [c1, c2, c3]], dtype=float)
        if s_notation == "6T2@LLB":
            ub = numpy.array([[a1, a2, a3],
                              [b1, b2, b3],
                              [c1, c2, c3]], dtype=float)
        elif s_notation == "5C1@LLB":
            ub = numpy.array([[-b1, -b2, -b3],
                              [a1, a2, a3],
                              [c1, c2, c3]], dtype=float)
            # ub = numpy.array([[-b1, -b2, -b3],
            #                  [ a1, a2, a3],
            #                  [ c1, c2, c3]], dtype=float)
        elif s_notation == "BusingLevy":
            ub = numpy.array([[b1, b2, b3],
                              [-a1, -a2, -a3],
                              [c1, c2, c3]], dtype=float)
        return ub

    @classmethod
    def calc_ub(cls, u, type="CCSL",
                abc_angles=(6.8, 6.8, 6.8, 90., 90., 90.)):
        """Calculate UB matrix by U matrix.

        U matrix is given in the notation defined by type
        variable and cell parameters (a, b, c, alpha, beta, gamma)
        """
        (a, b, c, al, be, ga) = abc_angles[:6]
        cell = Cell(length_a=a, length_b=b, length_c=c, angle_alpha=al,
                    angle_beta=be, angle_gamma=ga)
        m_b = cell.m_b
        obj = cls(type=type)

        # in general case it should not be correct as Op. x Op.^-1 != 1
        u_ccsl = obj.calc_matrix_from_ccsl(u)
        ub_ccsl = numpy.dot(u_ccsl, m_b)
        ub = obj.calc_matrix_from_ccsl(ub_ccsl)

        obj.ub_11 = ub[0, 0]
        obj.ub_12 = ub[0, 1]
        obj.ub_13 = ub[0, 2]
        obj.ub_21 = ub[1, 0]
        obj.ub_22 = ub[1, 1]
        obj.ub_23 = ub[1, 2]
        obj.ub_31 = ub[2, 0]
        obj.ub_32 = ub[2, 1]
        obj.ub_33 = ub[2, 2]
        obj.form_object()
        return obj

    def form_object(self) -> bool:
        """Redefine parent form object."""
        flag = False
        if self.is_defined():
            ub = numpy.array([[self.ub_11, self.ub_12, self.ub_13],
                              [self.ub_21, self.ub_22, self.ub_23],
                              [self.ub_31, self.ub_32, self.ub_33]],
                             dtype=float)
            self.__dict__["ub"] = ub
            self.__recalc_u_cell()

            u = numpy.array([[self.u_11, self.u_12, self.u_13],
                             [self.u_21, self.u_22, self.u_23],
                             [self.u_31, self.u_32, self.u_33]],
                            dtype=float)
            self.__dict__["u"] = u
            flag = True
        return flag

    def calc_e_up(self, phi=0., omega=0., chi=0.):
        """
        phi, omega, chi = angles of detector in radians.

        They are zeros by defaults.
        """
        orientation_ij = (self.u_11, self.u_12, self.u_13,
                          self.u_21, self.u_22, self.u_23,
                          self.u_31, self.u_32, self.u_33)

        phi_ij = (numpy.cos(phi), numpy.sin(phi), 0.*phi,
                  -numpy.sin(phi), numpy.cos(phi), 0.*phi,
                  0.*phi, 0.*phi, 0.*phi+1.)

        omega_ij = (numpy.cos(omega), numpy.sin(omega), 0.*omega,
                    -numpy.sin(omega), numpy.cos(omega), 0.*omega,
                    0.*omega, 0.*omega, 0.*omega+1.)

        chi_ij = (numpy.cos(chi), 0.*chi, numpy.sin(chi),
                  0.*chi, 0.*chi+1., 0.*chi,
                  -numpy.sin(chi), 0.*chi, numpy.cos(chi))

        u_11, u_12, u_13, u_21, u_22, u_23, u_31, u_32, u_33 = \
            calc_product_matrices(omega_ij, chi_ij, phi_ij, orientation_ij)
        ut_ij = (u_11, u_21, u_31, u_12, u_22, u_32, u_13, u_23, u_33)
        e_up_1, e_up_2, e_up_3 = calc_product_matrix_vector(
            ut_ij, (0., 0., 1.))
        return numpy.array([e_up_1, e_up_2, e_up_3], dtype=float)

    def calc_angle(self, index_h, index_k, index_l, wavelength: float = 1.4,
            diffracted_beam: str = "left", diffractometer_axes: str = "anticlockwise"):
        """Calculate scattering angles for given reflection hkl.

        Output
        ------
            gamma  is the azimuthal angle (degrees)
            nu is the elevation angle in the laboratory coordinate system
                (xyz), where x||k_i, z||H. (degrees)
            phi is ... (degrees)
        """
        ub = self.ub
        q_ub = numpy.dot(ub, [index_h, index_k, index_l])
        q = numpy.linalg.norm(q_ub)
        gamma_0 = 0.
        nu_0 = 0.

        q_final = [- q**2 * wavelength / 2.0, 0, q_ub[2]]
        q_y2 = q**2 - q_final[0]**2 - q_final[2]**2
        flag_hh_2 = True

        if q_y2 < 0:
            # print("Angles are not found.")
            flag_hh_2 = False
        if flag_hh_2:
            q_final[1] = numpy.sqrt(q_y2)
            k_f = [q_final[0] + 1.0/wavelength, q_final[1], q_final[2]]
            nu = numpy.arcsin(k_f[2]*wavelength)/numpy.pi*180. - nu_0

            gamma = numpy.arctan2(k_f[1], k_f[0]) / numpy.pi*180. - gamma_0
            phi = (numpy.arctan2(q_final[1], q_final[0]) -
                   numpy.arctan2(q_ub[1], q_ub[0]))/numpy.pi*180
            [phi, gamma, nu] = [phi if phi > 0. else phi + 360., gamma if
                                gamma > 0. else gamma + 360., nu]
            if ((diffracted_beam.strip().lower() == "left") and
                (diffractometer_axes.strip().lower() =="anticlockwise")): # left, anticlockwise
                pass
            elif ((diffracted_beam.strip().lower() != "left") and
                (diffractometer_axes.strip().lower() != "anticlockwise")): # right, clockwise
                pass
            elif ((diffracted_beam.strip().lower() == "left") and
                (diffractometer_axes.strip().lower() !="anticlockwise")): # right, anticlockwise
                gamma = -gamma
                phi = -phi
            elif ((diffracted_beam.strip().lower() != "left") and
                (diffractometer_axes.strip().lower() == "anticlockwise")): # left, clockwise
                gamma = -gamma
                phi = -phi

            return gamma, nu, phi

    def calc_q2(self, index_h: numpy.ndarray, index_k: numpy.ndarray,
                index_l: numpy.ndarray, cell: Cell) -> numpy.ndarray:
        r"""Calculate q2 = (sin \alpha)^2.

        \alpha is angle between magnetic field and scattering vector (?)
        """
        k_1, k_2, k_3 = cell.calc_k_loc(index_h, index_k, index_l)
        h_1_loc, h_2_loc = float(self.u_31), float(self.u_32)
        h_3_loc = float(self.u_33)
        res = 1.-(k_1*h_1_loc + k_2*h_2_loc + k_3*h_3_loc)**2
        return res

    def report(self):
        s_out = f"""Orientation matrix U:
          a* [c,  a*]        c
 X: {self.u_11: 8.5f} {self.u_12: 8.5f} {self.u_13: 8.5f}
 Y: {self.u_21: 8.5f} {self.u_22: 8.5f} {self.u_23: 8.5f}
 Z: {self.u_31: 8.5f} {self.u_32: 8.5f} {self.u_33: 8.5f}

axis 'X' is along incident beam;
axis 'Z' is vertical direction."""
        return s_out

    def report_html(self):
        s_out = f"""<b>Orientation matrix U:</b><br>
        <table>
<tr><th></th><th>a*</th><th>[c,  a*]</th><th>c</th></tr>
<tr><th>X:</th><td>{self.u_11: 8.5f}</td><td>{self.u_12: 8.5f}</td><td>{self.u_13: 8.5f}</td></tr>
<tr><th>Y:</th><td>{self.u_21: 8.5f}</td><td>{self.u_22: 8.5f}</td><td>{self.u_23: 8.5f}</td></tr>
<tr><th>Z:</th><td>{self.u_31: 8.5f}</td><td>{self.u_32: 8.5f}</td><td>{self.u_33: 8.5f}</td></tr>
        </table>
axis 'X' is along incident beam;<br>
axis 'Z' is vertical direction."""
        return s_out

    def get_dictionary(self):
        res = {}
        u_matrix = self.u
        e_up = numpy.array([u_matrix[2,0], u_matrix[2,1], u_matrix[2,2]], dtype=float)

        res["matrix_u"] = numpy.array([
            self.u_11, self.u_12, self.u_13,
            self.u_21, self.u_22, self.u_23,
            self.u_31, self.u_32, self.u_33], dtype = float)
        res["matrix_ub"] = numpy.array([
            self.ub_11, self.ub_12, self.ub_13,
            self.ub_21, self.ub_22, self.ub_23,
            self.ub_31, self.ub_32, self.ub_33], dtype = float)
        return res

class DiffrnOrientMatrixL(LoopN):
    """Description of DiffrnOrientMatrixL in loop."""

    ITEM_CLASS = DiffrnOrientMatrix
    ATTR_INDEX = "id"

    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(DiffrnOrientMatrixL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name

# s_cont = """
#  _diffrn_orient_matrix_UB_11           -0.04170
#  _diffrn_orient_matrix_UB_12           -0.01429
#  _diffrn_orient_matrix_UB_13           -0.02226
#  _diffrn_orient_matrix_UB_21           -0.00380
#  _diffrn_orient_matrix_UB_22           -0.05578
#  _diffrn_orient_matrix_UB_23           -0.05048
#  _diffrn_orient_matrix_UB_31            0.00587
#  _diffrn_orient_matrix_UB_32           -0.13766
#  _diffrn_orient_matrix_UB_33            0.02277
#   """

# obj = DiffrnOrientMatrix.from_cif(s_cont)
# print(obj, end="\n\n")
# print(obj.cell, end="\n\n")
