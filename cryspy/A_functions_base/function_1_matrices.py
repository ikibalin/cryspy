"""Mactices functions.

Functions
---------
    - calc_chi_sq
    - transform_string_to_r_b
    - transform_string_to_digits
    - transform_fraction_with_label_to_string
    - transform_digits_to_string
    - transform_r_b_to_string
    - calc_mRmCmRT
    - calc_rotation_matrix_ij_by_euler_angles
    - calc_euler_angles_by_rotation_matrix_ij
    - calc_determinant_matrix_ij
    - calc_inverse_matrix_ij
    - calc_rotation_matrix_ij_around_axis
    - calc_product_matrices
    - calc_product_matrix_vector
    - calc_vector_angle
    - calc_vector_product
    - scalar_product
    - calc_rotation_matrix_by_two_vectors
    - tri_linear_interpolation
    - ortogonalize_matrix
    - calc_moment_2d_by_susceptibility
    - calc_phase_3d
"""
import numpy
from fractions import Fraction
from typing import Tuple


def calc_chi_sq(y_exp, y_error, y_model, der_y=None):
    """Calc chi_sq."""
    y_diff = (y_exp - y_model)/y_error
    chi_sq = (numpy.square(numpy.abs(y_diff))).sum()
    if der_y is not None:
        der_chi_sq = (-2.*(y_diff/y_error)[:, numpy.newaxis]*der_y).sum(axis=0)
    else:
        der_chi_sq = None
    return chi_sq, der_chi_sq


def transform_string_to_r_b(name: str, labels=("x", "y", "z")) -> Tuple:
    """Transform string to rotation part and offset.

    Example
    -------
    x,y,-z   ->   0.0 1 0 0  0.0 0 1 0  0.0 0 0 -1
    """
    l_name = "".join(name.strip().split()).lstrip("(").rstrip(")").split(",")
    rij, bi = [], []
    for _name in l_name:
        coefficients, offset = transform_string_to_digits(_name, labels)
        rij.append(coefficients)
        bi.append(offset)
    res_rij = numpy.array(rij, dtype=object)
    res_bi = numpy.array(bi, dtype=object)
    return res_rij, res_bi


def transform_string_to_digits(name: str, labels: Tuple[str]):
    """Transform string to digits.

    Multiplication has to be implicit, division must be explicit.
    White space within the string is optional.

    >>> transform_string_to_digits('-a+2x/7-0.7t+4', ('a', 'x', 't'))
    """
    coefficients = []
    offset = Fraction(0, 1).limit_denominator(10)
    l_name = name.strip().replace(" ", "").replace("+", " +").replace(
        "-", " -").split()
    for _label in labels:
        res = Fraction(0, 1).limit_denominator(10)
        for _name in l_name:
            if _label in _name:
                s_1 = _name.replace(_label, "").replace("+/", "+1/").replace(
                    "-/", "-1/")
                if s_1 == "":
                    s_1 = "1"
                if s_1.startswith("/"):
                    s_1 = "1" + s_1
                if s_1.endswith("+"):
                    s_1 = s_1 + "1"
                if s_1.endswith("-"):
                    s_1 = s_1 + "1"
                res += Fraction(s_1).limit_denominator(10)
        coefficients.append(res)
    res = Fraction(0, 1).limit_denominator(10)
    for _name in l_name:
        flag = all([not (_label in _name) for _label in labels])
        if flag:
            res += Fraction(_name).limit_denominator(10)
    offset = res
    return coefficients, offset


def transform_fraction_with_label_to_string(number: Fraction, label: str) \
        -> str:
    """Transform fraction with label to string.

    Parameters
    ----------
    number : Fraction
        DESCRIPTION.
    label : str
        DESCRIPTION.

    Returns
    -------
    str
        DESCRIPTION.

    """
    if isinstance(number, Fraction):
        val = number
    else:
        val = Fraction(number).limit_denominator(10)
    if val == Fraction(0, 1):
        res = ""
    elif val == Fraction(1, 1):
        res = f"{label:}"
    elif val == Fraction(-1, 1):
        res = f"-{label:}"
    elif val.denominator == 1:
        res = f"{val.numerator}{label:}"
    elif val.numerator == 1:
        res = f"{label:}/{val.denominator}"
    else:
        res = f"{val.numerator}{label:}/{val.denominator}"
    return res


def transform_digits_to_string(labels: Tuple[str], coefficients,
                               offset: Fraction) -> str:
    """Form a string from digits.

    Arguments
    ---------
        labels: the tuple of lablels (ex.: ('x', 'y', 'z') or ('a', 'b', 'c')))
        coefficients: the parameters in front of label (ex.: (1.0, 0.5, 0.0))
        offset: the number (ex.: 2/3)

    Output
    ------
    string

    Example
    -------
    >>> transform_digits_to_string(('x', 'y', 'z'), (1.0, 0.5, 0.0), 0.6666667)
    x+1/2y+2/3
    """
    l_res = []
    for _coefficient, _label in zip(coefficients, labels):
        _name = transform_fraction_with_label_to_string(_coefficient, _label)
        if _name == "":
            pass
        elif _name.startswith("-"):
            l_res.append(_name)
        elif l_res == []:
            l_res.append(_name)
        else:
            l_res.append(f"+{_name:}")
    _name = str(Fraction(offset).limit_denominator(10))
    if _name == "0":
        if l_res == []:
            l_res.append(_name)
    elif ((l_res == []) | (_name.startswith("-"))):
        l_res.append(_name)
    else:
        l_res.append(f"+{_name:}")
    return "".join(l_res)


def transform_r_b_to_string(r, b, labels=("x", "y", "z")) -> str:
    """Transform r (matrix) b (vector) to string.

    Parameters
    ----------
    r : TYPE
        DESCRIPTION.
    b : TYPE
        DESCRIPTION.
    labels : TYPE, optional
        DESCRIPTION. The default is ("x", "y", "z").

    Returns
    -------
    str
        DESCRIPTION.

    """
    l_res = [transform_digits_to_string(labels, _ri, _bi) for _ri, _bi in
             zip(r, b)]
    return ",".join(l_res)


def calc_mRmCmRT(r_ij, c_ij):
    """Calculate matrix multiplication R*C*RT.

    Matrices are expressed through its component and can be expressed as
    nD-array.

    r_ij: r11, r12, r13, r21, r22, r23, r31, r32, r33
    c_ij: c11, c12, c13, c21, c22, c23, c31, c32, c33
    """
    r11, r12, r13, r21, r22, r23, r31, r32, r33 = r_ij
    rt_ij = r11, r21, r31, r12, r22, r32, r13, r23, r33
    rcrt_ij = calc_product_matrices(r_ij, c_ij, rt_ij)
    rcrt_11, rcrt_12, rcrt_13, rcrt_21, rcrt_22, rcrt_23, rcrt_31, rcrt_32, \
        rcrt_33 = rcrt_ij
    return rcrt_11, rcrt_12, rcrt_13, rcrt_21, rcrt_22, rcrt_23, rcrt_31, \
        rcrt_32, rcrt_33


def calc_rotation_matrix_ij_by_euler_angles(alpha, beta, gamma):
    """Calculate rotational matrix from Euler angle.

    Tait-Bryan convection for angles is used.
    psi - x-axis
    theta - y axis
    phi - z axis
        0 <= psi   <= 2 pi ???
    -pi/2 <= theta <= pi/2
        0 <= phi   <= 2 pi
    """
    psi, theta, phi = alpha, beta, gamma
    rm_11 = numpy.cos(theta)*numpy.cos(phi)
    rm_21 = numpy.cos(theta)*numpy.sin(phi)
    rm_31 = -numpy.sin(theta)
    rm_12 = numpy.sin(psi)*numpy.sin(theta)*numpy.cos(phi) - \
        numpy.cos(psi)*numpy.sin(phi)
    rm_22 = numpy.sin(psi)*numpy.sin(theta)*numpy.sin(phi) + \
        numpy.cos(psi)*numpy.cos(phi)
    rm_32 = numpy.sin(psi)*numpy.cos(theta)
    rm_13 = numpy.cos(psi)*numpy.sin(theta)*numpy.cos(phi) + \
        numpy.sin(psi)*numpy.sin(phi)
    rm_23 = numpy.cos(psi)*numpy.sin(theta)*numpy.sin(phi) - \
        numpy.sin(psi)*numpy.cos(phi)
    rm_33 = numpy.cos(psi)*numpy.cos(theta)
    return rm_11, rm_12, rm_13, rm_21, rm_22, rm_23, rm_31, rm_32, rm_33


def calc_euler_angles_by_rotation_matrix_ij(rm_ij):
    """Calculate Euler Angles from rotational matrix.

    Tait-Bryan convection for angles used.
    psi - x-axis
    theta - y axis
    phi - z axis.
    """
    rm_11, rm_12, rm_13, rm_21, rm_22, rm_23, rm_31, rm_32, rm_33 = rm_ij
    if (abs(rm_31) != 1.):
        theta1 = -1*numpy.arcsin(rm_31)
        theta2 = numpy.pi-theta1
        psi1 = numpy.arctan2((rm_32)/numpy.cos(theta1),
                             (rm_33)/numpy.cos(theta1))
        psi2 = numpy.arctan2((rm_32)/numpy.cos(theta2),
                             (rm_33)/numpy.cos(theta2))
        phi1 = numpy.arctan2((rm_21)/numpy.cos(theta1),
                             (rm_11)/numpy.cos(theta1))
        phi2 = numpy.arctan2((rm_21)/numpy.cos(theta2),
                             (rm_11)/numpy.cos(theta2))
        res = [[psi1, theta1, phi1], [psi2, theta2, phi2]]
    else:
        phi = 0
        if (rm_31 == -1):
            theta = numpy.pi/2
            psi = phi+numpy.arctan2(rm_12, rm_13)
            res = [psi, theta, phi]
        else:
            theta = -1*numpy.pi/2
            psi = -1*phi+numpy.arctan2(-1*rm_12, -1*rm_13)
            res = [psi, theta, phi]
    return res


def calc_determinant_matrix_ij(m_ij):
    """Calculate determinant of the matrix / matrices.

    Parameters
    ----------
    m_ij : TYPE
        m_ij = m_11, m_12, m_13, m_21, m_22, m_23, m_31, m_32, m_33

    Returns
    -------
    det : TYPE
        DESCRIPTION.

    """
    m_11, m_12, m_13, m_21, m_22, m_23, m_31, m_32, m_33 = m_ij
    det = m_11*m_22*m_33 - m_11*m_23*m_32  \
        - m_12*m_21*m_33 + m_12*m_23*m_31  \
        + m_13*m_21*m_32 - m_13*m_22*m_31
    return det


def calc_inverse_matrix_ij(m_ij):
    """Calculate inverse matrix / matrices.

    Parameters
    ----------
    m_ij : TYPE
        m_ij = m_11, m_12, m_13, m_21, m_22, m_23, m_31, m_32, m_33

    Returns
    -------
    m_i_11 : TYPE
        DESCRIPTION.
    m_i_12 : TYPE
        DESCRIPTION.
    m_i_13 : TYPE
        DESCRIPTION.
    m_i_21 : TYPE
        DESCRIPTION.
    m_i_22 : TYPE
        DESCRIPTION.
    m_i_23 : TYPE
        DESCRIPTION.
    m_i_31 : TYPE
        DESCRIPTION.
    m_i_32 : TYPE
        DESCRIPTION.
    m_i_33 : TYPE
        DESCRIPTION.

    """
    eps = 1.0e-12
    det = calc_determinant_matrix_ij(m_ij)
    m_11, m_12, m_13, m_21, m_22, m_23, m_31, m_32, m_33 = m_ij
    inv_det = numpy.where(det < eps, 0., 1./det)
    m_i_11 = +(m_22*m_33-m_23*m_32) * inv_det
    m_i_21 = -(m_21*m_33-m_23*m_31) * inv_det
    m_i_31 = +(m_21*m_32-m_22*m_31) * inv_det
    m_i_12 = -(m_12*m_33-m_13*m_32) * inv_det
    m_i_22 = +(m_11*m_33-m_13*m_31) * inv_det
    m_i_32 = -(m_11*m_32-m_12*m_31) * inv_det
    m_i_13 = +(m_12*m_23-m_13*m_22) * inv_det
    m_i_23 = -(m_11*m_23-m_13*m_21) * inv_det
    m_i_33 = +(m_11*m_22-m_12*m_21) * inv_det
    return m_i_11, m_i_12, m_i_13, m_i_21, m_i_22, m_i_23, m_i_31, m_i_32, \
        m_i_33


def calc_rotation_matrix_ij_around_axis(angle, axis="x"):
    """Calculate rotation matrix around given axis.

    Parameters
    ----------
    angle : TYPE
        DESCRIPTION.
    axis : TYPE, optional
        DESCRIPTION. The default is "x".

    Returns
    -------
    r_m_11 : TYPE
        DESCRIPTION.
    r_m_12 : TYPE
        DESCRIPTION.
    r_m_13 : TYPE
        DESCRIPTION.
    r_m_21 : TYPE
        DESCRIPTION.
    r_m_22 : TYPE
        DESCRIPTION.
    r_m_23 : TYPE
        DESCRIPTION.
    r_m_31 : TYPE
        DESCRIPTION.
    r_m_32 : TYPE
        DESCRIPTION.
    r_m_33 : TYPE
        DESCRIPTION.

    """
    np_zero = 0.*angle
    np_one = 1. + np_zero
    if axis.lower() == "y":
        r_m_11 = numpy.cos(angle)
        r_m_12 = np_zero
        r_m_13 = numpy.sin(angle)
        r_m_21 = np_zero
        r_m_22 = np_one
        r_m_23 = np_zero
        r_m_31 = -numpy.sin(angle)
        r_m_32 = np_zero
        r_m_33 = numpy.cos(angle)
    elif axis.lower() == "z":
        r_m_11 = numpy.cos(angle)
        r_m_12 = numpy.sin(angle)
        r_m_13 = np_zero
        r_m_21 = -numpy.sin(angle)
        r_m_22 = numpy.cos(angle)
        r_m_23 = np_zero
        r_m_31 = np_zero
        r_m_32 = np_zero
        r_m_33 = np_one
    else:
        # by default "x"
        r_m_11 = np_one
        r_m_12 = np_zero
        r_m_13 = np_zero
        r_m_21 = np_zero
        r_m_22 = numpy.cos(angle)
        r_m_23 = numpy.sin(angle)
        r_m_31 = np_zero
        r_m_32 = -numpy.sin(angle)
        r_m_33 = numpy.cos(angle)
    return r_m_11, r_m_12, r_m_13, r_m_21, r_m_22, r_m_23, r_m_31, r_m_32, \
        r_m_33


def calc_product_matrices(m_1_ij, m_2_ij, *argv):
    """Product of two or more matrices."""
    if len(argv) == 0:
        m_1_11, m_1_12, m_1_13, m_1_21, m_1_22, m_1_23, m_1_31, m_1_32, \
            m_1_33 = m_1_ij
        m_2_11, m_2_12, m_2_13, m_2_21, m_2_22, m_2_23, m_2_31, m_2_32, \
            m_2_33 = m_2_ij

        m_3_11 = m_1_11*m_2_11 + m_1_12*m_2_21 + m_1_13*m_2_31
        m_3_12 = m_1_11*m_2_12 + m_1_12*m_2_22 + m_1_13*m_2_32
        m_3_13 = m_1_11*m_2_13 + m_1_12*m_2_23 + m_1_13*m_2_33

        m_3_21 = m_1_21*m_2_11 + m_1_22*m_2_21 + m_1_23*m_2_31
        m_3_22 = m_1_21*m_2_12 + m_1_22*m_2_22 + m_1_23*m_2_32
        m_3_23 = m_1_21*m_2_13 + m_1_22*m_2_23 + m_1_23*m_2_33

        m_3_31 = m_1_31*m_2_11 + m_1_32*m_2_21 + m_1_33*m_2_31
        m_3_32 = m_1_31*m_2_12 + m_1_32*m_2_22 + m_1_33*m_2_32
        m_3_33 = m_1_31*m_2_13 + m_1_32*m_2_23 + m_1_33*m_2_33
    else:
        m_h_ij = calc_product_matrices(m_1_ij, m_2_ij)
        m_3_11, m_3_12, m_3_13, m_3_21, m_3_22, m_3_23, m_3_31, m_3_32, \
            m_3_33 = calc_product_matrices(m_h_ij, *argv)
    return m_3_11, m_3_12, m_3_13, m_3_21, m_3_22, m_3_23, m_3_31, m_3_32, \
        m_3_33


def calc_product_matrix_vector(m_ij, v_j):
    """Product matrix and vector."""
    m_11, m_12, m_13, m_21, m_22, m_23, m_31, m_32, m_33 = m_ij
    v_1, v_2, v_3 = v_j
    o_1 = m_11*v_1 + m_12*v_2 + m_13*v_3
    o_2 = m_21*v_1 + m_22*v_2 + m_23*v_3
    o_3 = m_31*v_1 + m_32*v_2 + m_33*v_3
    return o_1, o_2, o_3


def calc_vector_angle(v_1_i, v_2_i):
    """Angle between two vectors (in radians)."""
    eps = 1e-12
    v_1_1, v_1_2, v_1_3 = v_1_i
    v_2_1, v_2_2, v_2_3 = v_2_i
    prod_v_12 = v_1_1*v_2_1 + v_1_2*v_2_2 + v_1_3*v_2_3
    norm_v_1 = numpy.sqrt(numpy.square(v_1_1) + numpy.square(v_1_2) +
                          numpy.square(v_1_3))
    norm_v_2 = numpy.sqrt(numpy.square(v_2_1) + numpy.square(v_2_2) +
                          numpy.square(v_2_3))
    flag_12 = numpy.logical_or((norm_v_1 > eps), (norm_v_2 > eps))

    angle = numpy.arccos(numpy.where(flag_12, prod_v_12 / (norm_v_1*norm_v_2),
                                     0.99999))
    return angle


def calc_vector_product(v_1_i, v_2_i):
    """Vector product of vectors."""
    v_1_1, v_1_2, v_1_3 = v_1_i
    v_2_1, v_2_2, v_2_3 = v_2_i
    o_1 = v_1_2 * v_2_3 - v_1_3 * v_2_2
    o_2 = v_1_3 * v_2_1 - v_1_1 * v_2_3
    o_3 = v_1_1 * v_2_2 - v_1_2 * v_2_1
    return o_1, o_2, o_3


def scalar_product(v_1_i, v_2_i):
    """Scalar product of vectors."""
    v_1_1, v_1_2, v_1_3 = v_1_i
    v_2_1, v_2_2, v_2_3 = v_2_i
    return v_1_1*v_2_1 + v_1_2*v_2_2 + v_1_3*v_2_3


def calc_rotation_matrix_by_two_vectors(v_1_i, v_2_i, ax_1="x", ax_2="y"):
    """Calculate rotation matrix by two vectors."""
    norm_v_1 = numpy.sqrt(numpy.square(v_1_i[0]) + numpy.square(v_1_i[1]) +
                          numpy.square(v_1_i[2]))
    norm_v_2 = numpy.sqrt(numpy.square(v_2_i[0]) + numpy.square(v_2_i[1]) +
                          numpy.square(v_2_i[2]))
    e_x = (v_1_i[0]/norm_v_1, v_1_i[1]/norm_v_1, v_1_i[2]/norm_v_1)
    e_12 = (v_2_i[0]/norm_v_2, v_2_i[1]/norm_v_2, v_2_i[2]/norm_v_2)
    v_z = calc_vector_product(e_x, e_12)
    norm_v_z = numpy.sqrt(numpy.square(v_z[0]) + numpy.square(v_z[1]) +
                          numpy.square(v_z[2]))
    e_z = (v_z[0]/norm_v_z, v_z[1]/norm_v_z, v_z[2]/norm_v_z)
    e_y = calc_vector_product(e_z, e_x)
    rm_11, rm_12, rm_13 = e_x[0], e_y[0], e_z[0]
    rm_21, rm_22, rm_23 = e_x[1], e_y[1], e_z[1]
    rm_31, rm_32, rm_33 = e_x[2], e_y[2], e_z[2]
    return rm_11, rm_12, rm_13, rm_21, rm_22, rm_23, rm_31, rm_32, rm_33


def tri_linear_interpolation(ind_xyz, den_3d):
    """Linear interpolation by four points in 3d array.

    Parameters
    ----------
    ind_xyz : [:, 3] of float
        DESCRIPTION.
    den_3d : TYPE
        DESCRIPTION.

    Returns
    -------
    c : TYPE
        DESCRIPTION.

    """
    ind_x, ind_y, ind_z = ind_xyz[:, 0], ind_xyz[:, 1], ind_xyz[:, 2]
    (n_x, n_y, n_z) = den_3d.shape
    ind_x0 = (numpy.floor(ind_x)).astype(int)
    ind_x1 = ind_x0+1
    ind_y0 = (numpy.floor(ind_y)).astype(int)
    ind_y1 = ind_y0+1
    ind_z0 = (numpy.floor(ind_z)).astype(int)
    ind_z1 = ind_z0+1

    xd = (ind_x-ind_x0)/(ind_x1-ind_x0)
    yd = (ind_y-ind_y0)/(ind_y1-ind_y0)
    zd = (ind_z-ind_z0)/(ind_z1-ind_z0)

    Vx0y0z0 = den_3d[numpy.mod(ind_x0, n_x), numpy.mod(ind_y0, n_y),
                     numpy.mod(ind_z0, n_z)]
    Vx1y0z0 = den_3d[numpy.mod(ind_x1, n_x), numpy.mod(ind_y0, n_y),
                     numpy.mod(ind_z0, n_z)]
    Vx0y0z1 = den_3d[numpy.mod(ind_x0, n_x), numpy.mod(ind_y0, n_y),
                     numpy.mod(ind_z1, n_z)]
    Vx1y0z1 = den_3d[numpy.mod(ind_x1, n_x), numpy.mod(ind_y0, n_y),
                     numpy.mod(ind_z1, n_z)]
    Vx0y1z0 = den_3d[numpy.mod(ind_x0, n_x), numpy.mod(ind_y1, n_y),
                     numpy.mod(ind_z0, n_z)]
    Vx1y1z0 = den_3d[numpy.mod(ind_x1, n_x), numpy.mod(ind_y1, n_y),
                     numpy.mod(ind_z0, n_z)]
    Vx0y1z1 = den_3d[numpy.mod(ind_x0, n_x), numpy.mod(ind_y1, n_y),
                     numpy.mod(ind_z1, n_z)]
    Vx1y1z1 = den_3d[numpy.mod(ind_x1, n_x), numpy.mod(ind_y1, n_y),
                     numpy.mod(ind_z1, n_z)]

    c00 = Vx0y0z0*(1.-xd) + Vx1y0z0*xd
    c01 = Vx0y0z1*(1.-xd) + Vx1y0z1*xd
    c10 = Vx0y1z0*(1.-xd) + Vx1y1z0*xd
    c11 = Vx0y1z1*(1.-xd) + Vx1y1z1*xd

    c0 = c00*(1.-yd) + c10*yd
    c1 = c01*(1.-yd) + c11*yd

    c = c0*(1.-zd) + c1*zd

    return c


def ortogonalize_matrix(m_ij, m_norm_ij):
    """Ortogonalize matrix.

    matrix m_ij is defined in coordinate system (a, b, c).
    It is given as tuple

    Matrix m_norm  used to recalculate coordinates  from direct space
    (a1/|a1|, a2/|a2|, a3/|a3|) to Cartesian one (x||a*, z||c).

    x_cart = m_norm * x_direct

    m_norm  = [[(1 - cos**2 alpha1 - cos**2 alpha2 - cos**2 alpha3 +
             2 cos alpha1 cos alpha2 cos alpha3)**0.5/sin(alpha1),  0,  0],
     [  (cos alpha3 - cos alpha1 cos alpha2) / sin alpha1, sin alpha1,  0],
     [                                         cos alpha2, cos alpha1,  1]]

    Should be given as

    m_norm_ij = (_11, _12, _13, _21, _22, _23, _31, _32, _33)

    Output matrix s_ij is defined in Cartezian coordinate system where
    x||a*, z||c, y= [z x] (right handed)
    """
    s_11, s_12, s_13, s_21, s_22, s_23, s_31, s_32, s_33 = calc_mRmCmRT(
        m_norm_ij, m_ij)
    return s_11, s_12, s_13, s_21, s_22, s_23, s_31, s_32, s_33


def calc_moment_2d_by_susceptibility(r_ij, susc_i, m_norm_ij, h_loc):
    """Calculate  moment 2d by susceptibility.

    Recalculate chi_i given in reciprocal unit cell according to symmetry
    elements for each point.

    After susceptibility is multiplied in magnetic field defined

    r_ij:= r_11, r_12, r_13, r_21, r_22, r_23, r_31, r_32, r_33
    susc_i:= chi_11, chi_22, chi_33, chi_12, chi_13, chi_23

    Matrix m_norm  used to recalculate coordinates  from direct space
    (a1/|a1|, a2/|a2|, a3/|a3|) to Cartesian one (x||a*, z||c).

    x_cart = m_norm * x_direct

    m_norm  = [[(1 - cos**2 alpha1 - cos**2 alpha2 - cos**2 alpha3 +
                 2 cos alpha1 cos alpha2 cos alpha3)**0.5/sin(alpha1),  0,  0],
         [  (cos alpha3 - cos alpha1 cos alpha2) / sin alpha1, sin alpha1,  0],
         [                                         cos alpha2, cos alpha1,  1]]

    Matrix m_norm should be given as

        m_norm_ij = (_11, _12, _13, _21, _22, _23, _31, _32, _33)

    moment_2d: [points, symmetry]
    """
    r_11, r_12, r_13, r_21, r_22, r_23, r_31, r_32, r_33 = r_ij
    chi_11, chi_22, chi_33, chi_12, chi_13, chi_23 = susc_i
    # [ind, symm]
    chi_2d_ij = calc_mRmCmRT(
        (r_11[numpy.newaxis, :], r_12[numpy.newaxis, :],
         r_13[numpy.newaxis, :], r_21[numpy.newaxis, :],
         r_22[numpy.newaxis, :], r_23[numpy.newaxis, :],
         r_31[numpy.newaxis, :], r_32[numpy.newaxis, :],
         r_33[numpy.newaxis, :]),
        (chi_11[:, numpy.newaxis], chi_12[:, numpy.newaxis],
         chi_13[:, numpy.newaxis], chi_12[:, numpy.newaxis],
         chi_22[:, numpy.newaxis], chi_23[:, numpy.newaxis],
         chi_13[:, numpy.newaxis], chi_23[:, numpy.newaxis],
         chi_33[:, numpy.newaxis]))
    chi_orto_ij = ortogonalize_matrix(chi_2d_ij, m_norm_ij)
    moment_2d = calc_product_matrix_vector(chi_orto_ij, h_loc)
    return moment_2d


def calc_phase_3d(hkl, r_ij, b_i, fract_xyz):
    r"""Calculate phases over 3 dimensions (hkl, points, symmetry).

    fract_xyz.shape = (pointxs, 3)

    r_ij, b_i are elements of symmetry

    phase = exp(2\pi i (h*Rs*x + h*bs))

    phase_3d: [hkl, points, symmetry]
    """
    index_h, index_k, index_l = hkl[0], hkl[1], hkl[2]
    x, y, z = fract_xyz[:, 0], fract_xyz[:, 1], fract_xyz[:, 2]
    r_11, r_12, r_13, r_21, r_22, r_23, r_31, r_32, r_33 = r_ij
    b_1, b_2, b_3 = b_i

    # [ind, symm]
    mr_2d_1, mr_2d_2, mr_2d_3 = calc_product_matrix_vector(
        (r_11[numpy.newaxis, :], r_12[numpy.newaxis, :],
         r_13[numpy.newaxis, :], r_21[numpy.newaxis, :],
         r_22[numpy.newaxis, :], r_23[numpy.newaxis, :],
         r_31[numpy.newaxis, :], r_32[numpy.newaxis, :],
         r_33[numpy.newaxis, :]),
        (x[:, numpy.newaxis], y[:, numpy.newaxis], z[:, numpy.newaxis]))

    # [hkl, ind, symm]
    hrr_3d = scalar_product(
        (index_h[:, numpy.newaxis, numpy.newaxis],
         index_k[:, numpy.newaxis, numpy.newaxis],
         index_l[:, numpy.newaxis, numpy.newaxis]),
        (mr_2d_1[numpy.newaxis, :, :], mr_2d_2[numpy.newaxis, :, :],
         mr_2d_3[numpy.newaxis, :, :]))

    # [hkl, symm]
    hb_2d = scalar_product(
        (index_h[:, numpy.newaxis], index_k[:, numpy.newaxis],
         index_l[:, numpy.newaxis]),
        (b_1[numpy.newaxis, :], b_2[numpy.newaxis, :], b_3[numpy.newaxis, :]))

    phase_3d = numpy.exp(2.*numpy.pi*1j*(hrr_3d+hb_2d[:, numpy.newaxis, :]))
    return phase_3d
