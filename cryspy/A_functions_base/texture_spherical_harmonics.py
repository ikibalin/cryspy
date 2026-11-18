"""Texture of a powder by the generalized spherical harmonics (GSAS-II model).

The even part of the orientation distribution function gives the texture
correction of the integrated intensity of the reflection h measured along
the sample direction y (Bunge, 1982; Von Dreele, 1997):

    A(h, y) = 1 + sum_l 4 pi / (2l + 1) sum_m sum_n C_l^mn k_l^m(y) k_l^n(h),

l = 2, 4, ..., order_max. k_l^n(h) are the crystal harmonics symmetrized by
the Laue class, k_l^m(y) the sample harmonics symmetrized by the sample
symmetry (cylindrical, 2/m, mmm or none). The harmonics, their selection
rules, the orientation of the crystal and sample frames and the sample
angles follow GSAS-II (GSASIIlattice.py: OdfChk, GenSHCoeffT, CrsAng,
SamAng, GetKcl, GetKsl; GSASIIstrMath.py: SHTXcal), so the coefficients
C_l^mn refined in CrysPy and GSAS-II are the same numbers. The harmonics
are evaluated with the real spherical harmonics shared with the multipole
model (spherical_harmonics.py).

Frames
------
Crystal: z || c* and x || a* (x || a for the hexagonal axes), except the
monoclinic unique axis b (z || b*, x || c*), unique axis a (z || a*,
x || b*) and the rhombohedral axes (z || [111]*, x || a* - b*); x is made
normal to z. The frames of GSAS-II are the same, but GSAS-II takes the
azimuth of the reflections of oblique (monoclinic and triclinic) cells by
an approximate formula, so for them the coefficients agree with GSAS-II
only when the oblique angle is close to 90 degrees. For the Laue class -31m
GSAS-II uses the harmonics of -3m1 (not invariant under -31m for n = 3, 9,
...); here every class has harmonics invariant under its rotations.
The polar angle and the azimuth of the reflection are taken in this frame.
Sample: the sample direction of the scattering vector (psi, gamma) is
calculated from the Bragg angle, the goniometer angles omega, chi, phi, the
azimuth of the detector and the zeros of the sample orientation (GSAS-II
SamAng). For the cylindrical symmetry only psi, the angle between the
scattering vector and the cylinder axis, is used.

References
----------
H.-J. Bunge, Texture Analysis in Materials Science, Butterworths (1982);
Table 15.2.1 of the edition by Cuvillier Verlag (1993) for the cubic
harmonics.
R. B. Von Dreele, J. Appl. Cryst. 30, 517-525 (1997).
B. H. Toby, R. B. Von Dreele, J. Appl. Cryst. 46, 544-549 (2013) (GSAS-II).
"""
import numpy

from cryspy.A_functions_base.spherical_harmonics import (
    calc_real_spherical_harmonics_by_lm,
    cartesian_to_spherical,
)
from cryspy.A_functions_base.unit_cell import calc_eq_ccs_by_unit_cell_parameters

na = numpy.newaxis

SAMPLE_SYMMETRIES = ("cylindrical", "2/m", "mmm", "-1")

# Laue classes in the frames of GSAS-II. "-3R", "-3mR": rhombohedral axes;
# "-3m1" (2-fold axes along a, mirrors normal to x || a) and "-31m" (2-fold
# axes along a - b, mirrors containing x) differ by the orientation of the
# 2-fold axes in the hexagonal setting.
LAUE_CLASSES = ("-1", "2/m", "mmm", "4/m", "4/mmm", "-3", "-3R", "-3m1",
                "-31m", "-3mR", "6/m", "6/mmm", "m-3", "m-3m")

# Laue class candidates by the number of the rotations (with -1) of the class.
LAUE_CANDIDATES = {
    2: (("-1", ""),),
    4: (("2/m", "b"), ("2/m", "c"), ("2/m", "a")),
    8: (("mmm", ""), ("4/m", "")),
    16: (("4/mmm", ""),),
    6: (("-3", ""), ("-3R", "")),
    12: (("-3m1", ""), ("-31m", ""), ("-3mR", ""), ("6/m", "")),
    24: (("6/mmm", ""), ("m-3", "")),
    48: (("m-3m", ""),),
}

# Even order terms of Table 15.2.1 of Bunge (1993), p. 501: cubic
# harmonics with [001] parallel to z, coefficients of the normalized
# associated Legendre functions P_l^j cos(j beta), j = 0, 4, 8, ... (one row
# for each independent harmonic n = 1, 2, ...). Copied from GSAS-II (BOH).
CUBIC_COEFFICIENTS = {
    4: [[0.30469720, 0.36418281]],
    6: [[-0.14104740, 0.52775103]],
    8: [[0.28646862, 0.21545346, 0.32826995]],
    10: [[-0.16413497, 0.33078546, 0.39371345]],
    12: [[0.26141975, 0.27266871, 0.03277460, 0.32589402],
         [0.09298802, -0.23773812, 0.49446631, 0.0]],
    14: [[-0.17557309, 0.25821932, 0.27709173, 0.33645360]],
    16: [[0.24370673, 0.29873515, 0.06447688, 0.00377, 0.32574495],
         [0.12039646, -0.25330128, 0.23950998, 0.40962508, 0.0]],
    18: [[-0.16914245, 0.17017340, 0.34598142, 0.07433932, 0.32696037],
         [-0.06901768, 0.16006562, -0.24743528, 0.47110273, 0.0]],
    20: [[0.23067026, 0.31151832, 0.09287682, 0.01089683, 0.00037564, 0.32573563],
         [0.13615420, -0.25048007, 0.12882081, 0.28642879, 0.34620433, 0.0]],
    22: [[-0.16109560, 0.10244188, 0.36285175, 0.13377513, 0.01314399, 0.32585583],
         [-0.09620055, 0.20244115, -0.22389483, 0.17928946, 0.42017231, 0.0]],
    24: [[0.22050742, 0.31770654, 0.11661736, 0.02049853, 0.00150861, 0.00003426,
          0.32573505],
         [0.13651722, -0.21386648, 0.00522051, 0.33939435, 0.10837396, 0.32914497, 0.0],
         [0.05378596, -0.11945819, 0.16272298, -0.26449730, 0.44923956, 0.0, 0.0]],
    26: [[-0.15435003, 0.05261630, 0.35524646, 0.18578869, 0.03259103, 0.00186197,
          0.32574594],
         [-0.11306511, 0.22072681, -0.18706142, 0.05439948, 0.28122966, 0.35634355, 0.0]],
    28: [[0.21225019, 0.32031716, 0.13604702, 0.03132468, 0.00362703, 0.00018294,
          0.00000294, 0.32573501],
         [0.13219496, -0.17206256, -0.08742608, 0.32671661, 0.17973107, 0.02567515,
          0.32619598, 0.0],
         [0.07989184, -0.16735346, 0.18839770, -0.20705337, 0.12926808, 0.42715602,
          0.0, 0.0]],
    30: [[-0.14878368, 0.01524973, 0.33628434, 0.22632587, 0.05790047, 0.00609812,
          0.00022898, 0.32573594],
         [-0.11721726, 0.20915005, -0.11723436, -0.07815329, 0.31318947, 0.13655742,
          0.33241385, 0.0],
         [-0.04297703, 0.09317876, -0.11831248, 0.17355132, -0.28164031, 0.42719361,
          0.0, 0.0]],
    32: [[0.20533892, 0.32087437, 0.15187897, 0.04249238, 0.00670516, 0.00054977,
          0.00002018, 0.00000024, 0.32573501],
         [0.12775091, -0.13523423, -0.14935701, 0.28227378, 0.23670434, 0.05661270,
          0.00469819, 0.32578978, 0.0],
         [0.09703829, -0.19373733, 0.18610682, -0.14407046, 0.00220535, 0.26897090,
          0.36633402, 0.0, 0.0]],
    34: [[-0.14409234, -0.01343681, 0.31248977, 0.25557722, 0.08571889, 0.01351208,
          0.00095792, 0.00002550, 0.32573508],
         [-0.11527834, 0.18472133, -0.04403280, -0.16908618, 0.27227021, 0.21086614,
          0.04041752, 0.32688152, 0.0],
         [-0.06773139, 0.14120811, -0.15835721, 0.18357456, -0.19364673, 0.08377174,
          0.43116318, 0.0, 0.0]],
}
ORDER_MAX_CUBIC = max(CUBIC_COEFFICIENTS)


def is_allowed_harmonic(symmetry: str, l: int, m: int) -> bool:
    """Selection rule of the symmetrized harmonics k_l^m (GSAS-II OdfChk).

    Parameters
    ----------
    symmetry : str
        Laue class (one of LAUE_CLASSES) or sample symmetry (one of
        SAMPLE_SYMMETRIES).
    l : int
        Order; only even orders are used (the even part of the ODF).
    m : int
        Index, -l <= m <= l (for the cubic classes the number n >= 1 of the
        independent cubic harmonic of the order l).

    Returns
    -------
    bool
        True when the harmonic is allowed by the symmetry.
    """
    if l % 2 != 0 or abs(m) > l:
        return False
    if symmetry == "cylindrical":
        return m == 0
    if symmetry == "-1":
        return True
    if symmetry == "2/m":
        return abs(m) % 2 == 0
    if symmetry == "mmm":
        return abs(m) % 2 == 0 and m >= 0
    if symmetry == "4/m":
        return abs(m) % 4 == 0
    if symmetry == "4/mmm":
        return abs(m) % 4 == 0 and m >= 0
    if symmetry in ("-3", "-3R"):
        return abs(m) % 3 == 0
    if symmetry in ("-3m1", "-31m", "-3mR"):
        return abs(m) % 3 == 0 and m >= 0
    if symmetry == "6/m":
        return abs(m) % 6 == 0
    if symmetry == "6/mmm":
        return abs(m) % 6 == 0 and m >= 0
    if symmetry in ("m-3", "m-3m"):
        if m <= 0:
            return False
        if l % 12 == 2:
            return m <= l // 12
        return m <= l // 12 + 1
    raise ValueError(f"Unknown symmetry '{symmetry}' of the spherical harmonics.")


def form_texture_sh_lmn(laue_class: str, sample_symmetry: str, order_max: int):
    """Indexes (l, m, n) of the allowed coefficients C_l^mn (GSAS-II GenSHCoeffT).

    Parameters
    ----------
    laue_class : str
        Laue class of the crystal (one of LAUE_CLASSES).
    sample_symmetry : str
        Sample symmetry (one of SAMPLE_SYMMETRIES).
    order_max : int
        Maximal even order l.

    Returns
    -------
    numpy.ndarray
        Integer array (3, n_coefficient): l, m (sample), n (crystal).
    """
    if laue_class in ("m-3", "m-3m") and order_max > ORDER_MAX_CUBIC:
        raise ValueError(f"The cubic harmonics are tabulated up to the order {ORDER_MAX_CUBIC}.")
    l_lmn = []
    for l in range(2, int(order_max) + 1, 2):
        for m in range(-l, l + 1):
            if not is_allowed_harmonic(sample_symmetry, l, m):
                continue
            for n in range(-l, l + 1):
                if is_allowed_harmonic(laue_class, l, n):
                    l_lmn.append((l, m, n))
    return numpy.array(l_lmn, dtype=int).reshape(-1, 3).transpose()


def calc_crystal_frame(unit_cell_parameters, laue_class: str, unique_axis: str = ""):
    """Unit vectors x, y, z of the crystal frame of the harmonics.

    Parameters
    ----------
    unit_cell_parameters : array_like
        a, b, c (Å), alpha, beta, gamma (radians).
    laue_class : str
        Laue class (one of LAUE_CLASSES).
    unique_axis : str
        Unique axis "a", "b" or "c" of the monoclinic class 2/m.

    Returns
    -------
    numpy.ndarray
        Array (3, 3): rows are x, y, z in the Cartesian system of CrysPy
        (x || a*, z || c).
    """
    eq_reciprocal = calc_eq_ccs_by_unit_cell_parameters(
        numpy.array([[1, 0, 0, 1], [0, 1, 0, 1], [0, 0, 1, 1]], dtype=float),
        unit_cell_parameters)[0]
    a_star, b_star, c_star, star_111 = eq_reciprocal.transpose()
    if laue_class == "2/m" and unique_axis == "b":
        z_axis, x_axis = b_star, c_star
    elif laue_class == "2/m" and unique_axis == "a":
        z_axis, x_axis = a_star, b_star
    elif laue_class in ("-3R", "-3mR"):
        z_axis, x_axis = star_111, a_star - b_star
    elif laue_class in ("-3", "-3m1", "-31m", "6/m", "6/mmm"):
        z_axis, x_axis = c_star, numpy.cross(b_star, c_star)
    else:
        z_axis, x_axis = c_star, a_star
    z_axis = z_axis / numpy.linalg.norm(z_axis)
    x_axis = x_axis - numpy.dot(x_axis, z_axis) * z_axis
    x_axis = x_axis / numpy.linalg.norm(x_axis)
    y_axis = numpy.cross(z_axis, x_axis)
    return numpy.stack([x_axis, y_axis, z_axis], axis=0)


def calc_crystal_angles(index_hkl, unit_cell_parameters, laue_class: str,
                        unique_axis: str = ""):
    """Polar angle and azimuth of the reflections in the crystal frame.

    The analogue of GSAS-II CrsAng: the polar angle from c* (b*, a* or
    [111]*) and the azimuth from x of calc_crystal_frame.

    Parameters
    ----------
    index_hkl : array_like
        Miller indexes, shape (3, n_hkl).
    unit_cell_parameters : array_like
        a, b, c (Å), alpha, beta, gamma (radians).
    laue_class : str
        Laue class (one of LAUE_CLASSES).
    unique_axis : str
        Unique axis of the monoclinic class.

    Returns
    -------
    polar, azimuth : numpy.ndarray
        Angles in radians, shape (n_hkl,).
    """
    index_hkl = numpy.asarray(index_hkl, dtype=float)
    eq_hkl = calc_eq_ccs_by_unit_cell_parameters(index_hkl, unit_cell_parameters)[0]
    frame = calc_crystal_frame(unit_cell_parameters, laue_class, unique_axis)
    xyz = numpy.matmul(frame, eq_hkl)
    polar, azimuth = cartesian_to_spherical(xyz[0], xyz[1], xyz[2])[1:]
    return polar, azimuth


def calc_crystal_harmonics(index_ln, laue_class: str, polar, azimuth):
    """Symmetrized crystal harmonics k_l^n(h) (GSAS-II GetKcl).

    Parameters
    ----------
    index_ln : array_like of int
        Orders and indexes, shape (2, n_harmonic).
    laue_class : str
        Laue class (one of LAUE_CLASSES).
    polar, azimuth : array_like
        Angles of the reflections in the crystal frame, radians.

    Returns
    -------
    numpy.ndarray
        k_l^n of the shape (n_harmonic,) + polar.shape.
    """
    index_ln = numpy.asarray(index_ln, dtype=int).reshape(2, -1)
    polar = numpy.asarray(polar, dtype=float)
    # every harmonic is a combination of real spherical harmonics Y_l^j:
    # (l, j, weight) for each harmonic, evaluated together in one call
    l_terms = []
    for l, n in index_ln.transpose():
        l, n = int(l), int(n)
        if laue_class in ("m-3", "m-3m"):
            # Bunge's coefficients multiply P_l^j cos(j beta) without the
            # normalization of the azimuthal part
            l_terms.append([(l, 4 * i_j, coefficient * numpy.sqrt(
                (2. if i_j == 0 else 1.) * numpy.pi))
                for i_j, coefficient in enumerate(CUBIC_COEFFICIENTS[l][n - 1])])
        elif n == 0 or laue_class in ("mmm", "4/mmm", "6/mmm", "-31m"):
            # -31m: the mirrors contain x || a (beta -> -beta), so only
            # cos(n beta) is invariant; GSAS-II uses the rule of -3m1 for
            # -31m too, which is not invariant for n = 3, 9, ...
            l_terms.append([(l, abs(n), 1.)])
        elif laue_class in ("-3m1", "-3mR"):
            # mirrors normal to x (beta -> pi - beta): sin(n beta) for
            # n = 3, 9, ...
            l_terms.append([(l, -n, 1.)] if n % 6 == 3 else [(l, n, 1.)])
        else:
            # cos(n beta) + sin(n beta) of GSAS-II for both signs of n
            l_terms.append([(l, abs(n), 1.), (l, -abs(n), float(numpy.sign(n)))])
    return sum_real_spherical_harmonics(l_terms, polar, azimuth)


def sum_real_spherical_harmonics(l_terms, theta, phi):
    """Sums of weighted real spherical harmonics.

    Parameters
    ----------
    l_terms : list of list of (int, int, float)
        For every sum the terms (l, m, weight) of weight * Y_l^m.
    theta, phi : array_like
        Polar angle and azimuth, radians (broadcast together).

    Returns
    -------
    numpy.ndarray
        Sums of the shape (len(l_terms),) + broadcast shape of theta, phi;
        every distinct Y_l^m is evaluated once.
    """
    l_lm = sorted({(l, m) for terms in l_terms for l, m, weight in terms})
    d_index = {lm: i_lm for i_lm, lm in enumerate(l_lm)}
    theta, phi = numpy.broadcast_arrays(numpy.asarray(theta, dtype=float),
                                        numpy.asarray(phi, dtype=float))
    res = numpy.zeros((len(l_terms),) + theta.shape, dtype=float)
    if len(l_lm) == 0:
        return res
    y_lm = calc_real_spherical_harmonics_by_lm(numpy.array(l_lm, dtype=int).transpose(),
                                               theta, phi)
    for i_term, terms in enumerate(l_terms):
        for l, m, weight in terms:
            res[i_term] += weight * y_lm[d_index[(l, m)]]
    return res


def calc_sample_harmonics(index_lm, sample_symmetry: str, psi, gamma):
    """Symmetrized sample harmonics k_l^m(y) (GSAS-II GetKsl).

    Parameters
    ----------
    index_lm : array_like of int
        Orders and indexes, shape (2, n_harmonic).
    sample_symmetry : str
        Sample symmetry (one of SAMPLE_SYMMETRIES).
    psi, gamma : array_like
        Polar angle and azimuth of the scattering vector in the sample
        frame, radians.

    Returns
    -------
    numpy.ndarray
        k_l^m of the shape (n_harmonic,) + broadcast shape of psi, gamma.
    """
    index_lm = numpy.asarray(index_lm, dtype=int).reshape(2, -1)
    l_terms = []
    for l, m in index_lm.transpose():
        l, m = int(l), int(m)
        if m == 0 or sample_symmetry == "mmm":
            l_terms.append([(l, abs(m), 1.)])
        else:
            l_terms.append([(l, abs(m), 1.), (l, -abs(m), float(numpy.sign(m)))])
    return sum_real_spherical_harmonics(l_terms, psi, gamma)


def calc_sample_angles(theta, goniometer_omega, goniometer_chi, goniometer_phi,
                       azimuth, sample_omega=0., sample_chi=0., sample_phi=0.,
                       flag_coupled: bool = False):
    """Direction of the scattering vector in the sample frame (GSAS-II SamAng).

    Parameters
    ----------
    theta : array_like
        Bragg angle (half of 2theta), radians.
    goniometer_omega, goniometer_chi, goniometer_phi : float
        Goniometer angles of the sample, radians.
    azimuth : array_like
        Azimuth of the detector around the beam, radians (0 for the
        scattering in the horizontal plane).
    sample_omega, sample_chi, sample_phi : float
        Zeros of the sample orientation, radians.
    flag_coupled : bool
        True when omega follows theta (Bragg-Brentano geometry).

    Returns
    -------
    psi, gamma : numpy.ndarray
        Polar angle and azimuth of the scattering vector in the sample
        frame, radians (broadcast shape of theta and azimuth).
    """
    theta = numpy.asarray(theta, dtype=float)
    azimuth = numpy.asarray(azimuth, dtype=float)
    omega = goniometer_omega + theta if flag_coupled else goniometer_omega + 0. * theta
    s_omega, c_omega = numpy.sin(omega), numpy.cos(omega)
    s_theta, c_theta = numpy.sin(theta), numpy.cos(theta)
    s_azimuth, c_azimuth = numpy.sin(azimuth), numpy.cos(azimuth)
    s_chi, c_chi = numpy.sin(goniometer_chi), numpy.cos(goniometer_chi)
    s_phi = numpy.sin(goniometer_phi + sample_phi)
    c_phi = numpy.cos(goniometer_phi + sample_phi)
    s_s_omega, c_s_omega = numpy.sin(sample_omega), numpy.cos(sample_omega)
    s_s_chi, c_s_chi = numpy.sin(sample_chi), numpy.cos(sample_chi)

    a_t = -s_theta * c_omega + c_theta * c_azimuth * s_omega
    b_t = s_theta * s_omega + c_theta * c_azimuth * c_omega
    c_t = -c_theta * s_azimuth * s_chi
    d_t = -c_theta * s_azimuth * c_chi

    b_c1 = -a_t * s_phi + (c_t + b_t * c_chi) * c_phi
    b_c2 = d_t - b_t * s_chi
    b_c3 = a_t * c_phi + (c_t + b_t * c_chi) * s_phi

    b_c = b_c1 * c_s_omega * c_s_chi + b_c2 * c_s_omega * s_s_chi - b_c3 * s_s_omega
    psi = numpy.arccos(numpy.clip(b_c, -1., 1.))
    b_a = -b_c1 * s_s_chi + b_c2 * c_s_chi
    b_b = b_c1 * s_s_omega * c_s_chi + b_c2 * s_s_omega * s_s_chi + b_c3 * c_s_omega
    gamma = numpy.arctan2(b_b, b_a)
    return psi, gamma


def calc_texture_sh_harmonics(
        index_hkl, unit_cell_parameters, laue_class: str, unique_axis: str,
        sample_symmetry: str, texture_sh_lmn, psi, gamma):
    """Harmonics of the texture correction, kept apart to save memory.

    The correction is A = 1 + sum_lm k_l^m(y) sum_n 4 pi/(2l+1) C_l^mn
    k_l^n(h): the sample harmonics are evaluated once for every distinct
    (l, m) (one per order for the cylindrical symmetry), not for every
    coefficient, which matters for the points of a 2D detector.

    Parameters
    ----------
    index_hkl : array_like
        Miller indexes, shape (3, n_hkl).
    unit_cell_parameters : array_like
        a, b, c (Å), alpha, beta, gamma (radians).
    laue_class : str
        Laue class of the crystal (one of LAUE_CLASSES).
    unique_axis : str
        Unique axis of the monoclinic class ("" otherwise).
    sample_symmetry : str
        Sample symmetry (one of SAMPLE_SYMMETRIES).
    texture_sh_lmn : array_like of int
        Indexes l, m, n of the coefficients, shape (3, n_coefficient).
    psi, gamma : array_like
        Direction of the scattering vector in the sample frame, radians;
        the shape (..., n_hkl) (the last axis is the reflection).

    Returns
    -------
    k_crystal : numpy.ndarray
        4 pi/(2l+1) k_l^n(h) of the coefficients, shape (n_coefficient, n_hkl).
    index_coefficient_lm : numpy.ndarray of int
        Row of k_sample of every coefficient, shape (n_coefficient,).
    k_sample : numpy.ndarray
        k_l^m(y) of the distinct (l, m), shape (n_lm,) + psi.shape.
    """
    texture_sh_lmn = numpy.asarray(texture_sh_lmn, dtype=int).reshape(3, -1)
    polar, azimuth = calc_crystal_angles(index_hkl, unit_cell_parameters, laue_class,
                                         unique_axis)
    k_crystal = calc_crystal_harmonics(texture_sh_lmn[[0, 2]], laue_class, polar, azimuth)
    k_crystal = k_crystal * (4. * numpy.pi / (2. * texture_sh_lmn[0] + 1.))[:, na]
    index_lm, index_coefficient_lm = numpy.unique(
        texture_sh_lmn[[0, 1]], axis=1, return_inverse=True)
    k_sample = calc_sample_harmonics(index_lm, sample_symmetry, psi, gamma)
    return k_crystal, index_coefficient_lm.reshape(-1), k_sample


def calc_texture_sh_correction_by_harmonics(
        k_crystal, index_coefficient_lm, k_sample, texture_sh_coefficient):
    """Texture correction from the harmonics of calc_texture_sh_harmonics.

    Parameters
    ----------
    k_crystal : numpy.ndarray
        Shape (n_coefficient, n_hkl).
    index_coefficient_lm : numpy.ndarray of int
        Shape (n_coefficient,).
    k_sample : numpy.ndarray
        Shape (n_lm, ..., n_hkl).
    texture_sh_coefficient : array_like
        Coefficients C_l^mn, shape (n_coefficient,).

    Returns
    -------
    numpy.ndarray
        Correction A of the shape k_sample.shape[1:].
    """
    texture_sh_coefficient = numpy.asarray(texture_sh_coefficient, dtype=float)
    weight_lm = numpy.zeros((k_sample.shape[0], k_crystal.shape[1]), dtype=float)
    numpy.add.at(weight_lm, index_coefficient_lm, texture_sh_coefficient[:, na] * k_crystal)
    correction = numpy.ones(k_sample.shape[1:], dtype=float)
    for k_sample_lm, weight in zip(k_sample, weight_lm):
        correction += k_sample_lm * weight
    return correction


def calc_texture_sh_correction(
        index_hkl, unit_cell_parameters, laue_class: str, unique_axis: str,
        sample_symmetry: str, texture_sh_lmn, texture_sh_coefficient, psi, gamma,
        flag_texture_sh_coefficient: bool = False):
    """Texture correction A(h, y) of the integrated intensities.

    A = 1 + sum 4 pi / (2l + 1) C_l^mn k_l^m(y) k_l^n(h) (GSAS-II SHTXcal).

    Parameters
    ----------
    index_hkl : array_like
        Miller indexes, shape (3, n_hkl).
    unit_cell_parameters : array_like
        a, b, c (Å), alpha, beta, gamma (radians).
    laue_class : str
        Laue class of the crystal (one of LAUE_CLASSES).
    unique_axis : str
        Unique axis of the monoclinic class ("" otherwise).
    sample_symmetry : str
        Sample symmetry (one of SAMPLE_SYMMETRIES).
    texture_sh_lmn : array_like of int
        Indexes l, m, n of the coefficients, shape (3, n_coefficient).
    texture_sh_coefficient : array_like
        Coefficients C_l^mn, shape (n_coefficient,).
    psi, gamma : array_like
        Direction of the scattering vector in the sample frame, radians;
        the shape (..., n_hkl) (the last axis is the reflection).
    flag_texture_sh_coefficient : bool
        Calculate the derivatives over the coefficients.

    Returns
    -------
    correction : numpy.ndarray
        Texture correction of the shape of psi.
    dder : dict
        "texture_sh_coefficient": derivatives of the shape
        psi.shape + (n_coefficient,) when requested.
    """
    k_crystal, index_coefficient_lm, k_sample = calc_texture_sh_harmonics(
        index_hkl, unit_cell_parameters, laue_class, unique_axis, sample_symmetry,
        texture_sh_lmn, psi, gamma)
    correction = calc_texture_sh_correction_by_harmonics(
        k_crystal, index_coefficient_lm, k_sample, texture_sh_coefficient)
    dder = {}
    if flag_texture_sh_coefficient:
        # the correction is linear: the derivative is the term of the coefficient
        terms = k_sample[index_coefficient_lm] * k_crystal.reshape(
            (k_crystal.shape[0],) + (1,) * (k_sample.ndim - 2) + (k_crystal.shape[1],))
        dder["texture_sh_coefficient"] = numpy.moveaxis(terms, 0, -1)
    return correction, dder


def calc_texture_index(texture_sh_lmn, texture_sh_coefficient) -> float:
    """Texture index J = 1 + sum C_l^mn^2 / (2l + 1) (GSAS-II textureIndex).

    Parameters
    ----------
    texture_sh_lmn : array_like of int
        Indexes l, m, n, shape (3, n_coefficient).
    texture_sh_coefficient : array_like
        Coefficients C_l^mn, shape (n_coefficient,).

    Returns
    -------
    float
        1 for a random powder, larger for a texture.
    """
    texture_sh_lmn = numpy.asarray(texture_sh_lmn, dtype=int).reshape(3, -1)
    texture_sh_coefficient = numpy.asarray(texture_sh_coefficient, dtype=float)
    return float(1. + numpy.sum(numpy.square(texture_sh_coefficient)
                                / (2. * texture_sh_lmn[0] + 1.)))


def calc_laue_class_by_rotations(rotations, unit_cell_parameters):
    """Laue class and unique axis in the frames of the harmonics.

    The class is chosen among the candidates with the same number of
    rotations as the one whose crystal harmonics up to the order 8 are
    invariant under every rotation (this also finds the setting: unique axis
    of 2/m, -3m1 or -31m, hexagonal or rhombohedral axes).

    Parameters
    ----------
    rotations : array_like
        Rotations of the Laue class (with -1), shape (n, 3, 3), acting on hkl
        as h' = R^T h (rhochi_absorption.calc_laue_rotations).
    unit_cell_parameters : array_like
        a, b, c (Å), alpha, beta, gamma (radians).

    Returns
    -------
    laue_class, unique_axis : str
        A Laue class of LAUE_CLASSES and the unique axis ("" except for 2/m).
    """
    rotations = numpy.asarray(rotations, dtype=float)
    candidates = LAUE_CANDIDATES.get(rotations.shape[0], ())
    index_hkl = numpy.array([[0.31, 0.17, 0.83], [-0.52, 0.94, 0.11], [0.73, -0.28, 0.46],
                             [0.12, 0.65, -0.91], [1.0, 0.2, 0.4]], dtype=float).transpose()
    for laue_class, unique_axis in candidates:
        lmn = form_texture_sh_lmn(laue_class, "cylindrical", 8)
        index_ln = lmn[[0, 2]]
        polar, azimuth = calc_crystal_angles(index_hkl, unit_cell_parameters, laue_class,
                                             unique_axis)
        k_ref = calc_crystal_harmonics(index_ln, laue_class, polar, azimuth)
        scale = numpy.abs(k_ref).max()
        flag_invariant = True
        for rotation in rotations:
            index_hkl_rot = numpy.matmul(rotation.transpose(), index_hkl)
            polar, azimuth = calc_crystal_angles(index_hkl_rot, unit_cell_parameters,
                                                 laue_class, unique_axis)
            k_rot = calc_crystal_harmonics(index_ln, laue_class, polar, azimuth)
            if numpy.abs(k_rot - k_ref).max() > 1e-6 * scale:
                flag_invariant = False
                break
        if flag_invariant:
            return laue_class, unique_axis
    raise ValueError(
        "The Laue class of the phase is not found in the settings of the spherical "
        "harmonics texture (standard settings: 4-, 3- and 6-fold axes along c, "
        "rhombohedral 3-fold along [111]).")
