from typing import NoReturn
import numpy
import warnings

from cryspy.A_functions_base.function_1_atomic_vibrations import \
    apply_constraint_on_cell_by_type_cell
from cryspy.A_functions_base.function_2_crystallography_base import \
    calc_sthovl_by_hkl_abc_cosines, ortogonalize_matrix
from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN


class Cell(ItemN):
    """Unit cell description.

    Attributes
    ----------
        - length_a, length_b, length_c, angle_alpha, angle_beta, angle_gamma
          (mandatory)
        - m_g, m_g_reciprocal, m_b, m_b_norm, m_m, m_m_norm,
          cos_a, cos_b, cos_g, sin_a, sin_b, sin_g, cos_a_sq,
          cos_b_sq, cos_g_sq, sin_a_sq, sin_b_sq, sin_g_sq, cos_ia,
          cos_ib, cos_ig, sin_ia, sin_ib, sin_ig, cos_ia_sq,
          cos_ib_sq, cos_ig_sq, sin_ia_sq, sin_ib_sq, sin_ig_sq,
          reciprocal_length_a, reciprocal_length_b, reciprocal_length_c,
          reciprocal_angle_alpha, reciprocal_angle_beta,
          reciprocal_angle_gamma, volume (internal)
        - type_cell, it_coordinate_system_code, formula_units_z
          (internal, protected)
    """
    ATTR_MANDATORY_NAMES = ("length_a", "length_b", "length_c",
                            "angle_alpha", "angle_beta", "angle_gamma")
    ATTR_MANDATORY_TYPES = (float, float, float, float, float, float)
    ATTR_MANDATORY_CIF = ("length_a", "length_b", "length_c",
                          "angle_alpha", "angle_beta", "angle_gamma")

    ATTR_OPTIONAL_NAMES = ()
    ATTR_OPTIONAL_TYPES = ()
    ATTR_OPTIONAL_CIF = ()

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = (
        "m_g", "m_g_reciprocal", "m_b", "m_b_norm", "m_m", "m_m_norm",
        "cos_a", "cos_b", "cos_g", "sin_a", "sin_b", "sin_g", "cos_a_sq",
        "cos_b_sq", "cos_g_sq", "sin_a_sq", "sin_b_sq", "sin_g_sq", "cos_ia",
        "cos_ib", "cos_ig", "sin_ia", "sin_ib", "sin_ig", "cos_ia_sq",
        "cos_ib_sq", "cos_ig_sq", "sin_ia_sq", "sin_ib_sq", "sin_ig_sq",
        "reciprocal_length_a", "reciprocal_length_b", "reciprocal_length_c",
        "reciprocal_angle_alpha", "reciprocal_angle_beta",
        "reciprocal_angle_gamma", "volume")
    ATTR_INT_PROTECTED_NAMES = ("type_cell", "it_coordinate_system_code",
                                "formula_units_z")

    # parameters considered are refined parameters
    ATTR_REF = ('length_a', 'length_b', 'length_c',
                'angle_alpha', 'angle_beta', 'angle_gamma')
    ATTR_SIGMA = tuple([f"{_h:}_sigma" for _h in ATTR_REF])
    ATTR_CONSTR_FLAG = tuple([f"{_h:}_constraint" for _h in ATTR_REF])
    ATTR_REF_FLAG = tuple([f"{_h:}_refinement" for _h in ATTR_REF])
    ATTR_CONSTR_MARK = tuple([f"{_h:}_mark" for _h in ATTR_REF])

    # formats if cif format
    D_FORMATS = {"length_a": "{:.6f}", "length_b": "{:.6f}",
                 "length_c": "{:.6f}", "angle_alpha": "{:.6f}",
                 "angle_beta": "{:.6f}", "angle_gamma": "{:.6f}"}

    # constraints on the parameters
    D_CONSTRAINTS = {'type': ["c", "p"]}

    # default values for the parameters
    D_DEFAULT = {'length_a': 1., 'length_b': 1., 'length_c': 1.,
                 'angle_alpha': 90., 'angle_beta': 90., 'angle_gamma': 90.}
    for key in ATTR_SIGMA:
        D_DEFAULT[key] = 0.
    for key in (ATTR_CONSTR_FLAG + ATTR_REF_FLAG):
        D_DEFAULT[key] = False
    for key in ATTR_CONSTR_MARK:
        D_DEFAULT[key] = ""

    PREFIX = "cell"

    def __init__(self, **kwargs) -> NoReturn:
        super(Cell, self).__init__()

        # defined for any integer and float parameters
        D_MIN = {'legnth_a': 0, 'legnth_b': 0, 'legnth_c': 0,
                 'angle_alpha': 0., 'angle_beta': 0., 'angle_gamma': 0.}

        # defined for ani integer and float parameters
        D_MAX = {'angle_alpha': 180., 'angle_beta': 180., 'angle_gamma': 180.}

        self.__dict__["D_MIN"] = D_MIN
        self.__dict__["D_MAX"] = D_MAX
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)

    def apply_constraints_by_type_cell_it_coordinate_system_code(
            self, type_cell: str, it_coordinate_system_code: str) -> NoReturn:
        """Apply constraints."""
        self.type_cell = type_cell
        self.it_coordinate_system_code = it_coordinate_system_code
        self.form_object()

    def apply_constraints(self):
        """Apply constraints."""
        flag_1 = not(self.is_attribute("type_cell"))
        flag_2 = not(self.is_attribute("it_coordinate_system_code"))
        if (flag_1 | flag_2):
            return
        cell_p = (self.length_a, self.length_b, self.length_c,
                  self.angle_alpha, self.angle_beta, self.angle_gamma)
        cell_s = (self.length_a_sigma, self.length_b_sigma,
                  self.length_c_sigma, self.angle_alpha_sigma,
                  self.angle_beta_sigma, self.angle_gamma_sigma)
        cell_r = (self.length_a_refinement, self.length_b_refinement,
                  self.length_c_refinement, self.angle_alpha_refinement,
                  self.angle_beta_refinement, self.angle_gamma_refinement)
        cell_p, cell_s, cell_r, cell_constr = \
            apply_constraint_on_cell_by_type_cell(
                cell_p, cell_s, cell_r, self.type_cell,
                self.it_coordinate_system_code)
        self.__dict__['length_a'] = cell_p[0]
        self.__dict__['length_b'] = cell_p[1]
        self.__dict__['length_c'] = cell_p[2]
        self.__dict__['angle_alpha'] = cell_p[3]
        self.__dict__['angle_beta'] = cell_p[4]
        self.__dict__['angle_gamma'] = cell_p[5]
        self.__dict__['length_a_sigma'] = cell_s[0]
        self.__dict__['length_b_sigma'] = cell_s[1]
        self.__dict__['length_c_sigma'] = cell_s[2]
        self.__dict__['angle_alpha_sigma'] = cell_s[3]
        self.__dict__['angle_beta_sigma'] = cell_s[4]
        self.__dict__['angle_gamma_sigma'] = cell_s[5]
        self.__dict__['length_a_refinement'] = cell_r[0]
        self.__dict__['length_b_refinement'] = cell_r[1]
        self.__dict__['length_c_refinement'] = cell_r[2]
        self.__dict__['angle_alpha_refinement'] = cell_r[3]
        self.__dict__['angle_beta_refinement'] = cell_r[4]
        self.__dict__['angle_gamma_refinement'] = cell_r[5]
        self.__dict__['length_a_constraint'] = cell_constr[0]
        self.__dict__['length_b_constraint'] = cell_constr[1]
        self.__dict__['length_c_constraint'] = cell_constr[2]
        self.__dict__['angle_alpha_constraint'] = cell_constr[3]
        self.__dict__['angle_beta_constraint'] = cell_constr[4]
        self.__dict__['angle_gamma_constraint'] = cell_constr[5]
        self.form_object(flag_constraint=False)
        return

    def form_object(self, flag_constraint: bool = True):
        """Form object."""
        if flag_constraint:
            self.apply_constraints()

        rad = numpy.pi/180.

        c_a = numpy.cos(float(self.angle_alpha)*rad)
        c_b = numpy.cos(float(self.angle_beta)*rad)
        c_g = numpy.cos(float(self.angle_gamma)*rad)
        s_a = numpy.sin(float(self.angle_alpha)*rad)
        s_b = numpy.sin(float(self.angle_beta)*rad)
        s_g = numpy.sin(float(self.angle_gamma)*rad)
        self.cos_a = c_a
        self.cos_b = c_b
        self.cos_g = c_g
        self.sin_a = s_a
        self.sin_b = s_b
        self.sin_g = s_g

        c_a_sq, c_b_sq, c_g_sq = c_a**2, c_b**2, c_g**2
        s_a_sq, s_b_sq, s_g_sq = (1.-c_a_sq), (1.-c_b_sq), (1.-c_g_sq)
        self.cos_a_sq = c_a_sq
        self.cos_b_sq = c_b_sq
        self.cos_g_sq = c_g_sq
        self.sin_a_sq = s_a_sq
        self.sin_b_sq = s_b_sq
        self.sin_g_sq = s_g_sq

        a, b, c = self.length_a, self.length_b, self.length_c

        v_sqrt = (1.-c_a_sq-c_b_sq-c_g_sq+2.*c_a*c_b*c_g)**0.5
        vol = a*b*c*v_sqrt
        self.volume = vol

        # G matrix
        m_g = numpy.array([[a*a, a*b*c_g, a*c*c_b],
                           [a*b*c_g, b*b, b*c*c_a],
                           [a*c*c_b, b*c*c_a, c*c]], dtype=float)
        self.m_g = m_g

        # nM matrix
        m_m_norm = numpy.array([[v_sqrt/s_a, 0, 0],
                                [(c_g-c_a*c_b)/s_a, s_a, 0],
                                [c_b, c_a, 1]], dtype=float)
        self.m_m_norm = m_m_norm

        irad = 180./numpy.pi
        ialpha = numpy.arccos((c_b*c_g-c_a)/(s_b*s_g))*irad
        ibeta = numpy.arccos((c_g*c_a-c_b)/(s_g*s_a))*irad
        igamma = numpy.arccos((c_a*c_b-c_g)/(s_a*s_b))*irad
        ia, ib, ic = b*c*s_a/vol, c*a*s_b/vol, a*b*s_g/vol

        self.reciprocal_length_a = ia
        self.reciprocal_length_b = ib
        self.reciprocal_length_c = ic
        self.reciprocal_angle_alpha = ialpha
        self.reciprocal_angle_beta = ibeta
        self.reciprocal_angle_gamma = igamma

        c_ia, c_ib, c_ig = numpy.cos(ialpha*rad), numpy.cos(ibeta*rad), \
            numpy.cos(igamma*rad)
        s_ia, s_ib, s_ig = numpy.sin(ialpha*rad), numpy.sin(ibeta*rad), \
            numpy.sin(igamma*rad)
        self.cos_ia = c_ia
        self.cos_ib = c_ib
        self.cos_ig = c_ig
        self.sin_ia = s_ia
        self.sin_ib = s_ib
        self.sin_ig = s_ig

        c_ia_sq, c_ib_sq, c_ig_sq = c_ia**2, c_ib**2, c_ig**2
        s_ia_sq, s_ib_sq, s_ig_sq = (1.-c_ia_sq), (1.-c_ib_sq), (1.-c_ig_sq)
        self.cos_ia_sq = c_ia_sq
        self.cos_ib_sq = c_ib_sq
        self.cos_ig_sq = c_ig_sq
        self.sin_ia_sq = s_ia_sq
        self.sin_ib_sq = s_ib_sq
        self.sin_ig_sq = s_ig_sq

        # G matrix for reciprocal unit cell
        m_g_reciprocal = numpy.array([[ia*ia, ia*ib*c_ig, ia*ic*c_ib],
                                      [ia*ib*c_ig, ib*ib, ib*ic*c_ia],
                                      [ia*ic*c_ib, ib*ic*c_ia, ic*ic]],
                                     dtype=float)
        self.m_g_reciprocal = m_g_reciprocal

        # nB matrix
        m_b_norm = numpy.array([[ 1, c_ig, c_ib],
                                [ 0, s_ig, -s_ib*c_a],
                                [ 0, 0, v_sqrt/s_g]], dtype=float)
        self.m_b_norm = m_b_norm

        # B matrix
        m_b = numpy.array([[ia, ib*c_ig, ic*c_ib],
                           [0., ib*s_ig, -ic*s_ib*c_a],
                           [0., 0., 1./c]], dtype=float)
        self.m_b = m_b

        # M matrix j
        m_m = numpy.array([[ 1./ia, 0,0],
                           [-1*a*s_b*c_ig, b*s_a, 0],
                           [ a*c_b, b*c_a, c]], dtype=float)
        self.m_m = m_m

    def calc_sthovl(self, index_h: int, index_k: int, index_l: int):
        """
        Calculate sin(theta)/lambda for list of hkl reflections.

        Parameters
        ----------
        h : int
            Miller index h.
        k : int
            Miller index k.
        l : int
            Miller index l.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        return calc_sthovl_by_hkl_abc_cosines(
            index_h, index_k, index_l, self.length_a, self.length_b,
            self.length_c, self.cos_a, self.cos_b, self.cos_g)

    def calc_k_loc(self, index_h, index_k, index_l):
        """
        Calculate unity scattering vector in Cartesian coordinate system.
        Coordinate system is defined as (x||a*, z||c).

        Keyword arguments:

            h, k, l: Miller indices

        Output arguments:

            k_x, k_y, k_z: 1D numpy array of x, y, z components of unity
                           scattering vector
        """
        np_h = numpy.array(index_h, dtype=float)
        np_k = numpy.array(index_k, dtype=float)
        np_l = numpy.array(index_l, dtype=float)
        m_b = self.m_b
        k_x = m_b[0, 0]*np_h + m_b[0, 1]*np_k +m_b[0, 2]*np_l
        k_y = m_b[1, 0]*np_h + m_b[1, 1]*np_k +m_b[1, 2]*np_l
        k_z = m_b[2, 0]*np_h + m_b[2, 1]*np_k +m_b[2, 2]*np_l

        k_norm = (k_x**2 + k_y**2 + k_z**2)**0.5
        if not((type(index_h) is float) | (type(index_h) is int) |
               (type(index_h) is numpy.float64)):
            k_norm[k_norm == 0.] = 1.
        elif k_norm == 0.:
            k_norm = 1.

        k_x = k_x/k_norm
        k_y = k_y/k_norm
        k_z = k_z/k_norm

        return k_x, k_y, k_z

    def calc_matrix_t(self, index_h, index_k, index_l):
        """Determine rotation matrix to have new z axis along k_loc.

        Rotation matrix is defined by Euler angles
        Attention
        m_(x||a*, z||c) = T * M_(Z||hkl)
        """
        m_b = self.m_b
        k_x = m_b[0, 0] * index_h + m_b[0, 1] * index_k + m_b[0, 2] * index_l
        k_y = m_b[1, 0] * index_h + m_b[1, 1] * index_k + m_b[1, 2] * index_l
        k_z = m_b[2, 0] * index_h + m_b[2, 1] * index_k + m_b[2, 2] * index_l

        k_norm = (k_x**2 + k_y**2 + k_z**2)**0.5
        k_norm[k_norm == 0.] = 1.

        k_x = k_x/k_norm
        k_y = k_y/k_norm
        k_z = k_z/k_norm

        al = numpy.zeros(k_x.shape, dtype=float)

        be = numpy.arccos(k_z)
        sb = numpy.sin(be)
        flag = (sb != 0.)

        sa1 = k_x[flag]*1./sb[flag]
        ca2 = -1*k_y[flag]*1./sb[flag]
        sa1[sa1 > 1] = 1.
        sa1[sa1 < -1] = -1.

        ca2[ca2 > 1] = 1.
        ca2[ca2 < -1] = -1.

        al1 = numpy.arcsin(sa1)
        al2 = numpy.arccos(ca2)

        al_sh = numpy.copy(al1)
        al_sh[sa1 > 0.] = al2[sa1 > 0.]
        al_sh[sa1 <= 0.] = 2.*numpy.pi-al2[sa1 <= 0.]
        al_sh[numpy.abs(al2-al1) < 0.00001] = al1[numpy.abs(al2-al1) < 0.00001]

        al[flag] = al_sh

        ga = 0.
        ca, cb, cg = numpy.cos(al), numpy.cos(be), numpy.cos(ga)
        sa, sb, sg = numpy.sin(al), numpy.sin(be), numpy.sin(ga)

        # FIXME: I would like to recheck the expression for T.
        t_11, t_12, t_13 = ca*cg-sa*cb*sg, -ca*sg-sa*cb*cg,  sa*sb
        t_21, t_22, t_23 = sa*cg+ca*cb*sg, -sa*sg+ca*cb*cg, -ca*sb
        t_31, t_32, t_33 = sb*sg, sb*cg, cb

        flag = (((sa*sb-k_x)**2+(-ca*sb-k_y)**2+(cb-k_z)**2) > 0.0001)
        if any(flag):
            warnings.warn("Mistake with k_loc\nProgram is stopped",
                          UserWarning, stacklevel=2)
            quit()
        return t_11, t_12, t_13, t_21, t_22, t_23, t_31, t_32, t_33

    def calc_hkl(self, space_group, sthovl_min, sthovl_max):
        """
        A list of reflections hkl for cell in the range sthovl_min, sthovl_max
        taking into account the space group
        """
        lhkl, lmult = [], []
        l_hklres = []

        hmax = int(2.*self.length_a*sthovl_max)
        kmax = int(2.*self.length_b*sthovl_max)
        lmax = int(2.*self.length_c*sthovl_max)
        hmin, kmin, lmin = -1*hmax, -1*kmax, -1*lmax

        hmin = 0
        shift = space_group.shift
        r_s_g_s = space_group.reduced_space_group_symop
        orig_x, orig_y, orig_z = zip(*shift)
        np_orig_x = numpy.array(orig_x, dtype=float)
        np_orig_y = numpy.array(orig_y, dtype=float)
        np_orig_z = numpy.array(orig_z, dtype=float)

        r_11 = numpy.array(r_s_g_s.r_11, dtype=float)
        r_12 = numpy.array(r_s_g_s.r_12, dtype=float)
        r_13 = numpy.array(r_s_g_s.r_13, dtype=float)
        r_21 = numpy.array(r_s_g_s.r_21, dtype=float)
        r_22 = numpy.array(r_s_g_s.r_22, dtype=float)
        r_23 = numpy.array(r_s_g_s.r_23, dtype=float)
        r_31 = numpy.array(r_s_g_s.r_31, dtype=float)
        r_32 = numpy.array(r_s_g_s.r_32, dtype=float)
        r_33 = numpy.array(r_s_g_s.r_33, dtype=float)

        for h in range(hmin, hmax+1, 1):
            for k in range(kmin, kmax+1, 1):
                for l in range(lmin, lmax+1, 1):
                    flag=(abs(sum(numpy.exp(2.*numpy.pi*1j*(np_orig_x*h+np_orig_y*k+np_orig_z*l))))>0.00001)
                    if (flag):
                        lhkls = [(_1, _2, _3) for _1, _2, _3 in zip(h*r_11+k*r_21+l*r_31, h*r_12+k*r_22+l*r_32, h*r_13+k*r_23+l*r_33)]
                        lhkls.extend([(-hkl[0],-hkl[1],-hkl[2]) for hkl in lhkls])
                        lhkls.sort(key=lambda x:10000*x[0]+100*x[1]+x[2])
                        if (not(lhkls[-1] in lhkl)):
                            lhkl.append(lhkls[-1])
                            lmult.append(len(set(lhkls)))

        l_hklsthovl=[(hkl, self.calc_sthovl(hkl[0], hkl[1], hkl[2]), mult) for hkl, mult in zip(lhkl, lmult)]
        l_hklsthovl.sort(key=lambda x: x[1])
        l_hklres = [hklsthovl[0] for hklsthovl in l_hklsthovl if ((hklsthovl[1]>sthovl_min) & (hklsthovl[1]<sthovl_max))]
        l_multres = [hklsthovl[2] for hklsthovl in l_hklsthovl if ((hklsthovl[1]>sthovl_min) & (hklsthovl[1]<sthovl_max))]

        h = numpy.array([hh[0] for hh in l_hklres], dtype=int)
        k = numpy.array([hh[1] for hh in l_hklres], dtype=int)
        l = numpy.array([hh[2] for hh in l_hklres], dtype=int)
        mult = numpy.array(l_multres, dtype=int)
        return h, k, l, mult

    def calc_hkl_in_range(self, sthovl_min, sthovl_max):
        """
        Give a list of reflections hkl for cell in the range.
        sthovl_min, sthovl_max
        """
        h_max = int(2.*self.length_a*sthovl_max)
        k_max = int(2.*self.length_b*sthovl_max)
        l_max = int(2.*self.length_c*sthovl_max)
        h_min, k_min, l_min = -1*h_max, -1*k_max, -1*l_max

        np_h = numpy.array(range(h_min, h_max+1, 1), dtype=int)
        np_k = numpy.array(range(k_min, k_max+1, 1), dtype=int)
        np_l = numpy.array(range(l_min, l_max+1, 1), dtype=int)
        h_3d, k_3d, l_3d = numpy.meshgrid(np_h, np_k, np_l, indexing="ij")

        sthovl_3d = self.calc_sthovl(h_3d, k_3d, l_3d)
        flag_1 = sthovl_3d >= sthovl_min
        flag_2 = sthovl_3d <= sthovl_max
        flag_12 = numpy.logical_and(flag_1, flag_2)

        h = h_3d[flag_12]
        k = k_3d[flag_12]
        l = l_3d[flag_12]
        mult = numpy.ones(h.size, dtype=int)
        sthovl = sthovl_3d[flag_12]
        arg_sort = numpy.argsort(sthovl)
        return h[arg_sort], k[arg_sort], l[arg_sort], mult[arg_sort] 

    def calc_position_by_coordinate(self, x ,y, z):
        """
        Calculates position for coordinate :math:`(x,y,z)`
        in Carthezian coordinate system in which (x||a*, z||c)
        """
        m_m = self.m_m
        p_x = m_m[0, 0]*x + m_m[0, 1]*y + m_m[0, 2]*z
        p_y = m_m[1, 0]*x + m_m[1, 1]*y + m_m[1, 2]*z
        p_z = m_m[2, 0]*x + m_m[2, 1]*y + m_m[2, 2]*z
        return p_x, p_y, p_z

    def calc_coordinate_by_position(self, p_x, p_y, p_z):
        """
        Calculates coordinate :math:`(x,y,z)` for  position defined
        in Carthezian coordinate system in which (x||a*, z||c)
        """
        m_m = self.m_m
        m_im = numpy.linalg.inv(m_m)
        x = m_im[0, 0]*p_x + m_im[0, 1]*p_y + m_im[0, 2]*p_z
        y = m_im[1, 0]*p_x + m_im[1, 1]*p_y + m_im[1, 2]*p_z
        z = m_im[2, 0]*p_x + m_im[2, 1]*p_y + m_im[2, 2]*p_z
        return x, y, z

    def calc_length_sq(self, x, y, z):
        """
        According to IT_C Section 1.1.2 Lattice vectors, point rows and net planes

        .. math::
        t^{2} = x^{2} a^{2} + y^{2} b^{2} + z^{2} c^{2} 
            + 2 x y a b \\cos \\gamma  + 2 x z a c \\cos \\beta  + 2 y z b c \\cos \\alpha 
        """
        a, b, c = self.length_a, self.length_b, self.length_c
        c_a, c_b, c_g = self.cos_a, self.cos_b, self.cos_g
        t_sq = (x**2 * a**2 + y**2 * b**2 + z**2 * c**2 + 
                2.*x*y*a*b*c_g + 2.*x*z*a*c*c_b + 2.*y*z*b*c*c_a)
        return t_sq

    def closest_distance_between_fractions(self, fract_1: numpy.ndarray, fract_2: numpy.ndarray) -> numpy.ndarray:
        """
        Give closest distance between two fractions

        The shape of input arrays [3, n]
        """
        if len(fract_1.shape) == 1:
            fract_1_p = fract_1.reshape(fract_1.size, 1)
        else:
            fract_1_p = fract_1
        if len(fract_2.shape) == 1:
            fract_2_p = fract_2.reshape(fract_2.size, 1)
        else:
            fract_2_p = fract_2
        fract_1_p, fract_2_p = numpy.mod(fract_1_p, 1.), numpy.mod(fract_2_p, 1.)
        diff_fract = numpy.abs(fract_1_p - fract_2_p)
        diff_closest = numpy.where(diff_fract < 0.5,
                                   diff_fract, 1.0 - diff_fract)
        diff_pos = self.calc_position_by_coordinate(
            diff_closest[0, :], diff_closest[1, :], diff_closest[2, :])
        distance = numpy.sqrt((numpy.square(diff_pos)).sum(axis=0))
        return distance

    def calc_reciprocal_length_sq(self, h, k, l):
        """
        According to IT_C Section 1.1.2 Lattice vectors, point rows and net planes

        .. math::
            r^{*^{2}} = h^{2} a^{*^{2}} + k^{*^{2}} b^{*^{2}} + l^{2} c^{*^{2}} 
            + 2 h k a^{*} b^{*} \\cos \\gamma  + 2 h l a^{*} c^{*} \\cos \\beta  + 2 k l b^{*} c^{*} \\cos \\alpha 
        """
        a, b, c = self.reciprocal_length_a, self.reciprocal_length_b, self.reciprocal_length_c
        c_a, c_b, c_g = self.cos_ia, self.cos_ib, self.cos_ig
        x, y, z = h, k, l
        t_sq = (x**2 * a**2 + y**2 * b**2 + z**2 * c**2 + 
                2.*x*y*a*b*c_g + 2.*x*z*a*c*c_b + 2.*y*z*b*c*c_a)
        return t_sq

    def ortogonalize_matrix(self, m_ij):
        """
        matrix m_ij is defined in coordinate system (a, b, c)

        output matrix s_ij is defined in Cartezian coordinate system defined as x||a*, z||c, y= [z x] (right handed)
        """
        m_m_norm = self.m_m_norm

        m_norm_ij = (m_m_norm[0, 0], m_m_norm[0, 1], m_m_norm[0, 2],
                     m_m_norm[1, 0], m_m_norm[1, 1], m_m_norm[1, 2],
                     m_m_norm[2, 0], m_m_norm[2, 1], m_m_norm[2, 2])

        s_11, s_12, s_13, s_21, s_22, s_23, s_31, s_32, s_33 = \
            ortogonalize_matrix(m_ij, m_norm_ij)

        return s_11, s_12, s_13, s_21, s_22, s_23, s_31, s_32, s_33

    def report(self) -> str:
        """
        Make a report about cell object in string format.
        """
        ls_out = []
        ls_out.append(f"|Volume:| {self.volume: 9.2f} Ang.**3|")

        ls_out.append("\n## Reciprocal unit cell (1/Ang.)\n")
        ls_out.append(f"|{float(self.reciprocal_length_a):9.5f} |{float(self.reciprocal_length_b):9.5f} |{float(self.reciprocal_length_c):9.5f}|")
        ls_out.append(f"|{float(self.reciprocal_angle_alpha):9.3f} |{float(self.reciprocal_angle_beta):9.3f} |{float(self.reciprocal_angle_gamma):9.3f}|")

        ls_out.append("\n## B matrix (x||a*, z||c)\n")
        ls_out.append(f"|{self.m_b[0, 0]: 10.5f}|{self.m_b[0, 1]: 10.5f}|{self.m_b[0, 2]: 10.5f}|")
        ls_out.append(f"|{self.m_b[1, 0]: 10.5f}|{self.m_b[1, 1]: 10.5f}|{self.m_b[1, 2]: 10.5f}|")
        ls_out.append(f"|{self.m_b[2, 0]: 10.5f}|{self.m_b[2, 1]: 10.5f}|{self.m_b[2, 2]: 10.5f}|")

        ls_out.append("\n## M matrix (x||a*, z||c)\n")
        ls_out.append(f"|{self.m_m[0, 0]: 10.5f}|{self.m_m[0, 1]: 10.5f}|{self.m_m[0, 2]: 10.5f}|")
        ls_out.append(f"|{self.m_m[1, 0]: 10.5f}|{self.m_m[1, 1]: 10.5f}|{self.m_m[1, 2]: 10.5f}|")
        ls_out.append(f"|{self.m_m[2, 0]: 10.5f}|{self.m_m[2, 1]: 10.5f}|{self.m_m[2, 2]: 10.5f}|")

        ls_out.append("\n## Metric tensor for unit cell\n")
        ls_out.append(f"|{self.m_g[0, 0]: 10.5f}|{self.m_g[0, 1]: 10.5f}|{self.m_g[0, 2]: 10.5f}|")
        ls_out.append(f"|{self.m_g[1, 0]: 10.5f}|{self.m_g[1, 1]: 10.5f}|{self.m_g[1, 2]: 10.5f}|")
        ls_out.append(f"|{self.m_g[2, 0]: 10.5f}|{self.m_g[2, 1]: 10.5f}|{self.m_g[2, 2]: 10.5f}|")

        ls_out.append("\n## Metric tensor for reciprocal unit cell\n")
        ls_out.append(f"|{self.m_g_reciprocal[0, 0]: 10.5f}|{self.m_g_reciprocal[0, 1]: 10.5f}|{self.m_g_reciprocal[0, 2]: 10.5f}|")
        ls_out.append(f"|{self.m_g_reciprocal[1, 0]: 10.5f}|{self.m_g_reciprocal[1, 1]: 10.5f}|{self.m_g_reciprocal[1, 2]: 10.5f}|")
        ls_out.append(f"|{self.m_g_reciprocal[2, 0]: 10.5f}|{self.m_g_reciprocal[2, 1]: 10.5f}|{self.m_g_reciprocal[2, 2]: 10.5f}|")

        return "\n".join(ls_out)

    def get_unit_cell_parameters(self):
        unit_cell_parameters = numpy.array([
            self.length_a, self.length_b, self.length_c, 
            self.angle_alpha*numpy.pi/180., self.angle_beta*numpy.pi/180.,
            self.angle_gamma*numpy.pi/180.], dtype=float)
        return unit_cell_parameters

    def get_flags_unit_cell_parameters(self):
        flags_unit_cell_parameters = numpy.array([
            self.length_a_refinement, self.length_b_refinement, 
            self.length_c_refinement, self.angle_alpha_refinement,
            self.angle_beta_refinement, self.angle_gamma_refinement], dtype=bool)
        return flags_unit_cell_parameters

                

class CellL(LoopN):
    """
    Description of unit cell in loop.

    """
    ITEM_CLASS = Cell
    ATTR_INDEX = None
    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(CellL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name

# s_cont = """
# _cell_length_a 8.0
# _cell_length_b 3.76(5)
# _cell_length_c 4.0
# _cell_angle_alpha 35.0
# _cell_angle_beta 35.0
# _cell_angle_gamma 40.0
# """

# s_cont = """
# loop_
# _cell_length_a
# _cell_length_b
# _cell_length_c
# _cell_angle_alpha
# _cell_angle_beta
# _cell_angle_gamma
# 8.0  2.76546(5) 4.0 88.0 35.0() 87
# 3.0        3.0 1.0 92(2) 90.0 92
# """

# """
# val_1 = Cell.from_cif(s_cont)
# """

# val_2 = Cell(angle_alpha=750)
# val_2.length_a = 3.
# print(val_2, end="\n\n")

# cell_l = CellL.from_cif(s_cont)
# print(cell_l, end="\n\n")

# print(type(cell_l.length_a), end="\n\n")
# print(cell_l[0].calc_hkl_in_range(0.1, 0.2), end="\n\n")
# print(cell_l[0].report(), end="\n\n")
# print(cell_l[0].calc_sthovl(1,2,3), end="\n\n")
# print(cell_l[0], end="\n\n")
# print(cell_l[0].length_b_as_string, end="\n\n")
