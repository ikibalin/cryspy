"""Description of Crystal class."""
__author__ = 'ikibalin'
__version__ = "2020_08_19"
import math
import numpy

from typing import NoReturn

from cryspy.A_functions_base.function_1_strings import \
    value_error_to_string
from cryspy.A_functions_base.function_1_algebra import calc_m_sigma
from cryspy.A_functions_base.function_2_crystallography_base import \
    calc_phase_by_hkl_xyz_rb, calc_dwf, calc_form_factor_tensor_susceptibility

from cryspy.B_parent_classes.cl_3_data import DataN

from cryspy.C_item_loop_classes.cl_2_space_group import SpaceGroup
from cryspy.C_item_loop_classes.cl_1_cell import Cell
from cryspy.C_item_loop_classes.cl_1_atom_site import AtomSiteL
from cryspy.C_item_loop_classes.cl_1_atom_type import AtomTypeL
from cryspy.C_item_loop_classes.cl_1_atom_site_aniso import \
    AtomSiteAnisoL
from cryspy.C_item_loop_classes.cl_1_refln import ReflnL
from cryspy.C_item_loop_classes.cl_1_atom_site_susceptibility import \
    AtomSiteSusceptibilityL
from cryspy.C_item_loop_classes.cl_1_atom_type_scat import \
    AtomTypeScatL
from cryspy.C_item_loop_classes.cl_2_atom_site_scat import \
    AtomSiteScatL
from cryspy.C_item_loop_classes.cl_1_refln_susceptibility import \
    ReflnSusceptibilityL
from cryspy.C_item_loop_classes.cl_1_atom_local_axes import \
    AtomLocalAxesL
from cryspy.C_item_loop_classes.cl_1_atom_electron_configuration \
    import AtomElectronConfigurationL

from cryspy.D_functions_item_loop.function_1_report_magnetization_ellipsoid \
    import magnetization_ellipsoid_by_u_ij, report_main_axes_of_magnetization_ellipsoids

class Crystal(DataN):
    """
    Crystal structure description.

    Data items in the CRYSTAL category record details about
    crystal structure.

    Methods
    -------
        - calc_b_iso_beta
        - calc_f_nucl
        - calc_susceptibility_moment_tensor
        - calc_magnetic_moments_with_field_loc
        - report_main_axes_of_magnetization_ellipsoids
        - calc_main_axes_of_magnetization_ellipsoids
        - calc_magnetization_ellipsoid
        - calc_hkl_in_range
        - calc_hkl
        - calc_refln_susceptibility
        - calc_refln

    Attributes
    ----------
        - space_group, cell, atom_site (mandatory)
        - atom_type, atom_site_aniso, atom_site_susceptibility,
          atom_site_scat, atom_type_scat, atom_local_axes,
          atom_electron_confiduration (optional)
    """

    CLASSES_MANDATORY = (SpaceGroup, Cell, AtomSiteL)
    CLASSES_OPTIONAL = (
        AtomTypeL, AtomSiteAnisoL, AtomSiteSusceptibilityL, AtomSiteScatL,
        AtomTypeScatL, AtomLocalAxesL, AtomElectronConfigurationL)
    # CLASSES_INTERNAL = ()

    CLASSES = CLASSES_MANDATORY + CLASSES_OPTIONAL

    PREFIX = "crystal"

    # default values for the parameters
    D_DEFAULT = {}

    def __init__(self, data_name=None, **kwargs) -> NoReturn:
        super(Crystal, self).__init__()

        self.__dict__["items"] = []
        self.__dict__["data_name"] = data_name

        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)

    def form_object(self) -> NoReturn:
        """Redefined method of DataN."""
        self.apply_constraints()

    def apply_constraints(self) -> NoReturn:
        """
        Symmetry constraints on parameters.

        Returns
        -------
        NoReturn
            DESCRIPTION.

        """
        space_group = self.space_group
        space_group_wyckoff = space_group.space_group_wyckoff

        cell = self.cell
        cell.type_cell = space_group.bravais_type
        cell.it_coordinate_system_code = space_group.it_coordinate_system_code
        cell.apply_constraints()

        atom_site = self.atom_site
        atom_site.apply_constraints(space_group_wyckoff)

        if self.is_attribute("atom_site_aniso"):
            atom_site_aniso = self.atom_site_aniso
            atom_site_aniso.apply_space_group_constraint(
                atom_site, space_group)

        if self.is_attribute("atom_site_susceptibility"):
            atom_site_susceptibility = self.atom_site_susceptibility
            atom_site_susceptibility.apply_chi_iso_constraint(cell)
            atom_site_susceptibility.apply_moment_iso_constraint(cell)
            atom_site_susceptibility.apply_space_group_constraint(
                atom_site, space_group)

    def calc_b_iso_beta(self):
        """
        Calculate b_iso and beta_ij based on atom_site and atom_sites.

        For each atom defined in atom_site.
        """
        a_s = self.atom_site
        try:
            a_s_a = self.atom_site_aniso
        except AttributeError:
            a_s_a = None
        l_b_iso, l_beta = [], []
        coeff = float(8.*numpy.pi**2)
        cell = self.cell
        for item_a_s in a_s.items:
            label_atom = item_a_s.label
            try:
                adp_type = item_a_s.adp_type
            except AttributeError:
                adp_type = None
            b_iso = 0.
            beta = (0., 0., 0., 0., 0., 0.)
            if adp_type == "Uiso":
                u_iso = float(item_a_s.u_iso_or_equiv)
                b_iso = float(8.*numpy.pi**2*u_iso)
            elif adp_type == "Biso":
                b_iso = float(item_a_s.b_iso_or_equiv)
            elif adp_type == "Uovl":
                # FIXME: correct it
                u_iso = float(item_a_s.u_iso_or_equiv)
                b_iso = coeff*u_iso
            elif adp_type == "Umpe":
                # FIXME: correct it
                u_iso = float(item_a_s.u_iso_or_equiv)
                b_iso = float(8.*numpy.pi**2*u_iso)
            elif adp_type == "Uani":
                item_a_s_a = a_s_a[label_atom]
                beta = item_a_s_a.calc_beta(cell)
            elif adp_type == "Bovl":
                # FIXME: correct it
                b_iso = float(item_a_s.b_iso_or_equiv)
            elif adp_type == "Bani":
                item_a_s_a = a_s_a[label_atom]
                beta = (float(item_a_s_a.b_11), float(item_a_s_a.b_22),
                        float(item_a_s_a.b_33), float(item_a_s_a.b_12),
                        float(item_a_s_a.b_13), float(item_a_s_a.b_23))
            l_b_iso.append(b_iso)
            l_beta.append(beta)
        np_b_iso = numpy.array(l_b_iso, dtype=float)
        np_beta = numpy.array(l_beta, dtype=float)
        return np_b_iso, np_beta

    def calc_f_nucl(self, index_h, index_k, index_l):
        """
        Calculate nuclear structure factor.

        Keyword Arguments
        -----------------
            index_h, index_k, index_l: 1D numpy array of Miller indexes

        Output
        ------
            f_nucl: 1D numpy array of Nuclear structure factor

        Example
        -------
            >>> import numpy as np
            >>> h, k, l = np.array([1,2],dtype=int), np.array([1,0],dtype=int),
                          np.array([1,0],dtype=int)
            >>> f_nucl = crystal.calc_f_nucl(h, k, l)

        """
        if isinstance(index_h, (float, int)):
            index_h = numpy.array([index_h], dtype=float)
            index_k = numpy.array([index_k], dtype=float)
            index_l = numpy.array([index_l], dtype=float)
        elif isinstance(index_h, list):
            index_h = numpy.array(index_h, dtype=float)
            index_k = numpy.array(index_k, dtype=float)
            index_l = numpy.array(index_l, dtype=float)

        space_group = self.space_group
        r_s_g_s = space_group.reduced_space_group_symop

        cell = self.cell
        atom_site = self.atom_site
        occupancy = numpy.array(atom_site.occupancy, dtype=float)
        x = numpy.array(atom_site.fract_x, dtype=float)
        y = numpy.array(atom_site.fract_y, dtype=float)
        z = numpy.array(atom_site.fract_z, dtype=float)

        atom_multiplicity = numpy.array(atom_site.multiplicity, dtype=int)
        scat_length_neutron = numpy.array(atom_site.scat_length_neutron,
                                          dtype=complex)

        occ_mult = occupancy*atom_multiplicity

        r_11 = r_s_g_s.numpy_r_11.astype(float)
        r_12 = r_s_g_s.numpy_r_12.astype(float)
        r_13 = r_s_g_s.numpy_r_13.astype(float)
        r_21 = r_s_g_s.numpy_r_21.astype(float)
        r_22 = r_s_g_s.numpy_r_22.astype(float)
        r_23 = r_s_g_s.numpy_r_23.astype(float)
        r_31 = r_s_g_s.numpy_r_31.astype(float)
        r_32 = r_s_g_s.numpy_r_32.astype(float)
        r_33 = r_s_g_s.numpy_r_33.astype(float)
        b_1 = r_s_g_s.numpy_b_1.astype(float)
        b_2 = r_s_g_s.numpy_b_2.astype(float)
        b_3 = r_s_g_s.numpy_b_3.astype(float)

        phase_3d = calc_phase_by_hkl_xyz_rb(
            index_h, index_k, index_l, x, y, z, r_11, r_12, r_13, r_21, r_22,
            r_23, r_31, r_32, r_33, b_1, b_2, b_3)

        b_iso, beta = self.calc_b_iso_beta()

        dwf_3d = calc_dwf(cell, index_h, index_k, index_l, b_iso, beta, r_11,
                          r_12, r_13, r_21, r_22, r_23, r_31, r_32, r_33)

        hh = phase_3d*dwf_3d
        phase_2d = hh.sum(axis=2)  # sum over symmetry

        b_scat_2d = numpy.meshgrid(index_h, scat_length_neutron,
                                   indexing="ij")[1]
        occ_mult_2d = numpy.meshgrid(index_h, occ_mult, indexing="ij")[1]

        hh = phase_2d * b_scat_2d * occ_mult_2d
        # nuclear structure factor in assymetric unit cell
        f_hkl_as = hh.sum(axis=1)*1./r_11.size

        f_nucl = space_group.calc_f_hkl_by_f_hkl_as(index_h, index_k, index_l,
                                                    f_hkl_as)
        return f_nucl

    def calc_refln(self, index_h, index_k, index_l,
                   flag_internal: bool = True):
        """
        Calculate Refln cryspy object where nuclear structure factor is stored.

        Keyword Arguments:
        -----------------
            h, k, l: 1D numpy array of Miller indexes
            flag_internal: a flag to calculate or to use internal objects.
                           It should be True if user call the function.
                           It's True by default.

        Output:
        -------
            refln: object cryspy.Refln

        Example:
        -------
            >>> import numpy as np
            >>> h, k, l = np.array([1,2],dtype=int), np.array([1,0],dtype=int),
                          np.array([1,0],dtype=int)
            >>> refln = crystal.calc_refln(h, k, l)
            >>> print(refln.to_cif())
        """
        if isinstance(index_h, (float, int)):
            index_h = numpy.array([index_h], dtype=float)
            index_k = numpy.array([index_k], dtype=float)
            index_l = numpy.array([index_l], dtype=float)
        elif isinstance(index_h, list):
            index_h = numpy.array(index_h, dtype=float)
            index_k = numpy.array(index_k, dtype=float)
            index_l = numpy.array(index_l, dtype=float)
        
        f_nucl = self.calc_f_nucl(index_h, index_k, index_l)
        res = ReflnL(loop_name=self.data_name)
        res.numpy_index_h = index_h
        res.numpy_index_k = index_k
        res.numpy_index_l = index_l
        res.numpy_f_calc = f_nucl
        if flag_internal:
            res.numpy_to_items()
        return res

    def calc_susceptibility_moment_tensor(
            self, index_h, index_k, index_l, flag_only_orbital: bool = False):
        """
        Susceptibility tensor function.

        Calculate susceptibility tensor and moment tensor in Cartesian
        orthogonal system (x||a*, z||c)

        Keyword Arguments:
        -----------------
            index_h, index_k, index_l: 1D numpy array of Miller indexes
            flag_only_orbital default is False. When only orbital form-factor
            should be used put True

        Output:
        -------
            - CHI_11, CHI_12, CHI_13, CHI_21, CHI_22, CHI_23, CHI_31, CHI_32,
              CHI_33: 1D numpy array of susceptibility tensor

            - M_11, M_12, M_13, M_21, M_22, M_23, M_31, M_32, M_33: 1D numpy
              array of moment tensor

        Example:
        -------
            >>> import numpy as np
            >>> h, k, l = np.array([1,2],dtype=int), np.array([1,0],dtype=int),
                          np.array([1,0],dtype=int)
            >>> CHI_M = crystal.calc_susceptibility_moment_tensor(h, k, l)
            >>> CHI_11, CHI_12, CHI_13, CHI_21, CHI_22, CHI_23, CHI_31, CHI_32,
                CHI_33 = CHI_M[:9]
            >>> M_11, M_12, M_13, M_21, M_22, M_23, M_31, M_32, M_33 =
                CHI_M[9:]

        """
        if isinstance(index_h, (float, int)):
            index_h = numpy.array([index_h], dtype=float)
            index_k = numpy.array([index_k], dtype=float)
            index_l = numpy.array([index_l], dtype=float)
        elif isinstance(index_h, list):
            index_h = numpy.array(index_h, dtype=float)
            index_k = numpy.array(index_k, dtype=float)
            index_l = numpy.array(index_l, dtype=float)

        space_group = self.space_group
        r_s_g_s = space_group.reduced_space_group_symop

        cell = self.cell
        sthovl = cell.calc_sthovl(index_h, index_k, index_l)

        try:
            atom_site_scat = self.atom_site_scat
            atom_site_susceptibility = self.atom_site_susceptibility
        except AttributeError:
            s_11 = numpy.zeros(index_h.shape, dtype=float)
            s_12 = numpy.zeros(index_h.shape, dtype=float)
            s_13 = numpy.zeros(index_h.shape, dtype=float)
            s_21 = numpy.zeros(index_h.shape, dtype=float)
            s_22 = numpy.zeros(index_h.shape, dtype=float)
            s_23 = numpy.zeros(index_h.shape, dtype=float)
            s_31 = numpy.zeros(index_h.shape, dtype=float)
            s_32 = numpy.zeros(index_h.shape, dtype=float)
            s_33 = numpy.zeros(index_h.shape, dtype=float)
            sm_11 = numpy.zeros(index_h.shape, dtype=float)
            sm_12 = numpy.zeros(index_h.shape, dtype=float)
            sm_13 = numpy.zeros(index_h.shape, dtype=float)
            sm_21 = numpy.zeros(index_h.shape, dtype=float)
            sm_22 = numpy.zeros(index_h.shape, dtype=float)
            sm_23 = numpy.zeros(index_h.shape, dtype=float)
            sm_31 = numpy.zeros(index_h.shape, dtype=float)
            sm_32 = numpy.zeros(index_h.shape, dtype=float)
            sm_33 = numpy.zeros(index_h.shape, dtype=float)
            return s_11, s_12, s_13, s_21, s_22, s_23, s_31, s_32, s_33,\
                sm_11, sm_12, sm_13, sm_21, sm_22, sm_23, sm_31, sm_32, sm_33

        n_item = len(atom_site_susceptibility.items)

        try:
            chi_11 = numpy.array(atom_site_susceptibility.chi_11, dtype=float)
        except AttributeError:
            chi_11 = numpy.zeros(shape=(n_item, ), dtype=float)
        try:
            chi_22 = numpy.array(atom_site_susceptibility.chi_22, dtype=float)
        except AttributeError:
            chi_22 = numpy.zeros(shape=(n_item, ), dtype=float)
        try:
            chi_33 = numpy.array(atom_site_susceptibility.chi_33, dtype=float)
        except AttributeError:
            chi_33 = numpy.zeros(shape=(n_item, ), dtype=float)
        try:
            chi_12 = numpy.array(atom_site_susceptibility.chi_12, dtype=float)
        except AttributeError:
            chi_12 = numpy.zeros(shape=(n_item, ), dtype=float)
        try:
            chi_13 = numpy.array(atom_site_susceptibility.chi_13, dtype=float)
        except AttributeError:
            chi_13 = numpy.zeros(shape=(n_item, ), dtype=float)
        try:
            chi_23 = numpy.array(atom_site_susceptibility.chi_23, dtype=float)
        except AttributeError:
            chi_23 = numpy.zeros(shape=(n_item, ), dtype=float)

        try:
            moment_11 = numpy.array(atom_site_susceptibility.moment_11,
                                    dtype=float)
        except AttributeError:
            moment_11 = numpy.zeros(shape=(n_item, ), dtype=float)
        try:
            moment_22 = numpy.array(atom_site_susceptibility.moment_22,
                                    dtype=float)
        except AttributeError:
            moment_22 = numpy.zeros(shape=(n_item, ), dtype=float)
        try:
            moment_33 = numpy.array(atom_site_susceptibility.moment_33,
                                    dtype=float)
        except AttributeError:
            moment_33 = numpy.zeros(shape=(n_item, ), dtype=float)
        try:
            moment_12 = numpy.array(atom_site_susceptibility.moment_12,
                                    dtype=float)
        except AttributeError:
            moment_12 = numpy.zeros(shape=(n_item, ), dtype=float)
        try:
            moment_13 = numpy.array(atom_site_susceptibility.moment_13,
                                    dtype=float)
        except AttributeError:
            moment_13 = numpy.zeros(shape=(n_item, ), dtype=float)
        try:
            moment_23 = numpy.array(atom_site_susceptibility.moment_23,
                                    dtype=float)
        except AttributeError:
            moment_23 = numpy.zeros(shape=(n_item, ), dtype=float)

        atom_site = self.atom_site
        atom_site_scat.load_atom_type_scat_by_atom_site(atom_site)

        np_x_y_z_occ_mult = numpy.array([(
            atom_site[item.label].fract_x, atom_site[item.label].fract_y,
            atom_site[item.label].fract_z, atom_site[item.label].occupancy,
            atom_site[item.label].multiplicity)
            for item in atom_site_susceptibility.items], dtype=float)

        x = np_x_y_z_occ_mult[:, 0]
        y = np_x_y_z_occ_mult[:, 1]
        z = np_x_y_z_occ_mult[:, 2]
        occupancy = np_x_y_z_occ_mult[:, 3]
        atom_multiplicity = np_x_y_z_occ_mult[:, 4]

        occ_mult = occupancy*atom_multiplicity

        r_11 = r_s_g_s.numpy_r_11.astype(float)
        r_12 = r_s_g_s.numpy_r_12.astype(float)
        r_13 = r_s_g_s.numpy_r_13.astype(float)
        r_21 = r_s_g_s.numpy_r_21.astype(float)
        r_22 = r_s_g_s.numpy_r_22.astype(float)
        r_23 = r_s_g_s.numpy_r_23.astype(float)
        r_31 = r_s_g_s.numpy_r_31.astype(float)
        r_32 = r_s_g_s.numpy_r_32.astype(float)
        r_33 = r_s_g_s.numpy_r_33.astype(float)
        b_1 = r_s_g_s.numpy_b_1.astype(float)
        b_2 = r_s_g_s.numpy_b_2.astype(float)
        b_3 = r_s_g_s.numpy_b_3.astype(float)

        phase_3d = calc_phase_by_hkl_xyz_rb(
            index_h, index_k, index_l, x, y, z, r_11, r_12, r_13, r_21, r_22,
            r_23, r_31, r_32, r_33, b_1, b_2, b_3)

        try:
            atom_site_aniso = self.atom_site_aniso
            flag_adp = True
        except AttributeError:
            flag_adp = False

        # FIXME
        flag_adp = False
        if flag_adp:
            b_iso = atom_site.numpy_b_iso_or_equiv
            b_iso = numpy.where(numpy.isnan(b_iso), 0, b_iso)

            beta = atom_site_aniso.calc_beta(cell)

            dwf_3d = calc_dwf(cell, index_h, index_k, index_l, b_iso, beta,
                              r_11, r_12, r_13, r_21, r_22, r_23, r_31, r_32,
                              r_33)
        else:
            dwf_3d = numpy.ones(phase_3d.shape, dtype=float)

        hh = phase_3d*dwf_3d

        # phase_2d = hh.sum(axis=2)#sum over symmetry

        # b_scat_2d = numpy.meshgrid(h, scat_length_neutron, indexing="ij")[1]
        # form_factor = atom_site_scat.calc_form_factor(
        #     sthovl, flag_only_orbital=flag_only_orbital)
        
        form_factor = [atom_site_scat[item.label].calc_form_factor(
            sthovl, flag_only_orbital=flag_only_orbital)
            for item in atom_site_susceptibility.items]
        form_factor = numpy.array(list(zip(*form_factor)), dtype=float)
        

        # dimensions: hkl, magnetic atoms, reduced symmetry operators
        ff_11, ff_12, ff_13, ff_21, ff_22, ff_23, ff_31, ff_32, ff_33 = \
            calc_form_factor_tensor_susceptibility(
                chi_11, chi_22, chi_33, chi_12, chi_13, chi_23,
                r_s_g_s, form_factor, cell, index_h, index_k, index_l)

        ffm_11, ffm_12, ffm_13, ffm_21, ffm_22, ffm_23, ffm_31, ffm_32, \
            ffm_33 = calc_form_factor_tensor_susceptibility(
                moment_11, moment_22, moment_33, moment_12, moment_13,
                moment_23, r_s_g_s, form_factor, cell, index_h, index_k,
                index_l)

        occ_mult_2d = numpy.meshgrid(index_h, occ_mult, indexing="ij")[1]

        # dimensions: hkl, number of atoms, reduced symmetry operators,
        # 18 elements of susceptribility tensor
        hh = numpy.stack([ff_11, ff_12, ff_13, ff_21, ff_22, ff_23, ff_31,
                          ff_32, ff_33, ffm_11, ffm_12, ffm_13, ffm_21, ffm_22,
                          ffm_23, ffm_31, ffm_32, ffm_33], axis=-1)

        b_scat_4d = hh

        hh = (phase_3d * dwf_3d *
              occ_mult_2d[:, :, numpy.newaxis])[:, :, :, numpy.newaxis] *\
            b_scat_4d
        # nuclear structure factor in assymetric unit cell
        f_hkl_as_2d = (hh.sum(axis=2)*1./len(r_11)).sum(axis=1)

        f_2d = space_group.calc_f_hkl_by_f_hkl_as(index_h, index_k, index_l,
                                                  f_hkl_as_2d)
        hh = 0.2695*f_2d
        # the structure factor tensor in local coordinate system (ia, ib, ic)
        # chi in 10-12 cm; chim in muB (it is why here 0.2695)
        s_11, s_12, s_13, s_21, s_22, s_23, s_31, s_32, s_33 = \
            self._orto_matrix((hh[:, 0], hh[:, 1], hh[:, 2], hh[:, 3],
                               hh[:, 4], hh[:, 5], hh[:, 6], hh[:, 7],
                               hh[:, 8]))

        sm_11, sm_12, sm_13, sm_21, sm_22, sm_23, sm_31, sm_32, sm_33 = \
            self._orto_matrix((hh[:, 9], hh[:, 10], hh[:, 11], hh[:, 12],
                               hh[:, 13], hh[:, 14], hh[:, 15], hh[:, 16],
                               hh[:, 17]))

        return s_11, s_12, s_13, s_21, s_22, s_23, s_31, s_32, s_33, \
            sm_11, sm_12, sm_13, sm_21, sm_22, sm_23, sm_31, sm_32, sm_33

    def calc_refln_susceptibility(
            self, index_h, index_k, index_l, flag_internal: bool = True,
            flag_only_orbital: bool = False):
        """
        Calculate susceptibility tensor and moment tensor.

        They are given in Cartesian orthogonal system (x||a*, z||c).

        Keyword Arguments:
        -----------------
            h, k, l: 1D numpy array of Miller indexes

        Output:
        -------
            ReflnSusceptibilityL object of cryspy

        Example:
        -------
            >>> import numpy as np
            >>> h, k, l = np.array([1,2],dtype=int), np.array([1,0],dtype=int),
                          np.array([1,0],dtype=int)
            >>> refln_suscept = crystal.calc_refln_susceptibility(h, k, l)
        """
        if isinstance(index_h, (float, int)):
            index_h = numpy.array([index_h], dtype=float)
            index_k = numpy.array([index_k], dtype=float)
            index_l = numpy.array([index_l], dtype=float)
        elif isinstance(index_h, list):
            index_h = numpy.array(index_h, dtype=float)
            index_k = numpy.array(index_k, dtype=float)
            index_l = numpy.array(index_l, dtype=float)

        
        CHI_M = self.calc_susceptibility_moment_tensor(
            index_h, index_k, index_l, flag_only_orbital=flag_only_orbital)

        s_11, s_12, s_13, s_21, s_22, s_23, s_31, s_32, s_33 = CHI_M[:9]
        sm_11, sm_12, sm_13, sm_21, sm_22, sm_23, sm_31, sm_32, sm_33 = \
            CHI_M[9:]

        res = ReflnSusceptibilityL(loop_name=self.data_name)
        res.numpy_index_h = index_h
        res.numpy_index_k = index_k
        res.numpy_index_l = index_l
        res.numpy_chi_11_calc = s_11
        res.numpy_chi_12_calc = s_12
        res.numpy_chi_13_calc = s_13
        res.numpy_chi_21_calc = s_21
        res.numpy_chi_22_calc = s_22
        res.numpy_chi_23_calc = s_23
        res.numpy_chi_31_calc = s_31
        res.numpy_chi_32_calc = s_32
        res.numpy_chi_33_calc = s_33
        res.numpy_moment_11_calc = sm_11
        res.numpy_moment_12_calc = sm_12
        res.numpy_moment_13_calc = sm_13
        res.numpy_moment_21_calc = sm_21
        res.numpy_moment_22_calc = sm_22
        res.numpy_moment_23_calc = sm_23
        res.numpy_moment_31_calc = sm_31
        res.numpy_moment_32_calc = sm_32
        res.numpy_moment_33_calc = sm_33
        if flag_internal:
            res.numpy_to_items()
        return res

    def _orto_matrix(self, l_ij):
        """
        Matrix ortogonalization.

        matrix l_ij is defined in coordinate system (a, b, c)
        l_ij: = l_11, l_12, l_13, l_21, l_22, l_23, l_31, l_32, l_33
        output matrix s_ij is defined in Cartezian coordinate system defined
        as x||a*, z||c, y= [z x] (right handed)
        """
        cell = self.cell
        s_11, s_12, s_13, s_21, s_22, s_23, s_31, s_32, s_33 = \
            cell.ortogonalize_matrix(l_ij)
        return s_11, s_12, s_13, s_21, s_22, s_23, s_31, s_32, s_33

    def calc_hkl(self, sthol_min: float = 0., sthovl_max: float = 1.):
        """
        Calculate hkl and multiplicity taking into account the symmetry
        constraints.

        Parameters
        ----------
        sthol_min : float
            minimal sin(theta)/wavelength in inversed angstrems.
        sthovl_max : float
            maximal sin(theta)/wavelength in inversed angstrems.

        Returns
        -------
        h, k, l, mult : numpy.array[float]
            The Miller indeces: h, k, l, and its multiplicity: mult.

        """
        cell = self.cell
        space_group = self.space_group
        res = cell.calc_hkl(space_group, sthol_min, sthovl_max)
        return res

    def calc_hkl_in_range(self, sthol_min: float = 0., sthovl_max: float = 1.):
        """
        Give hkl without taking symmetry into account.

        Parameters
        ----------
        sthol_min : float
            DESCRIPTION.
        sthovl_max : float
            DESCRIPTION.

        Returns
        -------
        res : TYPE
            DESCRIPTION.

        """
        cell = self.cell
        res = cell.calc_hkl_in_range(sthol_min, sthovl_max)
        return res

    def calc_magnetization_ellipsoid(self):
        """
        Magnetization ellipsoids.

        The same coordinate system as U_ij (anisotropic Debye-Waller factor)

        Negtive eigenvalues of ellipsoid are replaced by positive.
        """
        try:
            cell = self.cell
            a_s_m_a = self.atom_site_susceptibility
        except AttributeError:
            return l_res

        l_res = []
        for item in a_s_m_a.items:
            chi_as_u_loc = magnetization_ellipsoid_by_u_ij(cell, item)
            l_res.append(chi_as_u_loc)
        return l_res

    def calc_main_axes_of_magnetization_ellipsoids(self):
        """Susceptibility along the main axes of magnetization ellipsoid.

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
        ll_moments = []
        ll_directions = []
        try:
            cell = self.cell
            a_s_s = self.atom_site_susceptibility
        except AttributeError:
            return ll_moments, ll_directions

        l_moments, l_moments_sigma, l_rot_matrix = \
            a_s_s.calc_main_axes_of_magnetization_ellipsoid(cell)
        return l_moments, l_moments_sigma, l_rot_matrix

    def report_main_axes_of_magnetization_ellipsoids(self):
        """
        Report about main axes of magnetization ellipsoids.

        Make a report about magnetization ellipsoids in string format.
        Calculations are performed by
        get_main_axes_of_magnetization_ellipsoids method
        and calc_magnetization_ellipsoid method of the crystal object.
        """
        # crystal is defined object of cryspy library;
        # type(crystal) is Crystal
        if self.is_attribute("atom_site_susceptibility"):
            a_s_m_a = self.atom_site_susceptibility
        else:
            return ""

        cell = self.cell
        s_out = report_main_axes_of_magnetization_ellipsoids(
            cell, a_s_m_a)

        return s_out

    def calc_magnetic_moments_with_field_loc(self, field_abc):
        """Orientation of magetic moment for atoms at applied magnetic field.

        The input magnetic field should be given in normalized unit cell
        (a/|a|, b/|b|, c/|c|)
        The output magnetic moment are given in normalized unit cell
        (a/|a|, b/|b|, c/|c|)
        
        Output
        ------
        l_lab_out - label of atoms
        l_xyz_out - position of atoms
        l_moment_out - moment of atoms

        """
        np_field = numpy.array(field_abc, dtype=float)
        l_lab_out, l_xyz_out, l_moment_out = [], [], []

        try:
            spgr = self.space_group
            cell = self.cell
            a_s = self.atom_site
            a_s_m_a = self.atom_site_susceptibility
        except AttributeError:
            return l_lab_out, l_xyz_out, l_moment_out

        m_m_norm = cell.m_m_norm
        m_mt_norm_m_norm_field = numpy.matmul(
            numpy.matmul(m_m_norm.transpose(), m_m_norm), np_field)

        for _l, _11, _22, _33, _12, _13, _23 in zip(
                a_s_m_a.label, a_s_m_a.chi_11, a_s_m_a.chi_22, a_s_m_a.chi_33,
                a_s_m_a.chi_12, a_s_m_a.chi_13, a_s_m_a.chi_23):
            m_chi = numpy.array([[_11, _12, _13],
                                 [_12, _22, _23],
                                 [_13, _23, _33]], dtype=float)
            x, y, z = a_s[_l].fract_x, a_s[_l].fract_y, a_s[_l].fract_z
            l_out = spgr.calc_rotated_matrix_for_position(m_chi, x, y, z)
            for _i_out, _out in enumerate(l_out):
                _xyz = _out[0]
                _chi = _out[1]
                _moment = numpy.matmul(_chi, m_mt_norm_m_norm_field)

                l_lab_out.append(f"{_l:}_{_i_out+1:}")
                l_xyz_out.append(_xyz)
                l_moment_out.append(_moment)
        return l_lab_out, l_xyz_out, l_moment_out

    def report(self):
        return self.report_main_axes_of_magnetization_ellipsoids()

    def plots(self):
        return []

# s_cont = """
#   data_Fe3O4
#   _cell_angle_alpha 90.0
#   _cell_angle_beta 90.0
#   _cell_angle_gamma 90.0
#   _cell_length_a 8.56212()
#   _cell_length_b 8.56212
#   _cell_length_c 8.56212
#   _space_group_it_coordinate_system_code 2
#   _space_group_IT_number    227

#   loop_
#   _atom_site_adp_type
#   _atom_site_B_iso_or_equiv
#   _atom_site_fract_x
#   _atom_site_fract_y
#   _atom_site_fract_z
#   _atom_site_label
#   _atom_site_occupancy
#   _atom_site_type_symbol
#   Uani 0.0 0.125 0.125 0.125 Fe3A 1.0 Fe3+
#   Uani 0.0 0.5 0.5 0.5 Fe3B 1.0 Fe3+
#   Uiso 0.0 0.25521 0.25521 0.25521 O1 1.0 O2-

#   loop_
#   _atom_type_scat_length_neutron
#   _atom_type_symbol
#     0.945 Fe3+
#   0.5803 O2-

#   loop_
#   _atom_site_aniso_U_11
#   _atom_site_aniso_U_12
#   _atom_site_aniso_U_13
#   _atom_site_aniso_U_22
#   _atom_site_aniso_U_23
#   _atom_site_aniso_U_33
#   _atom_site_aniso_label
#   0.0 0.0 0.0 0.0 0.0 0.0 Fe3A
#   0.0 0.0 0.0 0.0 0.0 0.0 Fe3B

#   loop_
#   _atom_site_scat_label
#   _atom_site_scat_lande
#   Fe3A 2.0
#   Fe3B 2.0

#   loop_
#   _atom_site_susceptibility_label
#   _atom_site_susceptibility_chi_type
#   _atom_site_susceptibility_chi_11
#   _atom_site_susceptibility_chi_12
#   _atom_site_susceptibility_chi_13
#   _atom_site_susceptibility_chi_22
#   _atom_site_susceptibility_chi_23
#   _atom_site_susceptibility_chi_33
#   Fe3A Cani -3.468(74) 0.0 0.0 -3.468 0.0 -3.468
#   Fe3B Cani 3.041      0.0 0.0  3.041(21) 0.0  3.041

# """

# obj = Crystal.from_cif(s_cont)
# for var_name in obj.get_variable_names():
#     print(var_name)
#     obj.set_variable_by_name(var_name, 7)
# print(obj.atom_site_susceptibility)
# obj.apply_constraints()
# print(obj.atom_site_susceptibility)
# print(obj.report_main_axes_of_magnetization_ellipsoids())
