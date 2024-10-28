"""Description of Crystal class."""
__author__ = 'ikibalin'
__version__ = "2020_08_19"
import numpy

from typing import NoReturn

from cryspy.A_functions_base.database import DATABASE
from cryspy.A_functions_base.charge_form_factor import calc_jl_for_ion, calc_scattering_amplitude_tabulated, get_atom_name_ion_charge_shell, get_atomic_symbol_ion_charge_isotope_number_by_ion_symbol
from cryspy.A_functions_base.debye_waller_factor import calc_param_iso_aniso_by_b_iso_beta, calc_u_ij_by_beta
from cryspy.A_functions_base.matrix_operations import calc_m1_m2_inv_m1, calc_m_v
from cryspy.A_functions_base.magnetic_form_factor import get_j0_j2_parameters
from cryspy.A_functions_base.unit_cell import calc_eq_ccs_by_unit_cell_parameters, calc_m_m_by_unit_cell_parameters, calc_reciprocal_by_unit_cell_parameters
from cryspy.A_functions_base.structure_factor import calc_f_nucl_by_dictionary, calc_sft_ccs_by_dictionary, calc_f_m_perp_by_sft, calc_bulk_susceptibility_by_dictionary
from cryspy.A_functions_base.symmetry_elements import calc_full_mag_elems, calc_symm_flags, define_bravais_type_by_symm_elems
from cryspy.A_functions_base.symmetry_constraints import calc_sc_beta, calc_sc_fract_sc_b, calc_sc_chi, calc_sc_chi_full

from cryspy.B_parent_classes.cl_3_data import DataN
from cryspy.B_parent_classes.preocedures import take_items_by_class

from cryspy.C_item_loop_classes.cl_1_atom_site import AtomSite, AtomSiteL
from cryspy.C_item_loop_classes.cl_1_atom_type import AtomTypeL
from cryspy.C_item_loop_classes.cl_1_atom_site_aniso import \
    AtomSiteAnisoL
from cryspy.C_item_loop_classes.cl_1_atom_site_susceptibility import \
    AtomSiteSusceptibilityL
from cryspy.C_item_loop_classes.cl_1_atom_site_exchange import AtomSiteExchangeL
from cryspy.C_item_loop_classes.cl_1_atom_type_scat import \
    AtomTypeScatL
from cryspy.C_item_loop_classes.cl_2_atom_site_scat import \
    AtomSiteScatL
from cryspy.C_item_loop_classes.cl_1_atom_site_moment import AtomSiteMoment, AtomSiteMomentL
from cryspy.C_item_loop_classes.cl_1_refln_susceptibility import \
    ReflnSusceptibilityL
from cryspy.C_item_loop_classes.cl_1_atom_local_axes import \
    AtomLocalAxesL
from cryspy.C_item_loop_classes.cl_1_atom_electron_configuration \
    import AtomElectronConfigurationL
from cryspy.C_item_loop_classes.cl_1_cell import Cell
from cryspy.C_item_loop_classes.cl_1_refln import ReflnL
from cryspy.C_item_loop_classes.cl_1_space_group_symop_magn_centering import \
    SpaceGroupSymopMagnCenteringL

from cryspy.C_item_loop_classes.cl_2_atom_rho_orbital_radial_slater \
    import AtomRhoOrbitalRadialSlaterL
from cryspy.C_item_loop_classes.cl_2_space_group_symop_magn_operation import \
    SpaceGroupSymopMagnOperationL
from cryspy.C_item_loop_classes.cl_2_space_group import SpaceGroup





from cryspy.D_functions_item_loop.function_1_report_magnetization_ellipsoid \
    import magnetization_ellipsoid_by_u_ij, report_main_axes_of_magnetization_ellipsoids

na = numpy.newaxis

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

    CLASSES_MANDATORY = (Cell, AtomSiteL,)
    CLASSES_OPTIONAL = (SpaceGroup,
        AtomTypeL, AtomSiteAnisoL, AtomSiteSusceptibilityL, AtomSiteExchangeL, AtomSiteScatL,
        AtomTypeScatL, AtomLocalAxesL, AtomElectronConfigurationL,
        AtomSiteMomentL, SpaceGroupSymopMagnOperationL, SpaceGroupSymopMagnCenteringL)
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
        if self.is_attribute("space_group"):
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
                atom_site_susceptibility.apply_space_group_constraint(
                    atom_site, space_group, cell)
            if self.is_attribute("atom_site_exchange"):
                atom_site_exchange = self.atom_site_exchange
                atom_site_exchange.apply_j_iso_constraint()
                atom_site_exchange.apply_space_group_constraint(
                    atom_site, space_group, cell)

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


    def calc_f_nucl(self, index_hkl):
        dict_crystal = self.get_dictionary()
        dict_in_out = {"index_hkl": index_hkl}
        f_nucl, dder = calc_f_nucl_by_dictionary(dict_crystal, dict_in_out, flag_use_precalculated_data=False)
        return f_nucl


    def calc_refln(self, index_hkl,
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
        # if isinstance(index_h, (float, int)):
        #     index_h = numpy.array([index_h], dtype=float)
        #     index_k = numpy.array([index_k], dtype=float)
        #     index_l = numpy.array([index_l], dtype=float)
        # elif isinstance(index_h, list):
        #     index_h = numpy.array(index_h, dtype=float)
        #     index_k = numpy.array(index_k, dtype=float)
        #     index_l = numpy.array(index_l, dtype=float)

        # index_hkl = numpy.stack([index_h, index_k, index_l], axis=0)
        f_nucl = self.calc_f_nucl(index_hkl)
        res = ReflnL(loop_name=self.data_name)
        res.numpy_index_h = index_hkl[0]
        res.numpy_index_k = index_hkl[1]
        res.numpy_index_l = index_hkl[2]
        res.numpy_f_calc = f_nucl
        if flag_internal:
            res.numpy_to_items()
        return res

    def calc_structure_factor_tensor_ccs(
            self, index_hkl, flag_only_orbital: bool = False, dict_in_out: dict = None):
        """Calculate structure factor tensor in CCS coordinate (X||a*, Z||c)."""
        
        dict_crystal = self.get_dictionary()
        if dict_in_out is None:
            dict_in_out = {"index_hkl": index_hkl, "flag_only_orbital": flag_only_orbital}
        else:
            dict_in_out["index_hkl"] = index_hkl
            dict_in_out["flag_only_orbital"] = flag_only_orbital

        sft_ccs, dder = calc_sft_ccs_by_dictionary(dict_crystal, dict_in_out, flag_use_precalculated_data=False)
        return sft_ccs, dder


    def calc_refln_susceptibility(
            self, index_hkl, flag_internal: bool = True,
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
        
        sft_ccs, dder = self.calc_structure_factor_tensor_ccs(
            index_hkl, flag_only_orbital=flag_only_orbital)

        res = ReflnSusceptibilityL(loop_name=self.data_name)
        res.numpy_index_h = index_hkl[0,:]
        res.numpy_index_k = index_hkl[1,:]
        res.numpy_index_l = index_hkl[2,:]
        res.numpy_chi_11_calc = sft_ccs[0, :]
        res.numpy_chi_12_calc = sft_ccs[1, :]
        res.numpy_chi_13_calc = sft_ccs[2, :]
        res.numpy_chi_21_calc = sft_ccs[3, :]
        res.numpy_chi_22_calc = sft_ccs[4, :]
        res.numpy_chi_23_calc = sft_ccs[5, :]
        res.numpy_chi_31_calc = sft_ccs[6, :]
        res.numpy_chi_32_calc = sft_ccs[7, :]
        res.numpy_chi_33_calc = sft_ccs[8, :]
        if flag_internal:
            res.numpy_to_items()
        return res

    def calc_f_m_perp(self, index_hkl, magnetic_field, dict_in_out: dict = None):
        """Calculate F_M_perpendicular
        """
        sft_ccs, dder_sft_ccs = self.calc_structure_factor_tensor_ccs(index_hkl, dict_in_out=dict_in_out)
        unit_cell_parameters = self.cell.get_unit_cell_parameters()
        eq_ccs, dder_eq_ccs = calc_eq_ccs_by_unit_cell_parameters(index_hkl, unit_cell_parameters)
        f_m_perp, dder_f_m_perp = calc_f_m_perp_by_sft(
            sft_ccs, magnetic_field, eq_ccs,
            flag_sft_ccs=False,
            flag_magnetic_field=False,
            flag_eq_ccs=False)
        return f_m_perp


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
        norm_field = numpy.sqrt(numpy.square(np_field).sum())
        if norm_field == 0.: 
            norm_field = 1.
        np_field = np_field/norm_field
        l_lab_out, l_xyz_out, l_moment_out = [], [], []

        try:
            spgr = self.space_group
            cell = self.cell
            a_s = self.atom_site
            a_s_m_a = self.atom_site_susceptibility
        except AttributeError:
            return None


        m_m_norm = cell.m_m_norm
        m_mt_norm_m_norm_field = numpy.matmul(
            numpy.matmul(m_m_norm.transpose(), m_m_norm), np_field)

        as_p1 = []
        asm_p1 = []
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
                label_p1 = f"{_l:}_{_i_out+1:}"
                l_lab_out.append(label_p1)
                l_xyz_out.append(_xyz)
                l_moment_out.append(_moment)
                as_p1.append(AtomSite(label=label_p1, fract_x=_xyz[0], fract_y=_xyz[1], fract_z=_xyz[2]))
                asm_p1.append(AtomSiteMoment(
                    label=label_p1,
                    crystalaxis_x=numpy.round(_moment[0], 5),
                    crystalaxis_y=numpy.round(_moment[1], 5),
                    crystalaxis_z=numpy.round(_moment[2], 5)))


        spgr_p1 = SpaceGroup(it_number=1)
        spgr_p1.form_object()
        atom_site_p1 = AtomSiteL()
        atom_site_p1.items = as_p1
        atom_site_moment_p1 = AtomSiteMomentL()
        atom_site_moment_p1.items = asm_p1

        cryspy_p1 = Crystal(space_group = spgr_p1, cell=cell, atom_site=atom_site_p1, atom_site_moment=atom_site_moment_p1)

        return cryspy_p1


    def report(self):
        return self.report_main_axes_of_magnetization_ellipsoids()


    def plots(self):
        return []


    def get_flags_atom_beta(self):
        atom_site = self.atom_site
        if self.is_attribute("atom_site_aniso"):
            atom_site_aniso = self.atom_site_aniso
            l_f = []
            for item_as in atom_site.items:
                if item_as.adp_type == "Bani":
                    item_asa = atom_site_aniso[item_as.label]
                    l_f.append(
                        numpy.array([
                            item_asa.b_11_refinement, item_asa.b_22_refinement, item_asa.b_33_refinement,
                            item_asa.b_12_refinement, item_asa.b_13_refinement, item_asa.b_23_refinement],
                            dtype=bool))
                elif item_as.adp_type == "Uani":
                    item_asa = atom_site_aniso[item_as.label]
                    l_f.append(
                        numpy.array([
                            item_asa.u_11_refinement, item_asa.u_22_refinement, item_asa.u_33_refinement,
                            item_asa.u_12_refinement, item_asa.u_13_refinement, item_asa.u_23_refinement],
                            dtype=bool))
                else:
                    l_f.append(numpy.zeros((6,), dtype=bool))
            res = numpy.stack(l_f, axis=0).transpose()
        else:
            res = numpy.zeros((6, len(atom_site.items)), dtype=bool)
        return res

    def get_dictionary(self):
        """Form dictionary. See documentation moduel CrysPy using Jupyter notebook.
        """
        self.form_object()
        ddict = {}
        space_group, cell, atom_site = None, None, None
        atom_site_susceptibility, atom_electron_configuration = None, None
        atom_site_exchange = None
        atom_type_scat = None
        atom_type = None
        atom_site_aniso = None
        atom_site_scat = None
        atom_site_moment = None
        space_group_symop_magn_centering = None
        space_group_symop_magn_operation = None

        l_obj = take_items_by_class(self, (SpaceGroup,))
        if len(l_obj) > 0:
            space_group = l_obj[0]

        l_obj = take_items_by_class(self, (Cell,))
        if len(l_obj) > 0:
            cell = l_obj[0]

        l_obj = take_items_by_class(self, (AtomSiteL,))
        if len(l_obj) > 0:
            atom_site = l_obj[0]

        l_obj = take_items_by_class(self, (AtomSiteAnisoL,))
        if len(l_obj) > 0:
            atom_site_aniso = l_obj[0]

        l_obj = take_items_by_class(self, (AtomSiteScatL,))
        if len(l_obj) > 0:
            atom_site_scat = l_obj[0]

        l_obj = take_items_by_class(self, (AtomTypeL,))
        if len(l_obj) > 0:
            atom_type = l_obj[0]

        l_obj = take_items_by_class(self, (AtomTypeScatL,))
        if len(l_obj) > 0:
            atom_type_scat = l_obj[0]

        l_obj = take_items_by_class(self, (AtomSiteSusceptibilityL,))
        if len(l_obj) > 0:
            atom_site_susceptibility = l_obj[0]

        l_obj = take_items_by_class(self, (AtomSiteExchangeL,))
        if len(l_obj) > 0:
            atom_site_exchange = l_obj[0]
            

        l_obj = take_items_by_class(self, (AtomElectronConfigurationL,))
        if len(l_obj) > 0:
            atom_electron_configuration = l_obj[0]

        l_obj = take_items_by_class(self, (AtomSiteMomentL,))
        if len(l_obj) > 0:
            atom_site_moment = l_obj[0]

        l_obj = take_items_by_class(self, (SpaceGroupSymopMagnCenteringL,))
        if len(l_obj) > 0:
            space_group_symop_magn_centering = l_obj[0]

        l_obj = take_items_by_class(self, (SpaceGroupSymopMagnOperationL,))
        if len(l_obj) > 0:
            space_group_symop_magn_operation = l_obj[0]

        ddict["name"] = self.data_name
        ddict["type_name"] = self.get_name()
        ind_mag = None
        rad = numpy.pi/180.
        if space_group is not None:
            full_symm_elems = space_group.full_space_group_symop.get_symm_elems()
            ddict["full_symm_elems"] = full_symm_elems
            reduced_symm_elems = space_group.reduced_space_group_symop.get_symm_elems()
            ddict["reduced_symm_elems"] = reduced_symm_elems
            ddict["centrosymmetry"] = space_group.centrosymmetry
            if ddict["centrosymmetry"]:
                p_centr = space_group.pcentr
                lcm = numpy.lcm.reduce([fr.denominator for fr in p_centr])
                vals = [int(fr.numerator * lcm / fr.denominator) for fr in p_centr]
                vals.append(lcm)
                ddict["centrosymmetry_position"] = numpy.array(vals, dtype=int)

            shift = space_group.shift
            l_vals = []
            for ind in range(len(shift)):
                lcm = numpy.lcm.reduce([fr.denominator for fr in shift[ind]])
                vals = [int(fr.numerator * lcm / fr.denominator) for fr in shift[ind]]
                vals.append(lcm)
                l_vals.append(vals)
            ddict["translation_elems"] = numpy.array(l_vals, dtype=int).transpose()

        if space_group_symop_magn_operation is not None:
            if space_group_symop_magn_centering is not None:
                se_sgsmo = self.space_group_symop_magn_operation.get_sym_elems()
                se_sgsmc = self.space_group_symop_magn_centering.get_sym_elems()
                full_mcif_elems = calc_full_mag_elems(se_sgsmo, se_sgsmc)
                ddict["full_mcif_elems"] = full_mcif_elems

        if cell is not None:
            if cell.is_attribute("type_cell"):
                type_cell, it_coordinate_system_code = cell.type_cell, cell.it_coordinate_system_code
            elif "full_mcif_elems" in ddict.keys(): 
                type_cell, it_coordinate_system_code = define_bravais_type_by_symm_elems(full_mcif_elems)
                cell.type_cell = type_cell
                cell.it_coordinate_system_code = it_coordinate_system_code
                cell.apply_constraints()

            unit_cell_parameters = cell.get_unit_cell_parameters()
            
            sc_uc, v_uc = calc_sc_v_unit_cell_parameters(type_cell, it_coordinate_system_code)
            ddict["sc_uc"] = sc_uc
            ddict["v_uc"] = v_uc
            
            ddict["unit_cell_parameters"] = numpy.dot(sc_uc, unit_cell_parameters)+v_uc
            ddict["flags_unit_cell_parameters"] = cell.get_flags_unit_cell_parameters()

        if atom_site is not None:
            atom_label = numpy.array(atom_site.label, dtype=str)
            atom_occupancy = numpy.array(atom_site.occupancy, dtype=float)
            if atom_site.is_attribute("multiplicity"):
                atom_multiplicity = numpy.array(atom_site.multiplicity, dtype=int)
                ddict["atom_multiplicity"] = atom_multiplicity
            atom_fract_xyz = numpy.array([atom_site.fract_x,
                                          atom_site.fract_y,
                                          atom_site.fract_z], dtype=float)
            atom_fract_xyz = numpy.mod(atom_fract_xyz, 1.)

            atom_type_symbol = numpy.array(atom_site.type_symbol, dtype=str)
            table_sthovl = numpy.linspace(0, 2, 501) # fro0 to 2 Å**-1
            d_dispersion = DATABASE["Dispersion"]
            table_wavelength = d_dispersion["table_wavelength"]
            l_table_atom_scattering_amplitude = []
            l_table_atom_dispersion = []
            flag_atom_scattering_amplitude = True
            for type_symbol in atom_type_symbol:
                try:
                    jl = calc_jl_for_ion(table_sthovl, type_symbol)
                    scattering_amplitude = jl[:,0]
                except UserWarning:
                    try:
                        scattering_amplitude = calc_scattering_amplitude_tabulated(type_symbol, table_sthovl)[0]
                    except KeyError:
                        flag_atom_scattering_amplitude = False
                if flag_atom_scattering_amplitude:
                    l_table_atom_scattering_amplitude.append(scattering_amplitude)

                type_atom = get_atom_name_ion_charge_shell(type_symbol)[0]
                try:
                    l_table_atom_dispersion.append(d_dispersion[type_atom.capitalize()])
                except KeyError:
                    l_table_atom_dispersion.append(numpy.zeros_like(table_wavelength))
            if flag_atom_scattering_amplitude:
                ddict["table_sthovl"] = table_sthovl
                ddict["table_atom_scattering_amplitude"] = numpy.stack(l_table_atom_scattering_amplitude, axis=0)
            ddict["table_wavelength"] = table_wavelength
            ddict["table_atom_dispersion"] = numpy.stack(l_table_atom_dispersion, axis=0)

            atom_b_scat = numpy.array(atom_site.scat_length_neutron, dtype=complex)

            b_d  = 3*10**6*numpy.ones(atom_fract_xyz.shape[1:], dtype=int) #10**6 it is precision
            b_elem = (atom_fract_xyz*b_d[na, :]).astype(int)
            gcd = numpy.gcd.reduce(numpy.gcd(b_elem, b_d[na,:]))

            atom_symm_elems = numpy.concatenate([
                numpy.floor_divide(b_elem, gcd[na, :]),
                numpy.floor_divide(b_d, gcd)[na, :],
                numpy.zeros((9, ) + atom_fract_xyz.shape[1:], dtype=int)],
                axis=0)


            if "full_mcif_elems" in ddict.keys():
                elems_fs = ddict["full_mcif_elems"]
            elif "full_symm_elems" in ddict.keys():
                elems_fs = ddict["full_symm_elems"]
            else:
                raise AttributeError("Symmetry elements are not found.")

            flag_2d = calc_symm_flags(elems_fs[:, :, na], atom_symm_elems[:, na, :])


            l_sc_fract, l_sc_b = [], []
            for ind in range(atom_label.size):
                symm_elems = elems_fs[:, flag_2d[:, ind]]
                sc_fract, sc_b = calc_sc_fract_sc_b(symm_elems, atom_fract_xyz[:, ind])
                l_sc_fract.append(sc_fract)
                l_sc_b.append(sc_b)
            atom_site_sc_fract = numpy.stack(l_sc_fract, axis=-1)

            ddict["atom_site_sc_fract"] = atom_site_sc_fract 
            ddict["atom_site_sc_b"] = numpy.stack(l_sc_b, axis=-1)

            flags_atom_fract_xyz = atom_site.get_flags_atom_fract_xyz()
            flags_atom_fract_xyz_constr = calc_m_v(atom_site_sc_fract, flags_atom_fract_xyz, flag_m=False, flag_v=False)[0].astype(bool)

            ddict["atom_label"] = atom_label
            ddict["atom_fract_xyz"] = atom_fract_xyz
            ddict["flags_atom_fract_xyz"] = flags_atom_fract_xyz_constr 
            ddict["atom_scat_length_neutron"] = atom_b_scat
            ddict["atom_occupancy"] = atom_occupancy
            ddict["flags_atom_occupancy"] = atom_site.get_flags_atom_occupancy()
            ddict["flags_atom_b_iso"] = atom_site.get_flags_atom_b_iso()
            ddict["atom_symm_elems_flag"] = flag_2d

        if atom_site_susceptibility is not None:
            atom_para_label = numpy.array(atom_site_susceptibility.label, dtype=str)
            atom_para_chi_type = numpy.array(atom_site_susceptibility.chi_type, dtype=str)
            atom_para_susceptibility = numpy.array([
                atom_site_susceptibility.chi_11, atom_site_susceptibility.chi_22,
                atom_site_susceptibility.chi_33, atom_site_susceptibility.chi_12,
                atom_site_susceptibility.chi_13, atom_site_susceptibility.chi_23],
                dtype=float)
            ddict["atom_para_label"] = atom_para_label
            ddict["atom_para_susceptibility"] = atom_para_susceptibility
            ddict["flags_atom_para_susceptibility"] = atom_site_susceptibility.get_flags_susceptibility()            
            if "atom_label" in ddict.keys():
                flag_2d = numpy.expand_dims(ddict["atom_label"], axis=1) == numpy.expand_dims(atom_para_label, axis=0)
                args = numpy.argwhere(flag_2d)
                atom_para_index = args[args[:,1].argsort()][:,0]
                ddict["atom_para_index"] = atom_para_index

                if "full_mcif_elems" in ddict.keys():
                    elems_fs = ddict["full_mcif_elems"]
                elif "full_symm_elems" in ddict.keys():
                    elems_fs = ddict["full_symm_elems"]
                else:
                    raise AttributeError("Symmetry elements are not found.")
                flag_2d = ddict["atom_symm_elems_flag"][:, atom_para_index]

                l_sc_chi = []
                for ind in range(atom_para_label.size):
                    if (atom_para_chi_type[ind]).lower() == "ciso":
                        sc_chi = numpy.array([
                            [1., 0., 0., 0., 0., 0.],
                            [1., 0., 0., 0., 0., 0.],
                            [1., 0., 0., 0., 0., 0.],
                            [0., 0., 0., 0., 0., 0.],
                            [0., 0., 0., 0., 0., 0.],
                            [0., 0., 0., 0., 0., 0.]], dtype=float)
                    else:
                        sc_chi = calc_sc_chi(
                            elems_fs[:, flag_2d[:, ind]],
                            unit_cell_parameters=unit_cell_parameters, flag_unit_cell_parameters=False)[0]
                    l_sc_chi.append(sc_chi)
                atom_para_sc_chi = numpy.stack(l_sc_chi, axis=-1)
                ddict["atom_para_sc_chi"] = atom_para_sc_chi
        
        if atom_site_exchange is not None:
            atom_exchange_label = numpy.array(atom_site_exchange.label, dtype=str)
            atom_exchange_type = numpy.array(atom_site_exchange.j_type, dtype=str)
            atom_exchange_tensor = numpy.array([
                atom_site_exchange.j_11, 
                atom_site_exchange.j_22, 
                atom_site_exchange.j_33, 
                atom_site_exchange.j_12, 
                atom_site_exchange.j_13, 
                atom_site_exchange.j_23, 
                ], dtype=float)
            flags_atom_exchange_tensor = numpy.array([
                atom_site_exchange.j_11_refinement,
                atom_site_exchange.j_22_refinement,
                atom_site_exchange.j_33_refinement,
                atom_site_exchange.j_12_refinement,
                atom_site_exchange.j_13_refinement,
                atom_site_exchange.j_23_refinement,
                ], dtype=bool)
            ddict["atom_exchange_label"] = atom_exchange_label
            ddict["atom_exchange_type"] = atom_exchange_type
            ddict["atom_exchange_tensor"] = atom_exchange_tensor
            ddict["flags_atom_exchange_tensor"] = flags_atom_exchange_tensor
            if "atom_label" in ddict.keys():

                if "full_mcif_elems" in ddict.keys():
                    elems_fs = ddict["full_mcif_elems"]
                elif "full_symm_elems" in ddict.keys():
                    elems_fs = ddict["full_symm_elems"]
                else:
                    raise AttributeError("Symmetry elements are not found.")

                flag_2d = numpy.expand_dims(ddict["atom_label"], axis=1) == numpy.expand_dims(atom_exchange_label, axis=0)
                args = numpy.argwhere(flag_2d)
                atom_exchange_index = args[args[:,1].argsort()][:,0]
                ddict["atom_exchange_index"] = atom_exchange_index

                l_sc_exchange = []
                flag_2d = ddict["atom_symm_elems_flag"][:, atom_exchange_index]
                for ind in range(atom_exchange_label.size):
                    if (atom_exchange_type[ind]).lower() == "jiso":
                        sc_exchange = numpy.array([
                            [1., 0., 0., 0., 0., 0.],
                            [1., 0., 0., 0., 0., 0.],
                            [1., 0., 0., 0., 0., 0.],
                            [0., 0., 0., 0., 0., 0.],
                            [0., 0., 0., 0., 0., 0.],
                            [0., 0., 0., 0., 0., 0.]], dtype=float)
                    else:
                        sc_exchange = calc_sc_chi(
                            elems_fs[:, flag_2d[:, ind]],
                            unit_cell_parameters=unit_cell_parameters, flag_unit_cell_parameters=False)[0]
                    l_sc_exchange.append(sc_exchange)
                atom_sc_exchange = numpy.stack(l_sc_exchange, axis=-1)
                ddict["atom_sc_exchange"] = atom_sc_exchange

        if atom_site_moment is not None:
            atom_ordered_label = numpy.array(atom_site_moment.label, dtype=str)
            atom_ordered_moment_crystalaxis_xyz = numpy.array([
                atom_site_moment.crystalaxis_x, 
                atom_site_moment.crystalaxis_y, 
                atom_site_moment.crystalaxis_z], dtype=float)
            flags_atom_ordered_moment_crystalaxis_xyz = numpy.array([
                atom_site_moment.crystalaxis_x_refinement, 
                atom_site_moment.crystalaxis_y_refinement, 
                atom_site_moment.crystalaxis_z_refinement], dtype=bool)
            ddict["atom_ordered_label"] = atom_ordered_label
            ddict["atom_ordered_moment_crystalaxis_xyz"] = atom_ordered_moment_crystalaxis_xyz
            ddict["flags_atom_ordered_moment_crystalaxis_xyz"] = flags_atom_ordered_moment_crystalaxis_xyz
            if "atom_label" in ddict.keys():
                flag_2d = numpy.expand_dims(ddict["atom_label"], axis=1) == numpy.expand_dims(atom_ordered_label, axis=0)
                args = numpy.argwhere(flag_2d)
                atom_ordered_index = args[args[:,1].argsort()][:,0]
                ddict["atom_ordered_index"] = atom_ordered_index

        if isinstance(self, Crystal):
            # FIXME: beta should be defined in atom_site_aniso_beta, b_iso for every atom
            atom_b_iso, atom_beta = self.calc_b_iso_beta()
            ddict["atom_b_iso"] = atom_b_iso
            ddict["atom_beta"] = atom_beta.transpose()
            ddict["flags_atom_beta"] = self.get_flags_atom_beta()

        if atom_site_scat is not None:
            mag_atom_label = numpy.array(atom_site_scat.label, dtype=str)
            mag_atom_lande_factor = numpy.array(atom_site_scat.lande, dtype=float)
            mag_atom_kappa = numpy.array(atom_site_scat.kappa, dtype=float)
            flags_mag_atom_lande_factor = atom_site_scat.get_flags_lande()
            flags_mag_atom_kappa = atom_site_scat.get_flags_kappa()

            ddict["mag_atom_label"] = mag_atom_label
            ddict["mag_atom_lande_factor"] = mag_atom_lande_factor
            ddict["mag_atom_kappa"] = mag_atom_kappa
            ddict["flags_mag_atom_lande_factor"] = flags_mag_atom_lande_factor
            ddict["flags_mag_atom_kappa"] = flags_mag_atom_kappa
            # l_ats = atom_site_scat.atom_type_scat
            if "atom_label" in ddict.keys():
                flag_2d = numpy.expand_dims(ddict["atom_label"], axis=1) == numpy.expand_dims(mag_atom_label, axis=0)
                args = numpy.argwhere(flag_2d)
                mag_atom_index = args[args[:,1].argsort()][:,0]
                ddict["mag_atom_index"] = mag_atom_index

                if atom_site_scat.is_attribute("type_symbol"):
                    mag_type_symbol = numpy.array(atom_site_scat.type_symbol, dtype=str)
                else:
                    as_tsymbol = numpy.array(atom_site.type_symbol, dtype=str)
                    mag_type_symbol = as_tsymbol[mag_atom_index] 

                j0_parameters, j2_parameters = get_j0_j2_parameters(mag_type_symbol)
                ddict["mag_atom_j0_parameters"] = j0_parameters
                ddict["mag_atom_j2_parameters"] = j2_parameters
            
            if "atom_para_label" in ddict.keys():
                flag_2d = numpy.expand_dims(ddict["atom_para_label"], axis=0) == numpy.expand_dims(mag_atom_label, axis=1)
                args = numpy.argwhere(flag_2d)
                mag_atom_para_index = args[args[:,1].argsort()][:,0]
                ddict["mag_atom_para_index"] = mag_atom_para_index

            if "atom_ordered_label" in ddict.keys():
                flag_2d = numpy.expand_dims(ddict["atom_ordered_label"], axis=0) == numpy.expand_dims(mag_atom_label, axis=1)
                args = numpy.argwhere(flag_2d)
                mag_atom_ordered_index = args[args[:,1].argsort()][:,0]
                ddict["mag_atom_ordered_index"] = mag_atom_ordered_index

        if atom_type is not None:
            pass

        if atom_type_scat is not None:
            pass

        if atom_site_aniso is not None:
            atom_site_aniso_label = numpy.array(atom_site_aniso.label, dtype=str)
            flag_2d = numpy.expand_dims(ddict["atom_label"], axis=1) == numpy.expand_dims(atom_site_aniso_label, axis=0)
            args = numpy.argwhere(flag_2d)
            atom_site_aniso_index = args[args[:,1].argsort()][:,0]
            ddict["atom_site_aniso_index"] = atom_site_aniso_index
            #FIXME: add transformation U_ani and others to beta

            if "full_mcif_elems" in ddict.keys():
                elems_fs = ddict["full_mcif_elems"]
            elif "full_symm_elems" in ddict.keys():
                elems_fs = ddict["full_symm_elems"]
            else:
                raise AttributeError("Symmetry elements are not found.")

            flag_2d = ddict["atom_symm_elems_flag"][:, atom_site_aniso_index]
            l_sc_beta = []
            atom_site_aniso_label = numpy.array(atom_site_aniso.label, dtype=str)

            for ind in range(atom_site_aniso_label.size):
                symm_elems = elems_fs[:, flag_2d[:, ind]]
                sc_beta = calc_sc_beta(symm_elems)
                l_sc_beta.append(sc_beta)
            ddict["atom_site_aniso_sc_beta"] = numpy.stack(l_sc_beta, axis=-1)


        if atom_electron_configuration is not None:
            d_slater_orbitals = DATABASE["Orbitals by radial Slater functions"]
            for aec_item in atom_electron_configuration.items:
                aec_label = aec_item.label
                at_symbol = atom_site[aec_label].type_symbol
                element_symbol, charge, isotope_number = get_atomic_symbol_ion_charge_isotope_number_by_ion_symbol(at_symbol)
                core_shells_populations = aec_item.get_core_shells_populations()
                l_population = []
                l_n = []
                l_zeta = []
                l_coeff = []
                for sh_pop in core_shells_populations:
                    shell, population = sh_pop
                    n, dzeta, coeff = d_slater_orbitals[(element_symbol, shell)]
                    l_population.append(population)
                    l_n.append(n)
                    l_zeta.append(dzeta)
                    l_coeff.append(coeff)
                if len(l_population) > 0:
                    dict_shell = {"core_population": l_population,
                        "core_n": l_n,
                        "core_zeta": l_zeta,
                        "core_coeff": l_coeff}
                    ddict[f"shell_{aec_label:}"] = dict_shell
        return ddict

    def take_parameters_from_dictionary(self, ddict_crystal, l_parameter_name: list=None, l_sigma: list=None):
        keys = ddict_crystal.keys()
        if l_parameter_name is None:
            l_parameter_name = []
            l_sigma = []

        if "atom_para_susceptibility" in keys:
            mag_atom_susceptibility = (ddict_crystal["atom_para_sc_chi"]  * ddict_crystal["atom_para_susceptibility"][na, :, :]).sum(axis=1)
            for i_item, item_ass in enumerate(self.atom_site_susceptibility.items):
                item_ass.chi_11 = mag_atom_susceptibility[0, i_item]
                item_ass.chi_22 = mag_atom_susceptibility[1, i_item]
                item_ass.chi_33 = mag_atom_susceptibility[2, i_item]
                item_ass.chi_12 = mag_atom_susceptibility[3, i_item]
                item_ass.chi_13 = mag_atom_susceptibility[4, i_item]
                item_ass.chi_23 = mag_atom_susceptibility[5, i_item]

        if "atom_exchange_tensor" in keys:
            j_tensor = ddict_crystal["atom_exchange_tensor"]
            atom_sc_exchange = ddict_crystal["atom_sc_exchange"]
            j_tensor = numpy.sum(numpy.expand_dims(j_tensor, axis=0)*atom_sc_exchange, axis=1)
            for i_item, item_ase in enumerate(self.atom_site_exchange.items):
                item_ase.j_11 = j_tensor[0, i_item]
                item_ase.j_22 = j_tensor[1, i_item]
                item_ase.j_33 = j_tensor[2, i_item]
                item_ase.j_12 = j_tensor[3, i_item]
                item_ase.j_13 = j_tensor[4, i_item]
                item_ase.j_23 = j_tensor[5, i_item]

        if "atom_ordered_moment_crystalaxis_xyz" in keys:
            hh = ddict_crystal["atom_ordered_moment_crystalaxis_xyz"]
            for i_item, item in enumerate(self.atom_site_moment.items):
                item.crystalaxis_x = hh[0, i_item]
                item.crystalaxis_y = hh[1, i_item]
                item.crystalaxis_z = hh[2, i_item]

        if "mag_atom_lande_factor" in keys:
            hh = ddict_crystal["mag_atom_lande_factor"]
            for i_item, item in enumerate(self.atom_site_scat.items):
                item.lande = hh[i_item]

        if "mag_atom_kappa" in keys:
            hh = ddict_crystal["mag_atom_kappa"]
            for i_item, item in enumerate(self.atom_site_scat.items):
                item.kappa = hh[i_item]

        if "atom_fract_xyz" in keys:
            hh = ddict_crystal["atom_fract_xyz"]
            atom_site_sc_fract = ddict_crystal["atom_site_sc_fract"] 
            atom_site_sc_b = ddict_crystal["atom_site_sc_b"] 
            hh = calc_m_v(atom_site_sc_fract, numpy.mod(hh, 1), flag_m=False, flag_v=False)[0] + atom_site_sc_b

            for i_item, item in enumerate(self.atom_site.items):
                item.fract_x = numpy.round(numpy.mod(hh[0, i_item], 1.), decimals=6)
                item.fract_y = numpy.round(numpy.mod(hh[1, i_item], 1.), decimals=6)
                item.fract_z = numpy.round(numpy.mod(hh[2, i_item], 1.), decimals=6)

        if "atom_occupancy" in keys:
            hh = ddict_crystal["atom_occupancy"]
            for i_item, item in enumerate(self.atom_site.items):
                item.occupancy = numpy.round(hh[i_item], decimals=6)

        if "atom_multiplicity" in keys:
            hh = ddict_crystal["atom_multiplicity"]
            for i_item, item in enumerate(self.atom_site.items):
                item.multiplicity = hh[i_item]

        if ("atom_b_iso" in keys):
            atom_b_iso = ddict_crystal["atom_b_iso"]
            if "atom_beta" in keys:
                atom_beta = ddict_crystal["atom_beta"]

                if "atom_site_aniso_sc_beta" in keys:
                    atom_site_aniso_sc_beta = ddict_crystal["atom_site_aniso_sc_beta"]
                    atom_site_aniso_index = ddict_crystal["atom_site_aniso_index"]
                    atom_sc_beta = numpy.zeros((6,)+atom_beta.shape, dtype=float)
                    atom_sc_beta[:, :, atom_site_aniso_index] = atom_site_aniso_sc_beta
                    atom_beta = (atom_sc_beta*numpy.expand_dims(atom_beta, axis=0)).sum(axis=1)

            else:
                atom_beta = numpy.zeros((6, )+ atom_b_iso.shape, dtype=float)
            unit_cell_parameters = ddict_crystal["unit_cell_parameters"]
            if not(self.atom_site.is_attribute("adp_type")):
                for item in self.atom_site:
                    item.adp_type = "Biso"
                    if not(item.is_attribute("b_iso_or_equiv")):
                        item.b_iso_or_equiv = 0.
            atom_adp_type = numpy.array(self.atom_site.adp_type, dtype=str)
            atom_iso_param, atom_aniso_param = calc_param_iso_aniso_by_b_iso_beta(unit_cell_parameters, atom_adp_type, atom_b_iso, atom_beta)


            if self.is_attribute("atom_site"):
                atom_site = self.atom_site
                for i_item, item in enumerate(self.atom_site.items):
                    if item.adp_type in ("Uiso", "Uovl", "Umpe"):
                        item.u_iso_or_equiv = atom_iso_param[i_item]
                    elif item.adp_type in ("Biso", "Bovl"):
                        item.b_iso_or_equiv = atom_iso_param[i_item]

                if self.is_attribute("atom_site_aniso"):
                    atom_site_aniso_index = ddict_crystal["atom_site_aniso_index"]

                    for item, index_as in zip(self.atom_site_aniso.items, atom_site_aniso_index):
                        aniso_param = atom_aniso_param[:, index_as]
                        item_as = self.atom_site.items[index_as]
                        if item_as.adp_type in ("Uani", ):
                            item.u_11 = aniso_param[0]
                            item.u_22 = aniso_param[1]
                            item.u_33 = aniso_param[2]
                            item.u_12 = aniso_param[3]
                            item.u_13 = aniso_param[4]
                            item.u_23 = aniso_param[5]
                        elif item_as.adp_type in ("Bani", ):
                            item.b_11 = aniso_param[0]
                            item.b_22 = aniso_param[1]
                            item.b_33 = aniso_param[2]
                            item.b_12 = aniso_param[3]
                            item.b_13 = aniso_param[4]
                            item.b_23 = aniso_param[5]

        if "unit_cell_parameters" in keys:
            hh = ddict_crystal["unit_cell_parameters"]
            sc_uc = ddict_crystal["sc_uc"]
            v_uc = ddict_crystal["v_uc"]
            hhh = numpy.dot(sc_uc, hh) + v_uc
            cell = self.cell
            cell.length_a = hhh[0]
            cell.length_b = hhh[1]
            cell.length_c = hhh[2]
            cell.angle_alpha = numpy.round(hhh[3]*180./numpy.pi, decimals=5)
            cell.angle_beta = numpy.round(hhh[4]*180./numpy.pi, decimals=5)
            cell.angle_gamma = numpy.round(hhh[5]*180./numpy.pi, decimals=5)

        for name, sigma in zip(l_parameter_name, l_sigma):
            parameter_label, ind_s = name
            if "atom_para_susceptibility" in parameter_label:
                ind_chi, ind_a = ind_s[0], ind_s[1]
                item_ass = self.atom_site_susceptibility.items[ind_a]
                if ind_chi == 0:
                    item_ass.chi_11_sigma = float(sigma)
                elif ind_chi == 1:
                    item_ass.chi_22_sigma = float(sigma)
                elif ind_chi == 2:
                    item_ass.chi_33_sigma = float(sigma)
                elif ind_chi == 3:
                    item_ass.chi_12_sigma = float(sigma)
                elif ind_chi == 4:
                    item_ass.chi_13_sigma = float(sigma)
                elif ind_chi == 5:
                    item_ass.chi_23_sigma = float(sigma)
            
            if "atom_exchange_tensor" in parameter_label:
                ind_j, ind_a = ind_s[0], ind_s[1]
                item = self.atom_site_exchange.items[ind_a]
                if ind_j == 0:
                    item.j_11_sigma = float(sigma)
                elif ind_j == 1:
                    item.j_22_sigma = float(sigma)
                elif ind_j == 2:
                    item.j_33_sigma = float(sigma)
                elif ind_j == 3:
                    item.j_12_sigma = float(sigma)
                elif ind_j == 4:
                    item.j_13_sigma = float(sigma)
                elif ind_j == 5:
                    item.j_23_sigma = float(sigma)
            
            if "atom_fract_xyz" in parameter_label:
                ind_p, ind_a = ind_s[0], ind_s[1]
                item = self.atom_site.items[ind_a]
                if ind_p == 0:
                    item.fract_x_sigma = float(sigma)
                elif ind_p == 1:
                    item.fract_y_sigma = float(sigma)
                elif ind_p == 2:
                    item.fract_z_sigma = float(sigma)
            if "atom_occupancy" in parameter_label:
                ind_a = ind_s[0]
                item = self.atom_site.items[ind_a]
                item.occupancy_sigma = float(sigma)

            if "atom_b_iso" in parameter_label:
                ind_a = ind_s[0]
                item = self.atom_site.items[ind_a]
                if item.adp_type in ("Uiso", "Uovl", "Umpe"):
                    item.u_iso_or_equiv_sigma = float(sigma)/float(8.*numpy.square(numpy.pi))
                elif item.adp_type in ("Biso", "Bovl", "Bmpe"):
                    item.b_iso_or_equiv_sigma = float(sigma)

            if "atom_beta" in parameter_label:
                ind_p, ind_a = ind_s[0], ind_s[1]
                item = self.atom_site_aniso[self.atom_site.items[ind_a].label]
                if self.atom_site[item.label].adp_type in ("Bani", ):
                    if ind_p == 0:
                        item.b_11_sigma = float(sigma)
                    elif ind_p == 1:
                        item.b_22_sigma = float(sigma)
                    elif ind_p == 2:
                        item.b_33_sigma = float(sigma)
                    elif ind_p == 3:
                        item.b_12_sigma = float(sigma)
                    elif ind_p == 4:
                        item.b_13_sigma = float(sigma)
                    elif ind_p == 5:
                        item.b_23_sigma = float(sigma)

                if self.atom_site[item.label].adp_type in ("Uani", ):
                    reciprocal_unit_cell_parameters = calc_reciprocal_by_unit_cell_parameters(unit_cell_parameters, flag_unit_cell_parameters=False)[0]
                    beta_hh = numpy.zeros((6,1), dtype=float)
                    coeff_hh = calc_u_ij_by_beta(beta_hh, reciprocal_unit_cell_parameters, flag_beta=True)[1]["beta"]
                    if ind_p == 0:
                        item.u_11_sigma = float(sigma*coeff_hh[0])
                    elif ind_p == 1:
                        item.u_22_sigma = float(sigma*coeff_hh[1])
                    elif ind_p == 2:
                        item.u_33_sigma = float(sigma*coeff_hh[2])
                    elif ind_p == 3:
                        item.u_12_sigma = float(sigma*coeff_hh[3])
                    elif ind_p == 4:
                        item.u_13_sigma = float(sigma*coeff_hh[4])
                    elif ind_p == 5:
                        item.u_23_sigma = float(sigma*coeff_hh[5])


            if "atom_ordered_moment_crystalaxis_xyz" in parameter_label:
                ind_p, ind_a = ind_s[0], ind_s[1]
                item = self.atom_site_moment.items[ind_a]
                if ind_p == 0:
                    item.crystalaxis_x_sigma = float(sigma)
                elif ind_p == 1:
                    item.crystalaxis_y_sigma = float(sigma)
                elif ind_p == 2:
                    item.crystalaxis_z_sigma = float(sigma)

            if "mag_atom_kappa" in parameter_label:
                ind_a = ind_s[0]
                item = self.atom_site_scat.items[ind_a]
                item.kappa_sigma = float(sigma)

            if "mag_atom_lande_factor" in parameter_label:
                ind_a = ind_s[0]
                item = self.atom_site_scat.items[ind_a]
                item.lande_sigma = float(sigma)

            if "unit_cell_parameters" in parameter_label:
                ind_p = ind_s[0]
                item = self.cell
                if ind_p == 0:
                    item.length_a_sigma = float(sigma)
                elif ind_p == 1:
                    item.length_b_sigma = float(sigma)
                elif ind_p == 2:
                    item.length_c_sigma = float(sigma)
                elif ind_p == 3:
                    item.angle_alpha_sigma = float(sigma)*180./numpy.pi
                elif ind_p == 4:
                    item.angle_beta_sigma = float(sigma)*180./numpy.pi
                elif ind_p == 5:
                    item.angle_gamma_sigma = float(sigma)*180./numpy.pi



def calc_sc_v_unit_cell_parameters(type_cell: str, it_coordinate_system_code: str):
    """
       Calculate symmetry constraints on unit_cell_parameters as
       unit_cell_parameters = sc_uc * unit_cell_parameters + v_uc
    """
    sc_uc = numpy.eye(6, dtype=float)
    v_uc = numpy.zeros((6,), dtype=float)

    if type_cell == "aP":
        pass
    elif type_cell.startswith("m"):
        sc_uc[3, 3] = 0.
        sc_uc[5, 5] = 0.
        v_uc[3] = numpy.pi/2.
        v_uc[5] = numpy.pi/2.
    elif type_cell.startswith("o"):
        sc_uc[3, 3] = 0.
        sc_uc[4, 4] = 0.
        sc_uc[5, 5] = 0.
        v_uc[3] = numpy.pi/2.
        v_uc[4] = numpy.pi/2.
        v_uc[5] = numpy.pi/2.
    elif ((type_cell.startswith("t"))):  # FIXME: check  | (type_cell == "hP")
        sc_uc[0, 0], sc_uc[0, 1] = 1., 0.
        sc_uc[1, 0], sc_uc[1, 1] = 1., 0.
        sc_uc[3, 3] = 0.
        sc_uc[4, 4] = 0.
        sc_uc[5, 5] = 0.
        v_uc[3] = numpy.pi/2.
        v_uc[4] = numpy.pi/2.
        v_uc[5] = numpy.pi/2.
    elif ((type_cell.startswith("hP"))):
        if it_coordinate_system_code.lower() == "h":
            sc_uc[0, 0], sc_uc[0, 1] = 1., 0.
            sc_uc[1, 0], sc_uc[1, 1] = 1., 0.
            sc_uc[3, 3] = 0.
            sc_uc[4, 4] = 0.
            sc_uc[5, 5] = 0.
            v_uc[3] = numpy.pi/2.
            v_uc[4] = numpy.pi/2.
            v_uc[5] = numpy.pi*2./3.
        else:
            sc_uc[0, 0], sc_uc[0, 1], sc_uc[0, 2] = 1., 0., 0.
            sc_uc[1, 0], sc_uc[1, 1], sc_uc[1, 2] = 1., 0., 0.
            sc_uc[2, 0], sc_uc[2, 1], sc_uc[2, 2] = 1., 0., 0.
            sc_uc[3, 3], sc_uc[3, 4], sc_uc[3, 5] = 1., 0., 0.
            sc_uc[4, 3], sc_uc[4, 4], sc_uc[4, 5] = 1., 0., 0.
            sc_uc[5, 3], sc_uc[5, 4], sc_uc[5, 5] = 1., 0., 0.
    elif (type_cell == "hR"):
        if it_coordinate_system_code.lower() == "h":
            sc_uc[0, 0], sc_uc[0, 1] = 1., 0.
            sc_uc[1, 0], sc_uc[1, 1] = 1., 0.
            sc_uc[3, 3] = 0.
            sc_uc[4, 4] = 0.
            sc_uc[5, 5] = 0.
            v_uc[3] = numpy.pi/2.
            v_uc[4] = numpy.pi/2.
            v_uc[5] = numpy.pi*2./3.
        else:
            sc_uc[0, 0], sc_uc[0, 1], sc_uc[0, 2] = 1., 0., 0.
            sc_uc[1, 0], sc_uc[1, 1], sc_uc[1, 2] = 1., 0., 0.
            sc_uc[2, 0], sc_uc[2, 1], sc_uc[2, 2] = 1., 0., 0.
            sc_uc[3, 3], sc_uc[3, 4], sc_uc[3, 5] = 1., 0., 0.
            sc_uc[4, 3], sc_uc[4, 4], sc_uc[4, 5] = 1., 0., 0.
            sc_uc[5, 3], sc_uc[5, 4], sc_uc[5, 5] = 1., 0., 0.
    elif type_cell.startswith("c"):
        sc_uc[0, 0], sc_uc[0, 1], sc_uc[0, 2] = 1., 0., 0.
        sc_uc[1, 0], sc_uc[1, 1], sc_uc[1, 2] = 1., 0., 0.
        sc_uc[2, 0], sc_uc[2, 1], sc_uc[2, 2] = 1., 0., 0.
        sc_uc[3, 3] = 0.
        sc_uc[4, 4] = 0.
        sc_uc[5, 5] = 0.
        v_uc[3] = numpy.pi/2.
        v_uc[4] = numpy.pi/2.
        v_uc[5] = numpy.pi/2.
    
    return sc_uc, v_uc