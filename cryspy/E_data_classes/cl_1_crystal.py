"""Description of Crystal class."""

__author__ = "ikibalin"
__version__ = "2020_08_19"
import dataclasses
import logging
import numpy

from typing import NoReturn

from cryspy.A_functions_base.matrix_operations import (
    calc_m_v,
)
from cryspy.A_functions_base.unit_cell import (
    calc_eq_ccs_by_unit_cell_parameters,
    calc_reciprocal_by_unit_cell_parameters,
    calc_m_m_by_unit_cell_parameters,
)
from cryspy.A_functions_base.structure_factor import (
    calc_f_nucl_by_dictionary,
    calc_sft_ccs_by_dictionary,
    calc_f_m_perp_by_sft,
    calc_f_m_perp_ordered_by_dictionary,
)

from cryspy.A_functions_base.structure_scene import (
    calc_atoms_in_unit_cell, calc_structure_scene, calc_tensors_xyz_by_symmetry)
from cryspy.D_functions_item_loop.function_1_view_structure import (
    ATOM_DISPLAY_ADP, ATOM_DISPLAY_SUSCEPTIBILITY,
    StructureViewSettings,
    ipywidgets,
    make_settings_widgets,
    make_structure_view,
    make_structure_view_html,
    take_element,
)


from cryspy.A_functions_base.symmetry_elements import (
    calc_full_mag_elems,
)

from cryspy.B_parent_classes.cl_3_data import DataN

from cryspy.C_item_loop_classes.cl_1_atom_local_axes import AtomLocalAxesL
from cryspy.C_item_loop_classes.cl_1_atom_rho_multipole import (
    AtomRhoMultipoleL,
)
from cryspy.C_item_loop_classes.cl_1_atom_site import AtomSite, AtomSiteL
from cryspy.C_item_loop_classes.cl_1_atom_type import AtomTypeL
from cryspy.C_item_loop_classes.cl_1_atom_site_aniso import AtomSiteAnisoL
from cryspy.C_item_loop_classes.cl_1_atom_site_susceptibility import (
    AtomSiteSusceptibilityL,
)
from cryspy.C_item_loop_classes.cl_1_atom_site_exchange import (
    AtomSiteExchangeL,
)
from cryspy.C_item_loop_classes.cl_1_atom_type_scat import AtomTypeScatL
from cryspy.C_item_loop_classes.cl_2_atom_site_scat import AtomSiteScatL
from cryspy.C_item_loop_classes.cl_1_atom_site_moment import (
    AtomSiteMoment,
    AtomSiteMomentL,
)
from cryspy.C_item_loop_classes.cl_1_refln_susceptibility import (
    ReflnSusceptibilityL,
)
from cryspy.C_item_loop_classes.cl_1_atom_electron_configuration import (
    AtomElectronConfigurationL,
)
from cryspy.C_item_loop_classes.cl_1_cell import Cell
from cryspy.C_item_loop_classes.cl_1_refln import ReflnL
from cryspy.C_item_loop_classes.cl_1_space_group_symop_magn_centering import (
    SpaceGroupSymopMagnCenteringL,
)

from cryspy.C_item_loop_classes.cl_2_space_group_symop_magn_operation import (
    SpaceGroupSymopMagnOperationL,
)
from cryspy.C_item_loop_classes.cl_2_space_group import SpaceGroup


from cryspy.D_functions_item_loop.function_1_report_magnetization_ellipsoid import (
    magnetization_ellipsoid_by_u_ij,
    report_main_axes_of_magnetization_ellipsoids,
)
from cryspy.A_functions_base.function_1_deprecation import deprecated
from cryspy.E_data_classes.cl_1_crystal_dictionary import (
    calc_sc_v_unit_cell_parameters,  # noqa: F401 (kept importable from here)
    form_crystal_dictionary,
    save_crystal_multipole_density_to_den,
    save_crystal_multipole_magnetization_density_to_den,
    take_crystal_parameters_from_dictionary,
)

logger = logging.getLogger(__name__)

na = numpy.newaxis


class Crystal(DataN):
    """The block ‘Crystal’ describes the crystal structure.

    - The atom_site loop provides the positions of atoms within the unit cell.
    - The cell defines the unit cell.

    The space group can be specified using:
    - space_group for crystal structures.
    - space_group_symop_magn_operation and space_group_symop_magn_centering for magnetic structures (MCIF format).

    The atom_site_susceptibility loop describes the paramagnetic properties of the crystal and can be defined in both nuclear and magnetic space groups.

    The atom_site_moment loop describes the ordered magnetic moments of the crystal and is only required if the magnetic space group (MCIF format) is defined.

    To define the magnetic form factor in the dipolar approximation, use the atom_site_scat loop.

    The anisotropic magnetic form factor can be defined using the multipole model.
     In this case, use the atom_rho_multipole loop and the atom_local_axes loop.
     The latter is necessary to define the local coordinate system in which the multipole functions are defined.

    Please note that the programme currently does not check the symmetry constraint on multipole parameters. Therefore, the user should verify this themselves.
    """

    CLASSES_MANDATORY = (
        Cell,
        AtomSiteL,
    )
    CLASSES_OPTIONAL = (
        SpaceGroup,
        AtomTypeL,
        AtomSiteAnisoL,
        AtomSiteSusceptibilityL,
        AtomSiteExchangeL,
        AtomSiteScatL,
        AtomTypeScatL,
        AtomLocalAxesL,
        AtomElectronConfigurationL,
        AtomSiteMomentL,
        SpaceGroupSymopMagnOperationL,
        SpaceGroupSymopMagnCenteringL,
        AtomRhoMultipoleL,
    )
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
            cell.it_coordinate_system_code = (
                space_group.it_coordinate_system_code
            )
            cell.apply_constraints()

            atom_site = self.atom_site
            atom_site.apply_constraints(space_group_wyckoff)

            if self.is_attribute("atom_site_aniso"):
                atom_site_aniso = self.atom_site_aniso
                atom_site_aniso.apply_space_group_constraint(
                    atom_site, space_group
                )

            if self.is_attribute("atom_site_susceptibility"):
                atom_site_susceptibility = self.atom_site_susceptibility
                atom_site_susceptibility.apply_chi_iso_constraint(cell)
                atom_site_susceptibility.apply_space_group_constraint(
                    atom_site, space_group, cell
                )
            if self.is_attribute("atom_site_exchange"):
                atom_site_exchange = self.atom_site_exchange
                atom_site_exchange.apply_j_iso_constraint()
                atom_site_exchange.apply_space_group_constraint(
                    atom_site, space_group, cell
                )

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
        coeff = float(8.0 * numpy.pi**2)
        cell = self.cell
        for item_a_s in a_s.items:
            label_atom = item_a_s.label
            try:
                adp_type = item_a_s.adp_type
            except AttributeError:
                adp_type = None
            b_iso = 0.0
            beta = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
            if adp_type == "Uiso":
                u_iso = float(item_a_s.u_iso_or_equiv)
                b_iso = float(8.0 * numpy.pi**2 * u_iso)
            elif adp_type == "Biso":
                b_iso = float(item_a_s.b_iso_or_equiv)
            elif adp_type == "Uovl":
                # FIXME: correct it
                u_iso = float(item_a_s.u_iso_or_equiv)
                b_iso = coeff * u_iso
            elif adp_type == "Umpe":
                # FIXME: correct it
                u_iso = float(item_a_s.u_iso_or_equiv)
                b_iso = float(8.0 * numpy.pi**2 * u_iso)
            elif adp_type == "Uani":
                item_a_s_a = a_s_a[label_atom]
                beta = item_a_s_a.calc_beta(cell)
            elif adp_type == "Bovl":
                # FIXME: correct it
                b_iso = float(item_a_s.b_iso_or_equiv)
            elif adp_type == "Bani":
                item_a_s_a = a_s_a[label_atom]
                beta = (
                    float(item_a_s_a.b_11),
                    float(item_a_s_a.b_22),
                    float(item_a_s_a.b_33),
                    float(item_a_s_a.b_12),
                    float(item_a_s_a.b_13),
                    float(item_a_s_a.b_23),
                )
            l_b_iso.append(b_iso)
            l_beta.append(beta)
        np_b_iso = numpy.array(l_b_iso, dtype=float)
        np_beta = numpy.array(l_beta, dtype=float)
        return np_b_iso, np_beta

    def calc_f_nucl(self, index_hkl):
        dict_crystal = self.get_dictionary()
        dict_in_out = {"index_hkl": index_hkl}
        f_nucl, dder = calc_f_nucl_by_dictionary(
            dict_crystal, dict_in_out, flag_use_precalculated_data=False
        )
        return f_nucl

    def calc_f_mag(
        self,
        index_hkl,
        magnetic_field_ccs: numpy.ndarray = numpy.array(
            [0.0, 0.0, 0.0], dtype=float
        ),
    ):
        flag_use_precalculated_data = False
        dict_crystal = self.get_dictionary()
        dict_crystal_keys = dict_crystal.keys()
        dict_in_out = {"index_hkl": index_hkl, "flag_only_orbital": False}

        f_mag = numpy.zeros(index_hkl.shape, dtype=numpy.complex128)
        if "atom_para_index" in dict_crystal_keys:
            sft_ccs = calc_sft_ccs_by_dictionary(
                dict_crystal,
                dict_in_out,
                flag_use_precalculated_data=flag_use_precalculated_data,
            )[0]
            f_mag_para = calc_m_v(
                sft_ccs, magnetic_field_ccs, flag_m=False, flag_v=False
            )[0]
            f_mag += f_mag_para
        if "atom_ordered_index" in dict_crystal_keys:
            calc_f_m_perp_ordered_by_dictionary(
                dict_crystal,
                dict_in_out,
                flag_use_precalculated_data=flag_use_precalculated_data,
            )
            if "f_m_o" in dict_in_out.keys():
                f_m_o = dict_in_out["f_m_o"]
                f_mag += f_m_o
        return f_mag

    def calc_refln(self, index_hkl, flag_internal: bool = True):
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
        self,
        index_hkl,
        flag_only_orbital: bool = False,
        dict_in_out: dict = None,
    ):
        """Calculate structure factor tensor in CCS coordinate (X||a*, Z||c)."""

        dict_crystal = self.get_dictionary()
        if dict_in_out is None:
            dict_in_out = {
                "index_hkl": index_hkl,
                "flag_only_orbital": flag_only_orbital,
            }
        else:
            dict_in_out["index_hkl"] = index_hkl
            dict_in_out["flag_only_orbital"] = flag_only_orbital

        sft_ccs, dder = calc_sft_ccs_by_dictionary(
            dict_crystal, dict_in_out, flag_use_precalculated_data=False
        )
        return sft_ccs, dder

    def calc_refln_susceptibility(
        self,
        index_hkl,
        flag_internal: bool = True,
        flag_only_orbital: bool = False,
    ):
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
            index_hkl, flag_only_orbital=flag_only_orbital
        )

        res = ReflnSusceptibilityL(loop_name=self.data_name)
        res.numpy_index_h = index_hkl[0, :]
        res.numpy_index_k = index_hkl[1, :]
        res.numpy_index_l = index_hkl[2, :]
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

    @deprecated
    def calc_f_m_perp(
        self, index_hkl, magnetic_field, dict_in_out: dict = None
    ):
        """Calculate F_M_perpendicular"""
        sft_ccs, dder_sft_ccs = self.calc_structure_factor_tensor_ccs(
            index_hkl, dict_in_out=dict_in_out
        )
        unit_cell_parameters = self.cell.get_unit_cell_parameters()
        eq_ccs, dder_eq_ccs = calc_eq_ccs_by_unit_cell_parameters(
            index_hkl, unit_cell_parameters
        )
        f_m_perp, dder_f_m_perp = calc_f_m_perp_by_sft(
            sft_ccs,
            magnetic_field,
            eq_ccs,
            flag_sft_ccs=False,
            flag_magnetic_field=False,
            flag_eq_ccs=False,
        )
        return f_m_perp

    def calc_hkl(self, sthol_min: float = 0.0, sthovl_max: float = 1.0):
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

    def calc_hkl_in_range(
        self, sthol_min: float = 0.0, sthovl_max: float = 1.0
    ):
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
            return []

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

        l_moments, l_moments_sigma, l_rot_matrix = (
            a_s_s.calc_main_axes_of_magnetization_ellipsoid(cell)
        )
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
        s_out = report_main_axes_of_magnetization_ellipsoids(cell, a_s_m_a)

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

        try:
            spgr = (
                self.space_group
            )  # FIXME: it should be also the case of magnetic space group
            cell = self.cell
            a_s = self.atom_site
            a_s_m_a = self.atom_site_susceptibility
        except AttributeError:
            return None

        m_m_norm = cell.m_m_norm
        np_field_xyz = numpy.matmul(m_m_norm, np_field)
        norm_field = numpy.sqrt(numpy.square(np_field_xyz).sum())
        if norm_field == 0.0:
            norm_field = 1.0
        np_field = np_field / norm_field
        l_lab_out, l_xyz_out, l_moment_out = [], [], []

        m_mt_norm_m_norm_field = numpy.matmul(
            numpy.matmul(m_m_norm.transpose(), m_m_norm), np_field
        )

        as_p1 = []
        asm_p1 = []
        for _l, _11, _22, _33, _12, _13, _23 in zip(
            a_s_m_a.label,
            a_s_m_a.chi_11,
            a_s_m_a.chi_22,
            a_s_m_a.chi_33,
            a_s_m_a.chi_12,
            a_s_m_a.chi_13,
            a_s_m_a.chi_23,
        ):
            m_chi = numpy.array(
                [[_11, _12, _13], [_12, _22, _23], [_13, _23, _33]],
                dtype=float,
            )
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
                as_p1.append(
                    AtomSite(
                        label=label_p1,
                        fract_x=_xyz[0],
                        fract_y=_xyz[1],
                        fract_z=_xyz[2],
                    )
                )
                asm_p1.append(
                    AtomSiteMoment(
                        label=label_p1,
                        crystalaxis_x=numpy.round(_moment[0], 5),
                        crystalaxis_y=numpy.round(_moment[1], 5),
                        crystalaxis_z=numpy.round(_moment[2], 5),
                    )
                )

        spgr_p1 = SpaceGroup(it_number=1)
        spgr_p1.form_object()
        atom_site_p1 = AtomSiteL()
        atom_site_p1.items = as_p1
        atom_site_moment_p1 = AtomSiteMomentL()
        atom_site_moment_p1.items = asm_p1

        cryspy_p1 = Crystal(
            space_group=spgr_p1,
            cell=cell,
            atom_site=atom_site_p1,
            atom_site_moment=atom_site_moment_p1,
        )

        return cryspy_p1

    def report(self):
        return self.report_main_axes_of_magnetization_ellipsoids()

    def plots(self):
        return []

    def get_structure_scene(self, settings=None):
        """Atoms, magnetic moments and cell edges for a 3D view.

        Symmetry comes from the magnetic space group (mcif operations and
        centrings, with time reversal) or from the space group; moments
        from ``atom_site_moment`` (μB along a/|a|, b/|b|, c/|c|).

        Parameters
        ----------
        settings : StructureViewSettings
            Number of cells and what to show; the defaults if None.

        Returns
        -------
        scene : dict
            See ``structure_scene.calc_structure_scene``; also
            ``atom_element`` and ``atom_label`` of the shown atoms of the
            asymmetric unit (``atom_index`` refers to them).
        """
        if settings is None:
            settings = StructureViewSettings()
        atom_site = self.atom_site
        unit_cell_parameters = self.cell.get_unit_cell_parameters()
        full_symm_elems = self.get_full_symm_elems()
        atom_label = numpy.array(atom_site.label, dtype=str)
        atom_fract_xyz = numpy.array(
            [atom_site.fract_x, atom_site.fract_y, atom_site.fract_z],
            dtype=float)
        atom_moment = self.calc_moments_crystalaxis()
        atom_tensor = None
        if settings.atom_display == ATOM_DISPLAY_ADP:
            atom_tensor = self.calc_adp_tensors_xyz()
        elif settings.atom_display == ATOM_DISPLAY_SUSCEPTIBILITY:
            atom_tensor = self.calc_susceptibility_tensors_xyz()
        if atom_site.is_attribute("type_symbol"):
            symbols = numpy.array(atom_site.type_symbol, dtype=str)
        else:
            symbols = atom_label
        # dummy atoms (type_symbol 'Dummy', e.g. for local axes) are hidden
        flags = numpy.ones(atom_label.size, dtype=bool)
        if settings.hide_dummy_atoms:
            flags = numpy.char.lower(symbols) != "dummy"
            if not numpy.any(flags):
                flags[:] = True
        scene = calc_structure_scene(
            unit_cell_parameters, full_symm_elems, atom_fract_xyz[:, flags],
            atom_moment_crystalaxis=(
                None if atom_moment is None else atom_moment[:, flags]),
            atom_tensor_xyz=(
                None if atom_tensor is None else atom_tensor[:, :, flags]),
            n_cells=settings.n_cells,
            flag_boundary=settings.show_boundary_atoms,
            flag_inner_edges=settings.show_inner_cells)
        scene["atom_element"] = [take_element(symbol) for symbol in symbols[flags]]
        scene["atom_label"] = list(atom_label[flags])
        return scene

    def get_full_symm_elems(self):
        """Symmetry elements with centring for the atoms of the cell.

        From the magnetic space group (mcif operations and centrings, 14
        rows with time reversal), else from the space group (13 rows), else
        the identity.
        """
        if (self.is_attribute("space_group_symop_magn_operation") and
                self.is_attribute("space_group_symop_magn_centering")):
            return calc_full_mag_elems(
                self.space_group_symop_magn_operation.get_sym_elems(),
                self.space_group_symop_magn_centering.get_sym_elems())
        if self.is_attribute("space_group"):
            return self.space_group.full_space_group_symop.get_symm_elems()
        return numpy.array(
            [[0], [0], [0], [1], [1], [0], [0], [0], [1], [0], [0], [0], [1]],
            dtype=int)

    def calc_moments_crystalaxis(self):
        """Ordered moments of the atoms from ``atom_site_moment``.

        Returns
        -------
        numpy.ndarray, shape (3, n_atom), or None without atom_site_moment
            μB along a/|a|, b/|b|, c/|c|; NaN for atoms without a moment.
        """
        if not self.is_attribute("atom_site_moment"):
            return None
        atom_label = numpy.array(self.atom_site.label, dtype=str)
        res = numpy.full((3, atom_label.size), numpy.nan)
        for moment in self.atom_site_moment.items:
            res[:, atom_label == moment.label] = numpy.array([
                moment.crystalaxis_x, moment.crystalaxis_y,
                moment.crystalaxis_z], dtype=float)[:, numpy.newaxis]
        return res

    def calc_atoms_in_unit_cell(self) -> dict:
        """All atoms of the unit cell, made by the symmetry elements.

        Positions, ordered moments (axial vectors, with time reversal of
        the magnetic space group), displacement tensors and susceptibility
        tensors of the atoms equivalent to every atom of ``atom_site``.

        Returns
        -------
        dict
            ``label`` (n,) new labels 'Fe1_1', 'Fe1_2', ...;
            ``atom_index`` (n,) index in ``atom_site``;
            ``fract_xyz`` (3, n) in [0, 1);
            ``moment_crystalaxis`` (3, n) μB along a/|a|, b/|b|, c/|c| or None;
            ``u_xyz`` (3, 3, n) U in Å², Cartesian (x || a*, z || c);
            ``chi_xyz`` (3, 3, n) χ in μB/T, Cartesian;
            NaN where the atom has no such value.
        """
        ucp = self.cell.get_unit_cell_parameters()
        full_symm_elems = self.get_full_symm_elems()
        atom_site = self.atom_site
        atom_fract_xyz = numpy.array(
            [atom_site.fract_x, atom_site.fract_y, atom_site.fract_z],
            dtype=float)
        moment = self.calc_moments_crystalaxis()
        moment_fract = None if moment is None else moment / ucp[:3, numpy.newaxis]
        fract_xyz, atom_index, moment_cell, symm_index = calc_atoms_in_unit_cell(
            full_symm_elems, atom_fract_xyz, moment_fract)
        l_label, d_count = [], {}
        for i_atom in atom_index:
            label = atom_site.label[i_atom]
            d_count[label] = d_count.get(label, 0) + 1
            l_label.append(f"{label}_{d_count[label]}")
        tensors = {}
        for key, tensor in (("u_xyz", self.calc_adp_tensors_xyz()),
                            ("chi_xyz", self.calc_susceptibility_tensors_xyz())):
            tensors[key] = calc_tensors_xyz_by_symmetry(
                ucp, full_symm_elems, symm_index, tensor[:, :, atom_index])
        return {
            "label": l_label,
            "atom_index": atom_index,
            "fract_xyz": fract_xyz,
            "moment_crystalaxis": (None if moment_cell is None else
                                   moment_cell * ucp[:3, numpy.newaxis]),
            **tensors}

    def to_p1(self, data_name: str = None):
        """The crystal with all atoms of the unit cell in space group P1.

        The atoms equivalent by symmetry become separate atoms 'Fe1_1',
        'Fe1_2', ... with their positions, occupancies, displacement
        parameters (anisotropic ones rotated), ordered moments (a magnetic
        crystal gets the magnetic group P1 with x,y,z,+1), susceptibility
        tensors (rotated, type Cani) and ``atom_site_scat``. Multipoles and
        local axes are not converted.

        Parameters
        ----------
        data_name : str
            Data name of the new crystal; the old one with '_p1' if None.

        Returns
        -------
        Crystal
        """
        from cryspy.H_functions_global.function_1_cryspy_objects import str_to_globaln
        if data_name is None:
            data_name = f"{self.data_name}_p1"
        cell, atom_site = self.cell, self.atom_site
        ucp = cell.get_unit_cell_parameters()
        atoms = self.calc_atoms_in_unit_cell()
        index, labels = atoms["atom_index"], atoms["label"]
        m_m = calc_m_m_by_unit_cell_parameters(ucp)[0].reshape(3, 3)
        inv_m = numpy.linalg.inv(m_m)
        rec = calc_reciprocal_by_unit_cell_parameters(ucp)[0][:3]

        ls_out = [f"data_{data_name}", "", cell.to_cif(), ""]
        magnetic = atoms["moment_crystalaxis"] is not None
        if magnetic:
            ls_out += ["loop_", "_space_group_symop_magn_operation.xyz",
                       "_space_group_symop_magn_operation.id", "x,y,z,+1 1", "",
                       "loop_", "_space_group_symop_magn_centering.xyz",
                       "_space_group_symop_magn_centering.id", "x,y,z,+1 1", ""]
        else:
            ls_out += ['_space_group_name_H-M_alt "P 1"', ""]

        ls_site = ["loop_", "_atom_site_label", "_atom_site_type_symbol",
                   "_atom_site_fract_x", "_atom_site_fract_y",
                   "_atom_site_fract_z", "_atom_site_occupancy",
                   "_atom_site_adp_type", "_atom_site_B_iso_or_equiv"]
        ls_aniso = ["loop_", "_atom_site_aniso_label"] + [
            f"_atom_site_aniso_U_{ij}" for ij in ("11", "22", "33", "12", "13", "23")]
        b_iso = self.calc_b_iso_beta()[0]
        for label, i_atom, xyz, u_xyz in zip(
                labels, index, atoms["fract_xyz"].transpose(),
                numpy.moveaxis(atoms["u_xyz"], 2, 0)):
            item = atom_site.items[i_atom]
            occupancy = item.occupancy if item.is_attribute("occupancy") else 1.
            type_symbol = item.type_symbol if item.is_attribute("type_symbol") else item.label
            adp_type = item.adp_type if item.is_attribute("adp_type") else "Biso"
            if adp_type in ("Uani", "Bani") and not numpy.any(numpy.isnan(u_xyz)):
                sigma = inv_m @ u_xyz @ inv_m.transpose()
                u_cif = sigma / (rec[:, numpy.newaxis] * rec[numpy.newaxis, :])
                ls_aniso.append(f"{label} " + " ".join(
                    f"{u_cif[i, j]:.6f}" for i, j in
                    ((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2))))
                adp_type, b_value = "Uani", 0.
            else:
                adp_type, b_value = "Biso", b_iso[i_atom]
            ls_site.append(
                f"{label} {type_symbol} {xyz[0]:.6f} {xyz[1]:.6f} {xyz[2]:.6f} "
                f"{occupancy:.6f} {adp_type} {b_value:.6f}")
        ls_out += ls_site + [""]
        if len(ls_aniso) > 8:
            ls_out += ls_aniso + [""]
        if magnetic:
            ls_moment = ["loop_", "_atom_site_moment_label",
                         "_atom_site_moment_crystalaxis_x",
                         "_atom_site_moment_crystalaxis_y",
                         "_atom_site_moment_crystalaxis_z"]
            for label, moment in zip(labels, atoms["moment_crystalaxis"].transpose()):
                if not numpy.any(numpy.isnan(moment)):
                    ls_moment.append(f"{label} " + " ".join(f"{m:.6f}" for m in moment))
            if len(ls_moment) > 5:
                ls_out += ls_moment + [""]
        ls_chi = ["loop_", "_atom_site_susceptibility_label",
                  "_atom_site_susceptibility_chi_type"] + [
            f"_atom_site_susceptibility_chi_{ij}" for ij in ("11", "22", "33", "12", "13", "23")]
        for label, chi in zip(labels, numpy.moveaxis(atoms["chi_xyz"], 2, 0)):
            if not numpy.any(numpy.isnan(chi)):
                ls_chi.append(f"{label} Cani " + " ".join(
                    f"{chi[i, j]:.6f}" for i, j in
                    ((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2))))
        if len(ls_chi) > 9:
            ls_out += ls_chi + [""]
        if self.is_attribute("atom_site_scat"):
            scat = self.atom_site_scat
            d_scat = {item.label: item for item in scat.items}
            l_item = []
            for label, i_atom in zip(labels, index):
                item = d_scat.get(atom_site.items[i_atom].label)
                if item is not None:
                    l_item.append(f"{label} " + " ".join(
                        str(getattr(item, name)) for name in item.ATTR_NAMES[1:]
                        if item.is_attribute(name)))
            if len(l_item) > 0:
                names = [name for name in scat.items[0].ATTR_NAMES[1:]
                         if scat.items[0].is_attribute(name)]
                ls_out += ["loop_", "_atom_site_scat_label"] + [
                    f"_atom_site_scat_{scat.items[0].ATTR_CIF[scat.items[0].ATTR_NAMES.index(name)]}"
                    for name in names] + l_item + [""]
        obj = str_to_globaln("\n".join(ls_out))
        return obj.items[0]

    def report_ellipsoids(self, kind: str = "ADP", flag_unit_cell: bool = False,
                          field: float = 1.) -> str:
        """Principal axes of the atomic ellipsoids in markdown.

        Parameters
        ----------
        kind : str
            'ADP': displacement ellipsoids (U_i in Å², rms displacements);
            'susceptibility': magnetization ellipsoids (|χ_i| in μB/T and
            the moments induced by ``field``).
        flag_unit_cell : bool
            For every atom of the unit cell (orientations of the equivalent
            atoms) instead of the atoms of ``atom_site``.
        field : float
            Magnetic field, T, for the induced moments.

        Returns
        -------
        str
            Directions are unit vectors in Cartesian axes (x || a*, z || c)
            and in the crystal axes [u v w] (fractions of a, b, c).
        """
        kind_low = kind.strip().lower()
        if kind_low.startswith("adp") or kind_low.startswith("u") or kind_low.startswith("therm"):
            flag_adp = True
        elif kind_low.startswith("sus") or kind_low.startswith("chi") or kind_low.startswith("magn"):
            flag_adp = False
        else:
            raise ValueError(f"kind is '{kind}', not 'ADP' or 'susceptibility'.")
        if flag_unit_cell:
            atoms = self.calc_atoms_in_unit_cell()
            labels = atoms["label"]
            tensors = atoms["u_xyz"] if flag_adp else atoms["chi_xyz"]
        else:
            labels = list(self.atom_site.label)
            tensors = (self.calc_adp_tensors_xyz() if flag_adp else
                       self.calc_susceptibility_tensors_xyz())
        m_m = calc_m_m_by_unit_cell_parameters(
            self.cell.get_unit_cell_parameters())[0].reshape(3, 3)
        inv_m = numpy.linalg.inv(m_m)
        if flag_adp:
            ls_out = ["Displacement ellipsoids: U_i (Å²) and rms displacement "
                      "sqrt(U_i) (Å) along the principal axes.", "",
                      "| atom | U_eq, Å² | U_i, Å² | rms, Å | axis xyz (x||a*, z||c) | axis [u v w] |",
                      "|---|---|---|---|---|---|"]
        else:
            ls_out = ["Magnetization ellipsoids: |χ_i| (μB/T) along the principal "
                      f"axes and the induced moment at {field:g} T (μB).", "",
                      "| atom | |χ_i|, μB/T | M_i, μB | axis xyz (x||a*, z||c) | axis [u v w] |",
                      "|---|---|---|---|---|"]
        n_atom = 0
        for label, tensor in zip(labels, numpy.moveaxis(tensors, 2, 0)):
            if numpy.any(numpy.isnan(tensor)):
                continue
            n_atom += 1
            if flag_adp:
                val, vec = numpy.linalg.eigh(0.5 * (tensor + tensor.transpose()))
            else:
                val, vec = numpy.linalg.eigh(tensor.transpose() @ tensor)
                val = numpy.sqrt(numpy.clip(val, 0., None))
            for i_axis in range(3):
                axis = vec[:, i_axis]
                uvw = inv_m @ axis
                uvw = uvw / numpy.abs(uvw).max()
                s_xyz = " ".join(f"{v:.3f}" for v in axis)
                s_uvw = " ".join(f"{v:.3f}" for v in uvw)
                s_label = label if i_axis == 0 else ""
                if flag_adp:
                    s_eq = f"{numpy.trace(tensor) / 3.:.5f}" if i_axis == 0 else ""
                    rms = numpy.sqrt(max(val[i_axis], 0.))
                    ls_out.append(f"| {s_label} | {s_eq} | {val[i_axis]:.5f} | "
                                  f"{rms:.4f} | {s_xyz} | {s_uvw} |")
                else:
                    ls_out.append(f"| {s_label} | {val[i_axis]:.4f} | "
                                  f"{val[i_axis] * field:.4f} | {s_xyz} | {s_uvw} |")
        if n_atom == 0:
            if flag_adp:
                return "No atom has displacement parameters."
            return "No atom has a susceptibility tensor (atom_site_susceptibility)."
        if flag_adp:
            ls_out.append("")
            ls_out.append("A negative U_i means non-positive definite displacement parameters.")
        return "\n".join(ls_out)

    def calc_adp_tensors_xyz(self):
        """Atomic displacement tensors U of the atoms in Cartesian axes.

        From ``adp_type`` of ``atom_site`` (Uiso, Biso, Uani, Bani, ...)
        and ``atom_site_aniso``, as the structure factor uses them.

        Returns
        -------
        numpy.ndarray, shape (3, 3, n_atom)
            U in Å², x || a*, z || c; NaN for atoms without displacement
            parameters.
        """
        b_iso, beta = self.calc_b_iso_beta()
        m_m = calc_m_m_by_unit_cell_parameters(
            self.cell.get_unit_cell_parameters())[0].reshape(3, 3)
        res = numpy.full((3, 3, b_iso.size), numpy.nan)
        for i_atom, (b, bb) in enumerate(zip(b_iso, beta)):
            if numpy.any(bb != 0.):
                beta_m = numpy.array([[bb[0], bb[3], bb[4]],
                                      [bb[3], bb[1], bb[5]],
                                      [bb[4], bb[5], bb[2]]], dtype=float)
                # exp(-h^T beta h): beta = 2 pi^2 Sigma_fract
                res[:, :, i_atom] = m_m @ beta_m @ m_m.transpose() / (2. * numpy.pi**2)
            elif b > 0.:
                res[:, :, i_atom] = numpy.eye(3) * b / (8. * numpy.pi**2)
        return res

    def calc_susceptibility_tensors_xyz(self):
        """Susceptibility tensors of the atoms in Cartesian axes.

        From ``atom_site_susceptibility`` (given in x || a*, z || c).

        Returns
        -------
        numpy.ndarray, shape (3, 3, n_atom)
            χ in μB/T; NaN for atoms without susceptibility.
        """
        atom_label = numpy.array(self.atom_site.label, dtype=str)
        res = numpy.full((3, 3, atom_label.size), numpy.nan)
        if not self.is_attribute("atom_site_susceptibility"):
            return res
        for item in self.atom_site_susceptibility.items:
            if (item.is_attribute("chi_type") and
                    item.chi_type.lower().startswith("ciso")):
                chi = numpy.eye(3) * float(item.chi_11)
            else:
                chi = numpy.array([
                    [item.chi_11, item.chi_12, item.chi_13],
                    [item.chi_12, item.chi_22, item.chi_23],
                    [item.chi_13, item.chi_23, item.chi_33]], dtype=float)
            res[:, :, atom_label == item.label] = chi[:, :, numpy.newaxis]
        return res

    def view_structure(
            self, settings=None, interactive: bool = False, js: str = None,
            **kwargs):
        """Interactive 3D view of the crystal structure (py3Dmol).

        Needs the optional extra ``pip install cryspy[view3d]``. In a
        Jupyter notebook the returned view shows itself; with
        ``interactive=True`` sliders and check boxes for the settings are
        shown above it.

        Parameters
        ----------
        settings : StructureViewSettings
            Settings of the view; the defaults if None.
        interactive : bool
            Show widgets to change the settings (Jupyter, ipywidgets).
        js : str
            URL of 3Dmol-min.js; the py3Dmol default (CDN) if None.
        **kwargs
            Fields of ``StructureViewSettings`` to change, e.g.
            ``n_cells_a=2, show_labels=True``.

        Returns
        -------
        py3Dmol.view, or ipywidgets.VBox when ``interactive`` is True.
        """
        if settings is None:
            settings = StructureViewSettings()
        settings = dataclasses.replace(settings, **kwargs)
        if not interactive:
            scene = self.get_structure_scene(settings)
            return make_structure_view(
                scene, scene["atom_element"], scene["atom_label"],
                settings=settings, js=js)

        widgets = make_settings_widgets(settings)
        output = ipywidgets.Output()

        def redraw(change=None):
            new = StructureViewSettings(
                **{name: widget.value for name, widget in widgets.items()})
            output.clear_output(wait=True)
            with output:
                self.view_structure(settings=new, js=js).show()

        for widget in widgets.values():
            widget.observe(redraw, names="value")
        redraw()
        return ipywidgets.VBox([
            ipywidgets.HBox(list(widgets.values()),
                            layout=ipywidgets.Layout(flex_flow="row wrap")),
            output])

    def get_structure_view_html(self, settings=None, js: str = None) -> str:
        """HTML page with the 3D view of the structure (for a web view).

        Parameters
        ----------
        settings : StructureViewSettings
            Settings of the view; the defaults if None.
        js : str
            URL of 3Dmol-min.js, e.g. a local file for work offline; the
            py3Dmol default (CDN) if None.
        """
        view = self.view_structure(settings=settings, js=js)
        return make_structure_view_html(view)

    def get_flags_atom_beta(self):
        atom_site = self.atom_site
        if self.is_attribute("atom_site_aniso"):
            atom_site_aniso = self.atom_site_aniso
            l_f = []
            for item_as in atom_site.items:
                if item_as.adp_type == "Bani":
                    item_asa = atom_site_aniso[item_as.label]
                    l_f.append(
                        numpy.array(
                            [
                                item_asa.b_11_refinement,
                                item_asa.b_22_refinement,
                                item_asa.b_33_refinement,
                                item_asa.b_12_refinement,
                                item_asa.b_13_refinement,
                                item_asa.b_23_refinement,
                            ],
                            dtype=bool,
                        )
                    )
                elif item_as.adp_type == "Uani":
                    item_asa = atom_site_aniso[item_as.label]
                    l_f.append(
                        numpy.array(
                            [
                                item_asa.u_11_refinement,
                                item_asa.u_22_refinement,
                                item_asa.u_33_refinement,
                                item_asa.u_12_refinement,
                                item_asa.u_13_refinement,
                                item_asa.u_23_refinement,
                            ],
                            dtype=bool,
                        )
                    )
                else:
                    l_f.append(numpy.zeros((6,), dtype=bool))
            res = numpy.stack(l_f, axis=0).transpose()
        else:
            res = numpy.zeros((6, len(atom_site.items)), dtype=bool)
        return res

    def get_dictionary(self):
        """Form dictionary. See documentation moduel CrysPy using Jupyter notebook."""
        self.form_object()
        ddict = super(Crystal, self).get_dictionary()
        return form_crystal_dictionary(self, ddict)

    def take_parameters_from_dictionary(
        self,
        ddict_crystal,
        l_parameter_name: list = None,
        l_sigma: list = None,
    ):
        return take_crystal_parameters_from_dictionary(
            self, ddict_crystal=ddict_crystal, l_parameter_name=l_parameter_name, l_sigma=l_sigma)

    def save_atom_rho_multipole_density_to_den(
        self,
        file_den: str = "file.den",
        Na: int = 48,
        Nb: int = 48,
        Nc: int = 48,
    ):
        """Save the multipole electron density of the atoms to a .den file.

        The density is calculated on a grid of the unit cell in the
        Cartesian coordinate system with Z along the axis c and X along the
        reciprocal axis a*.

        Parameters
        ----------
        file_den : str, optional
            The file to write.
        Na, Nb, Nc : int, optional
            Numbers of points along the axes a, b, c.
        """
        return save_crystal_multipole_density_to_den(
            self, file_den=file_den, Na=Na, Nb=Nb, Nc=Nc)

    def save_atom_rho_multipole_magnetization_density_to_den(
        self,
        file_den: str = "file.den",
        Na: int = 48,
        Nb: int = 48,
        Nc: int = 48,
    ):
        """Save the multipole magnetization density (mx, my, mz) to a .den file.

        The density is calculated on a grid of the unit cell in the
        Cartesian coordinate system with Z along the axis c and X along the
        reciprocal axis a*, from the ordered moments of the atoms.

        Parameters
        ----------
        file_den : str, optional
            The file to write.
        Na, Nb, Nc : int, optional
            Numbers of points along the axes a, b, c.
        """
        return save_crystal_multipole_magnetization_density_to_den(
            self, file_den=file_den, Na=Na, Nb=Nb, Nc=Nc)
