"""Description of Crystal class."""
__author__ = 'ikibalin'
__version__ = "2020_08_19"
from typing import NoReturn
import numpy
from warnings import warn

from cryspy.A_functions_base.function_3_mcif import calc_full_sym_elems, \
    calc_hkl_in_range

from cryspy.B_parent_classes.cl_3_data import DataN

from cryspy.C_item_loop_classes.cl_1_cell import Cell
from cryspy.C_item_loop_classes.cl_1_atom_site import AtomSiteL
from cryspy.C_item_loop_classes.cl_1_atom_type import AtomTypeL
from cryspy.C_item_loop_classes.cl_1_atom_site_aniso import \
    AtomSiteAnisoL
from cryspy.C_item_loop_classes.cl_1_refln import ReflnL
from cryspy.C_item_loop_classes.cl_1_atom_site_moment import AtomSiteMomentL
from cryspy.C_item_loop_classes.cl_1_space_group_symop_magn_centering import \
    SpaceGroupSymopMagnCenteringL
from cryspy.C_item_loop_classes.cl_2_space_group_symop_magn_operation import \
    SpaceGroupSymopMagnOperationL

from cryspy.C_item_loop_classes.cl_1_atom_type_scat import \
    AtomTypeScatL
from cryspy.C_item_loop_classes.cl_2_atom_site_scat import \
    AtomSiteScatL
from cryspy.C_item_loop_classes.cl_1_atom_local_axes import \
    AtomLocalAxesL
from cryspy.C_item_loop_classes.cl_1_atom_electron_configuration \
    import AtomElectronConfigurationL

from cryspy.D_functions_item_loop.function_1_calc_for_magcrystal import \
    calc_f_nucl, calc_f_mag
from cryspy.D_functions_item_loop.function_1_flip_ratio import \
    calc_fm_perp_for_fm_loc


class MagCrystal(DataN):
    """
    Crystal structure description with defined magnetic space group.

    Data items in the MAG_CRYSTAL category record details about
    magnetic crystal structure.

    """

    CLASSES_MANDATORY = (Cell, AtomSiteL, AtomSiteMomentL,
                         SpaceGroupSymopMagnOperationL)
    CLASSES_OPTIONAL = (
        AtomTypeL, AtomSiteAnisoL, SpaceGroupSymopMagnCenteringL,
        AtomSiteScatL, AtomTypeScatL, AtomLocalAxesL,
        AtomElectronConfigurationL)
    # CLASSES_INTERNAL = ()

    CLASSES = CLASSES_MANDATORY + CLASSES_OPTIONAL

    PREFIX = "mag_crystal"

    # default values for the parameters
    D_DEFAULT = {}

    def __init__(self, data_name=None, **kwargs) -> NoReturn:
        super(MagCrystal, self).__init__()

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
        try:
            atom_site = self.atom_site
            atom_site_scat = self.atom_site_scat
            atom_site_scat.load_atom_type_scat_by_atom_site(atom_site)
        except AttributeError:
            pass

        # space_group = self.space_group
        # space_group_wyckoff = space_group.space_group_wyckoff

        # cell = self.cell
        # cell.type_cell = space_group.bravais_type
        # cell.it_coordinate_system_code = space_group.it_coordinate_system_code
        # cell.apply_constraints()

        # atom_site = self.atom_site
        # atom_site.apply_constraints(space_group_wyckoff)
        # atom_site_aniso = self.atom_site_aniso
        # if atom_site_aniso is not None:
        #     atom_site_aniso.apply_space_group_constraint(atom_site,
        #                                                  space_group)
        # atom_site_susceptibility = self.atom_site_susceptibility
        # if atom_site_susceptibility is not None:
        #     atom_site_susceptibility.apply_chi_iso_constraint(cell)
        #     atom_site_susceptibility.apply_moment_iso_constraint(cell)
        #     atom_site_susceptibility.apply_space_group_constraint(atom_site,
        #                                                           space_group)

    def calc_b_iso_beta(self):
        """
        Calculate b_iso and beta_ij based on atom_site and atom_sites.

        For each atom defined in atom_site.
        """
        pass

    def calc_f_nucl(self, index_h, index_k, index_l):
        """
        Calculate nuclear structure factor.

        Arguments
        ---------
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
        cell = self.cell
        atom_site = self.atom_site
        try:
            atom_site_aniso = self.atom_site_aniso
        except AttributeError:
            atom_site_aniso = None

        space_group_symop_magn_operation = \
            self.space_group_symop_magn_operation
        space_group_symop_magn_centering = \
            self.space_group_symop_magn_centering

        f_hkl_as, dder = calc_f_nucl(
            index_h, index_k, index_l,
            space_group_symop_magn_operation,
            space_group_symop_magn_centering,
            cell, atom_site, atom_site_aniso, flag_derivatives=False)
        return f_hkl_as

    def calc_refln(self, index_h, index_k, index_l, flag_internal=True):
        """
        Calculate Refln cryspy object where nuclear structure factor is stored.

        Arguments
        ---------
            h, k, l: 1D numpy array of Miller indexes
            flag_internal: a flag to calculate or to use internal objects.
                           It should be True if user call the function.
                           It's True by default.

        Output
        ------
            refln: object cryspy.Refln

        """
        f_nucl = self.calc_f_nucl(index_h, index_k, index_l)
        res = ReflnL()
        res.numpy_index_h = index_h
        res.numpy_index_k = index_k
        res.numpy_index_l = index_l
        res.numpy_f_calc = f_nucl
        if flag_internal:
            res.numpy_to_items()
        return res

    def calc_f_mag(self, index_h, index_k, index_l):
        """
        Calculate magnetic structure factor.

        Keyword Arguments:
        -----------------
            index_h, index_k, index_l: 1D numpy array of Miller indexes

        Output:
        -------
            f_mag: 2D numpy array of Nuclear structure factor

        Example:
        -------
            >>> import numpy as np
            >>> h, k, l = np.array([1,2],dtype=int), np.array([1,0],dtype=int),
                          np.array([1,0],dtype=int)
            >>> f_nucl = crystal.calc_f_nucl(h, k, l)

        """
        cell = self.cell
        atom_site = self.atom_site
        try:
            atom_site_aniso = self.atom_site_aniso
        except AttributeError:
            atom_site_aniso = None

        atom_site_scat = self.atom_site_scat
        atom_site_moment = self.atom_site_moment

        space_group_symop_magn_operation = \
            self.space_group_symop_magn_operation
        space_group_symop_magn_centering = \
            self.space_group_symop_magn_centering

        flag_only_orbital = False
        flag_derivatives = False
        f_mag, dder = calc_f_mag(
            (index_h, index_k, index_l), space_group_symop_magn_operation,
            space_group_symop_magn_centering, cell, atom_site, atom_site_aniso,
            atom_site_scat, atom_site_moment,
            flag_derivatives=flag_derivatives,
            flag_only_orbital=flag_only_orbital)
        return f_mag

    def calc_f_mag_perp(self, index_h, index_k, index_l):
        """
        Calculate perpendicular component of magnetic structure factor.

        Keyword Arguments:
        -----------------
            index_h, index_k, index_l: 1D numpy array of Miller indexes

        Output:
        -------
            f_mag: 2D numpy array of Nuclear structure factor

        Example:
        -------
            >>> import numpy as np
            >>> h, k, l = np.array([1,2],dtype=int), np.array([1,0],dtype=int),
                          np.array([1,0],dtype=int)
            >>> f_nucl = crystal.calc_f_nucl(h, k, l)

        """
        f_mag = self.calc_f_mag(index_h, index_k, index_l)

        cell = self.cell
        k_loc = cell.calc_k_loc(index_h, index_k, index_l)
        f_mag_perp = calc_fm_perp_for_fm_loc(k_loc, f_mag)
        f_mag_perp = numpy.array(f_mag_perp, dtype=complex)
        return f_mag_perp

    def calc_hkl(self, sthovl_min: float, sthovl_max: float):
        """
        Give hkl taking symmetry into account.

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
        rad = numpy.pi / 180.
        unit_cell_parameters = (
            cell.length_a, cell.length_b, cell.length_c,
            cell.angle_alpha*rad, cell.angle_beta*rad, cell.angle_gamma*rad)

        sym_elems = self.space_group_symop_magn_operation.get_sym_elems()
        magn_centering = self.space_group_symop_magn_centering.get_sym_elems()
        full_sym_elems = calc_full_sym_elems(sym_elems, magn_centering)

        ind_hkl_mult = calc_hkl_in_range(
            full_sym_elems, unit_cell_parameters, sthovl_min, sthovl_max)

        return ind_hkl_mult


    def calc_hkl_in_range(self, sthovl_min: float, sthovl_max: float):
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
        rad = numpy.pi / 180.
        unit_cell_parameters = (
            cell.length_a, cell.length_b, cell.length_c,
            cell.angle_alpha*rad, cell.angle_beta*rad, cell.angle_gamma*rad)

        full_sym_elems = numpy.array([
            [0], [0], [0], [1], [1], [0], [0], [0], [1], [0], [0], [0], [1],
            [1], [0], [0], [0], [1], [0], [0], [0], [1]], dtype=int)

        ind_hkl_mult = calc_hkl_in_range(
            full_sym_elems, unit_cell_parameters, sthovl_min, sthovl_max)
        return ind_hkl_mult


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
