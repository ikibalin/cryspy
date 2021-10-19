"""DensityPoint and DensityPointL classes."""
from typing import NoReturn
from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN

import numpy

from cryspy.A_functions_base.function_1_matrices import calc_mRmCmRT, \
    calc_phase_3d, calc_moment_2d_by_susceptibility
from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary
from cryspy.A_functions_base.function_1_strings import transform_string_to_r_b

from cryspy.A_functions_base.function_2_mem import \
    calc_asymmetric_unit_cell_indexes, calc_factor_in_front_of_density_for_fm,\
    calc_index_atom_symmetry_closest_to_fract_xyz, calc_moment_perp

from cryspy.C_item_loop_classes.cl_1_cell import Cell
from cryspy.C_item_loop_classes.cl_1_space_group_symop import SpaceGroupSymopL
from cryspy.C_item_loop_classes.cl_1_atom_site import AtomSiteL
from cryspy.C_item_loop_classes.cl_1_atom_site_susceptibility import \
    AtomSiteSusceptibilityL
from cryspy.C_item_loop_classes.cl_1_atom_electron_configuration import \
    AtomElectronConfigurationL
from cryspy.C_item_loop_classes.cl_1_mem_parameters import MEMParameters
from cryspy.C_item_loop_classes.cl_2_space_group import SpaceGroup
from cryspy.C_item_loop_classes.cl_2_atom_rho_orbital_radial_slater import \
    AtomRhoOrbitalRadialSlaterL


class DensityPoint(ItemN):
    """
    Density point description.

    Attributes
    ----------
        - index_x, index_y, index_z (mandatory)
        - density, density_ferro, density_antiferro, multiplicity, fract_x,
          fract_y, fract_z, basin_atom_label, basin_atom_symop,
          basin_atom_distance (optional)

    Methods
    -------
        - calc_multiplicity
    """

    ATTR_MANDATORY_NAMES = ("index_x", "index_y", "index_z")
    ATTR_MANDATORY_TYPES = (int, int, int)
    ATTR_MANDATORY_CIF = ("index_x", "index_y", "index_z")

    ATTR_OPTIONAL_NAMES = ("density", "density_ferro", "density_antiferro",
                           "multiplicity", "fract_x", "fract_y", "fract_z",
                           "basin_atom_label", "basin_atom_symop",
                           "basin_atom_distance")
    ATTR_OPTIONAL_TYPES = (float, float, float, int, float, float, float,
                           str, str, float)
    ATTR_OPTIONAL_CIF = ("density", "density_ferro", "density_antiferro",
                         "multiplicity", "fract_x", "fract_y", "fract_z",
                         "basin_atom_label", "basin_atom_symop",
                         "basin_atom_distance")

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
    ATTR_CONSTR_MARK = tuple([f"{_h:}_mark" for _h in ATTR_REF])

    # formats if cif format
    D_FORMATS = {"density": "{:.5f}", "density_ferro": "{:.5f}",
                 "density_antiferro": "{:.5f}", "fract_x": "{:.5f}",
                 "fract_y": "{:.5f}", "fract_z": "{:.5f}",
                 "basin_atom_distance": "{:.5f}"}

    # constraints on the parameters
    D_CONSTRAINTS = {}

    # default values for the parameters
    D_DEFAULT = {}
    for key in ATTR_SIGMA:
        D_DEFAULT[key] = 0.
    for key in (ATTR_CONSTR_FLAG + ATTR_REF_FLAG):
        D_DEFAULT[key] = False
    for key in ATTR_CONSTR_MARK:
        D_DEFAULT[key] = ""

    PREFIX = "density_point"

    def __init__(self, **kwargs) -> NoReturn:
        super(DensityPoint, self).__init__()

        # defined for any integer and float parameters
        D_MIN = {"density": 0., "density_ferro": 0., "density_antiferro": 0.,
                 "multiplicity": 0, "basin_atom_distance": 0.}

        # defined for ani integer and float parameters
        D_MAX = {}

        self.__dict__["D_MIN"] = D_MIN
        self.__dict__["D_MAX"] = D_MAX
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)

    def calc_multiplicity(self, space_group: SpaceGroup):
        """Calculate the multiplicity of the point."""
        fract_x, fract_y, fract_z = self.fract_x, self.fract_y, self.fract_z
        mult = space_group.calc_xyz_mult(fract_x, fract_y, fract_z)[-1]
        # check of multiplicity
        el_2 = space_group.calc_el_symm_for_xyz(fract_x, fract_y,
                                                fract_z)[0].size
        mult_2 = int(len(space_group.full_space_group_symop.items)/el_2)
        if mult != mult_2:
            raise UserWarning(
                f"{fract_x:.5f} {fract_y:.5f} {fract_z:.5f} {mult:} {mult_2:}")
        self.multiplicity = mult


class DensityPointL(LoopN):
    """
    Description of Loop in loop.

    Internal Arguments
    ------------------
        - number_symmetry
        - number_unit_cell
        - volume_unit_cell

    Methods
    -------
        - form_asymmetric_unit_cell
    """

    ITEM_CLASS = DensityPoint
    ATTR_INDEX = None

    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(DensityPointL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name

    def form_asymmetric_unit_cell(self, space_group_symop: SpaceGroupSymopL,
                                  points_a: int = 48, points_b: int = 48,
                                  points_c: int = 48):
        """
        Form asymmetric unit cell.

        Parameters
        ----------
        space_group_symop : SpaceGroupSymop of cryspy library
            Full list of element symmetry.
        points_a : int, optional
            Number of points along axis a. The default is 48.
        points_b : int, optional
            Number of points along axis b. The default is 48.
        points_c : int, optional
            Number of points along axis c. The default is 48.

        Returns
        -------
        None.

        """
        r_11 = numpy.array(space_group_symop.r_11, dtype=float)
        r_12 = numpy.array(space_group_symop.r_12, dtype=float)
        r_13 = numpy.array(space_group_symop.r_13, dtype=float)

        r_21 = numpy.array(space_group_symop.r_21, dtype=float)
        r_22 = numpy.array(space_group_symop.r_22, dtype=float)
        r_23 = numpy.array(space_group_symop.r_23, dtype=float)

        r_31 = numpy.array(space_group_symop.r_31, dtype=float)
        r_32 = numpy.array(space_group_symop.r_32, dtype=float)
        r_33 = numpy.array(space_group_symop.r_33, dtype=float)

        b_1 = numpy.array(space_group_symop.b_1, dtype=float)
        b_2 = numpy.array(space_group_symop.b_2, dtype=float)
        b_3 = numpy.array(space_group_symop.b_3, dtype=float)

        r_ij = (r_11, r_12, r_13, r_21, r_22, r_23, r_31, r_32, r_33)
        b_i = (b_1, b_2, b_3)

        ind_x_a_u_c, ind_y_a_u_c, ind_z_a_u_c, counts_a_u_c = \
            calc_asymmetric_unit_cell_indexes(points_a, points_b, points_c,
                                              r_ij, b_i)

        fract_x_a_u_c = ind_x_a_u_c.astype(float)/float(points_a)
        fract_y_a_u_c = ind_y_a_u_c.astype(float)/float(points_b)
        fract_z_a_u_c = ind_z_a_u_c.astype(float)/float(points_c)

        l_item = []
        for _ix, _iy, _iz, _mult, _fx, _fy, _fz in zip(
                ind_x_a_u_c, ind_y_a_u_c, ind_z_a_u_c, counts_a_u_c,
                fract_x_a_u_c, fract_y_a_u_c, fract_z_a_u_c):

            item = DensityPoint(index_x=_ix, index_y=_iy, index_z=_iz,
                                multiplicity=_mult,
                                fract_x=_fx, fract_y=_fy, fract_z=_fz)
            l_item.append(item)
        self.items = l_item

        number_symmetry = r_11.size
        number_unit_cell = points_a*points_b*points_c

        self.number_symmetry = number_symmetry
        self.number_unit_cell = number_unit_cell

        self.calc_rbs_i(space_group_symop, points_a=points_a,
                        points_b=points_b, points_c=points_c)

    def separate_basins(self, space_group_symop: SpaceGroupSymopL, cell: Cell,
                        atom_site: AtomSiteL, points_a: int = 48,
                        points_b: int = 48, points_c: int = 48):
        """
        Separate unit cell on basins.

        Parameters
        ----------
        space_group_symop : SpaceGroupSymopL of cryspy library
            DESCRIPTION.
        cell : Cell of cryspy library
            DESCRIPTION.
        atom_site : AtomSiteL of cryspy library
            DESCRIPTION.
        points_a : int, optional
            DESCRIPTION. The default is 48.
        points_b : int, optional
            DESCRIPTION. The default is 48.
        points_c : int, optional
            DESCRIPTION. The default is 48.

        Returns
        -------
        None.

        """
        cell.form_object
        self.volume_unit_cell = cell.volume

        r_11 = numpy.array(space_group_symop.r_11, dtype=float)
        r_12 = numpy.array(space_group_symop.r_12, dtype=float)
        r_13 = numpy.array(space_group_symop.r_13, dtype=float)

        r_21 = numpy.array(space_group_symop.r_21, dtype=float)
        r_22 = numpy.array(space_group_symop.r_22, dtype=float)
        r_23 = numpy.array(space_group_symop.r_23, dtype=float)

        r_31 = numpy.array(space_group_symop.r_31, dtype=float)
        r_32 = numpy.array(space_group_symop.r_32, dtype=float)
        r_33 = numpy.array(space_group_symop.r_33, dtype=float)

        b_1 = numpy.array(space_group_symop.b_1, dtype=float)
        b_2 = numpy.array(space_group_symop.b_2, dtype=float)
        b_3 = numpy.array(space_group_symop.b_3, dtype=float)

        r_ij = (r_11, r_12, r_13, r_21, r_22, r_23, r_31, r_32, r_33)
        b_i = (b_1, b_2, b_3)

        fract_x = self.numpy_fract_x
        fract_y = self.numpy_fract_y
        fract_z = self.numpy_fract_z

        fract_xyz = numpy.transpose(numpy.vstack([fract_x, fract_y, fract_z]))

        fract_atom_x = numpy.array(atom_site.fract_x, dtype=float)
        fract_atom_y = numpy.array(atom_site.fract_y, dtype=float)
        fract_atom_z = numpy.array(atom_site.fract_z, dtype=float)

        fract_atom_xyz = numpy.transpose(numpy.vstack([
            fract_atom_x, fract_atom_y, fract_atom_z]))

        ind_at, ind_sym, numpy_basin_atom_distance = \
            calc_index_atom_symmetry_closest_to_fract_xyz(
                fract_xyz, fract_atom_xyz, r_ij, b_i, cell)

        numpy_atom_label = numpy.array(atom_site.label, str)
        numpy_basin_atom_label = numpy_atom_label[ind_at]

        numpy_symop_operation_xyz = numpy.array(
            space_group_symop.operation_xyz)
        numpy_basin_atom_symop = numpy_symop_operation_xyz[ind_sym]

        for item, basin_atom_label, basin_atom_symop, basin_atom_distance in \
            zip(self.items, numpy_basin_atom_label, numpy_basin_atom_symop,
                numpy_basin_atom_distance):
            item.basin_atom_label = str(basin_atom_label)
            item.basin_atom_symop = str(basin_atom_symop)
            item.basin_atom_distance = float(basin_atom_distance)

        numpy_structured_atom = numpy.array(
            list(zip(atom_site.label, atom_site.multiplicity)),
            dtype=[("label", "U10"), ("multiplicity", "i4")])

        self.numpy_structured_atom = numpy_structured_atom
        self.numpy_basin_atom_label = numpy_basin_atom_label
        self.numpy_basin_atom_symop = numpy_basin_atom_symop
        self.numpy_basin_atom_distance = numpy_basin_atom_distance

    def set_flat_density(self, l_magnetic_labes: list,
                         flag_two_channel: bool = False):
        r"""
        Define starting density as flat density by expression.

        p_{i} = \frac{N_{points} \Sum_{a} m_{a}} {V_{u.c.} \Sum_{i} m_{i}}
        """
        numpy_structured_atom = self.numpy_structured_atom
        n_u_c = self.number_unit_cell
        numpy_multiplicity = self.numpy_multiplicity

        volume_u_c = self.volume_unit_cell
        numpy_density = numpy.zeros(self.numpy_index_x.shape)
        numpy_basin_atom_label = self.numpy_basin_atom_label
        for structured_atom in numpy_structured_atom:
            mult_a = int(structured_atom["multiplicity"])
            label_a = str(structured_atom["label"])
            flag_atom = numpy_basin_atom_label == label_a
            sum_m_i_basin = (numpy_multiplicity[flag_atom]).sum()
            if sum_m_i_basin != 0:  # Case of two atoms in the same position
                den_i = float(n_u_c * mult_a) / float(sum_m_i_basin *
                                                      volume_u_c)
                numpy_density[flag_atom] = den_i

        self.numpy_density = numpy_density

        numpy_label = numpy.array(l_magnetic_labes, dtype=str)
        flag_2d = numpy_basin_atom_label[:, numpy.newaxis] == numpy_label[
            numpy.newaxis, :]
        flag_atom = numpy.any(flag_2d, axis=1)
        flag_background = numpy.logical_not(flag_atom)
        if flag_two_channel:
            flag_background = numpy.ones(flag_background.shape, bool)

        numpy_multiplicity = self.numpy_multiplicity

        sum_m_i_back = numpy_multiplicity[flag_background].sum()
        den_i = float(n_u_c * 1.) / float(sum_m_i_back * volume_u_c)
        numpy_density_ferro = numpy.zeros(self.numpy_index_x.shape)
        numpy_density_ferro[flag_background] += den_i

        numpy_density_antiferro = numpy.zeros(self.numpy_index_x.shape)
        numpy_density_antiferro[flag_background] += den_i

        self.numpy_density_ferro = numpy_density_ferro
        self.numpy_density_antiferro = numpy_density_antiferro

        self.renormalize_numpy_densities(flag_two_channel=flag_two_channel)
        self.numpy_to_items()

    def set_core_density(
            self, atom_site: AtomSiteL,
            atom_electron_configuration: AtomElectronConfigurationL,
            flag_two_channel: bool = False):
        r"""
        Define starting density as flat density by expression.

        p_{i} = \frac{N_{points} \Sum_{a} m_{a}} {V_{u.c.} \Sum_{i} m_{i}}
        """
        n_u_c = self.number_unit_cell
        volume_u_c = self.volume_unit_cell
        numpy_basin_atom_label = self.numpy_basin_atom_label
        numpy_basin_atom_distance = self.numpy_basin_atom_distance
        numpy_density = numpy.zeros(shape=numpy_basin_atom_distance.shape,
                                    dtype=float)

        for item_a_e_c in atom_electron_configuration.items:
            at_label = item_a_e_c.label
            at_type = atom_site[at_label].type_symbol

            flag_at_label = numpy_basin_atom_label == at_label

            l_sh_pop = item_a_e_c.get_core_shells_populations()

            l_arors = AtomRhoOrbitalRadialSlaterL().take_objects_for_atom_type(
                at_type)

            for sh_pop in l_sh_pop:
                shell, population = sh_pop
                obj = None
                for arors in l_arors:
                    try:
                        l_val = getattr(arors, f"coeff_{shell:}")
                        obj = arors
                        break
                    except AttributeError:
                        pass
                if obj is not None:
                    rho_sq = (obj.calc_normalized_rho(
                                numpy_basin_atom_distance[flag_at_label],
                                shell, kappa=1.))**2
                    rho_sq_1d = numpy.reshape(rho_sq, (rho_sq.size, ))

                    numpy_density[flag_at_label] += float(population) / \
                        (4.*numpy.pi) * rho_sq_1d

        self.numpy_density = numpy_density

        numpy_label = numpy.array([_h.label for _h in
                                   atom_electron_configuration.items],
                                  dtype=str)
        flag_2d = numpy_basin_atom_label[:, numpy.newaxis] == numpy_label[
            numpy.newaxis, :]
        flag_atom = numpy.any(flag_2d, axis=1)
        flag_background = numpy.logical_not(flag_atom)

        numpy_multiplicity = self.numpy_multiplicity
        sum_m_i_back = numpy_multiplicity[flag_background].sum()

        den_i = float(n_u_c * 1.) / float(sum_m_i_back * volume_u_c)
        numpy_density_ferro = numpy.zeros(numpy_basin_atom_distance.shape)
        numpy_density_ferro[flag_background] += den_i

        numpy_density_antiferro = numpy.zeros(numpy_basin_atom_distance.shape)
        numpy_density_antiferro[flag_background] += den_i

        self.numpy_density_ferro = numpy_density_ferro
        self.numpy_density_antiferro = numpy_density_antiferro

        self.renormalize_numpy_densities(flag_two_channel=flag_two_channel)
        self.numpy_to_items()


    def create_flat_density(self, space_group_symop: SpaceGroupSymopL,
                            cell: Cell, atom_site: AtomSiteL,
                            l_magnetic_labes: list, points_a: int = 48,
                            points_b: int = 48, points_c: int = 48,
                            flag_two_channel: bool = False):
        r"""Define starting density as flat density by expression.

        p_{i} = \frac{N_{points} \Sum_{a} m_{a}} {V_{u.c.} \Sum_{i} m_{i}}
        """
        self.form_asymmetric_unit_cell(
            space_group_symop, points_a=points_a, points_b=points_b,
            points_c=points_c)

        self.separate_basins(
            space_group_symop, cell, atom_site, points_a=points_a,
            points_b=points_b, points_c=points_c)

        self.set_flat_density(l_magnetic_labes,
                              flag_two_channel=flag_two_channel)

    def create_core_density(
            self, space_group_symop: SpaceGroupSymopL, cell: Cell,
            atom_site: AtomSiteL, atom_electron_configuration:
            AtomElectronConfigurationL, points_a: int = 48, points_b: int = 48,
            points_c: int = 48, flag_two_channel: bool = False):
        """Create core density.

        Parameters
        ----------
        space_group_symop : SpaceGroupSymopL
            DESCRIPTION.
        cell : Cell
            DESCRIPTION.
        atom_site : AtomSiteL
            DESCRIPTION.
        atom_electron_configuration : AtomElectronConfigurationL
            DESCRIPTION.
        l_magnetic_labes : TYPE
            DESCRIPTION.
        points_a : TYPE, optional
            DESCRIPTION. The default is 48.
        points_b : TYPE, optional
            DESCRIPTION. The default is 48.
        points_c : TYPE, optional
            DESCRIPTION. The default is 48.

        Returns
        -------
        None.

        """
        self.form_asymmetric_unit_cell(
            space_group_symop, points_a=points_a, points_b=points_b,
            points_c=points_c)

        self.separate_basins(
            space_group_symop, cell, atom_site, points_a=points_a,
            points_b=points_b, points_c=points_c)

        self.set_core_density(atom_site, atom_electron_configuration,
                              flag_two_channel=flag_two_channel)

    def calc_density_chi(self, cell: Cell,
                         atom_site_susceptibility: AtomSiteSusceptibilityL):
        """Calculate density chi."""
        density = self.density

        if atom_site_susceptibility is None:
            return numpy.zeros(len(density), dtype=float)

        label = self.basin_atom_label

        l_max_moment = [numpy.abs(
            item.calc_main_axes_of_magnetization_ellipsoid(cell)[0]).max()
            for item in atom_site_susceptibility.items]
        l_label_atom = atom_site_susceptibility.label

        l_den_chi = [den*l_max_moment[l_label_atom.index(lab)] if lab in
                     l_label_atom else 0. for lab, den in zip(label, density)]

        return numpy.array(l_den_chi, dtype=float)

    def save_to_file_den(
            self, mem_parameters: MEMParameters, space_group: SpaceGroup,
            cell: Cell, f_name: str = "file.den",
            atom_site_susceptibility: AtomSiteSusceptibilityL = None,
            f_background: str = "file_back.den"):
        """Save to file."""
        numpy_basin_atom_label = numpy.array(self.basin_atom_label, dtype=str)
        points_a = mem_parameters.points_a
        points_b = mem_parameters.points_b
        points_c = mem_parameters.points_c
        chi_iso_ferro = mem_parameters.chi_ferro
        chi_iso_antiferro = mem_parameters.chi_antiferro
        index_x = self.index_x
        index_y = self.index_y
        index_z = self.index_z
        # density = self.density
        density_chi = self.calc_density_chi(cell, atom_site_susceptibility)
        density_ferro = self.density_ferro
        density_antiferro = self.density_antiferro
        ls_out = []
        ls_out_b = ["Created by CrysPy, mu_B/Tesla"]
        ls_out.append("Created by CrysPy")
        ls_out.append("{:}".format(len(index_x)))
        ls_out_b.append("{:}".format(len(index_x)))

        for _x, _y, _z, _l, den_chi, den_f, den_a in \
            zip(index_x, index_y, index_z, numpy_basin_atom_label, density_chi,
                density_ferro, density_antiferro):
            ls_out.append(f"{_x:4} {_y:4} {_z:4} {den_chi:15.7f}")
            v1 = (chi_iso_ferro*den_f+chi_iso_antiferro*den_a)
            ls_out_b.append(f"{_x:4} {_y:4} {_z:4} {v1:15.7f}")

        a, b, c = cell.length_a, cell.length_b, cell.length_c
        al, be, ga = cell.angle_alpha, cell.angle_beta, cell.angle_gamma

        ls_out.append(
            f"{a:10.5f}{b:10.5f}{c:10.5f}{al:10.5f}{be:10.5f}{ga:10.5f}")
        ls_out_b.append(
            f"{a:10.5f}{b:10.5f}{c:10.5f}{al:10.5f}{be:10.5f}{ga:10.5f}")

        r_s_g_s = space_group.reduced_space_group_symop

        centr = int(space_group.centrosymmetry)
        l_orig = space_group.shift
        ls_out.append("{:5}{:5}{:5}{:5}{:5}{:5}".format(
            points_a, points_b, points_c, len(r_s_g_s.items), centr,
            len(l_orig)))
        ls_out_b.append("{:5}{:5}{:5}{:5}{:5}{:5}".format(
            points_a, points_b, points_c, len(r_s_g_s.items), centr,
            len(l_orig)))
        for _ in r_s_g_s.items:
            ls_out.append(
                f"{int(_.r_11):4}{int(_.r_21):4}{int(_.r_31):4}  \
{int(_.r_12):4}{int(_.r_22):4}{int(_.r_32):4}  {int(_.r_13):4}{int(_.r_23):4}\
{int(_.r_33):4}    {float(_.b_1):8.5f}{float(_.b_2):8.5f}{float(_.b_3):8.5f}")
            ls_out_b.append(
                f"{int(_.r_11):4}{int(_.r_21):4}{int(_.r_31):4}  \
{int(_.r_12):4}{int(_.r_22):4}{int(_.r_32):4}  {int(_.r_13):4}{int(_.r_23):4}\
{int(_.r_33):4}    {float(_.b_1):8.5f}{float(_.b_2):8.5f}{float(_.b_3):8.5f}")
        for orig in l_orig:
            ls_out.append(
                f"{float(orig[0]):8.4f}{float(orig[1]):8.4f}\
{float(orig[2]):8.4f}")
            ls_out_b.append(
                f"{float(orig[0]):8.4f}{float(orig[1]):8.4f}\
{float(orig[2]):8.4f}")

        with open(f_name, "w") as fid:
            fid.write("\n".join(ls_out))
        with open(f_background, "w") as fid:
            fid.write("\n".join(ls_out_b))
        return

    def calc_fm(self, factor_i, factor_ferro_i, factor_antiferro_i):
        """Calculate magnetic structure factor or its perpendicular component.

        Parameters
        ----------
        factor_i : TYPE [hkl, points]
            Factor is output of method calc_factor_in_front_of_density_for_fm
            for calculations of magnetic structure factor

            or

            it is output of method calc_factor_in_front_of_density_for_fm_perp
            for calculations of perpendicular component of magnetic structure
            factor.

        Returns
        -------
        msf_p_i : TYPE [hkl]
            Magnetic structure factor or perpendicular magnetic structure
            factor.
        delta_msf_i : TYPE [hkl, points]
            derivative over density

            d FM / d den

            or

            d FM_perp / d den.
        """
        np_den = self.numpy_density
        factor_x, factor_y, factor_z = factor_i

        f_m_perp_x = (factor_x * np_den[numpy.newaxis, :]).sum(axis=1)
        f_m_perp_y = (factor_y * np_den[numpy.newaxis, :]).sum(axis=1)
        f_m_perp_z = (factor_z * np_den[numpy.newaxis, :]).sum(axis=1)
        f_m_perp = (f_m_perp_x, f_m_perp_y, f_m_perp_z)

        delta_f_m_perp = (factor_x, factor_y, factor_z)

        np_den_f = self.numpy_density_ferro
        factor_f_x, factor_f_y, factor_f_z = factor_ferro_i

        f_m_perp_f_x = (factor_f_x * np_den_f[numpy.newaxis, :]).sum(axis=1)
        f_m_perp_f_y = (factor_f_y * np_den_f[numpy.newaxis, :]).sum(axis=1)
        f_m_perp_f_z = (factor_f_z * np_den_f[numpy.newaxis, :]).sum(axis=1)
        f_m_perp_f = (f_m_perp_f_x, f_m_perp_f_y, f_m_perp_f_z)

        delta_f_m_perp_f = (factor_f_x, factor_f_y, factor_f_z)

        np_den_a = self.numpy_density_antiferro

        factor_a_x, factor_a_y, factor_a_z = factor_antiferro_i

        f_m_perp_a_x = (factor_a_x * np_den_a[numpy.newaxis, :]).sum(axis=1)
        f_m_perp_a_y = (factor_a_y * np_den_a[numpy.newaxis, :]).sum(axis=1)
        f_m_perp_a_z = (factor_a_z * np_den_a[numpy.newaxis, :]).sum(axis=1)
        f_m_perp_a = (f_m_perp_a_x, f_m_perp_a_y, f_m_perp_a_z)

        delta_f_m_perp_a = (factor_a_x, factor_a_y, factor_a_z)

        f_m_p = (f_m_perp[0]+f_m_perp_f[0]+f_m_perp_a[0],
                 f_m_perp[1]+f_m_perp_f[1]+f_m_perp_a[1],
                 f_m_perp[2]+f_m_perp_f[2]+f_m_perp_a[2])
        return f_m_p, delta_f_m_perp, delta_f_m_perp_f, delta_f_m_perp_a

    def calc_phase_3d(self, hkl, space_group_symop: SpaceGroupSymopL):
        """Calculate phase 3d.

        Parameters
        ----------
        hkl : TYPE
            DESCRIPTION.
        space_group_symop : SpaceGroupSymopL
            DESCRIPTION.

        Returns
        -------
        phase_3d : TYPE
            DESCRIPTION.
        """
        fract_x = self.numpy_fract_x
        fract_y = self.numpy_fract_y
        fract_z = self.numpy_fract_z
        fract_xyz = numpy.transpose(numpy.vstack([fract_x, fract_y, fract_z]))

        r_11 = numpy.array(space_group_symop.r_11, dtype=float)
        r_12 = numpy.array(space_group_symop.r_12, dtype=float)
        r_13 = numpy.array(space_group_symop.r_13, dtype=float)

        r_21 = numpy.array(space_group_symop.r_21, dtype=float)
        r_22 = numpy.array(space_group_symop.r_22, dtype=float)
        r_23 = numpy.array(space_group_symop.r_23, dtype=float)

        r_31 = numpy.array(space_group_symop.r_31, dtype=float)
        r_32 = numpy.array(space_group_symop.r_32, dtype=float)
        r_33 = numpy.array(space_group_symop.r_33, dtype=float)

        b_1 = numpy.array(space_group_symop.b_1, dtype=float)
        b_2 = numpy.array(space_group_symop.b_2, dtype=float)
        b_3 = numpy.array(space_group_symop.b_3, dtype=float)

        r_ij = (r_11, r_12, r_13, r_21, r_22, r_23, r_31, r_32, r_33)
        b_i = (b_1, b_2, b_3)

        phase_3d = calc_phase_3d(hkl, r_ij, b_i, fract_xyz)
        return phase_3d

    def calc_rbs_i(self, space_group_symop: SpaceGroupSymopL,
                   points_a: int = 48, points_b: int = 48, points_c: int = 48):
        """
        Caclulate symmetry elements that transform point to itself.

        Parameters
        ----------
        space_group_symop : SpaceGroupSymopL
            DESCRIPTION.

        Returns
        -------
        None.

        """
        r_11 = numpy.array(space_group_symop.r_11, dtype=float)
        r_12 = numpy.array(space_group_symop.r_12, dtype=float)
        r_13 = numpy.array(space_group_symop.r_13, dtype=float)

        r_21 = numpy.array(space_group_symop.r_21, dtype=float)
        r_22 = numpy.array(space_group_symop.r_22, dtype=float)
        r_23 = numpy.array(space_group_symop.r_23, dtype=float)

        r_31 = numpy.array(space_group_symop.r_31, dtype=float)
        r_32 = numpy.array(space_group_symop.r_32, dtype=float)
        r_33 = numpy.array(space_group_symop.r_33, dtype=float)

        b_1 = numpy.array(space_group_symop.b_1, dtype=float)
        b_2 = numpy.array(space_group_symop.b_2, dtype=float)
        b_3 = numpy.array(space_group_symop.b_3, dtype=float)

        np_ind_x, np_ind_y = self.numpy_index_x, self.numpy_index_y
        np_ind_z = self.numpy_index_z
        if np_ind_x is None:
            np_ind_x = numpy.array(self.index_x, dtype=int)
            np_ind_y = numpy.array(self.index_y, dtype=int)
            np_ind_z = numpy.array(self.index_z, dtype=int)

        l_rs_i = []
        l_bs_i = []
        for i_x, i_y, i_z in zip(np_ind_x, np_ind_y, np_ind_z):
            np_x_n = numpy.mod((numpy.around((
                i_x*r_11 + i_y*r_12 + i_z*r_13 +
                points_a*b_1).astype(float), 0)).astype(int), points_a)
            np_y_n = numpy.mod((numpy.around((
                i_x*r_21 + i_y*r_22 + i_z*r_23 +
                points_b*b_2).astype(float), 0)).astype(int), points_b)
            np_z_n = numpy.mod((numpy.around((
                i_x*r_31 + i_y*r_32 + i_z*r_33 +
                points_c*b_3).astype(float), 0)).astype(int), points_c)
            flag_x = np_x_n == i_x
            flag_y = np_y_n == i_y
            flag_z = np_z_n == i_z
            flag_xy = numpy.logical_and(flag_x, flag_y)
            fl_xyz = numpy.logical_and(flag_xy, flag_z)
            r_11_c, r_12_c, r_13_c = r_11[fl_xyz], r_12[fl_xyz], r_13[fl_xyz]
            r_21_c, r_22_c, r_23_c = r_21[fl_xyz], r_22[fl_xyz], r_23[fl_xyz]
            r_31_c, r_32_c, r_33_c = r_31[fl_xyz], r_32[fl_xyz], r_33[fl_xyz]

            b_1_c, b_2_c, b_3_c = b_1[fl_xyz], b_2[fl_xyz], b_3[fl_xyz]
            l_rs_i.append((r_11_c, r_12_c, r_13_c, r_21_c, r_22_c, r_23_c,
                           r_31_c, r_32_c, r_33_c))
            l_bs_i.append((b_1_c, b_2_c, b_3_c))
        self.l_rs_i = l_rs_i
        self.l_bs_i = l_bs_i

    def calc_susc_i(self, atom_site_susceptibility: AtomSiteSusceptibilityL):
        """Calculate susceptibility of point i.

        Parameters
        ----------
        atom_site_susceptibility : AtomSiteSusceptibilityL
            DESCRIPTION.

        Returns
        -------
        susc_11 : TYPE
            DESCRIPTION.
        susc_22 : TYPE
            DESCRIPTION.
        susc_33 : TYPE
            DESCRIPTION.
        susc_12 : TYPE
            DESCRIPTION.
        susc_13 : TYPE
            DESCRIPTION.
        susc_23 : TYPE
            DESCRIPTION.
        """
        numpy_basin_atom_label = self.numpy_basin_atom_label
        numpy_basin_atom_symop = self.numpy_basin_atom_symop
        if numpy_basin_atom_label is None:
            numpy_basin_atom_label = numpy.array(self.basin_atom_label,
                                                 dtype=str)
            numpy_basin_atom_symop = numpy.array(self.basin_atom_symop,
                                                 dtype=str)
            self.numpy_basin_atom_label = numpy_basin_atom_label
            self.numpy_basin_atom_symop = numpy_basin_atom_symop

        np_label = numpy.array(atom_site_susceptibility.label, dtype=str)
        np_chi_11 = numpy.array(atom_site_susceptibility.chi_11, dtype=float)
        np_chi_22 = numpy.array(atom_site_susceptibility.chi_22, dtype=float)
        np_chi_33 = numpy.array(atom_site_susceptibility.chi_33, dtype=float)
        np_chi_12 = numpy.array(atom_site_susceptibility.chi_12, dtype=float)
        np_chi_13 = numpy.array(atom_site_susceptibility.chi_13, dtype=float)
        np_chi_23 = numpy.array(atom_site_susceptibility.chi_23, dtype=float)

        l_r = []
        for symop in numpy_basin_atom_symop:
            r, b = transform_string_to_r_b(symop)
            r = r.astype(float)
            l_r.append((r[0, 0], r[0, 1], r[0, 2], r[1, 0], r[1, 1], r[1, 2],
                        r[2, 0], r[2, 1], r[2, 2]))
        np_r = numpy.array(l_r, dtype=float)
        np_r_11, np_r_12, np_r_13 = np_r[:, 0], np_r[:, 1], np_r[:, 2]
        np_r_21, np_r_22, np_r_23 = np_r[:, 3], np_r[:, 4], np_r[:, 5]
        np_r_31, np_r_32, np_r_33 = np_r[:, 6], np_r[:, 7], np_r[:, 8]

        susc_11 = numpy.zeros(numpy_basin_atom_label.size, dtype=float)
        susc_22 = numpy.zeros(numpy_basin_atom_label.size, dtype=float)
        susc_33 = numpy.zeros(numpy_basin_atom_label.size, dtype=float)
        susc_12 = numpy.zeros(numpy_basin_atom_label.size, dtype=float)
        susc_13 = numpy.zeros(numpy_basin_atom_label.size, dtype=float)
        susc_23 = numpy.zeros(numpy_basin_atom_label.size, dtype=float)

        for _l, chi_11, chi_22, chi_33, chi_12, chi_13, chi_23 in \
                zip(np_label, np_chi_11, np_chi_22, np_chi_33, np_chi_12,
                    np_chi_13, np_chi_23):
            f_1 = numpy_basin_atom_label == _l
            chi_ij = (chi_11, chi_12, chi_13, chi_12, chi_22, chi_23,
                      chi_13, chi_23, chi_33)
            c_11, c_12, c_13, c_21, c_22, c_23, c_31, c_32, c_33 = \
                calc_mRmCmRT((np_r_11[f_1], np_r_12[f_1], np_r_13[f_1],
                              np_r_21[f_1], np_r_22[f_1], np_r_23[f_1],
                              np_r_31[f_1], np_r_32[f_1], np_r_33[f_1]),
                             chi_ij)

            susc_11[f_1] = c_11
            susc_22[f_1] = c_22
            susc_33[f_1] = c_33
            susc_12[f_1] = c_12
            susc_13[f_1] = c_13
            susc_23[f_1] = c_23
        l_rs_i = self.l_rs_i
        l_c_11, l_c_22, l_c_33, l_c_12, l_c_13, l_c_23 = [], [], [], [], [], []
        for s_11, s_22, s_33, s_12, s_13, s_23, rs_i in zip(
                susc_11, susc_22, susc_33, susc_12, susc_13, susc_23, l_rs_i):
            r_11, r_12, r_13, r_21, r_22, r_23, r_31, r_32, r_33 = rs_i
            c_11, c_12, c_13, c_21, c_22, c_23, c_31, c_32, c_33 = \
                calc_mRmCmRT(rs_i, (s_11, s_12, s_13,
                                    s_12, s_22, s_23,
                                    s_13, s_23, s_33))
            ns = c_11.size
            l_c_11.append(c_11.sum()/float(ns))
            l_c_22.append(c_22.sum()/float(ns))
            l_c_33.append(c_33.sum()/float(ns))
            l_c_12.append(c_12.sum()/float(ns))
            l_c_13.append(c_13.sum()/float(ns))
            l_c_23.append(c_23.sum()/float(ns))

        chi_11 = numpy.array(l_c_11, dtype=float)
        chi_22 = numpy.array(l_c_22, dtype=float)
        chi_33 = numpy.array(l_c_33, dtype=float)
        chi_12 = numpy.array(l_c_12, dtype=float)
        chi_13 = numpy.array(l_c_13, dtype=float)
        chi_23 = numpy.array(l_c_23, dtype=float)
        return chi_11, chi_22, chi_33, chi_12, chi_13, chi_23

    def calc_moment_2d(self, space_group_symop: SpaceGroupSymopL, cell: Cell,
                       atom_site_susceptibility: AtomSiteSusceptibilityL,
                       h_loc, chi_iso_ferro: float = 0.,
                       chi_iso_antiferro: float = 0.,
                       flag_two_channel: bool = False):
        """Calculate moments over points and symmetry elements.

        Parameters
        ----------
        space_group_symop : SpaceGroupSymopL
            DESCRIPTION.
        cell : Cell
            DESCRIPTION.
        atom_site_susceptibility : AtomSiteSusceptibilityL
            DESCRIPTION.
        h_loc : TYPE
            DESCRIPTION.
        chi_iso_ferro: TYPE
            Should be positive. Default is 0 mu_B per T per unit cell
        chi_iso_antiferro: TYPE
            Should be negative. Default is 0 mu_B per T per unit cell
        flag_back_all:bool
            True if background is defined for all points and
            False if background is defined for points where susceptibility is
            not defined

        Returns
        -------
        moment_2d : TYPE
            DESCRIPTION.
        moment_2d_ferro : TYPE
            DESCRIPTION.
        moment_2d_antiferro : TYPE
            DESCRIPTION.
        """
        cell.form_object
        m_norm = cell.m_m_norm
        m_norm_ij = (m_norm[0, 0], m_norm[0, 1], m_norm[0, 2],
                     m_norm[1, 0], m_norm[1, 1], m_norm[1, 2],
                     m_norm[2, 0], m_norm[2, 1], m_norm[2, 2])

        r_11 = numpy.array(space_group_symop.r_11, dtype=float)
        r_12 = numpy.array(space_group_symop.r_12, dtype=float)
        r_13 = numpy.array(space_group_symop.r_13, dtype=float)

        r_21 = numpy.array(space_group_symop.r_21, dtype=float)
        r_22 = numpy.array(space_group_symop.r_22, dtype=float)
        r_23 = numpy.array(space_group_symop.r_23, dtype=float)

        r_31 = numpy.array(space_group_symop.r_31, dtype=float)
        r_32 = numpy.array(space_group_symop.r_32, dtype=float)
        r_33 = numpy.array(space_group_symop.r_33, dtype=float)

        r_ij = (r_11, r_12, r_13, r_21, r_22, r_23, r_31, r_32, r_33)

        if not(flag_two_channel):
            s_11, s_22, s_33, s_12, s_13, s_23 = self.calc_susc_i(
                atom_site_susceptibility)

            susc_i = (0.2695*s_11, 0.2695*s_22, 0.2695*s_33,
                      0.2695*s_12, 0.2695*s_13, 0.2695*s_23)

            moment_2d = calc_moment_2d_by_susceptibility(
                r_ij, susc_i, m_norm_ij, h_loc)
        else:
            moment_2d = (
                numpy.zeros(shape=(self.numpy_index_x.size, r_ij[0].size)),
                numpy.zeros(shape=(self.numpy_index_x.size, r_ij[0].size)),
                numpy.zeros(shape=(self.numpy_index_x.size, r_ij[0].size)))

        moment_ferro = (0.2695*chi_iso_ferro*h_loc[0]*numpy.ones(
                            shape=(self.numpy_index_x.size, r_ij[0].size)),
                        0.2695*chi_iso_ferro*h_loc[1]*numpy.ones(
                            shape=(self.numpy_index_x.size, r_ij[0].size)),
                        0.2695*chi_iso_ferro*h_loc[2]*numpy.ones(
                            shape=(self.numpy_index_x.size, r_ij[0].size)))
        moment_antiferro = (0.2695*chi_iso_antiferro*h_loc[0]*numpy.ones(
                                shape=(self.numpy_index_x.size, r_ij[0].size)),
                            0.2695*chi_iso_antiferro*h_loc[1]*numpy.ones(
                                shape=(self.numpy_index_x.size, r_ij[0].size)),
                            0.2695*chi_iso_antiferro*h_loc[2]*numpy.ones(
                                shape=(self.numpy_index_x.size, r_ij[0].size)))

        return moment_2d, moment_ferro, moment_antiferro

    def calc_factor_in_front_of_density_for_fm(
            self, hkl, space_group_symop: SpaceGroupSymopL, cell: Cell,
            atom_site_susceptibility: AtomSiteSusceptibilityL, h_loc,
            chi_iso_ferro: float = 0., chi_iso_antiferro: float = 0.,
            flag_two_channel: bool = False):
        """Calculate factor in front of density for FM.

        Parameters
        ----------
        hkl : TYPE
            DESCRIPTION.
        space_group_symop : SpaceGroupSymopL
            DESCRIPTION.
        cell : Cell
            DESCRIPTION.
        atom_site_susceptibility : AtomSiteSusceptibilityL
            DESCRIPTION.
        h_loc : TYPE
            DESCRIPTION.
        chi_iso_ferro : float, optional
            DESCRIPTION. The default is 0..
        chi_iso_antiferro : float, optional
            DESCRIPTION. The default is 0..

        Returns
        -------
        v_hkl_2d_i : TYPE
            DESCRIPTION.
        v_ferro : TYPE
            DESCRIPTION.
        v_antiferro : TYPE
            DESCRIPTION.

        """
        mult_i = self.numpy_multiplicity
        np = self.number_unit_cell
        volume = self.volume_unit_cell
        moment_2d, moment_ferro, moment_antiferro = \
            self.calc_moment_2d(space_group_symop, cell,
                                atom_site_susceptibility, h_loc,
                                chi_iso_ferro=chi_iso_ferro,
                                chi_iso_antiferro=chi_iso_antiferro,
                                flag_two_channel=flag_two_channel)
        phase_3d = self.calc_phase_3d(hkl, space_group_symop)
        v_hkl_2d_i = \
            calc_factor_in_front_of_density_for_fm(mult_i, np, volume,
                                                   moment_2d, phase_3d, axis=2)
        v_ferro = \
            calc_factor_in_front_of_density_for_fm(
                mult_i, np, volume, moment_ferro, phase_3d, axis=2)

        v_antiferro = \
            calc_factor_in_front_of_density_for_fm(
                mult_i, np, volume, moment_antiferro, phase_3d, axis=2)

        return v_hkl_2d_i, v_ferro, v_antiferro

    def calc_factor_in_front_of_density_for_fm_perp(
            self, hkl, space_group_symop: SpaceGroupSymopL, cell: Cell,
            atom_site_susceptibility: AtomSiteSusceptibilityL, h_loc,
            chi_iso_ferro: float = 0., chi_iso_antiferro: float = 0.,
            flag_two_channel: bool = False):
        """Calculate factor in front of density for FM_perp.

        Arguments
        ---------
            - hkl
            - space_group_symop
            - ...
        """
        v_hkl_2d_i, v_ferro, v_antiferro = \
            self.calc_factor_in_front_of_density_for_fm(
                hkl, space_group_symop, cell, atom_site_susceptibility, h_loc,
                chi_iso_ferro=chi_iso_ferro,
                chi_iso_antiferro=chi_iso_antiferro,
                flag_two_channel=flag_two_channel)

        k_hkl_1, k_hkl_2, k_hkl_3 = cell.calc_k_loc(*hkl)
        na = numpy.newaxis
        v_hkl_perp_2d_i = calc_moment_perp(
            (k_hkl_1[:, na], k_hkl_2[:, na], k_hkl_3[:, na]), v_hkl_2d_i)
        v_perp_ferro = calc_moment_perp(
            (k_hkl_1[:, na], k_hkl_2[:, na], k_hkl_3[:, na]), v_ferro)
        v_perp_antiferro = calc_moment_perp(
            (k_hkl_1[:, na], k_hkl_2[:, na], k_hkl_3[:, na]), v_antiferro)
        return v_hkl_perp_2d_i, v_perp_ferro, v_perp_antiferro

    def calc_electrons_per_unit_cell(self):
        """
        Calculate electrons per unit cell.

        Returns
        -------
        numpy_atoms : TYPE
            DESCRIPTION.
        m_b_ferro : TYPE
            DESCRIPTION.
        m_b_aferro : TYPE
            DESCRIPTION.

        """
        numpy_density_ferro = self.numpy_density_ferro
        numpy_density_antiferro = self.numpy_density_antiferro
        numpy_multiplicity = self.numpy_multiplicity
        volume_unit_cell = self.volume_unit_cell
        number_unit_cell = self.number_unit_cell
        coeff = volume_unit_cell/float(number_unit_cell)
        m_b_ferro = (numpy_density_ferro*numpy_multiplicity).sum()*coeff
        m_b_aferro = (numpy_density_antiferro*numpy_multiplicity).sum()*coeff

        numpy_density = self.numpy_density
        numpy_basin_atom_label = self.numpy_basin_atom_label
        label_unique = numpy.unique(numpy_basin_atom_label)
        l_lab_m_at = []
        for atom_label in label_unique:
            f_at = (numpy_basin_atom_label == atom_label)
            m_at = (numpy_density[f_at]*numpy_multiplicity[f_at]).sum()*coeff
            l_lab_m_at.append((atom_label, m_at))

        numpy_atoms = numpy.array(l_lab_m_at, dtype=[("label", "U10"),
                                                     ("multiplicity", "f4")])

        return numpy_atoms, m_b_ferro, m_b_aferro

    def renormalize_numpy_densities(self, flag_two_channel: bool = False):
        """
        Renormalize numpy densitites.

        No input arguments. The internal arguments should be defiened.
        """
        numpy_density_ferro = self.numpy_density_ferro
        numpy_density_antiferro = self.numpy_density_antiferro
        numpy_multiplicity = self.numpy_multiplicity
        volume_unit_cell = self.volume_unit_cell
        number_unit_cell = self.number_unit_cell
        coeff = volume_unit_cell/float(number_unit_cell)
        m_b_ferro = (numpy_density_ferro*numpy_multiplicity).sum()*coeff
        m_b_aferro = (numpy_density_antiferro*numpy_multiplicity).sum()*coeff

        numpy_density_ferro /= m_b_ferro
        numpy_density_antiferro /= m_b_aferro

        if not(flag_two_channel):
            numpy_structured_atom = self.numpy_structured_atom

            numpy_density = self.numpy_density
            numpy_basin_atom_label = self.numpy_basin_atom_label

            for structured_atom in numpy_structured_atom:
                atom_label = structured_atom["label"]
                mult_a = structured_atom["multiplicity"]
                f_at = (numpy_basin_atom_label == atom_label)
                m_at = (numpy_density[f_at]*numpy_multiplicity[f_at]).sum() * \
                    coeff
                if m_at == 0.:
                    numpy_density[f_at] = float(mult_a)/(
                        coeff*float((numpy_multiplicity[f_at]).sum()))
                    m_at = float(mult_a)
                numpy_density[f_at] *= float(mult_a)/float(m_at)

# s_cont = """
#   loop_Ho1
#   _density_point_index_x
#   _density_point_index_y
#   _density_point_index_z
#   _density_point_density
#   _density_point_multiplicity
#   0  0  0  0.32  2
#   0  0  1  0.31  1
#   0  0  2  0.29  1
#   """

# obj = DensityPointL.from_cif(s_cont)
# print(obj, end="\n\n")
