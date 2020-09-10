"""Description of MEM class."""
__author__ = 'ikibalin'
__version__ = "2020_08_19"
from typing import NoReturn
import numpy

from cryspy.A_functions_base.function_2_mem import calc_moment_perp, \
    calc_fm_by_density

from cryspy.B_parent_classes.cl_4_global import GlobalN


from cryspy.C_item_loop_classes.cl_1_mem_parameters import MEMParameters
from cryspy.C_item_loop_classes.cl_1_refine_ls import RefineLs
from cryspy.C_item_loop_classes.cl_3_density_point import DensityPointL

from cryspy.E_data_classes.cl_1_crystal import Crystal
from cryspy.E_data_classes.cl_2_diffrn import Diffrn

from cryspy.F_functions_data.script_1_mem import maximize_entropy, \
    refine_susceptibility, make_cycle


class MEM(GlobalN):
    """
    Class to describe RhoChi class.

    Attributes
    ----------
        - crystal_#name (mandatory)
        - diffrn_#name
        - density_point_#name

    Methods
    -------
        - crystals()
        - experiments()
        - maximize_entropy()
        - refine_susceptibility()
        - make_cycle()
        - apply_constraints()
    """

    CLASSES_MANDATORY = (Crystal, Diffrn)
    CLASSES_OPTIONAL = (DensityPointL, MEMParameters)
    # CLASSES_INTERNAL = ()

    CLASSES = CLASSES_MANDATORY + CLASSES_OPTIONAL

    PREFIX = "mem"

    # default values for the parameters
    D_DEFAULT = {"mem_parameters": MEMParameters()}

    def __init__(self, global_name=None, **kwargs) -> NoReturn:
        super(MEM, self).__init__()

        self.__dict__["items"] = []
        self.__dict__["global_name"] = global_name

        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)

    def form_object(self):
        """Form object."""
        self.apply_constraint()

    def apply_constraint(self):
        """Apply constraints."""
        for item in self.items:
            if isinstance(item, Crystal):
                item.apply_constraint()

    def experiments(self):
        """List of expreiments."""
        return [item for item in self.items if isinstance(item, Diffrn)]

    def crystals(self):
        """List of crystals."""
        return [item for item in self.items if isinstance(item, Crystal)]

    def create_prior_density(self):
        """Create prior denisity."""
        crystal = self.crystals()[0]  # FIXME:
        atom_site_susceptibility = crystal.atom_site_susceptibility
        space_group = crystal.space_group
        space_group_symop = space_group.full_space_group_symop
        cell = crystal.cell
        atom_site = crystal.atom_site
        mem_parameters = self.mem_parameters
        prior_density = mem_parameters.prior_density
        points_a = mem_parameters.points_a
        points_b = mem_parameters.points_b
        points_c = mem_parameters.points_c
        flag_two_channel = mem_parameters.method == "2channel"

        l_magnetic_labes = atom_site_susceptibility.label
        density_point = DensityPointL()
        if prior_density == "core":
            atom_electron_configuration = crystal.atom_electron_configuration
            density_point.create_core_density(
                space_group_symop, cell, atom_site,
                atom_electron_configuration, points_a=points_a,
                points_b=points_b, points_c=points_c)
        else:
            density_point.create_flat_density(
                space_group_symop, cell, atom_site,
                l_magnetic_labes=l_magnetic_labes, points_a=points_a,
                points_b=points_b, points_c=points_c,
                flag_two_channel=flag_two_channel)
        self.add_items([density_point])

    def calc_fr(self):
        """Calculate Flip Ratios for diffraction experiments."""
        crystal = self.crystals()[0]  # FIXME:
        l_diffrn = self.experiments()  # FIXME:
        density_point = self.density_point
        mem_parameters = self.mem_parameters
        chi_iso_ferro = mem_parameters.chi_ferro
        chi_iso_antiferro = mem_parameters.chi_antiferro

        flag_two_channel = mem_parameters.method == "2channel"

        cell = crystal.cell
        space_group = crystal.space_group
        space_group_symop = space_group.full_space_group_symop
        atom_site_susceptibility = crystal.atom_site_susceptibility

        total_peaks = 0

        den_i = numpy.array(density_point.density, dtype=float)
        den_ferro_i = numpy.array(density_point.density_ferro, dtype=float)
        den_antiferro_i = numpy.array(density_point.density_antiferro,
                                      dtype=float)
        mult_i = numpy.array(density_point.multiplicity, dtype=int)

        volume = cell.volume
        n_points = mult_i.sum()

        for diffrn in l_diffrn:
            diffrn_orient_matrix = diffrn.diffrn_orient_matrix
            e_up = diffrn_orient_matrix.calc_e_up()
            setup = diffrn.setup
            field = float(setup.field)
            h_loc = (field*e_up[0], field*e_up[1], field*e_up[2])
            diffrn_refln = diffrn.diffrn_refln
            index_h = numpy.array(diffrn_refln.index_h, dtype=int)
            index_k = numpy.array(diffrn_refln.index_k, dtype=int)
            index_l = numpy.array(diffrn_refln.index_l, dtype=int)
            total_peaks += index_h.size
            hkl = (index_h, index_k, index_l)

            f_nucl = crystal.calc_f_nucl(*hkl)
            k_hkl = cell.calc_k_loc(*hkl)
            phase_3d = density_point.calc_phase_3d(hkl, space_group_symop)

            moment_2d, chi_2d_ferro, chi_2d_antiferro = \
                density_point.calc_moment_2d(
                    space_group_symop, cell, atom_site_susceptibility, h_loc,
                    chi_iso_ferro=1., chi_iso_antiferro=1.,
                    flag_two_channel=flag_two_channel)

            chi_ferro = calc_fm_by_density(mult_i, den_ferro_i, n_points,
                                           volume, chi_2d_ferro, phase_3d)
            chi_perp_ferro = calc_moment_perp(k_hkl, chi_ferro)

            chi_aferro = calc_fm_by_density(
                mult_i, den_antiferro_i, n_points, volume, chi_2d_antiferro,
                phase_3d)
            chi_perp_aferro = calc_moment_perp(k_hkl, chi_aferro)

            f_m = calc_fm_by_density(mult_i, den_i, n_points, volume,
                                     moment_2d, phase_3d)
            f_m_perp = calc_moment_perp(k_hkl, f_m)

            f_m_perp_sum = (
                f_m_perp[0] + chi_iso_ferro*chi_perp_ferro[0] +
                chi_iso_antiferro*chi_perp_aferro[0],
                f_m_perp[1] + chi_iso_ferro*chi_perp_ferro[1] +
                chi_iso_antiferro*chi_perp_aferro[1],
                f_m_perp[2] + chi_iso_ferro*chi_perp_ferro[2] +
                chi_iso_antiferro*chi_perp_aferro[2])

            fr_m, delta_fr_m = diffrn.calc_fr(cell, f_nucl, f_m_perp_sum,
                                              delta_f_m_perp=f_m_perp)
            diffrn.diffrn_refln.numpy_fr_calc = fr_m
        for diffrn in l_diffrn:
            # FIXME: not sure that input parameters should be modified
            diffrn.diffrn_refln.numpy_to_items()
            chi_sq, points = diffrn.diffrn_refln.calc_chi_sq_points()
            refine_ls = RefineLs(goodness_of_fit_all=chi_sq/points,
                                 number_reflns=points)
            diffrn.add_items([refine_ls])

    def maximize_entropy(self, c_lambda: float = 1e-3,
                         n_iterations: int = 3000, disp: bool = True):
        """Run entropy maximization.

        Arguments
        ---------
            - prior_density: "uniform" (default) or "core"
            - disp: True (default) or False
        """
        crystal = self.crystals()[0]  # FIXME:
        l_diffrn = self.experiments()  # FIXME:
        mem_parameters = self.mem_parameters
        chi_iso_ferro = mem_parameters.chi_ferro
        chi_iso_antiferro = mem_parameters.chi_antiferro
        points_a = mem_parameters.points_a
        points_b = mem_parameters.points_b
        points_c = mem_parameters.points_c
        prior_density = mem_parameters.prior_density
        flag_two_channel = mem_parameters.method == "2channel"
        gof_desired = mem_parameters.gof_desired

        density_point = maximize_entropy(
            crystal, l_diffrn, c_lambda=c_lambda, n_iterations=n_iterations,
            chi_iso_ferro=chi_iso_ferro, chi_iso_antiferro=chi_iso_antiferro,
            points_a=points_a, points_b=points_b, points_c=points_c,
            prior_density=prior_density, gof_desired=gof_desired,
            flag_two_channel=flag_two_channel, disp=disp)

        self.add_items([density_point])

    def refine_susceptibility(self, disp: bool = True) -> (float, float):
        """Refine susceptibility."""
        crystal = self.crystals()[0]  # FIXME:
        l_diffrn = self.experiments()  # FIXME:
        density_point = self.density_point
        mem_parameters = self.mem_parameters
        chi_iso_ferro = mem_parameters.chi_ferro
        chi_iso_antiferro = mem_parameters.chi_antiferro
        flag_ferro = mem_parameters.chi_ferro_refinement
        flag_antiferro = mem_parameters.chi_antiferro_refinement
        flag_two_channel = mem_parameters.method == "2channel"

        chi_iso_f, chi_iso_af = refine_susceptibility(
            crystal, l_diffrn, density_point, chi_iso_ferro=chi_iso_ferro,
            chi_iso_antiferro=chi_iso_antiferro, flag_ferro=flag_ferro,
            flag_antiferro=flag_antiferro, flag_two_channel=flag_two_channel,
            disp=disp)

        mem_parameters.chi_ferro = chi_iso_f
        mem_parameters.chi_antiferro = chi_iso_af

    def make_cycle(self, c_lambda: float = 1e-3, n_iterations: int = 3000,
                   n_cycle: int = 10, disp: bool = True):
        """Run Rho - Chi cycle refinement."""
        crystal = self.crystals()[0]  # FIXME:
        l_diffrn = self.experiments()  # FIXME:
        mem_parameters = self.mem_parameters
        chi_iso_ferro = mem_parameters.chi_ferro
        chi_iso_antiferro = mem_parameters.chi_antiferro
        flag_ferro = mem_parameters.chi_ferro_refinement
        flag_antiferro = mem_parameters.chi_antiferro_refinement
        points_a = mem_parameters.points_a
        points_b = mem_parameters.points_b
        points_c = mem_parameters.points_c
        prior_density = mem_parameters.prior_density
        flag_two_channel = mem_parameters.method == "2channel"
        gof_desired = mem_parameters.gof_desired

        density_point, chi_iso_f, chi_iso_af = \
            make_cycle(
                crystal, l_diffrn, chi_iso_ferro=chi_iso_ferro,
                chi_iso_antiferro=chi_iso_antiferro, flag_ferro=flag_ferro,
                flag_antiferro=flag_antiferro, points_a=points_a,
                points_b=points_b, points_c=points_c,
                prior_density=prior_density, c_lambda=c_lambda,
                n_iterations=n_iterations, n_cycle=n_cycle,
                gof_desired=gof_desired,
                flag_two_channel=flag_two_channel, disp=disp)
        self.add_items([density_point])
        mem_parameters.chi_ferro = chi_iso_f
        mem_parameters.chi_antiferro = chi_iso_af
