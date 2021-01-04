"""Description of MEM class."""
__author__ = 'ikibalin'
__version__ = "2020_08_19"
from typing import NoReturn
import numpy
import matplotlib.pyplot as plt

from cryspy.A_functions_base.function_2_crystallography_base import \
    calc_atoms_in_unit_cell
from cryspy.A_functions_base.function_2_mem import calc_moment_perp, \
    calc_fm_by_density, calc_index_atom_symmetry_closest_to_fract_xyz

from cryspy.B_parent_classes.cl_4_global import GlobalN

from cryspy.C_item_loop_classes.cl_3_density_point import DensityPointL
from cryspy.C_item_loop_classes.cl_1_mem_parameters import MEMParameters
from cryspy.C_item_loop_classes.cl_1_refine_ls import RefineLs
from cryspy.C_item_loop_classes.cl_2_section import SectionL

from cryspy.D_functions_item_loop.function_1_section_from_density_point \
    import calc_section_from_density_point
    
from cryspy.E_data_classes.cl_1_crystal import Crystal
from cryspy.E_data_classes.cl_2_diffrn import Diffrn

from cryspy.F_functions_data.script_1_mem import maximize_entropy, \
    refine_susceptibility, make_cycle, calc_moments_in_unit_cell, \
    calc_section_for_mem


class MEM(GlobalN):
    """
    Maximum Entropy Method: two-channel approach or chi approach.

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
    CLASSES_OPTIONAL = (DensityPointL, MEMParameters, SectionL)
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

    def save_to_file_den(self, f_name: str = "file.den", label_atom=None,
                         f_background: str = "file_back.den"):
        """Save density to ".den" files."""
        crystal = self.crystals()[0]
        space_group = crystal.space_group
        cell = crystal.cell
        mem_parameters = self.mem_parameters
        flag_two_channel = mem_parameters.method == "2channel"
        if flag_two_channel:
            atom_site_susceptibility = None
        else:
            atom_site_susceptibility = crystal.atom_site_susceptibility

        density_point = self.density_point
        density_point.save_to_file_den(
            mem_parameters, space_group, cell, f_name=f_name,
            atom_site_susceptibility=atom_site_susceptibility,
            f_background=f_background)

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
                points_b=points_b, points_c=points_c,
                flag_two_channel=flag_two_channel)
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

        # FIXME: temporary solution of calculation rbs_i if it's not defined
        density_point.calc_rbs_i(
            space_group_symop, points_a=mem_parameters.points_a,
            points_b=mem_parameters.points_b, points_c=mem_parameters.points_c)

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
            diffrn.diffrn_refln.numpy_to_items()
            chi_sq, points = diffrn.diffrn_refln.calc_chi_sq_points()
            refine_ls = RefineLs(goodness_of_fit_all=chi_sq/points,
                                 number_reflns=points)
            diffrn.add_items([refine_ls])

    def maximize_entropy(
            self, c_lambda: float = 1e-3, n_iterations: int = 3000,
            disp: bool = True, d_info: dict = None):
        """Run entropy maximization.

        Arguments
        ---------
            - prior_density: "uniform" (default) or "core"
            - disp: True (default) or False
        """
        crystal = self.crystals()[0]  # FIXME:
        l_diffrn = self.experiments()  # FIXME:
        mem_parameters = self.mem_parameters

        density_point = maximize_entropy(
            crystal, l_diffrn, mem_parameters, c_lambda=c_lambda,
            n_iterations=n_iterations, disp=disp, d_info=d_info)

        self.add_items([density_point])

    def refine_susceptibility(self, disp: bool = True, d_info: dict = None) \
            -> (float, float):
        """Refine susceptibility."""
        crystal = self.crystals()[0]  # FIXME:
        l_diffrn = self.experiments()  # FIXME:
        density_point = self.density_point
        mem_parameters = self.mem_parameters

        # FIXME: temporary solution
        density_point.volume_unit_cell = crystal.cell.volume
        density_point.number_unit_cell = \
            mem_parameters.points_a * mem_parameters.points_b * \
            mem_parameters.points_c

        refine_susceptibility(crystal, l_diffrn, density_point, mem_parameters,
                              disp=disp, d_info=d_info)

    def make_cycle(self, c_lambda: float = 1e-3, n_iterations: int = 3000,
                   n_cycle: int = 10, disp: bool = True, d_info: dict = None):
        """Run Rho - Chi cycle refinement."""
        crystal = self.crystals()[0]  # FIXME:
        l_diffrn = self.experiments()  # FIXME:
        mem_parameters = self.mem_parameters

        density_point = \
            make_cycle(crystal, l_diffrn, mem_parameters, c_lambda=c_lambda,
                       n_iterations=n_iterations, n_cycle=n_cycle, disp=disp,
                       d_info=d_info)
        self.add_items([density_point])

    def calc_moments_in_unit_cell(self, field_loc,
                                  f_name: str = "moments_in_unit_cell.dat"):
        """Calc moments in an unit cell.

        If f_name is given, the result is save in file.
        If it is None, the results are given in output.
        """
        density_point = self.density_point
        crystal = self.crystals()[0]
        mem_parameters = self.mem_parameters
        fract_x, fract_y, fract_z, moment_x, moment_y, moment_z = \
            calc_moments_in_unit_cell(density_point, crystal, mem_parameters,
                                      field_loc)
        if f_name is None:
            return fract_x, fract_y, fract_z, moment_x, moment_y, moment_z
        else:
            ls_out = ["# Moments in an unit cell"]
            ls_out.append("loop_\n_fract_x\n_fract_y\n_fract_z")
            ls_out.append("_moment_x\n_moment_y\n_moment_z")
            for f_x, f_y, f_z, m_x, m_y, m_z in zip(
                    fract_x, fract_y, fract_z, moment_x, moment_y, moment_z):
                ls_out.append(f"{f_x:9.5f} {f_y:9.5f} {f_z:9.5f} \
{m_x:10.3f} {m_y:10.3f} {m_z:10.3f}")
            with open(f_name, "w") as fid:
                fid.write("\n".join(ls_out))

    def calc_sections(self):
        """Calc moments in an unit cell.

        If f_name is given, the result is save in file.
        If it is None, the results are given in output.
        """
        density_point = self.density_point
        crystal = self.crystals()[0]
        loop_section = self.section
        mem_parameters = self.mem_parameters

        for section in loop_section.items:
            calc_section_for_mem(section, density_point, crystal,
                                 mem_parameters)

    def plots(self):
        return [self.plot_section(), ]
    
    def plot_section(self):
        """Plot first section
        """
        if not((self.is_attribute("section") &
                self.is_attribute("density_point") &
                self.is_attribute("mem_parameters") )):
            return
        if len(self.section.items) == 0:
            return

        section = self.section[0]
        crystal = self.crystals()[0]
        space_group = crystal.space_group
        f_s_g_s = space_group.full_space_group_symop
        r_11 = numpy.array(f_s_g_s.r_11, dtype=float)
        r_12 = numpy.array(f_s_g_s.r_12, dtype=float)
        r_13 = numpy.array(f_s_g_s.r_13, dtype=float)
        r_21 = numpy.array(f_s_g_s.r_21, dtype=float)
        r_22 = numpy.array(f_s_g_s.r_22, dtype=float)
        r_23 = numpy.array(f_s_g_s.r_23, dtype=float)
        r_31 = numpy.array(f_s_g_s.r_31, dtype=float)
        r_32 = numpy.array(f_s_g_s.r_32, dtype=float)
        r_33 = numpy.array(f_s_g_s.r_33, dtype=float)

        r_ij = (r_11, r_12, r_13, r_21, r_22, r_23, r_31, r_32, r_33)

        b_1 = numpy.array(f_s_g_s.b_1, dtype=float)
        b_2 = numpy.array(f_s_g_s.b_2, dtype=float)
        b_3 = numpy.array(f_s_g_s.b_3, dtype=float)

        b_i = (b_1, b_2, b_3)

        atom_site = crystal.atom_site
        fract_x = atom_site.numpy_fract_x
        fract_y = atom_site.numpy_fract_y
        fract_z = atom_site.numpy_fract_z
        fract_xyz = (fract_x, fract_y, fract_z)

        atom_label = atom_site.numpy_label

        fract_uc_x, fract_uc_y, fract_uc_z, label_uc = \
            calc_atoms_in_unit_cell(r_ij, b_i, fract_xyz, atom_label)
    
        cell = crystal.cell
        atom_site_susceptibility = crystal.atom_site_susceptibility

        density_point = self.density_point
        mem_parameters = self.mem_parameters

        atom_x, atom_y, atom_label = section.calc_atoms(
            cell, atom_site, f_s_g_s, distance_min=0.3)

        den_chi_section, den_b_section = \
            calc_section_from_density_point(
            section, density_point, mem_parameters, cell, f_s_g_s,
            atom_site, atom_site_susceptibility)

        fract_atom_xyz = numpy.array(fract_xyz, dtype=float
                                     ).transpose()

        fract_sec_xyz = section.calc_fractions(cell, atom_site)
        fract_sec_xyz = numpy.transpose(numpy.array(fract_sec_xyz,
                                                    dtype=float))

        n_atom_index, n_symmetry, distance = \
            calc_index_atom_symmetry_closest_to_fract_xyz(
                fract_sec_xyz, fract_atom_xyz, r_ij, b_i, cell)
        n_at_2d = numpy.transpose(n_atom_index.reshape(
            section.points_x, section.points_y))

        fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(4.2, 4.2),
                                       dpi=300)
        plt.set_cmap('Accent')
        ax1.imshow(n_at_2d, interpolation='bilinear',
                   extent=(-0.5*section.size_x, 0.5*section.size_x,
                           -0.5*section.size_y, 0.5*section.size_y),
                   alpha=0.1, origin="lower")

        den_x = numpy.linspace(-0.5*section.size_x, 0.5*section.size_x,
                               section.points_x)
        den_y = numpy.linspace(-0.5*section.size_y, 0.5*section.size_y,
                               section.points_y)

        blk = '#000000'
        ax1.contour(den_x, den_y, den_chi_section.transpose(),
                    levels=[0.1, 0.5, 1., 5., 10., 50.],
                    colors=[blk, blk, blk, blk, blk, blk],
                    linewidths=0.5)

        ax1.plot(atom_x, atom_y, 'ko', ms=3)
        for _1, _2, _3 in zip(atom_x, atom_y, atom_label):
            ax1.text(_1, _2, _3)
        ax1.set_title(
            f"Tensor. Max is {den_chi_section.max():.1f}")

        # plt.set_cmap('RdBu')
        ax2.imshow(n_at_2d, interpolation='bilinear',
                   extent=(-0.5*section.size_x, 0.5*section.size_x,
                           -0.5*section.size_y, 0.5*section.size_y),
                   alpha=0.1, origin="lower")

        hh = numpy.abs(den_b_section).max()
        rd = '#FF0000'
        ax2.contour(den_x, den_y, den_b_section.transpose(),
                    levels=[-50., -10., -5., -1., -0.5, -0.1, 
                            0.1, 0.5, 1., 5., 10., 50.],
                    colors=[rd, rd, rd, rd, rd, rd,
                            blk, blk, blk, blk, blk, blk],
                    linewidths=0.5)
        # ax2.imshow(den_b_section, interpolation='bilinear',
        #            extent=(-0.5*section.size_x, 0.5*section.size_x,
        #                    -0.5*section.size_y, 0.5*section.size_y),
        #            vmin=-hh, vmax=hh,
        #            alpha=1., origin="lower")
        ax2.set_title(f"2channel. Max is {hh:.1f}")
        ax2.plot(atom_x, atom_y, 'ko', ms=3)
        for _1, _2, _3 in zip(atom_x, atom_y, atom_label):
            ax2.text(_1, _2, _3)
    
        return (fig, (ax1, ax2))