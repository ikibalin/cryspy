# -*- coding: utf-8 -*-
"""The MEM algorithm is described."""
from typing import List
import numpy
import scipy
import scipy.optimize

from cryspy.A_functions_base.function_1_matrices import calc_chi_sq
from cryspy.A_functions_base.function_2_mem import calc_moment_perp, \
    calc_fm_by_density
from cryspy.C_item_loop_classes.cl_3_density_point import DensityPointL
from cryspy.E_data_classes.cl_1_crystal import Crystal
from cryspy.E_data_classes.cl_2_diffrn import Diffrn


def maximize_entropy(crystal: Crystal, l_diffrn: List[Diffrn],
                     c_lambda: float = 1e-7, n_cycle: int = 10000,
                     chi_iso_ferro: float = 0., chi_iso_antiferro: float = 0.,
                     n_x: int = 48, n_y: int = 48, n_z: int = 48,
                     prior_density: str = "uniform",
                     flag_two_channel: bool = False, disp: bool = True) -> \
        DensityPointL:
    """
    Collins algroritm.

    Parameters
    ----------
    crystal : Crystal
        DESCRIPTION.
    l_diffrn : TYPE
        DESCRIPTION.
    c_lambda : float, optional
        DESCRIPTION. The default is 1e-7.
    n_cycle : int, optional
        DESCRIPTION. The default is 10000.
    chi_iso_ferro : float, optional
        DESCRIPTION. The default is 0..
    chi_iso_antiferro : float, optional
        DESCRIPTION. The default is 0..
    n_x : int, optional
        DESCRIPTION. The default is 48.
    n_y : int, optional
        DESCRIPTION. The default is 48.
    n_z : int, optional
        DESCRIPTION. The default is 48.
    prior_density: str, optional
        Choose of prior density: "core" or "uniform". The default is "uniform".
    flag_two_channel: bool, optional
        DESCRIPTION. The default is False.
    disp : bool, optional
        DESCRIPTION. The default is True.

    Returns
    -------
    DensityPointL
        Magnetization density.

    """
    cell = crystal.cell
    space_group = crystal.space_group
    atom_site = crystal.atom_site
    space_group_symop = space_group.full_space_group_symop
    atom_site_susceptibility = crystal.atom_site_susceptibility
    l_magnetic_labes = atom_site_susceptibility.label
    density_point = DensityPointL()
    if prior_density == "core":
        print("The prior density is core's one.", end="\r")
        atom_electron_configuration = crystal.atom_electron_configuration
        density_point.create_core_density(
            space_group_symop, cell, atom_site, atom_electron_configuration,
            points_a=n_x, points_b=n_y, points_c=n_z)
    else:
        print("The prior density is uniform.", end="\r")
        density_point.create_flat_density(
            space_group_symop, cell, atom_site,
            l_magnetic_labes=l_magnetic_labes, points_a=n_x, points_b=n_y,
            points_c=n_z)

    l_f_nucl, l_v_2d_i, l_fr_e, l_fr_s = [], [], [], []
    total_peaks = 0
    for diffrn in l_diffrn:
        diffrn_orient_matrix = diffrn.diffrn_orient_matrix
        e_up = diffrn_orient_matrix.calc_e_up()
        setup = diffrn.setup
        field = float(setup.field)
        h_loc = (field*e_up[0], field*e_up[1], field*e_up[2])
        diffrn_refln = diffrn.diffrn_refln
        ind_h = numpy.array(diffrn_refln.index_h, dtype=int)
        ind_k = numpy.array(diffrn_refln.index_k, dtype=int)
        ind_l = numpy.array(diffrn_refln.index_l, dtype=int)
        total_peaks += ind_h.size
        hkl = (ind_h, ind_k, ind_l)
        fr_e = numpy.array(diffrn_refln.fr, dtype=float)
        fr_s = numpy.array(diffrn_refln.fr_sigma, dtype=float)
        v_hkl_perp_2d_i, v_b_ferro, v_b_antiferro = \
            density_point.calc_factor_in_front_of_density_for_fm_perp(
                hkl, space_group_symop, cell, atom_site_susceptibility, h_loc,
                chi_iso_ferro=chi_iso_ferro,
                chi_iso_antiferro=chi_iso_antiferro)
        f_nucl = crystal.calc_f_nucl(*hkl)
        l_f_nucl.append(f_nucl)
        l_v_2d_i.append((v_hkl_perp_2d_i, v_b_ferro, v_b_antiferro))
        l_fr_e.append(fr_e)
        l_fr_s.append(fr_s)

    def temp_func(numpy_den=None):
        l_chi_sq, l_der_chi_sq = [], []
        l_der_chi_sq_f, l_der_chi_sq_a = [], []
        for diffrn, f_nucl, v_2d_i, fr_e, fr_s in \
                zip(l_diffrn, l_f_nucl, l_v_2d_i, l_fr_e, l_fr_s):

            f_m_perp, delta_f_m_perp, delta_f_m_perp_f, delta_f_m_perp_a = \
                density_point.calc_fm(*v_2d_i)
            fr_m, delta_fr_m = diffrn.calc_fr(cell, f_nucl, f_m_perp,
                                              delta_f_nucl=None,
                                              delta_f_m_perp=delta_f_m_perp)
            delta_fr_m_f = diffrn.calc_fr(cell, f_nucl, f_m_perp,
                                          delta_f_nucl=None,
                                          delta_f_m_perp=delta_f_m_perp_f)[1]
            delta_fr_m_a = diffrn.calc_fr(cell, f_nucl, f_m_perp,
                                          delta_f_nucl=None,
                                          delta_f_m_perp=delta_f_m_perp_a)[1]

            diffrn.diffrn_refln.numpy_fr_calc = fr_m
            chi_sq, der_chi_sq = calc_chi_sq(fr_e, fr_s, fr_m, delta_fr_m)
            der_chi_sq_f = calc_chi_sq(fr_e, fr_s, fr_m, delta_fr_m_f)[1]
            der_chi_sq_a = calc_chi_sq(fr_e, fr_s, fr_m, delta_fr_m_a)[1]
            l_chi_sq.append(chi_sq)
            l_der_chi_sq.append(der_chi_sq)
            l_der_chi_sq_f.append(der_chi_sq_f)
            l_der_chi_sq_a.append(der_chi_sq_a)
        # print(" ".join([f" {val:10.2f}" for val in l_chi_sq]))
        return sum(l_chi_sq), sum(l_der_chi_sq), sum(l_der_chi_sq_f), \
            sum(l_der_chi_sq_a)

    chi_sq_diff_prev, chi_sq_diff = numpy.inf, numpy.inf
    chi_sq_prev, flag = numpy.inf, False

    i_cycle = 0
    while i_cycle <= n_cycle:
        i_cycle += 1
        numpy_density = density_point.numpy_density
        numpy_density_ferro = density_point.numpy_density_ferro
        numpy_density_antiferro = density_point.numpy_density_antiferro
        chi_sq, delta_chi_sq, delta_chi_sq_f, delta_chi_sq_a = temp_func()
        chi_sq_diff = chi_sq_prev - chi_sq
        if disp:
            chi_sq_n = chi_sq/float(total_peaks)
            print(f"cycle {i_cycle:5}. chi_sq/n: {chi_sq_n:.2f} {c_lambda*10**5:.3f}",
                  end="\r")
        if chi_sq > chi_sq_prev:
            if flag:
                print("U...")
                break
            flag = True
        else:
            flag = False

        if chi_sq_diff < chi_sq_diff_prev:  # FIXME: not sure about this
            c_lambda = 1.01 * c_lambda
        else:
            c_lambda = 0.99 * c_lambda

        chi_sq_diff_prev = chi_sq_diff
        chi_sq_prev = chi_sq

        numpy_density_new = numpy_density*numpy.exp(-c_lambda*delta_chi_sq)
        numpy_density_ferro_new = \
            numpy_density_ferro*numpy.exp(-c_lambda*delta_chi_sq_f)
        numpy_density_antiferro_new = \
            numpy_density_antiferro*numpy.exp(-c_lambda*delta_chi_sq_a)

        density_point.numpy_density = numpy_density_new
        density_point.numpy_density_ferro = numpy_density_ferro_new
        density_point.numpy_density_antiferro = numpy_density_antiferro_new
        density_point.renormalize_numpy_densities()
        if chi_sq_n < 1.:
            print(f"at cycle {i_cycle:5} chi_sq/n is less than 1.", end="\r")
            break
    density_point.numpy_to_items()
    for diffrn in l_diffrn:
        # FIXME: not sure that input parameters should be modified
        diffrn.diffrn_refln.numpy_to_items()
    return density_point


def refine_susceptibility(crystal: Crystal, l_diffrn: List[Diffrn],
                          density_point: DensityPointL,
                          chi_iso_ferro: float = 0.,
                          chi_iso_antiferro: float = 0.,
                          flag_ferro: bool = True, flag_antiferro: bool = True,
                          disp: bool = True) -> (float, float):
    """
    Refinement of susceptibility.

    Parameters
    ----------
    crystal : Crystal
        DESCRIPTION.
    l_diffrn : TYPE
        DESCRIPTION.
    density_point : DensityPointL
        DESCRIPTION.
    chi_iso_ferro : float, optional
        DESCRIPTION. The default is 0..
    chi_iso_antiferro : float, optional
        DESCRIPTION. The default is 0..
    flag_ferro : bool, optional
        DESCRIPTION. The default is True.
    flag_antiferro : bool, optional
        DESCRIPTION. The default is True.
    disp : bool, optional
        DESCRIPTION. The default is True.

    Returns
    -------
    float
        Chi_iso_ferro.
    float
        Chi_iso_antiferro.

    """
    cell = crystal.cell
    space_group = crystal.space_group
    atom_site = crystal.atom_site
    space_group_symop = space_group.full_space_group_symop
    atom_site_susceptibility = crystal.atom_site_susceptibility
    # l_magnetic_labes = atom_site_susceptibility.label

    volume = density_point.volume_unit_cell
    np = density_point.number_unit_cell

    den_i = density_point.numpy_density
    den_ferro_i = density_point.numpy_density_ferro
    den_antiferro_i = density_point.numpy_density_antiferro
    mult_i = density_point.numpy_multiplicity

    l_f_nucl, l_fr_e, l_fr_s, l_e_up, l_h_loc = [], [], [], [], []
    l_k_hkl, l_phase_3d = [], []
    l_chi_perp_ferro, l_chi_perp_antiferro = [], []
    total_peaks = 0
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
        fr_e = numpy.array(diffrn_refln.fr, dtype=float)
        fr_s = numpy.array(diffrn_refln.fr_sigma, dtype=float)
        f_nucl = crystal.calc_f_nucl(*hkl)
        k_hkl = cell.calc_k_loc(*hkl)
        phase_3d = density_point.calc_phase_3d(hkl, space_group_symop)

        moment_2d, chi_2d_ferro, chi_2d_antiferro = \
            density_point.calc_moment_2d(
                space_group_symop, cell, atom_site_susceptibility, h_loc, 1.,
                1.)

        chi_ferro = calc_fm_by_density(mult_i, den_ferro_i, np, volume,
                                       chi_2d_ferro, phase_3d)
        chi_perp_ferro = calc_moment_perp(k_hkl, chi_ferro)

        chi_aferro = calc_fm_by_density(mult_i, den_antiferro_i, np, volume,
                                        chi_2d_antiferro, phase_3d)
        chi_perp_aferro = calc_moment_perp(k_hkl, chi_aferro)

        l_f_nucl.append(f_nucl)
        l_fr_e.append(fr_e)
        l_fr_s.append(fr_s)
        l_e_up.append(e_up)
        l_h_loc.append(h_loc)
        l_k_hkl.append(k_hkl)
        l_phase_3d.append(phase_3d)
        l_chi_perp_ferro.append(chi_perp_ferro)
        l_chi_perp_antiferro.append(chi_perp_aferro)

    l_name = atom_site_susceptibility.get_variable_names()
    l_par_0 = [atom_site_susceptibility.get_variable_by_name(name) for name
               in l_name]
    if flag_ferro:
        l_par_0.append(chi_iso_ferro)
    if flag_antiferro:
        l_par_0.append(chi_iso_antiferro)

    def temp_func(l_par):
        if (flag_ferro & flag_antiferro):
            for name, parameter in zip(l_name, l_par[:-2]):
                atom_site_susceptibility.set_variable_by_name(name, parameter)
            chi_iso_f = float(l_par[-2])
            chi_iso_af = float(l_par[-1])
        elif (flag_ferro & (not(flag_antiferro))):
            for name, parameter in zip(l_name, l_par[:-1]):
                atom_site_susceptibility.set_variable_by_name(name, parameter)
            chi_iso_f = float(l_par[-1])
            chi_iso_af = chi_iso_antiferro
        elif ((not(flag_ferro)) & flag_antiferro):
            for name, parameter in zip(l_name, l_par[:-1]):
                atom_site_susceptibility.set_variable_by_name(name, parameter)
            chi_iso_f = chi_iso_ferro
            chi_iso_af = float(l_par[-1])
        else:
            for name, parameter in zip(l_name, l_par):
                atom_site_susceptibility.set_variable_by_name(name, parameter)
            chi_iso_f = chi_iso_ferro
            chi_iso_af = chi_iso_antiferro

        if chi_iso_f < 0.:
            chi_iso_f = 0
        if chi_iso_af > 0.:
            chi_iso_af = 0

        atom_site_susceptibility.apply_space_group_constraint(atom_site,
                                                              space_group)
        l_chi_sq = []
        l_der_chi_sq, l_der_chi_sq_f, l_der_chi_sq_a = [], [], []
        for diffrn, f_nucl, fr_e, fr_s, e_up, h_loc, k_hkl, phase_3d, \
            chi_perp_ferro, chi_perp_aferro in \
            zip(l_diffrn, l_f_nucl, l_fr_e, l_fr_s, l_e_up, l_h_loc, l_k_hkl,
                l_phase_3d, l_chi_perp_ferro, l_chi_perp_antiferro):

            moment_2d, moment_ferro, moment_antiferro = \
                density_point.calc_moment_2d(
                    space_group_symop, cell, atom_site_susceptibility, h_loc,
                    chi_iso_ferro, chi_iso_antiferro)

            f_m = calc_fm_by_density(mult_i, den_i, np, volume, moment_2d,
                                     phase_3d)
            f_m_perp = calc_moment_perp(k_hkl, f_m)

            # add ferro and anti_ferro

            f_m_perp_sum = (
                f_m_perp[0] + chi_iso_f*chi_perp_ferro[0] +
                chi_iso_af*chi_perp_aferro[0],
                f_m_perp[1] + chi_iso_f*chi_perp_ferro[1] +
                chi_iso_af*chi_perp_aferro[1],
                f_m_perp[2] + chi_iso_f*chi_perp_ferro[2] +
                chi_iso_af*chi_perp_aferro[2])

            fr_m, delta_fr_m = diffrn.calc_fr(cell, f_nucl, f_m_perp_sum,
                                              delta_f_m_perp=f_m_perp)
            delta_fr_m_f = delta_fr_m
            delta_fr_m_a = delta_fr_m

            diffrn.diffrn_refln.numpy_fr_calc = fr_m

            chi_sq, der_chi_sq = calc_chi_sq(fr_e, fr_s, fr_m, delta_fr_m)
            der_chi_sq_f = calc_chi_sq(fr_e, fr_s, fr_m, delta_fr_m_f)[1]
            der_chi_sq_a = calc_chi_sq(fr_e, fr_s, fr_m, delta_fr_m_a)[1]
            l_chi_sq.append(chi_sq)
            l_der_chi_sq.append(der_chi_sq)
            l_der_chi_sq_f.append(der_chi_sq_f)
            l_der_chi_sq_a.append(der_chi_sq_a)
        return sum(l_chi_sq)
    chi_sq = temp_func(l_par_0)
    print(f"Chi_sq before optimization {chi_sq/total_peaks:.5f}.           ",
          end="\r")
    res = scipy.optimize.minimize(temp_func, l_par_0, method="BFGS",
                                  callback=func_temp, options={"eps": 0.001})
    l_param = res.x
    chi_sq_new = res.fun

    hess_inv = res["hess_inv"]
    sigma = (abs(numpy.diag(hess_inv)))**0.5

    print(f"Chi_sq after optimization {chi_sq_new/total_peaks:.5f}.          ",
          end="\r")

    if (flag_ferro & flag_antiferro):
        for name, parameter, sig in zip(l_name, l_param[:-2], sigma[:-2]):
            atom_site_susceptibility.set_variable_by_name(name, parameter)
            name_sig = name[:-1] + ((f"{name[-1][0]:}_sigma", name[-1][1]), )
            atom_site_susceptibility.set_variable_by_name(name_sig, sig)
        chi_iso_f = float(l_param[-2])
        chi_iso_af = float(l_param[-1])
    elif (flag_ferro & (not(flag_antiferro))):
        for name, parameter, sig in zip(l_name, l_param[:-1], sigma[:-1]):
            atom_site_susceptibility.set_variable_by_name(name, parameter)
            name_sig = name[:-1] + ((f"{name[-1][0]:}_sigma", name[-1][1]), )
            atom_site_susceptibility.set_variable_by_name(name_sig, sig)
        chi_iso_f = float(l_param[-1])
        chi_iso_af = chi_iso_antiferro
    elif ((not(flag_ferro)) & flag_antiferro):
        for name, parameter, sig in zip(l_name, l_param[:-1], sigma[:-1]):
            atom_site_susceptibility.set_variable_by_name(name, parameter)
            name_sig = name[:-1] + ((f"{name[-1][0]:}_sigma", name[-1][1]), )
            atom_site_susceptibility.set_variable_by_name(name_sig, sig)
        chi_iso_f = chi_iso_ferro
        chi_iso_af = float(l_param[-1])
    else:
        for name, parameter, sig in zip(l_name, l_param, sigma):
            atom_site_susceptibility.set_variable_by_name(name, parameter)
            name_sig = name[:-1] + ((f"{name[-1][0]:}_sigma", name[-1][1]), )
            atom_site_susceptibility.set_variable_by_name(name_sig, sig)
        chi_iso_f = chi_iso_ferro
        chi_iso_af = chi_iso_antiferro

    if chi_iso_f < 0.:
        chi_iso_f = 0.
    if chi_iso_af > 0.:
        chi_iso_af = 0.
    for diffrn in l_diffrn:
        # FIXME: not sure that input parameters should be modified
        diffrn.diffrn_refln.numpy_to_items()

    return chi_iso_f, chi_iso_af


def func_temp(*argv):
    """Show data at optimization procedure."""
    l_param = argv[0]
    s_out = " ".join([f"{_:10.5f}" for _ in l_param])
    print(f"{s_out:}", end="\r")


def make_cycle(crystal: Crystal, l_diffrn: List[Diffrn],
               chi_iso_ferro: float = 0., chi_iso_antiferro: float = 0.,
               flag_ferro: bool = True, flag_antiferro: bool = True,
               disp: bool = True,
               points_a: int = 48, points_b: int = 48, points_c: int = 48,
               prior_density: str = "uniform",
               c_lambda: float = 1e-6, n_mem: int = 51, n_cycle: int = 10):
    """Rho - Chi cycle."""
    chi_iso_ferro_new, chi_iso_antiferro_new = chi_iso_ferro, chi_iso_antiferro
    # atom_site_susceptibility = crystal.atom_site_susceptibility
    # l_var = atom_site_susceptibility.get_variables()
    chi_iso_ferro_new = float(chi_iso_ferro)
    chi_iso_antiferro_new = float(chi_iso_antiferro)
    for i_cycle in range(n_cycle):
        print(f"Cycle {(i_cycle+1):4}. Entropy maximization                  ",
              end="\r")
        density_point = maximize_entropy(
            crystal, l_diffrn, c_lambda, n_mem,
            chi_iso_ferro=chi_iso_ferro_new,
            chi_iso_antiferro=chi_iso_antiferro_new,
            n_x=points_a, n_y=points_b, n_z=points_c,
            prior_density=prior_density, disp=disp)
        print(f"Cycle {(i_cycle+1):4}. Chi refinement                        ",
              end="\r")
        chi_iso_ferro_new, chi_iso_antiferro_new = refine_susceptibility(
            crystal, l_diffrn, density_point, chi_iso_ferro=chi_iso_ferro_new,
            chi_iso_antiferro=chi_iso_antiferro_new,
            flag_ferro=flag_ferro, flag_antiferro=flag_antiferro,
            disp=disp)
    return density_point, chi_iso_ferro_new, chi_iso_antiferro_new
