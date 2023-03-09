# -*- coding: utf-8 -*-
"""The MEM algorithm is described."""
from typing import List
import numpy
import scipy
import scipy.optimize
import copy

from cryspy.A_functions_base.function_1_matrices import calc_chi_sq, \
    tri_linear_interpolation, calc_product_matrix_vector
from cryspy.A_functions_base.function_2_mem import calc_moment_perp, \
    calc_fm_by_density, transfer_to_density_3d, transfer_to_chi_3d

from cryspy.C_item_loop_classes.cl_1_refine_ls import RefineLs
from cryspy.C_item_loop_classes.cl_1_mem_parameters import MEMParameters

from cryspy.C_item_loop_classes.cl_2_section import Section
from cryspy.C_item_loop_classes.cl_3_density_point import DensityPointL

from cryspy.E_data_classes.cl_1_crystal import Crystal
from cryspy.E_data_classes.cl_2_diffrn import Diffrn

from cryspy.D_functions_item_loop.function_1_flip_ratio import \
    calc_fm_perp_loc, calc_e_up_loc


def choose_max_clambda(c_lambda: float, den: numpy.ndarray, der: numpy.ndarray,
                       rel_diff: float = 0.1) -> float:
    """Search optimal c_lambda.

    den_new = den * exp(-c_lambda * der)
    c_lambda: (den_new - den) < rel_diff
    """
    exp_c_lambda = numpy.exp(-c_lambda*der)
    rel_max = (numpy.abs(1.-exp_c_lambda)).max()
    arg_max = (numpy.abs(1.-exp_c_lambda)).argmax()
    if (rel_max < rel_diff) & (rel_max > 0.0001 * rel_diff):
        return c_lambda
    elif (rel_max < 0.0001 * rel_diff):
        rel_diff = 0.01 * rel_diff

    # FIXME: check it.
    if exp_c_lambda[arg_max] > 1.:
        c_lambda_new = -1*numpy.log(1+rel_diff)/der[arg_max]
    else:  # <= 1
        c_lambda_new = -1*numpy.log(1-rel_diff)/der[arg_max]
    return c_lambda_new


def maximize_entropy(crystal: Crystal, l_diffrn: List[Diffrn],
                     mem_parameters: MEMParameters,
                     c_lambda: float = 1e-7, n_iterations: int = 10000,
                     disp: bool = True, d_info: dict = None) -> DensityPointL:
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
    n_iterations : int, optional
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
    flag_info = d_info is not None

    if flag_info:
        d_info_keys = d_info.keys()
        if "stop" not in d_info_keys:
            d_info["stop"] = False
        if "print" not in d_info_keys:
            d_info["print"] = ""

    chi_iso_ferro = mem_parameters.chi_ferro
    chi_iso_antiferro = mem_parameters.chi_antiferro
    points_a = mem_parameters.points_a
    points_b = mem_parameters.points_b
    points_c = mem_parameters.points_c
    prior_density = mem_parameters.prior_density
    flag_two_channel = mem_parameters.method == "2channel"
    gof_desired = mem_parameters.gof_desired

    cell = crystal.cell
    space_group = crystal.space_group
    atom_site = crystal.atom_site
    space_group_symop = space_group.full_space_group_symop
    if crystal.is_attribute("atom_site_susceptibility"):
        atom_site_susceptibility = crystal.atom_site_susceptibility
        l_magnetic_labes = atom_site_susceptibility.label
    else:
        atom_site_susceptibility = None
        l_magnetic_labes  = []
    density_point = DensityPointL()
    if prior_density == "core":
        print("The prior density is core's one.", end="\r")
        if flag_info:
            d_info["print"] = "The prior density is core's one."

        atom_electron_configuration = crystal.atom_electron_configuration
        density_point.create_core_density(
            space_group_symop, cell, atom_site, atom_electron_configuration,
            points_a=points_a, points_b=points_b, points_c=points_c,
            flag_two_channel=flag_two_channel)
    else:
        print("The prior density is uniform.", end="\r")
        if flag_info:
            d_info["print"] = "The prior density is uniform."
        density_point.create_flat_density(
            space_group_symop, cell, atom_site,
            l_magnetic_labes=l_magnetic_labes, points_a=points_a,
            points_b=points_b, points_c=points_c,
            flag_two_channel=flag_two_channel)

    l_f_nucl, l_v_2d_i, l_fr_e, l_fr_s, l_fm_orb_perp_loc = [], [], [], [], []
    total_peaks = 0

    for diffrn in l_diffrn:
        diffrn_orient_matrix = diffrn.diffrn_orient_matrix
        u_matrix = diffrn_orient_matrix.u
        e_up = calc_e_up_loc(0., 0., 0., u_matrix)

        setup = diffrn.setup
        field = float(setup.field)
        h_loc = (field*e_up[0], field*e_up[1], field*e_up[2])
        diffrn_refln = diffrn.diffrn_refln
        ind_h = numpy.array(diffrn_refln.index_h, dtype=int)
        ind_k = numpy.array(diffrn_refln.index_k, dtype=int)
        ind_l = numpy.array(diffrn_refln.index_l, dtype=int)
        total_peaks += ind_h.size

        chi_m = crystal.calc_susceptibility_moment_tensor(
            ind_h, ind_k, ind_l, flag_only_orbital=True)
        sft_ij = chi_m[:9]
        sftm_ij = chi_m[9:]

        k_loc_i = cell.calc_k_loc(ind_h, ind_k, ind_l)
        fm_orb_perp_loc = calc_fm_perp_loc(e_up, field, k_loc_i, sft_ij,
                                           sftm_ij)

        hkl = (ind_h, ind_k, ind_l)
        index_hkl = numpy.stack([ind_h, ind_k, ind_l], dtype=bool)
        fr_e = numpy.array(diffrn_refln.fr, dtype=float)
        fr_s = numpy.array(diffrn_refln.fr_sigma, dtype=float)
        v_hkl_perp_2d_i, v_b_ferro, v_b_antiferro = \
            density_point.calc_factor_in_front_of_density_for_fm_perp(
                hkl, space_group_symop, cell, atom_site_susceptibility, h_loc,
                chi_iso_ferro=chi_iso_ferro,
                chi_iso_antiferro=chi_iso_antiferro,
                flag_two_channel=flag_two_channel)
        f_nucl = crystal.calc_f_nucl(index_hkl)
        l_f_nucl.append(f_nucl)
        l_v_2d_i.append((v_hkl_perp_2d_i, v_b_ferro, v_b_antiferro))
        l_fr_e.append(fr_e)
        l_fr_s.append(fr_s)
        l_fm_orb_perp_loc.append(fm_orb_perp_loc)

    def temp_func(numpy_den=None):
        l_chi_sq, l_der_chi_sq = [], []
        l_der_chi_sq_f, l_der_chi_sq_a = [], []
        for diffrn, f_nucl, v_2d_i, fr_e, fr_s, fm_orb_perp_loc in \
                zip(l_diffrn, l_f_nucl, l_v_2d_i, l_fr_e, l_fr_s,
                    l_fm_orb_perp_loc):
            diffrn_refln = diffrn.diffrn_refln

            f_m_perp, delta_f_m_perp, delta_f_m_perp_f, delta_f_m_perp_a = \
                density_point.calc_fm(*v_2d_i)

            # FIXME: put condition
            f_m_perp = (f_m_perp[0] + fm_orb_perp_loc[0],
                        f_m_perp[1] + fm_orb_perp_loc[1],
                        f_m_perp[2] + fm_orb_perp_loc[2])

            fr_m, delta_fr_m = diffrn.calc_fr(cell, f_nucl, f_m_perp,
                                              delta_f_nucl=None,
                                              delta_f_m_perp=delta_f_m_perp)
            delta_fr_m_f = diffrn.calc_fr(cell, f_nucl, f_m_perp,
                                          delta_f_nucl=None,
                                          delta_f_m_perp=delta_f_m_perp_f)[1]
            delta_fr_m_a = diffrn.calc_fr(cell, f_nucl, f_m_perp,
                                          delta_f_nucl=None,
                                          delta_f_m_perp=delta_f_m_perp_a)[1]

            diffrn_refln.numpy_fr_calc = fr_m
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

    chi_sq_best, chi_sq_n_diff = numpy.inf, numpy.inf
    numpy_density_best, numpy_density_ferro_best = None, None
    numpy_density_antiferro_best = None
    delta_chi_sq_best, delta_chi_sq_f_best = None, None
    delta_chi_sq_a_best = None

    mult_i = density_point.numpy_multiplicity

    c_lambda_min = 1e-9  # min value
    i_cycle = 0
    while i_cycle <= n_iterations:
        i_cycle += 1
        numpy_density = density_point.numpy_density
        numpy_density_ferro = density_point.numpy_density_ferro
        numpy_density_antiferro = density_point.numpy_density_antiferro
        chi_sq, delta_chi_sq, delta_chi_sq_f, delta_chi_sq_a = temp_func()

        chi_sq_n = chi_sq/float(total_peaks)
        if disp:
            print(f"cycle {i_cycle:5}. chi_sq/n: {chi_sq_n:.2f} \
{c_lambda*10**5:.3f} {chi_sq_n_diff:.4f}",
                  end="\r")
        if flag_info:
            d_info["print"] = f"""Iteration {i_cycle:5}:
    chi_sq/n: {chi_sq_n:.2f};
    c_lambda: {c_lambda*10**5:.3f} * 10**-5;
    difference of chi_sq/n: {chi_sq_n_diff:.4f}."""

        if chi_sq > chi_sq_best:
            numpy_density = numpy_density_best
            numpy_density_ferro = numpy_density_ferro_best
            numpy_density_antiferro = numpy_density_antiferro_best
            delta_chi_sq = delta_chi_sq_best
            delta_chi_sq_f = delta_chi_sq_f_best
            delta_chi_sq_a = delta_chi_sq_a_best
            c_lambda = 0.5*c_lambda
        else:
            chi_sq_n_diff = (chi_sq_best - chi_sq)/float(total_peaks)
            chi_sq_best = chi_sq
            c_lambda = 1.03*c_lambda

        if not(flag_two_channel):
            numpy_density_best = copy.deepcopy(numpy_density)
            delta_chi_sq_best = copy.deepcopy(delta_chi_sq)
        numpy_density_ferro_best = copy.deepcopy(numpy_density_ferro)
        numpy_density_antiferro_best = copy.deepcopy(numpy_density_antiferro)
        delta_chi_sq_f_best = copy.deepcopy(delta_chi_sq_f)
        delta_chi_sq_a_best = copy.deepcopy(delta_chi_sq_a)

        if chi_sq_n < gof_desired:
            print(f"OUT: cycle {i_cycle:5} chi_sq/n is less than \
{gof_desired:.2f}",
                  end="\r")
            if flag_info:
                d_info["print"] = f"""OUT on iteration {i_cycle:5}:
    chi_sq/n is less than {gof_desired:.2f}."""
            break
        elif ((c_lambda < c_lambda_min) & False):
            print(f"OUT: cycle {i_cycle:5} chi_sq/n is {chi_sq_n:.2f}. \
c_lambda: {c_lambda:}",
                  end="\r")
            if flag_info:
                d_info["print"] = f"""OUT on iteration {i_cycle:5}:
    chi_sq/n is {chi_sq_n:.2f}. c_lambda less minimal"""
            if not(flag_two_channel):
                density_point.numpy_density = numpy_density_best
            density_point.numpy_density_ferro = numpy_density_ferro_best
            density_point.numpy_density_antiferro = \
                numpy_density_antiferro_best
            density_point.renormalize_numpy_densities(
                flag_two_channel=flag_two_channel)
            break
        elif ((chi_sq_n_diff < 0.001) and (i_cycle > 10)):
            print(f"OUT: cycle {i_cycle:5} chi_sq/n is {chi_sq_n:.2f} as diff \
of GoF is less than 0.001.",
                  end="\r")
            if flag_info:
                d_info["print"] = f"""OUT on iteration {i_cycle:5}:
    chi_sq/n is {chi_sq_n:.2f};
    difference of GoF is less than 0.001."""
            break
        elif flag_info:
            if d_info["stop"]:
                d_info["stop"] = False
                break

        delta_chi_sq_f_mult_i = delta_chi_sq_f/mult_i
        delta_chi_sq_a_mult_i = delta_chi_sq_a/mult_i

        rel_diff = 0.05
        if not(flag_two_channel):
            delta_chi_sq_mult_i = delta_chi_sq/mult_i
            c_lambda = choose_max_clambda(c_lambda, numpy_density,
                                          delta_chi_sq_mult_i, rel_diff)
            density_point.numpy_density = numpy_density*numpy.exp(
                -c_lambda*delta_chi_sq_mult_i)
        else:
            c_lambda = choose_max_clambda(c_lambda, numpy_density_ferro,
                                          delta_chi_sq_f_mult_i, rel_diff)

        density_point.numpy_density_ferro = \
            numpy_density_ferro*numpy.exp(-c_lambda*delta_chi_sq_f_mult_i)
        density_point.numpy_density_antiferro = \
            numpy_density_antiferro*numpy.exp(-c_lambda*delta_chi_sq_a_mult_i)

        density_point.renormalize_numpy_densities(
            flag_two_channel=flag_two_channel)

    density_point.numpy_to_items()
    for diffrn in l_diffrn:
        diffrn.diffrn_refln.numpy_to_items()
        chi_sq, points = diffrn.diffrn_refln.calc_chi_sq_points()
        refine_ls = RefineLs(goodness_of_fit_all=chi_sq/points,
                             number_reflns=points)
        diffrn.add_items([refine_ls])
    return density_point


def refine_susceptibility(crystal: Crystal, l_diffrn: List[Diffrn],
                          density_point: DensityPointL,
                          mem_parameters: MEMParameters, disp: bool = True,
                          d_info: dict = None) -> (float, float):
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
    flag_info = d_info is not None

    if flag_info:
        d_info_keys = d_info.keys()
        if "stop" not in d_info_keys:
            d_info["stop"] = False
        if "print" not in d_info_keys:
            d_info["print"] = ""

    crystal.apply_constraints()

    flag_two_channel = mem_parameters.method == "2channel"

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
        u_matrix = diffrn_orient_matrix.u
        e_up = calc_e_up_loc(0., 0., 0., u_matrix)

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
        index_hkl = numpy.stack([index_h, index_k, index_l], dtype=bool)
        f_nucl = crystal.calc_f_nucl(index_hkl)
        k_hkl = cell.calc_k_loc(*hkl)
        phase_3d = density_point.calc_phase_3d(hkl, space_group_symop)

        moment_2d, chi_2d_ferro, chi_2d_antiferro = \
            density_point.calc_moment_2d(
                space_group_symop, cell, atom_site_susceptibility, h_loc,
                chi_iso_ferro=1., chi_iso_antiferro=1.,
                flag_two_channel=flag_two_channel)

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

    l_name_1 = atom_site_susceptibility.get_variable_names()
    l_par_1_0 = [atom_site_susceptibility.get_variable_by_name(name) for name
                 in l_name_1]

    l_name_2 = mem_parameters.get_variable_names()
    l_par_2_0 = [mem_parameters.get_variable_by_name(name) for name
                 in l_name_2]
    l_name = l_name_1 + l_name_2
    l_par_0 = l_par_1_0 + l_par_2_0

    if len(l_name) == 0:
        return

    def temp_func(l_par):
        for name, parameter in zip(l_name, l_par):
            if name[0][0] == "atom_site_susceptibility":
                atom_site_susceptibility.set_variable_by_name(name, parameter)
            elif name[0][0] == "mem_parameters":
                mem_parameters.set_variable_by_name(name, parameter)

        chi_iso_f = mem_parameters.chi_ferro
        chi_iso_af = mem_parameters.chi_antiferro

        atom_site_susceptibility.apply_chi_iso_constraint(cell)
        atom_site_susceptibility.apply_space_group_constraint(
            atom_site, space_group, cell)

        l_chi_sq = []
        l_der_chi_sq, l_der_chi_sq_f, l_der_chi_sq_a = [], [], []

        # FIXME: add flag for orbital magnetic moment
        for diffrn, f_nucl, fr_e, fr_s, e_up, h_loc, k_hkl, phase_3d, \
            chi_perp_ferro, chi_perp_aferro in \
            zip(l_diffrn, l_f_nucl, l_fr_e, l_fr_s, l_e_up, l_h_loc, l_k_hkl,
                l_phase_3d, l_chi_perp_ferro, l_chi_perp_antiferro):

            moment_2d, moment_ferro, moment_antiferro = \
                density_point.calc_moment_2d(
                    space_group_symop, cell, atom_site_susceptibility, h_loc,
                    chi_iso_ferro=1., chi_iso_antiferro=1.,
                    flag_two_channel=flag_two_channel)

            f_m = calc_fm_by_density(mult_i, den_i, np, volume, moment_2d,
                                     phase_3d)
            f_m_perp = calc_moment_perp(k_hkl, f_m)

            # # correction on orbital moment
            # diffrn_refln = diffrn.diffrn_refln
            # setup = diffrn.setup
            # field = setup.field
            # ind_h = diffrn_refln.numpy_index_h
            # ind_k = diffrn_refln.numpy_index_k
            # ind_l = diffrn_refln.numpy_index_l
            # chi_m = crystal.calc_susceptibility_moment_tensor(
            #     ind_h, ind_k, ind_l, flag_only_orbital=True)
            # sft_ij = chi_m[:9]
            # sftm_ij = chi_m[9:]

            # fm_orb_perp_loc = calc_fm_perp_loc(e_up, field, k_hkl, sft_ij,
            #                                    sftm_ij)

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
    if flag_info:
        d_info["print"] = \
            f"Chi_sq/n before optimization {chi_sq/total_peaks:.5f}."

    res = scipy.optimize.minimize(
        temp_func, l_par_0, method="BFGS",
        callback=lambda x: func_temp(x, param_name=l_name, d_info=d_info),
        options={"eps": 0.001})

    l_param = res.x
    chi_sq_new = res.fun

    hess_inv = res["hess_inv"]
    sigma = (abs(numpy.diag(hess_inv)))**0.5

    print(f"Chi_sq after optimization {chi_sq_new/total_peaks:.5f}.          ",
          end="\r")
    if flag_info:
        d_info["print"] = \
            f"Chi_sq/n after optimization {chi_sq_new/total_peaks:.5f}."

    for name, parameter, sig in zip(l_name, l_param, sigma):
        name_sig = name[:-1] + ((f"{name[-1][0]:}_sigma", name[-1][1]), )
        if name[0][0] == "atom_site_susceptibility":
            atom_site_susceptibility.set_variable_by_name(name, parameter)
            atom_site_susceptibility.set_variable_by_name(name_sig, sig)
        elif name[0][0] == "mem_parameters":
            mem_parameters.set_variable_by_name(name, parameter)
            mem_parameters.set_variable_by_name(name_sig, sig)

    for diffrn in l_diffrn:
        diffrn.diffrn_refln.numpy_to_items()
        chi_sq, points = diffrn.diffrn_refln.calc_chi_sq_points()
        refine_ls = RefineLs(goodness_of_fit_all=chi_sq/points,
                             number_reflns=points)
        diffrn.add_items([refine_ls])


def func_temp(*argv, param_name: List[tuple] = None, d_info: dict = None) \
        -> bool:
    """Show data at optimization procedure."""
    flag_out = False
    l_param = argv[0]
    s_out = " ".join([f"{_:10.5f}" for _ in l_param])
    print(f"{s_out:}", end="\r")
    if d_info is not None:
        if param_name is not None:
            ls_out = [f"{name[-1][0]:}: {param: 10.5f}" for name, param in
                      zip(param_name, l_param)]
        else:
            ls_out = [f"{param: 10.5f}" for param in l_param]
        d_info["print"] = "Best solution:\n\n"+"\n".join(ls_out)
        flag_out = d_info["stop"]
    return flag_out


def make_cycle(crystal: Crystal, l_diffrn: List[Diffrn],
               mem_parameters: MEMParameters,
               c_lambda: float = 1e-6, n_iterations: int = 51,
               n_cycle: int = 10, disp: bool = True,
               d_info: dict = None):
    """Rho - Chi cycle."""
    flag_info = d_info is not None

    d_info_2 = None
    if flag_info:
        d_info_keys = d_info.keys()
        if "stop" not in d_info_keys:
            d_info["stop"] = False
        if "print" not in d_info_keys:
            d_info["print"] = ""
        if "d_info" not in d_info_keys:
            d_info_2 = {"stop": False, "print": ""}
            d_info["d_info"] = d_info_2

    flag_info = d_info is not None

    for i_cycle in range(n_cycle):
        print(f"Cycle {(i_cycle+1):4}. Entropy maximization                  ",
              end="\r")

        if flag_info:
            d_info["print"] = f"""Cycle {(i_cycle+1):4}/{n_cycle:}:
    Entropy maximization."""
            if d_info["stop"]:
                d_info["stop"] = False
                break
        density_point = maximize_entropy(
            crystal, l_diffrn, mem_parameters, c_lambda=c_lambda,
            n_iterations=n_iterations, disp=disp, d_info=d_info_2)

        if flag_info:
            d_info["print"] = f"""Cycle {(i_cycle+1):4}/{n_cycle:}:
    Chi refinement."""
            if d_info["stop"]:
                d_info["stop"] = False
                break
        print(f"Cycle {(i_cycle+1):4}. Chi refinement                        ",
              end="\r")
        refine_susceptibility(crystal, l_diffrn, density_point, mem_parameters,
                              disp=disp, d_info=d_info_2)
    return density_point


def calc_moments(fract_x, fract_y, fract_z, field_loc,
                 density_point: DensityPointL, crystal: Crystal,
                 mem_parameters: MEMParameters):
    """Calculate magnetic moments in points fract_x, fract_y, fract_z."""
    points_a = mem_parameters.points_a
    points_b = mem_parameters.points_b
    points_c = mem_parameters.points_c

    cell = crystal.cell
    atom_site = crystal.atom_site

    np_ind_xyz = numpy.transpose(numpy.array([fract_x*points_a, fract_y*points_b, fract_z*points_c], dtype=float))

    # np_ind_xyz = (fract_x*points_a, fract_y*points_b, fract_z*points_c)

    chi_3d_11, chi_3d_22, chi_3d_33, chi_3d_12, chi_3d_13, chi_3d_23, \
        chi_3d_iso = calc_densities_3d_for_mem(
        density_point, crystal, mem_parameters=mem_parameters)

    c_11 = tri_linear_interpolation(np_ind_xyz, chi_3d_11)
    c_22 = tri_linear_interpolation(np_ind_xyz, chi_3d_22)
    c_33 = tri_linear_interpolation(np_ind_xyz, chi_3d_33)
    c_12 = tri_linear_interpolation(np_ind_xyz, chi_3d_12)
    c_13 = tri_linear_interpolation(np_ind_xyz, chi_3d_13)
    c_23 = tri_linear_interpolation(np_ind_xyz, chi_3d_23)

    c_iso = tri_linear_interpolation(np_ind_xyz, chi_3d_iso)

    c_orto = cell.ortogonalize_matrix((c_11, c_12, c_13, c_12, c_22, c_23,
                                       c_13, c_23, c_33))
    m_x, m_y, m_z = calc_product_matrix_vector(c_orto, field_loc)
    m_b_x, m_b_y, m_b_z = c_iso*field_loc[0], c_iso*field_loc[1], \
        c_iso*field_loc[2]

    moment_x, moment_y, moment_z = m_x+m_b_x, m_y+m_b_y, m_z+m_b_z

    return moment_x, moment_y, moment_z


def calc_moments_in_unit_cell(density_point: DensityPointL, crystal: Crystal,
                              mem_parameters: MEMParameters, field_loc):
    """Calculate magnetization density in an unit cell."""
    points_a = mem_parameters.points_a
    points_b = mem_parameters.points_b
    points_c = mem_parameters.points_c

    h_x = numpy.linspace(0., 1., num=points_a, endpoint=False, dtype=float)
    h_y = numpy.linspace(0., 1., num=points_b, endpoint=False, dtype=float)
    h_z = numpy.linspace(0., 1., num=points_c, endpoint=False, dtype=float)

    fract_x, fract_y, fract_z = numpy.meshgrid(h_x, h_y, h_z, indexing="ij")
    fract_x, fract_y = fract_x.flatten(), fract_y.flatten()
    fract_z = fract_z.flatten()

    moment_x, moment_y, moment_z = calc_moments(
        fract_x, fract_y, fract_z, field_loc, density_point, crystal,
        mem_parameters)

    return fract_x, fract_y, fract_z, moment_x, moment_y, moment_z


def calc_section_for_mem(section: Section, density_point: DensityPointL,
                         crystal: Crystal, mem_parameters: MEMParameters):
    """
    Calculate magnetization density of paramagnetic compound at given field.

    The result is written into file.
    """
    field_loc = numpy.array([section.field_x, section.field_y,
                             section.field_z], dtype=float)
    cell = crystal.cell
    atom_site = crystal.atom_site
    np_fract_x, np_fract_y, np_fract_z = section.calc_fractions(
        cell, atom_site)

    moment_x, moment_y, moment_z = calc_moments(
        np_fract_x, np_fract_y, np_fract_z, field_loc, density_point, crystal,
        mem_parameters)

    size_x, size_y = float(section.size_x), float(section.size_y)

    np_x = (numpy.array(range(section.points_x), dtype=float) -
            0.5*section.points_x) * size_x / float(section.points_x)

    np_y = (numpy.array(range(section.points_y), dtype=float) -
            0.5*section.points_y) * size_y / float(section.points_y)

    np_x_2d, np_y_2d = numpy.meshgrid(np_x, np_y, indexing="ij")
    np_x, np_y = np_x_2d.flatten(), np_y_2d.flatten()

    v_pos_x, v_pos_y, v_pos_z = section.calc_axes_x_y_z(cell, atom_site)

    ls_out = ["loop_\n_section_cut_dist_x\n_section_cut_dist_y"]
    ls_out.append("_section_cut_moment_x\n_section_cut_moment_y")
    ls_out.append("_section_cut_moment_z")
    for _x, _y, v_m_x, v_m_y, v_m_z in zip(
            np_x, np_y, moment_x, moment_y, moment_z):
        val_x = v_pos_x[0]*v_m_x + v_pos_x[1]*v_m_y + v_pos_x[2]*v_m_z
        val_y = v_pos_y[0]*v_m_x + v_pos_y[1]*v_m_y + v_pos_y[2]*v_m_z
        val_z = v_pos_z[0]*v_m_x + v_pos_z[1]*v_m_y + v_pos_z[2]*v_m_z
        ls_out.append(f"{_x:9.5f} {_y:9.5f} {val_x:15.10f} {val_y:15.10f} \
{val_z:15.10f}")
    with open(section.url_out, "w") as fid:
        fid.write("\n".join(ls_out))


def calc_densities_3d_for_mem(density_point: DensityPointL, crystal: Crystal,
                              mem_parameters: MEMParameters):
    """
    Calculate 3d density for susceptibility tensor and isotropical background.

    Output:

        - chi_den_3d_11, chi_den_3d_22, chi_den_3d_33,
          chi_den_3d_12, chi_den_3d_13, chi_den_3d_23,
          chi_den_3d_iso
    """
    points_a = mem_parameters.points_a
    points_b = mem_parameters.points_b
    points_c = mem_parameters.points_c
    chi_iso_ferro = mem_parameters.chi_ferro
    chi_iso_antiferro = mem_parameters.chi_antiferro

    atom_site_susceptibility = crystal.atom_site_susceptibility

    space_group = crystal.space_group
    full_space_group_symop = space_group.full_space_group_symop

    points_abc = (points_a, points_b, points_c)

    density_point.calc_rbs_i(full_space_group_symop, points_a=points_a,
                             points_b=points_b, points_c=points_c)

    susc_11, susc_22, susc_33, susc_12, susc_13, susc_23 = \
        density_point.calc_susc_i(atom_site_susceptibility)

    index_xyz = (numpy.array(density_point.index_x, dtype=int),
                 numpy.array(density_point.index_y, dtype=int),
                 numpy.array(density_point.index_z, dtype=int))

    den_i = numpy.array(density_point.density, dtype=float)
    den_ferro_i = numpy.array(density_point.density_ferro, dtype=float)
    den_aferro_i = numpy.array(density_point.density_antiferro, dtype=float)

    r_11 = full_space_group_symop.numpy_r_11.astype(float)
    r_12 = full_space_group_symop.numpy_r_12.astype(float)
    r_13 = full_space_group_symop.numpy_r_13.astype(float)
    r_21 = full_space_group_symop.numpy_r_21.astype(float)
    r_22 = full_space_group_symop.numpy_r_22.astype(float)
    r_23 = full_space_group_symop.numpy_r_23.astype(float)
    r_31 = full_space_group_symop.numpy_r_31.astype(float)
    r_32 = full_space_group_symop.numpy_r_32.astype(float)
    r_33 = full_space_group_symop.numpy_r_33.astype(float)

    b_1 = full_space_group_symop.numpy_b_1.astype(float)
    b_2 = full_space_group_symop.numpy_b_2.astype(float)
    b_3 = full_space_group_symop.numpy_b_3.astype(float)

    r_ij = (r_11, r_12, r_13, r_21, r_22, r_23, r_31, r_32, r_33)
    b_i = (b_1, b_2, b_3)

    den_3d = transfer_to_density_3d(index_xyz, den_i, points_abc, r_ij, b_i)
    den_ferro_3d = transfer_to_density_3d(
        index_xyz, den_ferro_i, points_abc, r_ij, b_i)
    den_aferro_3d = transfer_to_density_3d(
        index_xyz, den_aferro_i, points_abc, r_ij, b_i)

    chi_3d_11, chi_3d_22, chi_3d_33, chi_3d_12, chi_3d_13, chi_3d_23 = \
        transfer_to_chi_3d(index_xyz, susc_11, susc_22, susc_33, susc_12,
                           susc_13, susc_23, points_abc, r_ij, b_i)

    chi_den_3d_11, chi_den_3d_22 = chi_3d_11*den_3d, chi_3d_22*den_3d
    chi_den_3d_33, chi_den_3d_12 = chi_3d_33*den_3d, chi_3d_12*den_3d
    chi_den_3d_13, chi_den_3d_23 = chi_3d_13*den_3d, chi_3d_23*den_3d

    chi_den_3d_iso = chi_iso_ferro*den_ferro_3d + \
        chi_iso_antiferro*den_aferro_3d

    return chi_den_3d_11, chi_den_3d_22, chi_den_3d_33,\
        chi_den_3d_12, chi_den_3d_13, chi_den_3d_23, chi_den_3d_iso
