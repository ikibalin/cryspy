# -*- coding: utf-8 -*-
"""The MEM algorithm is described."""
from typing import List
import numpy

from cryspy.A_functions_base.function_1_matrices import calc_chi_sq
from cryspy.C_item_loop_classes.cl_3_density_point import DensityPointL
from cryspy.E_data_classes.cl_1_crystal import Crystal
from cryspy.E_data_classes.cl_2_diffrn import Diffrn


def maximize_entropy(crystal: Crystal, l_diffrn: List[Diffrn],
                     c_lambda: float = 1e-7, n_cycle: int = 100,
                     chi_iso_ferro: float = 0., chi_iso_antiferro: float = 0.,
                     n_x: int = 48, n_y: int = 48, n_z: int = 48,
                     prior_density: str = "uniform", disp: bool = True) -> \
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
        DESCRIPTION. The default is 100.
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

    for i_cycle in range(n_cycle):
        numpy_density = density_point.numpy_density
        numpy_density_ferro = density_point.numpy_density_ferro
        numpy_density_antiferro = density_point.numpy_density_antiferro
        chi_sq, delta_chi_sq, delta_chi_sq_f, delta_chi_sq_a = temp_func()
        if disp:
            chi_sq_n = chi_sq/float(total_peaks)
            print(f"cycle {i_cycle:5}. chi_sq/n: {chi_sq_n:.2f}       ",
                  end="\r")
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
