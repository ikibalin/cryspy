# -*- coding: utf-8 -*-
"""
Magnetic structure factors: structure factor tensor, ordered magnetic moments
and bulk susceptibility.

Split from structure_factor.py, which still exports these names.
"""

import numpy

from .structure_factor import (
    get_atom_symm_elems_by_atom_fract_xyz,
    calc_pr1,
    calc_pr2,
    calc_pr3,
    calc_pr4,
    calc_pr5,
    calc_f_by_f_asym_a_pr,
    calc_sft_ccs_asym_a_by_pr,
)

from .matrix_operations import (
    calc_det_m,
    calc_m1_m2_inv_m1,
    calc_m_v,
    calc_vector_product_v1_v2_v1,
    calc_m_q_inv_m,
)
from .unit_cell import (
    calc_eq_ccs_by_unit_cell_parameters,
    calc_m_m_by_unit_cell_parameters,
    calc_m_m_norm_by_unit_cell_parameters,
    calc_sthovl_by_unit_cell_parameters,
)
from .debye_waller_factor import calc_dwf
from .symmetry_elements import (
    calc_multiplicity_by_atom_symm_elems,
    calc_full_symm_elems_by_reduced,
)
from .magnetic_form_factor import (
    calc_form_factor,
    calc_magnetic_form_factor_by_multipole_model,
)

na = numpy.newaxis


def calc_sft_ccs_by_dictionary(
    dict_crystal, dict_in_out, flag_use_precalculated_data: bool = False
):
    """Calculate structure factor tensor in CCS (X||a*, Z||c) based on the information given in dictionary.
    Output information is written in the same dictionary.
    """
    dict_crystal_keys = dict_crystal.keys()
    dict_in_out_keys = dict_in_out.keys()
    necessary_crystal_keys = set(
        [
            "unit_cell_parameters",
        ]
    )
    diff_set_crystal = necessary_crystal_keys.difference(
        set(dict_crystal_keys)
    )
    if len(diff_set_crystal) != 0:
        raise AttributeError(
            f"The following attributes have to be defined {diff_set_crystal:}"
        )

    flag_reduced_symm_elems = (
        len(
            set(
                ["reduced_symm_elems", "centrosymmetry", "translation_elems"]
            ).difference(set(dict_crystal_keys))
        )
        == 0
    )
    flag_full_symm_elems = (
        len(
            set(
                [
                    "full_symm_elems",
                ]
            ).difference(set(dict_crystal_keys))
        )
        == 0
    )
    flag_full_mcif_elems = (
        len(
            set(
                [
                    "full_mcif_elems",
                ]
            ).difference(set(dict_crystal_keys))
        )
        == 0
    )

    if not (
        flag_reduced_symm_elems or flag_full_symm_elems or flag_full_mcif_elems
    ):
        raise AttributeError("The symmetry elements have to be defined.")

    necessary_in_out_keys = set(
        [
            "index_hkl",
        ]
    )
    diff_set_in_out = necessary_in_out_keys.difference(set(dict_in_out_keys))
    if len(diff_set_in_out) != 0:
        raise AttributeError(
            f"The following attributes have to be defined {diff_set_in_out:}"
        )

    index_hkl = dict_in_out["index_hkl"]

    non_zero_keys = set(
        [
            "mag_atom_lande_factor",
            "mag_atom_kappa",
            "mag_atom_j0_parameters",
            "mag_atom_j2_parameters",
        ]
    )
    diff_set_crystal = non_zero_keys.difference(set(dict_crystal_keys))

    if len(diff_set_crystal) != 0:
        sft_ccs = numpy.zeros((9, index_hkl.shape[-1]), dtype=complex)
        dder = {}
        return sft_ccs, dder

    if "flag_only_orbital" in dict_in_out_keys:
        flag_only_orbital = dict_in_out["flag_only_orbital"]
    else:
        flag_only_orbital = False

    if flag_reduced_symm_elems:
        reduced_symm_elems = dict_crystal["reduced_symm_elems"]
        centrosymmetry = dict_crystal["centrosymmetry"]
        if centrosymmetry:
            centrosymmetry_position = dict_crystal["centrosymmetry_position"]
        else:
            centrosymmetry_position = None
        translation_elems = dict_crystal["translation_elems"]
    elif flag_full_symm_elems:
        full_symm_elems = dict_crystal["full_symm_elems"]
        reduced_symm_elems = full_symm_elems
        centrosymmetry = False
        centrosymmetry_position = None
        translation_elems = numpy.array([[0], [0], [0], [1]], dtype=int)
    elif flag_full_mcif_elems:
        full_mcif_elems = dict_crystal["full_mcif_elems"]
        reduced_symm_elems = full_mcif_elems
        centrosymmetry = False
        centrosymmetry_position = None
        translation_elems = numpy.array([[0], [0], [0], [1]], dtype=int)

    unit_cell_parameters = dict_crystal["unit_cell_parameters"]
    atom_para_index = dict_crystal["atom_para_index"]
    atom_para_fract_xyz = dict_crystal["atom_fract_xyz"][:, atom_para_index]
    atom_para_sc_fract = dict_crystal["atom_site_sc_fract"][:, atom_para_index]
    atom_para_sc_b = dict_crystal["atom_site_sc_b"][:, atom_para_index]
    atom_para_fract_xyz = (
        calc_m_v(
            atom_para_sc_fract,
            numpy.mod(atom_para_fract_xyz, 1),
            flag_m=False,
            flag_v=False,
        )[0]
        + atom_para_sc_b
    )

    atom_para_occupancy = dict_crystal["atom_occupancy"][atom_para_index]
    atom_para_b_iso = dict_crystal["atom_b_iso"][atom_para_index]
    atom_beta = dict_crystal["atom_beta"]
    if "atom_site_aniso_sc_beta" in dict_crystal_keys:
        atom_site_aniso_sc_beta = dict_crystal["atom_site_aniso_sc_beta"]
        atom_site_aniso_index = dict_crystal["atom_site_aniso_index"]
        atom_sc_beta = numpy.zeros((6,) + atom_beta.shape, dtype=float)
        atom_sc_beta[:, :, atom_site_aniso_index] = atom_site_aniso_sc_beta
        atom_beta = (atom_sc_beta * numpy.expand_dims(atom_beta, axis=0)).sum(
            axis=1
        )
    atom_para_beta = atom_beta[:, atom_para_index]

    mag_atom_para_index = dict_crystal["mag_atom_para_index"]
    atom_para_lande_factor = dict_crystal["mag_atom_lande_factor"][
        mag_atom_para_index
    ]
    atom_para_kappa = dict_crystal["mag_atom_kappa"][mag_atom_para_index]
    atom_para_j0_parameters = dict_crystal["mag_atom_j0_parameters"][
        :, mag_atom_para_index
    ]
    atom_para_j2_parameters = dict_crystal["mag_atom_j2_parameters"][
        :, mag_atom_para_index
    ]

    atom_para_susceptibility = dict_crystal["atom_para_susceptibility"]
    atom_para_sc_chi = dict_crystal["atom_para_sc_chi"]
    flag_unit_cell_parameters = numpy.any(
        dict_crystal["flags_unit_cell_parameters"]
    )
    flag_atom_para_fract_xyz = numpy.any(
        dict_crystal["flags_atom_fract_xyz"][:, atom_para_index]
    )
    flag_atom_para_occupancy = numpy.any(
        dict_crystal["flags_atom_occupancy"][atom_para_index]
    )
    flag_atom_para_b_iso = numpy.any(
        dict_crystal["flags_atom_b_iso"][atom_para_index]
    )
    flag_atom_para_beta = numpy.any(
        dict_crystal["flags_atom_beta"][:, atom_para_index]
    )
    flag_atom_para_susceptibility = numpy.any(
        dict_crystal["flags_atom_para_susceptibility"]
    )
    flag_atom_para_lande_factor = numpy.any(
        dict_crystal["flags_mag_atom_lande_factor"][mag_atom_para_index]
    )
    flag_atom_para_kappa = numpy.any(
        dict_crystal["flags_mag_atom_kappa"][mag_atom_para_index]
    )

    sft_ccs, dder = calc_sft_ccs(
        index_hkl,
        reduced_symm_elems,
        centrosymmetry,
        centrosymmetry_position,
        translation_elems,
        unit_cell_parameters,
        atom_para_fract_xyz,
        atom_para_occupancy,
        atom_para_susceptibility,
        atom_para_b_iso,
        atom_para_beta,
        atom_para_lande_factor,
        atom_para_kappa,
        atom_para_j0_parameters,
        atom_para_j2_parameters,
        atom_para_sc_chi,
        dict_in_out=dict_in_out,
        flag_only_orbital=flag_only_orbital,
        flag_unit_cell_parameters=flag_unit_cell_parameters,
        flag_atom_para_fract_xyz=flag_atom_para_fract_xyz,
        flag_atom_para_occupancy=flag_atom_para_occupancy,
        flag_atom_para_susceptibility=flag_atom_para_susceptibility,
        flag_atom_para_b_iso=flag_atom_para_b_iso,
        flag_atom_para_beta=flag_atom_para_beta,
        flag_atom_para_lande_factor=flag_atom_para_lande_factor,
        flag_atom_para_kappa=flag_atom_para_kappa,
        flag_use_precalculated_data=flag_use_precalculated_data,
    )
    return sft_ccs, dder


def calc_sft_ccs(
    index_hkl,
    reduced_symm_elems,
    centrosymmetry,
    centrosymmetry_position,
    translation_elems,
    unit_cell_parameters,
    atom_para_fract_xyz,
    atom_para_occupancy,
    atom_para_susceptibility,
    atom_para_b_iso,
    atom_para_beta,
    atom_para_lande_factor,
    atom_para_kappa,
    atom_para_j0_parameters,
    atom_para_j2_parameters,
    atom_para_sc_chi,
    dict_in_out: dict = None,
    flag_only_orbital: bool = False,
    flag_unit_cell_parameters: bool = False,
    flag_atom_para_fract_xyz: bool = False,
    flag_atom_para_occupancy: bool = False,
    flag_atom_para_susceptibility: bool = False,
    flag_atom_para_b_iso: bool = False,
    flag_atom_para_beta: bool = False,
    flag_atom_para_lande_factor: bool = False,
    flag_atom_para_kappa: bool = False,
    flag_use_precalculated_data: bool = False,
):
    """Calculate structure factor tensor in Cartesian coordinate system with X||a*, Z||c in 10**-12 cm.
    Note, that the susceptibility parameters are given in mu_B.
    """

    if dict_in_out is None:
        flag_dict = False
        dict_in_out_keys = []
    else:
        flag_dict = True
        dict_in_out_keys = dict_in_out.keys()
        if "index_hkl" in dict_in_out_keys:
            if numpy.any(dict_in_out["index_hkl"] != index_hkl):
                dict_in_out.clear()
                dict_in_out["index_hkl"] = index_hkl

    if flag_use_precalculated_data and (
        "atom_para_multiplicity" in dict_in_out_keys
    ):
        mag_atom_multiplicity = dict_in_out["atom_para_multiplicity"]
    else:
        atom_symm_elems = get_atom_symm_elems_by_atom_fract_xyz(
            atom_para_fract_xyz
        )

        if "full_symm_elems" in dict_in_out_keys:
            full_symm_elems = dict_in_out["full_symm_elems"]
        else:
            full_symm_elems = calc_full_symm_elems_by_reduced(
                reduced_symm_elems[:13],
                centrosymmetry,
                centrosymmetry_position,
                translation_elems,
            )
            if flag_dict:
                dict_in_out["full_symm_elems"] = full_symm_elems
        mag_atom_multiplicity = calc_multiplicity_by_atom_symm_elems(
            full_symm_elems, atom_symm_elems
        )
        if flag_dict:
            dict_in_out["atom_para_multiplicity"] = mag_atom_multiplicity

    flag_pr_1 = flag_atom_para_fract_xyz
    if (
        flag_use_precalculated_data
        and ("pr_1_atom_para" in dict_in_out_keys)
        and not (flag_atom_para_fract_xyz)
    ):
        pr_1 = dict_in_out["pr_1_atom_para"]
    else:
        pr_1, dder_pr_1 = calc_pr1(
            index_hkl,
            reduced_symm_elems[:13],
            atom_para_fract_xyz,
            flag_fract_xyz=flag_atom_para_fract_xyz,
        )
        if flag_dict:
            dict_in_out["pr_1_atom_para"] = pr_1

    if flag_use_precalculated_data and ("pr_2" in dict_in_out_keys):
        pr_2 = dict_in_out["pr_2"]
    else:
        pr_2 = calc_pr2(index_hkl, reduced_symm_elems[:13])
        if flag_dict:
            dict_in_out["pr_2"] = pr_2

    if flag_use_precalculated_data and ("pr_3" in dict_in_out_keys):
        pr_3 = dict_in_out["pr_3"]
    else:
        pr_3 = calc_pr3(index_hkl, translation_elems)
        if flag_dict:
            dict_in_out["pr_3"] = pr_3

    if flag_use_precalculated_data and ("pr_4" in dict_in_out_keys):
        pr_4 = dict_in_out["pr_4"]
    else:
        pr_4 = calc_pr4(index_hkl, centrosymmetry_position)
        if flag_dict:
            dict_in_out["pr_4"] = pr_4

    flag_sthovl = flag_unit_cell_parameters
    if (
        flag_use_precalculated_data
        and ("sthovl" in dict_in_out_keys)
        and not (flag_sthovl)
    ):
        sthovl = dict_in_out["sthovl"]
    else:
        sthovl, dder_sthovl = calc_sthovl_by_unit_cell_parameters(
            index_hkl,
            unit_cell_parameters,
            flag_unit_cell_parameters=flag_unit_cell_parameters,
        )
        if flag_dict:
            dict_in_out["sthovl"] = sthovl

    flag_pr_5 = flag_unit_cell_parameters
    if (
        flag_use_precalculated_data
        and ("pr_5" in dict_in_out_keys)
        and not (flag_pr_5)
    ):
        pr_5 = dict_in_out["pr_5"]
    else:
        pr_5, dder_pr_5 = calc_pr5(
            reduced_symm_elems[:13],
            unit_cell_parameters,
            flag_unit_cell_parameters=flag_unit_cell_parameters,
        )
        if flag_dict:
            dict_in_out["pr_5"] = pr_5

    flag_atom_para_form_factor = (
        flag_sthovl or flag_atom_para_lande_factor or flag_atom_para_kappa
    )
    flag_hh = True
    if "flag_only_orbital" in dict_in_out_keys:
        flag_hh = flag_only_orbital == dict_in_out["flag_only_orbital"]
    dict_in_out["flag_only_orbital"] = flag_only_orbital
    if (
        flag_use_precalculated_data
        and ("atom_para_form_factor" in dict_in_out_keys)
        and not (flag_atom_para_form_factor)
        and flag_hh
    ):
        atom_para_form_factor = dict_in_out["atom_para_form_factor"]
    else:
        atom_para_form_factor, dder_ff = calc_form_factor(
            sthovl[:, na],
            atom_para_lande_factor[na, :],
            atom_para_kappa[na, :],
            atom_para_j0_parameters[:, na, :],
            atom_para_j2_parameters[:, na, :],
            flag_lande_factor=flag_atom_para_lande_factor,
            flag_only_orbital=flag_only_orbital,
            flag_sthovl=flag_sthovl,
            flag_kappa=flag_atom_para_kappa,
        )
        if flag_dict:
            dict_in_out["atom_para_form_factor"] = atom_para_form_factor

    # dimensions ["hkl", "reduced symmetry", "atom"]
    flag_debye_waller_factor = (
        flag_sthovl or flag_atom_para_b_iso or flag_atom_para_beta
    )
    if (
        flag_use_precalculated_data
        and ("atom_para_debye_waller_factor" in dict_in_out_keys)
        and not (flag_debye_waller_factor)
    ):
        debye_waller_factor = dict_in_out["atom_para_debye_waller_factor"]
    else:
        debye_waller_factor, dder_dw = calc_dwf(
            index_hkl[:, :, na, na],
            sthovl[:, na, na],
            atom_para_b_iso[na, na, :],
            atom_para_beta[:, na, na, :],
            reduced_symm_elems[4:13, na, :, na],
            flag_sthovl=flag_sthovl,
            flag_b_iso=flag_atom_para_b_iso,
            flag_beta=flag_atom_para_beta,
        )
        if flag_dict:
            dict_in_out["atom_para_debye_waller_factor"] = debye_waller_factor

    flag_debye_waller = flag_atom_para_b_iso or flag_atom_para_beta

    flag_sft_ccs_asym = (
        flag_atom_para_form_factor
        or flag_debye_waller
        or flag_atom_para_occupancy
        or flag_atom_para_susceptibility
        or flag_pr_1
        or flag_pr_5
    )
    if (
        flag_use_precalculated_data
        and ("sft_ccs_asym" in dict_in_out_keys)
        and not (flag_sft_ccs_asym)
    ):
        sft_ccs_asym = dict_in_out["sft_ccs_asym"]
    else:
        sft_ccs_asym, dder_sft_ccs_asym = calc_sft_ccs_asym_a_by_pr(
            mag_atom_multiplicity,
            debye_waller_factor,
            atom_para_occupancy,
            atom_para_susceptibility,
            atom_para_sc_chi,
            pr_1,
            pr_2,
            pr_5,  # theta=theta,
            flag_debye_waller=flag_debye_waller,
            flag_atom_para_occupancy=flag_atom_para_occupancy,
            flag_atom_para_susceptibility=flag_atom_para_susceptibility,
            flag_pr_1=flag_pr_1,
            flag_pr_5=flag_pr_5,
        )
        if flag_dict:
            dict_in_out["sft_ccs_asym"] = sft_ccs_asym

    flag_sft_ccs = flag_sft_ccs_asym
    if (
        flag_use_precalculated_data
        and ("sft_ccs" in dict_in_out_keys)
        and not (flag_sft_ccs)
    ):
        sft_ccs = dict_in_out["sft_ccs"]
    else:
        sft_ccs, dder_sft_ccs = calc_f_by_f_asym_a_pr(
            sft_ccs_asym,
            atom_para_form_factor,
            pr_3,
            centrosymmetry,
            pr_4,
            flag_f_asym_a=flag_sft_ccs_asym,
            flag_scattering_length=flag_atom_para_form_factor,
        )
        if flag_dict:
            dict_in_out["sft_ccs"] = sft_ccs

    dder = {}
    if flag_unit_cell_parameters:
        dder["unit_cell_parameters"] = None

    if flag_atom_para_fract_xyz:
        dder["atom_para_fract_xyz"] = None

    if flag_atom_para_occupancy:
        dder["atom_para_occupancy"] = None

    if flag_atom_para_b_iso:
        dder["atom_para_b_iso"] = None

    if flag_atom_para_beta:
        dder["atom_para_beta"] = None
    if flag_atom_para_susceptibility:
        dder["atom_para_susceptibility"] = (
            dder_sft_ccs["f_asym_a_real"][:, na, :, :]
            * dder_sft_ccs_asym["atom_para_susceptibility"]
            + dder_sft_ccs["f_asym_a_imag"][:, na, :, :]
            * dder_sft_ccs_asym["atom_para_susceptibility"]
        )
    return sft_ccs, dder


def calc_f_m_perp_ordered_by_dictionary(
    dict_crystal, dict_in_out, flag_use_precalculated_data: bool = False
):
    dict_crystal_keys = dict_crystal.keys()
    dict_in_out_keys = dict_in_out.keys()
    necessary_crystal_keys = set(["unit_cell_parameters", "full_mcif_elems"])
    diff_set_crystal = necessary_crystal_keys.difference(
        set(dict_crystal_keys)
    )
    if len(diff_set_crystal) != 0:
        raise AttributeError(
            f"The following attributes have to be defined {diff_set_crystal:}"
        )

    necessary_in_out_keys = set(
        [
            "index_hkl",
        ]
    )
    diff_set_in_out = necessary_in_out_keys.difference(set(dict_in_out_keys))
    if len(diff_set_in_out) != 0:
        raise AttributeError(
            f"The following attributes have to be defined {diff_set_in_out:}"
        )

    index_hkl = dict_in_out["index_hkl"]

    non_zero_keys = set(
        [
            "mag_atom_lande_factor",
            "mag_atom_kappa",
            "mag_atom_j0_parameters",
            "mag_atom_j2_parameters",
        ]
    )
    diff_set_crystal = non_zero_keys.difference(set(dict_crystal_keys))

    if len(diff_set_crystal) != 0:
        sft_ccs = numpy.zeros((3, index_hkl.shape[-1]), dtype=complex)
        dder = {}
        return sft_ccs, dder

    if "flag_only_orbital" in dict_in_out_keys:
        flag_only_orbital = dict_in_out["flag_only_orbital"]
    else:
        flag_only_orbital = False

    full_mcif_elems = dict_crystal["full_mcif_elems"]

    unit_cell_parameters = dict_crystal["unit_cell_parameters"]
    atom_ordered_index = dict_crystal["atom_ordered_index"]
    atom_ordered_fract_xyz = dict_crystal["atom_fract_xyz"][
        :, atom_ordered_index
    ]
    atom_ordered_sc_fract = dict_crystal["atom_site_sc_fract"][
        :, atom_ordered_index
    ]
    atom_ordered_sc_b = dict_crystal["atom_site_sc_b"][:, atom_ordered_index]
    atom_ordered_fract_xyz = (
        calc_m_v(
            atom_ordered_sc_fract,
            numpy.mod(atom_ordered_fract_xyz, 1),
            flag_m=False,
            flag_v=False,
        )[0]
        + atom_ordered_sc_b
    )

    atom_ordered_occupancy = dict_crystal["atom_occupancy"][atom_ordered_index]
    atom_ordered_b_iso = dict_crystal["atom_b_iso"][atom_ordered_index]
    atom_beta = dict_crystal["atom_beta"]
    if "atom_site_aniso_sc_beta" in dict_crystal_keys:
        atom_site_aniso_sc_beta = dict_crystal["atom_site_aniso_sc_beta"]
        atom_site_aniso_index = dict_crystal["atom_site_aniso_index"]
        atom_sc_beta = numpy.zeros((6,) + atom_beta.shape, dtype=float)
        atom_sc_beta[:, :, atom_site_aniso_index] = atom_site_aniso_sc_beta
        atom_beta = (atom_sc_beta * numpy.expand_dims(atom_beta, axis=0)).sum(
            axis=1
        )
    atom_ordered_beta = atom_beta[:, atom_ordered_index]

    mag_atom_ordered_index = dict_crystal["mag_atom_ordered_index"]
    atom_ordered_lande_factor = dict_crystal["mag_atom_lande_factor"][
        mag_atom_ordered_index
    ]
    atom_ordered_kappa = dict_crystal["mag_atom_kappa"][mag_atom_ordered_index]
    atom_ordered_j0_parameters = dict_crystal["mag_atom_j0_parameters"][
        :, mag_atom_ordered_index
    ]
    atom_ordered_j2_parameters = dict_crystal["mag_atom_j2_parameters"][
        :, mag_atom_ordered_index
    ]

    atom_ordered_moment_crystalaxis_dn = dict_crystal[
        "atom_ordered_moment_crystalaxis_dn"
    ]
    flags_atom_ordered_moment_crystalaxis_dn = dict_crystal[
        "flags_atom_ordered_moment_crystalaxis_dn"
    ]
    flag_unit_cell_parameters = numpy.any(
        dict_crystal["flags_unit_cell_parameters"]
    )
    flag_atom_ordered_fract_xyz = numpy.any(
        dict_crystal["flags_atom_fract_xyz"][:, atom_ordered_index]
    )
    flag_atom_ordered_occupancy = numpy.any(
        dict_crystal["flags_atom_occupancy"][atom_ordered_index]
    )
    flag_atom_ordered_b_iso = numpy.any(
        dict_crystal["flags_atom_b_iso"][atom_ordered_index]
    )
    flag_atom_ordered_beta = numpy.any(
        dict_crystal["flags_atom_beta"][:, atom_ordered_index]
    )
    flag_atom_ordered_lande_factor = numpy.any(
        dict_crystal["flags_mag_atom_lande_factor"][mag_atom_ordered_index]
    )
    flag_atom_ordered_kappa = numpy.any(
        dict_crystal["flags_mag_atom_kappa"][mag_atom_ordered_index]
    )

    flag_multipole = all(
        [
            hh in dict_crystal.keys()
            for hh in ["mag_atom_rho_multipole_index", "atom_local_axes_label"]
        ]
    )
    if flag_multipole:
        mag_atom_rho_multipole_index = dict_crystal[
            "mag_atom_rho_multipole_index"
        ]
        atom_rho_multipole_transformation_matrix = dict_crystal[
            "atom_local_axes_matrix"
        ][:, mag_atom_rho_multipole_index]
        atom_rho_multipole_plm = dict_crystal["atom_rho_multipole_plm"][
            :, mag_atom_rho_multipole_index
        ]
        atom_rho_multipole_lm = dict_crystal["atom_rho_multipole_lm"]
        l_kappa = dict_crystal["atom_rho_multipole_kappa"][
            mag_atom_rho_multipole_index
        ]
        l_n_hh, l_zeta_hh, l_coeff_hh = dict_crystal[
            "atom_rho_multipole_n_zeta_coeff"
        ]
        # FIXME: not sure that it is needed
        l_n = [l_n_hh[hh] for hh in mag_atom_rho_multipole_index]
        l_zeta = [l_zeta_hh[hh] for hh in mag_atom_rho_multipole_index]
        l_coeff = [l_coeff_hh[hh] for hh in mag_atom_rho_multipole_index]

        magnetic_form_factor, dder = (
            calc_magnetic_form_factor_by_multipole_model(
                index_hkl,
                full_mcif_elems[4:-1, :],
                atom_rho_multipole_transformation_matrix,
                atom_rho_multipole_plm,
                atom_rho_multipole_lm,
                l_kappa,
                l_n,
                l_zeta,
                l_coeff,
                unit_cell_parameters,
            )
        )
    else:
        flag_sthovl = flag_unit_cell_parameters
        if (
            flag_use_precalculated_data
            and ("sthovl" in dict_in_out_keys)
            and not (flag_sthovl)
        ):
            sthovl = dict_in_out["sthovl"]
        else:
            sthovl, dder_sthovl = calc_sthovl_by_unit_cell_parameters(
                index_hkl,
                unit_cell_parameters,
                flag_unit_cell_parameters=flag_unit_cell_parameters,
            )
            dict_in_out["sthovl"] = sthovl

        atom_ordered_form_factor, dder_ff = calc_form_factor(
            sthovl[:, na],
            atom_ordered_lande_factor[na, :],
            atom_ordered_kappa[na, :],
            atom_ordered_j0_parameters[:, na, :],
            atom_ordered_j2_parameters[:, na, :],
            flag_lande_factor=flag_atom_ordered_lande_factor,
            flag_only_orbital=flag_only_orbital,
            flag_sthovl=flag_sthovl,
            flag_kappa=flag_atom_ordered_kappa,
        )
        magnetic_form_factor = numpy.expand_dims(
            atom_ordered_form_factor, axis=1
        )

    flag_atom_ordered_moment_crystalaxis_dn = numpy.any(
        flags_atom_ordered_moment_crystalaxis_dn
    )
    f_m_perp_o, dder = calc_f_m_perp_ordered(
        index_hkl,
        full_mcif_elems,
        unit_cell_parameters,
        atom_ordered_fract_xyz,
        atom_ordered_occupancy,
        atom_ordered_moment_crystalaxis_dn,
        atom_ordered_b_iso,
        atom_ordered_beta,
        atom_ordered_lande_factor,
        atom_ordered_kappa,
        atom_ordered_j0_parameters,
        atom_ordered_j2_parameters,
        magnetic_form_factor,
        dict_in_out=dict_in_out,
        flag_only_orbital=flag_only_orbital,
        flag_unit_cell_parameters=flag_unit_cell_parameters,
        flag_atom_ordered_fract_xyz=flag_atom_ordered_fract_xyz,
        flag_atom_ordered_occupancy=flag_atom_ordered_occupancy,
        flag_atom_ordered_moment_crystalaxis_dn=flag_atom_ordered_moment_crystalaxis_dn,
        flag_atom_ordered_b_iso=flag_atom_ordered_b_iso,
        flag_atom_ordered_beta=flag_atom_ordered_beta,
        flag_atom_ordered_lande_factor=flag_atom_ordered_lande_factor,
        flag_atom_ordered_kappa=flag_atom_ordered_kappa,
        flag_use_precalculated_data=flag_use_precalculated_data,
    )
    return f_m_perp_o, dder


def calc_f_m_perp_ordered(
    index_hkl,
    full_mcif_elems,
    unit_cell_parameters,
    atom_ordered_fract_xyz,
    atom_ordered_occupancy,
    atom_ordered_moment_crystalaxis_dn,
    atom_ordered_b_iso,
    atom_ordered_beta,
    atom_ordered_lande_factor,
    atom_ordered_kappa,
    atom_ordered_j0_parameters,
    atom_ordered_j2_parameters,
    magnetic_form_factor,
    dict_in_out: dict = None,
    flag_only_orbital: bool = False,
    flag_unit_cell_parameters: bool = False,
    flag_atom_ordered_fract_xyz: bool = False,
    flag_atom_ordered_occupancy: bool = False,
    flag_atom_ordered_moment_crystalaxis_dn: bool = False,
    flag_atom_ordered_b_iso: bool = False,
    flag_atom_ordered_beta: bool = False,
    flag_atom_ordered_lande_factor: bool = False,
    flag_atom_ordered_kappa: bool = False,
    flag_use_precalculated_data: bool = False,
):

    if dict_in_out is None:
        flag_dict = False
        dict_in_out_keys = []
    else:
        flag_dict = True
        dict_in_out_keys = dict_in_out.keys()
        if "index_hkl" in dict_in_out_keys:
            if numpy.any(dict_in_out["index_hkl"] != index_hkl):
                dict_in_out.clear()
                dict_in_out["index_hkl"] = index_hkl

    if flag_use_precalculated_data and (
        "atom_ordered_multiplicity" in dict_in_out_keys
    ):
        atom_ordered_multiplicity = dict_in_out["atom_ordered_multiplicity"]
    else:
        atom_symm_elems = get_atom_symm_elems_by_atom_fract_xyz(
            atom_ordered_fract_xyz
        )

        atom_ordered_multiplicity = calc_multiplicity_by_atom_symm_elems(
            full_mcif_elems[:13], atom_symm_elems
        )
        if flag_dict:
            dict_in_out["atom_ordered_multiplicity"] = (
                atom_ordered_multiplicity
            )

    if (
        flag_use_precalculated_data
        and ("pr_1_atom_ordered" in dict_in_out_keys)
        and not (flag_atom_ordered_fract_xyz)
    ):
        pr_1 = dict_in_out["pr_1_atom_ordered"]
    else:
        pr_1, dder_pr_1 = calc_pr1(
            index_hkl,
            full_mcif_elems[:13],
            atom_ordered_fract_xyz,
            flag_fract_xyz=flag_atom_ordered_fract_xyz,
        )
        if flag_dict:
            dict_in_out["pr_1_atom_ordered"] = pr_1

    if flag_use_precalculated_data and ("pr_2" in dict_in_out_keys):
        pr_2 = dict_in_out["pr_2"]
    else:
        pr_2 = calc_pr2(index_hkl, full_mcif_elems[:13])
        if flag_dict:
            dict_in_out["pr_2"] = pr_2

    flag_sthovl = flag_unit_cell_parameters
    if (
        flag_use_precalculated_data
        and ("sthovl" in dict_in_out_keys)
        and not (flag_sthovl)
    ):
        sthovl = dict_in_out["sthovl"]
    else:
        sthovl, dder_sthovl = calc_sthovl_by_unit_cell_parameters(
            index_hkl,
            unit_cell_parameters,
            flag_unit_cell_parameters=flag_unit_cell_parameters,
        )
        if flag_dict:
            dict_in_out["sthovl"] = sthovl

    flag_hh = True
    if "flag_only_orbital" in dict_in_out_keys:
        flag_hh = flag_only_orbital == dict_in_out["flag_only_orbital"]
    dict_in_out["flag_only_orbital"] = flag_only_orbital

    flag_atom_ordered_form_factor = (
        flag_atom_ordered_lande_factor
        or flag_sthovl
        or flag_atom_ordered_kappa
    )
    if (
        flag_use_precalculated_data
        and ("atom_ordered_form_factor" in dict_in_out_keys)
        and not (flag_atom_ordered_form_factor)
        and flag_hh
    ):
        atom_ordered_form_factor = dict_in_out["atom_ordered_form_factor"]
    else:
        atom_ordered_form_factor, dder_ff = calc_form_factor(
            sthovl[:, na],
            atom_ordered_lande_factor[na, :],
            atom_ordered_kappa[na, :],
            atom_ordered_j0_parameters[:, na, :],
            atom_ordered_j2_parameters[:, na, :],
            flag_lande_factor=flag_atom_ordered_lande_factor,
            flag_only_orbital=flag_only_orbital,
            flag_sthovl=flag_sthovl,
            flag_kappa=flag_atom_ordered_kappa,
        )
        if flag_dict:
            dict_in_out["atom_ordered_form_factor"] = atom_ordered_form_factor

    # dimensions ["hkl", "reduced symmetry", "atom"]
    flag_debye_waller_factor = (
        flag_sthovl or flag_atom_ordered_b_iso or flag_atom_ordered_beta
    )
    if (
        flag_use_precalculated_data
        and ("atom_ordered_debye_waller_factor" in dict_in_out_keys)
        and not (flag_debye_waller_factor)
    ):
        debye_waller_factor = dict_in_out["atom_ordered_debye_waller_factor"]
    else:
        debye_waller_factor, dder_dw = calc_dwf(
            index_hkl[:, :, na, na],
            sthovl[:, na, na],
            atom_ordered_b_iso[na, na, :],
            atom_ordered_beta[:, na, na, :],
            full_mcif_elems[4:13, na, :, na],
            flag_sthovl=flag_sthovl,
            flag_b_iso=flag_atom_ordered_b_iso,
            flag_beta=flag_atom_ordered_beta,
        )
        if flag_dict:
            dict_in_out["atom_ordered_debye_waller_factor"] = (
                debye_waller_factor
            )


    m_norm, der_m_norm = calc_m_m_norm_by_unit_cell_parameters(
        unit_cell_parameters,
        flag_unit_cell_parameters=flag_unit_cell_parameters,
    )
    atom_ordered_moment_crystalaxis_xyz = calc_m_v(
        m_norm, atom_ordered_moment_crystalaxis_dn, flag_m=False, flag_v=False
    )[0]
    flag_atom_ordered_moment_crystalaxis_xyz = (
        flag_atom_ordered_moment_crystalaxis_dn
    )
    r_direct = full_mcif_elems[4:13, :]
    rm, der_rm = calc_m_v(
        r_direct[:, :, na],
        atom_ordered_moment_crystalaxis_xyz[:, na, :],
        flag_m=False,
        flag_v=flag_atom_ordered_moment_crystalaxis_xyz,
    )
    rm_ccs, der_rm_ccs = calc_m_v(
        m_norm,
        rm,
        flag_m=flag_unit_cell_parameters,
        flag_v=flag_atom_ordered_moment_crystalaxis_xyz,
    )

    m_m = calc_m_m_by_unit_cell_parameters(
        unit_cell_parameters,
        flag_unit_cell_parameters=flag_unit_cell_parameters,
    )[0]
    r_ccs = calc_m1_m2_inv_m1(m_m, r_direct)[0]
    det_r, der_det_r = calc_det_m(r_ccs)
    theta_s = full_mcif_elems[13, :]
    # det_r is taken into account when multiplication on r_direct is performed
    # moment_ccs = 0.2695*rm_ccs*(theta_s)[na, :, na]
    moment_ccs = 0.2695 * rm_ccs * (theta_s * det_r)[na, :, na]

    # [N_atom]
    hh_1 = atom_ordered_multiplicity * atom_ordered_occupancy
    # [N_hkl, N_atom]
    # hh_2 = atom_ordered_form_factor*hh_1[na, :]
    hh_2 = hh_1[na, :]
    # [N_hkl, N_s, N_atom]
    hh_3 = pr_1 * debye_waller_factor * magnetic_form_factor * hh_2[:, na, :]
    # [3, N_hkl]
    f_m = (
        pr_2[na, :, :]
        * (hh_3[na, :, :, :] * moment_ccs[:, na, :, :]).sum(axis=3)
    ).sum(axis=2) / pr_2.shape[-1]
    eq_ccs, dder_eq_ccs = calc_eq_ccs_by_unit_cell_parameters(
        index_hkl,
        unit_cell_parameters,
        flag_unit_cell_parameters=flag_unit_cell_parameters,
    )
    dict_in_out["eq_ccs"] = eq_ccs

    dict_in_out["f_m_o"] = f_m

    flag_f_m = False
    f_m_perp_o, dder_f_m_perp_o = calc_vector_product_v1_v2_v1(
        eq_ccs, f_m, flag_v1=flag_unit_cell_parameters, flag_v2=flag_f_m
    )
    dict_in_out["f_m_perp_o"] = f_m_perp_o
    dder = {}
    if flag_atom_ordered_moment_crystalaxis_xyz:
        dder["atom_ordered_moment_crystalaxis_dn"] = None

    return f_m_perp_o, dder


def calc_bulk_susceptibility_by_dictionary(
    dict_crystal, dict_in_out, flag_use_precalculated_data: bool = False
):
    dict_crystal_keys = dict_crystal.keys()

    necessary_crystal_keys = set(
        [
            "unit_cell_parameters",
        ]
    )
    diff_set_crystal = necessary_crystal_keys.difference(
        set(dict_crystal_keys)
    )
    if len(diff_set_crystal) != 0:
        raise AttributeError(
            f"The following attributes have to be defined {diff_set_crystal:}"
        )

    flag_reduced_symm_elems = (
        len(
            set(
                ["reduced_symm_elems", "centrosymmetry", "translation_elems"]
            ).difference(set(dict_crystal_keys))
        )
        == 0
    )
    flag_full_symm_elems = (
        len(
            set(
                [
                    "full_symm_elems",
                ]
            ).difference(set(dict_crystal_keys))
        )
        == 0
    )
    flag_full_mcif_elems = (
        len(
            set(
                [
                    "full_mcif_elems",
                ]
            ).difference(set(dict_crystal_keys))
        )
        == 0
    )

    if not (
        flag_reduced_symm_elems or flag_full_symm_elems or flag_full_mcif_elems
    ):
        raise AttributeError("The symmetry elements have to be defined.")

    if flag_reduced_symm_elems:
        reduced_symm_elems = dict_crystal["reduced_symm_elems"]
        centrosymmetry = dict_crystal["centrosymmetry"]
        if centrosymmetry:
            centrosymmetry_position = dict_crystal["centrosymmetry_position"]
        else:
            centrosymmetry_position = None
        translation_elems = dict_crystal["translation_elems"]
    elif flag_full_symm_elems:
        full_symm_elems = dict_crystal["full_symm_elems"]
        reduced_symm_elems = full_symm_elems
        centrosymmetry = False
        centrosymmetry_position = None
        translation_elems = numpy.array([[0], [0], [0], [1]], dtype=int)
    elif flag_full_mcif_elems:
        full_mcif_elems = dict_crystal["full_mcif_elems"]
        reduced_symm_elems = full_mcif_elems[:13]
        centrosymmetry = False
        centrosymmetry_position = None
        translation_elems = numpy.array([[0], [0], [0], [1]], dtype=int)

    unit_cell_parameters = dict_crystal["unit_cell_parameters"]
    atom_para_index = dict_crystal["atom_para_index"]
    atom_para_fract_xyz = dict_crystal["atom_fract_xyz"][:, atom_para_index]
    atom_para_sc_fract = dict_crystal["atom_site_sc_fract"][:, atom_para_index]
    atom_para_sc_b = dict_crystal["atom_site_sc_b"][:, atom_para_index]
    atom_para_fract_xyz = (
        calc_m_v(
            atom_para_sc_fract,
            numpy.mod(atom_para_fract_xyz, 1),
            flag_m=False,
            flag_v=False,
        )[0]
        + atom_para_sc_b
    )
    atom_para_occupancy = dict_crystal["atom_occupancy"][atom_para_index]

    atom_para_susceptibility = dict_crystal["atom_para_susceptibility"]
    atom_para_sc_chi = dict_crystal[
        "atom_para_sc_chi"
    ]  # FIXME not sure that it should be in crystal dictionary but not dict_in_out
    flag_unit_cell_parameters = numpy.any(
        dict_crystal["flags_unit_cell_parameters"]
    )
    flag_atom_para_fract_xyz = numpy.any(
        dict_crystal["flags_atom_fract_xyz"][:, atom_para_index]
    )
    flag_atom_para_occupancy = numpy.any(
        dict_crystal["flags_atom_occupancy"][atom_para_index]
    )
    flag_atom_para_susceptibility = numpy.any(
        dict_crystal["flags_atom_para_susceptibility"]
    )

    bulk_susceptibility, dder = calc_bulk_susceptibility(
        reduced_symm_elems,
        centrosymmetry,
        centrosymmetry_position,
        translation_elems,
        unit_cell_parameters,
        atom_para_fract_xyz,
        atom_para_occupancy,
        atom_para_susceptibility,
        atom_para_sc_chi,
        dict_in_out=dict_in_out,
        flag_unit_cell_parameters=flag_unit_cell_parameters,
        flag_atom_para_fract_xyz=flag_atom_para_fract_xyz,
        flag_atom_para_occupancy=flag_atom_para_occupancy,
        flag_atom_para_susceptibility=flag_atom_para_susceptibility,
        flag_use_precalculated_data=flag_use_precalculated_data,
    )
    error_bars = numpy.zeros_like(bulk_susceptibility)
    if "atom_para_susceptibility" in dder.keys():
        dder_aps = dder["atom_para_susceptibility"]
        sigma_aps = dict_crystal["atom_para_susceptibility_sigma"]
        error_bars = numpy.sum(
            dder_aps * numpy.expand_dims(sigma_aps, axis=0), axis=1
        )

    return bulk_susceptibility, error_bars


def calc_bulk_susceptibility(
    reduced_symm_elems,
    centrosymmetry,
    centrosymmetry_position,
    translation_elems,
    unit_cell_parameters,
    atom_para_fract_xyz,
    atom_para_occupancy,
    atom_para_susceptibility,
    atom_para_sc_chi,
    dict_in_out: dict = None,
    flag_unit_cell_parameters: bool = False,
    flag_atom_para_fract_xyz: bool = False,
    flag_atom_para_occupancy: bool = False,
    flag_atom_para_susceptibility: bool = False,
    flag_use_precalculated_data: bool = False,
):

    if dict_in_out is None:
        flag_dict = True
        dict_in_out_keys = []
    else:
        dict_in_out_keys = dict_in_out.keys()
        flag_dict = True

    if flag_use_precalculated_data and (
        "atom_para_multiplicity" in dict_in_out_keys
    ):
        atom_para_multiplicity = dict_in_out["atom_para_multiplicity"]
    else:
        atom_symm_elems = get_atom_symm_elems_by_atom_fract_xyz(
            atom_para_fract_xyz
        )

        if "full_symm_elems" in dict_in_out_keys:
            full_symm_elems = dict_in_out["full_symm_elems"]
        else:
            full_symm_elems = calc_full_symm_elems_by_reduced(
                reduced_symm_elems,
                centrosymmetry,
                centrosymmetry_position,
                translation_elems,
            )
            if flag_dict:
                dict_in_out["full_symm_elems"] = full_symm_elems
        atom_para_multiplicity = calc_multiplicity_by_atom_symm_elems(
            full_symm_elems, atom_symm_elems
        )
        if flag_dict:
            dict_in_out["atom_para_multiplicity"] = atom_para_multiplicity

    flag_pr_5 = flag_unit_cell_parameters
    if (
        flag_use_precalculated_data
        and ("pr_5" in dict_in_out_keys)
        and not (flag_pr_5)
    ):
        pr_5 = dict_in_out["pr_5"]
    else:
        pr_5, dder_pr_5 = calc_pr5(
            reduced_symm_elems,
            unit_cell_parameters,
            flag_unit_cell_parameters=flag_unit_cell_parameters,
        )
        if flag_dict:
            dict_in_out["pr_5"] = pr_5

    mas_constr = (atom_para_sc_chi * atom_para_susceptibility[na, :, :]).sum(
        axis=1
    )

    hh, dder_hh = calc_m_q_inv_m(
        pr_5[:, :, na],
        mas_constr[:, na, :],
        flag_m=False,
        flag_q=flag_atom_para_susceptibility,
    )

    chi_bulk = (
        ((atom_para_multiplicity * atom_para_occupancy)[na, na, :] * hh).sum(
            axis=2
        )
    ).sum(axis=1) / pr_5.shape[-1]
    dder = {}
    return chi_bulk, dder
