"""Write the reconstructed MEM density into .den files (split from mempy_by_dictionary.py)."""

import logging
import numpy

from cryspy.A_functions_base.mempy import (
    save_spin_density_into_file,
    transform_to_spin_density_for_mcif,
    transform_to_density_auc_to_moments,
)







logger = logging.getLogger(__name__)


def save_density_to_den_file(
    dict_in_out, dict_mem_parameters, dict_crystal, l_dict_diffrn
):
    dict_in_out_keys = dict_in_out.keys()
    channel_col = "density_channel_col" in dict_in_out_keys
    channel_ani = "density_channel_ani" in dict_in_out_keys
    channel_nucl = "density_channel_nucl" in dict_in_out_keys
    file_spin_density = dict_mem_parameters["file_spin_density"]
    file_magnetization_density = dict_mem_parameters[
        "file_magnetization_density"
    ]
    file_nuclear_density = dict_mem_parameters["file_nuclear_density"]
    magnetization_plus = dict_in_out.get("magnetization_plus", None)
    magnetization_minus = dict_in_out.get("magnetization_minus", None)

    dict_crystal_keys = dict_crystal.keys()
    if "atom_para_label" in dict_crystal_keys:
        mag_atom_label = dict_crystal["atom_para_label"]
    elif "atom_ordered_label" in dict_crystal_keys:
        mag_atom_label = dict_crystal["atom_ordered_label"]
    else:
        mag_atom_label = []

    if "full_symm_elems" in dict_crystal.keys():
        flag_mcif = False

        reduced_symm_elems = dict_crystal["reduced_symm_elems"]
        centrosymmetry = dict_crystal["centrosymmetry"]
        if centrosymmetry:
            centrosymmetry_position = dict_crystal["centrosymmetry_position"]
        else:
            centrosymmetry_position = None
        translation_elems = dict_crystal["translation_elems"]
    elif "full_mcif_elems" in dict_crystal.keys():
        full_mcif_elems = dict_crystal["full_mcif_elems"]
        flag_mcif = True

    index_auc = dict_in_out["index_auc"]
    point_multiplicity = dict_in_out["point_multiplicity"]
    n_abc = dict_mem_parameters["points_abc"]
    unit_cell_parameters = dict_crystal["unit_cell_parameters"]

    point_atom_symm_elems_auc = dict_in_out.get(
        "point_atom_symm_elems_auc", None
    )
    atom_ordered_m_dn = dict_crystal.get(
        "atom_ordered_moment_crystalaxis_dn", None
    )
    if channel_ani and (atom_ordered_m_dn is not None):
        flag_moment_ordered = True
    else:
        flag_moment_ordered = False

    if channel_col and (file_spin_density is not None) and flag_mcif:
        density_col_best = dict_in_out.get("density_channel_col", None)
        spin_density_auc = density_col_best * numpy.array(
            [
                [
                    magnetization_plus,
                ],
                [
                    magnetization_minus,
                ],
            ],
            dtype=float,
        )
        for dict_diffrn in l_dict_diffrn:
            eh_ccs = dict_diffrn["matrix_u"][6:]
            index_full, spin_density = transform_to_spin_density_for_mcif(
                index_auc,
                point_multiplicity,
                spin_density_auc,
                n_abc,
                full_mcif_elems,
                eh_ccs,
            )
            centrosymmetry = False
            centrosymmetry_position = None
            translation_elems = numpy.array(
                [
                    [
                        0,
                    ],
                    [
                        0,
                    ],
                    [
                        0,
                    ],
                    [
                        1,
                    ],
                ],
                dtype=int,
            )
            one_symm_elems = full_mcif_elems[:13, :1]
            if len(l_dict_diffrn) <= 1:
                file_spin_density_diffrn = file_spin_density
            else:
                file_spin_density_diffrn = (
                    ".".join(file_spin_density.split(".")[:-1])
                    + dict_diffrn["type_name"]
                    + ".den"
                )

            save_spin_density_into_file(
                file_spin_density_diffrn,
                index_full,
                spin_density,
                n_abc,
                unit_cell_parameters,
                one_symm_elems,
                translation_elems,
                centrosymmetry,
                centrosymmetry_position,
            )
            logger.info(f"\nReconstructed spin density is written in file '{file_spin_density_diffrn:}'.")
    elif channel_col and (file_spin_density is not None):
        density_col_best = dict_in_out.get("density_channel_col", None)
        spin_density = density_col_best * numpy.array(
            [
                [
                    magnetization_plus,
                ],
                [
                    magnetization_minus,
                ],
            ],
            dtype=float,
        )
        save_spin_density_into_file(
            file_spin_density,
            index_auc,
            spin_density,
            n_abc,
            unit_cell_parameters,
            reduced_symm_elems,
            translation_elems,
            centrosymmetry,
            centrosymmetry_position,
        )
        logger.info(f"\nReconstructed spin density is written in file '{file_spin_density:}'.")
    if channel_nucl and (file_nuclear_density is not None):
        density_nucl_best = dict_in_out.get("density_channel_nucl", None)
        atom_scat_length_neutron = dict_crystal["atom_scat_length_neutron"]
        as_real = numpy.expand_dims(atom_scat_length_neutron.real, axis=0)
        as_pos_real = numpy.where(as_real > 0, as_real, 0)
        as_neg_real = numpy.where(as_real < 0, as_real, 0)

        as_imag = numpy.expand_dims(atom_scat_length_neutron.imag, axis=0)
        as_pos_imag = numpy.where(as_imag > 0, as_imag, 0)
        as_neg_imag = numpy.where(as_imag < 0, as_imag, 0)

        scattering_density_real = numpy.stack(
            [
                (density_nucl_best * as_pos_real).sum(axis=1),
                (density_nucl_best * as_neg_real).sum(axis=1),
            ],
            axis=0,
        )

        scattering_density_imag = numpy.stack(
            [
                (density_nucl_best * as_pos_imag).sum(axis=1),
                (density_nucl_best * as_neg_imag).sum(axis=1),
            ],
            axis=0,
        )

        if flag_mcif:
            reduced_symm_elems = full_mcif_elems[:13, :]
            translation_elems = numpy.array(
                [
                    [
                        0,
                    ],
                    [
                        0,
                    ],
                    [
                        0,
                    ],
                    [
                        1,
                    ],
                ],
                dtype=int,
            )
            centrosymmetry = False
            centrosymmetry_position = None
        file_real = file_nuclear_density.split(".den")[0] + "_real.den"
        file_imag = file_nuclear_density.split(".den")[0] + "_imag.den"
        save_spin_density_into_file(
            file_real,
            index_auc,
            scattering_density_real,
            n_abc,
            unit_cell_parameters,
            reduced_symm_elems,
            translation_elems,
            centrosymmetry,
            centrosymmetry_position,
        )
        save_spin_density_into_file(
            file_imag,
            index_auc,
            scattering_density_imag,
            n_abc,
            unit_cell_parameters,
            reduced_symm_elems,
            translation_elems,
            centrosymmetry,
            centrosymmetry_position,
        )
        logger.info(f"\nReconstructed nuclear density is written in files '{file_real:}', '{file_imag:}'.")

    if channel_ani and (file_magnetization_density is not None):
        density_ani_best = dict_in_out.get("density_channel_ani", None)

        # spin_density = numpy.stack([density_ani_best, numpy.zeros_like(density_ani_best)], axis=0)
        if flag_moment_ordered:
            index_full, moment_3d = transform_to_density_auc_to_moments(
                index_auc,
                point_atom_symm_elems_auc,
                point_multiplicity,
                density_ani_best,
                n_abc,
                full_mcif_elems,
                atom_ordered_m_dn,
                unit_cell_parameters,
            )

            centrosymmetry = False
            centrosymmetry_position = None
            translation_elems = numpy.array(
                [
                    [
                        0,
                    ],
                    [
                        0,
                    ],
                    [
                        0,
                    ],
                    [
                        1,
                    ],
                ],
                dtype=int,
            )
            one_symm_elems = full_mcif_elems[:13, :1]
            for ind, lab in enumerate(["mx", "my", "mz"]):
                moment_component = moment_3d[ind]
                spin_density_component = numpy.zeros(
                    (2,) + moment_component.shape, dtype=float
                )
                flag_positive = moment_component >= 0
                spin_density_component[0, flag_positive] = moment_component[
                    flag_positive
                ]
                spin_density_component[1, ~flag_positive] = moment_component[
                    ~flag_positive
                ]
                file_magnetization_density_component = (
                    ".".join(file_magnetization_density.split(".")[:-1])
                    + f"_{lab:}.den"
                )

                save_spin_density_into_file(
                    file_magnetization_density_component,
                    index_full,
                    spin_density_component,
                    n_abc,
                    unit_cell_parameters,
                    one_symm_elems,
                    translation_elems,
                    centrosymmetry,
                    centrosymmetry_position,
                )
                logger.info(f"\nReconstructed spin density is written in file '{file_magnetization_density_component:}'.")

        if flag_mcif:
            reduced_symm_elems = full_mcif_elems[:13, :]
            translation_elems = numpy.array(
                [
                    [
                        0,
                    ],
                    [
                        0,
                    ],
                    [
                        0,
                    ],
                    [
                        1,
                    ],
                ],
                dtype=int,
            )
            centrosymmetry = False
            centrosymmetry_position = None
        # save_spin_density_into_file(file_magnetization_density, index_auc, spin_density, n_abc, unit_cell_parameters,
        #       reduced_symm_elems, translation_elems, centrosymmetry, centrosymmetry_position)
        # print(f"\nReconstructed magnetization density is written in file '{file_magnetization_density:}'.")
        # np_atom_label_unique = numpy.unique(atom_label_auc_ani)
        for ind_atom, label in enumerate(mag_atom_label):
            spin_density_atom = numpy.stack(
                [
                    density_ani_best[:, ind_atom],
                    numpy.zeros_like(density_ani_best[:, ind_atom]),
                ],
                axis=0,
            )
            file_magnetization_density_atom = (
                ".".join(file_magnetization_density.split(".")[:-1])
                + "_"
                + label
                + ".den"
            )
            save_spin_density_into_file(
                file_magnetization_density_atom,
                index_auc,
                spin_density_atom,
                n_abc,
                unit_cell_parameters,
                reduced_symm_elems,
                translation_elems,
                centrosymmetry,
                centrosymmetry_position,
            )
            logger.info(f"Reconstructed magnetization density around atom '{label:}' is written in file '{file_magnetization_density_atom:}'.")
    return
