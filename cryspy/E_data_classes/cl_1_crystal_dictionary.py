"""Dictionary form of Crystal and the .den files of its multipole densities.

The methods of Crystal call these functions (split from cl_1_crystal.py).
"""

import logging
import numpy


from cryspy.A_functions_base.database import DATABASE
from cryspy.A_functions_base.charge_form_factor import (
    calc_jl_for_ion,
    calc_scattering_amplitude_tabulated,
    get_atom_name_ion_charge_shell,
    get_atomic_symbol_ion_charge_isotope_number_by_ion_symbol,
)
from cryspy.A_functions_base.debye_waller_factor import (
    calc_param_iso_aniso_by_b_iso_beta,
    calc_u_ij_by_beta,
)
from cryspy.A_functions_base.matrix_operations import (
    calc_m_v,
)
from cryspy.A_functions_base.magnetic_form_factor import get_j0_j2_parameters
from cryspy.A_functions_base.multipoles import (
    calc_transformation_matrix_from_global_to_local,
    calc_multipole_density,
)
from cryspy.A_functions_base.unit_cell import (
    calc_reciprocal_by_unit_cell_parameters,
    calc_m_m_by_unit_cell_parameters,
)

from cryspy.A_functions_base.function_3_mcif import calc_multiplicity
from cryspy.A_functions_base.mempy import (
    form_basins_2,
    calc_symm_elem_points_by_index_points,
    save_spin_density_into_file,
    transform_to_density_auc_to_moments,
)
from cryspy.A_functions_base.symmetry_elements import (
    calc_asymmetric_unit_cell_indexes,
)


from cryspy.A_functions_base.symmetry_elements import (
    calc_full_mag_elems,
    calc_symm_flags,
    define_bravais_type_by_symm_elems,
)
from cryspy.A_functions_base.symmetry_constraints import (
    calc_sc_beta,
    calc_sc_fract_sc_b,
    calc_sc_chi,
    calc_sc_ordered,
)

from cryspy.B_parent_classes.preocedures import take_items_by_class

from cryspy.C_item_loop_classes.cl_1_atom_local_axes import AtomLocalAxesL
from cryspy.C_item_loop_classes.cl_1_atom_rho_multipole import (
    AtomRhoMultipoleL,
)
from cryspy.C_item_loop_classes.cl_1_atom_site import AtomSiteL
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
    AtomSiteMomentL,
)
from cryspy.C_item_loop_classes.cl_1_atom_electron_configuration import (
    AtomElectronConfigurationL,
)
from cryspy.C_item_loop_classes.cl_1_cell import Cell
from cryspy.C_item_loop_classes.cl_1_space_group_symop_magn_centering import (
    SpaceGroupSymopMagnCenteringL,
)

from cryspy.C_item_loop_classes.cl_2_space_group_symop_magn_operation import (
    SpaceGroupSymopMagnOperationL,
)
from cryspy.C_item_loop_classes.cl_2_space_group import SpaceGroup
from cryspy.exceptions import RcifFormatError



logger = logging.getLogger(__name__)

na = numpy.newaxis


def form_crystal_dictionary(crystal, ddict):
    """Body of Crystal.get_dictionary, see its docstring."""
    space_group, cell, atom_site = None, None, None
    atom_site_susceptibility, atom_electron_configuration = None, None
    atom_site_exchange = None
    atom_type_scat = None
    atom_type = None
    atom_site_aniso = None
    atom_site_scat = None
    atom_site_moment = None
    space_group_symop_magn_centering = None
    space_group_symop_magn_operation = None
    atom_rho_multipole = None
    atom_local_axes = None

    l_obj = take_items_by_class(crystal, (SpaceGroup,))
    if len(l_obj) > 0:
        space_group = l_obj[0]

    l_obj = take_items_by_class(crystal, (Cell,))
    if len(l_obj) > 0:
        cell = l_obj[0]

    l_obj = take_items_by_class(crystal, (AtomSiteL,))
    if len(l_obj) > 0:
        atom_site = l_obj[0]

    l_obj = take_items_by_class(crystal, (AtomSiteAnisoL,))
    if len(l_obj) > 0:
        atom_site_aniso = l_obj[0]

    l_obj = take_items_by_class(crystal, (AtomSiteScatL,))
    if len(l_obj) > 0:
        atom_site_scat = l_obj[0]

    l_obj = take_items_by_class(crystal, (AtomTypeL,))
    if len(l_obj) > 0:
        atom_type = l_obj[0]

    l_obj = take_items_by_class(crystal, (AtomTypeScatL,))
    if len(l_obj) > 0:
        atom_type_scat = l_obj[0]

    l_obj = take_items_by_class(crystal, (AtomSiteSusceptibilityL,))
    if len(l_obj) > 0:
        atom_site_susceptibility = l_obj[0]

    l_obj = take_items_by_class(crystal, (AtomSiteExchangeL,))
    if len(l_obj) > 0:
        atom_site_exchange = l_obj[0]

    l_obj = take_items_by_class(crystal, (AtomElectronConfigurationL,))
    if len(l_obj) > 0:
        atom_electron_configuration = l_obj[0]

    l_obj = take_items_by_class(crystal, (AtomSiteMomentL,))
    if len(l_obj) > 0:
        atom_site_moment = l_obj[0]

    l_obj = take_items_by_class(crystal, (SpaceGroupSymopMagnCenteringL,))
    if len(l_obj) > 0:
        space_group_symop_magn_centering = l_obj[0]

    l_obj = take_items_by_class(crystal, (SpaceGroupSymopMagnOperationL,))
    if len(l_obj) > 0:
        space_group_symop_magn_operation = l_obj[0]

    l_obj = take_items_by_class(crystal, (AtomRhoMultipoleL,))
    if len(l_obj) > 0:
        atom_rho_multipole = l_obj[0]

    l_obj = take_items_by_class(crystal, (AtomLocalAxesL,))
    if len(l_obj) > 0:
        atom_local_axes = l_obj[0]

    ddict["name"] = crystal.data_name
    ddict["type_name"] = crystal.get_name()
    if space_group is not None:
        full_symm_elems = (
            space_group.full_space_group_symop.get_symm_elems()
        )
        ddict["full_symm_elems"] = full_symm_elems
        reduced_symm_elems = (
            space_group.reduced_space_group_symop.get_symm_elems()
        )
        ddict["reduced_symm_elems"] = reduced_symm_elems
        ddict["centrosymmetry"] = space_group.centrosymmetry
        if ddict["centrosymmetry"]:
            p_centr = space_group.pcentr
            lcm = numpy.lcm.reduce([fr.denominator for fr in p_centr])
            vals = [
                int(fr.numerator * lcm / fr.denominator) for fr in p_centr
            ]
            vals.append(lcm)
            ddict["centrosymmetry_position"] = numpy.array(vals, dtype=int)

        shift = space_group.shift
        l_vals = []
        for ind in range(len(shift)):
            lcm = numpy.lcm.reduce([fr.denominator for fr in shift[ind]])
            vals = [
                int(fr.numerator * lcm / fr.denominator)
                for fr in shift[ind]
            ]
            vals.append(lcm)
            l_vals.append(vals)
        ddict["translation_elems"] = numpy.array(
            l_vals, dtype=int
        ).transpose()

    if space_group_symop_magn_operation is not None:
        if space_group_symop_magn_centering is not None:
            se_sgsmo = (
                crystal.space_group_symop_magn_operation.get_sym_elems()
            )
            se_sgsmc = (
                crystal.space_group_symop_magn_centering.get_sym_elems()
            )
            full_mcif_elems = calc_full_mag_elems(se_sgsmo, se_sgsmc)
            ddict["full_mcif_elems"] = full_mcif_elems

    if cell is not None:
        if cell.is_attribute("type_cell"):
            type_cell, it_coordinate_system_code = (
                cell.type_cell,
                cell.it_coordinate_system_code,
            )
        elif "full_mcif_elems" in ddict.keys():
            type_cell, it_coordinate_system_code = (
                define_bravais_type_by_symm_elems(full_mcif_elems)
            )
            cell.type_cell = type_cell
            cell.it_coordinate_system_code = it_coordinate_system_code
            cell.apply_constraints()

        unit_cell_parameters = ddict["unit_cell_parameters"]

        sc_uc, v_uc = calc_sc_v_unit_cell_parameters(
            type_cell, it_coordinate_system_code
        )
        ddict["sc_uc"] = sc_uc
        ddict["v_uc"] = v_uc

        ddict["unit_cell_parameters"] = (
            numpy.dot(sc_uc, unit_cell_parameters) + v_uc
        )
        unit_cell_parameters = ddict["unit_cell_parameters"]
        # ddict["flags_unit_cell_parameters"] = cell.get_flags_unit_cell_parameters() : it should be already defined from cell.get_dictionary()

    if atom_site is not None:
        atom_label = numpy.array(atom_site.label, dtype=str)
        atom_occupancy = numpy.array(atom_site.occupancy, dtype=float)
        atom_fract_xyz = numpy.array(
            [atom_site.fract_x, atom_site.fract_y, atom_site.fract_z],
            dtype=float,
        )
        atom_fract_xyz = numpy.mod(atom_fract_xyz, 1.0)
        if "full_mcif_elems" in ddict.keys():
            atom_multiplicity = calc_multiplicity(
                ddict["full_mcif_elems"], atom_fract_xyz
            )
            ddict["atom_multiplicity"] = atom_multiplicity
        elif "full_symm_elems" in ddict.keys():
            atom_multiplicity = calc_multiplicity(
                ddict["full_symm_elems"], atom_fract_xyz
            )
            ddict["atom_multiplicity"] = atom_multiplicity
        elif atom_site.is_attribute("multiplicity"):
            atom_multiplicity = numpy.array(
                atom_site.multiplicity, dtype=int
            )
            ddict["atom_multiplicity"] = atom_multiplicity
        atom_type_symbol = numpy.array(atom_site.type_symbol, dtype=str)
        table_sthovl = numpy.linspace(0, 2, 501)  # fro0 to 2 Å**-1
        d_dispersion = DATABASE["Dispersion"]
        table_wavelength = d_dispersion["table_wavelength"]
        l_table_atom_scattering_amplitude = []
        l_table_atom_dispersion = []
        flag_atom_scattering_amplitude = True
        for type_symbol in atom_type_symbol:
            try:
                jl = calc_jl_for_ion(table_sthovl, type_symbol)
                scattering_amplitude = jl[:, 0]
            except RcifFormatError:
                try:
                    scattering_amplitude = (
                        calc_scattering_amplitude_tabulated(
                            type_symbol, table_sthovl
                        )[0]
                    )
                except KeyError:
                    flag_atom_scattering_amplitude = False
            if flag_atom_scattering_amplitude:
                l_table_atom_scattering_amplitude.append(
                    scattering_amplitude
                )

            type_atom = get_atom_name_ion_charge_shell(type_symbol)[0]
            try:
                l_table_atom_dispersion.append(
                    d_dispersion[type_atom.capitalize()]
                )
            except KeyError:
                l_table_atom_dispersion.append(
                    numpy.zeros_like(table_wavelength)
                )
        if flag_atom_scattering_amplitude:
            ddict["table_sthovl"] = table_sthovl
            ddict["table_atom_scattering_amplitude"] = numpy.stack(
                l_table_atom_scattering_amplitude, axis=0
            )
        ddict["table_wavelength"] = table_wavelength
        ddict["table_atom_dispersion"] = numpy.stack(
            l_table_atom_dispersion, axis=0
        )

        atom_b_scat = numpy.array(
            atom_site.scat_length_neutron, dtype=complex
        )

        b_d = (
            3 * 10**6 * numpy.ones(atom_fract_xyz.shape[1:], dtype=int)
        )  # 10**6 it is precision
        b_elem = (atom_fract_xyz * b_d[na, :]).astype(int)
        gcd = numpy.gcd.reduce(numpy.gcd(b_elem, b_d[na, :]))

        atom_symm_elems = numpy.concatenate(
            [
                numpy.floor_divide(b_elem, gcd[na, :]),
                numpy.floor_divide(b_d, gcd)[na, :],
                numpy.zeros((9,) + atom_fract_xyz.shape[1:], dtype=int),
            ],
            axis=0,
        )

        if "full_mcif_elems" in ddict.keys():
            elems_fs = ddict["full_mcif_elems"]
        elif "full_symm_elems" in ddict.keys():
            elems_fs = ddict["full_symm_elems"]
        else:
            raise AttributeError("Symmetry elements are not found.")

        flag_2d = calc_symm_flags(
            elems_fs[:, :, na], atom_symm_elems[:, na, :]
        )

        l_sc_fract, l_sc_b = [], []
        for ind in range(atom_label.size):
            symm_elems = elems_fs[:, flag_2d[:, ind]]
            sc_fract, sc_b = calc_sc_fract_sc_b(
                symm_elems, atom_fract_xyz[:, ind]
            )
            l_sc_fract.append(sc_fract)
            l_sc_b.append(sc_b)
        atom_site_sc_fract = numpy.stack(l_sc_fract, axis=-1)

        ddict["atom_site_sc_fract"] = atom_site_sc_fract
        ddict["atom_site_sc_b"] = numpy.stack(l_sc_b, axis=-1)

        flags_atom_fract_xyz = atom_site.get_flags_atom_fract_xyz()
        flags_atom_fract_xyz_constr = calc_m_v(
            atom_site_sc_fract,
            flags_atom_fract_xyz,
            flag_m=False,
            flag_v=False,
        )[0].astype(bool)

        ddict["atom_label"] = atom_label
        ddict["atom_fract_xyz"] = atom_fract_xyz
        ddict["flags_atom_fract_xyz"] = flags_atom_fract_xyz_constr
        ddict["atom_scat_length_neutron"] = atom_b_scat
        ddict["atom_occupancy"] = atom_occupancy
        ddict["flags_atom_occupancy"] = (
            atom_site.get_flags_atom_occupancy()
        )
        ddict["flags_atom_b_iso"] = atom_site.get_flags_atom_b_iso()
        ddict["atom_symm_elems_flag"] = flag_2d

    if atom_site_susceptibility is not None:
        atom_para_label = ddict["atom_para_label"]
        atom_para_chi_type = numpy.array(
            atom_site_susceptibility.chi_type, dtype=str
        )
        if "atom_label" in ddict.keys():
            flag_2d = numpy.expand_dims(
                ddict["atom_label"], axis=1
            ) == numpy.expand_dims(atom_para_label, axis=0)
            args = numpy.argwhere(flag_2d)
            atom_para_index = args[args[:, 1].argsort()][:, 0]
            ddict["atom_para_index"] = atom_para_index

            if "full_mcif_elems" in ddict.keys():
                elems_fs = ddict["full_mcif_elems"]
            elif "full_symm_elems" in ddict.keys():
                elems_fs = ddict["full_symm_elems"]
            else:
                raise AttributeError("Symmetry elements are not found.")
            flag_2d = ddict["atom_symm_elems_flag"][:, atom_para_index]

            l_sc_chi = []
            for ind in range(atom_para_label.size):
                if (atom_para_chi_type[ind]).lower() == "ciso":
                    sc_chi = numpy.array(
                        [
                            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        ],
                        dtype=float,
                    )
                else:
                    sc_chi = calc_sc_chi(
                        elems_fs[:, flag_2d[:, ind]],
                        unit_cell_parameters=unit_cell_parameters,
                        flag_unit_cell_parameters=False,
                    )[0]
                l_sc_chi.append(sc_chi)
            atom_para_sc_chi = numpy.stack(l_sc_chi, axis=-1)
            ddict["atom_para_sc_chi"] = atom_para_sc_chi

    if atom_site_exchange is not None:
        atom_exchange_label = numpy.array(
            atom_site_exchange.label, dtype=str
        )
        atom_exchange_type = numpy.array(
            atom_site_exchange.j_type, dtype=str
        )
        atom_exchange_tensor = numpy.array(
            [
                atom_site_exchange.j_11,
                atom_site_exchange.j_22,
                atom_site_exchange.j_33,
                atom_site_exchange.j_12,
                atom_site_exchange.j_13,
                atom_site_exchange.j_23,
            ],
            dtype=float,
        )
        flags_atom_exchange_tensor = numpy.array(
            [
                atom_site_exchange.j_11_refinement,
                atom_site_exchange.j_22_refinement,
                atom_site_exchange.j_33_refinement,
                atom_site_exchange.j_12_refinement,
                atom_site_exchange.j_13_refinement,
                atom_site_exchange.j_23_refinement,
            ],
            dtype=bool,
        )
        ddict["atom_exchange_label"] = atom_exchange_label
        ddict["atom_exchange_type"] = atom_exchange_type
        ddict["atom_exchange_tensor"] = atom_exchange_tensor
        ddict["flags_atom_exchange_tensor"] = flags_atom_exchange_tensor
        if "atom_label" in ddict.keys():

            if "full_mcif_elems" in ddict.keys():
                elems_fs = ddict["full_mcif_elems"]
            elif "full_symm_elems" in ddict.keys():
                elems_fs = ddict["full_symm_elems"]
            else:
                raise AttributeError("Symmetry elements are not found.")

            flag_2d = numpy.expand_dims(
                ddict["atom_label"], axis=1
            ) == numpy.expand_dims(atom_exchange_label, axis=0)
            args = numpy.argwhere(flag_2d)
            atom_exchange_index = args[args[:, 1].argsort()][:, 0]
            ddict["atom_exchange_index"] = atom_exchange_index

            l_sc_exchange = []
            flag_2d = ddict["atom_symm_elems_flag"][:, atom_exchange_index]
            for ind in range(atom_exchange_label.size):
                if (atom_exchange_type[ind]).lower() == "jiso":
                    sc_exchange = numpy.array(
                        [
                            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        ],
                        dtype=float,
                    )
                else:
                    sc_exchange = calc_sc_chi(
                        elems_fs[:, flag_2d[:, ind]],
                        unit_cell_parameters=unit_cell_parameters,
                        flag_unit_cell_parameters=False,
                    )[0]
                l_sc_exchange.append(sc_exchange)
            atom_sc_exchange = numpy.stack(l_sc_exchange, axis=-1)
            ddict["atom_sc_exchange"] = atom_sc_exchange

    if atom_site_moment is not None:
        atom_ordered_label = numpy.array(atom_site_moment.label, dtype=str)

        atom_ordered_moment_crystalaxis_dn = numpy.array(
            [
                atom_site_moment.crystalaxis_x,
                atom_site_moment.crystalaxis_y,
                atom_site_moment.crystalaxis_z,
            ],
            dtype=float,
        )

        flags_atom_ordered_moment_crystalaxis_dn = numpy.array(
            [
                atom_site_moment.crystalaxis_x_refinement,
                atom_site_moment.crystalaxis_y_refinement,
                atom_site_moment.crystalaxis_z_refinement,
            ],
            dtype=bool,
        )
        ddict["atom_ordered_label"] = atom_ordered_label
        ddict["atom_ordered_moment_crystalaxis_dn"] = (
            atom_ordered_moment_crystalaxis_dn
        )
        ddict["flags_atom_ordered_moment_crystalaxis_dn"] = (
            flags_atom_ordered_moment_crystalaxis_dn
        )
        if "atom_label" in ddict.keys():
            flag_2d = numpy.expand_dims(
                ddict["atom_label"], axis=1
            ) == numpy.expand_dims(atom_ordered_label, axis=0)
            args = numpy.argwhere(flag_2d)
            atom_ordered_index = args[args[:, 1].argsort()][:, 0]
            ddict["atom_ordered_index"] = atom_ordered_index

            if "full_mcif_elems" in ddict.keys():
                elems_fs = ddict["full_mcif_elems"]
            else:
                raise AttributeError(
                    "Magnetic symmetry elements are not found."
                )
            flag_2d = ddict["atom_symm_elems_flag"][:, atom_ordered_index]

            l_sc_ordered = []
            for ind in range(atom_ordered_label.size):
                sc_ordered = calc_sc_ordered(
                    elems_fs[:, flag_2d[:, ind]],
                    unit_cell_parameters=unit_cell_parameters,
                    flag_unit_cell_parameters=False,
                )[0]
                l_sc_ordered.append(sc_ordered)
            atom_ordered_sc_m_dn = numpy.stack(l_sc_ordered, axis=-1)
            ddict["atom_ordered_sc_m"] = atom_ordered_sc_m_dn  # in dn
            ddict["atom_ordered_moment_crystalaxis_dn"] = calc_m_v(
                ddict["atom_ordered_sc_m"],
                ddict["atom_ordered_moment_crystalaxis_dn"],
                flag_m=False,
                flag_v=False,
            )[0]

    # FIXME: beta should be defined in atom_site_aniso_beta, b_iso for every atom
    atom_b_iso, atom_beta = crystal.calc_b_iso_beta()
    ddict["atom_b_iso"] = atom_b_iso
    ddict["atom_beta"] = atom_beta.transpose()
    ddict["flags_atom_beta"] = crystal.get_flags_atom_beta()

    if atom_site_scat is not None:
        mag_atom_label = numpy.array(atom_site_scat.label, dtype=str)
        mag_atom_lande_factor = numpy.array(
            atom_site_scat.lande, dtype=float
        )
        mag_atom_kappa = numpy.array(atom_site_scat.kappa, dtype=float)
        flags_mag_atom_lande_factor = atom_site_scat.get_flags_lande()
        flags_mag_atom_kappa = atom_site_scat.get_flags_kappa()

        ddict["mag_atom_label"] = mag_atom_label
        ddict["mag_atom_lande_factor"] = mag_atom_lande_factor
        ddict["mag_atom_kappa"] = mag_atom_kappa
        ddict["flags_mag_atom_lande_factor"] = flags_mag_atom_lande_factor
        ddict["flags_mag_atom_kappa"] = flags_mag_atom_kappa
        # l_ats = atom_site_scat.atom_type_scat
        if "atom_label" in ddict.keys():
            flag_2d = numpy.expand_dims(
                ddict["atom_label"], axis=1
            ) == numpy.expand_dims(mag_atom_label, axis=0)
            args = numpy.argwhere(flag_2d)
            mag_atom_index = args[args[:, 1].argsort()][:, 0]
            ddict["mag_atom_index"] = mag_atom_index

            if atom_site_scat.is_attribute("type_symbol"):
                mag_type_symbol = numpy.array(
                    atom_site_scat.type_symbol, dtype=str
                )
            else:
                as_tsymbol = numpy.array(atom_site.type_symbol, dtype=str)
                mag_type_symbol = as_tsymbol[mag_atom_index]

            j0_parameters, j2_parameters = get_j0_j2_parameters(
                mag_type_symbol
            )
            ddict["mag_atom_j0_parameters"] = j0_parameters
            ddict["mag_atom_j2_parameters"] = j2_parameters

        if "atom_para_label" in ddict.keys():
            flag_2d = numpy.expand_dims(
                ddict["atom_para_label"], axis=0
            ) == numpy.expand_dims(mag_atom_label, axis=1)
            args = numpy.argwhere(flag_2d)
            mag_atom_para_index = args[args[:, 1].argsort()][:, 0]
            ddict["mag_atom_para_index"] = mag_atom_para_index

        if "atom_ordered_label" in ddict.keys():
            flag_2d = numpy.expand_dims(
                ddict["atom_ordered_label"], axis=0
            ) == numpy.expand_dims(mag_atom_label, axis=1)
            args = numpy.argwhere(flag_2d)
            mag_atom_ordered_index = args[args[:, 1].argsort()][:, 0]
            ddict["mag_atom_ordered_index"] = mag_atom_ordered_index

        if "atom_rho_multipole_label" in ddict.keys():
            flag_2d = numpy.expand_dims(
                ddict["atom_rho_multipole_label"], axis=0
            ) == numpy.expand_dims(mag_atom_label, axis=1)
            args = numpy.argwhere(flag_2d)
            mag_atom_rho_multipole_index = args[args[:, 1].argsort()][:, 0]
            ddict["mag_atom_rho_multipole_index"] = (
                mag_atom_rho_multipole_index
            )

    if atom_type is not None:
        pass

    if atom_type_scat is not None:
        pass

    if atom_site_aniso is not None:
        atom_site_aniso_label = numpy.array(
            atom_site_aniso.label, dtype=str
        )
        flag_2d = numpy.expand_dims(
            ddict["atom_label"], axis=1
        ) == numpy.expand_dims(atom_site_aniso_label, axis=0)
        args = numpy.argwhere(flag_2d)
        atom_site_aniso_index = args[args[:, 1].argsort()][:, 0]
        ddict["atom_site_aniso_index"] = atom_site_aniso_index
        # FIXME: add transformation U_ani and others to beta

        if "full_mcif_elems" in ddict.keys():
            elems_fs = ddict["full_mcif_elems"]
        elif "full_symm_elems" in ddict.keys():
            elems_fs = ddict["full_symm_elems"]
        else:
            raise AttributeError("Symmetry elements are not found.")

        flag_2d = ddict["atom_symm_elems_flag"][:, atom_site_aniso_index]
        l_sc_beta = []
        atom_site_aniso_label = numpy.array(
            atom_site_aniso.label, dtype=str
        )

        for ind in range(atom_site_aniso_label.size):
            symm_elems = elems_fs[:, flag_2d[:, ind]]
            sc_beta = calc_sc_beta(symm_elems)
            l_sc_beta.append(sc_beta)
        ddict["atom_site_aniso_sc_beta"] = numpy.stack(l_sc_beta, axis=-1)

    if atom_electron_configuration is not None:
        d_slater_orbitals = DATABASE["Orbitals by radial Slater functions"]
        for aec_item in atom_electron_configuration.items:
            aec_label = aec_item.label
            at_symbol = atom_site[aec_label].type_symbol
            element_symbol, charge, isotope_number = (
                get_atomic_symbol_ion_charge_isotope_number_by_ion_symbol(
                    at_symbol
                )
            )
            core_shells_populations = (
                aec_item.get_core_shells_populations()
            )
            l_population = []
            l_n = []
            l_zeta = []
            l_coeff = []
            for sh_pop in core_shells_populations:
                shell, population = sh_pop
                n, dzeta, coeff = d_slater_orbitals[
                    (element_symbol, shell)
                ]
                l_population.append(population)
                l_n.append(n)
                l_zeta.append(dzeta)
                l_coeff.append(coeff)
            if len(l_population) > 0:
                dict_shell = {
                    "core_population": l_population,
                    "core_n": l_n,
                    "core_zeta": l_zeta,
                    "core_coeff": l_coeff,
                }
                ddict[f"shell_{aec_label:}"] = dict_shell

    if atom_local_axes is not None:
        l_hh = []
        for item_ala in atom_local_axes.items:
            label_ax1 = item_ala.ax1.lower()
            label_ax2 = item_ala.ax2.lower()
            index_atom_origin_ax1 = numpy.squeeze(
                numpy.argwhere(ddict["atom_label"] == item_ala.atom_label)
            )
            index_atom_end_ax1 = numpy.squeeze(
                numpy.argwhere(ddict["atom_label"] == item_ala.atom0)
            )
            index_atom_origin_ax2 = numpy.squeeze(
                numpy.argwhere(ddict["atom_label"] == item_ala.atom1)
            )
            index_atom_end_ax2 = numpy.squeeze(
                numpy.argwhere(ddict["atom_label"] == item_ala.atom2)
            )
            matrix_m = calc_m_m_by_unit_cell_parameters(
                ddict["unit_cell_parameters"]
            )[0]
            atom_local_axes_matrix = (
                calc_transformation_matrix_from_global_to_local(
                    index_atom_origin_ax1=index_atom_origin_ax1,
                    index_atom_end_ax1=index_atom_end_ax1,
                    label_ax1=label_ax1,
                    index_atom_origin_ax2=index_atom_origin_ax2,
                    index_atom_end_ax2=index_atom_end_ax2,
                    label_ax2=label_ax2,
                    atom_fract_xyz=atom_fract_xyz,
                    matrix_m=matrix_m,
                )
            )
            l_hh.append(atom_local_axes_matrix)
        ddict["atom_local_axes_matrix"] = numpy.stack(l_hh, axis=1)

    if atom_rho_multipole is not None:
        pass

    return ddict



def take_crystal_parameters_from_dictionary(crystal, ddict_crystal, l_parameter_name: list = None, l_sigma: list = None):
    """Body of Crystal.take_parameters_from_dictionary, see its docstring."""
    keys = ddict_crystal.keys()
    if l_parameter_name is None:
        l_parameter_name = []
        l_sigma = []

    if "atom_para_susceptibility" in keys:
        mag_atom_susceptibility = (
            ddict_crystal["atom_para_sc_chi"]
            * ddict_crystal["atom_para_susceptibility"][na, :, :]
        ).sum(axis=1)
        for i_item, item_ass in enumerate(
            crystal.atom_site_susceptibility.items
        ):
            item_ass.chi_11 = mag_atom_susceptibility[0, i_item]
            item_ass.chi_22 = mag_atom_susceptibility[1, i_item]
            item_ass.chi_33 = mag_atom_susceptibility[2, i_item]
            item_ass.chi_12 = mag_atom_susceptibility[3, i_item]
            item_ass.chi_13 = mag_atom_susceptibility[4, i_item]
            item_ass.chi_23 = mag_atom_susceptibility[5, i_item]

    if "atom_exchange_tensor" in keys:
        j_tensor = ddict_crystal["atom_exchange_tensor"]
        atom_sc_exchange = ddict_crystal["atom_sc_exchange"]
        j_tensor = numpy.sum(
            numpy.expand_dims(j_tensor, axis=0) * atom_sc_exchange, axis=1
        )
        for i_item, item_ase in enumerate(crystal.atom_site_exchange.items):
            item_ase.j_11 = j_tensor[0, i_item]
            item_ase.j_22 = j_tensor[1, i_item]
            item_ase.j_33 = j_tensor[2, i_item]
            item_ase.j_12 = j_tensor[3, i_item]
            item_ase.j_13 = j_tensor[4, i_item]
            item_ase.j_23 = j_tensor[5, i_item]

    if "atom_ordered_moment_crystalaxis_dn" in keys:
        hh = calc_m_v(
            ddict_crystal["atom_ordered_sc_m"],
            ddict_crystal["atom_ordered_moment_crystalaxis_dn"],
        )[0]
        # hh = ddict_crystal["atom_ordered_moment_crystalaxis_xyz"]
        for i_item, item in enumerate(crystal.atom_site_moment.items):
            item.crystalaxis_x = hh[0, i_item]
            item.crystalaxis_y = hh[1, i_item]
            item.crystalaxis_z = hh[2, i_item]

    if "mag_atom_lande_factor" in keys:
        hh = ddict_crystal["mag_atom_lande_factor"]
        for i_item, item in enumerate(crystal.atom_site_scat.items):
            item.lande = hh[i_item]

    if "mag_atom_kappa" in keys:
        hh = ddict_crystal["mag_atom_kappa"]
        for i_item, item in enumerate(crystal.atom_site_scat.items):
            item.kappa = hh[i_item]

    if "atom_fract_xyz" in keys:
        hh = ddict_crystal["atom_fract_xyz"]
        atom_site_sc_fract = ddict_crystal["atom_site_sc_fract"]
        atom_site_sc_b = ddict_crystal["atom_site_sc_b"]
        hh = (
            calc_m_v(
                atom_site_sc_fract,
                numpy.mod(hh, 1),
                flag_m=False,
                flag_v=False,
            )[0]
            + atom_site_sc_b
        )

        for i_item, item in enumerate(crystal.atom_site.items):
            item.fract_x = numpy.round(
                numpy.mod(hh[0, i_item], 1.0), decimals=6
            )
            item.fract_y = numpy.round(
                numpy.mod(hh[1, i_item], 1.0), decimals=6
            )
            item.fract_z = numpy.round(
                numpy.mod(hh[2, i_item], 1.0), decimals=6
            )

    if "atom_occupancy" in keys:
        hh = ddict_crystal["atom_occupancy"]
        for i_item, item in enumerate(crystal.atom_site.items):
            item.occupancy = numpy.round(hh[i_item], decimals=6)

    if "atom_multiplicity" in keys:
        hh = ddict_crystal["atom_multiplicity"]
        for i_item, item in enumerate(crystal.atom_site.items):
            item.multiplicity = hh[i_item]

    if "atom_b_iso" in keys:
        atom_b_iso = ddict_crystal["atom_b_iso"]
        if "atom_beta" in keys:
            atom_beta = ddict_crystal["atom_beta"]

            if "atom_site_aniso_sc_beta" in keys:
                atom_site_aniso_sc_beta = ddict_crystal[
                    "atom_site_aniso_sc_beta"
                ]
                atom_site_aniso_index = ddict_crystal[
                    "atom_site_aniso_index"
                ]
                atom_sc_beta = numpy.zeros(
                    (6,) + atom_beta.shape, dtype=float
                )
                atom_sc_beta[:, :, atom_site_aniso_index] = (
                    atom_site_aniso_sc_beta
                )
                atom_beta = (
                    atom_sc_beta * numpy.expand_dims(atom_beta, axis=0)
                ).sum(axis=1)

        else:
            atom_beta = numpy.zeros((6,) + atom_b_iso.shape, dtype=float)
        unit_cell_parameters = ddict_crystal["unit_cell_parameters"]
        if not (crystal.atom_site.is_attribute("adp_type")):
            for item in crystal.atom_site:
                item.adp_type = "Biso"
                if not (item.is_attribute("b_iso_or_equiv")):
                    item.b_iso_or_equiv = 0.0
        atom_adp_type = numpy.array(crystal.atom_site.adp_type, dtype=str)
        atom_iso_param, atom_aniso_param = (
            calc_param_iso_aniso_by_b_iso_beta(
                unit_cell_parameters, atom_adp_type, atom_b_iso, atom_beta
            )
        )

        if crystal.is_attribute("atom_site"):
            for i_item, item in enumerate(crystal.atom_site.items):
                if item.adp_type in ("Uiso", "Uovl", "Umpe"):
                    item.u_iso_or_equiv = atom_iso_param[i_item]
                elif item.adp_type in ("Biso", "Bovl"):
                    item.b_iso_or_equiv = atom_iso_param[i_item]

            if crystal.is_attribute("atom_site_aniso"):
                atom_site_aniso_index = ddict_crystal[
                    "atom_site_aniso_index"
                ]

                for item, index_as in zip(
                    crystal.atom_site_aniso.items, atom_site_aniso_index
                ):
                    aniso_param = atom_aniso_param[:, index_as]
                    item_as = crystal.atom_site.items[index_as]
                    if item_as.adp_type in ("Uani",):
                        item.u_11 = aniso_param[0]
                        item.u_22 = aniso_param[1]
                        item.u_33 = aniso_param[2]
                        item.u_12 = aniso_param[3]
                        item.u_13 = aniso_param[4]
                        item.u_23 = aniso_param[5]
                    elif item_as.adp_type in ("Bani",):
                        item.b_11 = aniso_param[0]
                        item.b_22 = aniso_param[1]
                        item.b_33 = aniso_param[2]
                        item.b_12 = aniso_param[3]
                        item.b_13 = aniso_param[4]
                        item.b_23 = aniso_param[5]

    if "unit_cell_parameters" in keys:
        hh = ddict_crystal["unit_cell_parameters"]
        sc_uc = ddict_crystal["sc_uc"]
        v_uc = ddict_crystal["v_uc"]
        hhh = numpy.dot(sc_uc, hh) + v_uc
        cell = crystal.cell
        cell.length_a = hhh[0]
        cell.length_b = hhh[1]
        cell.length_c = hhh[2]
        cell.angle_alpha = numpy.round(
            hhh[3] * 180.0 / numpy.pi, decimals=5
        )
        cell.angle_beta = numpy.round(
            hhh[4] * 180.0 / numpy.pi, decimals=5
        )
        cell.angle_gamma = numpy.round(
            hhh[5] * 180.0 / numpy.pi, decimals=5
        )

    if "atom_rho_multipole_kappa" in keys:
        hh = ddict_crystal["atom_rho_multipole_kappa"]
        for i_item, item in enumerate(crystal.atom_rho_multipole.items):
            item.kappa = numpy.round(hh[i_item], decimals=6)

    if "atom_rho_multipole_plm" in keys:
        hh = ddict_crystal["atom_rho_multipole_plm"]
        for i_item, item in enumerate(crystal.atom_rho_multipole.items):
            item.load_plm_from_dictionary(
                numpy.round(hh[:, i_item], decimals=6)
            )

    for name, sigma in zip(l_parameter_name, l_sigma):
        parameter_label, ind_s = name
        if "atom_para_susceptibility" in parameter_label:
            ind_chi, ind_a = ind_s[0], ind_s[1]
            item_ass = crystal.atom_site_susceptibility.items[ind_a]
            if ind_chi == 0:
                item_ass.chi_11_sigma = float(sigma)
            elif ind_chi == 1:
                item_ass.chi_22_sigma = float(sigma)
            elif ind_chi == 2:
                item_ass.chi_33_sigma = float(sigma)
            elif ind_chi == 3:
                item_ass.chi_12_sigma = float(sigma)
            elif ind_chi == 4:
                item_ass.chi_13_sigma = float(sigma)
            elif ind_chi == 5:
                item_ass.chi_23_sigma = float(sigma)

        if "atom_exchange_tensor" in parameter_label:
            ind_j, ind_a = ind_s[0], ind_s[1]
            item = crystal.atom_site_exchange.items[ind_a]
            if ind_j == 0:
                item.j_11_sigma = float(sigma)
            elif ind_j == 1:
                item.j_22_sigma = float(sigma)
            elif ind_j == 2:
                item.j_33_sigma = float(sigma)
            elif ind_j == 3:
                item.j_12_sigma = float(sigma)
            elif ind_j == 4:
                item.j_13_sigma = float(sigma)
            elif ind_j == 5:
                item.j_23_sigma = float(sigma)

        if "atom_fract_xyz" in parameter_label:
            ind_p, ind_a = ind_s[0], ind_s[1]
            item = crystal.atom_site.items[ind_a]
            if ind_p == 0:
                item.fract_x_sigma = float(sigma)
            elif ind_p == 1:
                item.fract_y_sigma = float(sigma)
            elif ind_p == 2:
                item.fract_z_sigma = float(sigma)
        if "atom_occupancy" in parameter_label:
            ind_a = ind_s[0]
            item = crystal.atom_site.items[ind_a]
            item.occupancy_sigma = float(sigma)

        if "atom_b_iso" in parameter_label:
            ind_a = ind_s[0]
            item = crystal.atom_site.items[ind_a]
            if item.adp_type in ("Uiso", "Uovl", "Umpe"):
                item.u_iso_or_equiv_sigma = float(sigma) / float(
                    8.0 * numpy.square(numpy.pi)
                )
            elif item.adp_type in ("Biso", "Bovl", "Bmpe"):
                item.b_iso_or_equiv_sigma = float(sigma)

        if "atom_beta" in parameter_label:
            ind_p, ind_a = ind_s[0], ind_s[1]
            item = crystal.atom_site_aniso[crystal.atom_site.items[ind_a].label]
            if crystal.atom_site[item.label].adp_type in ("Bani",):
                if ind_p == 0:
                    item.b_11_sigma = float(sigma)
                elif ind_p == 1:
                    item.b_22_sigma = float(sigma)
                elif ind_p == 2:
                    item.b_33_sigma = float(sigma)
                elif ind_p == 3:
                    item.b_12_sigma = float(sigma)
                elif ind_p == 4:
                    item.b_13_sigma = float(sigma)
                elif ind_p == 5:
                    item.b_23_sigma = float(sigma)

            if crystal.atom_site[item.label].adp_type in ("Uani",):
                reciprocal_unit_cell_parameters = (
                    calc_reciprocal_by_unit_cell_parameters(
                        unit_cell_parameters,
                        flag_unit_cell_parameters=False,
                    )[0]
                )
                beta_hh = numpy.zeros((6, 1), dtype=float)
                coeff_hh = calc_u_ij_by_beta(
                    beta_hh,
                    reciprocal_unit_cell_parameters,
                    flag_beta=True,
                )[1]["beta"]
                if ind_p == 0:
                    item.u_11_sigma = float(sigma * coeff_hh[0].squeeze())
                elif ind_p == 1:
                    item.u_22_sigma = float(sigma * coeff_hh[1].squeeze())
                elif ind_p == 2:
                    item.u_33_sigma = float(sigma * coeff_hh[2].squeeze())
                elif ind_p == 3:
                    item.u_12_sigma = float(sigma * coeff_hh[3].squeeze())
                elif ind_p == 4:
                    item.u_13_sigma = float(sigma * coeff_hh[4].squeeze())
                elif ind_p == 5:
                    item.u_23_sigma = float(sigma * coeff_hh[5].squeeze())

        if "atom_ordered_moment_crystalaxis_dn" in parameter_label:
            ind_p, ind_a = ind_s[0], ind_s[1]
            item = crystal.atom_site_moment.items[ind_a]
            if ind_p == 0:
                item.crystalaxis_x_sigma = float(sigma)
            elif ind_p == 1:
                item.crystalaxis_y_sigma = float(sigma)
            elif ind_p == 2:
                item.crystalaxis_z_sigma = float(sigma)

        if "mag_atom_kappa" in parameter_label:
            ind_a = ind_s[0]
            item = crystal.atom_site_scat.items[ind_a]
            item.kappa_sigma = float(sigma)

        if "mag_atom_lande_factor" in parameter_label:
            ind_a = ind_s[0]
            item = crystal.atom_site_scat.items[ind_a]
            item.lande_sigma = float(sigma)

        if "unit_cell_parameters" in parameter_label:
            ind_p = ind_s[0]
            item = crystal.cell
            if ind_p == 0:
                item.length_a_sigma = float(sigma)
            elif ind_p == 1:
                item.length_b_sigma = float(sigma)
            elif ind_p == 2:
                item.length_c_sigma = float(sigma)
            elif ind_p == 3:
                item.angle_alpha_sigma = float(sigma) * 180.0 / numpy.pi
            elif ind_p == 4:
                item.angle_beta_sigma = float(sigma) * 180.0 / numpy.pi
            elif ind_p == 5:
                item.angle_gamma_sigma = float(sigma) * 180.0 / numpy.pi

        if "atom_rho_multipole_kappa" in parameter_label:
            ind_a = ind_s[0]
            item = crystal.atom_rho_multipole.items[ind_a]
            item.kappa_sigma = float(sigma)

        if "atom_rho_multipole_plm" in parameter_label:
            ind_p = ind_s[0]
            ind_a = ind_s[1]
            item = crystal.atom_rho_multipole.items[ind_a]
            item.set_plm_sigma_by_index(ind_p, float(sigma))



def save_crystal_multipole_density_to_den(crystal, file_den: str = 'file.den', Na: int = 48, Nb: int = 48, Nc: int = 48):
    """Body of Crystal.save_atom_rho_multipole_density_to_den, see its docstring."""
    crystal.form_object()
    n_abc = numpy.array([Na, Nb, Nc], dtype=int)
    d_crystal = crystal.get_dictionary()
    unit_cell_parameters = d_crystal["unit_cell_parameters"]

    if "full_mcif_elems" in d_crystal.keys():
        full_mcif_elems = d_crystal["full_mcif_elems"]
        full_symm_elems = full_mcif_elems  # [:13,:]
    else:
        full_symm_elems = d_crystal["full_symm_elems"]

    atom_fract_xyz = d_crystal["atom_fract_xyz"]
    atom_label = d_crystal["atom_label"]
    atom_rho_multipole_label = d_crystal["atom_rho_multipole_label"]

    atom_local_axes_matrix = d_crystal["atom_local_axes_matrix"]
    atom_local_axes_label = d_crystal["atom_local_axes_label"]

    l_hh_1, l_hh_2 = [], []
    for label in atom_rho_multipole_label:
        ind_1 = numpy.squeeze(numpy.argwhere(atom_label == label))
        ind_2 = numpy.squeeze(
            numpy.argwhere(atom_local_axes_label == label)
        )
        l_hh_1.append(atom_fract_xyz[..., ind_1])
        l_hh_2.append(atom_local_axes_matrix[..., ind_2])
    atom_rho_multipole_fract_xyz = numpy.stack(l_hh_1, axis=-1)
    atom_rho_multipole_transformation_matrix = numpy.stack(l_hh_2, axis=-1)

    index_auc, point_multiplicity = calc_asymmetric_unit_cell_indexes(
        n_abc, full_symm_elems
    )
    symm_elem_auc = calc_symm_elem_points_by_index_points(index_auc, n_abc)

    point_atom_distance_auc, point_atom_symm_elems_auc = form_basins_2(
        symm_elem_auc,
        full_symm_elems,
        unit_cell_parameters,
        atom_rho_multipole_fract_xyz,
    )

    atom_rho_multipole_plm = d_crystal["atom_rho_multipole_plm"]
    atom_rho_multipole_lm = d_crystal["atom_rho_multipole_lm"]
    l_kappa = d_crystal["atom_rho_multipole_kappa"]

    l_n, l_zeta, l_coeff = d_crystal["atom_rho_multipole_n_zeta_coeff"]
    multipole_density = calc_multipole_density(
        point_atom_symm_elems_auc,
        symm_elem_auc,
        unit_cell_parameters,
        atom_rho_multipole_fract_xyz,
        atom_rho_multipole_transformation_matrix,
        atom_rho_multipole_lm,
        atom_rho_multipole_plm,
        l_n,
        l_zeta,
        l_coeff,
        l_kappa,
    )
    logger.info(f"file den is '{file_den}'")
    if file_den != "":
        hh = numpy.sum(multipole_density, axis=1)
        flag_pos = hh >= 0.0
        spin_density = numpy.zeros(shape=(2,) + hh.shape, dtype=float)
        spin_density[0, flag_pos] = hh[flag_pos]
        spin_density[1, ~flag_pos] = hh[~flag_pos]
        if "reduced_symm_elems" in d_crystal.keys():
            reduced_symm_elems = d_crystal["reduced_symm_elems"]
            centrosymmetry = d_crystal["centrosymmetry"]
            centrosymmetry_position = d_crystal["centrosymmetry_position"]
            translation_elems = d_crystal["translation_elems"]
        else:
            reduced_symm_elems = full_symm_elems
            centrosymmetry = False
            centrosymmetry_position = numpy.array([0, 0, 0, 1])
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
                ]
            )
        save_spin_density_into_file(
            file_den,
            index_auc,
            spin_density,
            n_abc,
            unit_cell_parameters,
            reduced_symm_elems,
            translation_elems,
            centrosymmetry,
            centrosymmetry_position,
        )
    return (
        index_auc,
        point_atom_symm_elems_auc,
        point_multiplicity,
        multipole_density,
        unit_cell_parameters,
    )



def save_crystal_multipole_magnetization_density_to_den(crystal, file_den: str = 'file.den', Na: int = 48, Nb: int = 48, Nc: int = 48):
    """Body of Crystal.save_atom_rho_multipole_magnetization_density_to_den, see its docstring."""
    n_abc = numpy.array([Na, Nb, Nc], dtype=int)
    (
        index_auc,
        point_atom_symm_elems_auc,
        point_multiplicity,
        multipole_density,
        unit_cell_parameters,
    ) = crystal.save_atom_rho_multipole_density_to_den(
        file_den=file_den, Na=Na, Nb=Nb, Nc=Nc
    )
    d_crystal = crystal.get_dictionary()
    flag_moment_ordered = (
        "atom_ordered_moment_crystalaxis_dn" in d_crystal.keys()
    )
    if flag_moment_ordered:
        full_mcif_elems = d_crystal["full_mcif_elems"]
        atom_ordered_m_dn = d_crystal["atom_ordered_moment_crystalaxis_dn"]
        index_full, moment_3d = transform_to_density_auc_to_moments(
            index_auc,
            point_atom_symm_elems_auc,
            point_multiplicity,
            multipole_density,
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
            if file_den != "":
                file_magnetization_density_component = (
                    ".".join(file_den.split(".")[:-1]) + f"_{lab:}.den"
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


def calc_sc_v_unit_cell_parameters(
    type_cell: str, it_coordinate_system_code: str
):
    """
    Calculate symmetry constraints on unit_cell_parameters as
    unit_cell_parameters = sc_uc * unit_cell_parameters + v_uc
    """
    sc_uc = numpy.eye(6, dtype=float)
    v_uc = numpy.zeros((6,), dtype=float)

    if type_cell == "aP":
        pass
    elif type_cell is None:
        pass
    elif type_cell.startswith("m"):
        sc_uc[3, 3] = 0.0
        sc_uc[5, 5] = 0.0
        v_uc[3] = numpy.pi / 2.0
        v_uc[5] = numpy.pi / 2.0
    elif type_cell.startswith("o"):
        sc_uc[3, 3] = 0.0
        sc_uc[4, 4] = 0.0
        sc_uc[5, 5] = 0.0
        v_uc[3] = numpy.pi / 2.0
        v_uc[4] = numpy.pi / 2.0
        v_uc[5] = numpy.pi / 2.0
    elif type_cell.startswith("t"):  # FIXME: check  | (type_cell == "hP")
        sc_uc[0, 0], sc_uc[0, 1] = 1.0, 0.0
        sc_uc[1, 0], sc_uc[1, 1] = 1.0, 0.0
        sc_uc[3, 3] = 0.0
        sc_uc[4, 4] = 0.0
        sc_uc[5, 5] = 0.0
        v_uc[3] = numpy.pi / 2.0
        v_uc[4] = numpy.pi / 2.0
        v_uc[5] = numpy.pi / 2.0
    elif type_cell.startswith("hP"):
        if it_coordinate_system_code.lower() == "h":
            sc_uc[0, 0], sc_uc[0, 1] = 1.0, 0.0
            sc_uc[1, 0], sc_uc[1, 1] = 1.0, 0.0
            sc_uc[3, 3] = 0.0
            sc_uc[4, 4] = 0.0
            sc_uc[5, 5] = 0.0
            v_uc[3] = numpy.pi / 2.0
            v_uc[4] = numpy.pi / 2.0
            v_uc[5] = numpy.pi * 2.0 / 3.0
        else:
            sc_uc[0, 0], sc_uc[0, 1], sc_uc[0, 2] = 1.0, 0.0, 0.0
            sc_uc[1, 0], sc_uc[1, 1], sc_uc[1, 2] = 1.0, 0.0, 0.0
            sc_uc[2, 0], sc_uc[2, 1], sc_uc[2, 2] = 1.0, 0.0, 0.0
            sc_uc[3, 3], sc_uc[3, 4], sc_uc[3, 5] = 1.0, 0.0, 0.0
            sc_uc[4, 3], sc_uc[4, 4], sc_uc[4, 5] = 1.0, 0.0, 0.0
            sc_uc[5, 3], sc_uc[5, 4], sc_uc[5, 5] = 1.0, 0.0, 0.0
    elif type_cell == "hR":
        if it_coordinate_system_code.lower() == "h":
            sc_uc[0, 0], sc_uc[0, 1] = 1.0, 0.0
            sc_uc[1, 0], sc_uc[1, 1] = 1.0, 0.0
            sc_uc[3, 3] = 0.0
            sc_uc[4, 4] = 0.0
            sc_uc[5, 5] = 0.0
            v_uc[3] = numpy.pi / 2.0
            v_uc[4] = numpy.pi / 2.0
            v_uc[5] = numpy.pi * 2.0 / 3.0
        else:
            sc_uc[0, 0], sc_uc[0, 1], sc_uc[0, 2] = 1.0, 0.0, 0.0
            sc_uc[1, 0], sc_uc[1, 1], sc_uc[1, 2] = 1.0, 0.0, 0.0
            sc_uc[2, 0], sc_uc[2, 1], sc_uc[2, 2] = 1.0, 0.0, 0.0
            sc_uc[3, 3], sc_uc[3, 4], sc_uc[3, 5] = 1.0, 0.0, 0.0
            sc_uc[4, 3], sc_uc[4, 4], sc_uc[4, 5] = 1.0, 0.0, 0.0
            sc_uc[5, 3], sc_uc[5, 4], sc_uc[5, 5] = 1.0, 0.0, 0.0
    elif type_cell.startswith("c"):
        sc_uc[0, 0], sc_uc[0, 1], sc_uc[0, 2] = 1.0, 0.0, 0.0
        sc_uc[1, 0], sc_uc[1, 1], sc_uc[1, 2] = 1.0, 0.0, 0.0
        sc_uc[2, 0], sc_uc[2, 1], sc_uc[2, 2] = 1.0, 0.0, 0.0
        sc_uc[3, 3] = 0.0
        sc_uc[4, 4] = 0.0
        sc_uc[5, 5] = 0.0
        v_uc[3] = numpy.pi / 2.0
        v_uc[4] = numpy.pi / 2.0
        v_uc[5] = numpy.pi / 2.0

    return sc_uc, v_uc
