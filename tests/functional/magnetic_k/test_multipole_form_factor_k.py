"""Multipole (anisotropic) magnetic form factor of a structure with a
propagation vector (magnetic superspace group).

1. k = (1/2, 0, 0) with a two-fold axis: the powder pattern of the
   k-structure equals that of the magnetic supercell written explicitly in
   P1, where every generated atom has its own multipoles with local axes
   rotated by the operation that generates it (independent real-space
   construction, the mcif code for ordered moments).
2. MnF2 (flip ratios, multipoles of Mn and F): the magnetic space group and
   the same structure converted to the superspace (k = 0) give the same
   chi_sq.
"""
import itertools
import os

import numpy
import pytest

import cryspy
from cryspy.procedure_rhochi.rhochi_by_dictionary import (
    rhochi_calc_chi_sq_by_dictionary,
)
from cryspy.procedure_rhochi.rhochi_diffrn import (
    calc_chi_sq_for_diffrn_by_dictionary,
)
from cryspy.exceptions import CrysPyError
from cryspy.C_item_loop_classes.cl_1_atom_local_axes import AtomLocalAxesL
from cryspy.C_item_loop_classes.cl_1_atom_rho_multipole import AtomRhoMultipoleL

CELL = numpy.array([3.9, 5.1, 6.3])
K_VECTOR = numpy.array([0.5, 0.0, 0.0])
N_CELL = numpy.array([2, 1, 1])
# identity and the two-fold axis along a (keeps k); only proper rotations,
# so the local axes of the dummy atoms rotate as the density does
OPERATIONS = [("x1,x2,x3,x4,+1", numpy.eye(3)),
              ("x1,-x2,-x3,x4,+1", numpy.diag([1., -1., -1.]))]
# label, position, Fourier cos components along a, b, c (μB), multipoles
ATOMS = [
    ("Fe1", (0.11, 0.23, 0.37), (1.2, -0.7, 2.1), "P20 0.3, P21 0.15, P22 0.2, P2-1 -0.1"),
    ("Fe2", (0.61, 0.08, 0.19), (-0.4, 1.9, 0.5), None),
]
DUMMY_DISTANCE = 0.5  # Å, as in add_multipole_form_factor

EXPERIMENT = """
data_pd

_setup_wavelength 1.4
_setup_field 0.0
_setup_offset_2theta 0.0
_setup_radiation neutrons
_setup_K 0.0
_setup_cthm 0.0

_pd_instr_resolution_U 0.3
_pd_instr_resolution_V -0.2
_pd_instr_resolution_W 0.2
_pd_instr_resolution_X 0.0
_pd_instr_resolution_Y 0.0

_range_2theta_min 5.0
_range_2theta_max 120.0

loop_
_phase_label
_phase_scale
_phase_igsize
phase 1.0 0.0

loop_
_pd_background_2theta
_pd_background_intensity
5.0 0.0
120.0 0.0

loop_
_pd_meas_2theta
_pd_meas_intensity
_pd_meas_intensity_sigma
50.0 1.0 1.0
"""

ATOM_SITE_HEADER = (
    "loop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n"
    "_atom_site_fract_y\n_atom_site_fract_z\n_atom_site_adp_type\n"
    "_atom_site_B_iso_or_equiv\n_atom_site_occupancy")


def form_cell(lengths) -> list:
    return [f"_cell_length_{n} {v}" for n, v in zip("abc", lengths)] + [
        "_cell_angle_alpha 90.0", "_cell_angle_beta 90.0", "_cell_angle_gamma 90.0", ""]


def form_rcif_k() -> str:
    lines = ["global_\n\ndata_phase\n"] + form_cell(CELL) + [ATOM_SITE_HEADER]
    lines += [f"{lab} Fe {x} {y} {z} Biso 0.0 1.0" for lab, (x, y, z), _, _ in ATOMS]
    lines += ["", "loop_\n_atom_site_scat_label\n_atom_site_scat_type_symbol\n"
              "_atom_site_scat_lande"]
    lines += [f"{lab} Fe3+ 2.0" for lab, _, _, _ in ATOMS]
    lines += ["", "loop_\n_cell_wave_vector_seq_id\n_cell_wave_vector_x\n"
              "_cell_wave_vector_y\n_cell_wave_vector_z",
              "1 {} {} {}".format(*K_VECTOR), "",
              "loop_\n_space_group_symop_magn_ssg_operation.id\n"
              "_space_group_symop_magn_ssg_operation.algebraic"]
    lines += [f"{i + 1} {text}" for i, (text, _) in enumerate(OPERATIONS)]
    lines += ["", "loop_\n_space_group_symop_magn_ssg_centering.id\n"
              "_space_group_symop_magn_ssg_centering.algebraic", "1 x1,x2,x3,x4,+1", "",
              "loop_\n_atom_site_moment_Fourier.atom_site_label\n"
              "_atom_site_moment_Fourier.axis\n"
              "_atom_site_moment_Fourier.wave_vector_seq_id\n"
              "_atom_site_moment_Fourier_param.cos\n"
              "_atom_site_moment_Fourier_param.sin"]
    for lab, _, cos, _ in ATOMS:
        lines += [f"{lab} {axis} 1 {value} 0.0" for axis, value in zip("xyz", cos)]
    return "\n".join(lines) + "\n" + EXPERIMENT


def form_crystal_k():
    obj = cryspy.str_to_globaln(form_rcif_k())
    for lab, _, _, coefficients in ATOMS:
        if coefficients is not None:
            obj.crystal_phase.add_multipole_form_factor(lab, coefficients=coefficients)
    return obj


def generate_sites(position, cos):
    """Sites (small-cell units), moments and rotations of an atom in the
    supercell by the real-space rule M(r) = R M(source) cos(2 pi k.r_source)."""
    lattice = numpy.array(list(itertools.product(*[range(n) for n in N_CELL])))
    sites = {}
    for _, rotation in OPERATIONS:
        x = numpy.asarray(position)[numpy.newaxis, :] + lattice
        x_new = numpy.mod(x @ rotation.transpose(), N_CELL)
        moment = numpy.cos(2. * numpy.pi * (x @ K_VECTOR))[:, numpy.newaxis] * \
            numpy.asarray(cos)[numpy.newaxis, :] @ rotation.transpose()
        for site, site_moment in zip(x_new, moment):
            key = tuple(numpy.round(site, 6))
            if key not in sites:
                sites[key] = (site, site_moment, rotation)
    return list(sites.values())


def form_rcif_supercell() -> str:
    rows, scat, moments, axes, multipoles = [], [], [], [], []
    for lab, position, cos, coefficients in ATOMS:
        for i_site, (site, moment, rotation) in enumerate(generate_sites(position, cos)):
            name = f"{lab}_{i_site + 1}"
            xyz = site / N_CELL
            rows.append(f"{name} Fe {xyz[0]:.10f} {xyz[1]:.10f} {xyz[2]:.10f} Biso 0.0 1.0")
            scat.append(f"{name} Fe3+ 2.0")
            moments.append(f"{name} {moment[0]:.10f} {moment[1]:.10f} {moment[2]:.10f}")
            if coefficients is None:
                continue
            for axis, i_axis in (("c", 2), ("a", 0)):
                shift = numpy.zeros(3)
                shift[i_axis] = DUMMY_DISTANCE / CELL[i_axis]
                dummy = (site + rotation @ shift) / N_CELL
                rows.append(f"D{name}_{axis} Dummy {dummy[0]:.10f} {dummy[1]:.10f} "
                            f"{dummy[2]:.10f} Biso 0.0 1.0")
            axes.append(f"{name} D{name}_c Z {name} D{name}_a X")
            multipoles.append(name)
    obj_k = form_crystal_k()
    item = obj_k.crystal_phase.atom_rho_multipole["Fe1"]
    coeff = [name for name in item.ATTR_OPTIONAL_NAMES if item.is_attribute(name)]
    coeff_cif = [item.ATTR_OPTIONAL_CIF[item.ATTR_OPTIONAL_NAMES.index(name)] for name in coeff]
    lines = ["global_\n\ndata_phase\n"] + form_cell(CELL * N_CELL) + [
        "loop_\n_space_group_symop_magn_operation_xyz\n"
        "_space_group_symop_magn_operation_id", "x,y,z,+1 1", "",
        "loop_\n_space_group_symop_magn_centering_xyz\n"
        "_space_group_symop_magn_centering_id", "x,y,z,+1 1", "",
        ATOM_SITE_HEADER] + rows
    lines += ["", "loop_\n_atom_site_scat_label\n_atom_site_scat_type_symbol\n"
              "_atom_site_scat_lande"] + scat
    lines += ["", "loop_\n_atom_site_moment_label\n_atom_site_moment_crystalaxis_x\n"
              "_atom_site_moment_crystalaxis_y\n_atom_site_moment_crystalaxis_z"] + moments
    lines += ["", "loop_\n_atom_local_axes_atom_label\n_atom_local_axes_atom0\n"
              "_atom_local_axes_ax1\n_atom_local_axes_atom1\n_atom_local_axes_atom2\n"
              "_atom_local_axes_ax2"] + axes
    lines += ["", "loop_\n_atom_rho_multipole_atom_label\n_atom_rho_multipole_atom_type_symbol\n"
              "_atom_rho_multipole_orbital\n"
              + "\n".join(f"_atom_rho_multipole_{name}" for name in coeff_cif)]
    values = " ".join(f"{getattr(item, name)}" for name in coeff)
    lines += [f"{name} {item.atom_type_symbol} {item.orbital} {values}" for name in multipoles]
    return "\n".join(lines) + "\n" + EXPERIMENT


def calc_signal(obj, ttheta):
    dict_obj = obj.get_dictionary()
    dict_obj["pd_pd"]["ttheta"] = numpy.radians(ttheta)
    dict_obj["pd_pd"]["signal_exp"] = numpy.ones((2, ttheta.size))
    dict_obj["pd_pd"]["excluded_points"] = numpy.zeros(ttheta.size, dtype=bool)
    dict_in_out = {}
    rhochi_calc_chi_sq_by_dictionary(
        dict_obj, dict_in_out=dict_in_out, flag_use_precalculated_data=False,
        flag_calc_analytical_derivatives=False)
    pd = dict_in_out["pd_pd"]
    return pd["signal_plus"] + pd["signal_minus"]


def test_k_structure_with_multipoles_equals_supercell():
    ttheta = numpy.linspace(8., 115., 2000)
    signal_k = calc_signal(form_crystal_k(), ttheta)
    signal_s = calc_signal(cryspy.str_to_globaln(form_rcif_supercell()), ttheta)
    n_cell_sq = float(numpy.prod(N_CELL)) ** 2
    numpy.testing.assert_allclose(signal_s, n_cell_sq * signal_k,
                                  rtol=1e-6, atol=1e-6 * signal_s.max())


def test_multipoles_change_the_satellites():
    """The multipoles of Fe1 are used by the superspace path."""
    ttheta = numpy.linspace(8., 115., 2000)
    signal_m = calc_signal(form_crystal_k(), ttheta)
    signal_d = calc_signal(cryspy.str_to_globaln(form_rcif_k()), ttheta)
    assert numpy.abs(signal_m - signal_d).max() > 0.01 * signal_d.max()


MNF2 = os.path.join(os.path.dirname(__file__), "..", "MEM_test", "mnf2_moments",
                    "MnF2_multipole_prior.rcif")


def calc_chi_sq_mnf2(obj):
    dict_obj = obj.get_dictionary()
    return calc_chi_sq_for_diffrn_by_dictionary(
        dict_obj["diffrn_exp1"], dict_obj["crystal_mnf2"], {})[0]


def test_mnf2_multipoles_superspace_equals_magnetic_space_group():
    obj = cryspy.load_file(MNF2)
    chi_sq_msg = calc_chi_sq_mnf2(obj)
    obj.crystal_mnf2.convert_to_superspace()
    assert calc_chi_sq_mnf2(obj) == pytest.approx(chi_sq_msg, rel=1e-9)


def remove_multipoles(crystal):
    crystal.items[:] = [item for item in crystal.items if not isinstance(
        item, (AtomRhoMultipoleL, AtomLocalAxesL))]


def test_mnf2_multipoles_of_one_atom_only():
    """F keeps its dipole form factor, Mn gets multipoles by the method."""
    obj = cryspy.load_file(MNF2)
    crystal = obj.crystal_mnf2
    remove_multipoles(crystal)
    report = crystal.add_multipole_form_factor("Mn", coefficients="P20 0.05")
    assert "orbital 3d" in report
    chi_sq_msg = calc_chi_sq_mnf2(obj)
    crystal.convert_to_superspace()
    assert calc_chi_sq_mnf2(obj) == pytest.approx(chi_sq_msg, rel=1e-9)


def test_split_atom_with_multipoles_is_refused():
    """An atom whose orbit splits in the group of k cannot keep its
    multipoles silently: the error says to add them again afterwards."""
    obj = cryspy.load_file(MNF2)
    crystal = obj.crystal_mnf2
    remove_multipoles(crystal)
    crystal.add_multipole_form_factor("F")
    from cryspy.D_functions_item_loop.function_1_magnetic_k import split_atoms
    k_group = {"parent_rotation": numpy.array([[1, 0, 0, 0, 1, 0, 0, 0, 1],
                                               [-1, 0, 0, 0, -1, 0, 0, 0, 1]]).transpose(),
               "parent_translation": numpy.zeros((3, 2)),
               "rotation": numpy.eye(3).reshape(9, 1),
               "translation": numpy.zeros((3, 1)),
               "centring": numpy.zeros((3, 1))}
    with pytest.raises(CrysPyError, match="add_multipole_form_factor"):
        split_atoms(crystal, k_group)
