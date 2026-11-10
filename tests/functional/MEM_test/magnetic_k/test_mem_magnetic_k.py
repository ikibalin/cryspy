"""MEM of the magnetization density of structures with a propagation vector.

The density of every magnetic atom is reconstructed in the cell (the
shape of its moment); its moment in the cell L is given by the Fourier
amplitude, Re(A exp(-2 pi i k.L)). For a commensurate k the same
structure written in its magnetic supercell with the magnetic space group
(Crystal.make_magnetic_supercell(use_symmetry=True), reflections N H,
scale / det(N)^2) is calculated by the independent code for ordered
moments and must give the same MEM. An incommensurate k has no
supercell, MEM works all the same.
"""
import itertools

import numpy
import pytest

import cryspy
from cryspy.D_functions_item_loop.function_1_magnetic_k import take_supercell_multipliers
from cryspy.procedure_mempy.mempy_by_dictionary import (
    calc_preliminary_information_by_dictionary, mempy_reconstruction_by_dictionary)
from cryspy.procedure_rhochi.rhochi_by_dictionary import (
    rhochi_calc_chi_sq_by_dictionary)

CELL = (3.9, 5.1, 6.3)
CASES = {
    # I centring, time reversal with x4 + 1/2: supercell (2a, b, 2c)
    "k_half": dict(
        k_vector=(0.5, 0.0, 0.5),
        operations=["x1,x2,x3,x4,+1", "-x1,-x2,-x3,-x4,+1"],
        centrings=["x1,x2,x3,x4,+1", "x1+1/2,x2+1/2,x3+1/2,x4,+1",
                   "x1,x2,x3,x4+1/2,-1", "x1+1/2,x2+1/2,x3+1/2,x4+1/2,-1"],
        fourier=((1.2, 0.0), (-0.7, 0.0), (2.1, 0.0)), centring_i=True),
    # cos and sin components, phase 1/10: supercell (3a, b, c), group P1
    "k_third": dict(
        k_vector=(1. / 3., 0.0, 0.0),
        operations=["x1,x2,x3,x4,+1", "-x1,-x2,-x3,-x4+1/10,+1"],
        centrings=["x1,x2,x3,x4,+1", "x1,x2,x3,x4+1/2,-1"],
        fourier=((1.2, 0.3), (-0.7, 0.9), (2.1, -0.2)), centring_i=False),
    "incommensurate": dict(
        k_vector=(0.2817, 0.0, 0.0),
        operations=["x1,x2,x3,x4,+1", "-x1,-x2,-x3,-x4+1/10,+1"],
        centrings=["x1,x2,x3,x4,+1", "x1,x2,x3,x4+1/2,-1"],
        fourier=((1.2, 0.3), (-0.7, 0.9), (2.1, -0.2)), centring_i=False),
}

MEM = """global_

_mem_parameters_points_a 16
_mem_parameters_points_b 16
_mem_parameters_points_c 16
_mem_parameters_channel_col False
_mem_parameters_channel_ani True
_mem_parameters_magnetization_plus 1.0
_mem_parameters_magnetization_minus -4.0
_mem_parameters_prior_density core
_mem_parameters_use_asymmetry False
_mem_parameters_GoF_desired 1.0
"""

EXPERIMENT = """
data_exp1

_setup_wavelength 0.9
_setup_field 0.0
_setup_radiation neutrons
_setup_K 0.0
_setup_cthm 0.91

_diffrn_orient_matrix_ub_11 {ub_11}
_diffrn_orient_matrix_ub_12 0.0
_diffrn_orient_matrix_ub_13 0.0
_diffrn_orient_matrix_ub_21 0.0
_diffrn_orient_matrix_ub_22 {ub_22}
_diffrn_orient_matrix_ub_23 0.0
_diffrn_orient_matrix_ub_31 0.0
_diffrn_orient_matrix_ub_32 0.0
_diffrn_orient_matrix_ub_33 {ub_33}
_diffrn_orient_matrix_type CCSL

_extinction_model gauss
_extinction_mosaicity 0.0
_extinction_radius 2.0

_phase_label phase
_phase_scale {scale}

loop_
_diffrn_refln_index_h
_diffrn_refln_index_k
_diffrn_refln_index_l
_diffrn_refln_intensity
_diffrn_refln_intensity_sigma
"""


def form_crystal(case_name: str) -> str:
    case = CASES[case_name]
    lines = ["data_phase", ""]
    lines += [f"_cell_length_{n} {v}" for n, v in zip("abc", CELL)]
    lines += ["_cell_angle_alpha 90.0", "_cell_angle_beta 90.0", "_cell_angle_gamma 90.0", "",
              "loop_", "_atom_site_label", "_atom_site_type_symbol", "_atom_site_fract_x",
              "_atom_site_fract_y", "_atom_site_fract_z", "_atom_site_adp_type",
              "_atom_site_B_iso_or_equiv", "_atom_site_occupancy",
              "Fe1 Fe 0.11 0.23 0.37 Biso 0.3 1.0", "O1 O 0.32 0.41 0.05 Biso 0.5 1.0", "",
              "loop_", "_atom_electron_configuration_label",
              "_atom_electron_configuration_core", "Fe1 3d5", "",
              "loop_", "_atom_site_scat_label", "_atom_site_scat_type_symbol",
              "_atom_site_scat_lande", "Fe1 Fe3+ 2.0", "",
              "loop_", "_cell_wave_vector_seq_id", "_cell_wave_vector_x",
              "_cell_wave_vector_y", "_cell_wave_vector_z",
              "1 {:.10f} {:.10f} {:.10f}".format(*case["k_vector"]), "",
              "loop_", "_space_group_symop_magn_ssg_operation.id",
              "_space_group_symop_magn_ssg_operation.algebraic"]
    lines += [f"{i + 1} {text}" for i, text in enumerate(case["operations"])]
    lines += ["", "loop_", "_space_group_symop_magn_ssg_centering.id",
              "_space_group_symop_magn_ssg_centering.algebraic"]
    lines += [f"{i + 1} {text}" for i, text in enumerate(case["centrings"])]
    lines += ["", "loop_", "_atom_site_moment_Fourier.atom_site_label",
              "_atom_site_moment_Fourier.axis", "_atom_site_moment_Fourier.wave_vector_seq_id",
              "_atom_site_moment_Fourier_param.cos", "_atom_site_moment_Fourier_param.sin"]
    lines += [f"Fe1 {axis} 1 {cos} {sin}" for axis, (cos, sin) in zip("xyz", case["fourier"])]
    return "\n".join(lines) + "\n"


def form_reflections(case_name: str):
    """Nuclear reflections and satellites h ± k (fractional H)."""
    case = CASES[case_name]
    k_vector = numpy.array(case["k_vector"])
    l_hkl = []
    for hkl in itertools.product(range(-2, 3), range(0, 3), range(0, 3)):
        if not case["centring_i"] or sum(hkl) % 2 == 0:
            l_hkl.append(hkl)
            for sign in (1, -1):
                l_hkl.append(tuple(numpy.array(hkl) + sign * k_vector))
    hkl = numpy.unique(numpy.round(numpy.array(l_hkl), 10), axis=0)
    return hkl[numpy.any(hkl != 0., axis=1)]


def form_experiment(hkl, intensity, cell, scale):
    text = EXPERIMENT.format(ub_11=1. / cell[0], ub_22=1. / cell[1],
                             ub_33=1. / cell[2], scale=scale)
    rows = [f"{h:.10g} {k:.10g} {l:.10g} {value:.6f} {0.02 * value + 0.01:.6f}"
            for (h, k, l), value in zip(hkl, intensity)]
    return text + "\n".join(rows) + "\n"


def calc_chi_sq(rcif: str):
    obj = cryspy.str_to_globaln(rcif)
    dict_in_out = {}
    chi_sq = rhochi_calc_chi_sq_by_dictionary(
        obj.get_dictionary(), dict_in_out=dict_in_out,
        flag_use_precalculated_data=False, flag_calc_analytical_derivatives=False)[0]
    return chi_sq, dict_in_out


def form_rcif(case_name: str) -> str:
    """Intensities calculated with k, disturbed by 3 % (fixed seed)."""
    hkl = form_reflections(case_name)
    crystal = form_crystal(case_name)
    rcif = MEM + crystal + form_experiment(hkl, numpy.ones(len(hkl)), CELL, 10.)
    _, dict_in_out = calc_chi_sq(rcif)
    intensity = dict_in_out["diffrn_exp1"]["iint_plus"] + dict_in_out["diffrn_exp1"]["iint_minus"]
    intensity = 10. * numpy.asarray(intensity).flatten()
    random = numpy.random.default_rng(7)
    intensity = intensity * (1. + 0.03 * random.standard_normal(intensity.size))
    return MEM + crystal + form_experiment(hkl, intensity, CELL, 10.)


def take_crystal(obj):
    return [item for item in obj.items if isinstance(item, cryspy.Crystal)][0]


def take_supercell_dictionaries(obj):
    """MEM dictionaries of the structure written in its magnetic supercell."""
    dict_obj = obj.get_dictionary()
    n_cell = take_supercell_multipliers(dict_obj["crystal_phase"]["magnetic_k_vector"])
    supercell = take_crystal(obj).make_magnetic_supercell(data_name="phase", use_symmetry=True)
    dict_diffrn = dict(dict_obj["diffrn_exp1"])
    index_hkl = n_cell[:, numpy.newaxis] * dict_diffrn["index_hkl"]
    assert numpy.allclose(index_hkl, numpy.round(index_hkl), atol=1e-6)
    dict_diffrn["index_hkl"] = numpy.round(index_hkl)
    # F and the cell volume are det(N) times larger
    dict_diffrn["phase_scale"] = dict_diffrn["phase_scale"] / float(numpy.prod(n_cell)) ** 2
    dict_mem = dict(dict_obj["mem_parameters"])
    dict_mem["points_abc"] = numpy.array(dict_mem["points_abc"]) * n_cell
    return supercell, dict_mem, dict_diffrn, n_cell


def test_supercell_magnetic_space_group():
    crystal = take_crystal(cryspy.str_to_globaln(MEM + form_crystal("k_half")))
    supercell = crystal.make_magnetic_supercell(use_symmetry=True)
    elems = supercell.get_dictionary()["full_mcif_elems"]
    # -1 and I centring of the parent cell with the time reversal of the
    # half translations: 2 x 2 x 4 = 16 operations, one copy of every atom
    assert elems.shape[1] == 16
    assert list(supercell.atom_site.label) == ["Fe1", "O1"]
    assert sorted(set(elems[-1])) == [-1, 1]
    # phase 1/10 of k = 1/3: the section t = 0 keeps no rotation (P1)
    crystal = take_crystal(cryspy.str_to_globaln(MEM + form_crystal("k_third")))
    supercell = crystal.make_magnetic_supercell(use_symmetry=True)
    assert supercell.get_dictionary()["full_mcif_elems"].shape[1] == 1
    assert len(supercell.atom_site.label) == 12  # 3 cells x 2 x (Fe1, O1)


@pytest.mark.parametrize("case_name", ["k_half", "k_third"])
def test_supercell_same_chi_sq_as_k(case_name):
    rcif = form_rcif(case_name)
    chi_sq_k, _ = calc_chi_sq(rcif)
    obj = cryspy.str_to_globaln(rcif)
    supercell, _, dict_diffrn, n_cell = take_supercell_dictionaries(obj)
    hkl = dict_diffrn["index_hkl"].transpose()
    intensity = numpy.array(obj.diffrn_exp1.diffrn_refln.intensity)
    rcif_s = MEM + supercell.to_cif() + form_experiment(
        hkl, intensity, numpy.array(CELL) * n_cell, 10. / float(numpy.prod(n_cell)) ** 2)
    chi_sq_s, _ = calc_chi_sq(rcif_s)
    # the CIF text of the supercell keeps 5-6 digits of positions and moments
    assert chi_sq_s == pytest.approx(chi_sq_k, rel=1e-4)


def calc_prior_f_m_perp(dict_crystal, dict_mem, dict_diffrn):
    """F_M_perp of the prior density (spherical, 3d5) by the MEM kernel."""
    dict_in_out = {}
    calc_preliminary_information_by_dictionary(dict_crystal, dict_mem, [dict_diffrn], dict_in_out)
    return numpy.tensordot(
        dict_in_out["dict_in_out_diffrn_exp1"]["mem_ani"],
        dict_in_out["density_prior_channel_ani"], axes=([2, 3], [0, 1]))


@pytest.mark.parametrize("case_name", ["k_half", "k_third"])
def test_mem_kernel_k_equals_supercell(case_name):
    """The same density gives the same F_M_perp in the cell with k and in
    the supercell with the ordered moments (F is det(N) times larger)."""
    obj = cryspy.str_to_globaln(form_rcif(case_name))
    dict_obj = obj.get_dictionary()
    f_k = calc_prior_f_m_perp(
        dict_obj["crystal_phase"], dict_obj["mem_parameters"], dict_obj["diffrn_exp1"])
    supercell, dict_mem, dict_diffrn, n_cell = take_supercell_dictionaries(obj)
    f_s = calc_prior_f_m_perp(supercell.get_dictionary(), dict_mem, dict_diffrn)
    assert numpy.abs(f_k).max() > 0.05
    # In P1 (k_third) every atom of the supercell has its own channel: the
    # tails of its spherical prior near the other atoms differ slightly from
    # the channel of the cell (the nearest copy), < 1 % of F.
    numpy.testing.assert_allclose(f_s / float(numpy.prod(n_cell)), f_k, rtol=1e-2, atol=2e-3)


@pytest.mark.parametrize("case_name", ["k_half", "k_third"])
def test_mem_k_and_supercell_converge(case_name):
    rcif = form_rcif(case_name)
    res_k = cryspy.mempy_magnetization_density_reconstruction(cryspy.str_to_globaln(rcif))
    obj = cryspy.str_to_globaln(rcif)
    supercell, dict_mem, dict_diffrn, _ = take_supercell_dictionaries(obj)
    res_s = mempy_reconstruction_by_dictionary(
        supercell.get_dictionary(), dict_mem, [dict_diffrn], {})
    assert res_k["chi_sq"] < 1.
    assert res_s["chi_sq"] < 1.


def test_mem_incommensurate():
    obj = cryspy.str_to_globaln(form_rcif("incommensurate"))
    res = cryspy.mempy_magnetization_density_reconstruction(obj)
    assert res["chi_sq"] < 1.
    intensity_calc = numpy.array(obj.diffrn_exp1.diffrn_refln.intensity_calc)
    flag_satellite = numpy.abs(numpy.array(obj.diffrn_exp1.diffrn_refln.index_h) % 1) > 1e-4
    assert intensity_calc[flag_satellite].max() > 0.05
