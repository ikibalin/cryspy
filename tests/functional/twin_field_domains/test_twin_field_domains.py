"""Field-induced moments in twin domains (jana_twin, flip ratios).

A twin domain is the crystal rotated by the twin law, so the magnetic field
and the beam polarization must be taken in the Cartesian system of that
domain. The independent reference: the second domain alone (fraction 1) is
an untwinned crystal measured at the indices T h with the orientation
matrix UB T^-1 (the same scattering vectors in the laboratory).
"""
import warnings

import numpy
import pytest

import cryspy
from cryspy.procedure_rhochi.rhochi import rhochi_no_refinement

# tetragonal metric with anisotropic susceptibility (P 1): chi_11 != chi_22,
# so a rotation by 90 degrees about c changes the induced moment
CRYSTAL = """
data_fe
_cell_length_a 4.0
_cell_length_b 4.0
_cell_length_c 6.0
_cell_angle_alpha 90.0
_cell_angle_beta 90.0
_cell_angle_gamma 90.0
_space_group_name_H-M_alt "P 1"
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
_atom_site_adp_type
_atom_site_U_iso_or_equiv
Fe Fe3+ 0.10 0.20 0.30 1.0 Uiso 0.0
O  O2-  0.60 0.30 0.70 1.0 Uiso 0.0
loop_
_atom_site_susceptibility_label
_atom_site_susceptibility_chi_type
_atom_site_susceptibility_chi_11
_atom_site_susceptibility_chi_22
_atom_site_susceptibility_chi_33
_atom_site_susceptibility_chi_12
_atom_site_susceptibility_chi_13
_atom_site_susceptibility_chi_23
Fe Cani 1.0 0.1 0.4 0.2 0.0 0.0
loop_
_atom_site_scat_label
_atom_site_scat_Lande
_atom_site_scat_kappa
Fe 2.0 1.0
"""

# rotation by 90 degrees about c in the Jana convention h_2 = T h_1
TWIN_LAW = numpy.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])

INDEX_HKL = numpy.array([
    [1, 0, 1], [0, 1, 1], [1, 1, 0], [2, 1, 1], [1, 2, 1], [1, -1, 2],
    [2, 0, 1], [0, 2, 1], [1, 0, 3]], dtype=float).transpose()


def form_experiment(name: str, matrix_ub, index_hkl, twin_fraction=None,
                    twin_law=TWIN_LAW, time_reversal=None) -> str:
    ls_out = [
        f"data_{name:}", "_setup_wavelength 0.84", "_setup_field 5.0",
        "_setup_radiation neutrons", "_setup_K 0.0", "_setup_cthm 0.91",
        "_diffrn_radiation_polarization 0.95", "_diffrn_radiation_efficiency 1.0"]
    for ind_i in range(3):
        for ind_j in range(3):
            ls_out.append(f"_diffrn_orient_matrix_ub_{ind_i+1:}{ind_j+1:} {matrix_ub[ind_i, ind_j]:.10f}")
    ls_out.extend(["_diffrn_orient_matrix_type CCSL", "_phase_label fe",
                   "loop_", "_diffrn_refln_index_h", "_diffrn_refln_index_k",
                   "_diffrn_refln_index_l", "_diffrn_refln_fr", "_diffrn_refln_fr_sigma"])
    for hkl in index_hkl.transpose():
        ls_out.append(f"{int(round(hkl[0])):} {int(round(hkl[1])):} {int(round(hkl[2])):} 1.0 0.01")
    if twin_fraction is not None:
        ls_out.extend(["loop_", "_jana_twin_id"] + [
            f"_jana_twin_matrix[{ind_i+1:}][{ind_j+1:}]" for ind_i in range(3) for ind_j in range(3)] +
            ["_jana_twin_fraction"] +
            ([] if time_reversal is None else ["_jana_twin_time_reversal"]) +
            ["2 " + " ".join(f"{val:.1f}" for val in twin_law.flatten()) + f" {twin_fraction:}" +
             ("" if time_reversal is None else f" {time_reversal:}")])
    return "\n".join(ls_out)


def calc_flip_ratios(text_experiment: str, crystal: str = CRYSTAL) -> numpy.ndarray:
    global_object = cryspy.str_to_globaln(crystal + "\n" + text_experiment)
    rhochi_no_refinement(global_object)
    diffrn = [item for item in global_object.items if type(item).__name__ == "Diffrn"][0]
    return numpy.array(diffrn.diffrn_refln.fr_calc, dtype=float)


def matrix_ub_first_domain() -> numpy.ndarray:
    # a general orientation: the field is not along any axis of the crystal
    angle_1, angle_2 = 0.4, 0.9
    rot_1 = numpy.array([[numpy.cos(angle_1), -numpy.sin(angle_1), 0.],
                         [numpy.sin(angle_1), numpy.cos(angle_1), 0.], [0., 0., 1.]])
    rot_2 = numpy.array([[1., 0., 0.], [0., numpy.cos(angle_2), -numpy.sin(angle_2)],
                         [0., numpy.sin(angle_2), numpy.cos(angle_2)]])
    matrix_b = numpy.diag([1. / 4., 1. / 4., 1. / 6.])
    return rot_2.dot(rot_1).dot(matrix_b)


def test_second_domain_is_rotated_crystal():
    matrix_ub = matrix_ub_first_domain()
    fr_twin = calc_flip_ratios(form_experiment("twin", matrix_ub, INDEX_HKL, twin_fraction=1.0))
    fr_rotated = calc_flip_ratios(form_experiment(
        "rotated", matrix_ub.dot(numpy.linalg.inv(TWIN_LAW)), TWIN_LAW.dot(INDEX_HKL)))
    # the flip ratios differ from 1 by ~0.1-1; the model is exact
    assert numpy.allclose(fr_twin, fr_rotated, rtol=1e-8, atol=1e-10)
    # the check is not trivial: the second domain differs from the first one
    fr_first = calc_flip_ratios(form_experiment("first", matrix_ub, INDEX_HKL))
    assert not numpy.allclose(fr_twin, fr_first, rtol=1e-3)


def test_first_domain_only():
    matrix_ub = matrix_ub_first_domain()
    fr_twin = calc_flip_ratios(form_experiment("twin", matrix_ub, INDEX_HKL, twin_fraction=0.0))
    fr_first = calc_flip_ratios(form_experiment("first", matrix_ub, INDEX_HKL))
    assert numpy.allclose(fr_twin, fr_first, rtol=1e-10, atol=1e-12)


MEM_PARAMETERS = """
global_
_mem_parameters_points_a 12
_mem_parameters_points_b 12
_mem_parameters_points_c 18
_mem_parameters_channel_col False
_mem_parameters_channel_chi True
_mem_parameters_magnetization_plus 4.0
_mem_parameters_magnetization_minus -1.0
_mem_parameters_prior_density uniform
_mem_parameters_use_asymmetry False
_mem_parameters_GoF_desired 1.00
"""


def calc_mem_chi_sq(text_experiment: str) -> float:
    global_object = cryspy.str_to_globaln(MEM_PARAMETERS + CRYSTAL + "\n" + text_experiment)
    rhochi_no_refinement(global_object)
    diffrn = [item for item in global_object.items if type(item).__name__ == "Diffrn"][0]
    # "measured" flip ratios: the model ones, slightly changed
    diffrn.diffrn_refln.fr = [1.02 * float(value) for value in diffrn.diffrn_refln.fr_calc]
    return cryspy.mempy_magnetization_density_reconstruction(global_object)["chi_sq"]


def test_mem_second_domain_is_rotated_crystal():
    matrix_ub = matrix_ub_first_domain()
    chi_sq_twin = calc_mem_chi_sq(form_experiment("twin", matrix_ub, INDEX_HKL, twin_fraction=1.0))
    chi_sq_rotated = calc_mem_chi_sq(form_experiment(
        "rotated", matrix_ub.dot(numpy.linalg.inv(TWIN_LAW)), TWIN_LAW.dot(INDEX_HKL)))
    # the same MEM problem: equal to the convergence tolerance of MEM
    assert chi_sq_twin == pytest.approx(chi_sq_rotated, rel=1e-6)


def form_crystal_ordered(moment) -> str:
    """The same atoms with an ordered moment (in μB along a, b, c) of Fe."""
    text = CRYSTAL.split("loop_\n_atom_site_susceptibility_label")[0]
    return text + f"""
loop_
_space_group_symop_magn_operation_xyz
_space_group_symop_magn_operation_id
x,y,z,+1 1
loop_
_space_group_symop_magn_centering_xyz
_space_group_symop_magn_centering_id
x,y,z,+1 1
loop_
_atom_site_moment_label
_atom_site_moment_crystalaxis_x
_atom_site_moment_crystalaxis_y
_atom_site_moment_crystalaxis_z
_atom_site_moment_symmform
Fe {moment[0]:} {moment[1]:} {moment[2]:} mx,my,mz
loop_
_atom_site_scat_label
_atom_site_scat_Lande
_atom_site_scat_kappa
Fe 2.0 1.0
"""


@pytest.mark.filterwarnings("ignore:The net ordered moment")
def test_ordered_moments_magnetic_twin_operation():
    # the domain is the rotated crystal with the moments theta m
    matrix_ub = matrix_ub_first_domain()
    moment = (0.5, 1.5, -1.0)
    crystal = form_crystal_ordered(moment)
    text_rotated = form_experiment(
        "rotated", matrix_ub.dot(numpy.linalg.inv(TWIN_LAW)), TWIN_LAW.dot(INDEX_HKL))
    for time_reversal, moment_rotated in ((1, moment), (-1, tuple(-val for val in moment))):
        fr_twin = calc_flip_ratios(form_experiment(
            "twin", matrix_ub, INDEX_HKL, twin_fraction=1.0, time_reversal=time_reversal), crystal)
        fr_rotated = calc_flip_ratios(text_rotated, form_crystal_ordered(moment_rotated))
        assert numpy.allclose(fr_twin, fr_rotated, rtol=1e-8, atol=1e-10)


def test_warning_net_moment_against_field():
    # ferromagnet along c, field along c (U = 1); a twofold axis along a
    # turns the field of the second domain to -c
    matrix_ub = numpy.diag([1. / 4., 1. / 4., 1. / 6.])
    twin_law = numpy.diag([1., -1., -1.])
    crystal = form_crystal_ordered((0., 0., 2.))
    with pytest.warns(cryspy.CrysPyWarning, match="antiparallel"):
        calc_flip_ratios(form_experiment(
            "twin", matrix_ub, INDEX_HKL, twin_fraction=0.5, twin_law=twin_law), crystal)
    with warnings.catch_warnings():
        warnings.filterwarnings("error", message=".*antiparallel.*")
        calc_flip_ratios(form_experiment(
            "twin", matrix_ub, INDEX_HKL, twin_fraction=0.5, twin_law=twin_law,
            time_reversal=-1), crystal)


def test_assistant_tool_check_twin_time_reversal():
    # "measured" flip ratios of an antiferromagnet-like case decided by the
    # data: the second domain is time-reversed (theta = -1)
    matrix_ub = matrix_ub_first_domain()
    crystal = form_crystal_ordered((0.5, 1.5, -1.0))
    fr_measured = calc_flip_ratios(form_experiment(
        "twin", matrix_ub, INDEX_HKL, twin_fraction=0.4, time_reversal=-1), crystal)
    global_object = cryspy.str_to_globaln(crystal + "\n" + form_experiment(
        "twin", matrix_ub, INDEX_HKL, twin_fraction=0.4))
    diffrn = [item for item in global_object.items if type(item).__name__ == "Diffrn"][0]
    diffrn.diffrn_refln.fr = [float(value) for value in fr_measured]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = cryspy.call_assistant_tool(global_object, "check_twin_time_reversal")
    assert "-1 (was +1)" in result["text"]
    assert diffrn.jana_twin.items[0].time_reversal == -1
    assert result["metrics"]["chi_sq"] == pytest.approx(0., abs=1e-12)
