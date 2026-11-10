"""A structure with a propagation vector (msCIF superspace group) gives the
same powder pattern as the same structure written explicitly in the
magnetic supercell and calculated by the independent code for ordered
moments (mcif). The moments of the supercell are made by the real-space
msCIF rule of ``magnetic_k_reference`` (not by CrysPy).

k = (1/2, 0, 1/2) with I centring: supercell (2a, b, 2c).
k = (1/3, 0, 0), cos and sin components, tau = 1/10: supercell (3a, b, c).
"""
import numpy

import cryspy
from cryspy.procedure_rhochi.rhochi_by_dictionary import (
    rhochi_calc_chi_sq_by_dictionary,
)
from magnetic_k_reference import generate_moments, take_cells

CELL = (3.9, 5.1, 6.3)
# label, type, position, Fourier components (cos, sin) along a, b, c (μB)
ATOMS = [
    ("Fe1", "Fe", (0.11, 0.23, 0.37), ((1.2, 0.3), (-0.7, 0.9), (2.1, -0.2))),
    ("Fe2", "Fe", (0.61, 0.08, 0.19), ((-0.4, 0.0), (1.9, -0.5), (0.5, 1.1))),
    ("O1", "O", (0.32, 0.41, 0.05), None),
]
CASES = {
    "k_half": dict(
        k_vector=(0.5, 0.0, 0.5), n_cell=(2, 1, 2),
        operations=["x1,x2,x3,x4,+1", "-x1,-x2,-x3,-x4,+1"],
        centrings=["x1,x2,x3,x4,+1", "x1+1/2,x2+1/2,x3+1/2,x4,+1",
                   "x1,x2,x3,x4+1/2,-1", "x1+1/2,x2+1/2,x3+1/2,x4+1/2,-1"],
        sin=False),
    "k_third": dict(
        k_vector=(1. / 3., 0.0, 0.0), n_cell=(3, 1, 1),
        operations=["x1,x2,x3,x4,+1", "-x1,-x2,-x3,-x4+1/10,+1"],
        centrings=["x1,x2,x3,x4,+1", "x1,x2,x3,x4+1/2,-1"],
        sin=True),
}

EXPERIMENT = """
data_pd

_setup_wavelength 2.4
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
_range_2theta_max 100.0

loop_
_phase_label
_phase_scale
_phase_igsize
phase 1.0 0.0

loop_
_pd_background_2theta
_pd_background_intensity
5.0 0.0
100.0 0.0

loop_
_pd_meas_2theta
_pd_meas_intensity
_pd_meas_intensity_sigma
"""

def take_fourier(case, values):
    """(3, 2) cos and sin; the sin part is dropped when not used."""
    values = numpy.array(values, dtype=float)
    if not case["sin"]:
        values[:, 1] = 0.
    return values


def form_rcif_k(case_name: str, scale_moment: float = 1.) -> str:
    case = CASES[case_name]
    lines = ["global_\n\ndata_phase\n"]
    lines += [f"_cell_length_{n} {v}" for n, v in zip("abc", CELL)]
    lines += ["_cell_angle_alpha 90.0", "_cell_angle_beta 90.0",
              "_cell_angle_gamma 90.0", "",
              "loop_\n_atom_site_label\n_atom_site_type_symbol\n"
              "_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n"
              "_atom_site_adp_type\n_atom_site_B_iso_or_equiv\n"
              "_atom_site_occupancy"]
    lines += [f"{lab} {typ} {x} {y} {z} Biso 0.0 1.0"
              for lab, typ, (x, y, z), _ in ATOMS]
    lines += ["", "loop_\n_atom_site_scat_label\n_atom_site_scat_type_symbol\n"
              "_atom_site_scat_lande"]
    lines += [f"{lab} Fe3+ 2.0" for lab, _, _, m in ATOMS if m is not None]
    lines += ["", "loop_\n_cell_wave_vector_seq_id\n_cell_wave_vector_x\n"
              "_cell_wave_vector_y\n_cell_wave_vector_z",
              "1 {:.10f} {:.10f} {:.10f}".format(*case["k_vector"]), "",
              "loop_\n_space_group_symop_magn_ssg_operation.id\n"
              "_space_group_symop_magn_ssg_operation.algebraic"]
    lines += [f"{i + 1} {text}" for i, text in enumerate(case["operations"])]
    lines += ["", "loop_\n_space_group_symop_magn_ssg_centering.id\n"
              "_space_group_symop_magn_ssg_centering.algebraic"]
    lines += [f"{i + 1} {text}" for i, text in enumerate(case["centrings"])]
    lines += ["", "loop_\n_atom_site_moment_Fourier.atom_site_label\n"
              "_atom_site_moment_Fourier.axis\n"
              "_atom_site_moment_Fourier.wave_vector_seq_id\n"
              "_atom_site_moment_Fourier_param.cos\n"
              "_atom_site_moment_Fourier_param.sin"]
    for lab, _, _, fourier in ATOMS:
        if fourier is None:
            continue
        values = scale_moment * take_fourier(case, fourier)
        lines += [f"{lab} {axis} 1 {values[i, 0]} {values[i, 1]}"
                  for i, axis in enumerate("xyz")]
    return "\n".join(lines) + "\n" + EXPERIMENT + "50.0 1.0 1.0\n"


def form_rcif_supercell(case_name: str) -> str:
    """Every atom of the supercell written explicitly with its moment."""
    case = CASES[case_name]
    n_cell = numpy.array(case["n_cell"])
    assert numpy.all(n_cell == take_cells(case["k_vector"]))
    atoms = {lab: (xyz, take_fourier(case, fourier) if fourier is not None
                   else numpy.zeros((3, 2)))
             for lab, _, xyz, fourier in ATOMS}
    sites = generate_moments(case["operations"], case["centrings"],
                             case["k_vector"], atoms, n_cell, atol=1e-8)
    rows, moments = [], []
    for lab, typ, _, fourier in ATOMS:
        for i_site, (position, moment) in enumerate(sites[lab]):
            name = f"{lab}_{i_site + 1}"
            xyz = position / n_cell
            rows.append(f"{name} {typ} {xyz[0]:.10f} {xyz[1]:.10f} {xyz[2]:.10f} "
                        "Biso 0.0 1.0")
            if fourier is not None:
                moments.append((name, moment))
    lines = ["global_\n\ndata_phase\n"]
    lines += [f"_cell_length_{n} {m * v}" for n, m, v in zip("abc", n_cell, CELL)]
    lines += ["_cell_angle_alpha 90.0", "_cell_angle_beta 90.0",
              "_cell_angle_gamma 90.0", "",
              "loop_\n_space_group_symop_magn_operation_xyz\n"
              "_space_group_symop_magn_operation_id", "x,y,z,+1 1", "",
              "loop_\n_space_group_symop_magn_centering_xyz\n"
              "_space_group_symop_magn_centering_id", "x,y,z,+1 1", "",
              "loop_\n_atom_site_label\n_atom_site_type_symbol\n"
              "_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n"
              "_atom_site_adp_type\n_atom_site_B_iso_or_equiv\n"
              "_atom_site_occupancy"] + rows
    lines += ["", "loop_\n_atom_site_scat_label\n_atom_site_scat_type_symbol\n"
              "_atom_site_scat_lande"]
    lines += [f"{name} Fe3+ 2.0" for name, _ in moments]
    lines += ["", "loop_\n_atom_site_moment_label\n"
              "_atom_site_moment_crystalaxis_x\n_atom_site_moment_crystalaxis_y\n"
              "_atom_site_moment_crystalaxis_z"]
    lines += [f"{name} {m[0]:.10f} {m[1]:.10f} {m[2]:.10f}" for name, m in moments]
    return "\n".join(lines) + "\n" + EXPERIMENT + "50.0 1.0 1.0\n"


def calc_signal(rcif: str, ttheta):
    obj = cryspy.str_to_globaln(rcif)
    dict_obj = obj.get_dictionary()
    dict_obj["pd_pd"]["ttheta"] = numpy.radians(ttheta)
    dict_obj["pd_pd"]["signal_exp"] = numpy.ones((2, ttheta.size))
    dict_obj["pd_pd"]["excluded_points"] = numpy.zeros(ttheta.size, dtype=bool)
    dict_in_out = {}
    rhochi_calc_chi_sq_by_dictionary(
        dict_obj, dict_in_out=dict_in_out, flag_use_precalculated_data=False,
        flag_calc_analytical_derivatives=False)
    pd = dict_in_out["pd_pd"]
    return pd["signal_plus"] + pd["signal_minus"], pd["dict_in_out_phase"]


def check_powder(case_name: str):
    ttheta = numpy.linspace(8., 95., 1500)
    signal_k, phase_k = calc_signal(form_rcif_k(case_name), ttheta)
    signal_s, _ = calc_signal(form_rcif_supercell(case_name), ttheta)
    # magnetic satellites are present
    assert phase_k["iint_magnetic_k"].max() > 0.1
    # the supercell has N times more atoms: |F|^2 is N^2 times larger
    n_cell_sq = float(numpy.prod(CASES[case_name]["n_cell"])) ** 2
    numpy.testing.assert_allclose(signal_s, n_cell_sq * signal_k,
                                  rtol=1e-6, atol=1e-6 * signal_s.max())


def test_k_half_equals_supercell():
    check_powder("k_half")


def test_k_third_equals_supercell():
    """k = (1/3, 0, 0) is not equivalent to -k: cos and sin, both arms."""
    check_powder("k_third")


def test_magnetic_part_is_not_zero():
    ttheta = numpy.linspace(8., 95., 1500)
    signal_k, _ = calc_signal(form_rcif_k("k_half"), ttheta)
    signal_n, _ = calc_signal(form_rcif_k("k_half", scale_moment=0.), ttheta)
    assert numpy.abs(signal_k - signal_n).max() > 0.01 * signal_k.max()
