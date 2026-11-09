"""A structure with a propagation vector gives the same powder pattern as the
same structure written explicitly in the magnetic supercell and calculated by
the independent code for ordered moments (mcif).

k = (1/2, 0, 1/2): supercell (2a, b, 2c), m = S exp(-2 pi i k.R).
k = (1/3, 0, 0): supercell (3a, b, c), m = S exp(-2 pi i k.R) + c.c.
"""
import numpy

import cryspy
from cryspy.procedure_rhochi.rhochi_by_dictionary import (
    rhochi_calc_chi_sq_by_dictionary,
)

CELL = (3.9, 5.1, 6.3)
K_VECTOR = numpy.array([0.5, 0.0, 0.5])
# label, type, position, Fourier component along a, b, c (μB) or None
ATOMS = [
    ("Fe1", "Fe", (0.11, 0.23, 0.37), (1.2, -0.7, 2.1)),
    ("Fe2", "Fe", (0.61, 0.08, 0.19), (-0.4, 1.9, 0.5)),
    ("O1", "O", (0.32, 0.41, 0.05), None),
]
IMAG = {"Fe1": (0.0, 0.0, 0.0), "Fe2": (0.0, 0.0, 0.0)}
PHASE = {"Fe1": 0.0, "Fe2": 0.0}
LATTICE = ((0, 0, 0), (1, 0, 0), (0, 0, 1), (1, 0, 1))
SUPER = (2, 1, 2)
# operation, action on the moment, phase (2π)
OPERATIONS = [("x,y,z", "u,v,w", 0.0), ("-x,-y,-z", "-u,-v,w", 0.0)]
CENTRING = [(0.0, 0.0, 0.0), (0.5, 0.5, 0.5)]
MOMENT_MATRIX = {"u,v,w": numpy.eye(3), "-u,-v,w": numpy.diag([-1., -1., 1.])}
ROTATION = {"x,y,z": numpy.eye(3), "-x,-y,-z": -numpy.eye(3)}

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


def form_rcif_k():
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
    lines += ["", "loop_\n_magnetic_structure_k_number\n"
              "_magnetic_structure_k_vector_1\n_magnetic_structure_k_vector_2\n"
              "_magnetic_structure_k_vector_3",
              "1 {} {} {}".format(*K_VECTOR), "",
              "loop_\n_atom_site_moment_fourier_label\n"
              "_atom_site_moment_fourier_k_number\n"
              "_atom_site_moment_fourier_real_x\n"
              "_atom_site_moment_fourier_real_y\n"
              "_atom_site_moment_fourier_real_z\n"
              "_atom_site_moment_fourier_imag_x\n"
              "_atom_site_moment_fourier_imag_y\n"
              "_atom_site_moment_fourier_imag_z\n"
              "_atom_site_moment_fourier_phase"]
    lines += [f"{lab} 1 {m[0]} {m[1]} {m[2]} {IMAG[lab][0]} {IMAG[lab][1]} "
              f"{IMAG[lab][2]} {PHASE[lab]}"
              for lab, _, _, m in ATOMS if m is not None]
    lines += ["", "loop_\n_space_group_symop_magn_k_id\n"
              "_space_group_symop_magn_k_operation_xyz\n"
              "_space_group_symop_magn_k_moment_uvw\n"
              "_space_group_symop_magn_k_phase"]
    lines += [f"{i + 1} {op} {mom} {ph}"
              for i, (op, mom, ph) in enumerate(OPERATIONS)]
    lines += ["", "loop_\n_space_group_symop_magn_k_centering_id\n"
              "_space_group_symop_magn_k_centering_xyz",
              "1 x+1/2,y+1/2,z+1/2", ""]
    return "\n".join(lines) + EXPERIMENT + "50.0 1.0 1.0\n"


def form_rcif_supercell():
    """Every atom of the supercell written explicitly.

    m = S exp(-2 pi i k.R) for k equivalent to -k, otherwise
    m = S exp(-2 pi i k.R) + c.c. with S = (real + i imag) exp(-2 pi i phase).
    """
    k_half = numpy.allclose(numpy.mod(2 * K_VECTOR, 1.), 0.)
    rows, moments = [], []
    for lab, typ, xyz, fourier in ATOMS:
        n_atom = 0
        for op, mom, phase in OPERATIONS:
            for t_c in CENTRING:
                pos = ROTATION[op] @ numpy.array(xyz) + numpy.array(t_c)
                for lattice in LATTICE:
                    pos_l = pos + numpy.array(lattice)
                    pos_super = numpy.mod(pos_l / numpy.array(SUPER, dtype=float), 1.)
                    n_atom += 1
                    name = f"{lab}_{n_atom}"
                    rows.append(f"{name} {typ} {pos_super[0]:.8f} "
                                f"{pos_super[1]:.8f} {pos_super[2]:.8f} Biso 0.0 1.0")
                    if fourier is not None:
                        k_phase = numpy.dot(K_VECTOR, numpy.array(t_c) + lattice)
                        s_k = (numpy.array(fourier) + 1j * numpy.array(IMAG[lab])) * \
                            numpy.exp(-2j * numpy.pi * PHASE[lab])
                        m_c = MOMENT_MATRIX[mom] @ s_k * numpy.exp(
                            -2j * numpy.pi * (phase + k_phase))
                        m = m_c.real if k_half else 2. * m_c.real
                        moments.append((name, m))
    lines = ["global_\n\ndata_phase\n",
             f"_cell_length_a {SUPER[0] * CELL[0]}",
             f"_cell_length_b {SUPER[1] * CELL[1]}",
             f"_cell_length_c {SUPER[2] * CELL[2]}", "_cell_angle_alpha 90.0",
             "_cell_angle_beta 90.0", "_cell_angle_gamma 90.0", "",
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
    lines += [f"{name} {m[0]:.8f} {m[1]:.8f} {m[2]:.8f}" for name, m in moments]
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


def test_k_structure_equals_supercell():
    ttheta = numpy.linspace(8., 95., 1500)
    signal_k, phase_k = calc_signal(form_rcif_k(), ttheta)
    signal_s, _ = calc_signal(form_rcif_supercell(), ttheta)
    # magnetic satellites are present
    assert phase_k["iint_magnetic_k"].max() > 1.
    # the supercell has 4 times more atoms: |F|^2 is 16 times larger
    numpy.testing.assert_allclose(signal_s, 16. * signal_k,
                                  rtol=1e-6, atol=1e-6 * signal_s.max())


def test_incommensurate_like_k_equals_supercell(monkeypatch):
    """k = (1/3, 0, 0) is not equivalent to -k: complex S, phases, both arms."""
    module = globals()
    monkeypatch.setitem(module, "K_VECTOR", numpy.array([1. / 3., 0., 0.]))
    monkeypatch.setitem(module, "LATTICE", ((0, 0, 0), (1, 0, 0), (2, 0, 0)))
    monkeypatch.setitem(module, "SUPER", (3, 1, 1))
    monkeypatch.setitem(module, "CENTRING", [(0.0, 0.0, 0.0)])
    monkeypatch.setitem(module, "IMAG", {"Fe1": (0.3, 0.9, -0.2),
                                         "Fe2": (0.0, -0.5, 1.1)})
    monkeypatch.setitem(module, "PHASE", {"Fe1": 0.0, "Fe2": 0.17})
    monkeypatch.setitem(module, "OPERATIONS", [("x,y,z", "u,v,w", 0.0),
                                               ("-x,-y,-z", "-u,-v,w", 0.1)])
    rcif_k = form_rcif_k().replace(
        "\nloop_\n_space_group_symop_magn_k_centering_id\n"
        "_space_group_symop_magn_k_centering_xyz\n1 x+1/2,y+1/2,z+1/2\n", "")
    ttheta = numpy.linspace(8., 95., 1500)
    signal_k, phase_k = calc_signal(rcif_k, ttheta)
    signal_s, _ = calc_signal(form_rcif_supercell(), ttheta)
    assert phase_k["iint_magnetic_k"].max() > 1.
    # 3 times more atoms: |F|^2 is 9 times larger
    numpy.testing.assert_allclose(signal_s, 9. * signal_k,
                                  rtol=1e-6, atol=1e-6 * signal_s.max())


def test_magnetic_part_is_not_zero():
    ttheta = numpy.linspace(8., 95., 1500)
    signal_k, _ = calc_signal(form_rcif_k(), ttheta)
    rcif_nonmag = form_rcif_k().replace("1 1.2 -0.7 2.1", "1 0.0 0.0 0.0").replace(
        "1 -0.4 1.9 0.5", "1 0.0 0.0 0.0")
    signal_n, _ = calc_signal(rcif_nonmag, ttheta)
    assert numpy.abs(signal_k - signal_n).max() > 0.01 * signal_k.max()
