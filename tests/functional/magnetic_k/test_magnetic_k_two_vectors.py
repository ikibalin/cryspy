"""Magnetic structures with two propagation vectors (3+2 superspace).

The structure factor of the satellites h ± k1 and h ± k2 and the moments in
real space are compared with the direct calculation of
``magnetic_k_reference`` (independent of CrysPy), which generates every
atom of the supercell by the msCIF rule for two waves.
"""
import warnings

import numpy
import pytest

import cryspy
from cryspy.A_functions_base.structure_factor_magnetic_k import (
    calc_f_m_perp_k_by_dictionary,
    calc_moment_modulation_by_dictionary,
    calc_moments_at_positions,
)
from magnetic_k_reference import calc_f_m_perp, generate_moments_multi_k

CELL = (4.1, 5.3, 6.2, 90., 97., 90.)
K_VECTORS = ((1. / 3., 0., 0.), (0., 0.25, 0.))
N_CELL = numpy.array([3, 4, 1])
OPERATIONS = ["x1,x2,x3,x4,x5,+1", "-x1,x2,-x3,-x4,x5+1/2,+1",
              "-x1,-x2,-x3,-x4,-x5,+1", "x1,-x2,x3,x4,-x5-1/2,+1"]
CENTRINGS = ["x1,x2,x3,x4,x5,+1", "x1,x2,x3,x4+1/2,x5+1/2,-1"]
# label, position, Fourier components per wave (k1, k2): cos, sin along a, b, c (μB)
ATOMS = [
    ("Fe1", (0.13, 0.27, 0.39), [((1.2, 0.3), (-0.7, 0.9), (2.1, -0.2)),
                                 ((0.4, -1.1), (0.8, 0.2), (-0.3, 0.6))]),
    # site 2/m: no sin components (inversion, M(-x) = M(x)); the two-fold
    # axis keeps M along b for k1 and reverses it (x5 + 1/2) for k2
    ("Fe2", (0., 0., 0.), [((0., 0.), (1.3, 0.), (0., 0.)),
                           ((-1.4, 0.), (0., 0.), (0.7, 0.))]),
    ("Fe3", (0.61, 0.08, 0.19), [None, ((0.9, 0.4), (-0.2, 1.5), (0.3, -0.8))]),
]


def form_crystal():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return cryspy.str_to_globaln(form_crystal_text()).items[0]


def form_crystal_text() -> str:
    lines = ["global_", "data_k2", "",
             f"_cell_length_a {CELL[0]}", f"_cell_length_b {CELL[1]}",
             f"_cell_length_c {CELL[2]}", f"_cell_angle_alpha {CELL[3]}",
             f"_cell_angle_beta {CELL[4]}", f"_cell_angle_gamma {CELL[5]}", "",
             "loop_", "_cell_wave_vector_seq_id", "_cell_wave_vector_x",
             "_cell_wave_vector_y", "_cell_wave_vector_z"]
    lines += [f"{i + 1} {k[0]:.10f} {k[1]:.10f} {k[2]:.10f}" for i, k in enumerate(K_VECTORS)]
    lines += ["", "loop_", "_space_group_symop_magn_ssg_operation_id",
              "_space_group_symop_magn_ssg_operation_algebraic"]
    lines += [f"{i + 1} {text}" for i, text in enumerate(OPERATIONS)]
    lines += ["", "loop_", "_space_group_symop_magn_ssg_centering_id",
              "_space_group_symop_magn_ssg_centering_algebraic"]
    lines += [f"{i + 1} {text}" for i, text in enumerate(CENTRINGS)]
    lines += ["", "loop_", "_atom_site_Fourier_wave_vector_seq_id",
              "_atom_site_Fourier_wave_vector_q1_coeff",
              "_atom_site_Fourier_wave_vector_q2_coeff", "1 1 0", "2 0 1"]
    lines += ["", "loop_", "_atom_site_label", "_atom_site_type_symbol",
              "_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z",
              "_atom_site_occupancy", "_atom_site_B_iso_or_equiv"]
    lines += [f"{label} Fe {x} {y} {z} 1 0" for label, (x, y, z), _ in ATOMS]
    lines += ["O1 O 0.32 0.41 0.05 1 0"]
    lines += ["", "loop_", "_atom_site_scat_label", "_atom_site_scat_type_symbol",
              "_atom_site_scat_lande"]
    lines += [f"{label} Fe3+ 2.0" for label, _, _ in ATOMS]
    lines += ["", "loop_", "_atom_site_moment_Fourier_atom_site_label",
              "_atom_site_moment_Fourier_axis", "_atom_site_moment_Fourier_wave_vector_seq_id",
              "_atom_site_moment_Fourier_param.cos", "_atom_site_moment_Fourier_param.sin"]
    for label, _, waves in ATOMS:
        for i_wave, wave in enumerate(waves):
            if wave is None:
                continue
            for axis, (cos, sin) in zip("xyz", wave):
                lines.append(f"{label} {axis} {i_wave + 1} {cos} {sin}")
    return "\n".join(lines)


def take_reference_sites():
    atoms = {}
    for label, position, waves in ATOMS:
        atoms[label] = (position, numpy.array(
            [numpy.zeros((3, 2)) if wave is None else wave for wave in waves], dtype=float))
    return generate_moments_multi_k(OPERATIONS, CENTRINGS, K_VECTORS, atoms, N_CELL,
                                    atol=1e-8)


def test_structure_factor_of_two_propagation_vectors():
    crystal = form_crystal()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        dict_crystal = crystal.get_dictionary()
    assert dict_crystal["magnetic_k_vector"].shape == (3, 2)
    # form factor 1 and B = 0, as in the reference
    dict_crystal["mag_atom_j0_parameters"][:] = 0.
    dict_crystal["mag_atom_j0_parameters"][6] = 1.
    dict_crystal["mag_atom_j2_parameters"][:] = 0.
    dict_crystal["atom_b_iso"][:] = 0.
    rng = numpy.random.default_rng(3)
    index_h = rng.integers(-3, 4, size=(3, 8)).astype(float)
    index_hkl = numpy.concatenate(
        [index_h + sign * numpy.array(k)[:, numpy.newaxis]
         for k in K_VECTORS for sign in (1., -1.)], axis=1)
    f_m_perp = calc_f_m_perp_k_by_dictionary(
        dict_crystal, index_hkl, dict_crystal["unit_cell_parameters"])
    f_reference = calc_f_m_perp(take_reference_sites(), N_CELL, CELL[:3], CELL[3:], index_hkl)
    assert numpy.abs(f_reference).max() > 0.05
    numpy.testing.assert_allclose(f_m_perp, f_reference, atol=1e-8)


def test_moments_of_two_propagation_vectors():
    crystal = form_crystal()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        modulation = calc_moment_modulation_by_dictionary(crystal.get_dictionary())
    n_site = 0
    for l_site in take_reference_sites().values():
        position = numpy.stack([site[0] for site in l_site], axis=1)
        moment = numpy.stack([site[1] for site in l_site], axis=1)
        moment_cryspy = calc_moments_at_positions(
            position, modulation["fract_xyz"], modulation["amplitude"], modulation["k_vector"])
        numpy.testing.assert_allclose(moment_cryspy, moment, atol=1e-8)
        n_site += position.shape[1]
    assert n_site == 12 * (4 + 1 + 4)  # Fe1, Fe3 general, Fe2 on 2/m


def test_supercell_and_report_of_two_propagation_vectors():
    crystal = form_crystal()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        supercell = crystal.make_magnetic_supercell()
        report = crystal.report_magnetic_moments()
    assert "k1" in report and "k2" in report
    cell = supercell.cell
    numpy.testing.assert_allclose(
        [cell.length_a, cell.length_b, cell.length_c], numpy.array(CELL[:3]) * N_CELL)
    sites = take_reference_sites()
    moment_reference = {tuple(numpy.round(numpy.mod(position / N_CELL, 1.), 4)): moment
                        for l_site in sites.values() for position, moment in l_site}
    n_found = 0
    for item in supercell.atom_site_moment.items:
        atom = supercell.atom_site[item.label]
        key = tuple(numpy.round(numpy.mod([atom.fract_x, atom.fract_y, atom.fract_z], 1.), 4))
        moment = numpy.array([item.crystalaxis_x, item.crystalaxis_y, item.crystalaxis_z])
        numpy.testing.assert_allclose(moment, moment_reference[key], atol=2e-3)
        n_found += 1
    assert n_found == len(moment_reference)


def test_operation_exchanging_propagation_vectors_is_refused():
    from cryspy.A_functions_base.structure_factor_magnetic_k import parse_ssg_operation
    with pytest.raises(ValueError, match="exchange"):
        parse_ssg_operation("x2,x1,x3,x5,x4,+1")


def take_waves(crystal, n_wave, rng=None):
    """{label: (position, Fourier (n_wave, 3, 2))} of the crystal; with rng,
    random values of the free symmform symbols (per atom and wave)."""
    from cryspy.A_functions_base.magnetic_k_symmetry import parse_symmform
    d_symbol, atoms = {}, {}
    for atom in crystal.atom_site.items:
        fourier = numpy.zeros((n_wave, 3, 2))
        for item in crystal.atom_site_moment_fourier.items:
            if item.atom_site_label != atom.label:
                continue
            i_wave = int(item.wave_vector_seq_id) - 1
            for i_part, part in enumerate(("cos", "sin")):
                value = getattr(item, part)
                if rng is not None:
                    value = 0.
                    for symbol, coefficient in parse_symmform(
                            getattr(item, f"{part}_symmform")).items():
                        key = (atom.label, i_wave, symbol)
                        if key not in d_symbol:
                            d_symbol[key] = rng.uniform(-2., 2.)
                        value += coefficient * d_symbol[key]
                    setattr(item, part, value)
                fourier[i_wave, "xyz".index(item.axis), i_part] = value
        atoms[atom.label] = ((atom.fract_x, atom.fract_y, atom.fract_z), fourier)
    return atoms


def take_operations(crystal):
    return ([item.algebraic for item in crystal.space_group_symop_magn_ssg_operation.items],
            [item.algebraic for item in crystal.space_group_symop_magn_ssg_centering.items])


@pytest.mark.parametrize("irrep, mode", [("1", 1), ("2", 1), ("2", 4)])
def test_add_second_propagation_vector(irrep, mode):
    """Pnma, k1 = (1/3, 0, 0) then k2 = (0, 0, 1/2) (two-dimensional
    representations): the moments of k1 stay as they were, and the (3+2)
    group with random values of its free Fourier components gives one
    moment per site in the direct generation."""
    from test_magnetic_k_crystal_tools import form_crystal
    k_1, k_2 = (1. / 3., 0., 0.), (0., 0., 0.5)
    n_cell = numpy.array([3, 1, 2])
    rng = numpy.random.default_rng(5)
    crystal = form_crystal("P n m a", (5., 6., 7., 90, 90, 90),
                           [("Mn1", "Mn", 0.1, 0.25, 0.3), ("Mn2", "Mn", 0.0, 0.0, 0.0),
                            ("O1", "O", 0.2, 0.1, 0.4)], magnetic="Mn1")
    crystal.atom_site_scat.items.append(
        type(crystal.atom_site_scat.items[0])(label="Mn2", type_symbol="Mn3+", lande=2.))
    crystal.add_propagation_vector(*k_1, "Mn1,Mn2", irrep="1", moment=2.)
    operations, centrings = take_operations(crystal)
    atoms = take_waves(crystal, 1, rng)
    sites_1 = generate_moments_multi_k(operations, centrings, [k_1], atoms, n_cell, atol=1e-8)
    moment_1 = {tuple(numpy.round(position, 4)): moment
                for l_site in sites_1.values() for position, moment in l_site}

    report = crystal.add_propagation_vector(*k_2, "Mn1", irrep=irrep, mode=mode, moment=0.)
    assert "(3+2)" in report
    operations, centrings = take_operations(crystal)
    assert all(len(text.split(",")) == 6 for text in operations + centrings)
    sites = generate_moments_multi_k(operations, centrings, [k_1, k_2],
                                     take_waves(crystal, 2), n_cell, atol=1e-6)
    moment_2 = {tuple(numpy.round(position, 4)): moment
                for l_site in sites.values() for position, moment in l_site}
    assert moment_2.keys() == moment_1.keys()
    for key, moment in moment_1.items():
        numpy.testing.assert_allclose(moment_2[key], moment, atol=1e-6)

    # free Fourier components of both waves: the group gives one moment per site
    sites = generate_moments_multi_k(operations, centrings, [k_1, k_2],
                                     take_waves(crystal, 2, rng), n_cell, atol=1e-8)
    assert max(numpy.abs(moment).max() for l_site in sites.values()
               for _, moment in l_site) > 0.1
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        modulation = calc_moment_modulation_by_dictionary(crystal.get_dictionary())
    for label, l_site in sites.items():
        if label not in modulation["label"]:
            continue
        position = numpy.stack([site[0] for site in l_site], axis=1)
        moment = numpy.stack([site[1] for site in l_site], axis=1)
        numpy.testing.assert_allclose(calc_moments_at_positions(
            position, modulation["fract_xyz"], modulation["amplitude"],
            modulation["k_vector"]), moment, atol=1e-8)


def test_second_propagation_vector_equal_to_first_is_refused():
    from test_magnetic_k_crystal_tools import form_crystal
    crystal = form_crystal("P n m a", (5., 6., 7., 90, 90, 90),
                           [("Mn1", "Mn", 0.1, 0.25, 0.3)], magnetic="Mn1")
    crystal.add_propagation_vector(0.3, 0, 0, "Mn1", irrep="1", moment=2.)
    with pytest.raises(Exception, match="already"):
        crystal.add_propagation_vector(-0.3, 0, 0, "Mn1", irrep="1", moment=2.)


def form_rcif_supercell_of_reference() -> str:
    """The supercell (3a, 4b, c) with the moments of the direct generation."""
    from test_magnetic_k_supercell import EXPERIMENT
    sites = take_reference_sites()
    rows, moments = [], []
    for label, l_site in sites.items():
        for i_site, (position, moment) in enumerate(l_site):
            name = f"{label}_{i_site + 1}"
            xyz = position / N_CELL
            rows.append(f"{name} Fe {xyz[0]:.10f} {xyz[1]:.10f} {xyz[2]:.10f} Biso 0.0 1.0")
            moments.append((name, moment))
    lines = ["global_\n\ndata_phase\n"]
    lines += [f"_cell_length_{n} {m * v}" for n, m, v in zip("abc", N_CELL, CELL[:3])]
    lines += [f"_cell_angle_{n} {v}" for n, v in zip(("alpha", "beta", "gamma"), CELL[3:])]
    lines += ["", "loop_\n_space_group_symop_magn_operation_xyz\n"
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


def test_powder_pattern_of_two_propagation_vectors_equals_supercell():
    """Powder profile of the (3+2) description (without O1, B = 0) equals
    that of the supercell written atom by atom (ordered moments, mcif)."""
    from test_magnetic_k_supercell import EXPERIMENT, calc_signal
    # written by to_cif: k = 1/3 must keep its digits
    text = form_crystal().to_cif().replace("data_k2", "data_phase")
    text = "\n".join(line for line in text.split("\n") if not line.strip().startswith("O1 "))
    rcif_k = "global_\n\n" + text + "\n" + EXPERIMENT + "50.0 1.0 1.0\n"
    ttheta = numpy.linspace(8., 95., 1500)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        signal_k, phase_k = calc_signal(rcif_k, ttheta)
        signal_s, _ = calc_signal(form_rcif_supercell_of_reference(), ttheta)
    assert set(phase_k["k_index_magnetic_k"]) == {0, 1}
    n_cell_sq = float(numpy.prod(N_CELL)) ** 2
    numpy.testing.assert_allclose(signal_s, n_cell_sq * signal_k,
                                  rtol=1e-6, atol=1e-6 * signal_s.max())
