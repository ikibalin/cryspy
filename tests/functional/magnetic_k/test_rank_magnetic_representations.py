"""Choice of the representation of k by the data: a powder pattern made by
one representation of P 1 21/m 1 at k = (1/2, 0, 0) is fitted best by
that representation (zero chi_sq), the other three leave a misfit.

Γ2 and Γ3 allow the same Fourier components of Mn1 (moments in the ac
plane) and differ only by the sign of the moment of the atom at -r, so
the test checks that the procedure uses the full superspace group, not the
site symmetry only.
"""
import numpy
import pytest

import cryspy
from cryspy.D_functions_item_loop.function_1_magnetic_k import (
    take_irrep_candidates,
    take_k_group,
)
from cryspy.exceptions import CrysPyError
from test_magnetic_k_crystal_tools import form_crystal

CRYSTAL = """global_
data_phase
_space_group_name_H-M_alt 'P 1 21/m 1'
_cell_length_a 5.0
_cell_length_b 6.0
_cell_length_c 7.0
_cell_angle_alpha 90
_cell_angle_beta 100
_cell_angle_gamma 90
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_adp_type
_atom_site_B_iso_or_equiv
_atom_site_occupancy
Mn1 Mn 0.1 0.25 0.3 Biso 0.3 1.0
O1 O 0.35 0.05 0.2 Biso 0.5 1.0
loop_
_atom_site_scat_label
_atom_site_scat_type_symbol
_atom_site_scat_lande
Mn1 Mn3+ 2.0
"""
POWDER = """
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
_range_2theta_max 60.0
loop_
_phase_label
_phase_scale
_phase_igsize
phase 1.0 0.0
loop_
_pd_background_2theta
_pd_background_intensity
5.0 0.0
60.0 0.0
loop_
_pd_meas_2theta
_pd_meas_intensity
_pd_meas_intensity_sigma
"""
K_VECTOR = (0.5, 0., 0.)
TTHETA = numpy.arange(6., 59., 0.2)


def form_object(intensity=None, sigma=None):
    if intensity is None:
        intensity, sigma = numpy.ones_like(TTHETA), numpy.ones_like(TTHETA)
    rows = "".join(f"{t:.2f} {y:.6f} {s:.6f}\n" for t, y, s in zip(TTHETA, intensity, sigma))
    return cryspy.str_to_globaln(CRYSTAL + POWDER + rows)


@pytest.fixture(scope="module")
def measured():
    """Pattern of Γ2 with 3 μB, sigma of counting statistics."""
    truth = form_object()
    truth.crystal_phase.add_propagation_vector(*K_VECTOR, "Mn1", "2", 3.0)
    cryspy.rhochi_no_refinement(truth)
    proc = truth.pd_pd.pd_proc
    intensity = numpy.array(proc.intensity_plus_net) + numpy.array(proc.intensity_minus_net)
    return intensity, numpy.sqrt(intensity) + 0.5


def take_table(report):
    """[(representation, chi_sq/n)] of the table, in its order."""
    rows = [line.split("|") for line in report.splitlines()
            if line.startswith("| ") and line[2].isdigit()]
    return [(row[2].strip(), float(row[3])) for row in rows]


def test_generating_representation_is_best(measured):
    obj = form_object(*measured)
    report = cryspy.rhochi_rank_magnetic_representations(obj, *K_VECTOR, "Mn1", 1.0)
    table = take_table(report)
    assert [name for name, _ in table][0] == "Γ2"
    assert sorted(name for name, _ in table if "⊕" not in name) == ["Γ1", "Γ2", "Γ3", "Γ4"]
    # the sums of two one-dimensional representations are tried as well
    assert len([name for name, _ in table if "⊕" in name]) == 6
    assert "Recommended: Γ2" in report
    # the data are the model of Γ2: chi_sq is zero within the refinement
    # tolerance; the others miss the satellites by far more than the sigmas
    assert table[0][1] < 1e-2
    # the models without Γ2 leave a misfit
    assert min(value for name, value in table if "Γ2" not in name) > 0.5
    # the object is not changed
    assert obj.crystal_phase.is_attribute("space_group")
    assert not obj.crystal_phase.is_attribute("atom_site_moment_fourier")


def test_keep_best(measured):
    obj = form_object(*measured)
    cryspy.rhochi_rank_magnetic_representations(
        obj, *K_VECTOR, "Mn1", 1.0, keep_best=True)
    crystal = obj.crystal_phase
    assert crystal.is_attribute("atom_site_moment_fourier")
    assert cryspy.rhochi_no_refinement(obj)["chi_sq"] / TTHETA.size < 1e-2


def test_zero_start_moment_is_refused():
    with pytest.raises(CrysPyError, match="start moment"):
        cryspy.rhochi_rank_magnetic_representations(form_object(), *K_VECTOR, "Mn1", 0.)


def test_candidates_include_modes():
    """Pnma at k = (1/2, 0, 0): two representations of dimension 2 with
    four modes each (test_magnetic_k_multidimensional)."""
    crystal = form_crystal("P n m a", (5., 6., 7., 90, 90, 90),
                           [("Mn1", "Mn", 0.1, 0.25, 0.3)], magnetic="Mn1")
    k_group = take_k_group(crystal.space_group, numpy.array(K_VECTOR))
    candidates = take_irrep_candidates(k_group, numpy.array(K_VECTOR))
    assert candidates == [(1, mode) for mode in (1, 2, 3, 4)] + [
        (2, mode) for mode in (1, 2, 3, 4)]


TBTE3 = """global_
data_tbte3
_space_group_name_H-M_alt 'C m c m'
_cell_length_a 4.31
_cell_length_b 25.6
_cell_length_c 4.31
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_adp_type
_atom_site_B_iso_or_equiv
_atom_site_occupancy
Tb Tb 0.0 0.8306 0.25 Biso 0.3 1.0
Te1 Te 0.0 0.0705 0.25 Biso 0.3 1.0
"""


def test_magnetic_ion_is_added():
    """A crystal without atom_site_scat (a nuclear refinement): the ion is
    asked for, then atom_site_scat gets it with g_J = 3/2 of Tb3+ (4f8,
    7F6)."""
    crystal = cryspy.str_to_globaln(TBTE3).crystal_tbte3
    with pytest.raises(CrysPyError, match="magnetic_ion"):
        crystal.add_propagation_vector(0.5, 0.5, 0, "Tb", "1", 8.0)
    report = crystal.add_propagation_vector(0.5, 0.5, 0, "Tb", "1", 8.0,
                                            magnetic_ion="Tb3+")
    assert "atom_site_scat added" in report
    item = crystal.atom_site_scat["Tb"]
    assert item.type_symbol == "Tb3+"
    assert item.lande == pytest.approx(1.5)
    with pytest.raises(CrysPyError, match="form factor"):
        cryspy.str_to_globaln(TBTE3).crystal_tbte3.add_propagation_vector(
            0.5, 0.5, 0, "Tb", "1", 8.0, magnetic_ion="Tb9+")


def test_lande_factors_of_rare_earths():
    """g_J of the free 4f ions (Jensen, Mackintosh, table 1.1)."""
    from cryspy.D_functions_item_loop.function_1_magnetic_k import calc_lande_factor
    expected = {"Ce3+": 6 / 7, "Nd3+": 8 / 11, "Sm3+": 2 / 7, "Gd3+": 2., "Tb3+": 1.5,
                "Dy3+": 4 / 3, "Ho3+": 5 / 4, "Er3+": 6 / 5, "Tm3+": 7 / 6,
                "Yb3+": 8 / 7, "Fe3+": 2.}
    for ion, value in expected.items():
        assert calc_lande_factor(ion) == pytest.approx(value, abs=1e-6), ion
