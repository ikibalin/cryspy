"""Layers without correlations between them in powder experiments (loop layer_correlation).

The order of a phase that exists only inside the layers scatters along rods
(Warren line shape) instead of Bragg reflections. Checks: the rods keep the
integrated intensity of the Bragg reflections they replace (sum rule, an
independent check of the normalization, the form factor and the
orientation factor), the refinement finds the correlation length and the
scale, files without the loop are not changed.
"""
import os

import numpy
import pytest

import cryspy

DIR_TESTS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
F_MCIF_PD = os.path.join(DIR_TESTS, "functional", "mcif_test", "pd_test.rcif")

# square layers of Mn2+ (a = 4 A, layers 6 A apart), checkerboard
# antiferromagnet with the moments in the plane
SQUARE = """
global_

data_mn
_cell_length_a 4.0
_cell_length_b 4.0
_cell_length_c 6.0
_cell_angle_alpha 90.
_cell_angle_beta 90.
_cell_angle_gamma 90.

loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
_atom_site_adp_type
_atom_site_B_iso_or_equiv
Mn1 Mn2+ 0.0 0.0 0.0 1.0 Biso 0.3
Mn2 Mn2+ 0.5 0.5 0.0 1.0 Biso 0.3

loop_
_atom_site_scat_label
_atom_site_scat_type_symbol
_atom_site_scat_Lande
_atom_site_scat_kappa
Mn1 Mn2+ 2.0 1.0
Mn2 Mn2+ 2.0 1.0

loop_
_atom_site_moment_label
_atom_site_moment_crystalaxis_x
_atom_site_moment_crystalaxis_y
_atom_site_moment_crystalaxis_z
Mn1 4.0 0.0 1.0
Mn2 -4.0 0.0 -1.0

loop_
_space_group_symop_magn_operation_xyz
_space_group_symop_magn_operation_id
x,y,z,+1 1

loop_
_space_group_symop_magn_centering_xyz
_space_group_symop_magn_centering_id
x,y,z,+1 1

data_pd
_setup_wavelength 2.4
_setup_offset_2theta 0.0
_setup_radiation neutrons
_setup_K 0.0
_setup_cthm 0.91
_pd_instr_resolution_U 0.1
_pd_instr_resolution_V -0.1
_pd_instr_resolution_W 0.1
_pd_instr_resolution_X 0.0
_pd_instr_resolution_Y 0.0
_pd_instr_reflex_asymmetry_p1 0.1
_pd_instr_reflex_asymmetry_p2 0.0
_pd_instr_reflex_asymmetry_p3 0.0
_pd_instr_reflex_asymmetry_p4 0.0
_range_2theta_min 10.
_range_2theta_max 120.

loop_
_phase_label
_phase_scale
mn 1.0

loop_
_pd_background_2theta
_pd_background_intensity
5. 5.
130. 5.

loop_
_pd_meas_2theta
_pd_meas_intensity
_pd_meas_intensity_sigma
"""
STEP = 0.05


def form_square():
    ttheta = numpy.arange(10., 120. + 1e-6, STEP)
    lines = [f"{tth:.2f} 1.0 1.0" for tth in ttheta]
    return cryspy.str_to_globaln(SQUARE + "\n".join(lines) + "\n")


def take_pd(obj):
    return [item for item in obj.items if isinstance(item, cryspy.Pd)][0]


def calc_net(obj):
    """ttheta (degrees) and the signal without the background."""
    cryspy.rhochi_no_refinement(obj)
    pd_proc = take_pd(obj).pd_proc
    return (numpy.array(pd_proc.ttheta, dtype=float),
            numpy.array(pd_proc.intensity_plus_net, dtype=float)
            + numpy.array(pd_proc.intensity_minus_net, dtype=float))


@pytest.mark.parametrize("contribution", ["magnetic", "nuclear", "nuclear and magnetic"])
def test_rods_keep_the_intensity_of_the_bragg_reflections(contribution):
    """Sum rule: the rods spread the intensity of the Bragg reflections of
    the rods along them; over a range covering many reflections of a rod
    the integrated intensity is the same (moments with an out-of-plane
    component: the orientation factor changes along the rods)."""
    ttheta, net_bragg = calc_net(form_square())
    obj = form_square()
    cryspy.add_layer_correlation(obj, length=300., normal="0 0 1", contribution=contribution)
    net_rods = calc_net(obj)[1]
    assert numpy.abs(net_rods - net_bragg).max() > 0.1 * net_bragg.max()
    flag = (ttheta > 30.) & (ttheta < 118.)
    integral_bragg = net_bragg[flag].sum() * STEP
    integral_rods = net_rods[flag].sum() * STEP
    assert integral_rods == pytest.approx(integral_bragg, rel=0.02)
    # the front of the rod (1 0): nothing below, the line rises at 2theta_10
    tth_10 = numpy.degrees(2. * numpy.arcsin(2.4 / 8.))
    if contribution != "nuclear":
        below = (ttheta > tth_10 - 4.) & (ttheta < tth_10 - 2.)
        tail = (ttheta > tth_10 + 3.) & (ttheta < tth_10 + 5.)
        assert net_rods[below].max() < 0.02 * net_rods[tail].min()


def test_no_row_and_cif():
    """Without the loop nothing changes; the loop is written and read back."""
    obj = form_square()
    cryspy.add_layer_correlation(obj, length=120., normal="0 0 1")
    text = take_pd(obj).layer_correlation.to_cif()
    assert "_layer_correlation_length" in text and "magnetic" in text
    obj_2 = cryspy.str_to_globaln(obj.to_cif())
    row = take_pd(obj_2).layer_correlation.items[0]
    assert (row.normal_h, row.normal_k, row.normal_l) == (0, 0, 1)
    assert row.length == pytest.approx(120.)
    assert row.length_refinement and row.scale_refinement
    assert numpy.allclose(calc_net(obj_2)[1], calc_net(obj)[1], atol=2e-5)


def test_refine_length_and_scale():
    """Noise-free pattern of L = 80 A, scale 1: refined from L = 200, scale 0.7."""
    obj = form_square()
    cryspy.add_layer_correlation(obj, length=80., normal="0 0 1")
    cryspy.rhochi_no_refinement(obj)
    pd = take_pd(obj)
    total = numpy.array(pd.pd_proc.intensity_plus_net) + \
        numpy.array(pd.pd_proc.intensity_minus_net) + 5.
    for item, value in zip(pd.pd_meas.items, total):
        item.intensity = float(value)
        item.intensity_sigma = float(numpy.sqrt(value))
    row = pd.layer_correlation.items[0]
    row.length, row.scale = 200., 0.7
    cryspy.rhochi_rietveld_refinement(obj)
    assert row.length == pytest.approx(80., rel=0.02)
    assert row.scale == pytest.approx(1., rel=0.01)
    assert row.length_sigma > 0.
    report = cryspy.report_layer_correlation(obj)
    assert "correlation length L = 80" in report
    assert "34.92" in report


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_magnetic_space_group_with_centring():
    """A crystal with magnetic symmetry operations and centring: the layers
    suggested from the structure; the magnetic Bragg intensity becomes rods."""
    obj = cryspy.load_file(F_MCIF_PD)
    net_bragg = calc_net(obj)[1]
    cryspy.add_layer_correlation(obj, length=100.)
    pd = take_pd(obj)
    row = pd.layer_correlation.items[0]
    assert tuple(row.get_normal_index()) != (0, 0, 0)
    net_rods = calc_net(obj)[1]
    assert numpy.all(numpy.isfinite(net_rods))
    assert numpy.abs(net_rods - net_bragg).max() > 0.


def test_errors():
    obj = form_square()
    with pytest.raises(ValueError):
        cryspy.add_layer_correlation(obj, contribution="electronic")
    with pytest.raises(ValueError):
        cryspy.add_layer_correlation(obj, normal="1 0")
