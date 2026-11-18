"""Texture by spherical harmonics in the powder refinement."""
import os
import re

import numpy
import pytest

import cryspy

RCIF_PBSO4 = """
data_pbso4
_cell_length_a 8.4779
_cell_length_b 5.3968
_cell_length_c 6.9581
_cell_angle_alpha 90.0
_cell_angle_beta 90.0
_cell_angle_gamma 90.0
_space_group_name_H-M_alt "P n m a"
_space_group_IT_coordinate_system_code abc
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_adp_type
_atom_site_B_iso_or_equiv
_atom_site_occupancy
Pb Pb 0.1872 0.25 0.1662 Biso 0.8 1.0
S S 0.0643 0.25 0.6826 Biso 0.5 1.0
O1 O 0.9079 0.25 0.5960 Biso 1.5 1.0
O2 O 0.1926 0.25 0.5417 Biso 1.0 1.0
O3 O 0.0804 0.0289 0.8073 Biso 0.8 1.0

data_pd
_setup_wavelength 1.912
_setup_field 0.0
_setup_offset_2theta 0.0
_setup_radiation neutrons
_pd_instr_resolution_U 0.12
_pd_instr_resolution_V -0.33
_pd_instr_resolution_W 0.28
_pd_instr_resolution_X 0.0
_pd_instr_resolution_Y 0.15
_range_2theta_min 15.0
_range_2theta_max 95.0
loop_
_phase_label
_phase_scale
pbso4 1.0
loop_
_pd_background_2theta
_pd_background_intensity
10.0 50.0
100.0 50.0
loop_
_pd_meas_2theta
_pd_meas_intensity
_pd_meas_intensity_sigma
{points}
"""

# mmm, cylindrical, order 4: C_2^00, C_2^02, C_4^00, C_4^02, C_4^04
COEFFICIENTS = {(2, 0, 0): 0.9, (2, 0, 2): -0.6, (4, 0, 0): 0.5, (4, 0, 2): 0.3, (4, 0, 4): -0.4}


def form_pbso4(intensity=None):
    ttheta = numpy.arange(15., 95., 0.1)
    if intensity is None:
        intensity = numpy.full(ttheta.shape, 100.)
    sigma = numpy.sqrt(numpy.abs(intensity)) + 1.
    points = "\n".join(f"{tth:.2f} {i:.4f} {s:.4f}" for tth, i, s in zip(ttheta, intensity, sigma))
    return cryspy.str_to_globaln(RCIF_PBSO4.format(points=points))


def take_intensity_calc(obj) -> numpy.ndarray:
    cryspy.rhochi_no_refinement(obj)
    pd_proc = obj.pd_pd.pd_proc
    return numpy.array(pd_proc.intensity_plus_net) + numpy.array(pd_proc.intensity_minus_net) \
        + numpy.array(pd_proc.intensity_bkg_calc)


def set_coefficients(obj, d_value: dict, flag_refinement: bool = False):
    for item in obj.pd_pd.texture_sh_coefficient.items:
        item.value = d_value[(item.l, item.m, item.n)]
        item.value_refinement = flag_refinement


def test_add_texture_and_random_powder():
    obj = form_pbso4()
    intensity_random = take_intensity_calc(obj)
    text = cryspy.add_texture_spherical_harmonics(obj, order_max=4)
    assert "Laue class mmm" in text and "5 coefficients" in text
    lmn = {(item.l, item.m, item.n) for item in obj.pd_pd.texture_sh_coefficient.items}
    assert lmn == set(COEFFICIENTS)
    # all coefficients 0: the pattern of the random powder
    assert numpy.allclose(take_intensity_calc(obj), intensity_random, rtol=1e-10)
    report = cryspy.report_texture_spherical_harmonics(obj)
    assert "J = 1.000" in report
    # order 6 keeps the values of order 4
    set_coefficients(obj, COEFFICIENTS)
    cryspy.add_texture_spherical_harmonics(obj, order_max=6)
    d_value = {(item.l, item.m, item.n): item.value
               for item in obj.pd_pd.texture_sh_coefficient.items}
    assert len(d_value) == 5 + 4
    assert all(d_value[key] == value for key, value in COEFFICIENTS.items())


def test_refinement_recovers_the_coefficients():
    obj = form_pbso4()
    cryspy.add_texture_spherical_harmonics(obj, order_max=4)
    set_coefficients(obj, COEFFICIENTS)
    intensity_texture = take_intensity_calc(obj)

    obj = form_pbso4(intensity_texture)
    cryspy.add_texture_spherical_harmonics(obj, order_max=4)
    chi_sq_start = cryspy.rhochi_no_refinement(obj)["chi_sq"]
    cryspy.rhochi_rietveld_refinement(obj)
    chi_sq = cryspy.rhochi_no_refinement(obj)["chi_sq"]
    assert chi_sq < 1e-4 * chi_sq_start
    for item in obj.pd_pd.texture_sh_coefficient.items:
        assert item.value == pytest.approx(COEFFICIENTS[(item.l, item.m, item.n)], abs=1e-3)
        assert item.value_sigma > 0.
    report = cryspy.report_texture_spherical_harmonics(obj)
    j_index = 1. + sum(value ** 2 / (2 * key[0] + 1) for key, value in COEFFICIENTS.items())
    assert f"J = {j_index:.3f}" in report

    # the CIF keeps the texture
    obj_2 = cryspy.str_to_globaln(obj.to_cif())
    assert numpy.allclose(take_intensity_calc(obj_2), take_intensity_calc(obj))


def insert_in_pd(obj, text: str):
    text_cif, n_sub = re.subn(r"(\ndata_pd\s*\n)", lambda match: match.group(1) + text + "\n",
                              obj.to_cif(), count=1)
    assert n_sub == 1
    return cryspy.str_to_globaln(text_cif)


def test_goniometer_changes_the_correction():
    obj = form_pbso4()
    cryspy.add_texture_spherical_harmonics(obj, order_max=4)
    set_coefficients(obj, COEFFICIENTS)
    intensity_vertical = take_intensity_calc(obj)
    # the default orientation is chi = 90 (vertical cylinder axis)
    obj_2 = insert_in_pd(obj, "_sample_goniometer_chi 90.0")
    assert obj_2.pd_pd.is_attribute("sample_goniometer")
    assert numpy.allclose(take_intensity_calc(obj_2), intensity_vertical)
    obj_3 = insert_in_pd(obj, "_sample_goniometer_chi 0.0")
    assert not numpy.allclose(take_intensity_calc(obj_3), intensity_vertical, rtol=1e-3)


def test_march_dollase_and_spherical_harmonics_are_exclusive():
    obj = form_pbso4()
    cryspy.add_texture_spherical_harmonics(obj, order_max=2)
    obj_2 = insert_in_pd(obj, "loop_\n_texture_label\n_texture_g_1\n_texture_g_2\n"
                         "_texture_h_ax\n_texture_k_ax\n_texture_l_ax\npbso4 0.9 0.0 0 0 1")
    with pytest.raises(ValueError):
        cryspy.rhochi_no_refinement(obj_2)
    with pytest.raises(ValueError):
        cryspy.add_texture_spherical_harmonics(obj_2, order_max=4)


def test_wrong_arguments():
    obj = form_pbso4()
    with pytest.raises(ValueError):
        cryspy.add_texture_spherical_harmonics(obj, sample_symmetry="cubic")
    with pytest.raises(ValueError):
        cryspy.add_texture_spherical_harmonics(obj, order_max=40)
    with pytest.raises(ValueError):
        cryspy.add_texture_spherical_harmonics(obj, phase_label="no_phase")


DIR = os.path.dirname(__file__)
F_TOF = os.path.join(DIR, "test_tof", "model_gauss.rcif")
F_PD2D = os.path.join(DIR, "..", "..", "examples", "rhochi_polarized_powder_2d_DyAl_5K_5T",
                      "rhochi_polarized_powder_2d_DyAl_5K_5T.rcif")


def test_tof_texture_refinement():
    obj = cryspy.file_to_globaln(F_TOF)
    chi_sq_random = cryspy.rhochi_no_refinement(obj)["chi_sq"]
    cryspy.add_texture_spherical_harmonics(obj, order_max=2)
    tof = obj.tof_tof1
    assert cryspy.rhochi_no_refinement(obj)["chi_sq"] == pytest.approx(chi_sq_random, rel=1e-10)
    d_value = {}
    for item, value in zip(tof.texture_sh_coefficient.items, (0.6,)):
        item.value, item.value_refinement = value, False
        d_value[(item.l, item.m, item.n)] = value
    # only the texture is refined (the profile is that of the simulation)
    cryspy.call_assistant_tool(obj, "fix_parameters")
    cryspy.rhochi_no_refinement(obj)
    tof_proc = tof.tof_proc
    intensity = numpy.array(tof_proc.intensity_plus_net) + numpy.array(
        tof_proc.intensity_minus_net) + numpy.array(tof_proc.intensity_bkg_calc)
    for item, value in zip(tof.tof_meas.items, intensity):
        item.intensity = float(value)
    for item in tof.texture_sh_coefficient.items:
        item.value, item.value_refinement = 0., True
    cryspy.rhochi_rietveld_refinement(obj)
    for item in tof.texture_sh_coefficient.items:
        assert item.value == pytest.approx(d_value[(item.l, item.m, item.n)], abs=1e-3)


def test_pd2d_texture():
    obj = cryspy.file_to_globaln(F_PD2D)
    pd2d = obj.pd2d_powder2
    pd2d.items.remove(pd2d.texture)
    # the first calculation of this example differs from the next ones by
    # 4e-7 (also without the texture); compare with a repeated one
    cryspy.rhochi_no_refinement(obj)
    chi_sq_random = cryspy.rhochi_no_refinement(obj)["chi_sq"]
    cryspy.add_texture_spherical_harmonics(obj, order_max=6)
    assert cryspy.rhochi_no_refinement(obj)["chi_sq"] == pytest.approx(chi_sq_random, rel=1e-10)
    for item in pd2d.texture_sh_coefficient.items:
        item.value = 0.5
    chi_sq = cryspy.rhochi_no_refinement(obj)["chi_sq"]
    assert numpy.isfinite(chi_sq) and not numpy.isclose(chi_sq, chi_sq_random, rtol=1e-4)


def form_two_experiments():
    obj = form_pbso4()
    text = obj.to_cif()
    i_pd = text.index("data_pd")
    i_end = text.index("data_pbso4") if text.index("data_pbso4") > i_pd else len(text)
    block = text[i_pd:i_end]
    return cryspy.str_to_globaln(text + "\n" + block.replace("data_pd", "data_pd2", 1))


def test_coefficients_are_shared_by_experiments():
    obj = form_two_experiments()
    text = cryspy.add_texture_spherical_harmonics(obj, order_max=4)
    assert "shared" in text
    l_mark = [[item.value_mark for item in pd.texture_sh_coefficient.items]
              for pd in (obj.pd_pd, obj.pd_pd2)]
    assert l_mark[0] == l_mark[1] and len(set(l_mark[0])) == 5
    # one experiment: no marks
    obj = form_two_experiments()
    cryspy.add_texture_spherical_harmonics(obj, order_max=4, experiment="pd2")
    assert not obj.pd_pd.is_attribute("texture_sh")
    assert all(not getattr(item, "value_mark", "") for item in obj.pd_pd2.texture_sh_coefficient.items)


def test_new_sample_symmetry_keeps_only_m_0():
    obj = form_pbso4()
    cryspy.add_texture_spherical_harmonics(obj, order_max=2, sample_symmetry="mmm")
    for item in obj.pd_pd.texture_sh_coefficient.items:
        item.value = 0.5
    cryspy.add_texture_spherical_harmonics(obj, order_max=2, sample_symmetry="2/m")
    for item in obj.pd_pd.texture_sh_coefficient.items:
        assert item.value == (0.5 if item.m == 0 else 0.)


def take_net_intensity(obj) -> numpy.ndarray:
    cryspy.rhochi_no_refinement(obj)
    pd_proc = obj.pd_pd.pd_proc
    return numpy.array(pd_proc.intensity_plus_net) + numpy.array(pd_proc.intensity_minus_net)


def test_march_dollase_texture_of_one_of_two_phases():
    """A phase without texture next to a phase with the March-Dollase
    texture keeps the intensities of a random powder."""
    march = ("loop_\n_texture_label\n_texture_g_1\n_texture_g_2\n_texture_h_ax\n"
             "_texture_k_ax\n_texture_l_ax\npbso4 0.8 0.0 0 0 1")
    obj = form_pbso4()
    intensity_random = take_net_intensity(obj)
    intensity_march = take_net_intensity(insert_in_pd(obj, march))
    text = obj.to_cif()
    i_crystal = text.index("data_pbso4")
    crystal = text[i_crystal:].split("\ndata_")[0]
    text, n_sub = re.subn(r"\n[ \t]*pbso4[ \t]+1\.0[ \t]+0\.0[ \t]*\n", "\npbso4 0.5 0.0\npbso4b 0.5 0.0\n", text)
    assert n_sub == 1
    text += "\n" + crystal.replace("data_pbso4", "data_pbso4b", 1)
    obj_2 = insert_in_pd(cryspy.str_to_globaln(text), march)
    # pd_proc keeps 5 decimals
    assert numpy.allclose(take_net_intensity(obj_2), 0.5 * (intensity_random + intensity_march),
                          rtol=1e-6, atol=2e-5)
