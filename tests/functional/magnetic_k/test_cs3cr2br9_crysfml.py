"""Single-crystal magnetic intensities of Cs3Cr2Br9 (k = (1/4, 1/4, 0)) from
the CrysFML example CCB.cfl / CCB.int (see the header of the rcif file).

The crystal is given with its paramagnetic space group P 63/m m c and is
converted by ``Crystal.add_propagation_vector``; the representation of
CrysFML (basis functions BASR) is Γ2 of the group of k, the Fourier
component of Cr is (0, 0, 1.12064) μB.
"""
import os

import numpy
import pytest

import cryspy

F_RCIF = os.path.join(os.path.dirname(__file__), "cs3cr2br9_crysfml.rcif")
S_CR = 1.12064


def load_converted(irrep: int = 2):
    obj = cryspy.load_file(F_RCIF)
    crystal = obj.crystal_ccb
    crystal.add_propagation_vector(0.25, 0.25, 0., "Cr", irrep=irrep)
    crystal.atom_site_moment_fourier["Cr"].real_z = S_CR
    # CrysFML MagOptim uses B = 0 for the magnetic atom
    crystal.atom_site["Cr"].b_iso_or_equiv = 0.
    return obj


def calc_agreement(obj):
    cryspy.rhochi_rietveld_refinement(obj)
    refln = obj.diffrn_ccbmag.diffrn_refln
    i_obs = numpy.array(refln.intensity)
    i_calc = numpy.array(refln.intensity_calc)
    sigma = numpy.array(refln.intensity_sigma)
    chi_sq_n = numpy.mean(numpy.square((i_obs - i_calc) / sigma))
    r_f2 = numpy.abs(i_obs - i_calc).sum() / i_obs.sum()
    return chi_sq_n, r_f2, i_calc


def test_irreps_of_group_of_k():
    obj = cryspy.load_file(F_RCIF)
    report = obj.crystal_ccb.report_irreps_of_k(0.25, 0.25, 0., "Cr")
    assert "4 of 24 operations" in report
    assert "Γ4" in report and "Γ5" not in report


def test_conversion_splits_atoms_and_keeps_parent_group():
    obj = load_converted()
    crystal = obj.crystal_ccb
    assert not crystal.is_attribute("space_group")
    assert crystal.parent_space_group.name_hm_alt == "P 63/m m c"
    labels = list(crystal.atom_site.label)
    # 6h and 12k split into three orbits of the group of k, 2b and 4f do not
    assert labels == ["Cs1", "Cs2", "Cr", "Br1", "Br1_2", "Br1_3",
                      "Br2", "Br2_2", "Br2_3"]
    assert list(crystal.space_group_symop_magn_k.phase) == [0., 0., 0.5, 0.5]
    # the converted crystal is written and read back
    text = "global_\n" + crystal.to_cif() + "\n" + obj.diffrn_ccbmag.to_cif()
    obj_2 = cryspy.str_to_globaln(text)
    assert obj_2.crystal_ccb.is_attribute("atom_site_moment_fourier")


def test_crysfml_data_are_described():
    chi_sq_n, r_f2, i_calc = calc_agreement(load_converted())
    # CrysFML: the same model; 62 reflections, 12 of them zero
    assert chi_sq_n < 3.
    assert r_f2 < 0.09
    i_obs = numpy.array(cryspy.load_file(F_RCIF).diffrn_ccbmag.diffrn_refln.intensity)
    assert numpy.all(i_calc[i_obs == 0.] < 0.5)


def test_other_representation_does_not_describe_data():
    _, r_f2, _ = calc_agreement(load_converted(irrep=1))
    assert r_f2 > 0.5


def test_moment_direction_is_refined_along_c():
    obj = load_converted()
    fourier = obj.crystal_ccb.atom_site_moment_fourier["Cr"]
    fourier.real_x, fourier.real_y = 0.3, 0.1
    fourier.real_x_refinement, fourier.real_y_refinement = True, True
    cryspy.rhochi_rietveld_refinement(obj)
    assert fourier.real_x == pytest.approx(0., abs=0.15)
    assert fourier.real_y == pytest.approx(0., abs=0.15)


def test_magnetic_domains():
    obj = load_converted()
    assert "Magnetic domains: 3" in obj.crystal_ccb.report_magnetic_domains()
    chi_sq_one = calc_agreement(obj)[0]
    report = obj.diffrn_ccbmag.define_magnetic_domains(obj.crystal_ccb)
    assert "3 magnetic domains" in report
    assert len(obj.diffrn_ccbmag.jana_twin.items) == 2
    # the measured satellites belong to the arm k of the first domain; the
    # other two domains scatter elsewhere, so only the scale changes
    chi_sq_three = calc_agreement(obj)[0]
    assert chi_sq_three == pytest.approx(chi_sq_one, rel=1e-4)


def test_simulation_has_satellites():
    obj = load_converted()
    obj.items[:] = [obj.crystal_ccb]
    cryspy.simulation_single_crystal(obj)
    refln = obj.diffrn_sc.diffrn_refln
    index_h = numpy.array(refln.index_h)
    flag_satellite = numpy.abs(index_h - numpy.round(index_h)) > 1e-3
    intensity = numpy.array(refln.intensity)
    assert flag_satellite.sum() > 100
    assert intensity[flag_satellite].max() > 0.5
