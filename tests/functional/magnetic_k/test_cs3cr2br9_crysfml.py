"""Single-crystal magnetic intensities of Cs3Cr2Br9 (k = (1/4, 1/4, 0)) from
the CrysFML example CCB.cfl / CCB.int (see the header of the rcif file).

The crystal is given with its paramagnetic space group P 63/m m c and is
converted by ``Crystal.add_propagation_vector``; the representation of
CrysFML (basis functions BASR) is Γ2 of the group of k, the Fourier
component of Cr is (0, 0, 1.12064) μB in the FullProf convention
``m = S exp(-2 pi i k.L) + c.c.``, so the amplitude of the moment is
``|Mcos + i Msin| = 2 |S|`` (the global phase does not change the
intensities of one arm).
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
    crystal.atom_site_moment_fourier.take_rows("Cr")["z"].cos = 2. * S_CR
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
    # 6h and 12k split into two orbits of the superspace group, 2b and 4f
    # do not
    assert labels == ["Cs1", "Cs2", "Cr", "Br1", "Br1_2", "Br2", "Br2_2"]
    # group of k (4) and the operations mapping k to -k (4)
    assert len(crystal.space_group_symop_magn_ssg_operation.items) == 8
    assert crystal.cell_wave_vector.items[0].x == pytest.approx(0.25)
    rows = crystal.atom_site_moment_fourier.take_rows("Cr")
    assert (rows["y"].cos_symmform, rows["z"].sin_symmform) == ("-mxc1", "0")
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
    rows = obj.crystal_ccb.atom_site_moment_fourier.take_rows("Cr")
    rows["x"].cos, rows["x"].sin = 0.3, 0.1
    rows["y"].cos, rows["y"].sin = -0.3, 0.1
    rows["x"].cos_refinement, rows["x"].sin_refinement = True, True
    cryspy.rhochi_rietveld_refinement(obj)
    # y follows x by the symmetry (symmform)
    assert rows["y"].cos == pytest.approx(-rows["x"].cos)
    assert rows["x"].cos == pytest.approx(0., abs=0.3)
    assert rows["x"].sin == pytest.approx(0., abs=0.3)


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
    # a second call does not take the domains for structural twins
    obj.diffrn_ccbmag.define_magnetic_domains(obj.crystal_ccb)
    assert len(obj.diffrn_ccbmag.jana_twin.items) == 2


def test_magnetic_domains_keep_structural_twins():
    obj = load_converted()
    diffrn = obj.diffrn_ccbmag
    diffrn.define_magnetic_domains(obj.crystal_ccb)
    twin = diffrn.jana_twin.items[0]
    # a structural twin (c -> -c, a -> -a) with a fifth of the volume
    for i in (1, 2, 3):
        for j in (1, 2, 3):
            setattr(twin, f"matrix_{i}{j}", {1: -1., 2: 1., 3: -1.}[i] if i == j else 0.)
    twin.fraction = 0.2
    diffrn.jana_twin.items[:] = [twin]
    report = diffrn.define_magnetic_domains(obj.crystal_ccb)
    assert "combined with 1 structural twins" in report
    fractions = numpy.array(diffrn.jana_twin.fraction)
    assert len(fractions) == 5
    # the twin keeps its volume, shared by its three magnetic domains
    assert fractions.sum() == pytest.approx(1. - 0.8 / 3., abs=1e-5)


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


def test_simulation_of_satellites_only():
    obj = load_converted()
    obj.items[:] = [obj.crystal_ccb]
    cryspy.simulation_single_crystal(obj, wavelength=2.4, ttheta_max=60., satellites_only=True)
    refln = obj.diffrn_sc.diffrn_refln
    index_hkl = numpy.array([refln.index_h, refln.index_k, refln.index_l], dtype=float)
    assert index_hkl.shape[1] > 0
    assert numpy.all(numpy.any(numpy.abs(index_hkl - numpy.round(index_hkl)) > 1e-3, axis=0))
