"""From the paramagnetic crystal to the magnetic structure with single-crystal
satellites (the path of the recipe magnetic_structure_k, as for TbTe3 with
k = (0, 0, 1/2) and (1/2, 1/2, 0) measured on D10):

- a nuclear refinement with the satellites already in the object works;
- find_propagation_vectors gives k from the fractional indices;
- the ranking gives the experiments the magnetic domains (two arms of the
  star of k = (1/2, 1/2, 0) in C m c m) and tries the sums Γi⊕Γj; the
  generating model is found with its domains, R_F2 of the satellites zero.
"""
import warnings

import numpy
import pytest

import cryspy
from cryspy.procedure_rhochi.rhochi_magnetic_irreps import is_satellite_of_k

CRYSTAL = """global_
data_tbte3
_space_group_name_H-M_alt 'C m c m'
_cell_length_a 4.28
_cell_length_b 25.29
_cell_length_c 4.29
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
Tb Tb 0.0 0.1692 0.25 Biso 0.3 1.0
Te1 Te 0.0 0.2948 0.25 Biso 0.3 1.0
Te2 Te 0.0 0.5720 0.25 Biso 0.3 1.0
Te3 Te 0.0 0.9282 0.25 Biso 0.3 1.0
"""
DIFFRN = """
data_{name}
_setup_wavelength 2.369
_setup_field 0.0
_diffrn_orient_matrix_ub_11 0.2336
_diffrn_orient_matrix_ub_12 0.0
_diffrn_orient_matrix_ub_13 0.0
_diffrn_orient_matrix_ub_21 0.0
_diffrn_orient_matrix_ub_22 0.0395
_diffrn_orient_matrix_ub_23 0.0
_diffrn_orient_matrix_ub_31 0.0
_diffrn_orient_matrix_ub_32 0.0
_diffrn_orient_matrix_ub_33 0.2331
_diffrn_orient_matrix_type CCSL
_phase_label tbte3
_phase_scale 1.0
loop_
_diffrn_refln_index_h
_diffrn_refln_index_k
_diffrn_refln_index_l
_diffrn_refln_intensity
_diffrn_refln_intensity_sigma
{refln}
"""
K_VECTOR = numpy.array([0.5, 0.5, 0.])
NUCLEAR = numpy.array([[0, 2, 0], [0, 4, 0], [1, 1, 0], [1, 3, 1], [0, 0, 2], [2, 0, 0],
                       [1, 5, 1], [0, 6, 1], [2, 2, 1], [1, 7, 0]], dtype=float).transpose()
INTEGER = numpy.array([[0, 0, 0], [0, 1, 0], [0, 2, 1], [1, 0, 1], [0, 3, 0], [1, 1, 1],
                       [0, 4, 1], [-1, 2, 0], [0, 5, 1], [1, 3, 0], [-1, 4, 1], [0, 7, 0],
                       [1, -2, 1], [-1, -3, 0]], dtype=float).transpose()
SATELLITES = numpy.concatenate([INTEGER + K_VECTOR[:, None], INTEGER - K_VECTOR[:, None]], axis=1)


def form_diffrn(name, index_hkl, intensity=None) -> str:
    if intensity is None:
        intensity = numpy.ones(index_hkl.shape[1])
    refln = "\n".join(f"{h:.4f} {k:.4f} {l:.4f} {i:.6f} {0.01 + 0.01 * abs(i):.6f}"
                      for (h, k, l), i in zip(index_hkl.transpose(), intensity))
    return DIFFRN.format(name=name, refln=refln)


def calc_measured():
    """Intensities of Γ1 of k = (1/2, 1/2, 0) on Tb (5 μB) with its two
    magnetic domains in equal parts."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        obj = cryspy.str_to_globaln(CRYSTAL + form_diffrn("nuc", NUCLEAR) +
                                    form_diffrn("sat", SATELLITES))
        crystal = obj.crystal_tbte3
        crystal.add_propagation_vector(*K_VECTOR, "Tb", "1", 5.0, magnetic_ion="Tb3+")
        for diffrn in (obj.diffrn_nuc, obj.diffrn_sat):
            assert "magnetic domains" in diffrn.define_magnetic_domains(crystal)
        cryspy.rhochi_no_refinement(obj)
    return (numpy.array(obj.diffrn_nuc.diffrn_refln.intensity_calc),
            numpy.array(obj.diffrn_sat.diffrn_refln.intensity_calc))


@pytest.fixture(scope="module")
def paramagnetic():
    intensity_nuc, intensity_sat = calc_measured()
    return CRYSTAL + form_diffrn("nuc", NUCLEAR, intensity_nuc) + \
        form_diffrn("sat", SATELLITES, intensity_sat)


def test_nuclear_refinement_with_satellites_in_the_object(paramagnetic):
    """Refined nuclear parameters and satellites of a crystal without k:
    the satellites are calculated as zero (no crash on the derivatives)."""
    obj = cryspy.str_to_globaln(paramagnetic.replace("Te1 Te 0.0 0.2948", "Te1 Te 0.0 0.2948()"))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cryspy.rhochi_rietveld_refinement(obj)
    assert numpy.allclose(obj.diffrn_sat.diffrn_refln.intensity_calc, 0.)
    assert obj.crystal_tbte3.atom_site["Te1"].fract_y == pytest.approx(0.2948, abs=1e-3)


def test_find_propagation_vectors(paramagnetic):
    report = cryspy.find_propagation_vectors(cryspy.str_to_globaln(paramagnetic))
    assert "| (0.5, 0.5, 0) | 28 |" in report
    assert "10 reflections with integer indices" in report


def test_is_satellite_of_k():
    index_hkl = numpy.array([[0.5, 0.5, 0.], [1., 0., 0.], [0., 0., 0.5], [-0.5, 1.5, 1.]]).T
    assert list(is_satellite_of_k(index_hkl, K_VECTOR)) == [True, False, False, True]


def test_rank_with_domains_and_superpositions(paramagnetic):
    obj = cryspy.str_to_globaln(paramagnetic)
    report = cryspy.rhochi_rank_magnetic_representations(
        obj, *K_VECTOR, "Tb", 3.0, magnetic_ion="Tb3+")
    rows = [line.split("|") for line in report.splitlines()
            if line.startswith("| ") and line[2].isdigit()]
    names = [row[2].strip() for row in rows]
    assert "Γ1⊕Γ2" in names
    # Γ1 and the sums that contain it reproduce the data with the two
    # domains; the one with the fewest parameters is recommended
    assert "Recommended: Γ1 (1 Fourier parameters" in report or "Recommended: Γ1," in report
    row_gamma_1 = rows[names.index("Γ1")]
    assert float(row_gamma_1[4].strip(" %")) < 0.5
    assert row_gamma_1[6].strip() == "2"
    # without the domains the generating model misses the satellites of the second arm
    report = cryspy.rhochi_rank_magnetic_representations(
        obj, *K_VECTOR, "Tb", 3.0, magnetic_ion="Tb3+", with_domains=False,
        try_superpositions=False)
    rows = [line.split("|") for line in report.splitlines()
            if line.startswith("| ") and line[2].isdigit()]
    assert min(float(row[4].strip(" %")) for row in rows) > 1.


def test_keep_best_keeps_the_domains(paramagnetic):
    obj = cryspy.str_to_globaln(paramagnetic)
    report = cryspy.rhochi_rank_magnetic_representations(
        obj, *K_VECTOR, "Tb", 3.0, magnetic_ion="Tb3+", keep_best=True,
        try_superpositions=False)
    assert "2 magnetic domains" in report
    assert obj.diffrn_sat.is_attribute("jana_twin")
    assert obj.crystal_tbte3.is_attribute("atom_site_moment_fourier")
