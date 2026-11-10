"""Group of k of a centred lattice: R^T k - k must be a vector of the
reciprocal lattice of the centred cell (h + k even for C), not only integer."""
import numpy

import cryspy
from cryspy.A_functions_base.magnetic_k_symmetry import calc_k_group_flags
from cryspy.A_functions_base.structure_factor_magnetic_k import (
    calc_satellite_index_hkl, is_k_equivalent_to_minus_k)

CRYSTAL = """global_
data_tbte3
_cell_length_a 4.31
_cell_length_b 25.4
_cell_length_c 4.31
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_space_group_name_H-M_alt "C m c m"
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Tb Tb3+ 0.0 0.169455 0.25
"""
C_CENTRING = numpy.array([[0., 0.5], [0., 0.5], [0., 0.]])


def test_k_group_flags_c_centring():
    mirror_x = numpy.diag([-1., 1., 1.]).flatten()[:, numpy.newaxis]
    k_vector = numpy.array([0.5, 0.5, 0.])
    # R^T k - k = (-1, 0, 0): integer, but not a vector of the C lattice
    assert calc_k_group_flags(mirror_x, k_vector)[0]
    assert not calc_k_group_flags(mirror_x, k_vector, centring_translations=C_CENTRING)[0]


def test_k_minus_k_c_centring():
    assert is_k_equivalent_to_minus_k([0.5, 0.5, 0.], centring_translations=C_CENTRING)
    # 2k = (1, 0, 0) is forbidden by the C centring
    assert not is_k_equivalent_to_minus_k([0.5, 0., 0.], centring_translations=C_CENTRING)
    index_hkl, sign_k = calc_satellite_index_hkl(
        0., 0.3, numpy.array([4., 4., 4., numpy.pi / 2, numpy.pi / 2, numpy.pi / 2]),
        numpy.array([0.5, 0., 0.]), centring_translations=C_CENTRING)
    assert set(sign_k.tolist()) == {-1, 1}


def test_report_irreps_cmcm_k_half_half():
    crystal = cryspy.str_to_globaln(CRYSTAL).crystal_tbte3
    report = crystal.report_irreps_of_k(0.5, 0.5, 0., "Tb")
    assert "4 of 8 operations" in report
    assert "## Γ4" in report
