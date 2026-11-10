"""Magnetic structure factor of ordered moments (mCIF atom_site_moment) in
oblique cells, compared with the direct sum of ``magnetic_k_reference``
(independent of CrysPy) over all atoms of the cell.

The moments are components along the unit vectors of a, b, c (mCIF
crystalaxis). CrysPy applied the normalisation matrix twice and the
rotation of the operations to Cartesian components, so F_M was wrong in
monoclinic, triclinic and hexagonal cells (up to 36 % of |F|^2 for
gamma = 120°).
"""
import warnings

import numpy
import pytest

import cryspy
from cryspy.A_functions_base.structure_factor_magnetic import (
    calc_f_m_perp_ordered_by_dictionary,
)
from cryspy.A_functions_base.structure_factor_magnetic_k import parse_operation_xyz
from magnetic_k_reference import calc_f_m_perp

ATOMS = [("Fe1", (0.1, 0.2, 0.3), (1.5, -0.7, 2.1)),
         ("Fe2", (0.6, 0.7, 0.1), (-0.4, 1.9, 0.5))]
CASES = {
    "monoclinic": ((4.1, 5.3, 6.2, 90., 97., 90.), ["x,y,z,+1", "-x,y,-z,+1"]),
    "triclinic": ((4.1, 5.3, 6.2, 80., 97., 105.), ["x,y,z,+1", "-x,-y,-z,+1"]),
    "hexagonal": ((4.1, 4.1, 6.2, 90., 90., 120.),
                  ["x,y,z,+1", "-y,x-y,z,+1", "-x+y,-x,z+1/2,-1"]),
}


def form_crystal(cell, operations):
    lines = ["global_", "data_c"]
    lines += [f"_cell_length_{name} {value}" for name, value in zip("abc", cell[:3])]
    lines += [f"_cell_angle_{name} {value}"
              for name, value in zip(("alpha", "beta", "gamma"), cell[3:])]
    lines += ["loop_", "_space_group_symop_magn_operation_xyz",
              "_space_group_symop_magn_operation_id"]
    lines += [f"{text} {i + 1}" for i, text in enumerate(operations)]
    lines += ["loop_", "_space_group_symop_magn_centering_xyz",
              "_space_group_symop_magn_centering_id", "x,y,z,+1 1"]
    lines += ["loop_", "_atom_site_label", "_atom_site_type_symbol", "_atom_site_fract_x",
              "_atom_site_fract_y", "_atom_site_fract_z", "_atom_site_adp_type",
              "_atom_site_B_iso_or_equiv", "_atom_site_occupancy"]
    lines += [f"{label} Fe {x} {y} {z} Biso 0 1" for label, (x, y, z), _ in ATOMS]
    lines += ["loop_", "_atom_site_scat_label", "_atom_site_scat_type_symbol",
              "_atom_site_scat_lande"]
    lines += [f"{label} Fe3+ 2.0" for label, _, _ in ATOMS]
    lines += ["loop_", "_atom_site_moment_label", "_atom_site_moment_crystalaxis_x",
              "_atom_site_moment_crystalaxis_y", "_atom_site_moment_crystalaxis_z"]
    lines += [f"{label} {m[0]} {m[1]} {m[2]}" for label, _, m in ATOMS]
    return "\n".join(lines)


def take_sites(cell, operations):
    """All atoms of the cell with their moments along the unit vectors of
    a, b, c: M' = theta det(R) D R D^-1 M, D = diag(|a|, |b|, |c|)."""
    lengths = numpy.diag(cell[:3])
    sites = {}
    for label, position, moment in ATOMS:
        sites[label] = []
        for text in operations:
            xyz, theta = text.rsplit(",", 1)
            rotation, translation = parse_operation_xyz(xyz)
            rotation = rotation.reshape(3, 3)
            rotation_unit = lengths @ rotation @ numpy.linalg.inv(lengths)
            sites[label].append((
                numpy.mod(rotation @ numpy.array(position) + translation, 1.),
                int(theta) * numpy.linalg.det(rotation) * rotation_unit @ numpy.array(moment)))
    return sites


@pytest.mark.parametrize("case", sorted(CASES))
def test_ordered_moments_in_oblique_cell(case):
    cell, operations = CASES[case]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        dict_crystal = cryspy.str_to_globaln(form_crystal(cell, operations)).get_dictionary()["crystal_c"]
    # form factor 1 and B = 0, as in the reference
    for key in dict_crystal:
        if key.endswith("j0_parameters") or key.endswith("j2_parameters"):
            dict_crystal[key][:] = 0.
            if key.endswith("j0_parameters"):
                dict_crystal[key][6] = 1.
    index_hkl = numpy.random.default_rng(2).integers(-3, 4, size=(3, 20)).astype(float)
    f_m_perp = calc_f_m_perp_ordered_by_dictionary(dict_crystal, {"index_hkl": index_hkl})[0]
    f_reference = calc_f_m_perp(take_sites(cell, operations), numpy.ones(3, dtype=int),
                                cell[:3], cell[3:], index_hkl)
    assert numpy.abs(f_reference).max() > 0.05
    # F_M_perp is defined up to the phase of the origin: compare |F|^2 per component sum
    numpy.testing.assert_allclose((numpy.abs(f_m_perp) ** 2).sum(axis=0),
                                  (numpy.abs(f_reference) ** 2).sum(axis=0), atol=1e-10)
