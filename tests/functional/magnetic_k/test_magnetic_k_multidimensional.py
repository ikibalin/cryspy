"""Magnetic structures of multidimensional small representations: every
isotropy subgroup (mode) of the two-dimensional representations of Pnma
at k = (1/2, 0, 0) gives moments that belong to the chosen representation.

The moments of the magnetic supercell (2a, b, c) are made from the
superspace group written by ``add_propagation_vector`` with the real-space
msCIF rule of ``magnetic_k_reference`` (not by CrysPy) and projected on the
representations by the projection operator of the space group (Izyumov,
Naish, Ozerov, "Neutron diffraction of magnetic materials" (1991)) with the
characters of the X1 and X2 representations of Pnma (Bilbao REPRES:
2 for the identity, -2 and 2 for the mirror x,-y+1/2,z, zero otherwise).
"""
import numpy
import pytest

from cryspy.A_functions_base.magnetic_k_symmetry import parse_symmform
from cryspy.A_functions_base.structure_factor_magnetic_k import parse_operation_xyz
from magnetic_k_reference import generate_moments
from test_magnetic_k_crystal_tools import form_crystal

K_VECTOR = (0.5, 0., 0.)
N_CELL = numpy.array([2, 1, 1])
PNMA = ["x,y,z", "-x+1/2,-y,z+1/2", "-x,y+1/2,-z", "x+1/2,-y+1/2,-z+1/2",
        "-x,-y,-z", "x+1/2,y,-z+1/2", "x,-y+1/2,z", "-x+1/2,y+1/2,z+1/2"]
CHARACTERS = {"1": [2, 0, 0, 0, 0, 0, -2, 0], "2": [2, 0, 0, 0, 0, 0, 2, 0]}


def project(sites, characters):
    """Projection of the supercell moments on a representation.

    The group of the supercell has the 8 operations times the translations
    0 and a; the character of the translation a is exp(-2 pi i k.a) = -1.
    """
    position = numpy.array([site[0] for site in sites])
    moment = numpy.array([site[1] for site in sites])
    out = numpy.zeros_like(moment)
    n_group = 2 * len(PNMA)
    for xyz, character in zip(PNMA, characters):
        rotation, translation = parse_operation_xyz(xyz)
        rotation = rotation.reshape(3, 3)
        for shift, character_shift in ((0, 1.), (1, -1.)):
            new_position = numpy.mod(position @ rotation.T + translation + [shift, 0, 0], N_CELL)
            for site_position, site_moment in zip(new_position, moment):
                distance = numpy.mod(position - site_position + 0.5 * N_CELL, N_CELL) - 0.5 * N_CELL
                j = numpy.flatnonzero(numpy.all(numpy.abs(distance) < 1e-5, axis=1))[0]
                # axial vector: det(R) R m
                out[j] += (2. / n_group) * character * character_shift * \
                    numpy.linalg.det(rotation) * rotation @ site_moment
    return moment, out


def take_atoms(crystal, rng):
    """Random values of the free Fourier parameters (symmform symbols per atom)."""
    d_symbol = {}
    atoms = {}
    for atom in crystal.atom_site.items:
        fourier = numpy.zeros((3, 2))
        for item in crystal.atom_site_moment_fourier.items:
            if item.atom_site_label != atom.label:
                continue
            for i_part, symmform in enumerate((item.cos_symmform, item.sin_symmform)):
                value = 0.
                for symbol, coefficient in parse_symmform(symmform).items():
                    key = (atom.label, symbol)
                    if key not in d_symbol:
                        d_symbol[key] = rng.uniform(-2., 2.)
                    value += coefficient * d_symbol[key]
                fourier["xyz".index(item.axis), i_part] = value
        atoms[atom.label] = ((atom.fract_x, atom.fract_y, atom.fract_z), fourier)
    return atoms


@pytest.mark.parametrize("position", [(0.1, 0.25, 0.3), (0., 0., 0.), (0.1, 0.2, 0.3)])
@pytest.mark.parametrize("irrep", ["1", "2"])
def test_modes_belong_to_representation(position, irrep):
    other = "2" if irrep == "1" else "1"
    rng = numpy.random.default_rng(7)
    report = form_crystal("P n m a", (5., 6., 7., 90, 90, 90), [("Mn1", "Mn", *position)],
                          magnetic="Mn1").report_irreps_of_k(*K_VECTOR, "Mn1")
    assert "Γ1 (dim 2)" in report and "Γ2 (dim 2)" in report
    n_nonzero = 0
    for mode in (1, 2, 3, 4):
        assert f"### Γ{irrep} mode {mode}" in report
        crystal = form_crystal("P n m a", (5., 6., 7., 90, 90, 90),
                               [("Mn1", "Mn", *position)], magnetic="Mn1")
        crystal.add_propagation_vector(*K_VECTOR, "Mn1", irrep=irrep, mode=mode, moment=2.)
        operations = [item.algebraic for item in crystal.space_group_symop_magn_ssg_operation.items]
        centrings = [item.algebraic for item in crystal.space_group_symop_magn_ssg_centering.items]
        assert len(operations) < len(PNMA)  # a subgroup of the group of k
        # generate_moments checks that the moments at one site agree
        sites = sum(generate_moments(operations, centrings, K_VECTOR,
                                     take_atoms(crystal, rng), N_CELL).values(), [])
        moment, projected = project(sites, CHARACTERS[irrep])
        _, projected_other = project(sites, CHARACTERS[other])
        numpy.testing.assert_allclose(projected, moment, atol=1e-6)
        numpy.testing.assert_allclose(projected_other, 0., atol=1e-6)
        n_nonzero += numpy.abs(moment).max() > 1e-3
    assert n_nonzero >= 3


def test_mode_of_one_dimensional_representation_is_refused():
    crystal = form_crystal("P m m m", (5., 6., 7., 90, 90, 90), [("Mn1", "Mn", 0.1, 0.2, 0.3)],
                           magnetic="Mn1")
    with pytest.raises(Exception, match="mode"):
        crystal.add_propagation_vector(*K_VECTOR, "Mn1", irrep="1", mode=2, moment=2.)
