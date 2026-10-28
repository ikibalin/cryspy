"""3D view of crystal structures: the scene (pure numpy) and the py3Dmol view."""
import dataclasses
import os

import numpy
import pytest

import cryspy
from cryspy.A_functions_base.structure_scene import (
    calc_atoms_in_cells, calc_atoms_in_unit_cell, calc_cell_edges,
    calc_structure_scene)

DIR_EXAMPLES = os.path.join(
    os.path.dirname(__file__), os.pardir, os.pardir, "examples")
F_MCIF = os.path.join(
    DIR_EXAMPLES, "rhochi_unpolarized_powder_1d_mcif_URu0.96Rh0.04Si2",
    "rhochi_unpolarized_powder_1d_mcif_URu0.96Rh0.04Si2.rcif")
F_PBSO4 = os.path.join(
    DIR_EXAMPLES, "rhochi_unpolarized_powder_1d_PbSO4",
    "rhochi_unpolarized_powder_1d_PbSO4.rcif")

IDENTITY = [0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1]
INVERSION = [0, 0, 0, 1, -1, 0, 0, 0, -1, 0, 0, 0, -1]
BODY_CENTRING = [1, 1, 1, 2, 1, 0, 0, 0, 1, 0, 0, 0, 1]


def take_crystal(f_name):
    rcif = cryspy.load_file(f_name)
    return [item for item in rcif.items if type(item).__name__ == "Crystal"][0]


def test_atoms_in_unit_cell_merges_special_positions():
    symm = numpy.array([IDENTITY, INVERSION, BODY_CENTRING]).transpose()
    fract = numpy.array([[0., 0.1], [0., 0.2], [0., 0.3]])
    xyz, index, moment = calc_atoms_in_unit_cell(symm, fract)
    # (0,0,0): identity and inversion give one site; centring a second one
    assert numpy.sum(index == 0) == 2
    assert numpy.sum(index == 1) == 3
    assert moment is None


def test_moments_are_axial_vectors_with_time_reversal():
    inversion_tr = INVERSION + [-1]
    symm = numpy.array([IDENTITY + [1], inversion_tr]).transpose()
    fract = numpy.array([[0.1], [0.2], [0.3]])
    moment = numpy.array([[1.], [2.], [3.]])
    xyz, index, moment_out = calc_atoms_in_unit_cell(symm, fract, moment)
    # inversion keeps an axial vector, time reversal flips it
    numpy.testing.assert_allclose(moment_out[:, 1], [-1., -2., -3.])
    numpy.testing.assert_allclose(xyz[:, 1], [0.9, 0.8, 0.7])


def test_atoms_in_cells_with_boundary():
    fract = numpy.array([[0.], [0.], [0.]])
    xyz, index = calc_atoms_in_cells(fract, (2, 1, 1), flag_boundary=True)
    assert xyz.shape[1] == 3 * 2 * 2   # corners of a 2x1x1 block
    xyz, index = calc_atoms_in_cells(fract, (2, 1, 1), flag_boundary=False)
    assert xyz.shape[1] == 2


def test_cell_edges():
    assert calc_cell_edges((1, 1, 1)).shape == (2, 3, 12)
    assert calc_cell_edges((2, 1, 1), flag_inner=True).shape == (2, 3, 4 + 6 + 6)


def test_scene_cartesian_cubic():
    ucp = numpy.array([4., 4., 4., numpy.pi/2, numpy.pi/2, numpy.pi/2])
    symm = numpy.array([IDENTITY]).transpose()
    scene = calc_structure_scene(
        ucp, symm, numpy.array([[0.5], [0.5], [0.5]]),
        atom_moment_crystalaxis=numpy.array([[0.], [0.], [5.]]))
    numpy.testing.assert_allclose(scene["atom_xyz"][:, 0], [2., 2., 2.])
    numpy.testing.assert_allclose(scene["moment_xyz"][:, 0], [0., 0., 5.],
                                  atol=1e-12)


def test_crystal_scene_from_space_group():
    crystal = take_crystal(F_PBSO4)
    scene = crystal.get_structure_scene(
        cryspy.StructureViewSettings(show_boundary_atoms=False))
    multiplicity = crystal.get_dictionary()["atom_multiplicity"]
    assert scene["atom_index"].size == multiplicity.sum()
    assert scene["moment_xyz"] is None
    assert set(scene["atom_element"]) == {"Pb", "S", "O"}


def test_crystal_scene_from_mcif():
    crystal = take_crystal(F_MCIF)
    scene = crystal.get_structure_scene(
        cryspy.StructureViewSettings(show_boundary_atoms=False))
    moment = scene["moment_xyz"]
    magnetic = ~numpy.isnan(moment[0])
    # U1_1 at the corner and the body centre (moment -0.63 along c),
    # U1_2 at x = 1/6, 1/3, 2/3, 5/6 (+0.505 along c)
    assert magnetic.sum() == 6
    numpy.testing.assert_allclose(numpy.sort(moment[2, magnetic]),
                                  [-0.63, -0.63] + [0.505] * 4)
    numpy.testing.assert_allclose(moment[:2, magnetic], 0., atol=1e-12)


def test_view_structure_html():
    pytest.importorskip("py3Dmol")
    crystal = take_crystal(F_MCIF)
    html = crystal.get_structure_view_html(
        cryspy.StructureViewSettings(n_cells_a=2, show_labels=True),
        js="3Dmol-min.js")
    assert "3Dmol-min.js" in html
    assert html.count("addArrow") >= 3 + 6 * 2  # axes + moments
    assert "U1_1" in html


def test_settings_describe_every_field():
    for field in dataclasses.fields(cryspy.StructureViewSettings):
        assert field.metadata["description"].endswith(".")
