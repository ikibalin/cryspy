"""3D view of crystal structures: the scene (pure numpy) and the py3Dmol view."""
import dataclasses
import os

import numpy
import pytest

import cryspy
from cryspy.A_functions_base.structure_scene import (
    calc_adp_probability_factor, calc_atoms_in_cells, calc_atoms_in_unit_cell,
    calc_cell_edges, calc_ellipsoid_matrix_by_adp,
    calc_ellipsoid_matrix_by_susceptibility, calc_ellipsoid_mesh,
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
    xyz, index, moment, symm_index = calc_atoms_in_unit_cell(symm, fract)
    # (0,0,0): identity and inversion give one site; centring a second one
    assert numpy.sum(index == 0) == 2
    assert numpy.sum(index == 1) == 3
    assert moment is None


def test_moments_are_axial_vectors_with_time_reversal():
    inversion_tr = INVERSION + [-1]
    symm = numpy.array([IDENTITY + [1], inversion_tr]).transpose()
    fract = numpy.array([[0.1], [0.2], [0.3]])
    moment = numpy.array([[1.], [2.], [3.]])
    xyz, index, moment_out, symm_index = calc_atoms_in_unit_cell(symm, fract, moment)
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


F_EMSO = os.path.join(DIR_EXAMPLES, "only_for_test", "emso.rcif")


def test_tensors_rotate_with_symmetry():
    # a fourfold axis along z turns U_xx into U_yy
    ucp = numpy.array([4., 4., 5., numpy.pi/2, numpy.pi/2, numpy.pi/2])
    four = [0, 0, 0, 1, 0, -1, 0, 1, 0, 0, 0, 0, 1]
    symm = numpy.array([IDENTITY, four]).transpose()
    u = numpy.diag([0.01, 0.02, 0.03])[:, :, numpy.newaxis]
    scene = calc_structure_scene(
        ucp, symm, numpy.array([[0.1], [0.2], [0.3]]), atom_tensor_xyz=u,
        flag_boundary=False)
    tensors = numpy.moveaxis(scene["tensor_xyz"], 2, 0)
    numpy.testing.assert_allclose(numpy.diag(tensors[0]), [0.01, 0.02, 0.03])
    numpy.testing.assert_allclose(numpy.diag(tensors[1]), [0.02, 0.01, 0.03])


def test_adp_probability_factor():
    assert calc_adp_probability_factor(0.5) == pytest.approx(1.5382, abs=1e-4)
    assert calc_adp_probability_factor(0.99) == pytest.approx(3.3682, abs=1e-3)


def test_ellipsoid_matrices():
    u = numpy.diag([0.04, 0.01, 0.09])[:, :, numpy.newaxis]
    matrix = calc_ellipsoid_matrix_by_adp(u, 0.5)[:, :, 0]
    numpy.testing.assert_allclose(
        numpy.diag(matrix), 1.53817 * numpy.array([0.2, 0.1, 0.3]), rtol=1e-4)
    chi = numpy.array([[0., 2., 0.], [2., 0., 0.], [0., 0., -1.]])[:, :, numpy.newaxis]
    matrix = calc_ellipsoid_matrix_by_susceptibility(chi, 3.)[:, :, 0]
    # the image of the sphere |B| = 3 T: semi-axes 6, 6, 3 μB
    numpy.testing.assert_allclose(numpy.linalg.eigvalsh(matrix), [3., 6., 6.])
    vertices, normals, faces = calc_ellipsoid_mesh(numpy.zeros(3), matrix)
    assert faces.max() < vertices.shape[1]
    numpy.testing.assert_allclose(numpy.sqrt(numpy.square(normals).sum(axis=0)), 1.)


def test_susceptibility_ellipsoids_view():
    pytest.importorskip("py3Dmol")
    crystal = take_crystal(F_EMSO)
    settings = cryspy.StructureViewSettings(
        atom_display="susceptibility ellipsoids", show_boundary_atoms=False)
    scene = crystal.get_structure_scene(settings)
    magnetic = ~numpy.isnan(scene["tensor_xyz"][0, 0])
    assert magnetic.sum() == 9  # Er1 at 9d of R-3m
    html = crystal.get_structure_view_html(settings, js="3Dmol-min.js")
    assert html.count("addCustom") == 9
    with pytest.raises(ValueError):
        cryspy.StructureViewSettings(atom_display="cubes")


def test_adp_ellipsoids_and_dummy_atoms():
    crystal = take_crystal(F_PBSO4)
    settings = cryspy.StructureViewSettings(
        atom_display="ADP ellipsoids", show_boundary_atoms=False)
    scene = crystal.get_structure_scene(settings)
    u = scene["tensor_xyz"]
    # Uiso: isotropic tensors U * 1; S has a negative Uiso and stays a sphere
    flags = numpy.isfinite(u[0, 0])
    assert (~flags).sum() == 4
    numpy.testing.assert_allclose(u[0, 1, flags], 0., atol=1e-12)
    numpy.testing.assert_allclose(u[0, 0, flags], u[2, 2, flags])
    crystal.atom_site.items[1].type_symbol = "Dummy"
    n_all = scene["atom_index"].size
    scene = crystal.get_structure_scene(settings)
    assert "S" not in scene["atom_element"]
    assert scene["atom_index"].size == n_all - 4
    scene = crystal.get_structure_scene(
        dataclasses.replace(settings, hide_dummy_atoms=False))
    assert scene["atom_index"].size == n_all


def test_to_p1_keeps_structure_factors():
    crystal = take_crystal(F_PBSO4)
    crystal_p1 = crystal.to_p1()
    assert len(crystal_p1.atom_site.items) == crystal.get_dictionary()["atom_multiplicity"].sum()
    index_hkl = numpy.array([[1, 0, 2, 3], [0, 1, 1, 2], [1, 2, 5, 1]], dtype=float)
    numpy.testing.assert_allclose(numpy.abs(crystal_p1.calc_f_nucl(index_hkl)),
                                  numpy.abs(crystal.calc_f_nucl(index_hkl)), rtol=1e-4, atol=1e-8)


def test_to_p1_magnetic():
    crystal = take_crystal(F_MCIF)
    crystal_p1 = crystal.to_p1()
    moments = [(item.label, item.crystalaxis_z) for item in crystal_p1.atom_site_moment.items]
    assert len(moments) == 6
    assert crystal_p1.is_attribute("space_group_symop_magn_operation")


def test_assistant_tools_unit_cell_and_ellipsoids():
    rcif = cryspy.load_file(F_EMSO)
    text = cryspy.call_assistant_tool(rcif, "report_ellipsoids",
                                      {"kind": "susceptibility", "unit_cell": True})["text"]
    assert "Er1_9" in text and "3.8939" in text
    text = cryspy.call_assistant_tool(rcif, "list_atoms_in_unit_cell", {})["text"]
    assert "66 atoms" in text
    result = cryspy.call_assistant_tool(rcif, "convert_to_p1", {})
    assert result["flag_object_changed"]
    assert any(item.data_name == "EMSO_p1" for item in rcif.items)
