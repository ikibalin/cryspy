"""Interatomic distances, angles, torsion angles and hydrogen bonds."""
import math
import os

import numpy
import pytest

import cryspy
from cryspy.A_functions_base.structure_geometry import (
    calc_angle_xyz, calc_geometry_sigma, calc_images_near_position,
    calc_torsion_xyz)

F_CH30 = os.path.join(os.path.dirname(__file__), os.pardir, "functional", "ch30", "ch30.rcif")

IDENTITY = [0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1]

S_CSCL = """global_
data_cscl
_cell_length_a 4.0
_cell_length_b 4.0
_cell_length_c 4.0
_cell_angle_alpha 90.0
_cell_angle_beta 90.0
_cell_angle_gamma 90.0
_space_group_name_H-M_alt "P m -3 m"
_space_group_IT_coordinate_system_code 1
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
Cs Cs 0.0 0.0 0.0 1.0
Cl Cl 0.5 0.5 0.5 1.0
"""

# a water dimer O1-H1...O2 along x in a large P1 cell
S_HBOND = """global_
data_water
_cell_length_a 10.0
_cell_length_b 10.0
_cell_length_c 10.0
_cell_angle_alpha 90.0
_cell_angle_beta 90.0
_cell_angle_gamma 90.0
_space_group_name_H-M_alt "P 1"
_space_group_IT_coordinate_system_code 1
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
O1 O 0.100 0.5 0.5 1.0
H1 H 0.196 0.5 0.5 1.0
O2 O2- 0.380 0.5 0.5 1.0
"""


def take_crystal(rcif):
    return [item for item in rcif.items if type(item).__name__ == "Crystal"][0]


def test_cscl_neighbours_and_angles():
    crystal = take_crystal(cryspy.str_to_globaln(S_CSCL))
    result = crystal.calc_distance_or_angle("Cs Cl")
    assert result["value"] == pytest.approx(4. * math.sqrt(3.) / 2., abs=1e-6)
    assert result["sigma"] == 0.
    # angle Cl-Cs-Cl between two of the eight nearest Cl (cube corners)
    result = crystal.calc_distance_or_angle("Cl Cs Cl")
    angle_cube = math.degrees(math.acos(1. / 3.))
    assert min(abs(result["value"] - angle) for angle in (
        angle_cube, 180. - angle_cube, 180.)) < 1e-4
    text = crystal.report_distances_and_angles("Cs", max_distance=3.5)
    assert "8 distances" in text and "28 angles" in text
    # the next Cs is a lattice translation away
    assert crystal.calc_distance_or_angle("Cs Cs")["value"] == pytest.approx(4.)


def test_images_cover_oblique_cell():
    # monoclinic cell: brute force over a large block of translations
    ucp = numpy.array([5., 6., 7., math.pi / 2, math.radians(115.), math.pi / 2])
    symm = numpy.array([IDENTITY]).transpose()
    atom, center = numpy.array([0.9, 0.2, 0.95]), numpy.array([0.1, 0.1, 0.1])
    images = calc_images_near_position(ucp, symm, atom, center, 6.)
    from cryspy.A_functions_base.structure_geometry import calc_m_m
    m_m = calc_m_m(ucp)
    grid = numpy.stack(numpy.meshgrid(*[numpy.arange(-5, 6)] * 3, indexing="ij")).reshape(3, -1)
    distance = numpy.linalg.norm(m_m @ (atom[:, None] + grid - center[:, None]), axis=0)
    expected = numpy.sort(distance[distance <= 6.])
    numpy.testing.assert_allclose(images["distance"], expected, atol=1e-10)


def test_torsion_sign_and_angle():
    p_1, p_2 = numpy.array([1., 0., 0.]), numpy.array([0., 0., 0.])
    p_3 = numpy.array([0., 0., 1.])
    # looking from 2 to 3 (along +z) bond 1-2 turned clockwise by 90° to 3-4
    assert calc_torsion_xyz(p_1, p_2, p_3, numpy.array([0., -1., 1.])) == pytest.approx(90.)
    assert calc_torsion_xyz(p_1, p_2, p_3, numpy.array([0., 1., 1.])) == pytest.approx(-90.)
    assert calc_torsion_xyz(p_1, p_2, p_3, numpy.array([1., 0., 1.])) == pytest.approx(0.)
    assert calc_angle_xyz(p_1, p_2, p_3) == pytest.approx(90.)


def test_distance_sigma_analytical():
    ucp = numpy.array([4., 5., 6., math.pi / 2, math.pi / 2, math.pi / 2])
    fract = numpy.array([[0.1, 0.4], [0.2, 0.2], [0.3, 0.3]])
    fract_sigma = numpy.array([[0.001, 0.002], [0., 0.], [0., 0.]])
    ucp_sigma = numpy.array([0.01, 0., 0., 0., 0., 0.])
    eye = numpy.eye(3)
    l_position = [(0, eye, numpy.zeros(3)), (1, eye, numpy.zeros(3))]
    sigma = calc_geometry_sigma(ucp, ucp_sigma, fract, fract_sigma, l_position)
    # d = a (x2 - x1): sigma^2 = (a s1)^2 + (a s2)^2 + ((x2 - x1) s_a)^2
    expected = math.sqrt((4. * 0.001) ** 2 + (4. * 0.002) ** 2 + (0.3 * 0.01) ** 2)
    assert sigma == pytest.approx(expected, rel=1e-5)


def test_ch30_distance_by_hand():
    crystal = take_crystal(cryspy.load_file(F_CH30))
    result = crystal.calc_distance_or_angle("N1 C3")
    # orthogonal cell a = b = 7.405, c = 8.6492, both atoms of atom_site
    d_xyz = numpy.array([(0.4238 - 0.39331) * 7.405, (0.9238 - 0.89331) * 7.405,
                         (0.240632 - 0.063492) * 8.6492])
    assert result["value"] == pytest.approx(numpy.linalg.norm(d_xyz), abs=1e-5)
    # only z is refined: dd/dz = c^2 dz / d
    derivative = 8.6492 ** 2 * (0.240632 - 0.063492) / numpy.linalg.norm(d_xyz)
    assert result["sigma"] == pytest.approx(derivative * math.hypot(80e-6, 88e-6), rel=1e-3)
    assert result["symmetry"] == ["x,y,z", "x,y,z"]
    text = crystal.report_distances_and_angles("N1 C1")
    assert text.startswith("Distance N1–C1 = 2.35")
    with pytest.raises(ValueError, match="X1"):
        crystal.calc_distance_or_angle("X1 C1")


def test_hydrogen_bonds():
    crystal = take_crystal(cryspy.str_to_globaln(S_HBOND))
    text = crystal.report_hydrogen_bonds()
    assert "| O1 | H1 | O2 | 0.960 | 1.840 | 2.800 | 180.00 |" in text
    assert "No hydrogen bonds" in crystal.report_hydrogen_bonds(max_distance_ha=1.5)


def test_geometry_loops_saved_and_read(tmp_path):
    rcif = cryspy.load_file(F_CH30)
    crystal = take_crystal(rcif)
    text = crystal.add_distances_and_angles(max_distance=1.7, torsions="N1 C3 C4 H4A")
    assert "12 distances" in text and "1 torsion angles" in text
    f_name = os.path.join(tmp_path, "ch30_geom.rcif")
    with open(f_name, "w") as fid:
        fid.write(rcif.to_cif())
    crystal_2 = take_crystal(cryspy.load_file(f_name))
    l_bond = [(item.atom_site_label_1, item.atom_site_label_2, item.distance)
              for item in crystal_2.geom_bond.items]
    assert ("N1", "C3", "1.5650(10)") in l_bond
    assert len(crystal_2.geom_angle.items) == 21
    assert crystal_2.geom_torsion.items[0].atom_site_label_4 == "H4A"
    # called again without torsions: the old loop is removed
    crystal_2.add_distances_and_angles(max_distance=1.7)
    assert not crystal_2.is_attribute("geom_torsion")


def test_assistant_tools():
    rcif = cryspy.load_file(F_CH30)
    text = cryspy.call_assistant_tool(rcif, "report_distances_and_angles",
                                      {"atoms": "N1 C1"})["text"]
    assert "Distance N1–C1" in text
    text = cryspy.call_assistant_tool(rcif, "report_distances_and_angles",
                                      {"atoms": "N1 C3", "cif": True})["text"]
    assert "_geom_bond.distance" in text and "_space_group_symop.operation_xyz" in text
    text = cryspy.call_assistant_tool(rcif, "list_atoms_in_unit_cell",
                                      {"atoms": "N1 C1"})["text"]
    assert "6 atoms of N1 C1" in text and "C3_1" not in text
    from cryspy.procedure_assistant.methods import find_matching_methods
    l_method = find_matching_methods(rcif, "What is the distance between N1 and C1?")
    assert l_method[0]["name"] == "report_distances_and_angles"


def test_restraint_difference_of_torsion_wraps():
    from cryspy.A_functions_base.structure_geometry import calc_restraint_difference
    assert calc_restraint_difference(179., -179., 4) == pytest.approx(-2.)
    assert calc_restraint_difference(179., -179., 3) == pytest.approx(358.)


def test_restraints_in_dictionary_and_rcif(tmp_path):
    from cryspy.procedure_rhochi.rhochi_by_dictionary import calc_restraints_chi_sq
    rcif = cryspy.load_file(F_CH30)
    crystal = take_crystal(rcif)
    text = crystal.add_restraint("N1 C3")
    assert "differs by 0.0000" in text
    assert crystal.restr_distance.items[0].target == pytest.approx(1.5650, abs=1e-4)
    crystal.add_restraint("N1 C3 C4", target="110", sigma=0.5)
    crystal.add_restraint("N1 C3 C4 H4A", target="-50")
    # the same atoms again replace the old row
    crystal.add_restraint("N1 C3", target="1.50", sigma=0.01)
    assert len(crystal.restr_distance.items) == 1
    value_angle = crystal.calc_distance_or_angle("N1 C3 C4")["value"]
    value_torsion = crystal.calc_distance_or_angle("N1 C3 C4 H4A")["value"]
    expected = (((1.5650 - 1.50) / 0.01) ** 2 + ((value_angle - 110.) / 0.5) ** 2
                + ((value_torsion + 50.) / 2.) ** 2)
    chi_sq = calc_restraints_chi_sq(rcif.get_dictionary())
    assert chi_sq == pytest.approx(expected, rel=1e-3)
    assert "Sum (diff/sigma)^2" in crystal.update_restraints()
    # a coordinate written outside [0, 1) gives the same restraint value
    crystal.atom_site["C4"].fract_x -= 1.
    assert calc_restraints_chi_sq(rcif.get_dictionary()) == pytest.approx(chi_sq, rel=1e-6)

    f_name = os.path.join(tmp_path, "ch30_restr.rcif")
    with open(f_name, "w") as fid:
        fid.write(rcif.to_cif())
    rcif_2 = cryspy.load_file(f_name)
    crystal_2 = take_crystal(rcif_2)
    assert crystal_2.restr_torsion.items[0].target == pytest.approx(-50.)
    assert crystal_2.restr_distance.items[0].get_site_symmetry() == [".", "."]
    assert calc_restraints_chi_sq(rcif_2.get_dictionary()) == pytest.approx(chi_sq, rel=1e-3)
    assert len(rcif_2.get_variable_names()) == len(rcif.get_variable_names())


def test_restraint_method_hint():
    from cryspy.procedure_assistant.methods import find_matching_methods
    rcif = cryspy.load_file(F_CH30)
    l_method = find_matching_methods(rcif, "Restrain the distance N1 C3 to 1.5 A")
    assert l_method[0]["name"] == "add_restraint"
