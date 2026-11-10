"""Crystal methods for structures with a propagation vector: the parameters
the data cannot define, the atoms of the cell, the parent cell, rounded
satellite indices, fractions and superpositions of representations, the
real-space moments."""
import os
import warnings

import numpy
import pytest

import cryspy
from cryspy.A_functions_base.structure_factor_magnetic_k import (
    calc_f_m_perp_k_by_dictionary,
    warn_not_satellites,
)
from cryspy.exceptions import CrysPyError

DIR = os.path.join(os.path.dirname(__file__), "magndata")


def form_crystal(space_group: str, cell, atoms, magnetic: str = "Mn"):
    text = f"""global_
data_cr
_space_group_name_H-M_alt '{space_group}'
_cell_length_a {cell[0]}
_cell_length_b {cell[1]}
_cell_length_c {cell[2]}
_cell_angle_alpha {cell[3]}
_cell_angle_beta {cell[4]}
_cell_angle_gamma {cell[5]}
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
_atom_site_B_iso_or_equiv
""" + "".join(f"{label} {symbol} {x} {y} {z} 1 0.3\n"
              for label, symbol, x, y, z in atoms) + f"""
loop_
_atom_site_scat_label
_atom_site_scat_type_symbol
_atom_site_scat_lande
{magnetic} Mn3+ 2.0
"""
    return cryspy.str_to_globaln(text).items[0]


def take_fixed(crystal):
    """Refine all Fourier components and return the flags kept."""
    for item in crystal.atom_site_moment_fourier.items:
        item.cos_refinement, item.sin_refinement = True, True
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ddict = crystal.get_dictionary()
    return ddict["flags_atom_fourier_cos"], ddict["flags_atom_fourier_sin"]


def test_fourier_gauge():
    """P1, k = (0, 1/2, 0): the atom sees only x4 = k.p + n/2, so 3 of 6
    parameters are defined. P 1 21 1 with the characters ±i: both atoms
    see different x4, all 6 are defined. Incommensurate P1: the global
    phase leaves 5 of 6."""
    crystal = form_crystal("P 1", (5, 6, 7, 90, 90, 90), [("Mn", "Mn", 0.11, 0.2, 0.33)])
    crystal.add_propagation_vector(0, 0.5, 0, "Mn", "1", 2.0)
    flags_cos, flags_sin = take_fixed(crystal)
    assert flags_cos.sum() + flags_sin.sum() == 3

    crystal = form_crystal("P 1 21 1", (5, 6, 7, 90, 100, 90), [("Mn", "Mn", 0.11, 0.2, 0.33)])
    crystal.add_propagation_vector(0, 0.5, 0, "Mn", "1", 2.0)
    flags_cos, flags_sin = take_fixed(crystal)
    assert flags_cos.sum() + flags_sin.sum() == 6

    crystal = form_crystal("P 1", (5, 6, 7, 90, 90, 90), [("Mn", "Mn", 0.11, 0.2, 0.33)])
    crystal.add_propagation_vector(0, 0.3, 0, "Mn", "1", 2.0)
    flags_cos, flags_sin = take_fixed(crystal)
    assert flags_cos.sum() + flags_sin.sum() == 5


def test_atoms_of_cell_from_superspace_group():
    crystal = form_crystal("P n m a", (5, 6, 7, 90, 90, 90), [("Mn", "Mn", 0, 0, 0.5)])
    crystal.add_propagation_vector(0, 0.28, 0, "Mn", "1", 2.0)
    elems = crystal.get_full_symm_elems()
    assert elems.shape[0] == 14
    obj = cryspy.GlobalN.make_container((), (cryspy.Crystal,), "global")
    obj.add_items([crystal])
    text = cryspy.call_assistant_tool(obj, "list_atoms_in_unit_cell", {})["text"]
    assert "4 atoms in the unit cell" in text


def test_parent_group_in_another_cell():
    """MAGNDATA 1.1.84 written in a non-standard cell of the parent C2/m:
    the parent cannot constrain the cell (gamma = 122.45 is kept)."""
    with open(os.path.join(DIR, "magndata_1.1.84.mcif")) as fid:
        text = fid.read()
    text = text.replace('"a,b,c;0,0,0"', '"a,c,-b;0,0,0"').replace(
        "_cell_angle_gamma", "_cell_angle_gamma 122.45\n#")
    label = cryspy.str_to_globaln(text).items[0].atom_site_moment_fourier.items[0].atom_site_label
    text += ("\nloop_\n_atom_site_scat_label\n_atom_site_scat_type_symbol\n"
             f"_atom_site_scat_lande\n{label} Fe3+ 2.0\n")
    crystal = cryspy.str_to_globaln(text).items[0]
    with pytest.warns(UserWarning, match="another cell"):
        ddict = crystal.get_dictionary()
    assert numpy.degrees(ddict["unit_cell_parameters"][5]) == pytest.approx(122.45)
    with pytest.raises(CrysPyError, match="another cell"):
        crystal.report_magnetic_domains()


def test_rounded_satellite_index_warns():
    crystal = form_crystal("P n m a", (5, 6, 7, 90, 90, 90), [("Mn", "Mn", 0, 0, 0.5)])
    crystal.add_propagation_vector(0, 0.2817, 0, "Mn", "1", 2.0)
    ddict = crystal.get_dictionary()
    index_hkl = numpy.array([[1., 1.], [0.2817, 0.282], [0., 0.]])
    with pytest.warns(UserWarning, match="not h ± k"):
        warn_not_satellites(index_hkl, ddict["magnetic_k_vector"])
    f_m_perp = calc_f_m_perp_k_by_dictionary(
        ddict, index_hkl, ddict["unit_cell_parameters"])
    assert numpy.abs(f_m_perp[:, 0]).max() > 0.01
    assert numpy.abs(f_m_perp[:, 1]).max() == 0.


def test_fraction_of_k_and_superposition():
    """'1/3' as k by the assistant, 0.3333 snapped to 1/3, Γ1⊕Γ2 keeps the
    operations with equal characters (half of the group of k of Pnma)."""
    from cryspy.procedure_assistant.tools import convert_argument
    assert convert_argument("1/3", {"type": "number"}) == pytest.approx(1. / 3.)
    crystal = form_crystal("P n m a", (5, 6, 7, 90, 90, 90), [("Mn", "Mn", 0.1, 0.2, 0.3)])
    report = crystal.report_irreps_of_k(0, 0.3333, 0, "Mn")
    assert "0.3333 → 1/3" in report
    report = crystal.add_propagation_vector(0, 0.28, 0, "Mn", "1,2", 2.0)
    assert "Γ1⊕Γ2" in report
    assert len(crystal.space_group_symop_magn_ssg_operation.items) == 4


def test_report_magnetic_moments_of_hexagonal_cell():
    """P 6, k || c, Γ4 (complex character of the 6-fold axis): a circular
    helix; the start moment 2 μB is the real length of the moment in
    hexagonal axes (components 2 ∓ 1.1547i along a and b)."""
    crystal = form_crystal("P 6", (3, 3, 5, 90, 90, 120), [("Mn", "Mn", 0, 0, 0)])
    crystal.add_propagation_vector(0, 0, 0.2137, "Mn", "4", 2.0)
    report = crystal.report_magnetic_moments()
    rows = [line for line in report.splitlines() if line.startswith("| 1 ")]
    assert len(rows) == 1
    assert "| 2.000–2.000 | circular helix" in rows[0]


def test_convert_space_group_to_superspace_keeps_nuclear_structure_factor():
    """A crystal with only space_group (no magnetic group) is written as the
    superspace group with k = 0; the nuclear structure factor is unchanged."""
    crystal = form_crystal("P n m a", (5.5, 7.6, 5.4, 90, 90, 90),
                           [("Mn", "Mn", 0.0, 0.0, 0.5), ("O1", "O", 0.48, 0.25, 0.08),
                            ("La", "La", 0.05, 0.25, 0.99)])
    index_hkl = numpy.array([[1, 0, 2, 1, 3], [0, 1, 1, 2, 2], [1, 1, 0, 3, 1]])
    f_nucl = crystal.calc_f_nucl(index_hkl)
    crystal.convert_to_superspace()
    assert crystal.is_attribute("space_group_symop_magn_ssg_operation")
    assert not crystal.is_attribute("space_group")
    assert len(crystal.space_group_symop_magn_ssg_operation.items) == 8
    assert numpy.allclose(crystal.calc_f_nucl(index_hkl), f_nucl, atol=1e-10)
    with pytest.raises(CrysPyError):
        crystal.convert_to_superspace()


def test_assistant_finds_convert_to_superspace_by_words():
    """Requests of several words find methods by their docstrings
    (the assistant asked 'symmetry operation' and found nothing)."""
    crystal = form_crystal("P n m a", (5.5, 7.6, 5.4, 90, 90, 90),
                           [("Mn", "Mn", 0.0, 0.0, 0.5)])
    for pattern in ("symmetry operation", "superspace", "convert to superspace",
                    "ssg", "Convert this crystal to superspace"):
        result = cryspy.call_assistant_tool(crystal, "list_methods", {"pattern": pattern})
        assert not result["flag_error"], pattern
        assert "convert_to_superspace" in result["text"], pattern


def test_3d_view_keeps_moments_after_conversion_to_superspace():
    """The 3D view draws the same moments at the same atoms before and
    after convert_to_superspace (k = 0)."""
    f_name = os.path.join(os.path.dirname(__file__), "..", "..", "..", "examples",
                          "rhochi_unpolarized_powder_1d_mcif_URu0.96Rh0.04Si2",
                          "rhochi_unpolarized_powder_1d_mcif_URu0.96Rh0.04Si2.rcif")
    crystal = [item for item in cryspy.load_file(f_name).items
               if isinstance(item, cryspy.Crystal)][0]
    settings = cryspy.StructureViewSettings(n_cells_a=2, n_cells_b=1, n_cells_c=2)

    def take_atoms(scene):
        moment = numpy.nan_to_num(scene["moment_xyz"])
        return sorted(map(tuple, numpy.round(numpy.vstack([scene["atom_fract_xyz"], moment]), 4).T))

    l_before = take_atoms(crystal.get_structure_scene(settings))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        crystal.convert_to_superspace()
    l_after = take_atoms(crystal.get_structure_scene(settings))
    assert l_after == l_before
    assert max(abs(value) for row in l_after for value in row[3:]) > 0.5
