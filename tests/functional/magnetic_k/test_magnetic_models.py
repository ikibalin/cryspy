"""Magnetic models of a propagation vector without data
(Crystal.report_magnetic_models) and the choice of one by its number
(Crystal.add_magnetic_model)."""
import pytest

import cryspy
from cryspy.exceptions import CrysPyError
from test_rank_magnetic_representations import CRYSTAL, K_VECTOR, TBTE3, take_rows


def take_models(report) -> dict:
    """{name: (model, free Fourier parameters, domains)}"""
    return {row[2].strip(): (int(row[1]), row[4].strip(), row[5].strip())
            for row in take_rows(report)}


def test_models_of_p21m():
    """P 1 21/m 1, k = (1/2, 0, 0): every operation keeps k, so a
    one-dimensional representation has one domain. A sum Γi⊕Γj keeps the
    operations where the characters are equal; the inversion keeps Friedel
    pairs together, so only the sums that lose the inversion (Γ1⊕Γ3,
    Γ2⊕Γ4) have a second domain, the 2_1 axis."""
    crystal = cryspy.str_to_globaln(CRYSTAL).crystal_phase
    d_model = take_models(crystal.report_magnetic_models(*K_VECTOR, "Mn1"))
    assert [d_model[f"Γ{i}"][0] for i in (1, 2, 3, 4)] == [1, 2, 3, 4]
    assert len(d_model) == 10
    # Γ1 allows the moment along b only (cos and sin), Γ2 in the ac plane
    assert d_model["Γ1"][1] == "Mn1 2" and d_model["Γ2"][1] == "Mn1 4"
    assert d_model["Γ1"][2] == "1"
    assert d_model["Γ1⊕Γ3"][2] == "2: -x,y+1/2,-z,+1"
    assert d_model["Γ2⊕Γ4"][2] == "2: -x,y+1/2,-z,+1"
    assert d_model["Γ1⊕Γ2"][2] == "1"
    # the crystal is not changed
    assert crystal.is_attribute("space_group")
    assert not crystal.is_attribute("atom_site_moment_fourier")
    # without the superpositions the numbers of the others are the same
    d_single = take_models(crystal.report_magnetic_models(
        *K_VECTOR, "Mn1", try_superpositions=False))
    assert {name: value[0] for name, value in d_single.items()} == {
        f"Γ{i}": i for i in (1, 2, 3, 4)}


def test_add_magnetic_model_is_add_propagation_vector():
    crystal = cryspy.str_to_globaln(CRYSTAL).crystal_phase
    model = take_models(crystal.report_magnetic_models(*K_VECTOR, "Mn1"))["Γ1⊕Γ3"][0]
    report = crystal.add_magnetic_model(*K_VECTOR, "Mn1", model, 2.0)
    assert report.startswith(f"Model {model} (Γ1⊕Γ3)")
    reference = cryspy.str_to_globaln(CRYSTAL).crystal_phase
    reference.add_propagation_vector(*K_VECTOR, "Mn1", "1,3", 2.0)
    assert crystal.to_cif() == reference.to_cif()
    with pytest.raises(CrysPyError, match="model must be from 1 to 10"):
        cryspy.str_to_globaln(CRYSTAL).crystal_phase.add_magnetic_model(*K_VECTOR, "Mn1", 11)


def test_models_without_form_factor():
    """TbTe3 after a nuclear refinement (no atom_site_scat): the models are
    listed all the same; the two domains are the two arms of the star of
    k = (1/2, 1/2, 0) in C m c m."""
    crystal = cryspy.str_to_globaln(TBTE3).crystal_tbte3
    d_model = take_models(crystal.report_magnetic_models(0.5, 0.5, 0., "Tb"))
    assert all(value[2] == "2: -x,y,-z+1/2,+1" for name, value in d_model.items()
               if "⊕" not in name)
    assert not crystal.is_attribute("atom_site_scat")
