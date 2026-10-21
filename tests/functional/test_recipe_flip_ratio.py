"""Recipes for flip ratios of a paramagnetic crystal with a nuclear space group (Ho2Ti2O7)."""
import os

import cryspy

F_HOTI = os.path.join(os.path.dirname(__file__), "HoTi_single_test", "main.rcif")


def load_without_magnetic_atoms():
    obj = cryspy.load_file(F_HOTI)
    crystal = obj.crystal_ho2ti2o7
    crystal.items.remove(crystal.atom_site_susceptibility)
    crystal.items.remove(crystal.atom_site_scat)
    return obj


def test_suggest_recipes_for_flip_ratios():
    obj = cryspy.load_file(F_HOTI)
    text = cryspy.call_assistant_tool(obj, "suggest_recipes")["text"]
    assert "expHo2Ti2O7 (single_crystal, polarized, flip_ratio, nuclear): refine_flip_ratio_nuclear" in text
    assert "Run: refine_flip_ratio_nuclear for expHo2Ti2O7" in text


def test_check_asks_for_field_and_magnetic_atoms():
    obj = load_without_magnetic_atoms()
    obj.diffrn_expho2ti2o7.setup.field = 0.
    result = cryspy.call_assistant_tool(obj, "check_flip_ratio_model")
    assert [d_missing["name"] for d_missing in result["missing"]] == ["field", "magnetic_atoms", "lande"]
    assert "magnetic field of expHo2Ti2O7 (setup.field) is zero" in result["text"]
    d_atoms = result["missing"][1]
    assert cryspy.convert_assistant_answer(d_atoms, "Ho1") == ["Ho1"]
    result = cryspy.call_assistant_tool(
        obj, "check_flip_ratio_model", {"field": 1., "magnetic_atoms": ["Ho1"], "lande": 1.25})
    assert not result["flag_error"], result["text"]
    assert result["missing"] == []
    assert "the model gives magnetic scattering" in result["text"]
    assert obj.diffrn_expho2ti2o7.setup.field == 1.
    assert obj.crystal_ho2ti2o7.atom_site_scat.items[0].lande == 1.25


def test_recipe_refines_susceptibility_from_scratch():
    """A new magnetic atom: the start value of chi is scanned, then Ciso,
    Cani, kappa, polarization, positions and extinction are refined."""
    obj = load_without_magnetic_atoms()
    cryspy.call_assistant_tool(obj, "check_flip_ratio_model", {"magnetic_atoms": ["Ho1"], "lande": 1.25})
    result = cryspy.run_recipe(obj, "refine_flip_ratio_nuclear")
    d_status = {name: status for name, status, _ in result["steps"]}
    assert d_status["start"] == "done"
    assert d_status["chi_iso"] == "done"
    assert d_status["chi_aniso"] == "done"
    item = obj.crystal_ho2ti2o7.atom_site_susceptibility.items[0]
    assert item.chi_type == "Cani"
    # the values of the refined model in main.rcif: chi_11 = 3.329(29), chi_12 = 3.254(26)
    assert 2.5 < item.chi_11 < 4.5 and 2.5 < item.chi_12 < 4.5
    assert "at the end" in result["text"]
