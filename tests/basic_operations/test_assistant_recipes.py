"""Recipes and the history of refinements for AI assistants."""
import cryspy
from cryspy.procedure_assistant import history
from cryspy.procedure_assistant.parameters import take_parameter_record
from cryspy.procedure_assistant.recipes import L_RECIPES, take_recipe
from cryspy.procedure_assistant.tools import take_assistant_tool

RCIF = """global_

data_fe
_cell_length_a 5.0
_cell_length_b 5.0
_cell_length_c 5.0
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_space_group_name_H-M_alt "P m -3 m"

loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
_atom_site_adp_type
_atom_site_B_iso_or_equiv
  Fe  Fe3+  0.0  0.0  0.0  1.0  Biso  0.3
  O   O2-   0.5  0.0  0.0  1.0  Biso  0.5
"""


def test_recipe_actions_are_tools():
    for recipe in L_RECIPES:
        assert len(recipe.steps) > 0
        for step in recipe.steps:
            for name, arguments in step.actions:
                tool = take_assistant_tool(name)
                assert tool is not None, f"{recipe.name}.{step.name}: unknown tool {name}"
                assert set(arguments) <= set(tool.parameters)


def test_get_recipe_lists_steps():
    result = cryspy.call_assistant_tool(None, "get_recipe", {"name": "refine_structure"})
    assert not result["flag_error"]
    for step in take_recipe("refine_structure").steps:
        assert step.name in result["text"]
    assert "refine_structure" in cryspy.call_assistant_tool(None, "list_recipes")["text"]


def test_get_unknown_recipe_is_error_with_names():
    result = cryspy.call_assistant_tool(None, "get_recipe", {"name": "magic"})
    assert result["flag_error"]
    assert "refine_structure" in result["text"]


def test_run_unknown_recipe_raises():
    obj = cryspy.str_to_globaln(RCIF)
    try:
        cryspy.run_recipe(obj, "magic")
    except ValueError as error:
        assert "refine_structure" in str(error)
    else:
        raise AssertionError("ValueError is expected")


def test_snapshot_restores_values_sigmas_and_flags():
    obj = cryspy.str_to_globaln(RCIF)
    path = "crystal_fe.atom_site[Fe].b_iso_or_equiv"
    record = take_parameter_record(obj, path)
    history.keep_checkpoint(obj, 10.)
    record.set_refinement(True)  # marked after the checkpoint (refine_parameters)
    history.keep_snapshot(obj)
    record.set_value(-1.)  # as a refinement would do
    setattr(record.item, "b_iso_or_equiv_sigma", 0.5)
    history.keep_checkpoint(obj, 20.)
    snapshot = history.restore_snapshot(obj)
    assert snapshot is not None
    assert record.value == 0.3
    assert record.sigma == 0.
    assert not record.flag_refinement
    assert history.take_history(obj).chi_sq == 10.
    assert history.restore_snapshot(obj) is None


def test_undo_without_refinement_is_error():
    obj = cryspy.str_to_globaln(RCIF)
    result = cryspy.call_assistant_tool(obj, "undo_refinement")
    assert result["flag_error"]


def test_warnings_negative_adp_and_worse_chi_sq():
    obj = cryspy.str_to_globaln(RCIF)
    take_parameter_record(obj, "crystal_fe.atom_site[O].b_iso_or_equiv").set_value(-0.2)
    d_metrics = history.form_metrics({"chi_sq": 110., "n_point": 10}, obj, chi_sq_previous=100.)
    assert d_metrics["chi_sq_per_n"] == 11.
    assert abs(d_metrics["chi_sq_change"] - 0.1) < 1e-12
    s_warning = " ".join(d_metrics["warnings"])
    assert "worse" in s_warning
    assert "atom_site[O].b_iso_or_equiv" in s_warning
    assert "atom_site[Fe]" not in s_warning
    text = history.format_metrics(d_metrics)
    assert "worse" in text and "undo_refinement" in text


def test_no_warnings_for_better_chi_sq():
    obj = cryspy.str_to_globaln(RCIF)
    d_metrics = history.form_metrics({"chi_sq": 90.}, obj, chi_sq_previous=100.)
    assert d_metrics["warnings"] == []
    assert "better" in history.format_metrics(d_metrics)


def test_tools_without_numbers_have_empty_metrics():
    obj = cryspy.str_to_globaln(RCIF)
    result = cryspy.call_assistant_tool(obj, "summarize_object")
    assert result["metrics"] == {}
