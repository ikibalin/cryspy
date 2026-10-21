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
    history.keep_checkpoint(obj, {"chi_sq_per_n": 10.})
    record.set_refinement(True)  # marked after the checkpoint (refine_parameters)
    history.keep_snapshot(obj)
    record.set_value(-1.)  # as a refinement would do
    setattr(record.item, "b_iso_or_equiv_sigma", 0.5)
    history.keep_checkpoint(obj, {"chi_sq_per_n": 20.})
    snapshot = history.restore_snapshot(obj)
    assert snapshot is not None
    assert record.value == 0.3
    assert record.sigma == 0.
    assert not record.flag_refinement
    assert history.take_history(obj).d_metrics == {"chi_sq_per_n": 10.}
    assert history.restore_snapshot(obj) is None


def test_undo_without_refinement_is_error():
    obj = cryspy.str_to_globaln(RCIF)
    result = cryspy.call_assistant_tool(obj, "undo_refinement")
    assert result["flag_error"]


def test_chi_sq_per_point_is_compared():
    """The number of points may change: chi_sq/n is compared, not chi_sq."""
    obj = cryspy.str_to_globaln(RCIF)
    d_metrics = history.form_metrics({"chi_sq": 110., "n_point": 10}, obj,
                                     d_previous={"chi_sq_per_n": 10., "n_point": 5})
    assert d_metrics["chi_sq_per_n"] == 11.
    assert abs(d_metrics["chi_sq_per_n_change"] - 0.1) < 1e-12
    assert "chi_sq/n got worse" in " ".join(d_metrics["warnings"])
    text = history.format_metrics(d_metrics)
    assert "10 -> 11" in text and "worse" in text and "undo_refinement" in text
    # more points, a larger chi_sq, the same fit
    d_metrics = history.form_metrics({"chi_sq": 200.}, obj,
                                     d_previous={"chi_sq_per_n": 10., "n_point": 20})
    assert d_metrics["warnings"] == []
    assert "unchanged" in history.format_metrics(d_metrics)


def test_warnings_negative_adp_of_refined_parameters_only():
    obj = cryspy.str_to_globaln(RCIF)
    record_o = take_parameter_record(obj, "crystal_fe.atom_site[O].b_iso_or_equiv")
    record_o.set_value(-0.2)
    assert history.form_metrics({"chi_sq": 90.}, obj)["warnings"] == []
    record_o.set_refinement(True)
    s_warning = " ".join(history.form_metrics({"chi_sq": 90.}, obj)["warnings"])
    assert "atom_site[O].b_iso_or_equiv" in s_warning
    assert "atom_site[Fe]" not in s_warning


class Reflection:
    def __init__(self, intensity, intensity_calc, excluded=False):
        self.intensity = intensity
        self.intensity_calc = intensity_calc
        self.excluded = excluded

    def is_attribute(self, name):
        return True


class Reflections:
    def __init__(self, items):
        self.items = items


def test_r_f2():
    """R_F2 = sum |I_obs - I_calc| / sum I_obs without excluded reflections."""
    loop = Reflections([Reflection(10., 8.), Reflection(30., 33.), Reflection(5., 50., excluded=True)])
    assert abs(history.calc_r_f2(loop) - 5. / 40.) < 1e-12
    assert history.calc_r_f2(Reflections([])) is None


def test_r_f2_change_is_a_warning():
    obj = cryspy.str_to_globaln(RCIF)
    d_metrics = {"chi_sq": 1., "r_f2": {"exp": 0.12}, "r_f2_previous": {"exp": 0.10}}
    s_warning = " ".join(history.find_warnings(d_metrics, obj))
    assert "R_F2 of exp got worse" in s_warning
    text = history.format_metrics(dict(d_metrics, warnings=[]))
    assert "R_F2 of exp: 0.1000 -> 0.1200 (+20 %, worse)" in text


def test_tools_without_numbers_have_empty_metrics():
    obj = cryspy.str_to_globaln(RCIF)
    result = cryspy.call_assistant_tool(obj, "summarize_object")
    assert result["metrics"] == {}
