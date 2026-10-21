"""Recipes as text files: the format, the folder of the user, the kinds of experiments."""
import os

import pytest

import cryspy
from cryspy.procedure_assistant import history, recipes
from cryspy.procedure_assistant.parameters import take_parameter_record

RCIF_MAGNETIC_ATOM = """global_

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

loop_
_atom_site_susceptibility_label
_atom_site_susceptibility_chi_type
_atom_site_susceptibility_chi_11
_atom_site_susceptibility_chi_22
_atom_site_susceptibility_chi_33
_atom_site_susceptibility_chi_12
_atom_site_susceptibility_chi_13
_atom_site_susceptibility_chi_23
  Fe  Cani  1.0  1.0  1.0  0.0  0.0  0.0

loop_
_atom_site_scat_label
_atom_site_scat_lande
_atom_site_scat_kappa
  Fe  2.0  1.0
"""


@pytest.fixture
def user_folder(tmp_path):
    recipes.set_user_recipe_folder(str(tmp_path))
    yield str(tmp_path)
    recipes.set_user_recipe_folder(None)


def test_template_is_a_correct_recipe():
    recipe = recipes.parse_recipe(recipes.S_TEMPLATE)
    assert recipe.name == "my_recipe"
    assert recipe.kinds == ("single_crystal", "unpolarized")
    assert [step.name for step in recipe.steps] == ["start", "scale"]
    assert recipe.steps[0].actions == (("fix_parameters", {}), ("calculate", {}))
    assert recipe.steps[1].actions == (("refine_parameters", {"pattern": "scale"}), ("refine", {}))
    assert recipe.steps[1].on_fail == recipes.S_UNDO  # the default of a step without 'Otherwise'
    assert recipes.check_recipe(recipe) == []
    assert "2. scale: Scale factor." in recipe.to_text()


def test_continuation_lines_and_comments():
    recipe = recipes.parse_recipe("""Recipe two_lines: Title
# a comment
1. first: A goal written
   on two lines.
   Actions: calculate
""")
    assert recipe.steps[0].goal == "A goal written on two lines."


@pytest.mark.parametrize("text, message", [
    ("", "empty"),
    ("1. start: goal", "line 1: the recipe starts with"),
    ("Recipe a: b\n", "has no steps"),
    ("Recipe a: b\n1. s: g\n   Actions: refine {pattern}\n", "line 2: step s: the arguments of refine are not JSON"),
    ("Recipe a: b\n1. s: g\nActions: calculate\n", "line 3: unexpected line"),
    ("Recipe a: b\n   Actions: calculate\n", "line 2: 'Actions' outside of a step"),
])
def test_wrong_recipe_has_the_number_of_the_line(text, message):
    with pytest.raises(ValueError, match=message):
        recipes.parse_recipe(text)


def test_check_recipe_finds_unknown_tools_arguments_and_kinds():
    recipe = recipes.parse_recipe("""Recipe a: b
For: single_crystal, liquid
1. s: g
   Actions: magic; refine_parameters {"patern": "scale"}
""")
    s_problem = "; ".join(recipes.check_recipe(recipe))
    assert "unknown kinds in 'For': liquid" in s_problem
    assert "unknown tool magic" in s_problem
    assert "refine_parameters has no arguments patern" in s_problem


def test_user_recipe_replaces_the_one_of_cryspy(user_folder):
    text = recipes.S_TEMPLATE.replace("my_recipe", "refine_structure")
    path = recipes.save_user_recipe(text)
    assert path == os.path.join(user_folder, "refine_structure.txt")
    recipe = recipes.take_recipe("refine_structure")
    assert recipe.flag_user()
    assert [step.name for step in recipe.steps] == ["start", "scale"]
    d_file = {(d_out["name"], d_out["flag_user"]): d_out for d_out in recipes.list_recipe_files()}
    assert d_file[("refine_structure", False)]["flag_replaced"]
    assert not d_file[("refine_structure", True)]["flag_replaced"]
    assert recipes.delete_user_recipe("refine_structure")
    assert not recipes.take_recipe("refine_structure").flag_user()
    assert not recipes.delete_user_recipe("refine_structure")


def test_wrong_user_file_is_listed_with_its_error_and_not_used(user_folder):
    with open(os.path.join(user_folder, "broken.txt"), "w") as fid:
        fid.write("Recipe broken: b\n1. s: g\n   Actions: magic\n")
    d_file = {d_out["name"]: d_out for d_out in recipes.list_recipe_files()}
    assert "unknown tool magic" in d_file["broken"]["error"]
    assert recipes.take_recipe("broken") is None
    with pytest.raises(ValueError, match="unknown tool magic"):
        recipes.save_user_recipe("Recipe broken: b\n1. s: g\n   Actions: magic\n")


def test_save_without_user_folder_is_error():
    recipes.set_user_recipe_folder(None)
    with pytest.raises(ValueError, match="folder"):
        recipes.save_user_recipe(recipes.S_TEMPLATE)


def test_suggest_by_kinds_most_specific_first():
    l_name = [recipe.name for recipe in recipes.suggest_recipes(
        ["single_crystal", "polarized", "flip_ratio", "nuclear"])]
    assert l_name[0] == "refine_flip_ratio_nuclear"
    assert "refine_structure" not in l_name  # it is for unpolarized data
    assert "choose_recipes" not in l_name  # 'any': never suggested for an experiment
    l_name = [recipe.name for recipe in recipes.suggest_recipes(["powder", "powder_1d", "unpolarized", "nuclear"])]
    assert l_name == ["refine_structure"]


def test_set_chi_type_and_kappa_warning():
    obj = cryspy.str_to_globaln(RCIF_MAGNETIC_ATOM)
    result = cryspy.call_assistant_tool(obj, "set_chi_type", {"chi_type": "Ciso"})
    assert "Ciso is set for:\nfe: Fe" in result["text"]
    assert obj.crystal_fe.atom_site_susceptibility.items[0].chi_type == "Ciso"
    assert "already" in cryspy.call_assistant_tool(obj, "set_chi_type", {"chi_type": "Ciso"})["text"]
    record = take_parameter_record(obj, "crystal_fe.atom_site_scat[Fe].kappa")
    record.set_value(0.2)
    assert history.find_warnings({}, obj) == []  # not refined
    record.set_refinement(True)
    assert "kappa outside" in " ".join(history.find_warnings({}, obj))
