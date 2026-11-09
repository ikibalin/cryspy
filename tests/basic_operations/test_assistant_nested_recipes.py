"""Recipes built from recipes: a step with use_recipe does the steps of
another recipe of the CookBook."""
import pytest

import cryspy
from cryspy.procedure_assistant import recipes

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
"""

BRICK = """Recipe brick: Fix and list
Requests: fix all
For: any
1. fix: Fix all parameters.
   Actions: fix_parameters
2. look: List the parameters.
   Actions: list_parameters
"""

BLOCK = """Recipe block: The brick, then the cell
Requests: fix all and look at the cell
For: any
1. first: The brick.
   Actions: use_recipe {"name": "brick"}
2. cell: Show the cell.
   Actions: show_item {"path": "fe.cell"}
"""


@pytest.fixture
def user_folder(tmp_path):
    recipes.set_user_recipe_folder(str(tmp_path))
    yield tmp_path
    recipes.set_user_recipe_folder(None)


def test_recipe_built_from_a_recipe_runs_its_steps(user_folder):
    recipes.save_user_recipe(BRICK)
    recipes.save_user_recipe(BLOCK)
    obj = cryspy.str_to_globaln(RCIF)
    l_called = []
    result = cryspy.run_recipe(obj, "block", callback=lambda *args: l_called.append(args[:2]))
    assert l_called == [("first/fix", "done"), ("first/look", "done"), ("first", "done"), ("cell", "done")]
    assert "- first/look: done" in result["text"]


def test_use_recipe_gives_the_steps_to_the_model(user_folder):
    recipes.save_user_recipe(BRICK)
    result = cryspy.call_assistant_tool(None, "use_recipe", {"name": "brick"})
    assert not result["flag_error"]
    assert "next step of the recipe that uses it" in result["text"] and "1. fix: Fix all" in result["text"]
    assert cryspy.call_assistant_tool(None, "use_recipe", {"name": "nothing"})["flag_error"]


def test_unknown_used_recipe_and_circle_are_errors(user_folder):
    text = recipes.check_recipe_text(BLOCK)
    assert "use_recipe: recipe 'brick' is not in the CookBook" in text
    recipes.save_user_recipe(BRICK)
    assert "is correct: 2 steps" in recipes.check_recipe_text(BLOCK)
    # brick now uses block, which uses brick
    circle = BRICK.replace("Actions: list_parameters", 'Actions: use_recipe {"name": "block"}')
    recipes.save_user_recipe(BLOCK)
    text = recipes.check_recipe_text(circle)
    assert "recipes use each other in a circle: brick -> block -> brick" in text
    itself = BRICK.replace("Actions: list_parameters", 'Actions: use_recipe {"name": "brick"}')
    assert "circle: brick -> brick" in recipes.check_recipe_text(itself)


def test_guide_tells_about_use_recipe():
    assert 'use_recipe {"name"' in recipes.form_recipe_guide()
