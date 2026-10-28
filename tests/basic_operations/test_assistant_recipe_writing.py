"""Tools of an AI assistant writing a recipe from the words of the user."""
import pytest

import cryspy
from cryspy.procedure_assistant import recipes

RECIPE = """Recipe scale_then_adp: Scale factor, then thermal parameters
Requests: refine the scale and then the adp
For: single_crystal, unpolarized

1. start: Fix all parameters and calculate the start model.
   Actions: fix_parameters; calculate
2. scale: Scale factor.
   Actions: refine_parameters {"pattern": "scale"}; refine
3. adp: Isotropic thermal parameters.
   Actions: refine_parameters {"pattern": "adp"}; refine

Final answer: chi_sq/n at the start and at the end.
"""


@pytest.fixture
def user_folder(tmp_path):
    recipes.set_user_recipe_folder(str(tmp_path))
    yield tmp_path
    recipes.set_user_recipe_folder(None)


def call(name, arguments=None):
    return cryspy.call_assistant_tool(None, name, arguments)


def test_guide_has_format_kinds_examples_and_way_of_work():
    result = call("get_recipe_guide")
    assert not result["flag_error"]
    for text in ("Recipe <name>: <title>", "flip_ratio", "refine_structure", "ask_user",
                 "check_recipe_text", "save_recipe"):
        assert text in result["text"]


def test_check_correct_recipe_gives_steps():
    result = call("check_recipe_text", {"text": RECIPE})
    assert not result["flag_error"]
    assert "scale_then_adp is correct: 3 steps (start, scale, adp)" in result["text"]
    assert "Remarks" not in result["text"]


def test_check_wrong_recipe_gives_line_and_unknown_tool():
    result = call("check_recipe_text", {"text": "Recipe a: b\n1. one: x\nfoo"})
    assert result["flag_error"] and "line 3" in result["text"]
    text = RECIPE.replace("refine_parameters {\"pattern\": \"adp\"}", "refine_adp")
    result = call("check_recipe_text", {"text": text})
    assert result["flag_error"] and "unknown tool refine_adp" in result["text"]


def test_check_gives_remarks():
    text = RECIPE.replace("Actions: refine_parameters {\"pattern\": \"adp\"}; refine",
                          "Actions: refine_parameters {\"pattern\": \"adp\"}")
    text = text.replace("scale_then_adp", "refine_structure").replace("For: single_crystal, unpolarized\n", "")
    result = call("check_recipe_text", {"text": text})
    assert not result["flag_error"]
    assert "no step refines them" in result["text"]
    assert "recipe of CrysPy" in result["text"]
    assert "no 'For:'" in result["text"]


def test_ask_user_asks_as_missing_and_returns_answers():
    arguments = {"question": "Which parameters are refined after the scale",
                 "options": ["adp", "atomic positions"]}
    result = call("ask_user", arguments)
    assert not result["flag_error"]
    [d_missing] = result["missing"]
    assert d_missing["question"] == "Which parameters are refined after the scale?"
    assert "1: adp; 2: atomic positions" in d_missing["hint"]
    answer = cryspy.convert_assistant_answer(d_missing, "2")
    result = call("ask_user", dict(arguments, **{d_missing["name"]: answer}))
    assert result["missing"] == []
    assert "-> atomic positions" in result["text"]


def test_ask_user_several_questions():
    result = call("ask_user", {"questions": ["Which experiment?", "Undo when R_F2 is worse?"]})
    assert [d["name"] for d in result["missing"]] == ["answer_1", "answer_2"]
    result = call("ask_user", {"questions": ["Which experiment?", "Undo when R_F2 is worse?"],
                               "answer_1": "exp1", "answer_2": "yes"})
    assert "Which experiment? -> exp1" in result["text"]
    assert "Undo when R_F2 is worse? -> yes" in result["text"]
    assert call("ask_user", {})["flag_error"]


def test_save_recipe_and_ask_before_replacing(user_folder):
    result = call("save_recipe", {"text": RECIPE})
    assert not result["flag_error"], result["text"]
    path = user_folder / "scale_then_adp.txt"
    assert path.is_file()
    assert recipes.take_recipe("scale_then_adp") is not None
    new_text = RECIPE.replace("Scale factor, then", "The scale, then")
    result = call("save_recipe", {"text": new_text})
    [d_missing] = result["missing"]
    assert d_missing["name"] == "replace"
    assert "scale_then_adp" in d_missing["question"]
    no = cryspy.convert_assistant_answer(d_missing, "no")
    result = call("save_recipe", {"text": new_text, "replace": no})
    assert result["flag_error"]
    assert "Scale factor, then" in path.read_text()
    yes = cryspy.convert_assistant_answer(d_missing, "yes")
    result = call("save_recipe", {"text": new_text, "replace": yes})
    assert not result["flag_error"]
    assert "The scale, then" in path.read_text()


def test_save_wrong_recipe_or_without_folder_is_error(user_folder):
    result = call("save_recipe", {"text": "Recipe a: b"})
    assert result["flag_error"] and "not saved" in result["text"]
    recipes.set_user_recipe_folder(None)
    assert call("save_recipe", {"text": RECIPE})["flag_error"]


def test_recipe_step_may_ask_the_user():
    text = RECIPE.replace("Actions: fix_parameters; calculate",
                          "Actions: ask_user {\"question\": \"Start from the current values?\"}; calculate")
    assert not call("check_recipe_text", {"text": text})["flag_error"]
