"""Methods of the data blocks, loops and items for AI assistants: found in
the class definitions, described from docstrings and type hints."""
import os
import warnings

import pytest

import cryspy
from cryspy.procedure_assistant import recipes
from cryspy.procedure_assistant.methods import take_methods
from cryspy.procedure_assistant.procedures import parse_docstring

DIR_EXAMPLE = os.path.join(os.path.dirname(__file__), "..", "..", "examples",
                           "rhochi_polarized_powder_1d_Ho2Ti2O7_5K_1T")
FILE_RCIF = os.path.join(DIR_EXAMPLE, "rhochi_polarized_powder_1d_Ho2Ti2O7_5K_1T.rcif")


@pytest.fixture
def powder():
    return cryspy.load_file(FILE_RCIF)


def call(obj, name, arguments=None):
    return cryspy.call_assistant_tool(obj, name, arguments)


def test_methods_from_classes_without_internal_and_deprecated():
    d_method = take_methods(cryspy.Pd)
    assert "define_background_points" in d_method and "estimate_background" in d_method
    assert not any(name.startswith(("plot", "get_", "_", "to_")) for name in d_method)
    assert "add_items" not in d_method  # a method of DataN
    assert "calc_asymmetric_cell" not in take_methods(cryspy.SpaceGroup)  # deprecated
    assert "calc_sthovl" in take_methods(cryspy.Cell)


def test_list_methods_short_and_with_pattern(powder):
    result = call(powder, "list_methods")
    assert not result["flag_error"]
    assert "pd_powder1 (Pd): define_background_points, estimate_background" in result["text"]
    assert "pd_powder1.pd_background (PdBackgroundL): define_points" in result["text"]
    result = call(powder, "list_methods", {"pattern": "background"})
    assert ("define_background_points(step_ttheta: number=1.0): Define the background points "
            "from the measured profile.") in result["text"]
    assert "define_points(pd_meas: PdMeasL, step_ttheta: number=10.0)" in result["text"]
    assert call(powder, "list_methods", {"path": "nothing_like_this"})["flag_error"]


def test_call_method_of_data_block_by_its_class(powder):
    with warnings.catch_warnings():
        warnings.simplefilter("error", DeprecationWarning)
        result = call(powder, "call_method", {"method": "define_background_points", "step_ttheta": 3})
    assert not result["flag_error"] and result["flag_object_changed"]
    assert result["text"] == "pd_powder1.define_background_points finished."
    n_point = len(powder.pd_powder1.pd_meas.items)
    assert len(powder.pd_powder1.pd_background.items) == n_point


def test_call_method_finds_items_of_the_object(powder):
    # pd_meas of define_points is taken from the same data block
    result = call(powder, "call_method", {"path": "pd_powder1.pd_background", "method": "define_points",
                                          "arguments": {"step_ttheta": 5}})
    assert not result["flag_error"], result["text"]
    assert len(powder.pd_powder1.pd_background.items) == len(powder.pd_powder1.pd_meas.items)


def test_call_method_asks_values_without_default():
    cell = cryspy.Cell(length_a=2., length_b=2., length_c=2.,
                       angle_alpha=90., angle_beta=90., angle_gamma=90.)
    result = call(cell, "call_method", {"method": "calc_sthovl", "index_h": 1})
    assert [d_missing["name"] for d_missing in result["missing"]] == ["index_k", "index_l"]
    result = call(cell, "call_method", {"method": "calc_sthovl", "index_h": 1, "index_k": 0, "index_l": 0})
    assert not result["flag_error"] and "0.25" in result["text"]


def test_call_unknown_method_names_where_it_is(powder):
    result = call(powder, "call_method", {"path": "pd_powder1.pd_background",
                                          "method": "define_background_points"})
    assert result["flag_error"] and "it is a method of pd_powder1" in result["text"]


def test_docstring_with_several_parameters_in_one_line():
    d_doc = parse_docstring("""Save.

    Parameters
    ----------
    Na, Nb, Nc : int, optional
        Numbers of points along the axes a, b, c.
    """)
    assert d_doc["parameters"]["Nb"]["description"] == "Numbers of points along the axes a, b, c."


def test_recipe_step_name_in_several_words_and_method_as_action():
    text = """Recipe polarized_powder: Polarized powder 1d
Requests: refine polarized powder data
For: powder, polarized
1. Define background points: describe the background
   Actions: call_method {"method": "define_background_points", "step_ttheta": 2}
2. scale: Scale factor.
   Actions: refine_parameters {"pattern": "scale"}; refine
"""
    result = recipes.check_recipe_text(text)
    assert "is correct: 2 steps (define_background_points, scale)" in result
    result = recipes.check_recipe_text(text.replace(
        'call_method {"method": "define_background_points", "step_ttheta": 2}', "define_background_points"))
    assert "it is a method of Pd: write call_method" in result
    result = recipes.check_recipe_text(text.replace("define_background_points\", ", "no_such_method\", "))
    assert "no CrysPy class has the method 'no_such_method'" in result
