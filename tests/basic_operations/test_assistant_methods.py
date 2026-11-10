"""Methods of the data blocks, loops and items for AI assistants: found in
the class definitions, described from docstrings and type hints."""
import os
import warnings

import pytest

import cryspy
from cryspy.procedure_assistant import recipes
from cryspy.procedure_assistant.methods import take_methods
from cryspy.procedure_assistant.procedures import parse_docstring
from cryspy.A_functions_base.function_1_deprecation import deprecated

DIR_EXAMPLE = os.path.join(os.path.dirname(__file__), "..", "..", "examples",
                           "rhochi_polarized_powder_1d_Ho2Ti2O7_5K_1T")
FILE_RCIF = os.path.join(DIR_EXAMPLE, "rhochi_polarized_powder_1d_Ho2Ti2O7_5K_1T.rcif")


@pytest.fixture
def powder():
    return cryspy.load_file(FILE_RCIF)


def call(obj, name, arguments=None):
    return cryspy.call_assistant_tool(obj, name, arguments)


def test_methods_from_classes_without_internal_and_deprecated(monkeypatch):
    d_method = take_methods(cryspy.Pd)
    assert "define_background_points" in d_method and "estimate_background" in d_method
    assert not any(name.startswith(("plot", "get_", "_", "to_")) for name in d_method)
    assert "add_items" not in d_method  # a method of DataN
    def calc_twice(self, h: int = 1):
        """Deprecated method: not offered to assistants."""
        return 2 * h

    monkeypatch.setattr(cryspy.Cell, "calc_twice", deprecated(calc_twice), raising=False)
    assert "calc_twice" not in take_methods(cryspy.Cell)
    assert "calc_sthovl" in take_methods(cryspy.Cell)


def test_list_methods_short_and_with_pattern(powder):
    result = call(powder, "list_methods")
    assert not result["flag_error"]
    # all methods, one line each: the model chooses by their meaning
    assert ("pd_powder1 (Pd):\n  define_background_points: Define the background points "
            "from the measured profile.") in result["text"]
    assert "pd_powder1.pd_background (PdBackgroundL):\n  define_points: " in result["text"]
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


def test_list_methods_never_empty_and_asks_when_unsure(powder):
    """Without a matching word all methods are still listed (the model
    chooses by their descriptions); equal candidates suggest ask_user."""
    text = call(powder, "list_methods", {"pattern": "zebra giraffe"})["text"]
    assert "No method has the words" in text and "ask_user" in text
    assert "define_background_points: Define the background points" in text
    text = call(powder, "list_methods", {"pattern": "ellipsoids of zebra"})["text"]
    assert text.startswith("Methods with only some words of the request")
    text = call(powder, "list_methods", {"pattern": "atoms ellipsoids"})["text"]
    assert text.startswith("Best matching methods:")
    assert "report_ellipsoids(" in text and "Other methods:" in text


def test_methods_matching_a_request_are_found(powder):
    """A small model answered "cannot convert to superspace" without
    looking at the methods; they are given with the request."""
    from cryspy.procedure_assistant.methods import form_matching_methods_text
    text = form_matching_methods_text(powder, "Define the background points")
    assert "define_background_points" in text
    assert "call_method" in text
    crystal = cryspy.str_to_globaln("""global_
data_mn
_space_group_name_H-M_alt 'P n m a'
_cell_length_a 5.5
_cell_length_b 7.6
_cell_length_c 5.4
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Mn Mn 0 0 0.5
""")
    text = form_matching_methods_text(crystal, "Convert the structure to superspace")
    assert text.split("\n")[1].startswith("- convert_to_superspace of crystal_mn")
    assert form_matching_methods_text(crystal, "hello") == ""
    assert form_matching_methods_text(crystal, "what is the chi square") == ""
