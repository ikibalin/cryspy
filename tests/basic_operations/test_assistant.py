"""Help pages and assistant tools (cryspy.procedure_assistant)."""
import json
import os
import re

import numpy
import pytest

import cryspy
from cryspy.procedure_assistant import help as cryspy_help
from cryspy.procedure_assistant.parameters import (
    find_parameter_records, take_parameter_record)

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
  Fe  Fe3+  0.0  0.0  0.0  1.0  Biso  0.3()
  O   O2-   0.5  0.0  0.0  1.0  Biso  0.5
"""


def take_object():
    return cryspy.str_to_globaln(RCIF)


def take_all_cif_tags() -> set:
    from cryspy.C_item_loop_classes.cl_1_atom_rho_multipole import \
        AtomRhoMultipole
    set_tag = set()
    for item_class in cryspy.L_ITEM_CLASS + [AtomRhoMultipole]:
        for name in item_class.ATTR_CIF:
            set_tag.add(f"_{item_class.PREFIX}_{name}".lower())
    return set_tag


def test_guides_use_existing_cif_tags():
    """A guide must not invent rcif tags."""
    set_tag = take_all_cif_tags()
    n_tag = 0
    for page in cryspy_help.load_guide_pages():
        l_tag = re.findall(r"`(_[A-Za-z][\w\-\[\]/]*)", page.text)
        # tags in the rcif examples
        l_tag += [tag for tag in re.findall(
            r"^\s*(_[A-Za-z][\w\-\[\]/]*)", page.text, flags=re.MULTILINE)]
        n_tag += len(l_tag)
        for tag in l_tag:
            assert tag.lower() in set_tag, f"{page.name}: unknown tag {tag}"
    assert n_tag > 30


def test_guides_have_titles_and_keywords():
    l_page = cryspy_help.load_guide_pages()
    assert len(l_page) >= 10
    for page in l_page:
        assert page.title != page.name
        assert len(page.keywords) > 0
        assert not page.text.startswith("---")


@pytest.mark.parametrize("query, page_name", [
    ("how to describe thermal vibrations", "atomic_displacement"),
    ("что необходимо использовать для описания тепловых колебаний",
     "atomic_displacement"),
    ("как описать мультиполи", "multipoles"),
    ("multipole model of form factor", "multipoles"),
    ("twins", "twinning"),
    ("background of powder pattern", "powder"),
    ("atom_site_aniso", "atom_site_aniso"),
])
def test_search_help(query, page_name):
    l_page = cryspy_help.search_help(query)
    assert page_name in [page.name for page in l_page[:2]]


def test_get_help_by_name_and_unknown():
    assert "atom_rho_multipole" in cryspy.get_help("multipoles")
    assert "_atom_site_aniso_U_11" in cryspy.get_help("atom_site_aniso")
    assert "Nothing is found" in cryspy.get_help("zzzz qqqq")


def test_category_pages_for_all_classes():
    l_name = [page.name for page in cryspy_help.load_category_pages()]
    for item_class in cryspy.L_ITEM_CLASS:
        assert item_class.PREFIX.lower() in l_name
    assert "atom_rho_multipole" in l_name
    assert "crystal" in l_name


def test_parameter_paths():
    obj = take_object()
    record = take_parameter_record(obj, "crystal_fe.atom_site[Fe].b_iso_or_equiv")
    assert record is not None
    assert record.flag_refinement
    assert record.value == pytest.approx(0.3)
    # brackets are literal, '*' is a wildcard
    assert len(find_parameter_records(obj, "atom_site[O]")) == 5
    assert len(find_parameter_records(obj, "*.atom_site[*].fract_x")) == 2
    assert len(find_parameter_records(obj, "", flag_refined_only=True)) == 1


def test_tools_change_object():
    obj = take_object()
    call = cryspy.call_assistant_tool
    result = call(obj, "refine_parameters", {"pattern": "atom_site[O].fract"})
    assert result["flag_object_changed"]
    # O on 3d site (1/2, 0, 0) of Pm-3m: all coordinates are fixed
    assert "Not refinable" in result["text"]
    call(obj, "set_parameter", {"path": "cell.length_a", "value": "5.1"})
    assert obj.crystal_fe.cell.length_a == pytest.approx(5.1)
    call(obj, "fix_parameters", {})
    assert len(obj.get_variable_names()) == 0
    result = call(obj, "add_atom", {
        "crystal": "fe", "label": "Fe2", "type_symbol": "Fe2+",
        "fract_x": 0.5, "fract_y": 0.5, "fract_z": 0.5})
    assert "added" in result["text"]
    assert len(obj.crystal_fe.atom_site.items) == 3
    obj_2 = cryspy.str_to_globaln(obj.to_cif())
    assert len(obj_2.crystal_fe.atom_site.items) == 3


def test_tool_errors_are_text():
    obj = take_object()
    call = cryspy.call_assistant_tool
    assert call(obj, "unknown", {})["text"].startswith("Error")
    assert "needs arguments" in call(obj, "set_parameter", {"path": "a"})["text"]
    assert "Error" in call(obj, "set_parameter",
                           {"path": "cell.length_a", "value": "abc"})["text"]
    assert "no CrysPy object" in call(None, "list_parameters", {})["text"]
    assert "Atomic displacement" in call(None, "get_help",
                                         {"query": "thermal vibrations"})["text"]


def test_tool_descriptions_are_json():
    l_tool = cryspy.assistant_tools()
    json.dumps(l_tool)
    l_name = [tool["function"]["name"] for tool in l_tool]
    assert "get_help" in l_name and "refine" in l_name
    assert "refine_parameters(pattern: string)" in cryspy.assistant_tools_as_text()


def test_calculate_powder_example():
    dir_cryspy = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
    f_name = os.path.join(
        dir_cryspy, "examples", "rhochi_unpolarized_powder_1d_PbSO4",
        "rhochi_unpolarized_powder_1d_PbSO4.rcif")
    obj = cryspy.load_file(f_name)
    result = cryspy.call_assistant_tool(obj, "calculate")
    assert "chi_sq" in result["text"]


@pytest.mark.parametrize("pattern, n_record", [
    ("atom_site.fract_x", 2),  # labels left out, as models write it
    ("atomic positions", 6),
    ("fract_x, fract_y", 4),
    ("b_iso | cell.length_a", 3),
    ("unit cell", 6),
])
def test_patterns_of_requests(pattern, n_record):
    assert len(find_parameter_records(take_object(), pattern)) == n_record


def test_no_match_is_error_with_names():
    obj = take_object()
    result = cryspy.call_assistant_tool(
        obj, "refine_parameters", {"pattern": "scale_factor"})
    assert result["flag_error"]
    assert result["text"].startswith("Error: no parameters match")
    assert "atom_site.fract_x" in result["text"]
    result = cryspy.call_assistant_tool(
        obj, "refine_parameters", {"pattern": "atom_site[O].fract"})
    # all coordinates of O are fixed by symmetry: nothing is done
    assert result["flag_error"]
    result = cryspy.call_assistant_tool(
        obj, "refine_parameters", {"pattern": "cell.length_a"})
    assert not result["flag_error"]


def test_fix_all_parameters():
    obj = take_object()
    call = cryspy.call_assistant_tool
    call(obj, "refine_parameters", {"pattern": "cell.length_a, b_iso"})
    assert len(obj.get_variable_names()) == 3
    result = call(obj, "fix_parameters", {"pattern": ""})
    assert "All 3 refined parameters" in result["text"]
    assert len(obj.get_variable_names()) == 0
    assert not call(obj, "fix_parameters", {})["flag_error"]
    assert call(obj, "fix_parameters", {"pattern": "polarization"})["flag_error"]


def test_structure_factor_of_one_reflection():
    obj = take_object()
    call = cryspy.call_assistant_tool
    # Pm-3m, Fe at 0 0 0 and O on 3d (1/2 0 0), (0 1/2 0), (0 0 1/2):
    # F(110) = b(Fe) T(Fe) + b(O) T(O) (-1 - 1 + 1)
    text_1 = call(obj, "calculate_structure_factor",
                  {"h": 1, "k": 1, "l": 0})["text"]
    text_2 = call(obj, "calculate_structure_factor", {"hkl": "(1 1 0)"})["text"]
    assert text_1 == text_2
    b_fe, b_o = 0.9450, 0.5803  # 10^-12 cm
    sthovl_sq = 2 / (2 * 5.0)**2
    f_110 = (b_fe * numpy.exp(-0.3 * sthovl_sq)
             - b_o * numpy.exp(-0.5 * sthovl_sq))
    f_text = float(re.search(r"F_N = (-?[\d.]+)", text_1).group(1))
    assert f_text == pytest.approx(f_110, abs=2e-3)
    assert call(obj, "calculate_structure_factor", {"hkl": "1 0"})["flag_error"]


# Arguments asked from the user ("missing", as MCP elicitation)

def tool_plan(cryspy_object, cell: list, geometry: str,
              wavelength: float = 0.) -> str:
    if wavelength == 0.:
        raise cryspy.ArgumentsNeeded(
            [{"name": "wavelength", "question": "the wavelength",
              "type": "number", "hint": "Å", "proposal": 1.26}],
            message="The wavelength is not given in the data.")
    return (f"Plan for {' '.join(f'{value:g}' for value in cell):} in "
            f"{geometry:} geometry at {wavelength:g} Å.")


@pytest.fixture
def plan_tool(monkeypatch):
    from cryspy.procedure_assistant import tools
    tool = tools.AssistantTool(
        "plan_measurement", "Plan of the measurement.", tool_plan,
        {"cell": {"type": "array", "items": {"type": "number"},
                  "minItems": 6, "maxItems": 6,
                  "question": "the unit cell parameters",
                  "hint": "a b c alpha beta gamma"},
         "geometry": {"type": "string", "enum": ["normal beam", "four circle"],
                      "question": "the geometry of measurements"},
         "wavelength": {"type": "number"}},
        required=("cell", "geometry"), flag_object_needed=False)
    monkeypatch.setattr(tools, "L_ASSISTANT_TOOLS",
                        tools.L_ASSISTANT_TOOLS + [tool])
    return tool


def test_missing_arguments_are_asked(plan_tool):
    result = cryspy.call_assistant_tool(None, "plan_measurement", {})
    assert not result["flag_error"] and not result["flag_object_changed"]
    assert result["text"] == ("plan_measurement needs to know: the unit cell "
                              "parameters; the geometry of measurements.")
    assert [d_missing["name"] for d_missing in result["missing"]] == [
        "cell", "geometry"]
    d_cell = result["missing"][0]
    assert d_cell["hint"] == "a b c alpha beta gamma" and d_cell["minItems"] == 6
    json.dumps(result["missing"])  # plain data for any program

    arguments = {"cell": cryspy.convert_assistant_answer(d_cell, "5 6 4 90 90 90"),
                 "geometry": cryspy.convert_assistant_answer(
                     result["missing"][1], "Normal Beam")}
    assert arguments == {"cell": [5., 6., 4., 90., 90., 90.],
                         "geometry": "normal beam"}
    # the tool asks for a value it finds missing itself
    result = cryspy.call_assistant_tool(None, "plan_measurement", arguments)
    assert result["text"] == ("The wavelength is not given in the data.\n"
                              "plan_measurement needs to know: the wavelength.")
    assert result["missing"] == [{"name": "wavelength", "question": "the wavelength",
                                  "type": "number", "hint": "Å", "proposal": 1.26}]
    arguments["wavelength"] = 1.26
    result = cryspy.call_assistant_tool(None, "plan_measurement", arguments)
    assert result["text"] == "Plan for 5 6 4 90 90 90 in normal beam geometry at 1.26 Å."
    assert result["missing"] == []


@pytest.mark.parametrize("text, message", [
    ("5 6 4", "6 values are expected, 3 are given."),
    ("5 6 4 90 90 x", "'x' is not a number."),
])
def test_wrong_answer_is_explained(plan_tool, text, message):
    d_cell = cryspy.call_assistant_tool(None, "plan_measurement", {})["missing"][0]
    with pytest.raises(ValueError, match=re.escape(message)):
        cryspy.convert_assistant_answer(d_cell, text)


def test_answer_types():
    convert = cryspy.convert_assistant_answer
    assert convert({"type": "integer"}, " 63 ") == 63
    assert convert({"type": "array", "items": {"type": "number"}}, "[5, 6; 4]") == [5., 6., 4.]
    with pytest.raises(ValueError, match="'Cmcm' is not an integer. Expected: number of the space group."):
        convert({"type": "integer", "hint": "number of the space group"}, "Cmcm")
    with pytest.raises(ValueError, match="is not one of: Biso, Uiso"):
        convert({"type": "string", "enum": ["Biso", "Uiso"]}, "Bani")


def test_required_argument_without_question_stays_error():
    result = cryspy.call_assistant_tool(take_object(), "set_parameter", {"value": 1})
    assert result["flag_error"] and result["missing"] == []
