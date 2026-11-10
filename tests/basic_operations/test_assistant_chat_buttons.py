"""Data for the buttons of the chat with the assistant: options of
questions, requests proposed from the object, progress of long tools."""
import pytest

import cryspy

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


@pytest.fixture
def single_crystal():
    """The crystal with simulated single crystal intensities, its scale
    changed so that a refinement has something to do."""
    cryspy_object = cryspy.str_to_globaln(RCIF)
    cryspy.simulation_single_crystal(cryspy_object)
    cryspy_object.diffrn_sc.phase.scale = 0.7
    return cryspy_object


def test_requests_for_empty_object_and_crystal():
    l_request = cryspy.suggest_assistant_requests(None)
    assert {"label": "What can you do?", "request": "What can you do?"} in l_request
    assert any(item["request"] == "@add_crystal_from_cod" for item in l_request)
    l_request = cryspy.suggest_assistant_requests(cryspy.str_to_globaln(RCIF))
    l_text = [item["request"] for item in l_request]
    assert "@run_procedure simulation_unpolarized_neutron_powder_diffraction" in l_text
    assert "@summarize_object" in l_text
    assert len(l_request) <= 5


def test_requests_for_experiments_are_recipes(single_crystal):
    l_request = cryspy.suggest_assistant_requests(single_crystal)
    assert l_request[0]["request"].startswith("@refine")
    assert "sc" in l_request[0]["label"]
    # no simulation is proposed when there are experiments
    assert not any("simulation" in item["request"] for item in l_request)


def test_refine_without_marked_parameters_asks_with_options(single_crystal):
    result = cryspy.call_assistant_tool(single_crystal, "refine")
    assert not result["flag_error"]
    [d_missing] = result["missing"]
    assert d_missing["name"] == "parameters" and "scale" in d_missing["options"]
    result = cryspy.call_assistant_tool(single_crystal, "refine", {"parameters": "scale"})
    assert not result["flag_error"], result["text"]
    assert "Marked for refinement" in result["text"]
    assert single_crystal.diffrn_sc.phase.scale == pytest.approx(1., rel=1e-3)
    assert "chi_sq_per_n" in result["metrics"]


def test_refine_reports_progress_and_stops(single_crystal):
    single_crystal.diffrn_sc.phase.scale_refinement = True
    l_info = []
    result = cryspy.call_assistant_tool(single_crystal, "refine", progress=l_info.append)
    assert not result["flag_error"] and len(l_info) > 0 and "chi_sq_per_n" in l_info[0]
    single_crystal.diffrn_sc.phase.scale = 0.5
    result = cryspy.call_assistant_tool(single_crystal, "refine", progress=lambda info: False)
    assert result["text"].startswith("Refinement stopped by the user.")


def test_other_keys_of_tool_results_pass(monkeypatch):
    from cryspy.procedure_assistant import tools
    tool = tools.take_assistant_tool("summarize_object")
    monkeypatch.setattr(tool, "function", lambda obj: {"text": "done", "extra": [1, 2]})
    result = cryspy.call_assistant_tool(cryspy.str_to_globaln(RCIF), "summarize_object")
    assert result["text"] == "done" and result["extra"] == [1, 2] and result["metrics"] == {}
