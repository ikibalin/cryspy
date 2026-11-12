"""Selection of parameters one by one (tool refine_one_by_one)."""
import os

import pytest

import cryspy
from cryspy.procedure_assistant.parameters import take_parameter_record
from cryspy.procedure_assistant.refinement import refine_one_by_one

DIR = os.path.dirname(__file__)
F_MNF2 = os.path.join(DIR, "MEM_test", "mnf2_moments", "mnf2_unpol_SC.rcif")
PATTERN = "atom_site_moment[Mn]"


def take_mnf2():
    obj = cryspy.file_to_globaln(F_MNF2)
    cryspy.call_assistant_tool(obj, "fix_parameters")
    cryspy.call_assistant_tool(obj, "refine_parameters", {"pattern": "diffrn_exp1.phase.scale"})
    return obj


def test_keeps_only_the_moment_that_improves_the_fit():
    """MnF2 (unpolarized single crystal): the moments of Mn start at 0, where
    their gradient vanishes; with start_value only m_z decreases chi_sq/n."""
    obj = take_mnf2()
    cryspy.call_assistant_tool(obj, "calculate")
    l_info = []
    result = cryspy.call_assistant_tool(
        obj, "refine_one_by_one", {"pattern": PATTERN, "rejected_value": 0, "start_value": 1},
        progress=lambda info: l_info.append(info))
    assert not result["flag_error"], result["text"]
    d_decision = {step["path"]: step["decision"] for step in result["steps"]}
    path_z = "crystal_mnf2.atom_site_moment[Mn].crystalaxis_z"
    assert d_decision[path_z] == "kept", result["text"]
    assert list(d_decision.values()).count("kept") == 1, result["text"]
    for path, decision in d_decision.items():
        record = take_parameter_record(obj, path)
        assert record.flag_refinement == (decision == "kept")
        if decision == "fixed":
            assert record.value == 0.
    assert take_parameter_record(obj, "diffrn_exp1.phase.scale").flag_refinement
    assert abs(take_parameter_record(obj, path_z).value) > 3.
    assert result["metrics"]["chi_sq_per_n"] < result["metrics"]["chi_sq_per_n_previous"]
    assert {info["path"] for info in l_info} >= {path_z}
    assert all("step" in info and "n_step" in info for info in l_info)
    assert "kept" in result["text"] and "fixed" in result["text"]

    cryspy.call_assistant_tool(obj, "undo_refinement")
    record = take_parameter_record(obj, path_z)
    assert record.value == pytest.approx(5.1) and not record.flag_refinement


def test_stop_returns_the_values_of_the_step():
    obj = take_mnf2()
    path_z = "crystal_mnf2.atom_site_moment[Mn].crystalaxis_z"
    d_out = refine_one_by_one(obj, path_z, progress=lambda info: False if info["path"] else None)
    assert d_out["stopped"] and d_out["steps"] == []
    record = take_parameter_record(obj, path_z)
    assert record.value == pytest.approx(5.1) and not record.flag_refinement


def test_high_correlation_rejects_the_parameter():
    """With a correlation limit of 0 every new pair of refined parameters
    is 'high': the second candidate is fixed whatever chi_sq/n does."""
    obj = take_mnf2()
    d_out = refine_one_by_one(obj, "crystal_mnf2.atom_site_moment[Mn].crystalaxis_z, "
                              "crystal_mnf2.atom_site[Mn].b_iso_or_equiv", correlation_limit=0.)
    d_step = {step["path"]: step for step in d_out["steps"]}
    step_b = d_step["crystal_mnf2.atom_site[Mn].b_iso_or_equiv"]
    assert step_b["decision"] == "fixed" and "correlation" in step_b["reason"], d_out
    assert not take_parameter_record(obj, "crystal_mnf2.atom_site[Mn].b_iso_or_equiv").flag_refinement


def test_unknown_pattern_is_an_error():
    obj = take_mnf2()
    result = cryspy.call_assistant_tool(obj, "refine_one_by_one", {"pattern": "no_such_parameter"})
    assert result["flag_error"]
