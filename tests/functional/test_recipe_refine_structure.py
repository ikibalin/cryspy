"""The recipe of structure refinement on real single crystal data sets."""
import os

import cryspy
from cryspy.procedure_assistant.parameters import take_parameter_record

DIR = os.path.dirname(__file__)
F_MNF2 = os.path.join(DIR, "MEM_test", "mnf2_moments", "mnf2_unpol_SC.rcif")
F_MNBITE = os.path.join(DIR, "mnbi4te7_single_lambdaover2", "main.rcif")


def test_refine_tool_compares_chi_sq_per_point_and_r_f2():
    obj = cryspy.file_to_globaln(F_MNF2)
    cryspy.call_assistant_tool(obj, "fix_parameters")
    start = cryspy.call_assistant_tool(obj, "calculate")["metrics"]
    assert set(start["r_f2"]) == {"exp1"}
    cryspy.call_assistant_tool(obj, "refine_parameters", {"pattern": "scale"})
    result = cryspy.call_assistant_tool(obj, "refine")
    d_metrics = result["metrics"]
    assert d_metrics["chi_sq_per_n_previous"] == start["chi_sq_per_n"]
    assert d_metrics["chi_sq_per_n"] < start["chi_sq_per_n"]
    assert d_metrics["r_f2_previous"] == start["r_f2"]
    assert d_metrics["r_f2"]["exp1"] < start["r_f2"]["exp1"]
    assert d_metrics["n_refined"] == 1
    assert d_metrics["warnings"] == []
    assert "R_F2 of exp1:" in result["text"] and "better" in result["text"]


def test_run_recipe_undoes_negative_adp():
    """MnBi4Te7: free U_iso of Te become negative, the recipe takes the step
    back; single crystal: no background step.

    The scale of this data set is badly defined, so whether the positions
    step passes (correlations, a tiny change of chi_sq/n) depends on the
    platform; the parameters must follow its status either way."""
    obj = cryspy.file_to_globaln(F_MNBITE)
    result = cryspy.run_recipe(obj, "refine_structure")
    d_status = {name: status for name, status, _ in result["steps"]}
    d_text = {name: text for name, _, text in result["steps"]}
    assert d_status["background"] == "skipped"
    assert d_status["positions"] in ("done", "undone"), d_text["positions"]
    assert d_status["adp"] == "undone", d_text["adp"]
    record = take_parameter_record(obj, "crystal_mnbite.atom_site[Te1].u_iso_or_equiv")
    assert record.value >= 0.
    assert not record.flag_refinement
    flag_positions = take_parameter_record(obj, "crystal_mnbite.atom_site[Bi1].fract_z").flag_refinement
    assert flag_positions == (d_status["positions"] == "done"), d_text["positions"]
