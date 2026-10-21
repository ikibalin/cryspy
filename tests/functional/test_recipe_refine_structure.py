"""The recipe of structure refinement on a real single crystal data set."""
import os

import numpy

import cryspy
from cryspy.procedure_assistant.parameters import take_parameter_record

F_NAME = os.path.join(os.path.dirname(__file__), "twin_tbte3", "main.rcif")


def test_refine_tool_compares_with_previous_calculation():
    obj = cryspy.file_to_globaln(F_NAME)
    cryspy.call_assistant_tool(obj, "fix_parameters")
    start = cryspy.call_assistant_tool(obj, "calculate")
    cryspy.call_assistant_tool(obj, "refine_parameters", {"pattern": "scale"})
    result = cryspy.call_assistant_tool(obj, "refine")
    d_metrics = result["metrics"]
    assert d_metrics["chi_sq_previous"] == start["metrics"]["chi_sq"]
    assert d_metrics["chi_sq"] < d_metrics["chi_sq_previous"]
    assert d_metrics["n_point"] == start["metrics"]["n_point"]
    assert d_metrics["n_refined"] == 1
    assert "better" in result["text"]


def test_run_recipe_undoes_negative_adp():
    """TbTe3: free U_iso become negative, the recipe takes the step back."""
    obj = cryspy.file_to_globaln(F_NAME)
    result = cryspy.run_recipe(obj, "refine_structure")
    d_status = {name: status for name, status, _ in result["steps"]}
    assert d_status == {"start": "done", "scale": "done", "background": "skipped",
                        "positions": "done", "adp": "undone", "extinction": "done"}
    record = take_parameter_record(obj, "crystal_tbte3.atom_site[Tb].u_iso_or_equiv")
    assert numpy.isclose(record.value, 0.)
    assert not record.flag_refinement
    assert take_parameter_record(obj, "crystal_tbte3.atom_site[Tb].fract_y").flag_refinement
