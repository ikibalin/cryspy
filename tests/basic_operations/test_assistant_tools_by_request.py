"""General (core) tools are always described to a model; the special
tools are found by the words of the request (find_matching_tools)."""
import pytest

from cryspy.procedure_assistant.tools import (
    CORE_TOOLS, L_ASSISTANT_TOOLS, find_matching_tools, form_matching_tools_text)


def names(request):
    return [tool.name for tool in find_matching_tools(request)]


def test_core_tools_are_tools():
    assert CORE_TOOLS <= {tool.name for tool in L_ASSISTANT_TOOLS}
    assert any(not tool.flag_core for tool in L_ASSISTANT_TOOLS)


@pytest.mark.parametrize("request_text, name", [
    ("Is an absorption correction needed?", "estimate_absorption"),
    ("Calculate the distances between atoms", "report_distances_and_angles"),
    ("Load the crystal from COD 1000041", "add_crystal_from_cod"),
    ("Convert the crystal to P1", "convert_to_p1"),
    ("Check the time reversal of the twins", "check_twin_time_reversal"),
    ("Add a constraint x+y=1", "add_constraint"),
    ("What do my documents say about extinction?", "search_documents"),
])
def test_special_tool_found(request_text, name):
    assert name in names(request_text)


@pytest.mark.parametrize("request_text", [
    "refine the cell", "What is the chi square?", "Show the structure", "refine",
])
def test_general_requests_get_no_special_tool(request_text):
    assert names(request_text) == []


def test_text_with_arguments_and_description():
    text = form_matching_tools_text("Load the crystal from COD 1000041")
    assert text.startswith("Actions that may do this:")
    assert "- add_crystal_from_cod(" in text and "Crystallography Open Database" in text
    assert form_matching_tools_text("refine the cell") == ""
