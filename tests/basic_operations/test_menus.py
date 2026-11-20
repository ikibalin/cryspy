"""Menus of programs for the procedures and the methods: grouped by the
first word of the names, or by a line "Menu:" in the docstring; "Menu: -"
hides a function from the menus."""
import cryspy
from cryspy.procedure_assistant.menus import (
    MIN_METHODS_IN_SUBMENU, UNSORTED_PROCEDURES, group_methods_for_menu,
    group_procedures_for_menu, is_hidden_in_menu, take_first_word_menu, take_menu)
from cryspy.procedure_assistant.procedures import describe_procedure, parse_docstring


def make(name, doc="Do something."):
    def function(self):
        pass
    function.__name__ = name
    function.__doc__ = doc
    return function


def test_first_word():
    assert take_first_word_menu("calc_sthovl") == "Calculate"
    assert take_first_word_menu("report") == "Report"
    assert take_first_word_menu("plots_html") == "Plot"
    assert take_first_word_menu("rhochi_no_refinement") == "RhoChi"


def test_menu_line_of_the_docstring():
    doc = """Calculate the extinction.

    Menu: Extinction

    Parameters
    ----------
    x : float
        Value.
    """
    function = make("calc_extinction", doc)
    assert take_menu(function) == "Extinction"
    d_doc = parse_docstring(doc)
    assert d_doc["menu"] == "Extinction"
    assert "Menu" not in d_doc["summary"] + d_doc["description"]
    assert d_doc["parameters"]["x"]["description"] == "Value."
    assert take_menu(make("a", "Do.\n\nMenu: -")) == ""
    assert is_hidden_in_menu(make("a", "Do.\n\nMenu: -"))
    assert take_menu(make("a")) is None
    assert parse_docstring("Do.")["menu"] is None


def test_methods_by_first_word_with_threshold():
    d_function = {name: make(name) for name in
                  ["calc_a", "calc_b", "calc_c", "calc_d", "report_a", "report_b", "get_x"]}
    assert group_methods_for_menu(d_function, min_methods=4) == [
        ("Calculate", ["calc_a", "calc_b", "calc_c", "calc_d"]),
        (None, ["report_a", "report_b", "get_x"])]
    assert group_methods_for_menu(d_function, min_methods=5) == [(None, list(d_function))]


def test_methods_menu_line_and_hidden():
    d_function = {
        "report_a": make("report_a"),
        "estimate_b": make("estimate_b", "Estimate.\n\nMenu: Report"),
        "internal": make("internal", "Internal step.\n\nMenu: -"),
        "fit": make("fit", "Fit.\n\nMenu: Refinement"),
    }
    # an explicit menu is a submenu even for one method, and takes the
    # methods with its first word
    assert group_methods_for_menu(d_function) == [
        ("Refinement", ["fit"]), ("Report", ["report_a", "estimate_b"])]


def test_real_classes_lose_no_method():
    for cls in (cryspy.Crystal, cryspy.Cell, cryspy.Pd, cryspy.Diffrn, cryspy.AtomSiteL):
        groups = group_methods_for_menu(cls)
        l_name = [name for _, names in groups for name in names]
        l_visible = [name for name, function in vars(cls).items()
                     if callable(function) and not name.startswith("_")
                     and not is_hidden_in_menu(function)]
        assert sorted(l_name) == sorted(set(l_name)) and set(l_name) <= set(l_visible)
    d_crystal = dict(group_methods_for_menu(cryspy.Crystal))
    assert len(d_crystal["Calculate"]) >= MIN_METHODS_IN_SUBMENU


def test_procedures_menus():
    d_menu = group_procedures_for_menu()
    assert "RhoChi" in d_menu and "MEMPy" in d_menu
    l_all = [function for functions in d_menu.values() for function in functions]
    assert len(l_all) == len({function.__name__ for function in l_all})
    for menu, functions in d_menu.items():
        if menu != UNSORTED_PROCEDURES:
            assert len(functions) >= 2 or take_menu(functions[0]) == menu


def test_procedure_menu_line():
    def rhochi_extra(obj: cryspy.GlobalN):
        """Do the extra step.

        Menu: Extra
        """
    def mempy_hidden(obj: cryspy.GlobalN):
        """Hidden step.

        Menu: -
        """
    d_menu = group_procedures_for_menu([rhochi_extra, mempy_hidden])
    assert d_menu == {"Extra": [rhochi_extra]}
    d_procedure = describe_procedure(rhochi_extra)
    assert d_procedure["group"] == "Extra" and not d_procedure["hidden_in_menu"]
    assert d_procedure["summary"] == "Do the extra step."
    assert describe_procedure(mempy_hidden)["hidden_in_menu"]
    assert describe_procedure(mempy_hidden)["group"] == "MEMPy"


def test_list_methods_by_menus(monkeypatch):
    """A long list of methods for the assistant is shown by menus (no
    method is left out), and a menu name as pattern opens the menu."""
    import os
    from cryspy.procedure_assistant import tools
    f_rcif = os.path.join(os.path.dirname(__file__), "..", "..", "examples",
                          "rhochi_polarized_powder_1d_Ho2Ti2O7_5K_1T",
                          "rhochi_polarized_powder_1d_Ho2Ti2O7_5K_1T.rcif")
    obj = cryspy.load_file(f_rcif)
    text = cryspy.call_assistant_tool(obj, "list_methods", {"pattern": "Calculate"})["text"]
    assert text.startswith("Methods of the menu Calculate:")
    assert "calc_hkl: " in text and "report_ellipsoids" not in text
    monkeypatch.setattr(tools, "N_CHAR_METHOD_LINES", 500)
    text = cryspy.call_assistant_tool(obj, "list_methods", {})["text"]
    assert "menu Calculate: " in text and "calc_hkl" in text
    assert "estimate_background: " in text  # out of menus: with its summary
    assert "a menu name" in text
