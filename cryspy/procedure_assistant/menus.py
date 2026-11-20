"""Menus of programs (CrysPy Editor) for the procedures and the methods.

The procedures (procedures.take_procedures) and the methods of the data
blocks, loops and items are grouped into menus by the first word of their
names, so that a new function finds its place without a list of menus
anywhere. A line "Menu:" in the docstring sets the menu by hand:

    def calc_extinction(...):
        '''Calculate the extinction of the reflections.

        Menu: Extinction
        '''

"Menu: -" hides the function from the menus (it is still a procedure or a
method for AI assistants and for scripts). The line is not a part of the
summary and description of the docstring (procedures.parse_docstring).

    >>> import cryspy
    >>> cryspy.group_methods_for_menu(cryspy.Crystal)
    [('Add', [...]), ('Calculate', [...]), ..., (None, [...])]
    >>> cryspy.group_procedures_for_menu()
    {'MEMPy': [...], 'RhoChi': [...], ..., 'Unsorted procedures': [...]}
"""
import inspect
import re
from types import FunctionType
from typing import Callable, Dict, List, Optional, Tuple, Union

# a submenu of the methods appears when at least so many methods share the
# first word (a menu set by "Menu:" in the docstring is always a submenu)
MIN_METHODS_IN_SUBMENU = 4
# a menu of the procedures appears when at least so many procedures share
# the first word; the others go to UNSORTED_PROCEDURES
MIN_PROCEDURES_IN_MENU = 2
UNSORTED_PROCEDURES = "Unsorted procedures"
# "Menu: -" (or "Menu: none") hides the function from the menus
HIDDEN = ("-", "none", "hidden")

# titles of the menus for the first words of the names; other words are
# capitalized ("refine" -> "Refine")
D_MENU_WORD = {
    "calc": "Calculate", "mempy": "MEMPy", "rhochi": "RhoChi",
    "plot": "Plot", "plots": "Plot", "del": "Delete", "rm": "Remove",
    "is": "Check",
}

MENU_LINE = re.compile(r"^\s*Menu\s*:\s*(.*?)\s*$", re.IGNORECASE)


def split_menu_line(doc: Optional[str]) -> Tuple[str, Optional[str]]:
    """(docstring without the line "Menu: ...", the text of that line or
    None when there is none)."""
    if not doc:
        return doc or "", None
    l_line = inspect.cleandoc(doc).splitlines()
    menu = None
    l_keep = []
    for line in l_line:
        match = MENU_LINE.match(line)
        if match is not None and menu is None:
            menu = match.group(1)
            continue
        l_keep.append(line)
    if menu is None:
        return doc, None
    return "\n".join(l_keep), menu


def take_menu(function: Callable) -> Optional[str]:
    """The menu set in the docstring ("Menu: Calculate" -> "Calculate"),
    "" when the function is hidden from the menus ("Menu: -"), None when
    the docstring sets no menu."""
    menu = split_menu_line(getattr(function, "__doc__", None))[1]
    if menu is None:
        return None
    if menu.lower() in HIDDEN or menu == "":
        return ""
    return menu


def is_hidden_in_menu(function: Callable) -> bool:
    """True for "Menu: -" in the docstring."""
    return take_menu(function) == ""


def take_first_word_menu(name: str) -> str:
    """Menu of a function by the first word of its name.

    >>> take_first_word_menu("calc_sthovl")
    'Calculate'
    >>> take_first_word_menu("report")
    'Report'
    >>> take_first_word_menu("rhochi_no_refinement")
    'RhoChi'
    """
    word = name.split("_")[0].lower()
    return D_MENU_WORD.get(word, word.capitalize())


def group_functions(
        d_function: Dict[str, Callable], min_functions: int,
) -> Tuple[Dict[str, List[str]], List[str]]:
    """({menu: names}, names out of the menus) in the order of d_function,
    without the hidden functions. A menu set in the docstring is always
    there; a menu by the first word of the names needs min_functions."""
    d_menu, d_word = {}, {}
    for name, function in d_function.items():
        menu = take_menu(function)
        if menu == "":
            continue
        if menu is not None:
            d_menu.setdefault(menu, []).append(name)
        else:
            d_word.setdefault(take_first_word_menu(name), []).append(name)
    l_alone = []
    for menu, names in d_word.items():
        if menu in d_menu or len(names) >= min_functions:
            d_menu.setdefault(menu, []).extend(names)
        else:
            l_alone.extend(names)
    l_order = [name for name in d_function]
    d_menu = {menu: sorted(names, key=l_order.index) for menu, names in d_menu.items()}
    l_alone = sorted(l_alone, key=l_order.index)
    return d_menu, l_alone


def take_own_methods(cls: type) -> Dict[str, Callable]:
    """{name: function} of the public methods defined in the class itself
    (the methods of the context menu of an object in CrysPy Editor)."""
    return {name: function for name, function in vars(cls).items()
            if isinstance(function, FunctionType) and not name.startswith("_")}


def group_methods_for_menu(
        methods: Union[type, object, Dict[str, Callable]],
        min_methods: int = MIN_METHODS_IN_SUBMENU,
) -> List[Tuple[Optional[str], List[str]]]:
    """Methods of a class (or an object, or {name: function}) for a menu.

    Returns [(submenu, names), ...] in alphabetical order of the submenus,
    then (None, names) of the methods directly in the menu. The methods
    with "Menu: -" in the docstring are left out.
    """
    if isinstance(methods, dict):
        d_function = methods
    else:
        d_function = take_own_methods(methods if isinstance(methods, type) else type(methods))
    d_menu, l_alone = group_functions(d_function, min_methods)
    l_out = [(menu, d_menu[menu]) for menu in sorted(d_menu)]
    if l_alone:
        l_out.append((None, l_alone))
    return l_out


def group_procedures_for_menu(
        functions: Optional[List[Callable]] = None,
        min_procedures: int = MIN_PROCEDURES_IN_MENU,
) -> Dict[str, List[Callable]]:
    """{menu: procedures} of the menus of CrysPy Editor: the procedures
    (procedures.take_procedures, or the given functions, e.g. with the
    functions of the user) by the first word of their names or by
    "Menu:" in the docstring; single ones go to "Unsorted procedures";
    "Menu: -" hides a procedure."""
    if functions is None:
        from cryspy.procedure_assistant.procedures import take_procedures
        functions = list(take_procedures().values())
    d_function = {}
    for function in functions:
        d_function.setdefault(function.__name__, function)
    d_menu, l_alone = group_functions(d_function, min_procedures)
    d_out = {menu: [d_function[name] for name in d_menu[menu]] for menu in sorted(d_menu)}
    if l_alone:
        d_out.setdefault(UNSORTED_PROCEDURES, []).extend(d_function[name] for name in l_alone)
    return d_out
