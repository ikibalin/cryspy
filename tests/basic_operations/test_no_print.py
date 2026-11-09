"""The library reports through logging, not print (see CLAUDE.md)."""
import ast
import os

import cryspy

DIR_CRYSPY = os.path.dirname(cryspy.__file__)

# print is the purpose of these: console entry point and helpers that show
# lists to the user
ALLOWED = {
    ("__main__.py", None),
    ("__init__.py", "functions"),
    ("__init__.py", "classes"),
    (os.path.join("A_functions_base", "function_2_space_group.py"), "print_long_list"),
}


def print_lines(node, rel_name):
    """Lines of print calls in node, outside the allowed functions."""
    if isinstance(node, ast.FunctionDef) and (rel_name, node.name) in ALLOWED:
        return []
    l_line = []
    if (isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
            and node.func.id == "print"):
        l_line.append(node.lineno)
    for child in ast.iter_child_nodes(node):
        l_line.extend(print_lines(child, rel_name))
    return l_line


def find_print_calls():
    l_found = []
    for dir_path, dir_names, file_names in os.walk(DIR_CRYSPY):
        for file_name in file_names:
            if not file_name.endswith(".py"):
                continue
            f_name = os.path.join(dir_path, file_name)
            rel_name = os.path.relpath(f_name, DIR_CRYSPY)
            if (rel_name, None) in ALLOWED:
                continue
            with open(f_name, encoding="utf-8") as fid:
                tree = ast.parse(fid.read())
            l_found.extend(
                f"{rel_name}:{lineno}" for lineno in print_lines(tree, rel_name))
    return l_found


def test_no_print_in_library():
    assert find_print_calls() == []


def test_cryspy_logger_has_null_handler():
    import logging
    handlers = logging.getLogger("cryspy").handlers
    assert any(isinstance(handler, logging.NullHandler) for handler in handlers)
