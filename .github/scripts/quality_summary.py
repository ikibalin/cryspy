"""Markdown of the code duplication (jscpd report) and of the docstring coverage of a package.

Usage: python quality_summary.py <jscpd-report.json> <package directory>

Only a report for the pull request comment: nothing here fails the CI.
"""
import ast
import collections
import json
import os
import sys

N_CLONES = 5


def take_package(path: str, root: str) -> str:
    """First sub-package of a file inside root ("(top level)" for the files of root)."""
    parts = os.path.relpath(path, root).split(os.sep)
    return parts[0] if len(parts) > 1 else "(top level)"


def form_duplication_lines(file_report: str, root: str) -> list:
    """Table of the duplicated lines per package and the largest clones."""
    if not os.path.isfile(file_report):
        return ["### Code duplication", "", "No jscpd report.", ""]
    with open(file_report, encoding="utf-8") as fid:
        report = json.load(fid)
    total = report["statistics"]["total"]
    by_package = collections.defaultdict(lambda: [0, 0])
    for name, stat in report["statistics"].get("formats", {}).get("python", {}).get(
            "sources", {}).items():
        package = by_package[take_package(name, root)]
        package[0] += stat["lines"]
        package[1] += stat["duplicatedLines"]
    ls_out = ["### Code duplication: {:.1f} %".format(total["percentage"]), "",
              "Lines of code in clones of 6 lines or more (jscpd); {:} clones.".format(
                  total["clones"]), "",
              "| Package | Lines | Duplicated |", "|---|---:|---:|"]
    for name in sorted(by_package):
        lines, duplicated = by_package[name]
        ls_out.append("| `{:}` | {:} | {:.0f} % |".format(name, lines, 100. * duplicated / max(lines, 1)))
    clones = sorted(report["duplicates"], key=lambda item: -item["lines"])[:N_CLONES]
    if clones:
        ls_out += ["", "Largest clones:", ""]
        for clone in clones:
            first, second = clone["firstFile"], clone["secondFile"]
            ls_out.append("- {:} lines: `{:}:{:}` and `{:}:{:}`".format(
                clone["lines"], os.path.relpath(first["name"], os.path.dirname(root)),
                first["start"], os.path.relpath(second["name"], os.path.dirname(root)),
                second["start"]))
    return ls_out + [""]


def count_docstrings(file_name: str) -> tuple:
    """Public functions, classes and methods of a file and how many have a docstring."""
    with open(file_name, encoding="utf-8") as fid:
        tree = ast.parse(fid.read())
    n_all, n_doc = 0, 0
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)) and \
                not node.name.startswith("_"):
            n_all += 1
            n_doc += ast.get_docstring(node) is not None
    return n_all, n_doc


def form_docstring_lines(root: str) -> list:
    """Table of the public objects with a docstring per package."""
    by_package = collections.defaultdict(lambda: [0, 0])
    for directory, _, files in os.walk(root):
        for name in files:
            if name.endswith(".py"):
                n_all, n_doc = count_docstrings(os.path.join(directory, name))
                package = by_package[take_package(os.path.join(directory, name), root)]
                package[0] += n_all
                package[1] += n_doc
    n_all = sum(value[0] for value in by_package.values())
    n_doc = sum(value[1] for value in by_package.values())
    ls_out = ["### Docstrings: {:.1f} %".format(100. * n_doc / max(n_all, 1)), "",
              "Public functions, classes and methods with a docstring (the menus and the "
              "assistant read them).", "",
              "| Package | Public | With docstring |", "|---|---:|---:|"]
    for name in sorted(by_package):
        total, documented = by_package[name]
        ls_out.append("| `{:}` | {:} | {:.0f} % |".format(name, total, 100. * documented / max(total, 1)))
    return ls_out + [""]


def main(file_report: str, root: str) -> str:
    root = os.path.abspath(root)
    return "\n".join(form_duplication_lines(file_report, root) + form_docstring_lines(root))


if __name__ == "__main__":
    print(main(sys.argv[1], sys.argv[2]))
