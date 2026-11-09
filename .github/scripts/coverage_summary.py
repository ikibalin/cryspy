"""Markdown table of test coverage per subpackage, from `coverage json`.

Usage: python coverage_summary.py coverage.json > coverage.md
"""
import json
import sys
from collections import defaultdict

with open(sys.argv[1], encoding="utf-8") as fid:
    report = json.load(fid)

statements = defaultdict(int)
covered = defaultdict(int)
for file_name, data in report["files"].items():
    # file names may be relative or absolute: group by the folder right
    # below the last "cryspy" folder
    parts = file_name.replace("\\", "/").split("/")
    i_cryspy = max(i for i, part in enumerate(parts[:-1]) if part == "cryspy")
    sub_parts = parts[i_cryspy + 1:]
    package = sub_parts[0] if len(sub_parts) > 1 else "(top level)"
    statements[package] += data["summary"]["num_statements"]
    covered[package] += data["summary"]["covered_lines"]

total = report["totals"]
lines = [
    "<!-- cryspy-coverage -->",
    f"### Test coverage: {total['percent_covered']:.1f} %",
    "",
    "Unit, functional and examples tests together. "
    "The HTML report is in the `coverage-html` artifact of the run.",
    "",
    "| Package | Statements | Covered |",
    "|---|---:|---:|",
]
for package in sorted(statements):
    n_statement = statements[package]
    percent = 100. * covered[package] / n_statement if n_statement else 100.
    lines.append(f"| `{package}` | {n_statement} | {percent:.0f} % |")
lines.append(
    f"| **total** | **{total['num_statements']}** | "
    f"**{total['percent_covered']:.1f} %** |")
print("\n".join(lines))
