"""Tables of the internal knowledge base of CrysPy (cryspy/knowledge/*.md).

The knowledge base is a folder of Markdown files written for people and
for AI assistants (searched as documents, see
procedure_assistant/knowledge.py). Some of them hold Markdown tables with
the data that CrysPy uses itself, so the values are written once, in a
file that people read and correct, and not in the code:

- the magnetic ions (magnetic_ions_*.md: "Ion", "Usual", "Lande g",
  "Moment (μB)", ...): the Lande factor of an ion, the free-ion moment, the
  usual ion of an element;
- the usual charges of the non-magnetic ions (ion_charges.md: "Element",
  "Charge") for the charge balance;
- the special points of the Brillouin zone (special_k_points.md:
  "Lattice", "Point", "k") to name a propagation vector.

neutron_absorbers.md is for reading only: its cross sections are those of
the CrysPy database (a test checks it).
"""
import os
import re
from fractions import Fraction
from functools import lru_cache
from typing import Dict, List, Optional

DIR_KNOWLEDGE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                             "knowledge")
MAGNETIC_ION_FILES = ("magnetic_ions_4f.md", "magnetic_ions_3d.md", "magnetic_ions_4d_5d_5f.md")
ION_CHARGE_FILE = "ion_charges.md"
RE_ION = re.compile(r"^([A-Z][a-z]?)(\d*)([+-])$")


def split_table_row(line: str) -> List[str]:
    return [cell.strip() for cell in line.strip().strip("|").split("|")]


def read_markdown_tables(path: str) -> List[List[Dict[str, str]]]:
    """Tables of a Markdown file: a list of tables, each a list of rows
    {column title: cell text}."""
    with open(path, "r", encoding="utf-8") as fid:
        l_line = fid.read().splitlines()
    l_table, ind = [], 0
    while ind < len(l_line) - 1:
        line, line_next = l_line[ind].strip(), l_line[ind + 1].strip()
        if line.startswith("|") and re.match(r"^\|[\s\-:|]+\|$", line_next):
            l_title = split_table_row(line)
            l_row, ind = [], ind + 2
            while ind < len(l_line) and l_line[ind].strip().startswith("|"):
                l_cell = split_table_row(l_line[ind])
                l_row.append(dict(zip(l_title, l_cell + [""] * (len(l_title) - len(l_cell)))))
                ind += 1
            l_table.append(l_row)
        else:
            ind += 1
    return l_table


def read_number(text: str) -> Optional[float]:
    """'6/7', '9.0', '2' -> float; '–', '' -> None."""
    text = text.strip()
    try:
        return float(Fraction(text))
    except (ValueError, ZeroDivisionError):
        return None


def split_ion(ion: str):
    """'Tb3+' -> ('Tb', 3); 'O2-' -> ('O', -2); 'Fe' -> ('Fe', None)."""
    match = RE_ION.match(ion.strip())
    if match is None:
        element = "".join(char for char in ion.strip() if char.isalpha())
        return element[:1].upper() + element[1:].lower(), None
    element, digits, sign = match.groups()
    charge = int(digits or 1)
    return element, charge if sign == "+" else -charge


@lru_cache(maxsize=None)
def take_magnetic_ions() -> Dict[str, dict]:
    """{ion: {"ion", "element", "charge", "usual" (bool), "configuration",
    "lande" (None when the ion has no moment), "moment" (free-ion moment in
    μB: g_J·J of a 4f ion, g·S of a 3d ion), "notes", "file"}} of the
    tables of magnetic_ions_*.md."""
    d_ion = {}
    for file_name in MAGNETIC_ION_FILES:
        for l_row in read_markdown_tables(os.path.join(DIR_KNOWLEDGE, file_name)):
            for row in l_row:
                if "Ion" not in row:
                    continue
                ion = row["Ion"]
                element, charge = split_ion(ion)
                d_ion[ion] = {
                    "ion": ion, "element": element, "charge": charge,
                    "usual": row.get("Usual", "").lower() == "yes",
                    "configuration": row.get("Configuration", ""),
                    "ground_term": row.get("Ground term", ""),
                    "lande": read_number(row.get("Lande g", "")),
                    "moment": read_number(row.get("Moment (μB)", "")) or 0.,
                    "notes": row.get("Notes", ""), "file": file_name}
    return d_ion


def take_ions_of_element(element: str) -> List[dict]:
    """Ions of the element in the tables, the usual one first."""
    l_ion = [d_ion for d_ion in take_magnetic_ions().values() if d_ion["element"] == element]
    return sorted(l_ion, key=lambda d_ion: not d_ion["usual"])


def take_magnetic_elements() -> List[str]:
    """Elements that have an ion with a moment in the tables."""
    l_element = []
    for d_ion in take_magnetic_ions().values():
        if d_ion["moment"] > 0. and d_ion["element"] not in l_element:
            l_element.append(d_ion["element"])
    return l_element


def take_usual_ion(element: str) -> Optional[dict]:
    """The usual magnetic ion of the element (Tb -> Tb3+), None if unknown."""
    for d_ion in take_ions_of_element(element):
        if d_ion["usual"]:
            return d_ion
    return None


def take_lande_factor(ion: str) -> float:
    """Lande factor of the ion from the tables (Tb3+: 1.5, Fe3+: 2); 2 (spin
    only) for an ion that is not there or has no moment."""
    d_ion = take_magnetic_ions().get(ion.strip())
    if d_ion is None or d_ion["lande"] is None:
        return 2.
    return d_ion["lande"]


@lru_cache(maxsize=None)
def take_ion_charges() -> Dict[str, int]:
    """{element: usual charge} of the non-magnetic ions (ion_charges.md)."""
    d_charge = {}
    for l_row in read_markdown_tables(os.path.join(DIR_KNOWLEDGE, ION_CHARGE_FILE)):
        for row in l_row:
            if "Element" in row and "Charge" in row:
                d_charge[row["Element"]] = int(row["Charge"].replace("+", ""))
    return d_charge


SPECIAL_K_POINT_FILE = "special_k_points.md"


@lru_cache(maxsize=None)
def take_special_k_points() -> Dict[str, List[tuple]]:
    """{lattice ('orthorhombic C'): [(label, k (3,) tuple of floats)]} of
    special_k_points.md (conventional reciprocal basis)."""
    d_point = {}
    for l_row in read_markdown_tables(os.path.join(DIR_KNOWLEDGE, SPECIAL_K_POINT_FILE)):
        for row in l_row:
            if not {"Lattice", "Point", "k"} <= set(row):
                continue
            k_vector = tuple(read_number(word) for word in row["k"].split())
            d_point.setdefault(row["Lattice"], []).append((row["Point"], k_vector))
    return d_point
