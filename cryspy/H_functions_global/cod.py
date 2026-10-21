"""Crystal structures from the Crystallography Open Database (COD).

    >>> import cryspy
    >>> cif_text = cryspy.load_cif_from_cod(formula="Al")
    >>> crystal = cryspy.cif_to_crystal(cif_text)
    >>> obj = cryspy.str_to_globaln("global_")
    >>> cryspy.add_crystal_from_cod(obj, formula="Al")

Only the standard library is used (urllib). The search goes through the
REST interface of COD, https://wiki.crystallography.net/RESTful_API/ ,
with the mirror qiserver.ugr.es as a fallback.
"""
import json
import logging
import re
import urllib.error
import urllib.parse
import urllib.request
from typing import List

from cryspy.B_parent_classes.cl_3_data import DataN
from cryspy.B_parent_classes.cl_4_global import GlobalN
from cryspy.C_item_loop_classes.cl_1_atom_site import AtomSiteL
from cryspy.C_item_loop_classes.cl_1_atom_site_aniso import AtomSiteAnisoL
from cryspy.C_item_loop_classes.cl_1_cell import Cell
from cryspy.C_item_loop_classes.cl_2_space_group import SpaceGroup
from cryspy.E_data_classes.cl_1_crystal import Crystal
from cryspy.H_functions_global.function_1_cryspy_objects import str_to_globaln

logger = logging.getLogger(__name__)

URL_COD = "https://www.crystallography.net/cod"
URL_COD_MIRROR = "https://qiserver.ugr.es/cod"
TIMEOUT = 30  # seconds


def form_hill_formula(formula: str) -> str:
    """Formula in the Hill notation with spaces, as COD searches it:
    'TiO2' -> 'O2 Ti', 'Al' -> 'Al', 'C2H6O' -> 'C2 H6 O'.

    Elements are C, then H, then the others alphabetically (all
    alphabetically without carbon); the count 1 is omitted.
    """
    l_element = re.findall(r"([A-Z][a-z]?)\s*(\d*\.?\d*)", formula)
    if len(l_element) == 0:
        raise ValueError(f"'{formula:}' is not a chemical formula.")
    d_count = {}
    for element, count in l_element:
        d_count[element] = d_count.get(element, 0.) + (
            float(count) if count != "" else 1.)
    l_order = sorted(d_count)
    if "C" in d_count:
        l_order = ["C"] + (["H"] if "H" in d_count else []) + [
            element for element in l_order if element not in ("C", "H")]
    ls_part = []
    for element in l_order:
        count = d_count[element]
        s_count = "" if count == 1. else f"{count:g}"
        ls_part.append(f"{element:}{s_count:}")
    return " ".join(ls_part)


def read_url(url: str, data: bytes = None) -> str:
    request = urllib.request.Request(
        url, data=data, headers={"User-Agent": "cryspy"})
    with urllib.request.urlopen(request, timeout=TIMEOUT) as response:
        return response.read().decode("utf-8", errors="replace")


def search_cod_ids(formula: str) -> List[str]:
    """COD numbers of the structures with this chemical formula.

    Raises ValueError if COD cannot be reached.
    """
    hill = form_hill_formula(formula)
    query = urllib.parse.urlencode({"formula": hill, "format": "json"})
    try:
        l_entry = json.loads(read_url(f"{URL_COD:}/result?{query:}"))
        return [str(entry["file"]) for entry in l_entry if "file" in entry]
    except (urllib.error.URLError, OSError, ValueError, KeyError,
            TypeError) as error:
        logger.info("COD search failed (%s), the mirror is tried.", error)
    data = urllib.parse.urlencode({"formula": hill}).encode("utf-8")
    try:
        text = read_url(f"{URL_COD_MIRROR:}/result.php", data=data)
    except (urllib.error.URLError, OSError) as error:
        raise ValueError(
            f"COD cannot be reached to search for '{hill:}': {error:}"
        ) from error
    return list(dict.fromkeys(re.findall(r"\"(\d+)\.cif\"", text)))


def load_cif_from_cod(formula: str = "", cod_id: str = "") -> str:
    """Text of a CIF file from COD: by its COD number, or the first
    structure found for the chemical formula ('Al', 'TiO2').

    Raises ValueError if nothing is found or COD cannot be reached.
    """
    cod_id = str(cod_id).strip()
    if cod_id == "":
        if formula.strip() == "":
            raise ValueError("Give a chemical formula or a COD number.")
        l_cod_id = search_cod_ids(formula)
        if len(l_cod_id) == 0:
            raise ValueError(
                f"COD has no structure with the formula "
                f"'{form_hill_formula(formula):}'.")
        cod_id = l_cod_id[0]
    if not cod_id.isdigit():
        raise ValueError(f"'{cod_id:}' is not a COD number.")
    try:
        return read_url(f"{URL_COD:}/{cod_id:}.cif")
    except (urllib.error.URLError, OSError) as error:
        raise ValueError(
            f"CIF file {cod_id:} cannot be loaded from COD: {error:}"
        ) from error


def take_type_symbol(label: str) -> str:
    """Element of an atom label: 'Al1' -> 'Al', 'O2a' -> 'O'."""
    match = re.match(r"([A-Z][a-z]?)", label.strip().capitalize())
    return match.group(1) if match is not None else label


# old tags of CIF 1.0 (still in many COD files) -> tags read by CrysPy
D_OLD_SPACE_GROUP_TAGS = {
    "_space_group_name_H-M_alt": ("_symmetry_space_group_name_H-M", ),
    "_space_group_IT_number": ("_symmetry_Int_Tables_number", ),
    "_space_group_name_Hall": ("_symmetry_space_group_name_Hall", ),
}


def take_modern_space_group_tags(cif_text: str) -> str:
    """CIF text where the space group is given by the tags of CrysPy.

    A file with only the old tags (_symmetry_space_group_name_H-M) gets the
    new ones; a file with no space group is taken as P 1 (logged).
    """
    l_line = []
    for tag, l_old_tag in D_OLD_SPACE_GROUP_TAGS.items():
        if re.search(rf"^\s*{re.escape(tag):}\s", cif_text,
                     flags=re.MULTILINE | re.IGNORECASE):
            continue
        for old_tag in l_old_tag:
            match = re.search(rf"^\s*{re.escape(old_tag):}\s+(.+?)\s*$",
                              cif_text, flags=re.MULTILINE | re.IGNORECASE)
            if match is not None:
                l_line.append(f"{tag:} {match.group(1):}")
                break
    if not re.search(r"^\s*_(space_group|symmetry)_", cif_text,
                     flags=re.MULTILINE | re.IGNORECASE):
        logger.info("No space group in the CIF text: P 1 is used.")
        l_line.append('_space_group_name_H-M_alt "P 1"')
    if len(l_line) == 0:
        return cif_text
    # the new tags go right after the first data_ line
    match = re.search(r"^\s*data_\S*\s*$", cif_text,
                      flags=re.MULTILINE | re.IGNORECASE)
    if match is None:
        return cif_text
    return (cif_text[:match.end():] + "\n" + "\n".join(l_line) + "\n"
            + cif_text[match.end():])


def cif_to_crystal(cif_text: str, data_name: str = "") -> Crystal:
    """Crystal of the first data block of a CIF file (from COD or other
    programs): cell, space group, atom_site and atom_site_aniso.

    Missing type symbols are taken from the labels; atoms without
    displacement parameters get adp_type 'Biso' with B = 0.
    """
    globaln = str_to_globaln(take_modern_space_group_tags(cif_text))
    l_data = [item for item in globaln.items if isinstance(item, DataN)]
    if len(l_data) == 0:
        raise ValueError("The CIF text has no data block.")
    data = l_data[0]
    d_item = {}
    for item in data.items:
        for item_class in (Cell, SpaceGroup, AtomSiteL, AtomSiteAnisoL):
            if isinstance(item, item_class) and item_class not in d_item:
                d_item[item_class] = item
    if Cell not in d_item or AtomSiteL not in d_item:
        raise ValueError("The CIF text has no cell or no atom_site loop.")
    for atom in d_item[AtomSiteL].items:
        if not atom.is_attribute("type_symbol"):
            atom.type_symbol = take_type_symbol(atom.label)
        if not atom.is_attribute("adp_type"):
            if atom.is_attribute("u_iso_or_equiv"):
                atom.adp_type = "Uiso"
            else:
                atom.adp_type = "Biso"
                if not atom.is_attribute("b_iso_or_equiv"):
                    atom.b_iso_or_equiv = 0.
        if not atom.is_attribute("occupancy"):
            atom.occupancy = 1.
    if SpaceGroup not in d_item:
        raise ValueError("The CIF text has no space group.")
    if data_name.strip() == "":
        data_name = re.sub(r"\W", "_", str(data.data_name)) or "cod"
    crystal = Crystal(data_name=data_name, items=list(d_item.values()))
    crystal.apply_constraints()
    return crystal


def add_crystal_from_cod(globaln: GlobalN, formula: str = "",
                         cod_id: str = "", data_name: str = "") -> Crystal:
    """Load a structure from COD and add it to the object as a crystal.

    The name is the given one, otherwise the formula or 'cod_<number>';
    a name already used in the object gets a number.
    """
    cif_text = load_cif_from_cod(formula=formula, cod_id=cod_id)
    if data_name.strip() == "":
        data_name = formula.strip() if cod_id.strip() == "" else f"cod_{cod_id:}"
    data_name = re.sub(r"\W", "_", data_name.strip()).lower() or "cod"
    l_used = [str(getattr(item, "data_name", "")).lower()
              for item in globaln.items]
    name, index = data_name, 1
    while name in l_used:
        index += 1
        name = f"{data_name:}_{index:}"
    crystal = cif_to_crystal(cif_text, data_name=name)
    globaln.add_items([crystal])
    return crystal
