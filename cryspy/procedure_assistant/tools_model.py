"""Assistant tools that read and change the model: attributes, items, atoms and crystals.

They are registered in tools.py (L_ASSISTANT_TOOLS).
"""
import re
from typing import Optional

import numpy

from cryspy.B_parent_classes.attributes import take_allowed_values
from cryspy.B_parent_classes.cl_3_data import DataN
from cryspy.C_item_loop_classes.cl_1_atom_site import AtomSite
from cryspy.E_data_classes.cl_1_crystal import Crystal

from cryspy.procedure_assistant.attributes import (
    describe_item_class, explain_allowed_value, find_attribute_records, find_item_classes,
    take_attribute_record, take_block, take_object_item_classes)
from cryspy.procedure_assistant import items
from cryspy.procedure_assistant.procedures import (
    D_BASIC_TYPE)

from cryspy.procedure_assistant.tool_base import (
    MAX_PARAMETERS_IN_LIST, ArgumentsNeeded, format_value)


def take_data_block(cryspy_object, name: str):
    name_sh = name.strip().lower()
    for item in cryspy_object.items:
        if isinstance(item, DataN):
            if name_sh in (str(item.data_name).lower(), item.get_name()):
                return item
    return None


def tool_show_model(cryspy_object, block: str = "",
                    max_length: int = 8000) -> str:
    if block.strip() == "":
        text = cryspy_object.to_cif()
    else:
        data = take_data_block(cryspy_object, block)
        if data is None:
            return f"Error: data block '{block:}' is not found."
        text = data.to_cif()
    if len(text) > max_length:
        text = text[:max_length] + "\n[... the rest is cut ...]"
    return text


def tool_list_attributes(cryspy_object, pattern: str = "") -> str:
    l_record = find_attribute_records(cryspy_object, pattern)
    if len(l_record) == 0:
        s_explanation = explain_allowed_value(pattern)
        if s_explanation:
            return f"Error: no attributes match '{pattern:}'. {s_explanation:}"
        return (f"Error: no attributes match '{pattern:}'. Use summarize_object to see "
                "the data blocks and describe_category for the attributes of a category.")
    ls_out = ["path | value | allowed values"]
    for record in l_record[:MAX_PARAMETERS_IN_LIST]:
        allowed = record.allowed_values
        s_allowed = ""
        if allowed is not None and len(allowed) <= 30:
            s_allowed = ", ".join(str(_) for _ in allowed)
        elif allowed is not None:
            s_allowed = f"{len(allowed):} values"
        ls_out.append(f"{record.path:} | {format_value(record.value):} | {s_allowed:}")
    if len(l_record) > MAX_PARAMETERS_IN_LIST:
        ls_out.append(f"... {len(l_record) - MAX_PARAMETERS_IN_LIST:} more, "
                      "give a narrower pattern.")
    return "\n".join(ls_out)


def tool_set_attribute(cryspy_object, path: str, value: str) -> str:
    record = take_attribute_record(cryspy_object, path)
    if record is None:
        l_record = find_attribute_records(cryspy_object, path, flag_defined_only=False)
        if len(l_record) == 1:
            record = l_record[0]
        elif len(l_record) > 1:
            return ("Error: the path is ambiguous, give one of:\n"
                    + "\n".join(record.path for record in l_record[:20]))
        else:
            return (f"Error: no attribute '{path:}'. Use list_attributes to find the path.")
    if record.flag_constraint:
        return f"Error: {record.path:} is constrained by symmetry and cannot be set."
    old_value = record.value
    try:
        record.set_value(value)
    except ValueError as error:
        return f"Error: {error:}"
    return f"{record.path:}: {format_value(old_value):} -> {format_value(record.value):}"


def tool_show_item(cryspy_object, path: str, max_length: int = 8000) -> str:
    block = take_block(cryspy_object, path)
    s_value = ""
    if block is None:  # the path of one attribute: its value and its item
        record = take_attribute_record(cryspy_object, path)
        if record is None:
            l_record = find_attribute_records(cryspy_object, path)
            record = l_record[0] if len(l_record) == 1 else None
        if record is None:
            l_name = [str(data.get_name()) for data in getattr(cryspy_object, "items", [])]
            return (f"Error: '{path:}' is not found. Give a data block, a category, an item of a "
                    "loop or an attribute: 'crystal_fe', 'crystal_fe.cell', "
                    "'crystal_fe.atom_site[Fe]', 'diffrn_exp1.setup.field'. The object has: "
                    f"{', '.join(l_name):}. list_attributes finds the paths of attributes.")
        block = record.item
        s_value = f"{record.path:} = {format_value(record.value):}\n\n"
    text = s_value + block.to_cif()
    if len(text) > max_length:
        text = text[:max_length] + "\n[... the rest is cut ...]"
    return text


def tool_describe_category(cryspy_object, name: str = "") -> str:
    if name.strip() == "":
        l_class = take_object_item_classes(cryspy_object) if cryspy_object is not None else []
        if len(l_class) == 0:
            return "Error: give the name of a category, e.g. 'cell' or 'atom_site'."
    else:
        l_class = find_item_classes(name)
        if len(l_class) == 0:
            return f"Error: there is no category '{name:}'."
        if len(l_class) > 10:
            return ("Several categories match, give one of: "
                    + ", ".join(item_class.PREFIX for item_class in l_class))
    return "\n\n".join(describe_item_class(item_class) for item_class in l_class)


def tool_add_atom(cryspy_object, crystal: str, label: str, type_symbol: str,
                  fract_x: float, fract_y: float, fract_z: float,
                  occupancy: float = 1., adp_type: str = "",
                  b_iso_or_equiv: float = 0.) -> str:
    data = take_data_block(cryspy_object, crystal)
    if not isinstance(data, Crystal):
        return f"Error: crystal '{crystal:}' is not found."
    atom_site = data.atom_site
    if label in [str(item.label) for item in atom_site.items]:
        return f"Error: atom '{label:}' already exists in '{crystal:}'."
    if adp_type == "":
        # the same kind of isotropic parameter as the other atoms
        l_adp_type = [item.adp_type for item in atom_site.items
                      if item.is_attribute("adp_type")]
        adp_type = "Uiso" if "Uiso" in l_adp_type else "Biso"
    if adp_type not in ("Biso", "Uiso"):
        return ("Error: adp_type should be 'Biso' or 'Uiso'; "
                "anisotropic parameters are set in atom_site_aniso.")
    atom = AtomSite(label=label, type_symbol=type_symbol, fract_x=fract_x,
                    fract_y=fract_y, fract_z=fract_z, occupancy=occupancy,
                    adp_type=adp_type)
    if adp_type == "Biso":
        atom.b_iso_or_equiv = b_iso_or_equiv
    else:
        atom.u_iso_or_equiv = b_iso_or_equiv
    # keep the same columns as the other atoms, so the loop stays printable
    for name in ("b_iso_or_equiv", "u_iso_or_equiv"):
        if not atom.is_attribute(name) and all(
                item.is_attribute(name) for item in atom_site.items):
            setattr(atom, name, 0.)
    atom_site.items.append(atom)
    data.apply_constraints()
    return (f"Atom {label:} ({type_symbol:}) is added to {crystal:} at "
            f"({fract_x:}, {fract_y:}, {fract_z:}).")


def tool_add_item(cryspy_object, category: str, block: str = "",
                  values: Optional[dict] = None, **extra_values) -> str:
    # the values may come as a dict or as separate arguments
    d_value = dict(values or {}, **extra_values)
    try:
        path, container = items.add_item(cryspy_object, category, d_value, block)
    except items.BlockNeeded as error:
        raise ArgumentsNeeded([{
            "name": "block", "type": "string", "enum": list(error.args[0]),
            "question": f"the data block where {category.strip():} is added"}]) from error
    except items.ValuesNeeded as error:
        item_class = items.take_category_class(category)
        l_missing = []
        for name in error.args:
            d_missing = {"name": name, "question": f"the {name:} of {item_class.PREFIX:}",
                         "type": D_BASIC_TYPE.get(
                             item_class.ATTR_TYPES[item_class.ATTR_NAMES.index(name)], "string")}
            allowed = take_allowed_values(item_class, name)
            if allowed is not None and len(allowed) <= 30:
                d_missing["enum"] = [str(value) for value in allowed]
            l_missing.append(d_missing)
        raise ArgumentsNeeded(l_missing, f"{item_class.PREFIX:} needs these values.") from error
    except ValueError as error:
        return f"Error: {error:}"
    text = container.to_cif()
    if len(text) > 3000:
        text = text[:3000] + "\n[... the rest is cut ...]"
    return f"{path:} is added:\n{text:}"


def tool_remove_item(cryspy_object, path: str) -> str:
    try:
        removed = items.remove_item(cryspy_object, path)
    except ValueError as error:
        return f"Error: {error:}"
    return f"{removed:} is removed."


def take_hkl(h, k, l, hkl: str):
    """numpy array (3, 1) of the indices given as h, k, l or as '5 0 3'."""
    if hkl.strip() != "":
        l_index = re.findall(r"-?\d+(?:\.\d*)?", hkl)
        if len(l_index) != 3:
            raise ValueError(f"hkl '{hkl:}' should have three indices.")
        h, k, l = (float(index) for index in l_index)
    if h is None or k is None or l is None:
        raise ValueError("give the indices h, k, l (or hkl as '5 0 3').")
    return numpy.array([[h], [k], [l]], dtype=float)


def tool_calculate_structure_factor(
        cryspy_object, h: float = None, k: float = None, l: float = None,
        hkl: str = "", crystal: str = "") -> str:
    try:
        index_hkl = take_hkl(h, k, l, hkl)
    except ValueError as error:
        return f"Error: {error:}"
    if crystal.strip() == "":
        l_crystal = [item for item in cryspy_object.items
                     if isinstance(item, Crystal)]
    else:
        l_crystal = [take_data_block(cryspy_object, crystal)]
        if not isinstance(l_crystal[0], Crystal):
            return f"Error: crystal '{crystal:}' is not found."
    if len(l_crystal) == 0:
        return "Error: the object has no crystal."
    s_hkl = " ".join(f"{index:g}" for index in index_hkl[:, 0])
    ls_out = [f"Nuclear structure factor F_N of the reflection ({s_hkl:}), "
              "in 10^-12 cm (per unit cell, with the Debye-Waller factor):"]
    for item in l_crystal:
        f_nucl = complex(numpy.asarray(item.calc_f_nucl(index_hkl)).ravel()[0])
        sthovl = float(numpy.asarray(item.cell.calc_sthovl(
            *index_hkl[:, 0])).ravel()[0])
        ls_out.append(
            f"{item.data_name:}: F_N = {f_nucl.real:.5f} {f_nucl.imag:+.5f}i, "
            f"|F_N| = {abs(f_nucl):.5f}, |F_N|^2 = {abs(f_nucl)**2:.5f}, "
            f"sin(theta)/lambda = {sthovl:.5f} 1/Å")
    return "\n".join(ls_out)


def take_named_crystals(cryspy_object, crystal: str = "") -> list:
    """Crystals of the object, or the one named; ValueError if not found."""
    if crystal.strip() == "":
        l_crystal = [item for item in cryspy_object.items if isinstance(item, Crystal)]
        if len(l_crystal) == 0:
            raise ValueError("the object has no crystal.")
        return l_crystal
    item = take_data_block(cryspy_object, crystal)
    if not isinstance(item, Crystal):
        raise ValueError(f"crystal '{crystal:}' is not found.")
    return [item]


def tool_list_atoms_in_unit_cell(cryspy_object, crystal: str = "",
                                 max_atoms: int = 500) -> str:
    try:
        l_crystal = take_named_crystals(cryspy_object, crystal)
    except ValueError as error:
        return f"Error: {error:}"
    ls_out = []
    for item in l_crystal:
        atoms = item.calc_atoms_in_unit_cell()
        moment = atoms["moment_crystalaxis"]
        n_atom = len(atoms["label"])
        ls_out.append(f"Crystal {item.data_name:}: {n_atom:} atoms in the unit cell "
                      "(generated by the symmetry elements from atom_site).")
        s_head = "| atom | from | x | y | z |"
        s_line = "|---|---|---|---|---|"
        if moment is not None:
            s_head += " M_a, M_b, M_c (μB, along a/|a|, b/|b|, c/|c|) |"
            s_line += "---|"
        ls_out += ["", s_head, s_line]
        for i_out in range(min(n_atom, max_atoms)):
            xyz = atoms["fract_xyz"][:, i_out]
            line = (f"| {atoms['label'][i_out]:} | "
                    f"{item.atom_site.items[atoms['atom_index'][i_out]].label:} | "
                    f"{xyz[0]:.5f} | {xyz[1]:.5f} | {xyz[2]:.5f} |")
            if moment is not None:
                m = moment[:, i_out]
                line += (" |" if numpy.any(numpy.isnan(m)) else
                         " " + ", ".join(f"{v:.4f}" for v in m) + " |")
            ls_out.append(line)
        if n_atom > max_atoms:
            ls_out.append(f"[... {n_atom - max_atoms:} more atoms are cut ...]")
        ls_out.append("")
    return "\n".join(ls_out)


def tool_convert_to_p1(cryspy_object, crystal: str = "", name: str = "") -> str:
    try:
        l_crystal = take_named_crystals(cryspy_object, crystal)
    except ValueError as error:
        return f"Error: {error:}"
    if len(l_crystal) > 1:
        return ("Error: the object has several crystals ("
                + ", ".join(item.data_name for item in l_crystal)
                + "); give the crystal to convert.")
    item = l_crystal[0]
    data_name = name.strip() if name.strip() != "" else f"{item.data_name}_p1"
    if take_data_block(cryspy_object, data_name) is not None:
        return f"Error: the data block '{data_name:}' exists already; give another name."
    crystal_p1 = item.to_p1(data_name)
    cryspy_object.add_items([crystal_p1])
    ls_out = [f"Crystal {data_name:} in P1 is added: {len(crystal_p1.atom_site.items):} atoms "
              f"(from {len(item.atom_site.items):} atoms of {item.data_name:})."]
    for loop_name in ("atom_site_aniso", "atom_site_moment", "atom_site_susceptibility",
                      "atom_site_scat"):
        if crystal_p1.is_attribute(loop_name):
            ls_out.append(f"{loop_name:}: {len(getattr(crystal_p1, loop_name).items):} atoms.")
    ls_out.append("Multipoles and local axes are not converted. Experiments still use "
                  f"{item.data_name:}; the new crystal is not in their phases.")
    return "\n".join(ls_out)


def tool_report_ellipsoids(cryspy_object, kind: str = "ADP", crystal: str = "",
                           unit_cell: bool = False, field: float = 1.) -> str:
    try:
        l_crystal = take_named_crystals(cryspy_object, crystal)
        ls_out = []
        for item in l_crystal:
            ls_out.append(f"Crystal {item.data_name:}.")
            ls_out.append(item.report_ellipsoids(kind, flag_unit_cell=unit_cell, field=field))
            ls_out.append("")
    except ValueError as error:
        return f"Error: {error:}"
    return "\n".join(ls_out)


def tool_add_crystal_from_cod(cryspy_object, formula: str = "",
                              cod_id: str = "", name: str = "") -> str:
    from cryspy.H_functions_global.cod import add_crystal_from_cod
    if formula.strip() == "" and str(cod_id).strip() == "":
        raise ArgumentsNeeded([{
            "name": "formula", "type": "string",
            "question": "the chemical formula of the crystal or its COD number",
            "hint": "e.g. TiO2 or 9008460"}])
    if formula.strip().isdigit() and str(cod_id).strip() == "":  # the COD number as the answer
        formula, cod_id = "", formula.strip()
    try:
        crystal = add_crystal_from_cod(
            cryspy_object, formula=formula, cod_id=str(cod_id),
            data_name=name)
    except ValueError as error:
        return f"Error: {error:}"
    cell = crystal.cell
    return (f"Crystal {crystal.data_name:} is added from COD: space group "
            f"{crystal.space_group.name_hm_alt:}, a, b, c = {cell.length_a:g}, "
            f"{cell.length_b:g}, {cell.length_c:g} Å, "
            f"{len(crystal.atom_site.items):} atoms.")
