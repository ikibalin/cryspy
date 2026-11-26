"""Multipole (Hansen-Coppens) description of the magnetic form factor of an
atom: rows of atom_rho_multipole and atom_local_axes, with dummy atoms for
the local axes.

The same rows serve magnetic space groups and magnetic superspace groups
(propagation vectors): the form factor of a generated atom is rotated by the
operation in both cases
(``calc_magnetic_form_factor_with_multipoles_by_dictionary``).
"""
import re
from typing import Optional, Tuple

import numpy

from cryspy.A_functions_base.function_1_strings import split_list_text
from cryspy.A_functions_base.knowledge_tables import take_magnetic_ions
from cryspy.C_item_loop_classes.cl_1_atom_local_axes import (
    AtomLocalAxes, AtomLocalAxesL)
from cryspy.C_item_loop_classes.cl_1_atom_rho_multipole import (
    AtomRhoMultipole, AtomRhoMultipoleL)
from cryspy.C_item_loop_classes.cl_1_atom_site import AtomSite
from cryspy.exceptions import ArgumentsNeeded, CrysPyError

# distance (Å) of the dummy atoms which give the directions of the local axes
DUMMY_DISTANCE = 0.5
# a dummy atom put before by the user nearer than this (Å) gives a local axis
MAX_DUMMY_DISTANCE = 1.0
# cos of the angle below which a direction is taken as parallel to c (~5°)
COS_PARALLEL = 0.996


def take_element_of(type_symbol: str) -> str:
    match = re.match(r"\s*\d*([A-Za-z][a-z]?)", str(type_symbol))
    return match.group(1).capitalize() if match else ""


def parse_values_of_atoms(text: str, l_label: list) -> dict:
    """{label: value} of 'Mn: Mn3+, F: O2-'; a value without a label
    ('Mn3+') is for every atom. Labels and values have no spaces."""
    d_value = {}
    for part in re.split(r"[,;]", str(text)):
        part = part.strip()
        if part == "":
            continue
        match = re.fullmatch(r"(\S+)\s*[:=]\s*(\S+)", part)
        if match is None:
            d_value.update({label: part for label in l_label if label not in d_value})
            continue
        label, value = match.groups()
        if label not in l_label:
            raise CrysPyError(f"'{part}': no atom {label} among {', '.join(l_label)}.")
        d_value[label] = value
    return d_value


def take_magnetic_ion_of_atom(crystal, label: str, magnetic_ion: str = "") -> Tuple[str, Optional[dict]]:
    """Ion of the atom and, when it is not sure, the choice to ask.

    The ion is ``magnetic_ion``, else the one of atom_site_scat, else the
    type_symbol of atom_site when it has a charge (Cu2+). The ion of
    atom_site_scat is not sure when its element differs from the element
    of the atom (O2- for F): then the choice {"label", "options",
    "proposal", "reason"} is returned with it."""
    ion = magnetic_ion.strip()
    if ion != "":
        return ion, None
    type_symbol = str(crystal.atom_site[label].type_symbol)
    site_ion = type_symbol if ("+" in type_symbol or "-" in type_symbol) else ""
    if crystal.is_attribute("atom_site_scat") and \
            label in list(crystal.atom_site_scat.label):
        ion = str(crystal.atom_site_scat[label].type_symbol).strip()
    if ion != "" and take_element_of(ion) != take_element_of(type_symbol):
        l_option = [ion] + ([site_ion] if site_ion else [])
        return ion, {"label": label, "options": l_option, "proposal": ion,
                     "reason": f"{ion} of atom_site_scat (its dipole form factor), "
                               f"but the atom is {type_symbol}"}
    if ion == "":
        ion = site_ion
    if ion == "":
        raise CrysPyError(
            f"The ion of the atom {label} is unknown: give magnetic_ion "
            "(as Cu2+ or Fe3+).")
    return ion, None


def take_orbitals_of_element(element: str) -> list:
    """Orbitals of the element with radial Slater functions in CrysPy."""
    from cryspy.A_functions_base.charge_form_factor import DATABASE
    return [orbital for name, orbital in DATABASE["Orbitals by radial Slater functions"]
            if name == element]


def take_orbital_of_ion(ion: str) -> Tuple[str, Optional[dict]]:
    """Orbital of the radial function of an ion and, when it is not sure,
    the choice to ask.

    The open shell of a magnetic ion of the knowledge base (Cu2+ 3d9 ->
    3d) is sure. For another ion (an anion with a transferred spin, as
    O2- or F1-) the orbitals with Slater functions in CrysPy are the
    options and the outer shell of its configuration (2p) the proposal."""
    from cryspy.A_functions_base.charge_form_factor import get_atom_name_ion_charge_shell
    d_ion = take_magnetic_ions().get(ion)
    if d_ion is not None:
        match = re.match(r"\s*(\d[spdf])", d_ion["configuration"])
        if match is not None:
            return match.group(1), None
    try:
        element, _, configuration = get_atom_name_ion_charge_shell(ion)
    except Exception:  # noqa: BLE001 (any unknown symbol: the orbital is asked)
        element, configuration = take_element_of(ion), ""
    l_orbital = take_orbitals_of_element(element)
    if len(l_orbital) == 0:
        raise CrysPyError(f"No radial (Slater) functions of {element} in CrysPy: "
                          "the multipoles of this atom cannot be described.")
    l_shell = [re.match(r"(\d[spdf])", item).group(1) for item in configuration.split()
               if re.match(r"(\d[spdf])", item)]
    outer = l_shell[-1] if l_shell and l_shell[-1] in l_orbital else ""
    reason = f"{ion} is not a magnetic ion of the knowledge base"
    if outer:
        reason += f"; {outer} is its outer shell ({configuration})"
    return outer, {"options": l_orbital, "proposal": outer, "reason": reason}


def find_dummy_atom_near(crystal, label: str) -> Optional[Tuple[str, bool]]:
    """(label, along c) of the dummy atom of atom_site nearest to the atom
    (within MAX_DUMMY_DISTANCE, lattice translations taken into account),
    None when there is none. A dummy atom put there before gives the
    direction of a local axis."""
    from cryspy.A_functions_base.unit_cell import calc_m_m_by_unit_cell_parameters
    m_m = calc_m_m_by_unit_cell_parameters(crystal.cell.get_unit_cell_parameters())[0]
    m_m = m_m.reshape(3, 3)
    atom = crystal.atom_site[label]
    xyz = numpy.array([atom.fract_x, atom.fract_y, atom.fract_z], dtype=float)
    l_found = []
    for item in crystal.atom_site.items:
        if str(item.type_symbol).strip().lower() != "dummy":
            continue
        diff = numpy.array([item.fract_x, item.fract_y, item.fract_z], dtype=float) - xyz
        diff -= numpy.round(diff)
        vector = m_m @ diff
        distance = float(numpy.linalg.norm(vector))
        if 1e-3 < distance <= MAX_DUMMY_DISTANCE:
            along_c = abs(vector[2]) / distance > COS_PARALLEL
            l_found.append((distance, str(item.label), along_c))
    if len(l_found) == 0:
        return None
    _, dummy_label, along_c = min(l_found)
    return dummy_label, along_c


def parse_multipole_coefficients(coefficients: str) -> dict:
    """{'coeff_P20': 0.1, 'coeff_P2m2': -0.05} from 'P20 0.1, P2-2=-0.05'."""
    d_coeff = {}
    for part in re.split(r"[,;]", coefficients):
        if part.strip() == "":
            continue
        match = re.fullmatch(
            r"\s*(?:coeff_)?P(\d)(-|m)?(\d)\s*[= ]\s*([-+0-9.eE]+)\s*", part)
        if match is None:
            raise CrysPyError(
                f"Cannot read the multipole coefficient '{part.strip()}': write "
                "it as 'P20 0.1' or 'P2-2=0.05'.")
        l_index, minus, m_index, value = match.groups()
        if int(m_index) > int(l_index) or int(l_index) > 4:
            raise CrysPyError(
                f"There is no multipole P{l_index}{minus or ''}{m_index}: "
                "l = 0..4, |m| <= l.")
        name = f"coeff_P{l_index}{'m' if minus else ''}{m_index}"
        d_coeff[name] = float(value)
    return d_coeff


def form_dummy_atoms(crystal, label: str, axes: tuple = ("c", "a")):
    """Dummy atoms near the atom along the given cell axes ("c", "a").

    With ``ax1 = Z`` along c and ``ax2 = X`` along a (orthogonalised to c),
    the local axes are the Cartesian axes of CrysPy (x || a*, z || c).
    """
    atom = crystal.atom_site[label]
    unit_cell_parameters = crystal.cell.get_unit_cell_parameters()
    xyz = numpy.array([atom.fract_x, atom.fract_y, atom.fract_z], dtype=float)
    l_dummy = []
    for axis in axes:
        i_axis = "abc".index(axis)
        position = xyz.copy()
        position[i_axis] += DUMMY_DISTANCE / unit_cell_parameters[i_axis]
        l_dummy.append(AtomSite(
            label=f"D{label}_{axis}", type_symbol="Dummy",
            fract_x=float(numpy.round(position[0], 8)),
            fract_y=float(numpy.round(position[1], 8)),
            fract_z=float(numpy.round(position[2], 8)),
            occupancy=1., adp_type="Biso", b_iso_or_equiv=0.))
    return l_dummy


def add_multipole_form_factor(
        crystal, atom_label: str, magnetic_ion: str = "", orbital: str = "",
        coefficients: str = "", kappa: float = 1.) -> str:
    """Body of Crystal.add_multipole_form_factor, see its docstring."""
    l_label = split_list_text(atom_label)
    l_site = list(crystal.atom_site.label)
    l_unknown = [label for label in l_label if label not in l_site]
    if len(l_label) == 0 or l_unknown:
        raise CrysPyError(
            f"There is no atom '{', '.join(l_unknown) or atom_label}' in atom_site; "
            "the atoms are: " + ", ".join(l_site) + ".")
    l_with = list(crystal.atom_rho_multipole.atom_label) \
        if crystal.is_attribute("atom_rho_multipole") else []
    l_skip = [label for label in l_label if label in l_with]
    l_label = [label for label in l_label if label not in l_with]
    if len(l_label) == 0:
        raise CrysPyError(
            f"The atom{'s' if len(l_skip) > 1 else ''} {', '.join(l_skip)} "
            "has multipoles already; change its row of atom_rho_multipole "
            "or remove it first.")
    d_coeff = {"coeff_P00": 1.}
    d_coeff.update(parse_multipole_coefficients(coefficients))

    # every value is found before anything is changed: one round of
    # questions when an ion or an orbital is not sure
    d_ion_given = parse_values_of_atoms(magnetic_ion, l_label)
    d_orbital_given = parse_values_of_atoms(orbital, l_label)
    d_ion, l_ion_doubt = {}, []
    for label in l_label:
        d_ion[label], d_doubt = take_magnetic_ion_of_atom(
            crystal, label, d_ion_given.get(label, ""))
        if d_doubt is not None:
            l_ion_doubt.append(d_doubt)
    if l_ion_doubt:
        raise ArgumentsNeeded([form_question_of_atoms(
            "magnetic_ion", "the ion of the radial function of", l_ion_doubt,
            d_ion_given, l_label)])
    d_orbital, l_orbital_doubt = {}, []
    for label in l_label:
        if d_orbital_given.get(label, "").strip():
            d_orbital[label] = d_orbital_given[label].strip()
            continue
        d_orbital[label], d_doubt = take_orbital_of_ion(d_ion[label])
        if d_doubt is not None:
            l_orbital_doubt.append(dict(d_doubt, label=label))
    if l_orbital_doubt:
        raise ArgumentsNeeded([form_question_of_atoms(
            "orbital", "the orbital of the radial function of", l_orbital_doubt,
            d_orbital_given, l_label)])

    ls_out = []
    for label in l_label:
        ls_out.append(add_multipole_of_atom(
            crystal, label, d_ion[label], d_orbital[label], d_coeff, float(kappa)))
    if l_skip:
        ls_out.append(f"{', '.join(l_skip)}: multipoles already there, not changed.")
    ls_out.append(
        "The multipoles replace the dipole form factor of the atom for magnetic "
        "space groups and propagation vectors (superspace groups) alike. "
        "Give only the coefficients allowed by the site symmetry: CrysPy "
        "does not check them; refine kappa and the coefficients against flip ratios.")
    return "\n\n".join(ls_out)


def form_question_of_atoms(name: str, question: str, l_doubt: list, d_given: dict,
                           l_label: list) -> dict:
    """One question of a parameter for the atoms in doubt; the answer is
    'label: value' for each of them (values known already are kept)."""
    l_doubt_label = [d_doubt["label"] for d_doubt in l_doubt]
    d_known = {label: value for label, value in d_given.items() if label not in l_doubt_label}

    def form_answer(d_value: dict) -> str:
        d_all = dict(d_known, **d_value)
        return ", ".join(f"{label}: {d_all[label]}" for label in l_label if label in d_all)

    d_question = {
        "name": name, "type": "string",
        "question": f"{question} {', '.join(l_doubt_label)}",
        "hint": "; ".join(f"{d_doubt['label']}: {', '.join(d_doubt['options'])} "
                          f"({d_doubt['reason']})" for d_doubt in l_doubt)
        + "; answer as 'label: value'"}
    if all(d_doubt["proposal"] for d_doubt in l_doubt):
        d_question["proposal"] = form_answer(
            {d_doubt["label"]: d_doubt["proposal"] for d_doubt in l_doubt})
        d_question["s_proposal"] = "; ".join(
            f"{d_doubt['label']}: {d_doubt['proposal']} ({d_doubt['reason']})"
            for d_doubt in l_doubt)
    if len(l_doubt) == 1:
        d_question["options"] = [form_answer({l_doubt[0]["label"]: option})
                                 for option in l_doubt[0]["options"]]
    return d_question


def add_multipole_of_atom(crystal, label: str, ion: str, orbital: str, d_coeff: dict,
                          kappa: float) -> str:
    """Rows of atom_rho_multipole and, when missing, atom_local_axes of one
    atom; a dummy atom near the atom (put there before) is used for an axis."""
    item_multipole = AtomRhoMultipole(
        atom_label=label, atom_type_symbol=ion, orbital=orbital,
        kappa=kappa, **d_coeff)
    try:
        item_multipole.get_n_zeta_coeff()
    except KeyError:
        raise CrysPyError(
            f"No radial (Slater) function of the orbital {orbital} of {ion} "
            "in the CrysPy database; give another orbital.") from None

    ls_out = [f"# Multipole form factor of {label}", "",
              f"Ion {ion}, orbital {orbital}, kappa {kappa:g}; "
              + ", ".join(f"{name[6:].replace('m', '-')} = {value:g}"
                          for name, value in d_coeff.items()) + "."]
    flag_axes = crystal.is_attribute("atom_local_axes") and \
        label in list(crystal.atom_local_axes.atom_label)
    if flag_axes:
        ls_out.append(f"The local axes of {label} are kept.")
    else:
        ls_out.append(add_local_axes(crystal, label))
    if crystal.is_attribute("atom_rho_multipole"):
        crystal.atom_rho_multipole.items.append(item_multipole)
    else:
        loop = AtomRhoMultipoleL()
        loop.items = [item_multipole]
        crystal.add_items([loop])

    flag_magnetic = crystal.is_attribute("atom_site_scat") and \
        label in list(crystal.atom_site_scat.label)
    if not flag_magnetic:
        ls_out.append(
            f"{label} has no atom_site_scat row: the multipoles are used for "
            "magnetic atoms only (add the magnetic ion first).")
    return "\n".join(ls_out)


def add_local_axes(crystal, label: str) -> str:
    """Row of atom_local_axes of the atom: z || c, x || a* by new dummy
    atoms, or the direction of a dummy atom near the atom for one axis."""
    found = find_dummy_atom_near(crystal, label)
    if found is None:
        l_dummy = form_dummy_atoms(crystal, label, ("c", "a"))
        crystal.atom_site.items.extend(l_dummy)
        item_axes = AtomLocalAxes(
            atom_label=label, atom0=l_dummy[0].label, ax1="Z",
            atom1=label, atom2=l_dummy[1].label, ax2="X")
        text = (f"Local axes: z || c, x || a* (dummy atoms {l_dummy[0].label}, "
                f"{l_dummy[1].label}); change atom_local_axes to align them with "
                "the bonds.")
    else:
        dummy, along_c = found
        new_axis = "a" if along_c else "c"
        new_dummy = form_dummy_atoms(crystal, label, (new_axis, ))[0]
        crystal.atom_site.items.append(new_dummy)
        if along_c:
            item_axes = AtomLocalAxes(atom_label=label, atom0=dummy, ax1="Z",
                                      atom1=label, atom2=new_dummy.label, ax2="X")
            text = (f"Local axes: z to the dummy atom {dummy} (it was there, along c), "
                    f"x || a* (new dummy atom {new_dummy.label}).")
        else:
            item_axes = AtomLocalAxes(atom_label=label, atom0=dummy, ax1="X",
                                      atom1=label, atom2=new_dummy.label, ax2="Z")
            text = (f"Local axes: x to the dummy atom {dummy} (it was there), z || c "
                    f"orthogonalised to x (new dummy atom {new_dummy.label}); change "
                    "atom_local_axes if another choice is meant.")
    if crystal.is_attribute("atom_local_axes"):
        crystal.atom_local_axes.items.append(item_axes)
    else:
        loop = AtomLocalAxesL()
        loop.items = [item_axes]
        crystal.add_items([loop])
    return text
