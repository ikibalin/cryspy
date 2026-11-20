"""Distances, bond angles and torsion angles of the atoms of ``atom_site``.

Reports in Markdown and the CIF loops ``_geom_bond``, ``_geom_angle``,
``_geom_torsion`` (International Union of Crystallography core CIF
dictionary) made from the cell, the atoms and the symmetry elements of a
crystal. The numbers come from
:mod:`cryspy.A_functions_base.structure_geometry`; nothing is stored.
"""
import itertools
import re
from typing import List

import numpy

from cryspy.A_functions_base.function_1_strings import value_error_mark_to_string
from cryspy.A_functions_base.structure_geometry import (
    calc_geometry_sigma, calc_geometry_value, calc_images_near_position,
    calc_restraint_difference, calc_restraint_value, find_chain_positions,
    form_position, form_symmetry_code, identity_index)
from cryspy.C_item_loop_classes.cl_1_geom import (
    GeomAngle, GeomAngleL, GeomBond, GeomBondL, GeomHbond, GeomHbondL,
    GeomTorsion, GeomTorsionL)
from cryspy.C_item_loop_classes.cl_1_restr import (
    RestrAngle, RestrAngleL, RestrDistance, RestrDistanceL, RestrTorsion,
    RestrTorsionL)

GEOMETRY_KINDS = {2: ("distance", "Å"), 3: ("angle", "°"), 4: ("torsion angle", "°")}


class GeometryInput:
    """Cell, atoms and symmetry elements of a crystal for the geometry."""

    def __init__(self, cell, atom_site, full_symm_elems):
        self.unit_cell_parameters = cell.get_unit_cell_parameters()
        sigma = [getattr(cell, name, None) for name in (
            "length_a_sigma", "length_b_sigma", "length_c_sigma",
            "angle_alpha_sigma", "angle_beta_sigma", "angle_gamma_sigma")]
        sigma = numpy.array([0. if value is None else float(value) for value in sigma])
        sigma[3:] = numpy.radians(sigma[3:])
        self.unit_cell_parameters_sigma = numpy.nan_to_num(sigma)
        self.labels = [str(item.label) for item in atom_site.items]
        self.fract_xyz = numpy.array(
            [atom_site.fract_x, atom_site.fract_y, atom_site.fract_z], dtype=float)
        self.fract_xyz_sigma = numpy.nan_to_num(numpy.array(
            [[0. if getattr(item, f"fract_{axis}_sigma", None) is None
              else float(getattr(item, f"fract_{axis}_sigma"))
              for item in atom_site.items] for axis in "xyz"], dtype=float))
        self.full_symm_elems = full_symm_elems

    def take_atom_index(self, label: str) -> int:
        if label in self.labels:
            return self.labels.index(label)
        l_lower = [name.lower() for name in self.labels]
        if label.lower() in l_lower:
            return l_lower.index(label.lower())
        raise ValueError(f"atom '{label:}' is not in atom_site; the atoms are "
                         + ", ".join(self.labels) + ".")

    def split_labels(self, atoms) -> List[str]:
        """Labels from 'N1 C1', 'N1, C1, C2', 'N1-C1-C2' or a list."""
        if isinstance(atoms, (list, tuple)):
            return [str(label).strip() for label in atoms if str(label).strip()]
        l_label = [label for label in re.split(r"[\s,;]+", str(atoms).strip()) if label]
        if len(l_label) == 1 and "-" in l_label[0] and l_label[0] not in self.labels:
            l_label = [label for label in l_label[0].split("-") if label]
        return l_label

    def value_sigma(self, l_position) -> tuple:
        value = calc_geometry_value(self.unit_cell_parameters, self.fract_xyz, l_position)
        sigma = calc_geometry_sigma(
            self.unit_cell_parameters, self.unit_cell_parameters_sigma,
            self.fract_xyz, self.fract_xyz_sigma, l_position)
        return value, sigma

    def position(self, i_atom, symm_index, translation):
        return form_position(self.full_symm_elems, i_atom, symm_index, translation)

    def symmetry_code(self, symm_index, translation) -> str:
        return form_symmetry_code(self.full_symm_elems, symm_index, translation)

    def cif_symmetry_code(self, symm_index, translation) -> str:
        """'.' for the atom itself, otherwise 'n_klm' of the core CIF
        dictionary (n is the number of the operation in the loop
        _space_group_symop written with the geometry, k = 5 + translation)."""
        if symm_index == identity_index(self.full_symm_elems) and not numpy.any(translation):
            return "."
        return f"{symm_index + 1:}_" + "".join(str(5 + int(t)) for t in translation)


def format_value_sigma(value: float, sigma: float, n_digits: int) -> str:
    if sigma > 0.:
        return value_error_mark_to_string(value, sigma, "")
    return f"{value:.{n_digits:}f}"


def calc_geometry(geometry: GeometryInput, atoms) -> dict:
    """Distance, angle or torsion angle of 2, 3 or 4 atoms of atom_site
    (see :meth:`Crystal.calc_distance_or_angle`)."""
    l_label = geometry.split_labels(atoms)
    if len(l_label) not in GEOMETRY_KINDS:
        raise ValueError("give 2 atoms for a distance, 3 for an angle (the vertex "
                         f"in the middle) or 4 for a torsion angle; got {len(l_label):}.")
    l_index = [geometry.take_atom_index(label) for label in l_label]
    l_chain = find_chain_positions(geometry.unit_cell_parameters, geometry.full_symm_elems,
                                   geometry.fract_xyz, l_index)
    value, sigma = geometry.value_sigma([geometry.position(*item) for item in l_chain])
    kind, unit = GEOMETRY_KINDS[len(l_label)]
    return {"kind": kind, "value": value, "sigma": sigma, "unit": unit,
            "labels": [geometry.labels[i_atom] for i_atom in l_index],
            "symmetry": [geometry.symmetry_code(symm_index, translation)
                         for _, symm_index, translation in l_chain],
            "cif_symmetry": [geometry.cif_symmetry_code(symm_index, translation)
                             for _, symm_index, translation in l_chain]}


def find_neighbours(geometry: GeometryInput, i_center: int, max_distance: float) -> List[dict]:
    """Atoms (equivalent positions of all atoms of atom_site) within
    max_distance of the atom i_center, nearest first."""
    l_neighbour = []
    center = geometry.fract_xyz[:, i_center]
    for i_atom in range(len(geometry.labels)):
        images = calc_images_near_position(
            geometry.unit_cell_parameters, geometry.full_symm_elems,
            geometry.fract_xyz[:, i_atom], center, max_distance)
        for i_pos in range(images["distance"].size):
            l_neighbour.append({
                "atom_index": i_atom, "symm_index": int(images["symm_index"][i_pos]),
                "translation": images["translation"][:, i_pos],
                "distance": float(images["distance"][i_pos])})
    l_neighbour.sort(key=lambda item: (item["distance"], item["atom_index"]))
    return l_neighbour


def calc_bonds_and_angles(geometry: GeometryInput, l_center: List[int],
                          max_distance: float, flag_angles: bool = True) -> tuple:
    """Bonds (distances up to max_distance) and the angles between them at
    the atoms l_center; each bond between two different atoms of
    atom_site once."""
    i_identity = identity_index(geometry.full_symm_elems)
    zero = numpy.zeros(3, dtype=int)
    l_bond, l_angle = [], []
    flag_all = len(l_center) == len(geometry.labels)
    for i_center in l_center:
        l_neighbour = find_neighbours(geometry, i_center, max_distance)
        position_center = geometry.position(i_center, i_identity, zero)
        for neighbour in l_neighbour:
            if flag_all and neighbour["atom_index"] < i_center:
                continue  # given at that atom
            position = geometry.position(neighbour["atom_index"], neighbour["symm_index"],
                                         neighbour["translation"])
            value, sigma = geometry.value_sigma([position_center, position])
            l_bond.append({"labels": [geometry.labels[i_center],
                                      geometry.labels[neighbour["atom_index"]]],
                           "value": value, "sigma": sigma,
                           "neighbours": [neighbour]})
        if not flag_angles:
            continue
        for neighbour_1, neighbour_3 in itertools.combinations(l_neighbour, 2):
            l_position = [
                geometry.position(neighbour_1["atom_index"], neighbour_1["symm_index"],
                                  neighbour_1["translation"]),
                position_center,
                geometry.position(neighbour_3["atom_index"], neighbour_3["symm_index"],
                                  neighbour_3["translation"])]
            value, sigma = geometry.value_sigma(l_position)
            l_angle.append({"labels": [geometry.labels[neighbour_1["atom_index"]],
                                       geometry.labels[i_center],
                                       geometry.labels[neighbour_3["atom_index"]]],
                            "value": value, "sigma": sigma,
                            "neighbours": [neighbour_1, neighbour_3]})
    return l_bond, l_angle


def form_code(geometry: GeometryInput, neighbour: dict, flag_cif: bool) -> str:
    if flag_cif:
        return geometry.cif_symmetry_code(neighbour["symm_index"], neighbour["translation"])
    return geometry.symmetry_code(neighbour["symm_index"], neighbour["translation"])


def report_chain(geometry: GeometryInput, atoms, max_distance: float) -> str:
    """Markdown text of one distance, angle or torsion angle; for two atoms
    also the other distances between them up to max_distance."""
    result = calc_geometry(geometry, atoms)
    n_digits = 4 if result["kind"] == "distance" else 2
    s_value = format_value_sigma(result["value"], result["sigma"], n_digits)
    s_atoms = "–".join(result["labels"])
    ls_out = [f"{result['kind'].capitalize():} {s_atoms:} = {s_value:} {result['unit']:}."]
    ls_out.append("Atoms: " + "; ".join(
        f"{label:} at {code:}" for label, code in zip(result["labels"], result["symmetry"])) + ".")
    if result["kind"] == "distance":
        i_1, i_2 = (geometry.take_atom_index(label) for label in result["labels"])
        images = calc_images_near_position(
            geometry.unit_cell_parameters, geometry.full_symm_elems,
            geometry.fract_xyz[:, i_2], geometry.fract_xyz[:, i_1], max_distance)
        if images["distance"].size > 1:
            ls_out.append(f"All {result['labels'][1]:} atoms within {max_distance:g} Å "
                          f"of {result['labels'][0]:}:")
            for i_pos in range(images["distance"].size):
                position = geometry.position(i_2, images["symm_index"][i_pos],
                                             images["translation"][:, i_pos])
                value, sigma = geometry.value_sigma(
                    [geometry.position(i_1, identity_index(geometry.full_symm_elems),
                                       numpy.zeros(3, dtype=int)), position])
                code = geometry.symmetry_code(images["symm_index"][i_pos],
                                              images["translation"][:, i_pos])
                ls_out.append(f"- {format_value_sigma(value, sigma, 4):} Å ({code:})")
    if numpy.all(geometry.fract_xyz_sigma == 0.) and numpy.all(
            geometry.unit_cell_parameters_sigma == 0.):
        ls_out.append("No standard uncertainties: the coordinates and the cell have no sigma.")
    return "\n".join(ls_out)


def report_bonds_and_angles(geometry: GeometryInput, l_center: List[int],
                            max_distance: float, flag_angles: bool = True,
                            flag_cif: bool = False) -> str:
    """Markdown tables (or CIF loops) of the bonds and angles at the atoms."""
    l_bond, l_angle = calc_bonds_and_angles(geometry, l_center, max_distance, flag_angles)
    if flag_cif:
        return form_geometry_cif(geometry, l_bond, l_angle)
    s_atoms = ", ".join(geometry.labels[i_atom] for i_atom in l_center)
    ls_out = [f"Bonds up to {max_distance:g} Å of {s_atoms:} "
              f"({len(l_bond):} distances, {len(l_angle):} angles).", "",
              "| atom 1 | atom 2 | distance, Å | atom 2 at |", "|---|---|---|---|"]
    for bond in l_bond:
        ls_out.append(f"| {bond['labels'][0]:} | {bond['labels'][1]:} | "
                      f"{format_value_sigma(bond['value'], bond['sigma'], 4):} | "
                      f"{form_code(geometry, bond['neighbours'][0], False):} |")
    if flag_angles and l_angle:
        ls_out += ["", "| atom 1 | vertex | atom 3 | angle, ° | atom 1 at | atom 3 at |",
                   "|---|---|---|---|---|---|"]
        for angle in l_angle:
            ls_out.append(
                f"| {angle['labels'][0]:} | {angle['labels'][1]:} | {angle['labels'][2]:} | "
                f"{format_value_sigma(angle['value'], angle['sigma'], 2):} | "
                f"{form_code(geometry, angle['neighbours'][0], False):} | "
                f"{form_code(geometry, angle['neighbours'][1], False):} |")
    return "\n".join(ls_out)


def form_symop_cif(geometry: GeometryInput) -> List[str]:
    ls_out = ["loop_", "_space_group_symop.id", "_space_group_symop.operation_xyz"]
    for symm_index in range(numpy.asarray(geometry.full_symm_elems).shape[1]):
        ls_out.append(f"{symm_index + 1:} {geometry.symmetry_code(symm_index, (0, 0, 0)):}")
    return ls_out


def form_geometry_cif(geometry: GeometryInput, l_bond: List[dict], l_angle: List[dict]) -> str:
    """CIF loops _geom_bond and _geom_angle with the operations they refer
    to (the symmetry codes 'n_klm' number the operations of this
    _space_group_symop loop, centring included)."""
    ls_out = form_symop_cif(geometry)
    ls_out += ["", "loop_", "_geom_bond.atom_site_label_1", "_geom_bond.atom_site_label_2",
               "_geom_bond.distance", "_geom_bond.site_symmetry_2"]
    for bond in l_bond:
        ls_out.append(f"{bond['labels'][0]:} {bond['labels'][1]:} "
                      f"{format_value_sigma(bond['value'], bond['sigma'], 4):} "
                      f"{form_code(geometry, bond['neighbours'][0], True):}")
    if l_angle:
        ls_out += ["", "loop_", "_geom_angle.atom_site_label_1",
                   "_geom_angle.atom_site_label_2", "_geom_angle.atom_site_label_3",
                   "_geom_angle.value", "_geom_angle.site_symmetry_1",
                   "_geom_angle.site_symmetry_3"]
        for angle in l_angle:
            ls_out.append(
                f"{angle['labels'][0]:} {angle['labels'][1]:} {angle['labels'][2]:} "
                f"{format_value_sigma(angle['value'], angle['sigma'], 2):} "
                f"{form_code(geometry, angle['neighbours'][0], True):} "
                f"{form_code(geometry, angle['neighbours'][1], True):}")
    return "\n".join(ls_out)


def form_chain_cif(geometry: GeometryInput, result: dict) -> str:
    """CIF loop of one distance, angle or torsion angle of calc_geometry."""
    labels, codes = result["labels"], result["cif_symmetry"]
    n_digits = 4 if result["kind"] == "distance" else 2
    s_value = format_value_sigma(result["value"], result["sigma"], n_digits)
    if result["kind"] == "distance":
        names = ["_geom_bond.atom_site_label_1", "_geom_bond.atom_site_label_2",
                 "_geom_bond.distance", "_geom_bond.site_symmetry_1",
                 "_geom_bond.site_symmetry_2"]
    elif result["kind"] == "angle":
        names = [f"_geom_angle.atom_site_label_{i:}" for i in (1, 2, 3)] + [
            "_geom_angle.value"] + [f"_geom_angle.site_symmetry_{i:}" for i in (1, 2, 3)]
    else:
        names = [f"_geom_torsion.atom_site_label_{i:}" for i in (1, 2, 3, 4)] + [
            "_geom_torsion.value"] + [f"_geom_torsion.site_symmetry_{i:}" for i in (1, 2, 3, 4)]
    ls_out = form_symop_cif(geometry) + ["", "loop_"] + names
    ls_out.append(" ".join(labels + [s_value] + codes))
    return "\n".join(ls_out)



def take_element(type_symbol: str) -> str:
    """Element of a type symbol: 'O2-' -> 'O', 'Fe3+' -> 'Fe', 'D' -> 'D'."""
    match = re.match(r"[A-Za-z][a-z]?", str(type_symbol).strip())
    return "" if match is None else match.group(0).capitalize()


def find_hbonds(geometry: GeometryInput, type_symbols: List[str],
                acceptors: str = "N O F Cl S", max_distance_dh: float = 1.2,
                max_distance_ha: float = 2.6, min_angle_dha: float = 120.) -> List[dict]:
    """Hydrogen bonds D–H···A of the hydrogen (deuterium) atoms of atom_site.

    The donor D is the atom (not H) closest to H within max_distance_dh;
    the acceptors are the atoms of the elements ``acceptors`` within
    max_distance_ha of H, apart from the donor, with the angle D–H···A not
    smaller than min_angle_dha (geometric criteria of Jeffrey, An
    Introduction to Hydrogen Bonding, 1997; give other limits for weak
    bonds).

    Returns
    -------
    list of dict
        ``labels`` (D, H, A), ``neighbours`` (D and A positions as in
        :func:`find_neighbours`), ``values`` and ``sigmas`` of D–H, H···A,
        D···A (Å) and D–H···A (degrees).
    """
    l_element = [take_element(symbol) for symbol in type_symbols]
    s_acceptor = {take_element(symbol) for symbol in acceptors.replace(",", " ").split()}
    i_identity = identity_index(geometry.full_symm_elems)
    zero = numpy.zeros(3, dtype=int)
    l_hbond = []
    for i_h, element in enumerate(l_element):
        if element not in ("H", "D"):
            continue
        l_neighbour = find_neighbours(geometry, i_h, max(max_distance_dh, max_distance_ha))
        l_donor = [item for item in l_neighbour if l_element[item["atom_index"]] not in ("H", "D")
                   and item["distance"] <= max_distance_dh]
        if len(l_donor) == 0:
            continue
        donor = l_donor[0]
        position_d = geometry.position(donor["atom_index"], donor["symm_index"], donor["translation"])
        position_h = geometry.position(i_h, i_identity, zero)
        for acceptor in l_neighbour:
            if (acceptor is donor or l_element[acceptor["atom_index"]] not in s_acceptor
                    or acceptor["distance"] > max_distance_ha):
                continue
            position_a = geometry.position(acceptor["atom_index"], acceptor["symm_index"],
                                           acceptor["translation"])
            l_value_sigma = [geometry.value_sigma(l_position) for l_position in (
                [position_d, position_h], [position_h, position_a], [position_d, position_a],
                [position_d, position_h, position_a])]
            if l_value_sigma[3][0] < min_angle_dha:
                continue
            l_hbond.append({
                "labels": [geometry.labels[donor["atom_index"]], geometry.labels[i_h],
                           geometry.labels[acceptor["atom_index"]]],
                "neighbours": [donor, acceptor],
                "values": [value for value, _ in l_value_sigma],
                "sigmas": [sigma for _, sigma in l_value_sigma]})
    return l_hbond


def report_hbonds(geometry: GeometryInput, l_hbond: List[dict]) -> str:
    """Markdown table of the hydrogen bonds of :func:`find_hbonds`."""
    if len(l_hbond) == 0:
        return "No hydrogen bonds found."
    ls_out = [f"{len(l_hbond):} hydrogen bonds D–H···A.", "",
              "| D | H | A | D–H, Å | H···A, Å | D···A, Å | D–H···A, ° | D at | A at |",
              "|---|---|---|---|---|---|---|---|---|"]
    for hbond in l_hbond:
        s_values = " | ".join(format_value_sigma(value, sigma, 2 if i_value == 3 else 3)
                              for i_value, (value, sigma) in enumerate(
                                  zip(hbond["values"], hbond["sigmas"])))
        ls_out.append(f"| {' | '.join(hbond['labels']):} | {s_values:} | "
                      f"{form_code(geometry, hbond['neighbours'][0], False):} | "
                      f"{form_code(geometry, hbond['neighbours'][1], False):} |")
    return "\n".join(ls_out)


def form_geometry_loops(geometry: GeometryInput, l_bond: List[dict], l_angle: List[dict],
                        l_torsion: List[dict], l_hbond: List[dict]) -> dict:
    """Loop classes geom_bond, geom_angle, geom_torsion, geom_hbond (only
    the non-empty ones) of the calculated geometry; {attribute name: loop}."""
    d_loop = {}
    if l_bond:
        loop = GeomBondL()
        loop.items = [GeomBond(
            atom_site_label_1=bond["labels"][0], atom_site_label_2=bond["labels"][1],
            distance=format_value_sigma(bond["value"], bond["sigma"], 4),
            site_symmetry_1=".", site_symmetry_2=form_code(geometry, bond["neighbours"][0], True))
            for bond in l_bond]
        d_loop["geom_bond"] = loop
    if l_angle:
        loop = GeomAngleL()
        loop.items = [GeomAngle(
            atom_site_label_1=angle["labels"][0], atom_site_label_2=angle["labels"][1],
            atom_site_label_3=angle["labels"][2],
            value=format_value_sigma(angle["value"], angle["sigma"], 2),
            site_symmetry_1=form_code(geometry, angle["neighbours"][0], True),
            site_symmetry_2=".",
            site_symmetry_3=form_code(geometry, angle["neighbours"][1], True))
            for angle in l_angle]
        d_loop["geom_angle"] = loop
    if l_torsion:
        loop = GeomTorsionL()
        loop.items = [GeomTorsion(
            **{f"atom_site_label_{i + 1:}": label for i, label in enumerate(torsion["labels"])},
            **{f"site_symmetry_{i + 1:}": code for i, code in enumerate(torsion["cif_symmetry"])},
            value=format_value_sigma(torsion["value"], torsion["sigma"], 2))
            for torsion in l_torsion]
        d_loop["geom_torsion"] = loop
    if l_hbond:
        loop = GeomHbondL()
        l_item = []
        for hbond in l_hbond:
            s_values = [format_value_sigma(value, sigma, 2 if i_value == 3 else 3)
                        for i_value, (value, sigma) in enumerate(
                            zip(hbond["values"], hbond["sigmas"]))]
            l_item.append(GeomHbond(
                atom_site_label_d=hbond["labels"][0], atom_site_label_h=hbond["labels"][1],
                atom_site_label_a=hbond["labels"][2],
                distance_dh=s_values[0], distance_ha=s_values[1], distance_da=s_values[2],
                angle_dha=s_values[3],
                site_symmetry_d=form_code(geometry, hbond["neighbours"][0], True),
                site_symmetry_h=".",
                site_symmetry_a=form_code(geometry, hbond["neighbours"][1], True)))
        loop.items = l_item
        d_loop["geom_hbond"] = loop
    return d_loop


RESTRAINT_CLASSES = {2: (RestrDistance, RestrDistanceL, "restr_distance"),
                     3: (RestrAngle, RestrAngleL, "restr_angle"),
                     4: (RestrTorsion, RestrTorsionL, "restr_torsion")}
# sigma of a restraint when the user gives none: Å for distances, degrees
DEFAULT_RESTRAINT_SIGMA = {2: 0.01, 3: 1., 4: 2.}


def take_symmetry_index(geometry: GeometryInput, code: str) -> int:
    """Index of the symmetry element of a site symmetry code 'n_klm' ('.'
    or empty is the atom itself)."""
    code = "." if code is None else str(code).strip()
    if code in (".", "?", "", "None"):
        return identity_index(geometry.full_symm_elems)
    n_symm = numpy.asarray(geometry.full_symm_elems).shape[1]
    try:
        symm_index = int(code.split("_")[0]) - 1
    except ValueError:
        raise ValueError(f"site symmetry '{code:}' is not 'n_klm'.") from None
    if not 0 <= symm_index < n_symm:
        raise ValueError(f"site symmetry '{code:}': the crystal has {n_symm:} "
                         "symmetry operations.")
    return symm_index


def take_restraint_arrays(geometry: GeometryInput, item) -> tuple:
    """Atom indices, rotations and offsets of the atoms of a restraint item
    for :func:`calc_restraint_value`."""
    l_atom_index = [geometry.take_atom_index(label) for label in item.get_labels()]
    symm = numpy.asarray(geometry.full_symm_elems, dtype=float)
    l_rotation, l_offset = [], []
    for code in item.get_site_symmetry():
        symm_index = take_symmetry_index(geometry, code)
        l_rotation.append(symm[4:13, symm_index].reshape(3, 3))
        l_offset.append(symm[0:3, symm_index] / symm[3, symm_index])
    return l_atom_index, l_rotation, l_offset


def calc_restraint_item_value(geometry: GeometryInput, item) -> float:
    return calc_restraint_value(geometry.unit_cell_parameters, geometry.fract_xyz,
                                *take_restraint_arrays(geometry, item))


def form_restraint_item(geometry: GeometryInput, atoms, target=None, sigma: float = 0.):
    """Restraint item (restr_distance, restr_angle or restr_torsion) of 2–4
    atoms: the equivalent atoms of calc_geometry, the target is the
    present value when not given."""
    result = calc_geometry(geometry, atoms)
    n_atom = len(result["labels"])
    item_class = RESTRAINT_CLASSES[n_atom][0]
    if sigma is None or sigma <= 0.:
        sigma = DEFAULT_RESTRAINT_SIGMA[n_atom]
    target = result["value"] if target is None else float(target)
    item = item_class(
        **{f"atom_site_label_{i + 1:}": label for i, label in enumerate(result["labels"])},
        **{f"site_symmetry_{i + 1:}": code for i, code in enumerate(result["cif_symmetry"])},
        target=target, target_weight_param=float(sigma))
    item.diff = calc_restraint_difference(result["value"], target, n_atom)
    return item


def report_restraints(geometry: GeometryInput, l_loop) -> str:
    """Markdown table of the restraints with their present values."""
    ls_out = ["| restraint | atoms | target | sigma | value | diff | (diff/sigma)^2 |",
              "|---|---|---|---|---|---|---|"]
    chi_sq = 0.
    for loop in l_loop:
        for item in loop.items:
            labels = item.get_labels()
            value = calc_restraint_item_value(geometry, item)
            diff = calc_restraint_difference(value, item.target, len(labels))
            weighted = (diff / item.target_weight_param) ** 2
            chi_sq += weighted
            n_digits = 4 if len(labels) == 2 else 2
            ls_out.append(
                f"| {GEOMETRY_KINDS[len(labels)][0]:} | {'–'.join(labels):} | "
                f"{item.target:.{n_digits:}f} | {item.target_weight_param:g} | "
                f"{value:.{n_digits:}f} | {diff:.{n_digits:}f} | {weighted:.2f} |")
    ls_out.append("")
    ls_out.append(f"Sum (diff/sigma)^2 = {chi_sq:.2f} (added to chi_sq in the refinement).")
    return "\n".join(ls_out)
