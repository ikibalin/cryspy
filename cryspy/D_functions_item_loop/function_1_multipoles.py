"""Multipole (Hansen-Coppens) description of the magnetic form factor of an
atom: rows of atom_rho_multipole and atom_local_axes, with dummy atoms for
the local axes.

The same rows serve magnetic space groups and magnetic superspace groups
(propagation vectors): the form factor of a generated atom is rotated by the
operation in both cases
(``calc_magnetic_form_factor_with_multipoles_by_dictionary``).
"""
import re

import numpy

from cryspy.A_functions_base.knowledge_tables import take_magnetic_ions
from cryspy.C_item_loop_classes.cl_1_atom_local_axes import (
    AtomLocalAxes, AtomLocalAxesL)
from cryspy.C_item_loop_classes.cl_1_atom_rho_multipole import (
    AtomRhoMultipole, AtomRhoMultipoleL)
from cryspy.C_item_loop_classes.cl_1_atom_site import AtomSite
from cryspy.exceptions import CrysPyError

# distance (Å) of the dummy atoms which give the directions of the local axes
DUMMY_DISTANCE = 0.5


def take_magnetic_ion_of_atom(crystal, label: str, magnetic_ion: str = "") -> str:
    """Ion of the atom: ``magnetic_ion``, else atom_site_scat, else the
    type_symbol of atom_site when it has a charge (Cu2+)."""
    ion = magnetic_ion.strip()
    if ion == "" and crystal.is_attribute("atom_site_scat") and \
            label in list(crystal.atom_site_scat.label):
        ion = str(crystal.atom_site_scat[label].type_symbol).strip()
    if ion == "":
        type_symbol = str(crystal.atom_site[label].type_symbol)
        ion = type_symbol if ("+" in type_symbol or "-" in type_symbol) else ""
    if ion == "":
        raise CrysPyError(
            f"The ion of the atom {label} is unknown: give magnetic_ion "
            "(as Cu2+ or Fe3+).")
    return ion


def take_orbital_of_ion(ion: str) -> str:
    """Open shell of a magnetic ion from the knowledge base (Cu2+ 3d9 -> 3d)."""
    d_ion = take_magnetic_ions().get(ion)
    if d_ion is None:
        raise CrysPyError(
            f"The ion '{ion}' is not in the tables of magnetic ions: give "
            "the orbital (as 3d or 4f).")
    match = re.match(r"\s*(\d[spdf])", d_ion["configuration"])
    if match is None:
        raise CrysPyError(
            f"No open shell in the configuration '{d_ion['configuration']}' "
            f"of {ion}: give the orbital (as 3d or 4f).")
    return match.group(1)


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


def form_dummy_atoms(crystal, label: str):
    """Two dummy atoms near the atom: along c and along a.

    With ``ax1 = Z`` along c and ``ax2 = X`` along a (orthogonalised to c),
    the local axes are the Cartesian axes of CrysPy (x || a*, z || c).
    """
    atom = crystal.atom_site[label]
    unit_cell_parameters = crystal.cell.get_unit_cell_parameters()
    xyz = numpy.array([atom.fract_x, atom.fract_y, atom.fract_z], dtype=float)
    l_dummy = []
    for axis, i_axis in (("c", 2), ("a", 0)):
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
    label = atom_label.strip()
    if label not in list(crystal.atom_site.label):
        raise CrysPyError(
            f"There is no atom '{label}' in atom_site; the atoms are: "
            + ", ".join(crystal.atom_site.label) + ".")
    if crystal.is_attribute("atom_rho_multipole") and \
            label in list(crystal.atom_rho_multipole.atom_label):
        raise CrysPyError(
            f"The atom {label} has multipoles already; change its row of "
            "atom_rho_multipole or remove it first.")
    ion = take_magnetic_ion_of_atom(crystal, label, magnetic_ion)
    orbital = orbital.strip() or take_orbital_of_ion(ion)
    d_coeff = {"coeff_P00": 1.}
    d_coeff.update(parse_multipole_coefficients(coefficients))
    item_multipole = AtomRhoMultipole(
        atom_label=label, atom_type_symbol=ion, orbital=orbital,
        kappa=float(kappa), **d_coeff)
    try:
        item_multipole.get_n_zeta_coeff()
    except KeyError:
        raise CrysPyError(
            f"No radial (Slater) function of the orbital {orbital} of {ion} "
            "in the CrysPy database; give another orbital.") from None

    ls_out = [f"# Multipole form factor of {label}", "",
              f"Ion {ion}, orbital {orbital}, kappa {float(kappa):g}; "
              + ", ".join(f"{name[6:].replace('m', '-')} = {value:g}"
                          for name, value in d_coeff.items()) + "."]
    flag_axes = crystal.is_attribute("atom_local_axes") and \
        label in list(crystal.atom_local_axes.atom_label)
    if flag_axes:
        ls_out.append(f"The local axes of {label} are kept.")
    else:
        l_dummy = form_dummy_atoms(crystal, label)
        crystal.atom_site.items.extend(l_dummy)
        item_axes = AtomLocalAxes(
            atom_label=label, atom0=l_dummy[0].label, ax1="Z",
            atom1=label, atom2=l_dummy[1].label, ax2="X")
        if crystal.is_attribute("atom_local_axes"):
            crystal.atom_local_axes.items.append(item_axes)
        else:
            loop = AtomLocalAxesL()
            loop.items = [item_axes]
            crystal.add_items([loop])
        ls_out.append(
            f"Local axes: z || c, x || a* (dummy atoms {l_dummy[0].label}, "
            f"{l_dummy[1].label}); change atom_local_axes to align them with "
            "the bonds.")
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
    ls_out.append(
        "The multipoles replace the dipole form factor of the atom for magnetic "
        "space groups and propagation vectors (superspace groups) alike. "
        "Give only the coefficients allowed by the site symmetry: CrysPy "
        "does not check them; refine kappa and the coefficients against flip ratios.")
    return "\n".join(ls_out)
