"""Magnetic atoms and magnetic ions of a crystal guessed from its
composition and the knowledge base of CrysPy (cryspy/knowledge/*.md).

Used to propose values that the user confirms: an AI assistant asks "Am I
right that the magnetic ion is Tb3+?" instead of an open question.
"""
from typing import Dict, List, Optional

import numpy

from cryspy.A_functions_base.knowledge_tables import (
    split_ion, take_ion_charges, take_ions_of_element, take_magnetic_elements,
    take_magnetic_ions, take_usual_ion)
from cryspy.A_functions_base.magnetic_form_factor import get_j0_j2_parameters

# anions with a reliable formal charge: the charge balance is tried only
# when one of them is in the compound (not in metals and intermetallics)
BALANCE_ANIONS = ("O", "F", "Cl", "Br", "I", "S", "Se")
# a valence closer than this to an integer is taken as that integer
VALENCE_TOLERANCE = 0.05


def take_element(type_symbol: str) -> str:
    return split_ion(str(type_symbol))[0]


def has_form_factor(ion: str) -> bool:
    """True when the CrysPy table has the <j0> coefficients of the ion."""
    return not numpy.allclose(get_j0_j2_parameters(numpy.array([ion]))[0], 0.)


def take_unit_cell_composition(crystal) -> Dict[str, float]:
    """{element: number of atoms in the unit cell} (occupancies included)."""
    atoms = crystal.calc_atoms_in_unit_cell()
    l_site = crystal.atom_site.items
    d_count = {}
    for i_atom in atoms["atom_index"]:
        item = l_site[int(i_atom)]
        occupancy = float(item.occupancy) if item.is_attribute("occupancy") else 1.
        element = take_element(item.type_symbol)
        d_count[element] = d_count.get(element, 0.) + occupancy
    return d_count


def guess_magnetic_atoms(crystal) -> List[str]:
    """Labels of the atoms of atom_site that are magnetic: those of
    atom_site_scat when the crystal has it, otherwise the atoms of the
    elements with magnetic ions in the knowledge base (Tb, Tb_2 of TbTe3)."""
    l_label = [str(item.label) for item in crystal.atom_site.items]
    if crystal.is_attribute("atom_site_scat"):
        l_scat = [label for label in crystal.atom_site_scat.label if label in l_label]
        if l_scat:
            return l_scat
    l_magnetic = take_magnetic_elements()
    return [str(item.label) for item in crystal.atom_site.items
            if take_element(item.type_symbol) in l_magnetic]


def calc_balance_valence(crystal, element: str) -> Optional[float]:
    """Valence of the element from the charge balance of the unit cell,
    None when it does not apply: no anion of BALANCE_ANIONS, another
    element without a fixed charge (a second magnetic element, Te, ...)."""
    d_count = take_unit_cell_composition(crystal)
    d_charge = take_ion_charges()
    if element not in d_count or not any(anion in d_count for anion in BALANCE_ANIONS):
        return None
    total = 0.
    for other, number in d_count.items():
        if other == element:
            continue
        if other not in d_charge:
            return None
        total += d_charge[other] * number
    return -total / d_count[element]


def describe_ion(d_ion: dict) -> str:
    """'Tb3+ (4f8, g_J 3/2, free-ion moment 9 μB)'."""
    l_text = [d_ion["configuration"]]
    if d_ion["lande"] is not None and d_ion["moment"] > 0.:
        l_text.append(f"Lande g {d_ion['lande']:.4g}, free-ion moment {d_ion['moment']:g} μB")
    else:
        l_text.append("no moment")
    return f"{d_ion['ion']} ({', '.join(text for text in l_text if text)})"


def guess_magnetic_ion(crystal, element: str) -> dict:
    """The magnetic ion of the element in the crystal.

    The charge balance of the unit cell gives it in ionic compounds (oxides,
    fluorides, ...; Fe in Fe2O3: Fe3+); otherwise the usual ion of the
    element in the knowledge base (rare earths: R3+).

    Returns
    -------
    dict
        "ion" (as Tb3+; "" when nothing is known), "reason" (text),
        "form_factor" (True when CrysPy has its form factor), "ion_data"
        (row of the knowledge base or None).
    """
    d_ion_all = take_magnetic_ions()
    d_usual = take_usual_ion(element)
    valence = calc_balance_valence(crystal, element)
    d_out, reason = None, ""
    if valence is not None:
        n_valence = int(round(valence))
        ion = f"{element}{n_valence}+"
        if abs(valence - n_valence) > VALENCE_TOLERANCE:
            l_near = [f"{element}{int(numpy.floor(valence))}+",
                      f"{element}{int(numpy.ceil(valence))}+"]
            reason = (f"the charge balance gives {valence:.2f}+ (mixed valence of "
                      f"{' and '.join(l_near)})")
            d_out = d_ion_all.get(f"{element}{n_valence}+")
        elif ion in d_ion_all and d_ion_all[ion]["moment"] > 0.:
            d_out = d_ion_all[ion]
            reason = "from the charge balance of the unit cell"
        else:
            reason = f"the charge balance gives {ion}, which has no moment"
    if d_out is None and d_usual is not None:
        d_out = d_usual
        text_usual = f"the usual ion of {element} in compounds"
        reason = f"{text_usual} ({reason})" if reason else text_usual
    if d_out is None:
        l_known = [d_ion["ion"] for d_ion in take_ions_of_element(element)]
        return {"ion": "", "reason": f"{element} is not in the tables of magnetic ions"
                + (f" (known: {', '.join(l_known)})" if l_known else ""),
                "form_factor": False, "ion_data": None}
    return {"ion": d_out["ion"], "reason": reason, "form_factor": has_form_factor(d_out["ion"]),
            "ion_data": d_out}
