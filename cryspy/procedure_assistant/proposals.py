"""Values proposed to the user before a question is asked.

When a method or a procedure needs the propagation vector, the magnetic
atoms or the magnetic ion, the value is estimated from the object first:

- k from the fractional indices of the measured satellites
  (take_propagation_vectors), the vectors not yet in the crystal first;
  the three components k_1, k_2, k_3 are asked as one vector "k";
- the magnetic atoms: those with atom_site_scat, or the atoms of the
  elements with magnetic ions in the knowledge base (cryspy/knowledge);
- the magnetic ion (when the magnetic atoms have no atom_site_scat and no
  charge in type_symbol): the charge balance of the unit cell or the usual
  ion of the element (Tb -> Tb3+).

The question then carries "proposal" and "s_proposal" (see
call_assistant_tool): the user answers "yes" or gives another value, and is
asked an open question only after "no".
"""
from typing import List, Optional

import numpy

from cryspy.B_parent_classes.cl_4_global import GlobalN

K_NAMES = ("k_1", "k_2", "k_3")
D_K_QUESTION = {
    "name": "k", "question": "the propagation vector k (reciprocal lattice units)",
    "type": "array", "items": {"type": "number"}, "minItems": 3, "maxItems": 3,
    "hint": "three numbers or fractions, as 0 0 1/2"}


def expand_k_argument(arguments: dict, l_param: list) -> dict:
    """arguments with "k" (a vector: the answer of the question of
    D_K_QUESTION, or of a model) written as k_1, k_2, k_3 when the function
    has these parameters and no parameter k."""
    if "k" not in arguments or "k" in l_param or not all(name in l_param for name in K_NAMES):
        return arguments
    from cryspy.procedure_assistant.tools import convert_argument
    arguments = dict(arguments)
    l_value = convert_argument(arguments.pop("k"), D_K_QUESTION)
    for name, value in zip(K_NAMES, l_value):
        arguments.setdefault(name, value)
    return arguments


def take_crystal_for(cryspy_object, block=None, crystal_name: str = ""):
    """The crystal the call is about: the block itself, the named one or the
    only crystal of the object; None when it is not clear."""
    from cryspy.E_data_classes.cl_1_crystal import Crystal
    if isinstance(block, Crystal):
        return block
    if not isinstance(cryspy_object, GlobalN):
        return None
    l_crystal = [item for item in cryspy_object.items if isinstance(item, Crystal)]
    if crystal_name:
        l_crystal = [item for item in l_crystal
                     if item.data_name.lower() == str(crystal_name).strip().lower()]
    return l_crystal[0] if len(l_crystal) == 1 else None


def take_crystal_k_vectors(crystal) -> list:
    """Canonical propagation vectors already in the crystal."""
    from cryspy.procedure_rhochi.rhochi_magnetic_irreps import take_canonical_k
    if crystal is None or not crystal.is_attribute("cell_wave_vector"):
        return []
    return [take_canonical_k([item.x, item.y, item.z]) for item in crystal.cell_wave_vector.items]


def propose_k(cryspy_object, crystal) -> dict:
    """The question of the propagation vector with the proposal from the
    measured satellites (the vector with most satellites not yet in the
    crystal)."""
    from cryspy.exceptions import CrysPyError
    from cryspy.procedure_rhochi.rhochi_magnetic_irreps import (
        form_k_fraction_text, take_canonical_k, take_propagation_vectors)
    d_question = dict(D_K_QUESTION)
    try:
        l_vector = take_propagation_vectors(cryspy_object)["vectors"] \
            if isinstance(cryspy_object, GlobalN) else []
    except CrysPyError:  # no single-crystal experiments
        l_vector = []
    if len(l_vector) == 0:
        return d_question
    l_present = take_crystal_k_vectors(crystal)
    l_new = [d_vector for d_vector in l_vector
             if not any(numpy.allclose(take_canonical_k(d_vector["k"]), present, atol=1e-4)
                        for present in l_present)]
    l_text = [f"{form_k_fraction_text(d_vector['k'])} with {d_vector['n_satellite']} satellites"
              for d_vector in l_vector]
    d_question["hint"] = (f"{d_question['hint']}; the measured satellites give "
                          + "; ".join(l_text))
    if len(l_new) == 0:
        return d_question
    d_best = l_new[0]
    s_proposal = (f"{form_k_fraction_text(d_best['k'])} (from {d_best['n_satellite']} "
                  f"measured satellites in {', '.join(d_best['experiments'])})")
    l_other = [form_k_fraction_text(d_vector["k"]) for d_vector in l_new[1:]]
    if l_other:
        s_proposal += f"; the data have also {', '.join(l_other)}, the next k"
    d_question["proposal"] = [float(value) for value in d_best["k"]]
    d_question["s_proposal"] = s_proposal
    return d_question


def propose_magnetic_atoms(crystal, d_missing: dict) -> dict:
    from cryspy.D_functions_item_loop.function_1_magnetic_ions import (
        guess_magnetic_atoms, guess_magnetic_ion, take_element)
    l_label = guess_magnetic_atoms(crystal)
    if len(l_label) == 0:
        return d_missing
    l_element = []
    for item in crystal.atom_site.items:
        if str(item.label) in l_label and take_element(item.type_symbol) not in l_element:
            l_element.append(take_element(item.type_symbol))
    l_ion = [guess_magnetic_ion(crystal, element)["ion"] or element for element in l_element]
    proposal = ", ".join(l_label)
    return dict(d_missing, proposal=proposal,
                s_proposal=f"{proposal} ({', '.join(l_ion)}: the magnetic element"
                + ("s" if len(l_ion) > 1 else "") + " of the crystal)")


def take_atoms_without_ion(crystal, magnetic_atoms: str) -> list:
    """Magnetic atoms without atom_site_scat whose type_symbol has no
    charge: add_magnetic_scattering needs magnetic_ion for them."""
    from cryspy.D_functions_item_loop.function_1_magnetic_k import split_labels
    l_scat = list(crystal.atom_site_scat.label) if crystal.is_attribute("atom_site_scat") else []
    d_site = {str(item.label): item for item in crystal.atom_site.items}
    return [label for label in split_labels(str(magnetic_atoms))
            if label in d_site and label not in l_scat
            and not any(sign in str(d_site[label].type_symbol) for sign in "+-")]


def propose_magnetic_ion(crystal, l_label: list) -> Optional[dict]:
    """Question of the magnetic ion of the atoms l_label (one element), with
    the proposal of guess_magnetic_ion; None for atoms of several elements
    (the error of the method names them)."""
    from cryspy.D_functions_item_loop.function_1_magnetic_ions import (
        describe_ion, guess_magnetic_ion, take_element)
    l_element = sorted({take_element(crystal.atom_site[label].type_symbol) for label in l_label})
    if len(l_element) != 1:
        return None
    element = l_element[0]
    d_guess = guess_magnetic_ion(crystal, element)
    d_question = {"name": "magnetic_ion", "type": "string",
                  "question": f"the magnetic ion of {', '.join(l_label)} (its form factor)",
                  "hint": "an ion as Tb3+, Fe3+, Mn2+"}
    if d_guess["ion"] == "":
        d_question["hint"] += f"; {d_guess['reason']}"
    elif not d_guess["form_factor"]:
        d_question["hint"] += (f"; {d_guess['ion']} is {d_guess['reason']}, but CrysPy has no "
                               "form factor of it: give another ion or atom_type_scat")
    else:
        d_question["proposal"] = d_guess["ion"]
        d_question["s_proposal"] = f"{describe_ion(d_guess['ion_data'])}, {d_guess['reason']}"
    return d_question


def propose_arguments(cryspy_object, block, l_param: list, kwargs: dict,
                      l_missing: List[dict]) -> List[dict]:
    """The questions l_missing of a call with the values proposed from the
    object; also the question of the magnetic ion when it has a default but
    the call would fail without it.

    l_param: names of the parameters of the function; kwargs: the values
    known already."""
    crystal = take_crystal_for(cryspy_object, block, kwargs.get("crystal_name", ""))
    l_name = [d_missing["name"] for d_missing in l_missing]
    l_out = []
    if all(name in l_name for name in K_NAMES):
        for d_missing in l_missing:
            if d_missing["name"] == K_NAMES[0]:
                l_out.append(propose_k(cryspy_object, crystal))
            elif d_missing["name"] not in K_NAMES:
                l_out.append(d_missing)
    else:
        l_out = list(l_missing)
    if crystal is None:
        return l_out
    for ind, d_missing in enumerate(l_out):
        if d_missing["name"] == "magnetic_atoms" and "proposal" not in d_missing:
            l_out[ind] = propose_magnetic_atoms(crystal, d_missing)
    # the ion is asked when the magnetic atoms are known (the next round after their question)
    if "magnetic_ion" in l_param and not str(kwargs.get("magnetic_ion", "")).strip() \
            and "magnetic_atoms" in kwargs and "magnetic_ion" not in l_name:
        l_label = take_atoms_without_ion(crystal, kwargs["magnetic_atoms"])
        if l_label:
            d_question = propose_magnetic_ion(crystal, l_label)
            if d_question is not None:
                l_out.append(d_question)
    return l_out
