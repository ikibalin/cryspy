"""Knowledge base of CrysPy (cryspy/knowledge/*.md): the tables read by
CrysPy (magnetic ions, charges), the guess of the magnetic ion of a crystal
and the search of the knowledge base."""
import pytest

import cryspy
from cryspy.A_functions_base.knowledge_tables import (
    take_ion_charges, take_lande_factor, take_magnetic_elements, take_magnetic_ions,
    take_usual_ion)
from cryspy.D_functions_item_loop.function_1_magnetic_ions import (
    guess_magnetic_atoms, guess_magnetic_ion, has_form_factor)
from cryspy.procedure_assistant.knowledge import split_knowledge_chunks

LANTHANIDES = ("Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb")


def calc_hund_4f(n_f: int):
    """(g_J, g_J·J) of the Hund's rules ground multiplet of n_f 4f
    electrons, (None, 0) for J = 0: an independent check of the table."""
    m_l = [3, 2, 1, 0, -1, -2, -3]
    n_up = min(n_f, 7)
    spin = 0.5 * (n_up - (n_f - n_up))
    orbital = abs(sum(m_l[:n_up]) + sum(m_l[:n_f - n_up]))
    total = abs(orbital - spin) if n_f < 7 else orbital + spin
    if total == 0.:
        return None, 0.
    j_sq = total * (total + 1.)
    lande = 1. + (j_sq + spin * (spin + 1.) - orbital * (orbital + 1.)) / (2. * j_sq)
    return lande, lande * total


def test_4f_table_follows_hunds_rules():
    for d_ion in take_magnetic_ions().values():
        if d_ion["element"] not in LANTHANIDES:
            continue
        n_f = LANTHANIDES.index(d_ion["element"]) + 4 - d_ion["charge"]
        assert d_ion["configuration"] == f"4f{n_f}", d_ion["ion"]
        lande, moment = calc_hund_4f(n_f)
        if lande is None:
            assert d_ion["lande"] is None and d_ion["moment"] == 0., d_ion["ion"]
            continue
        assert d_ion["lande"] == pytest.approx(lande, abs=1e-6), d_ion["ion"]
        assert d_ion["moment"] == pytest.approx(moment, abs=0.006), d_ion["ion"]


def test_3d_moments_are_spin_only():
    for d_ion in take_magnetic_ions().values():
        if d_ion["file"] != "magnetic_ions_3d.md" or d_ion["moment"] == 0.:
            continue
        assert d_ion["lande"] == 2., d_ion["ion"]


def test_every_element_has_one_usual_ion():
    for element in take_magnetic_elements():
        l_usual = [d_ion for d_ion in take_magnetic_ions().values()
                   if d_ion["element"] == element and d_ion["usual"]]
        assert len(l_usual) == 1, element


def test_missing_form_factors_are_noted():
    """An ion without a CrysPy form factor says so in its notes (the
    usual ions of the 4d, 5d metals and Ce3+), so the knowledge base does
    not promise what CrysPy cannot calculate."""
    for d_ion in take_magnetic_ions().values():
        if d_ion["usual"] and d_ion["moment"] > 0. and not has_form_factor(d_ion["ion"]):
            assert "No" in d_ion["notes"] and "form factor" in d_ion["notes"], d_ion["ion"]


def test_lande_factor_from_the_table():
    assert take_lande_factor("Tb3+") == pytest.approx(1.5)
    assert take_lande_factor("Ce3+") == pytest.approx(6. / 7.)
    assert take_lande_factor("Fe3+") == 2.
    assert take_lande_factor("Eu3+") == 2.  # J = 0: no g_J, spin only as before
    assert take_lande_factor("Xx5+") == 2.
    assert take_usual_ion("Tb")["ion"] == "Tb3+"
    assert take_ion_charges()["O"] == -2


def form_p1_crystal(atoms: str) -> object:
    return cryspy.str_to_globaln(f"""global_
data_test
_space_group_name_H-M_alt 'P 1'
_cell_length_a 5.0
_cell_length_b 5.0
_cell_length_c 5.0
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
{atoms}
""").crystal_test


@pytest.mark.parametrize("atoms, element, ion, word", [
    ("Fe1 Fe 0 0 0 1\nFe2 Fe 0.5 0 0 1\nO1 O 0 0.3 0 1\nO2 O 0 0.6 0 1\nO3 O 0 0 0.5 1",
     "Fe", "Fe3+", "charge balance"),
    ("Mn1 Mn 0 0 0 1\nO1 O 0.5 0.5 0.5 1", "Mn", "Mn2+", "charge balance"),
    ("La1 La 0 0 0 1\nMn1 Mn 0.5 0.5 0.5 1\nO1 O 0.5 0 0 1\nO2 O 0 0.5 0 1\nO3 O 0 0 0.5 1",
     "Mn", "Mn3+", "charge balance"),
    ("Tb1 Tb 0 0 0 1\nTe1 Te 0.5 0 0 1\nTe2 Te 0 0.5 0 1\nTe3 Te 0 0 0.5 1",
     "Tb", "Tb3+", "usual ion"),
    ("Fe1 Fe 0 0 0 1\nFe2 Fe 0.5 0 0 1\nFe3 Fe 0 0.5 0 1\nO1 O 0 0 0.5 1\nO2 O 0.5 0.5 0 1"
     "\nO3 O 0.5 0 0.5 1\nO4 O 0 0.5 0.5 1", "Fe", "Fe3+", "mixed valence"),
])
def test_guess_magnetic_ion(atoms, element, ion, word):
    crystal = form_p1_crystal(atoms)
    d_guess = guess_magnetic_ion(crystal, element)
    assert d_guess["ion"] == ion
    assert word in d_guess["reason"]
    assert d_guess["form_factor"]


def test_guess_magnetic_atoms():
    crystal = form_p1_crystal("Tb1 Tb 0 0 0 1\nTe1 Te 0.5 0 0 1\nO1 O 0 0.5 0 1")
    assert guess_magnetic_atoms(crystal) == ["Tb1"]


def test_table_rows_are_fragments():
    content = "---\ntitle: Ions\n---\nText before.\n\n| Ion | Usual |\n|---|---|\n| Tb3+ | yes |\n"
    l_chunk = split_knowledge_chunks(content)
    assert (8, "Ions: Ion Tb3+ (Tb); Usual yes") in l_chunk
    assert any("Text before." in text for _, text in l_chunk)


@pytest.mark.parametrize("query, text", [
    ("usual ion of Tb", "Ion Tb3+ (Tb); Usual yes"),
    ("Lande factor Ho3+", "Ion Ho3+ (Ho)"),
    ("moment of Co2+", "Ion Co2+ (Co)"),
])
def test_search_knowledge(query, text):
    l_found = cryspy.search_knowledge(query, n_results=2)
    assert text in l_found[0]["text"]
    result = cryspy.call_assistant_tool(None, "search_knowledge", {"query": query})
    assert text in result["text"] and "knowledge/" in result["text"]


def test_lande_written_with_its_source():
    """lande 0: the value of the knowledge base is written to atom_site_scat
    and the report says where it comes from."""
    from cryspy.D_functions_item_loop.function_1_magnetic_k import add_magnetic_scattering
    crystal = form_p1_crystal("U1 U 0 0 0 1\nFe1 Fe 0.5 0.5 0.5 1")
    l_text = add_magnetic_scattering(crystal, ["U1"], "U4+", 0.)
    assert crystal.atom_site_scat["U1"].lande == pytest.approx(0.8)
    assert "Lande 0.8: g_J of Hund's rules for U4+" in l_text[0]
    l_text = add_magnetic_scattering(crystal, ["Fe1"], "Fe3+", 0.)
    assert "Lande 2: spin-only value of Fe3+" in l_text[0]
    crystal_given = form_p1_crystal("Fe1 Fe 0 0 0 1")
    assert "Lande 2.1: given" in add_magnetic_scattering(crystal_given, ["Fe1"], "Fe3+", 2.1)[0]
