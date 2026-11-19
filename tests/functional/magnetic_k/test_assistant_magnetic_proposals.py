"""The AI assistant proposes k, the magnetic atoms and the magnetic ion
from the data (TbTe3-like crystal with the satellites of two propagation
vectors), the user confirms them; a simulation adds the satellites of a k
that the crystal does not have (zero intensities, a list for planning)."""
import warnings

import numpy
import pytest

import cryspy
from test_magnetic_k_path_from_paramagnetic import CRYSTAL, INTEGER, NUCLEAR, form_diffrn

K_00H = numpy.array([0., 0., 0.5])
K_HH0 = numpy.array([0.5, 0.5, 0.])


def form_object():
    """Paramagnetic crystal, nuclear reflections, 28 satellites of (1/2, 1/2, 0)
    and 14 of (0, 0, 1/2) (+ k only)."""
    sat_hh0 = numpy.concatenate([INTEGER + K_HH0[:, None], INTEGER - K_HH0[:, None]], axis=1)
    sat_00h = INTEGER + K_00H[:, None]
    return cryspy.str_to_globaln(CRYSTAL + form_diffrn("nuc", NUCLEAR) + form_diffrn(
        "sat_hh0", sat_hh0) + form_diffrn("sat_00h", sat_00h))


def answer_yes(obj, arguments: dict) -> dict:
    """Call add_propagation_vector, answering 'yes' to every proposal (as
    the Dialog of CrysPy Editor does), until it runs."""
    l_question = []
    for _ in range(4):
        result = cryspy.call_assistant_tool(obj, "call_method", arguments)
        if not result.get("missing"):
            return {"result": result, "questions": l_question}
        for d_missing in result["missing"]:
            assert "proposal" in d_missing, d_missing
            l_question.append(d_missing)
            arguments[d_missing["name"]] = d_missing["proposal"]
    raise AssertionError("questions do not end")


def test_k_atoms_and_ion_are_proposed():
    obj = form_object()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        d_out = answer_yes(obj, {"method": "add_propagation_vector"})
    l_question = d_out["questions"]
    assert [d_missing["name"] for d_missing in l_question] == ["k", "magnetic_atoms",
                                                               "magnetic_ion"]
    d_k, d_atoms, d_ion = l_question
    assert d_k["proposal"] == pytest.approx(list(K_HH0))  # more satellites first
    assert "(1/2, 1/2, 0)" in d_k["s_proposal"] and "(0, 0, 1/2)" in d_k["s_proposal"]
    assert d_atoms["proposal"] == "Tb"
    assert d_ion["proposal"] == "Tb3+" and "usual ion" in d_ion["s_proposal"]
    crystal = obj.crystal_tbte3
    assert crystal.atom_site_scat["Tb"].lande == pytest.approx(1.5)
    # the second call proposes the other vector
    result = cryspy.call_assistant_tool(obj, "call_method", {"method": "add_propagation_vector"})
    d_k = result["missing"][0]
    assert d_k["proposal"] == pytest.approx(list(K_00H))
    # the magnetic atoms are those of atom_site_scat now: no question of the ion
    assert [d_missing["name"] for d_missing in result["missing"]] == ["k", "magnetic_atoms"]


def test_refused_proposal_takes_the_answer():
    """The user answers another k: it is written as k_1, k_2, k_3."""
    obj = form_object()
    result = cryspy.call_assistant_tool(obj, "call_method", {
        "method": "report_irreps_of_k", "k": cryspy.convert_assistant_answer(
            {"type": "array", "items": {"type": "number"}, "minItems": 3, "maxItems": 3},
            "0 0 1/2")})
    assert not result.get("missing")
    assert "0.5" in result["text"]


def test_rank_procedure_gets_the_proposals():
    obj = form_object()
    result = cryspy.call_assistant_tool(obj, "run_procedure", {
        "name": "rhochi_rank_magnetic_representations"})
    l_name = [d_missing["name"] for d_missing in result["missing"]]
    assert l_name == ["k", "magnetic_atoms"]
    assert all("proposal" in d_missing for d_missing in result["missing"])


def test_simulation_with_a_k_of_a_paramagnetic_crystal():
    obj = cryspy.str_to_globaln(CRYSTAL)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = cryspy.call_assistant_tool(obj, "run_procedure", {
            "name": "simulation_single_crystal", "wavelength": 2.36, "ttheta_max": 60.,
            "satellites_only": True, "k": [0, 0, 0.5]})
    assert "zero intensity" in result["text"]
    refln = obj.diffrn_sc.diffrn_refln
    index_l = numpy.array(refln.index_l, dtype=float)
    assert len(refln.items) > 5
    assert numpy.allclose(numpy.abs(numpy.mod(index_l, 1.)), 0.5)
    assert numpy.allclose(refln.intensity_calc, 0.)


def test_simulation_adds_a_second_k():
    obj = cryspy.str_to_globaln(CRYSTAL)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        obj.crystal_tbte3.add_propagation_vector(*K_00H, "Tb", "1", 5.0, magnetic_ion="Tb3+")
        text = cryspy.simulation_single_crystal(obj, 2.36, 1., 60., True, *K_HH0)
    refln = obj.diffrn_sc.diffrn_refln
    index = numpy.array([refln.index_h, refln.index_k, refln.index_l], dtype=float)
    intensity = numpy.array(refln.intensity_calc)
    flag_00h = numpy.abs(numpy.mod(index[2], 1.) - 0.5) < 1e-4
    assert intensity[flag_00h].max() > 0.  # satellites of the crystal's k
    assert numpy.allclose(intensity[~flag_00h], 0.)  # satellites of the new k
    assert (~flag_00h).sum() > 0 and "not a propagation vector" in text
    # a k the crystal has is not added twice
    obj_same = cryspy.str_to_globaln(CRYSTAL)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        obj_same.crystal_tbte3.add_propagation_vector(*K_00H, "Tb", "1", 5.0,
                                                      magnetic_ion="Tb3+")
        cryspy.simulation_single_crystal(obj_same, 2.36, 1., 60., True)
        n_own = len(obj_same.diffrn_sc.diffrn_refln.items)
        cryspy.simulation_single_crystal(obj_same, 2.36, 1., 60., True, 0., 0., -0.5)
    assert len(obj_same.items[-1].diffrn_refln.items) == n_own


def test_points_of_the_brillouin_zone():
    """k of TbTe3 are the points S and Z of the C-centred orthorhombic zone."""
    from cryspy.procedure_rhochi.rhochi_magnetic_irreps import name_special_k_point
    crystal = cryspy.str_to_globaln(CRYSTAL).crystal_tbte3
    assert name_special_k_point(crystal, (0.5, -0.5, 0.)) == "S"
    assert name_special_k_point(crystal, (0., 0., 0.5)) == "Z"
    assert name_special_k_point(crystal, (1., 0., 0.)) == "Y"  # not Γ in a C cell
    assert name_special_k_point(crystal, (0.5, 0., 0.)) == ""
    result = cryspy.call_assistant_tool(form_object(), "call_method",
                                        {"method": "add_propagation_vector"})
    assert "(1/2, 1/2, 0) (point S)" in result["missing"][0]["s_proposal"]
