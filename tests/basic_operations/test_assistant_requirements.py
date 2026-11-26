"""Categories a function needs and adds (sections Needs and Adds of its
docstring): the assistant adds only what is missing and asks the possible
values (procedure_assistant/requirements.py)."""
import cryspy
from cryspy.procedure_assistant.methods import form_matching_methods_text
from cryspy.procedure_assistant.procedures import parse_docstring
from cryspy.procedure_assistant.requirements import check_requirements, take_requirements

# MnF2 as in the file of the user: msCIF superspace group with k = 0, the
# ion of atom_site_scat of F is O2- (its dipole form factor), two dummy
# atoms D1 (along [110] from F) and D2 (along c from Mn) put there before
MNF2 = """global_
data_mnf2
_cell_length_a 4.8736
_cell_length_b 4.8736
_cell_length_c 3.3
_cell_angle_alpha 90.0
_cell_angle_beta 90.0
_cell_angle_gamma 90.0
_parent_space_group_name_H-M_alt "P 42/m n m"
_parent_space_group_IT_number 136
_parent_space_group_IT_coordinate_system_code 1
loop_
_cell_wave_vector_seq_id
_cell_wave_vector_x
_cell_wave_vector_y
_cell_wave_vector_z
1 0.0 0.0 0.0
loop_
_space_group_symop_magn_ssg_operation_algebraic
_space_group_symop_magn_ssg_operation_id
x1,x2,x3,x4,+1 1
-x1,-x2,x3,x4,+1 2
-x2,-x1,-x3,x4,-1 3
x2,x1,-x3,x4,-1 4
-x1+1/2,x2+1/2,x3+1/2,x4,+1 5
x1+1/2,-x2+1/2,x3+1/2,x4,+1 6
-x2+1/2,x1+1/2,-x3+1/2,x4,-1 7
x2+1/2,-x1+1/2,-x3+1/2,x4,-1 8
loop_
_space_group_symop_magn_ssg_centering_algebraic
_space_group_symop_magn_ssg_centering_id
x1,x2,x3,x4,+1 1
x1,x2,x3,x4+1/2,-1 2
loop_
_atom_site_fourier_wave_vector_seq_id
_atom_site_fourier_wave_vector_q1_coeff
1 1
loop_
_atom_site_moment_fourier_atom_site_label
_atom_site_moment_fourier_axis
_atom_site_moment_fourier_wave_vector_seq_id
_atom_site_moment_fourier_param.cos
_atom_site_moment_fourier_param.sin
_atom_site_moment_fourier_param.cos_symmform
_atom_site_moment_fourier_param.sin_symmform
Mn x 1 0.0 0.0 mxc1 mxs1
Mn y 1 0.0 0.0 -mxc1 -mxs1
Mn z 1 -4.368 0.0 mzc1 mzs1
F x 1 -0.060 0.0 mxc1 mxs1
F y 1 0.060 0.0 -mxc1 -mxs1
F z 1 0.032 0.0 mzc1 mzs1
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
_atom_site_adp_type
_atom_site_B_iso_or_equiv
Mn Mn2+ 0.0 0.0 0.0 1.0 Biso 0.22
F F1- 0.303691 0.303691 0.0 1.0 Biso 0.46
D1 Dummy 0.4 0.4 0.0 1.0 Biso 0.0
D2 Dummy 0.0 0.0 0.1 1.0 Biso 0.0
loop_
_atom_site_scat_label
_atom_site_scat_type_symbol
_atom_site_scat_Lande
_atom_site_scat_kappa
Mn Mn3+ 2.0 1.0
F O2- 2.0 1.0
"""

REQUEST = "Add to mnf2 crystal description of magnetic form-factors for Mn and F through multipoles"


def call_method(obj, **arguments):
    return cryspy.call_assistant_tool(
        obj, "call_method", dict(method="add_multipole_form_factor", **arguments))


def test_parse_needs_and_adds():
    d_doc = parse_docstring(cryspy.Crystal.add_multipole_form_factor.__doc__)
    assert [item["category"] for item in d_doc["adds"]] == ["atom_rho_multipole", "atom_local_axes"]
    # the needed categories come from the class of the added rows
    d_requirement = take_requirements(cryspy.Crystal.add_multipole_form_factor)
    assert [(item["category"], item["column"], item["parameter"])
            for item in d_requirement["needs"]] == [
        ("atom_site", "label", "atom_label"), ("atom_site_scat", "label", "atom_label")]
    # the sections are not in the description (method search, menus)
    assert "Needs" not in d_doc["description"]


def test_check_requirements_of_multipoles():
    obj = cryspy.str_to_globaln(MNF2)
    crystal = obj.crystal_mnf2
    d_check = check_requirements(obj, type(crystal).add_multipole_form_factor, crystal)
    d_choice = d_check["choices"]["atom_label"]
    assert d_choice["values"] == ["Mn", "F", "D1", "D2"]
    assert d_choice["eligible"] == ["Mn", "F"]
    assert d_choice["lacking"] == {"D1": ["atom_site_scat"], "D2": ["atom_site_scat"]}
    assert d_choice["several"]
    assert all(len(d_item["blocks"]) == 0 for d_item in d_check["adds"])
    d_check = check_requirements(obj, type(crystal).add_multipole_form_factor, crystal,
                                 {"atom_label": "Mn and D1"})
    assert len(d_check["errors"]) == 1 and "D1 has no row of atom_site_scat" in d_check["errors"][0]


def test_request_of_the_user_gives_the_multipole_method_with_its_requirements():
    """The request gave add_magnetic_model and the assistant asked for the
    propagation vector k (add_multipole_form_factor was not in develop):
    the method comes first with what the crystal has and misses."""
    obj = cryspy.str_to_globaln(MNF2)
    text = form_matching_methods_text(obj, REQUEST)
    assert text.split("\n")[1].startswith("- add_multipole_form_factor of crystal_mnf2:")
    assert "add_magnetic_model" not in text and "add_propagation_vector" not in text
    assert "Requirements of add_multipole_form_factor:" in text
    assert "- atom_label: Mn, F (one or several); D1, D2 lack atom_site_scat" in text
    assert "adds atom_rho_multipole (atom_label = atom_label): not in the object" in text


def test_missing_atoms_are_asked_with_the_possible_ones():
    obj = cryspy.str_to_globaln(MNF2)
    d_out = call_method(obj)
    (d_missing, ) = d_out["missing"]
    assert d_missing["name"] == "atom_label"
    assert d_missing["options"] == ["Mn, F", "Mn", "F"]
    assert "D1, D2 lack atom_site_scat" in d_missing["hint"]
    # a dummy atom cannot get multipoles: the needed category is named
    d_out = call_method(obj, atom_label="D1")
    assert not d_out.get("missing")
    assert 'D1 has no row of atom_site_scat' in d_out["text"]
    assert '(add_item {"category": "atom_site_scat"})' in d_out["text"]
    assert not obj.crystal_mnf2.is_attribute("atom_rho_multipole")


def test_multipoles_of_mn_and_f_step_by_step():
    """Mn: ion and orbital are sure (Mn3+ of atom_site_scat, 3d). F: the
    ion of atom_site_scat is O2- (another element) and O2- is not a
    magnetic ion: the ion and then the orbital are asked with the possible
    values; nothing is changed before the answers."""
    obj = cryspy.str_to_globaln(MNF2)
    crystal = obj.crystal_mnf2
    d_out = call_method(obj, atom_label="Mn and F")
    (d_ion, ) = d_out["missing"]
    assert d_ion["name"] == "magnetic_ion" and d_ion["proposal"] == "F: O2-"
    assert d_ion["options"] == ["F: O2-", "F: F1-"]
    assert not crystal.is_attribute("atom_rho_multipole")
    d_out = call_method(obj, atom_label="Mn and F", magnetic_ion=d_ion["proposal"])
    (d_orbital, ) = d_out["missing"]
    assert d_orbital["name"] == "orbital" and d_orbital["proposal"] == "F: 2p"
    assert d_orbital["options"] == ["F: 1s", "F: 2s", "F: 2p"]
    assert not crystal.is_attribute("atom_rho_multipole")
    d_out = call_method(obj, atom_label="Mn and F", magnetic_ion="F: O2-",
                        orbital=d_orbital["proposal"])
    assert not d_out.get("missing"), d_out["text"]
    d_row = {item.atom_label: (item.atom_type_symbol, item.orbital)
             for item in crystal.atom_rho_multipole.items}
    assert d_row == {"Mn": ("Mn3+", "3d"), "F": ("O2-", "2p")}
    # the dummy atoms put there before give an axis each
    d_axes = {item.atom_label: (item.atom0, item.ax1, item.atom2, item.ax2)
              for item in crystal.atom_local_axes.items}
    assert d_axes == {"Mn": ("D2", "Z", "DMn_a", "X"), "F": ("D1", "X", "DF_c", "Z")}
    # nothing is left to add: the reason is given, no question
    d_out = call_method(obj)
    assert not d_out.get("missing")
    assert "Mn, F have atom_rho_multipole, atom_local_axes already" in d_out["text"]


def test_atoms_with_multipoles_are_kept():
    obj = cryspy.str_to_globaln(MNF2)
    crystal = obj.crystal_mnf2
    crystal.add_multipole_form_factor("Mn", coefficients="P20 0.1")
    text = crystal.add_multipole_form_factor("Mn, F", magnetic_ion="F: F1-", orbital="F: 2p")
    assert "Mn: multipoles already there, not changed." in text
    assert [item.atom_label for item in crystal.atom_rho_multipole.items] == ["Mn", "F"]
    assert crystal.atom_rho_multipole["Mn"].coeff_P20 == 0.1


def test_tool_check_requirements_of_a_procedure():
    obj = cryspy.str_to_globaln("""global_
data_pbso4
_cell_length_a 8.4779
_cell_length_b 5.3968
_cell_length_c 6.9581
_cell_angle_alpha 90.0
_cell_angle_beta 90.0
_cell_angle_gamma 90.0
_space_group_name_H-M_alt "P n m a"
_space_group_IT_coordinate_system_code abc
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Pb Pb 0.1872 0.25 0.1662

data_pd
_setup_wavelength 1.912
_setup_field 0.0
_setup_offset_2theta 0.0
_setup_radiation neutrons
_pd_instr_resolution_U 0.12
_pd_instr_resolution_V -0.33
_pd_instr_resolution_W 0.28
_pd_instr_resolution_X 0.0
_pd_instr_resolution_Y 0.15
_range_2theta_min 15.0
_range_2theta_max 95.0
loop_
_pd_background_2theta
_pd_background_intensity
10.0 50.0
100.0 50.0
loop_
_phase_label
_phase_scale
pbso4 1.0
loop_
_pd_meas_2theta
_pd_meas_intensity
_pd_meas_intensity_sigma
20.0 100.0 10.0
30.0 100.0 10.0
""")
    text = cryspy.call_assistant_tool(
        obj, "check_requirements", {"name": "add_texture_spherical_harmonics"})["text"]
    assert "- needs phase (label = phase_label): rows pbso4" in text
    assert "- adds texture_sh_coefficient (phase_label = phase_label): not in the object" in text
    assert "- phase_label: pbso4" in text
    d_out = cryspy.call_assistant_tool(obj, "run_procedure", {
        "name": "add_texture_spherical_harmonics", "phase_label": "pbso5"})
    assert "no pbso5 in phase.label (phase_label); the values are pbso4" in d_out["text"]
    text = cryspy.call_assistant_tool(
        obj, "check_requirements", {"name": "texture_sh_coefficient"})["text"]
    assert "needs a row of texture_sh for the same phase_label; present: none" in text
    assert "new rows possible for: pbso4" in text
    assert "procedure add_texture_spherical_harmonics" in text
    text = cryspy.call_assistant_tool(
        obj, "check_requirements", {"name": "add_restraint"})["text"]
    assert "does not name the categories" in text


# --- any category by add_item, from the definitions of the classes ---------

def add_item(obj, **arguments):
    return cryspy.call_assistant_tool(obj, "add_item", dict(category="atom_rho_multipole", **arguments))


def test_class_dependencies_of_multipoles():
    from cryspy.C_item_loop_classes.cl_1_atom_rho_multipole import AtomRhoMultipole
    from cryspy.procedure_assistant.dependencies import take_key_attribute
    assert take_key_attribute(AtomRhoMultipole) == "atom_label"
    assert AtomRhoMultipole.REQUIRED_CATEGORIES == ("atom_site_scat", "atom_local_axes")
    obj = cryspy.str_to_globaln(MNF2)
    d_proposal = AtomRhoMultipole.propose_values(
        obj.crystal_mnf2, {"atom_label": "F", "atom_type_symbol": "O2-"})
    assert d_proposal["atom_type_symbol"]["options"] == ["O2-", "F1-"]
    assert d_proposal["orbital"]["proposal"] == "2p" and not d_proposal["orbital"]["sure"]
    d_proposal = AtomRhoMultipole.propose_values(
        obj.crystal_mnf2, {"atom_label": "Mn", "atom_type_symbol": "Mn3+"})
    assert list(d_proposal) == ["orbital"] and d_proposal["orbital"]["sure"]


def test_add_item_asks_the_atoms_with_the_possible_ones():
    obj = cryspy.str_to_globaln(MNF2)
    (d_question, ) = add_item(obj)["missing"]
    assert d_question["name"] == "atom_label"
    assert d_question["options"] == ["Mn, F", "Mn", "F"]
    assert "D1, D2: no atom_site_scat" in d_question["hint"]


def test_add_item_adds_multipoles_with_the_rows_they_need():
    """No special method: atom_rho_multipole of Mn and F is added with
    atom_local_axes (required by its class); atom_site_scat is there. The
    sure values are taken (Mn3+ of atom_site_scat, 3d), the others are
    asked in one round with proposals: the ion and orbital of F, the
    local axes from the nearest atoms (D2 along c for Mn, D1 for F)."""
    obj = cryspy.str_to_globaln(MNF2)
    crystal = obj.crystal_mnf2
    d_out = add_item(obj, atom_label="Mn and F")
    d_question = {item["name"]: item for item in d_out["missing"]}
    assert set(d_question) == {
        "atom_type_symbol", "orbital", "atom_local_axes.atom0", "atom_local_axes.ax1",
        "atom_local_axes.atom2", "atom_local_axes.ax2"}
    assert d_question["atom_type_symbol"]["proposal"] == "F: O2-"
    assert d_question["orbital"]["proposal"] == "F: 2p"
    assert d_question["atom_local_axes.atom0"]["proposal"] == "Mn: D2, F: D1"
    assert d_question["atom_local_axes.ax1"]["proposal"] == "Mn: Z, F: X"
    # no second direction near F that is not along D1: no proposal, asked
    assert "proposal" not in d_question["atom_local_axes.atom2"]
    assert "add_multipole_form_factor of crystal_mnf2" in d_question["atom_local_axes.atom0"]["hint"]
    assert not crystal.is_attribute("atom_rho_multipole")
    d_answer = {name: item.get("proposal", "") for name, item in d_question.items()}
    d_answer["atom_local_axes.atom2"] = "Mn: F, F: D2"
    d_out = add_item(obj, atom_label="Mn and F", **d_answer)
    assert not d_out.get("missing"), d_out["text"]
    d_row = {item.atom_label: (item.atom_type_symbol, item.orbital, item.coeff_P00)
             for item in crystal.atom_rho_multipole.items}
    assert d_row == {"Mn": ("Mn3+", "3d", 1.), "F": ("O2-", "2p", 1.)}
    d_axes = {item.atom_label: (item.atom0, item.ax1, item.atom1, item.atom2, item.ax2)
              for item in crystal.atom_local_axes.items}
    assert d_axes == {"Mn": ("D2", "Z", "Mn", "F", "X"), "F": ("D1", "X", "F", "D2", "Z")}
    # the rows are there: nothing more to add
    assert "nothing is added" in add_item(obj, atom_label="Mn")["text"]
    text = add_item(obj)["text"]
    assert "no new row of atom_rho_multipole is possible" in text
    assert "D1, D2: no atom_site_scat" in text


def test_plain_answer_does_not_replace_a_sure_value():
    """'2p', the answer to the question about F, must not become the
    orbital of Mn (its 3d is sure)."""
    obj = cryspy.str_to_globaln(MNF2)
    crystal = obj.crystal_mnf2
    crystal.add_multipole_form_factor("Mn, F", magnetic_ion="F: O2-", orbital="F: 2p")
    obj.crystal_mnf2.atom_rho_multipole.items.clear()
    d_out = add_item(obj, atom_label="Mn, F", atom_type_symbol="O2-", orbital="2p")
    assert not d_out.get("missing"), d_out["text"]
    d_row = {item.atom_label: (item.atom_type_symbol, item.orbital)
             for item in crystal.atom_rho_multipole.items}
    assert d_row == {"Mn": ("Mn3+", "3d"), "F": ("O2-", "2p")}


def test_reference_values_are_checked_for_any_category():
    obj = cryspy.str_to_globaln(MNF2)
    text = cryspy.call_assistant_tool(obj, "add_item", {
        "category": "atom_site_scat", "label": "Fe"})["text"]
    assert "no Fe in atom_site.label" in text
    text = cryspy.call_assistant_tool(obj, "add_item", {
        "category": "atom_site_scat", "label": "D1"})["text"]
    # the type symbol comes from the row of atom_site with the same label
    assert "type_symbol Dummy" in text


def test_category_report_names_the_methods_that_compute_it():
    obj = cryspy.str_to_globaln(MNF2)
    text = cryspy.call_assistant_tool(obj, "check_requirements", {"name": "atom_local_axes"})["text"]
    assert "one row per atom_site.label (atom_label); present: none" in text
    assert "also made by add_multipole_form_factor of crystal_mnf2" in text
    text = form_matching_methods_text(obj, REQUEST)
    assert 'Without the method: add_item {"category": "atom_rho_multipole"}' in text


def test_help_says_what_multipoles_need():
    """Iurii: the assistant must know that multipoles need local axes and
    whether the core part of the electron shell is needed."""
    from cryspy.procedure_assistant.help import get_help
    from cryspy.C_item_loop_classes.cl_1_atom_rho_multipole import AtomRhoMultipole
    text = get_help("magnetic form factor through multipoles core electron configuration")
    assert "local axes" in text and "atom_local_axes" in text
    assert "core electrons" in text and "not needed for the magnetic form factor" in text
    assert "local axes of the atom" in AtomRhoMultipole.__doc__
