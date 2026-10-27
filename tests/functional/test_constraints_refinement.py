"""Refinement with a linear constraint between two parameters (x+y=b)."""
import logging
import os

import numpy

import cryspy
from cryspy.procedure_assistant.parameters import take_dict_paths, take_parameter_records
from cryspy.procedure_rhochi.rhochi import log_refined_parameters

F_NAME = os.path.join(os.path.dirname(__file__), "twin_tbte3", "main.rcif")


def test_refinement_keeps_linear_constraint(caplog):
    obj = cryspy.load_file(F_NAME)
    # only the two coordinates: the extinction radius 0 of the file makes
    # the Hessian of all parameters singular
    obj.fix_variables()
    cryspy.add_constraint(obj, "y1 + y2 = 1.5", {
        "y1": "tbte3.atom_site[Te1].fract_y", "y2": "tbte3.atom_site[Te2].fract_y"})
    atom_site = obj.crystal_tbte3.atom_site
    # the start values do not satisfy the constraint
    assert not numpy.isclose(atom_site["Te1"].fract_y + atom_site["Te2"].fract_y, 1.5)

    with caplog.at_level(logging.INFO, logger="cryspy"):
        cryspy.rhochi_rietveld_refinement(obj)
    l_message = [record.getMessage() for record in caplog.records]
    # fitting and constrained parameters are logged by their full paths
    i_fit = l_message.index("Number of fitting parameters 1")
    assert l_message[i_fit + 1].startswith(" - crystal_tbte3.atom_site[Te1].fract_y ")
    assert " - crystal_tbte3.atom_site[Te2].fract_y" in l_message
    # optimal parameters are logged by their full paths and marks
    assert any(record.getMessage().startswith(" - crystal_tbte3.atom_site[Te2].fract_y ")
               and record.getMessage().endswith("[y2]") for record in caplog.records)

    y_1, y_2 = atom_site["Te1"].fract_y, atom_site["Te2"].fract_y
    assert numpy.isclose(y_1 + y_2, 1.5, atol=1e-10)
    assert abs(y_1 - 0.92893) < 0.01
    # dependent parameter: sigma from the covariance of the free one
    assert atom_site["Te2"].fract_y_sigma > 0.
    assert numpy.isclose(atom_site["Te1"].fract_y_sigma, atom_site["Te2"].fract_y_sigma)


def test_refinement_with_soft_constraint():
    def refine(sigma):
        obj = cryspy.load_file(F_NAME)
        obj.fix_variables()
        cryspy.add_constraint(obj, "y1 + y2 = 1.5", {
            "y1": "tbte3.atom_site[Te1].fract_y", "y2": "tbte3.atom_site[Te2].fract_y"},
            sigma=sigma)
        cryspy.rhochi_rietveld_refinement(obj)
        atom_site = obj.crystal_tbte3.atom_site
        return atom_site["Te1"].fract_y + atom_site["Te2"].fract_y

    # the data alone move the sum to 1.569 (the other parameters are fixed,
    # chi_sq ~ 1e7): a tight restraint keeps the start sum 1.5, a loose one
    # lets the data decide
    assert abs(refine(1e-5) - 1.5) < 1e-3
    assert refine(1e-2) > 1.55


def test_dict_paths_of_refined_parameters():
    obj = cryspy.load_file(F_NAME)
    obj.fix_variables()
    crystal = obj.crystal_tbte3
    crystal.atom_site["Te1"].occupancy_refinement = True
    crystal.atom_site["Te2"].fract_y_refinement = True
    crystal.cell.length_c_refinement = True
    value_before = crystal.atom_site["Te1"].occupancy
    d_path = take_dict_paths(obj, obj.get_dictionary())
    assert sorted(d_path.values()) == [
        "crystal_tbte3.atom_site[Te1].occupancy", "crystal_tbte3.atom_site[Te2].fract_y",
        "crystal_tbte3.cell.length_c"]
    # the ways are those of the refinement
    assert d_path[("crystal_tbte3", "unit_cell_parameters", (2, ))] == "crystal_tbte3.cell.length_c"
    # values are restored after the search
    assert numpy.isclose(crystal.atom_site["Te1"].occupancy, value_before, rtol=0., atol=1e-12)


def test_log_with_text_values(caplog):
    # loops unknown to the data block (atom_site_susceptibility of a diffrn
    # block) keep the values written by a refinement as text
    f_name = os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "examples",
        "rhochi_single_susceptibility_Ho2Ti2O7", "rhochi_single_susceptibility_Ho2Ti2O7.rcif")
    obj = cryspy.load_file(f_name)
    dict_global = obj.get_dictionary()
    record = [record for record in take_parameter_records(obj)
              if record.path == "diffrn_ho2ti2o7.atom_site_susceptibility[Ho1].chi_11"][0]
    setattr(record.item, record.attribute, "3.329")
    d_path = take_dict_paths(obj, dict_global)
    assert "crystal_ho2ti2o7.atom_site_susceptibility[Ho1].chi_11" in d_path.values()
    assert record.value == "3.329"
    with caplog.at_level(logging.INFO, logger="cryspy"):
        log_refined_parameters(obj)
    assert " - diffrn_ho2ti2o7.atom_site_susceptibility[Ho1].chi_11 3.329" in \
        [record.getMessage() for record in caplog.records]
