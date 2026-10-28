"""Labels of the inversed Hessian are the paths of the refinement log.

The labels were made from the ways in the dictionary
('t1p8nuc,scale,(np.int64(0),)') and did not match the names of the
optimal parameters ('diffrn_t1p8nuc.phase.scale').
"""
import os

import cryspy
from cryspy.procedure_rhochi.rhochi import rhochi_inversed_hessian
from cryspy.procedure_assistant.parameters import take_parameter_records

F_NAME = os.path.join(os.path.dirname(__file__), "twin_tbte3", "main.rcif")


def test_inversed_hessian_labels_are_paths():
    cryspy_object = cryspy.load_file(F_NAME)
    l_path = [record.path for record in take_parameter_records(cryspy_object)
              if record.flag_refinement]
    assert len(l_path) > 1

    inv_hessian = rhochi_inversed_hessian(cryspy_object)

    # the names and the order of the log of the refinement
    assert list(inv_hessian.label) == l_path
    # the labels survive the rcif text of the matrix
    inv_hessian_2 = cryspy.InversedHessian(with_labels=inv_hessian.with_labels)
    inv_hessian_2.form_object()
    assert list(inv_hessian_2.label) == list(inv_hessian.label)


def test_inversed_hessian_report_html():
    inv_hessian = cryspy.InversedHessian(
        with_labels="crystal_fe.cell.length_a 0.0001 0.00009\n"
                    "crystal_fe.atom_site[Fe1].occupancy 0.00009 0.0001")
    s_html = inv_hessian.report_html()
    assert "crystal_fe.atom_site[Fe1].occupancy" in s_html
    # each row of the tables is opened and closed
    assert s_html.count("<tr>") == s_html.count("</tr>") == 3 + 3
    # strong correlation (0.9) is marked
    assert "bgcolor" in s_html
    s_md = inv_hessian.report()
    assert "| 2 | crystal_fe.atom_site[Fe1].occupancy |" in s_md


def test_numerical_correlation_follows_the_result_order():
    from cryspy.procedure_assistant.tools import take_numerical_correlation
    inv_hessian = cryspy.InversedHessian(
        with_labels="b.y 4.0 1.0 0.0\na.x 1.0 1.0 0.0\nc.z 0.0 0.0 9.0")
    inv_hessian.form_object()
    cryspy_object = cryspy.load_file(F_NAME)
    cryspy_object.add_items([inv_hessian])
    result = {"parameter_name": [("a", "x", (0,)), ("b", "y", (0,)), ("c", "z", (0,))],
              "parameter_path": ["a.x", "b.y", "c.z"]}
    matrix = take_numerical_correlation(cryspy_object, result)["correlation_matrix"]
    assert abs(matrix[0, 1] - 0.5) < 1e-12
    assert abs(matrix[0, 2]) < 1e-12
