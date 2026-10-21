"""Progress callback of the RhoChi refinement: report and stop."""
import os

import numpy

import cryspy
from cryspy.procedure_rhochi.rhochi_by_dictionary import \
    rhochi_rietveld_refinement_by_dictionary

F_NAME = os.path.join(
    os.path.dirname(__file__), "HoTi_single_test", "main.rcif")

INFO_KEYS = {"iteration", "chi_sq", "n_point", "chi_sq_per_n",
             "parameter_names", "parameter_values"}


def test_callback_reports_every_iteration():
    rhochi = cryspy.file_to_globaln(F_NAME)
    l_info = []
    dict_out = cryspy.rhochi_rietveld_refinement(
        rhochi, callback=lambda info: l_info.append(info))
    assert not dict_out["stopped"]
    assert len(l_info) > 2
    assert [info["iteration"] for info in l_info] == list(range(1, len(l_info) + 1))
    for info in l_info:
        assert set(info.keys()) == INFO_KEYS
        assert len(info["parameter_values"]) == len(info["parameter_names"])
        assert numpy.isclose(info["chi_sq_per_n"], info["chi_sq"] / info["n_point"])
    # BFGS accepts only steps which do not increase chi_sq
    l_chi_sq = [info["chi_sq"] for info in l_info]
    assert l_chi_sq[-1] <= l_chi_sq[0]


def test_callback_returning_false_stops_refinement():
    rhochi = cryspy.file_to_globaln(F_NAME)
    n_stop = 2
    l_info = []

    def callback(info):
        l_info.append(info)
        return info["iteration"] < n_stop

    dict_out = cryspy.rhochi_rietveld_refinement(rhochi, callback=callback)
    assert dict_out["stopped"]
    assert len(l_info) == n_stop

    # the object keeps the parameters of the last completed iteration:
    # compare with an object which takes the same values from its dictionary
    info = l_info[-1]
    reference = cryspy.file_to_globaln(F_NAME)
    reference_dict = reference.get_dictionary()
    for way, value in zip(info["parameter_names"], info["parameter_values"]):
        reference_dict[way[0]][way[1]][way[2]] = value
    reference.take_parameters_from_dictionary(reference_dict)
    global_dict = rhochi.get_dictionary()
    reference_dict = reference.get_dictionary()
    for way in info["parameter_names"]:
        assert numpy.allclose(
            global_dict[way[0]][way[1]], reference_dict[way[0]][way[1]],
            rtol=1e-12, atol=0.)

    chi_sq = cryspy.rhochi_no_refinement(rhochi)["chi_sq"]
    chi_sq_reference = cryspy.rhochi_no_refinement(reference)["chi_sq"]
    assert numpy.isclose(chi_sq, chi_sq_reference, rtol=1e-10)


def test_callback_stops_basinhopping():
    global_dict = cryspy.file_to_globaln(F_NAME).get_dictionary()
    chi_sq, parameter_names, dict_in_out, res = rhochi_rietveld_refinement_by_dictionary(
        global_dict, method="basinhopping niter=2",
        progress_callback=lambda info: False)
    assert res["stopped"]
    assert res["nit"] == 1
    values = [global_dict[way[0]][way[1]][way[2]] for way in parameter_names]
    assert numpy.allclose(values, res.x, rtol=1e-12, atol=0.)


def test_signature_stays_one_click_in_cryspy_editor():
    # CrysPy Editor runs a procedure without a dialog when its positional
    # arguments are only annotated GlobalN objects; callback must not change it
    for func in (cryspy.rhochi_rietveld_refinement, ):
        code = func.__code__
        assert code.co_argcount == 1
        annotations = {key for key in func.__annotations__ if key != "return"}
        assert annotations == {"cryspy_object"}
