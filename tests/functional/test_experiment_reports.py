"""Reports of the experiments shown in the View tab of CrysPy Editor."""
import os
import io
import contextlib

import cryspy

DIR_EXAMPLES = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "examples")


def _load_and_calc(*l_path):
    rcif_object = cryspy.load_file(os.path.join(DIR_EXAMPLES, *l_path))
    with contextlib.redirect_stdout(io.StringIO()):
        dict_res = cryspy.rhochi_no_refinement(rcif_object)
    return rcif_object, dict_res


def test_report_powder_1d_unpolarized():
    rcif_object, dict_res = _load_and_calc(
        "rhochi_unpolarized_powder_1d_PbSO4",
        "rhochi_unpolarized_powder_1d_PbSO4.rcif")
    pd = rcif_object.pd_pd
    s_md = pd.report()
    assert "# Agreement with the model" in s_md
    assert "# Phases" in s_md
    chi_sq_per_n = float(dict_res["chi_sq"]) / float(dict_res["n_point"])
    assert f"| intensity | yes | {int(dict_res['n_point']):} | {chi_sq_per_n:.2f} |" in s_md
    s_html = pd.report_html()
    assert s_html.count("<table") == s_html.count("</table>") == 3


def test_report_tof():
    rcif_object, dict_res = _load_and_calc(
        "rhochi_unpolarized_tof_powder_CeCuAl",
        "rhochi_unpolarized_tof_powder_CeCuAl.rcif")
    s_md = rcif_object.tof_tof1.report()
    chi_sq_per_n = float(dict_res["chi_sq"]) / float(dict_res["n_point"])
    assert f"{chi_sq_per_n:.2f}" in s_md
    assert "# Refinement" in s_md


def test_report_single_crystal():
    rcif_object, dict_res = _load_and_calc(
        "rhochi_single_susceptibility_Ho2Ti2O7",
        "rhochi_single_susceptibility_Ho2Ti2O7.rcif")
    diffrn = rcif_object.diffrn_ho2ti2o7
    s_html = diffrn.report_html()
    assert "Agreement with the model" in s_html
    assert s_html.count("<table") == s_html.count("</table>")
    assert s_html.count("<tr>") == s_html.count("</tr>")
    # inner tags of the html are closed: no '<' from the text is left
    assert " < " not in s_html
