import os
import re

import numpy

import cryspy

F_NAME = os.path.join(os.path.dirname(__file__), "model_gauss.rcif")


def _chi_sq_with_texture(g_1: float):
    with open(F_NAME, "r") as fid:
        text = fid.read()
    texture = (
        "\nloop_\n_texture_label\n_texture_g_1\n_texture_g_2\n"
        "_texture_h_ax\n_texture_k_ax\n_texture_l_ax\n"
        f"cecual {g_1:.5f} 0.00000 0.000 0.000 1.000\n")
    text, n_sub = re.subn(r"(\ndata_tof\w*\n)", lambda match: match.group(1) + texture, text, count=1)
    assert n_sub == 1
    rhochi = cryspy.str_to_globaln(text)
    res = cryspy.rhochi_no_refinement(rhochi)
    return res["chi_sq"]


def test_texture_in_tof_refinement():
    rhochi = cryspy.file_to_globaln(F_NAME)
    chi_sq_no_texture = cryspy.rhochi_no_refinement(rhochi)["chi_sq"]

    # g_1 = 1 is a random powder: the March-Dollase factor is 1 for every hkl.
    assert numpy.isclose(_chi_sq_with_texture(1.), chi_sq_no_texture, rtol=1e-10)

    chi_sq_texture = _chi_sq_with_texture(0.8)
    assert numpy.isfinite(chi_sq_texture)
    assert not numpy.isclose(chi_sq_texture, chi_sq_no_texture, rtol=1e-6)
