"""Plots of single-crystal experiments."""
import os

import matplotlib

import cryspy

matplotlib.use("Agg")

F_NAME_INTENSITY = os.path.join(os.path.dirname(__file__), "twin_tbte3", "main.rcif")
F_NAME_FLIP_RATIO = os.path.join(os.path.dirname(__file__), "HoTi_single_test", "main.rcif")


def _diffrn_plots(f_name):
    rcif_obj = cryspy.load_file(f_name)
    cryspy.rhochi_no_refinement(rcif_obj)
    diffrn = [item for item in rcif_obj.items if type(item).__name__ == "Diffrn"][0]
    return diffrn.plots()


def test_plots_unpolarised_intensity():
    # data given by intensity (no flip ratio) must not ask for 'fr'
    l_fig_ax = _diffrn_plots(F_NAME_INTENSITY)
    assert len(l_fig_ax) == 1


def test_plots_flip_ratio():
    l_fig_ax = _diffrn_plots(F_NAME_FLIP_RATIO)
    assert len(l_fig_ax) == 2
