"""hkl labels of peak positions on interactive powder plots."""
import os

import matplotlib

import cryspy

matplotlib.use("Agg")

F_NAME_PD = os.path.join(
    os.path.dirname(__file__), "Ho2Ti2O7_5K_1T_powder_1d_test", "main.rcif")
F_NAME_TOF = os.path.join(os.path.dirname(__file__), "test_tof", "model_pV.rcif")


def _plot(f_name, class_name):
    rcif_obj = cryspy.load_file(f_name)
    cryspy.rhochi_no_refinement(rcif_obj)
    data = [item for item in rcif_obj.items if type(item).__name__ == class_name][0]
    fig, ax = data.plots()[0]
    return data, fig, ax


def _labels(ax):
    return [t for t in ax.texts if t.get_text().startswith("(")]


def _check_zoom(data, fig, ax, peak_class_name, attr_x):
    peaks = [item for item in data.items if type(item).__name__ == peak_class_name][0]
    x = _sorted(getattr(peaks, attr_x))
    assert len(x) > 5
    # the full pattern has many peaks: no labels
    assert _labels(ax) == []

    # zoom on the first peak only: its label is shown
    dx = 0.1 * min(abs(x[1] - x[0]), 1.0) + 1e-6
    ax.set_xlim(x[0] - dx, x[0] + dx)
    l_text = [t.get_text() for t in _labels(ax)]
    assert len(l_text) >= 1
    i = list(getattr(peaks, attr_x)).index(x[0])
    label = f"({peaks.index_h[i]} {peaks.index_k[i]} {peaks.index_l[i]})"
    assert any(label in text for text in l_text)

    # zoom on the lower axes (shared x) also updates the labels
    ax_hkl = fig.axes[1]
    ax_hkl.set_xlim(x[0] - 1e-3, x[-1] + 1e-3)
    assert _labels(ax) == []
    fig.canvas.draw()


def _sorted(values):
    return sorted(float(v) for v in values)


def test_pd_hkl_labels_on_zoom():
    data, fig, ax = _plot(F_NAME_PD, "Pd")
    _check_zoom(data, fig, ax, "PdPeakL", "ttheta")


def test_tof_hkl_labels_on_zoom():
    data, fig, ax = _plot(F_NAME_TOF, "TOF")
    _check_zoom(data, fig, ax, "TOFPeakL", "time")


def test_group_peak_labels():
    from cryspy.A_functions_base.plot_peak_labels import group_peak_labels
    groups = group_peak_labels(
        [20.167, 24.171, 24.171, 20.17], ["(3 3 1)", "(3 3 3)", "(5 1 1)", "(3 3 1)"],
        ["C0", "C0", "C0", "C1"])
    assert groups == [(20.167, ["(3 3 1)"], "C0"), (24.171, ["(3 3 3)", "(5 1 1)"], "C0")]
