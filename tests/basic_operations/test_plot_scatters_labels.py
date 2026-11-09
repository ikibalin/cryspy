"""hkl labels of the scatter plots follow the zoom at once."""
import matplotlib
matplotlib.use("Agg")

from cryspy.C_item_loop_classes.cl_1_diffrn_refln import \
    plot_scatters_with_labels


def texts(fig):
    return [t.get_text().strip() for t in fig.axes[0].texts]


def test_label_appears_at_the_first_zoom():
    x = [0.1 * i for i in range(20)]
    labels = [f"({i} 0 0)" for i in range(20)]
    fig = plot_scatters_with_labels(x, x, sy=[0.01] * 20, labels=labels,
                                    n_visible=10)
    ax = fig.axes[0]
    ax.set_xlim(-1, 3)
    ax.set_ylim(-1, 3)
    assert texts(fig) == []
    # zoom of the toolbar to one point: both limits are set, then drawn
    ax.set_xlim(0.45, 0.55)
    ax.set_ylim(0.45, 0.55)
    assert texts(fig) == ["(5 0 0)"]
    fig.canvas.draw()
    assert texts(fig) == ["(5 0 0)"]
    ax.set_xlim(-1, 3)
    ax.set_ylim(-1, 3)
    assert texts(fig) == []
