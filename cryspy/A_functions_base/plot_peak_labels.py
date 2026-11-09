"""Interactive hkl labels of peak positions on matplotlib plots.

The functions work on matplotlib objects given as arguments and do not
import matplotlib, so `import cryspy` works without the `cryspy[plot]`
extra.
"""
import numpy


def form_hkl_labels(index_h, index_k, index_l) -> list:
    """Form labels "(h k l)" of reflections.

    Parameters
    ----------
    index_h, index_k, index_l : array-like of int, shape (n_hkl,)
        Miller indices.

    Returns
    -------
    list of str
        Labels as in the single-crystal plots, for example "(1 0 -2)".
    """
    return [f"({h:} {k:} {l:})" for h, k, l in zip(index_h, index_k, index_l)]


def group_peak_labels(positions, labels, colors=None, decimals: int = 2) -> list:
    """Group labels of peaks that have the same position.

    Parameters
    ----------
    positions : array-like of float, shape (n_peak,)
        Peak positions (2theta in degrees, time in microseconds, ...).
    labels : sequence of str, shape (n_peak,)
        Label of every peak.
    colors : sequence, shape (n_peak,), optional
        Color of every peak (a matplotlib color or None).
    decimals : int
        Positions equal after rounding to `decimals` form one group.

    Returns
    -------
    list of tuple
        `(position, labels, color)` per group, in the order of `positions`:
        the position and the color of the first peak of the group and the
        list of its distinct labels.
    """
    if colors is None:
        colors = [None] * len(labels)
    np_x = numpy.asarray(positions, dtype=float)
    groups = {}
    for key, position, label, color in zip(numpy.round(np_x, decimals), np_x, labels, colors):
        group = groups.setdefault(float(key), (float(position), [], color))
        if label not in group[1]:
            group[1].append(label)
    return list(groups.values())


def add_interactive_peak_labels(
        ax, l_peaks, x_curve=None, y_curve=None, n_visible: int = 5):
    """Show hkl labels of peaks when few peaks are in the visible x-range.

    The labels follow zoom and pan of the interactive plot: they are shown
    when at most `n_visible` peak positions are inside the current x-range
    and hidden when there are more, as for the single-crystal plots
    (`plot_scatters_with_labels`). Every Axes sharing the x-axis with `ax`
    triggers the update.

    Peaks of all tuples with the same position (rounded to 0.01) count as
    one and give one label listing their distinct labels: a data block may
    hold the same peaks twice (an old unnamed loop read from the file and
    the loop calculated for the phase).

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes where the labels are drawn.
    l_peaks : list of tuple
        One tuple `(positions, labels, color)` per phase: positions of the
        peaks along x (array of float, shape (n_peak,)), their labels
        (sequence of str, shape (n_peak,)) and the text color (a matplotlib
        color or None for black).
    x_curve, y_curve : array-like of float, shape (n_point,), optional
        Curve above which a label is put at the peak position (for example
        the upper envelope of the experiment and the model). Without a
        curve the labels are put at the bottom of `ax`.
    n_visible : int
        Maximal number of visible peak positions for which labels are shown.

    Returns
    -------
    callable
        The update function, called once here; call it again after
        changing the limits without emitting events.
    """
    positions, labels, colors = [], [], []
    for np_x, l_label, color in l_peaks:
        positions.extend(numpy.asarray(np_x, dtype=float))
        labels.extend(l_label)
        colors.extend([color] * len(l_label))
    l_groups = []
    for position, labs, color in group_peak_labels(positions, labels, colors):
        chunks = [labs[i:i+5] for i in range(0, len(labs), 5)]
        text = "\n".join(", ".join(chunk) for chunk in chunks)
        l_groups.append((position, text, color))

    if x_curve is not None and y_curve is not None:
        x_curve = numpy.asarray(x_curve, dtype=float)
        y_curve = numpy.asarray(y_curve, dtype=float)
        order = numpy.argsort(x_curve)
        x_curve, y_curve = x_curve[order], y_curve[order]
        flag_curve = x_curve.size > 0
    else:
        flag_curve = False

    active_texts = []

    def update_labels(ax_changed=None):
        """Recreate labels for the current limits.

        Called when the axis limits change (zoom, pan, set_xlim), before
        the figure is drawn, so the labels appear in the same redraw.
        """
        nonlocal active_texts
        for t in active_texts:
            t.remove()
        active_texts = []

        xmin, xmax = sorted(ax.get_xlim())
        visible = [group for group in l_groups if xmin <= group[0] <= xmax]
        if len(visible) == 0 or len(visible) > n_visible:
            return

        ymin, ymax = ax.get_ylim()
        for position, text, color in visible:
            # above the curve; for a high peak the label hangs from the top
            # of the axes, right of the peak position
            xy, xycoords = (position, 0.02), ("data", "axes fraction")
            xytext, ha, va = (0, 4), "center", "bottom"
            if flag_curve:
                y = float(numpy.interp(position, x_curve, y_curve))
                if y <= ymin + 0.6 * (ymax - ymin):
                    xy, xycoords = (position, max(y, ymin)), "data"
                else:
                    xy = (position, 1.0)
                    xytext, ha, va = (2, -4), "left", "top"
            t = ax.annotate(
                text, xy=xy, xycoords=xycoords, xytext=xytext,
                textcoords="offset points", rotation=90, ha=ha,
                va=va, fontsize=9,
                color="black" if color is None else color,
                annotation_clip=False, clip_on=True)
            active_texts.append(t)

    l_ax = ax.get_shared_x_axes().get_siblings(ax)
    for ax_shared in l_ax:
        ax_shared.callbacks.connect("xlim_changed", update_labels)
    ax.callbacks.connect("ylim_changed", update_labels)

    update_labels()
    return update_labels
