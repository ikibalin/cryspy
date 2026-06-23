import numpy

na = numpy.newaxis


def cutoff_select(delta_2d, cutoff_fwhm, half_width, params):
    """Restrict each peak to its cutoff window (FullProf "WDT").

    Returns ``(None, params, delta_2d)`` when ``cutoff_fwhm`` is ``inf`` (no
    cutoff). Otherwise returns the boolean mask, the per-point ``params``
    gathered onto the in-window points, and the in-window deltas as a
    column, so the kernels run only there (faster for a tighter cutoff).
    """
    if numpy.isinf(cutoff_fwhm):
        return None, params, delta_2d
    keep = numpy.abs(delta_2d) <= half_width[:, na]
    row = numpy.nonzero(keep)[0]
    return keep, tuple(p[row] for p in params), delta_2d[keep][:, na]


def cutoff_place(keep, res):
    """Scatter column results back onto the full grid (pass-through if no cutoff)."""
    if keep is None:
        return res
    out = numpy.zeros(keep.shape)
    out[keep] = res[:, 0]
    return out
