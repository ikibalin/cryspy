"""Merging equivalent magnetic satellites and their labels."""
import numpy

from cryspy.A_functions_base.plot_peak_labels import form_fractional_hkl_labels
from cryspy.A_functions_base.structure_factor_magnetic_k import (
    merge_equivalent_satellites)


def test_merge_equivalent_satellites():
    # (±1/2 0 ±1/2) of k1: equal 2theta and intensity in pairs; (1/2 0 1/2)
    # of k2 at the same position stays separate
    index_hkl = numpy.array([
        [0.5, -0.5, 0.5, -0.5, 0.5],
        [0.0, 0.0, 0.0, 0.0, 0.0],
        [0.5, -0.5, -0.5, 0.5, 0.5]])
    ttheta = numpy.radians([20.0, 20.0, 20.01, 20.01, 20.0])
    iint = numpy.array([3.0, 3.0 * (1. + 1e-9), 0.0, 0.0, 1.0])
    k_index = numpy.array([0, 0, 0, 0, 1])

    index_m, ttheta_m, iint_m, k_m, mult = merge_equivalent_satellites(
        index_hkl, ttheta, iint, k_index)

    assert index_m.shape == (3, 3)
    rows = {(tuple(index_m[:, i]), int(k_m[i])): (iint_m[i], mult[i])
            for i in range(3)}
    assert numpy.allclose(rows[((0.5, 0.0, 0.5), 0)], (6.0, 2))
    assert numpy.allclose(rows[((0.5, 0.0, -0.5), 0)], (0.0, 2))
    assert numpy.allclose(rows[((0.5, 0.0, 0.5), 1)], (1.0, 1))
    assert numpy.isclose(iint_m.sum(), iint.sum())


def test_merge_equivalent_satellites_empty():
    res = merge_equivalent_satellites(
        numpy.zeros((3, 0)), numpy.zeros(0), numpy.zeros(0),
        numpy.zeros(0, dtype=int))
    assert res[0].shape == (3, 0)
    assert all(item.size == 0 for item in res[1:])


def test_form_fractional_hkl_labels():
    labels = form_fractional_hkl_labels(
        [0.5, -0.0, 0.2341], [0.0, 1.0, 0.0], [-1.25, 2.0, 1.0])
    assert labels == ["(0.5 0 -1.25)", "(0 1 2)", "(0.234 0 1)"]
