"""Magnetic satellites of a powder pattern in the output and on the plot.

Ho2BaNiO5 (k = (1/2, 0, 1/2), D1B): after a calculation the
experiment holds the loop pd_peak_magnetic_k with one row per orbit of
satellites h + k, and the zoomed plot labels the satellite peaks.
"""
import os

import matplotlib
import numpy

import cryspy

matplotlib.use("Agg")

F_RCIF = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))))),
    "examples", "rhochi_unpolarized_powder_1d_k_Ho2BaNiO5",
    "rhochi_unpolarized_powder_1d_k_Ho2BaNiO5.rcif")
K_VECTOR = numpy.array([0.5, 0.0, 0.5])


def _calculate():
    obj = cryspy.load_file(F_RCIF)
    cryspy.rhochi_no_refinement(obj)
    pd = obj.pd_d1b
    satellites = [item for item in pd.items
                  if type(item).__name__ == "PdPeakMagneticKL"]
    return obj, pd, satellites


def test_satellite_loop_positions_and_indices():
    obj, pd, satellites = _calculate()
    assert len(satellites) == 1
    loop = satellites[0]
    assert loop.loop_name == "Ho2BaNiO5"
    index_hkl = numpy.stack([loop.numpy_index_h, loop.numpy_index_k,
                             loop.numpy_index_l])
    # every row is h + k or h - k with integer h
    integer_part = numpy.minimum(
        numpy.abs(index_hkl - K_VECTOR[:, None] - numpy.round(index_hkl - K_VECTOR[:, None])),
        numpy.abs(index_hkl + K_VECTOR[:, None] - numpy.round(index_hkl + K_VECTOR[:, None])))
    assert numpy.allclose(integer_part, 0., atol=1e-4)
    assert numpy.all(numpy.diff(loop.numpy_ttheta) >= 0.)
    assert numpy.all(loop.numpy_intensity > 0.)
    assert numpy.all(loop.numpy_index_mult >= 1)
    assert set(loop.numpy_wave_vector_seq_id) == {1}

    # 2theta of (1/2 0 1/2) from Bragg's law; the cell of the example is
    # monoclinic (beta = 89.95 deg), 1/d^2 = H G* H with the metric G = A A^T
    cell = obj.crystal_ho2banio5.cell
    beta = numpy.radians(cell.angle_beta)
    vectors = numpy.array([  # a, b, c in Cartesian coordinates (rows)
        [cell.length_a, 0., 0.],
        [0., cell.length_b, 0.],
        [cell.length_c * numpy.cos(beta), 0., cell.length_c * numpy.sin(beta)]])
    metric_reciprocal = numpy.linalg.inv(vectors @ vectors.transpose())
    index = numpy.array([0.5, 0., 0.5])
    inv_d = numpy.sqrt(index @ metric_reciprocal @ index)
    ttheta = numpy.degrees(2. * numpy.arcsin(0.5 * pd.setup.wavelength * inv_d))
    ttheta += pd.setup.offset_ttheta
    i_row = int(numpy.argmin(numpy.abs(index_hkl - numpy.array(
        [[0.5], [0.], [0.5]])).sum(axis=0)))
    assert numpy.allclose(index_hkl[:, i_row], [0.5, 0., 0.5])
    assert abs(loop.numpy_ttheta[i_row] - ttheta) < 2e-3


def test_satellite_loop_survives_cif():
    obj, pd, satellites = _calculate()
    obj_2 = cryspy.str_to_globaln(obj.to_cif())
    loops = [item for item in obj_2.pd_d1b.items
             if type(item).__name__ == "PdPeakMagneticKL"]
    assert len(loops) == 1
    assert numpy.allclose(loops[0].numpy_ttheta, satellites[0].numpy_ttheta)
    assert numpy.allclose(loops[0].numpy_index_h, satellites[0].numpy_index_h)


def test_satellite_labels_on_zoom():
    obj, pd, satellites = _calculate()
    fig, ax = pd.plots()[0]
    ttheta = satellites[0].numpy_ttheta[0]
    ax.set_xlim(ttheta - 0.2, ttheta + 0.2)
    labels = [t.get_text() for t in ax.texts]
    assert "(0.5 0 0.5)" in labels
