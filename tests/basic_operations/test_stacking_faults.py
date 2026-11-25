"""Stacking faults of layers in powder (A_functions_base/stacking_faults.py,
Hendricks-Teller model) and the diagnosis of the peak widths
(A_functions_base/peak_width_diagnostics.py)."""
import numpy
import pytest

from cryspy.A_functions_base.stacking_faults import (
    calc_atoms_of_cell, calc_fault_variants, calc_flag_faulted, calc_layer_cut,
    calc_layer_structure_factor, calc_orientations, calc_reciprocal_metric,
    calc_rod_intensity, calc_rod_width, find_stacking_operation, select_layer, split_normal)
from cryspy.A_functions_base.peak_width_diagnostics import (
    calc_lorentzian_columns, describe_shift, find_fault_rule, fit_sharp_and_broad,
    group_reflections)

UCP = numpy.array([10.36, 8.34, 13.39, 0.5 * numpy.pi, 0.5 * numpy.pi, 0.5 * numpy.pi])
# Pbcn: x,y,z; -x+1/2,-y+1/2,z+1/2; -x,y,-z+1/2; x+1/2,-y+1/2,-z and the inversion
PBCN = [(numpy.diag([1, 1, 1]), (0, 0, 0)), (numpy.diag([-1, -1, 1]), (.5, .5, .5)),
        (numpy.diag([-1, 1, -1]), (0, 0, .5)), (numpy.diag([1, -1, -1]), (.5, .5, 0)),
        (numpy.diag([-1, -1, -1]), (0, 0, 0)), (numpy.diag([1, 1, -1]), (.5, .5, .5)),
        (numpy.diag([1, -1, 1]), (0, 0, .5)), (numpy.diag([-1, 1, 1]), (.5, .5, 0))]
ROTATIONS = numpy.array([op[0] for op in PBCN])
TRANSLATIONS = numpy.array([op[1] for op in PBCN], dtype=float)


def form_layer(xyz, metric):
    n = xyz.shape[0]
    return {"xyz": xyz, "amplitude": numpy.ones(n, dtype=complex), "b_iso": numpy.zeros(n),
            "beta": numpy.zeros((n, 3, 3)), "metric": metric}


def form_rod(layer, rotation, shifts, orientations, h, k, tau):
    q = numpy.stack([numpy.full_like(tau, h), numpy.full_like(tau, k), tau], axis=1)
    f_layer = numpy.array([calc_layer_structure_factor(q @ r, layer) for r in orientations])
    phase = numpy.array([[numpy.exp(2j * numpy.pi * ((q @ r) @ t)) for t in shifts]
                         for r in orientations])
    return f_layer, phase


def test_stacking_operation_layer_and_faulted_reflections():
    xyz, index_atom, _ = calc_atoms_of_cell(
        numpy.array([[0.], [0.1232], [0.25]]), ROTATIONS, TRANSLATIONS)
    assert xyz.shape == (4, 3) and numpy.all(index_atom == 0)
    assert split_normal((0, 0, 2))[1] == 2
    rotation, translation = find_stacking_operation((0, 0, 2), ROTATIONS, TRANSLATIONS)
    assert (0, 0, 2) @ translation == pytest.approx(1.)
    assert calc_orientations(rotation).shape[0] == 2
    rotation_1, translation_1 = find_stacking_operation((0, 0, 1), ROTATIONS, TRANSLATIONS)
    assert numpy.array_equal(rotation_1, numpy.eye(3)) and translation_1[2] == pytest.approx(1.)
    # Fe at z = 1/4 and 3/4: the cut is in the gap at z = 0 (N.r = 0) or 1/2
    cut = calc_layer_cut(xyz, (0, 0, 2))
    assert min(abs(cut), abs(cut - 1.)) < 1e-6
    xyz_layer = select_layer(xyz, (0, 0, 2), cut)[0]
    assert xyz_layer.shape[0] == 2
    variants = calc_fault_variants((0., 0.5, 0.), (0, 0, 2), ROTATIONS, TRANSLATIONS)
    assert variants.shape[0] == 1
    index_hkl = numpy.array([[0, 0, 2], [1, 1, 0], [1, 2, 1], [1, 1, 1]]).T
    flag = calc_flag_faulted(index_hkl, calc_orientations(rotation), variants)
    assert list(flag) == [False, True, False, True]


def test_hexagonal_shift_has_two_variants():
    rotations = numpy.array([numpy.eye(3, dtype=int), [[0, -1, 0], [1, -1, 0], [0, 0, 1]],
                             [[-1, 1, 0], [-1, 0, 0], [0, 0, 1]]])
    translations = numpy.zeros((3, 3))
    variants = calc_fault_variants((1. / 3., 2. / 3., 0.), (0, 0, 1), rotations, translations)
    assert variants.shape[0] == 1  # (1/3, 2/3) is kept by the 3-fold axis
    variants = calc_fault_variants((1. / 3., 0., 0.), (0, 0, 1), rotations, translations)
    assert variants.shape[0] == 3


def test_random_faults_of_a_translation_stacking_are_the_lorentzian_formula():
    """Identical layers, S0 = translation c: I = |F|^2 (1 - |b|^2) / |1 - b e^(2 pi i l)|^2."""
    metric = calc_reciprocal_metric(UCP)[0]
    layer = form_layer(numpy.array([[0.1, 0.2, 0.3], [0.4, 0.1, 0.6]]), metric)
    alpha, shift = 0.2, numpy.array([0., 0.5, 0.])
    orientations = numpy.array([numpy.eye(3, dtype=int)])
    shifts = [numpy.array([0., 0., 1.]), numpy.array([0., 0.5, 1.])]
    tau = numpy.linspace(-0.5, 0.5, 201)
    f_layer, phase = form_rod(layer, None, shifts, orientations, 1, 1, tau)
    intensity = calc_rod_intensity(f_layer, phase, alpha, alpha, 1)
    b = 1. - alpha + alpha * numpy.exp(2j * numpy.pi * (shift @ [1, 1, 0]))
    expected = numpy.abs(f_layer[0]) ** 2 * (1. - abs(b) ** 2) / numpy.abs(
        1. - b * numpy.exp(2j * numpy.pi * tau)) ** 2
    assert numpy.allclose(intensity, expected, rtol=1e-8)
    # the width of the rod
    width = calc_rod_width(phase[:, :, 100:101], alpha, alpha, 1)[0]
    assert width == pytest.approx(-numpy.log(abs(b)) / (2. * numpy.pi))


@pytest.mark.parametrize("alpha, beta", [(0.2, 0.2), (0.1, 0.6), (0.3, 0.05)])
def test_rod_keeps_the_intensity_of_the_bragg_reflections(alpha, beta):
    """Layers related by the c-glide (x, -y, z+1/2), two per cell."""
    metric = calc_reciprocal_metric(UCP)[0]
    xyz_cell, _, _ = calc_atoms_of_cell(
        numpy.array([[0.0, 0.2633, 0.73], [0.1232, 0.1269, 0.1259], [0.25, 0.3509, 0.5013]]),
        ROTATIONS, TRANSLATIONS)
    cut = calc_layer_cut(xyz_cell, (0, 0, 2))
    xyz_layer = select_layer(xyz_cell, (0, 0, 2), cut)[0]
    rotation, translation = find_stacking_operation((0, 0, 2), ROTATIONS, TRANSLATIONS)
    orientations = calc_orientations(rotation)
    shifts = [translation, translation + numpy.array([0., 0.5, 0.])]
    # the intensity moves between the reflections of a rod, the total over
    # many of them is kept (no displacements: no decay along the rod)
    tau = numpy.linspace(-20., 20., 200001)[:-1]
    f_layer, phase = form_rod(form_layer(xyz_layer, metric), rotation, shifts, orientations,
                              1, 1, tau)
    intensity = calc_rod_intensity(f_layer, phase, alpha, beta, 2)
    assert numpy.all(intensity > -1e-9)
    hkl = numpy.array([[1., 1., float(l)] for l in range(-20, 20)])
    f2_bragg = numpy.abs(calc_layer_structure_factor(hkl, form_layer(xyz_cell, metric))) ** 2
    assert intensity.mean() == pytest.approx(f2_bragg.mean(), rel=0.03)


def test_small_probability_gives_bragg_peaks():
    metric = calc_reciprocal_metric(UCP)[0]
    layer = form_layer(numpy.array([[0.1, 0.2, 0.3]]), metric)
    orientations = numpy.array([numpy.eye(3, dtype=int)])
    shifts = [numpy.array([0., 0., 1.]), numpy.array([0., 0.5, 1.])]
    tau = numpy.array([0., 0.05, 0.5])
    f_layer, phase = form_rod(layer, None, shifts, orientations, 1, 1, tau)
    intensity = calc_rod_intensity(f_layer, phase, 1e-4, 1e-4, 1)
    assert intensity[0] > 1e3 * intensity[1] > 1e3 * intensity[2] * 0.


def test_describe_shift():
    assert describe_shift(numpy.array([0., 0.5, 0.])) == "k odd"
    assert describe_shift(numpy.array([0.5, 0.5, 0.])) == "h + k odd"
    assert describe_shift(numpy.array([1. / 3., 2. / 3., 0.])) == "h - k not 3n"


def test_fit_sharp_and_broad_and_the_rule():
    x = numpy.linspace(10., 60., 2001)
    index_hkl = numpy.array([[0, 0, 2], [1, 1, 0], [1, 1, 1], [2, 0, 0], [1, 1, 2], [0, 2, 1],
                             [1, 2, 1], [1, 3, 0], [2, 0, 2]]).T
    position = numpy.linspace(15., 55., index_hkl.shape[1])
    fwhm = numpy.full(position.size, 0.3)
    flag_broad = index_hkl[1] % 2 == 1
    width_true = numpy.where(flag_broad, 4. * fwhm, fwhm)
    y = 1000. * calc_lorentzian_columns(x, position, width_true).sum(axis=1) + 50.
    sigma = numpy.sqrt(y)
    y = y + numpy.random.default_rng(2).normal(0., 1., x.size) * sigma - 50.
    column_sharp = calc_lorentzian_columns(x, position, fwhm)
    column_broad = [calc_lorentzian_columns(x, position, f * fwhm) for f in (2.5, 5.)]
    i_sharp, i_broad, s_int, _ = fit_sharp_and_broad(y, sigma, column_sharp, column_broad, x)
    part = i_broad / (i_sharp + i_broad)
    assert numpy.all(part[flag_broad] > 0.8) and numpy.all(part[~flag_broad] < 0.2)
    l_rule = find_fault_rule([index_hkl[:, [i]] for i in range(position.size)], part,
                             numpy.ones(position.size))
    assert l_rule[0]["rule"] == "k odd"
    assert numpy.allclose(l_rule[0]["shift"], (0., 0.5, 0.))


def test_group_reflections():
    group = group_reflections(numpy.array([10., 10.01, 12., 20.]), numpy.full(4, 0.2))
    assert list(group) == [0, 0, 1, 2]
