"""Small parts of the conversion of single crystal FullProf files: Fortran
formats of .int files, magnetic form factors and the extinction."""
import numpy
import pytest

from cryspy.H_functions_global.fullprof import FullProfFormatError
from cryspy.H_functions_global.fullprof_single_crystal import (
    FullProfSingleCrystalPhase, calc_metric_reciprocal, form_extinction_rcif,
    form_fullprof_magnetic_type, parse_fortran_format, read_fixed_line, take_magnetic_ion)


def test_parse_fortran_format():
    assert parse_fortran_format("(3i4,i5,2F14.4,i5,4f8.2)") == (
        [("i", 4)] * 3 + [("i", 5)] + [("f", 14)] * 2 + [("i", 5)] + [("f", 8)] * 4)


def test_read_fixed_line():
    """Fixed columns: numbers without spaces between them are separated."""
    l_field = parse_fortran_format("(3i4,i5,2F14.4,i5)")
    line = "   1  -2  10    1      288.8600        1.4100    1"
    assert read_fixed_line(line, l_field) == [1., -2., 10., 1., 288.86, 1.41, 1.]
    assert read_fixed_line("-10-11-12", parse_fortran_format("(3i3)")) == [-10., -11., -12.]


def test_magnetic_ion():
    assert take_magnetic_ion("MHO3") == "Ho3+"
    assert form_fullprof_magnetic_type("Ho3+") == "MHO3"
    assert form_fullprof_magnetic_type("Mn2+") == "MMN2"
    with pytest.raises(FullProfFormatError):
        take_magnetic_ion("JHO3")


def test_extinction_isotropic_as_qani():
    """FullProf model 1, y = (1 + 0.001 E F² λ³/sin 2θ)^-1/2, equals qani,
    y = (1 + 0.001/4 F² λ⁵/(sin³θ cos θ) |q_ij h_i h_j|)^-1/2, with
    q = E/2 g*_ii and E g*_ij: (2 sinθ/λ)² = g*_ij h_i h_j."""
    phase = FullProfSingleCrystalPhase()
    phase.cell = [4.26, 25.16, 4.28, 90., 97., 90.]
    phase.extinction = [0.5, 0., 0., 0., 0., 0., 0.]
    phase.extinction_codes = [0.] * 7
    phase.extinction_model = 1
    d_q = dict(line.split()[0:2] for line in form_extinction_rcif(phase)[0][1:] if line.strip())
    g_r = calc_metric_reciprocal(phase.cell)
    hkl = numpy.array([1., 2., 3.])
    sum_q = sum(float(d_q[f"_extinction_q_{name}"]) * hkl[i] * hkl[j] for name, i, j in (
        ("hh", 0, 0), ("kk", 1, 1), ("ll", 2, 2), ("hk", 0, 1), ("hl", 0, 2), ("kl", 1, 2)))
    wavelength = 2.36
    sthovl = 0.5 * numpy.sqrt(hkl @ g_r @ hkl)
    sin_theta = sthovl * wavelength
    cos_theta = numpy.sqrt(1. - sin_theta ** 2)
    y_fullprof = 0.001 * 0.5 * wavelength ** 3 / (2. * sin_theta * cos_theta)
    y_cryspy = 0.001 / 4. * wavelength ** 5 / (sin_theta ** 3 * cos_theta) * sum_q
    assert y_cryspy == pytest.approx(y_fullprof, rel=1e-6)
