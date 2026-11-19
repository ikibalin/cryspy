"""Magnetic correlation length of powder phases (loop magnetic_correlation).

The magnetic and interference intensities of a phase get an extra
Lorentzian (half width 1/xi in Q); the nuclear intensity keeps its profile.
Old files (without the loop) give the same patterns as before.
"""
import os

import numpy
import pytest

import cryspy

DIR_TESTS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
F_MCIF_PD = os.path.join(DIR_TESTS, "functional", "mcif_test", "pd_test.rcif")
F_HO2TI2O7 = os.path.join(DIR_TESTS, "functional", "Ho2Ti2O7_5K_1T_powder_1d_test", "main.rcif")
F_HO2BANIO5 = os.path.join(
    os.path.dirname(DIR_TESTS), "examples", "rhochi_unpolarized_powder_1d_k_Ho2BaNiO5",
    "rhochi_unpolarized_powder_1d_k_Ho2BaNiO5.rcif")


def take_pd(obj):
    return [item for item in obj.items if isinstance(item, cryspy.Pd)][0]


def calc_pattern(obj):
    """ttheta (degrees), signal plus, signal minus without the background."""
    cryspy.rhochi_no_refinement(obj)
    pd_proc = take_pd(obj).pd_proc
    return (numpy.array(pd_proc.ttheta, dtype=float),
            numpy.array(pd_proc.intensity_plus_net, dtype=float),
            numpy.array(pd_proc.intensity_minus_net, dtype=float))


@pytest.mark.filterwarnings("ignore::UserWarning")
@pytest.mark.parametrize("f_name", [F_MCIF_PD, F_HO2TI2O7, F_HO2BANIO5])
def test_infinite_length_and_ellipsoid_of_one_length(f_name):
    """An infinite length changes nothing; the ellipsoid made from one
    length gives the pattern of this length (k = 0, polarized
    paramagnetic, k != 0 satellites)."""
    ttheta, plus_0, minus_0 = calc_pattern(cryspy.load_file(f_name))
    obj = cryspy.load_file(f_name)
    cryspy.add_magnetic_correlation_length(obj, 1e12)
    plus_inf, minus_inf = calc_pattern(obj)[1:]
    # pd_proc keeps 5 decimals
    assert numpy.allclose(plus_inf, plus_0, rtol=1e-6, atol=2e-5)
    assert numpy.allclose(minus_inf, minus_0, rtol=1e-6, atol=2e-5)

    obj = cryspy.load_file(f_name)
    cryspy.add_magnetic_correlation_length(obj, 30.)
    plus_30, minus_30 = calc_pattern(obj)[1:]
    assert numpy.abs(plus_30 + minus_30 - plus_0 - minus_0).max() > \
        1e-3 * (plus_0 + minus_0).max()
    obj = cryspy.load_file(f_name)
    cryspy.add_magnetic_correlation_length(obj, 30., anisotropic=True)
    plus_a, minus_a = calc_pattern(obj)[1:]
    assert numpy.allclose(plus_a, plus_30, rtol=1e-8, atol=2e-5)
    assert numpy.allclose(minus_a, minus_30, rtol=1e-8, atol=2e-5)


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_satellites_keep_their_integrated_intensity():
    """Ho2BaNiO5, k = (1/2, 0, 1/2): purely magnetic satellites get broader
    and lower, the integrated intensity of a satellite well inside the
    pattern stays."""
    ttheta, plus_0, minus_0 = calc_pattern(cryspy.load_file(F_HO2BANIO5))
    obj = cryspy.load_file(F_HO2BANIO5)
    cryspy.add_magnetic_correlation_length(obj, 40.)
    plus_40, minus_40 = calc_pattern(obj)[1:]
    signal_0, signal_40 = plus_0 + minus_0, plus_40 + minus_40
    assert signal_40.max() < signal_0.max()
    # first satellite, at about 20.4 degrees (lambda 2.524 A)
    peak = numpy.argmax(numpy.where(ttheta < 25., signal_0, 0.))
    window = numpy.abs(ttheta - ttheta[peak]) < 3.
    ratio = signal_40[window].sum() / signal_0[window].sum()
    # the Lorentzian tails outside +-3 degrees and of the neighbours: 10 %
    assert 0.85 < ratio < 1.1
    assert signal_40[peak] < 0.8 * signal_0[peak]


def test_round_trip_of_the_loop():
    obj = cryspy.load_file(F_MCIF_PD)
    cryspy.add_magnetic_correlation_length(obj, 25., anisotropic=True)
    pd = take_pd(obj)
    pd.magnetic_correlation.items[0].aniso_33 = 100.
    text = obj.to_cif()
    assert "_magnetic_correlation_aniso_33" in text
    obj_2 = cryspy.str_to_globaln(text)
    row = take_pd(obj_2).magnetic_correlation.items[0]
    assert row.aniso_33 == pytest.approx(100.)
    assert row.aniso_11_refinement and not row.aniso_12_refinement
    report = cryspy.report_magnetic_correlation(obj_2)
    assert "c* 10.0 A" in report
    assert "2.5 times shorter" in report


def test_refinement_recovers_the_correlation_length():
    """Pattern calculated with xi = 40 A, refined from 80 A."""
    obj = cryspy.load_file(F_MCIF_PD)
    obj.fix_variables()
    cryspy.add_magnetic_correlation_length(obj, 40.)
    obj.fix_variables()
    ttheta, plus, minus = calc_pattern(obj)
    signal = plus + minus + numpy.array(take_pd(obj).pd_proc.intensity_bkg_calc, dtype=float)

    obj = cryspy.load_file(F_MCIF_PD)
    obj.fix_variables()
    pd = take_pd(obj)
    ttheta_meas = numpy.array(pd.pd_meas.ttheta, dtype=float)
    signal_meas = numpy.interp(ttheta_meas, ttheta, signal)
    for item, value in zip(pd.pd_meas.items, signal_meas):
        item.intensity = float(value)
        item.intensity_sigma = float(numpy.sqrt(max(value, 1.)))
    cryspy.add_magnetic_correlation_length(obj, 80.)
    cryspy.rhochi_rietveld_refinement(obj)
    row = take_pd(obj).magnetic_correlation.items[0]
    assert row.length == pytest.approx(40., rel=1e-3)
    assert row.length_sigma > 0.


def test_procedures_for_assistants():
    from cryspy.procedure_assistant.methods import find_matching_procedures
    for request in ("Set magnetic domain size", "magnetic correlation length",
                    "magnetic peaks are broader than nuclear"):
        names = [item["name"] for item in find_matching_procedures(request)]
        assert names[0] == "add_magnetic_correlation_length"
    description = cryspy.describe_procedure(cryspy.add_magnetic_correlation_length)
    assert [param["name"] for param in description["parameters"]] == [
        "rcif_object", "length", "phase_label", "anisotropic", "experiment"]
    obj = cryspy.load_file(F_MCIF_PD)
    with pytest.raises(ValueError):
        cryspy.add_magnetic_correlation_length(obj, -1.)
    assert cryspy.report_magnetic_correlation(obj).startswith("No pd experiment")
