"""Extra broadening of the magnetic reflections of powder phases (loop magnetic_broadening).

The magnetic and interference intensities of a phase get an extra
Lorentzian (half width `width` in Q); the nuclear intensity keeps its
profile. width = 0 and old files (without the loop) give the same patterns
as before.
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
def test_zero_width_and_tensor_of_one_width(f_name):
    """Width 0 (refined, so the extra profile is calculated) changes
    nothing; the tensor made from one width gives the pattern of this
    width (k = 0, polarized paramagnetic, k != 0 satellites)."""
    ttheta, plus_0, minus_0 = calc_pattern(cryspy.load_file(f_name))
    for anisotropic in (False, True):
        obj = cryspy.load_file(f_name)
        cryspy.add_magnetic_broadening(obj, 0., anisotropic=anisotropic)
        plus_w, minus_w = calc_pattern(obj)[1:]
        # pd_proc keeps 5 decimals
        assert numpy.allclose(plus_w, plus_0, rtol=1e-10, atol=2e-5)
        assert numpy.allclose(minus_w, minus_0, rtol=1e-10, atol=2e-5)

    obj = cryspy.load_file(f_name)
    cryspy.add_magnetic_broadening(obj, 0.03)
    plus_3, minus_3 = calc_pattern(obj)[1:]
    assert numpy.abs(plus_3 + minus_3 - plus_0 - minus_0).max() > \
        1e-3 * (plus_0 + minus_0).max()
    obj = cryspy.load_file(f_name)
    cryspy.add_magnetic_broadening(obj, 0.03, anisotropic=True)
    plus_a, minus_a = calc_pattern(obj)[1:]
    assert numpy.allclose(plus_a, plus_3, rtol=1e-8, atol=2e-5)
    assert numpy.allclose(minus_a, minus_3, rtol=1e-8, atol=2e-5)


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_satellites_keep_their_integrated_intensity():
    """Ho2BaNiO5, k = (1/2, 0, 1/2): purely magnetic satellites get broader
    and lower, the integrated intensity of a satellite stays."""
    ttheta, plus_0, minus_0 = calc_pattern(cryspy.load_file(F_HO2BANIO5))
    obj = cryspy.load_file(F_HO2BANIO5)
    cryspy.add_magnetic_broadening(obj, 0.025)
    plus_w, minus_w = calc_pattern(obj)[1:]
    signal_0, signal_w = plus_0 + minus_0, plus_w + minus_w
    assert signal_w.max() < signal_0.max()
    # first satellite, at about 20.4 degrees (lambda 2.524 A)
    peak = numpy.argmax(numpy.where(ttheta < 25., signal_0, 0.))
    window = numpy.abs(ttheta - ttheta[peak]) < 3.
    ratio = signal_w[window].sum() / signal_0[window].sum()
    # the Lorentzian tails outside +-3 degrees and of the neighbours: 10 %
    assert 0.85 < ratio < 1.1
    assert signal_w[peak] < 0.8 * signal_0[peak]


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_negative_width_makes_magnetic_peaks_narrower():
    """pd_test has the Lorentzian term Y = 0.030 deg of the instrument: a
    negative width takes it away from the magnetic intensity only."""
    ttheta, plus_0, minus_0 = calc_pattern(cryspy.load_file(F_MCIF_PD))
    obj = cryspy.load_file(F_MCIF_PD)
    cryspy.add_magnetic_broadening(obj, -0.005)
    plus_n, minus_n = calc_pattern(obj)[1:]
    signal_0, signal_n = plus_0 + minus_0, plus_n + minus_n
    assert signal_n.max() >= signal_0.max() - 1e-4
    assert numpy.abs(signal_n - signal_0).max() > 1e-4 * signal_0.max()
    report = cryspy.report_magnetic_broadening(obj)
    assert "longer than the nuclear ones" in report
    assert "Y = 0.0301 deg" in report


def test_round_trip_of_the_loop():
    obj = cryspy.load_file(F_MCIF_PD)
    cryspy.add_magnetic_broadening(obj, 0.01, anisotropic=True)
    pd = take_pd(obj)
    pd.magnetic_broadening.items[0].width_33 = 0.04
    text = obj.to_cif()
    assert "_magnetic_broadening_width_33" in text
    obj_2 = cryspy.str_to_globaln(text)
    row = take_pd(obj_2).magnetic_broadening.items[0]
    assert row.width_33 == pytest.approx(0.04)
    assert row.width_11_refinement and not row.width_12_refinement
    report = cryspy.report_magnetic_broadening(obj_2)
    assert "c* 0.04000" in report


def test_refinement_recovers_the_width():
    """Pattern calculated with width 0.025 1/A, refined from 0."""
    obj = cryspy.load_file(F_MCIF_PD)
    obj.fix_variables()
    cryspy.add_magnetic_broadening(obj, 0.025)
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
    cryspy.add_magnetic_broadening(obj, 0.)
    cryspy.rhochi_rietveld_refinement(obj)
    row = take_pd(obj).magnetic_broadening.items[0]
    assert row.width == pytest.approx(0.025, rel=2e-3)
    assert row.width_sigma > 0.


def test_procedures_for_assistants():
    from cryspy.procedure_assistant.methods import find_matching_procedures
    for request in ("Set magnetic domain size", "magnetic broadening",
                    "magnetic peaks are broader than nuclear"):
        names = [item["name"] for item in find_matching_procedures(request)]
        assert names[0] == "add_magnetic_broadening", request
    description = cryspy.describe_procedure(cryspy.add_magnetic_broadening)
    assert [param["name"] for param in description["parameters"]] == [
        "rcif_object", "width", "phase_label", "anisotropic", "experiment"]
    obj = cryspy.load_file(F_MCIF_PD)
    assert cryspy.report_magnetic_broadening(obj).startswith("No pd experiment")
