"""Stacking faults in 1D powder (loop stacking_fault, Hendricks-Teller model, diagnose_peak_widths)."""
import os

import numpy
import pytest

import cryspy

DIR = os.path.dirname(__file__)
F_PBSO4 = os.path.join(DIR, "PbSO4_unpol_powder_test", "main.rcif")


def take_block(rcif, class_name: str):
    return [block for block in rcif.items if type(block).__name__ == class_name][0]


def load_pbso4(ttheta_max: float = 70.):
    rcif = cryspy.file_to_globaln(F_PBSO4)
    experiment = take_block(rcif, "Pd")
    experiment.range.ttheta_max = ttheta_max
    rcif.fix_variables()
    return rcif, experiment


def calc_signal(experiment):
    proc = experiment.pd_proc
    return numpy.array(proc.intensity_plus_net) + numpy.array(proc.intensity_minus_net)


def simulate_faulted(probability: float = 0.15, probability_repeat: float = -1.):
    """PbSO4 pattern calculated with faults, written as the measured one."""
    rcif, experiment = load_pbso4()
    cryspy.add_stacking_fault(rcif, normal="0 0 1", shift="0 1/2 0", probability=probability,
                              probability_repeat=probability_repeat)
    rcif.fix_variables()
    cryspy.rhochi_no_refinement(rcif)
    proc = experiment.pd_proc
    y_calc = calc_signal(experiment) + numpy.array(proc.intensity_bkg_calc)
    ttheta = numpy.array(proc.ttheta)
    meas = experiment.pd_meas
    tth_meas = numpy.array(meas.ttheta)
    y_meas = numpy.interp(tth_meas, ttheta, y_calc)
    for item, value in zip(meas.items, y_meas):
        item.intensity = float(value)
        item.intensity_sigma = float(numpy.sqrt(max(value, 1.)))
    return rcif, experiment


def test_zero_probability_is_the_bragg_pattern():
    rcif, experiment = load_pbso4()
    cryspy.rhochi_no_refinement(rcif)
    signal_bragg = calc_signal(experiment)
    cryspy.add_stacking_fault(rcif, normal="0 0 1", shift="0 1/2 0", probability=0.)
    rcif.fix_variables()
    cryspy.rhochi_no_refinement(rcif)
    signal_zero = calc_signal(experiment)
    # sampled on a grid of a fifth of the FWHM
    assert numpy.abs(signal_zero - signal_bragg).max() < 0.01 * signal_bragg.max()


def test_faults_keep_the_integrated_intensity_and_spread_k_odd():
    rcif, experiment = load_pbso4()
    cryspy.rhochi_no_refinement(rcif)
    signal_bragg = calc_signal(experiment)
    cryspy.add_stacking_fault(rcif, normal="0 0 1", shift="0 1/2 0", probability=0.2)
    rcif.fix_variables()
    cryspy.rhochi_no_refinement(rcif)
    signal_fault = calc_signal(experiment)
    # the intensity moves between the reflections, the total is kept
    assert signal_fault.sum() == pytest.approx(signal_bragg.sum(), rel=0.05)
    assert signal_fault.max() < signal_bragg.max() * 1.001
    assert numpy.abs(signal_fault - signal_bragg).max() > 0.05 * signal_bragg.max()
    s_report = cryspy.report_stacking_fault(rcif)
    assert "k odd" in s_report and "0.2000" in s_report
    # rcif round trip
    rcif_2 = cryspy.str_to_globaln(rcif.to_cif())
    row = take_block(rcif_2, "Pd").stacking_fault.items[0]
    assert row.shift_y == pytest.approx(0.5) and row.probability == pytest.approx(0.2)


def test_probability_is_refined_in_the_profile_matching():
    rcif, experiment = simulate_faulted(0.15)
    row = experiment.stacking_fault.items[0]
    row.probability = 0.05
    row.probability_refinement = True
    cryspy.rhochi_profile_matching(rcif)
    assert row.probability == pytest.approx(0.15, abs=0.01)


def test_diagnosis_finds_the_rule_and_the_layers():
    rcif, experiment = simulate_faulted(0.2)
    experiment.items.remove(experiment.stacking_fault)
    s_diagnosis = cryspy.diagnose_peak_widths(rcif)
    assert "SELECTIVE_WIDTH" in s_diagnosis and "k odd" in s_diagnosis
    assert "STACKING_FAULTS" in s_diagnosis and "best layers (" in s_diagnosis
    from cryspy.procedure_assistant.tools import tool_diagnose_fit
    assert "SELECTIVE_WIDTH" in tool_diagnose_fit(rcif)


def test_correlated_faults_are_refined():
    rcif, experiment = simulate_faulted(0.15, 0.5)
    row = experiment.stacking_fault.items[0]
    row.probability, row.probability_repeat = 0.15, 0.25
    row.probability_refinement, row.probability_repeat_refinement = False, True
    cryspy.rhochi_profile_matching(rcif)
    assert row.probability_repeat == pytest.approx(0.5, abs=0.05)
    assert "after a fault" in cryspy.report_stacking_fault(rcif)


def test_wrong_layers_are_refused():
    rcif, _ = load_pbso4()
    with pytest.raises(ValueError, match="plane of the layers"):
        cryspy.add_stacking_fault(rcif, normal="0 0 1", shift="0 1/2 1/2")
    with pytest.raises(ValueError, match="symmetry"):
        cryspy.add_stacking_fault(rcif, normal="0 0 3", shift="0 1/2 0")


F_HO2TI2O7 = os.path.join(DIR, "Ho2Ti2O7_5K_1T_powder_1d_test", "main.rcif")


def calc_polarized_signals(polarization: float, probability: float):
    """I+ and I- of Ho2Ti2O7 (paramagnetic, 1 T) with stacking faults of the layers (0 0 4)."""
    from cryspy.procedure_rhochi.rhochi_by_dictionary import rhochi_calc_chi_sq_by_dictionary
    rcif = cryspy.file_to_globaln(F_HO2TI2O7)
    experiment = take_block(rcif, "Pd")
    experiment.range.ttheta_max = 40.
    experiment.diffrn_radiation.polarization = polarization
    if probability > 0.:
        cryspy.add_stacking_fault(rcif, normal="0 0 4", shift="1/4 1/4 0",
                                  probability=probability, probability_repeat=0.4)
    dict_in_out = {}
    rhochi_calc_chi_sq_by_dictionary(
        rcif.get_dictionary(), dict_in_out=dict_in_out, flag_use_precalculated_data=False,
        flag_calc_analytical_derivatives=False)
    dict_pd = dict_in_out[f"pd_{experiment.data_name.lower():}"]
    return dict_pd["signal_plus"], dict_pd["signal_minus"]


def test_polarized_faults_sum_equals_unpolarized():
    """I+ + I- of a polarized pattern is the unpolarized one, with the same line shapes.

    Fd-3m with layers (0 0 4): some waves do not see the fault shift
    (coherent, |lambda| = 1), they stay Bragg peaks instead of a division
    by zero.
    """
    plus, minus = calc_polarized_signals(0.88, 0.2)
    plus_0, minus_0 = calc_polarized_signals(0., 0.2)
    plus_b, minus_b = calc_polarized_signals(0.88, 0.)
    total = plus + minus
    assert numpy.all(numpy.isfinite(total))
    assert numpy.allclose(total, plus_0 + minus_0, rtol=1e-10, atol=1e-8 * total.max())
    # the faults spread the lines but keep the scale and the polarized difference
    assert total.max() < 3. * (plus_b + minus_b).max()
    assert total.min() > -1e-6 * total.max()
    assert numpy.abs(plus - minus).max() > 0.2 * total.max()
