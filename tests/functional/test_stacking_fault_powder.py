"""Random stacking faults in 1D powder (loop stacking_fault, diagnose_peak_widths)."""
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


def simulate_faulted(probability: float = 0.15):
    """PbSO4 pattern calculated with faults, written as the measured one."""
    rcif, experiment = load_pbso4()
    cryspy.add_stacking_fault(rcif, normal="0 0 1", shift="0 1/2 0", probability=probability)
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
    assert signal_fault.sum() == pytest.approx(signal_bragg.sum(), rel=0.03)
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
