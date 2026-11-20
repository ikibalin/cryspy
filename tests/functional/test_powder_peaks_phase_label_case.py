"""Peak (pd_peak, tof_peak) and refln loops after a calculation, phase labels in any case.

The refinement keys the calculated data of a phase by its label in lower
case; the loops of the output must be found for a label such as "Fe3O4" too.
"""
import os

import pytest

import cryspy

DIR = os.path.dirname(__file__)
CASES = (
    (os.path.join(DIR, "PbSO4_unpol_powder_test", "main.rcif"), "Pd", "PdPeakL"),
    (os.path.join(DIR, "test_tof", "model_pV.rcif"), "TOF", "TOFPeakL"),
    (os.path.join(DIR, "Fe3O4_150K_6T_powder_2d_test", "Fe3O4_150K_6.0T.rcif"), "Pd2d", "PdPeakL"),
)


@pytest.mark.parametrize("f_name, class_name, peak_class_name", CASES)
def test_peak_loop_with_mixed_case_phase_label(f_name, class_name, peak_class_name):
    rcif_obj = cryspy.load_file(f_name)
    experiment = [item for item in rcif_obj.items if type(item).__name__ == class_name][0]
    experiment.items[:] = [
        item for item in experiment.items
        if type(item).__name__ not in (peak_class_name, "ReflnL")]
    for item_phase in experiment.phase.items:
        item_phase.label = item_phase.label.capitalize()  # "pbso4" -> "Pbso4"

    cryspy.rhochi_no_refinement(rcif_obj)

    peaks = [item for item in experiment.items if type(item).__name__ == peak_class_name]
    assert len(peaks) == len(experiment.phase.items)
    assert len(peaks[0].items) > 0
    assert peaks[0].loop_name == experiment.phase.items[0].label
