"""Asymmetry of powder peaks from the axial divergence (Finger-Cox-Jephcoat).

Item pd_instr_axial_divergence (S/L, H/L): S/L = 0 or H/L = 0 gives the
symmetric pattern; a pattern calculated with known S/L and H/L is fitted
back from other values.
"""
import os

import numpy
import pytest

import cryspy

DIR_TESTS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
F_PBSO4 = os.path.join(DIR_TESTS, "functional", "PbSO4_unpol_powder_test", "main.rcif")


def take_pd(obj):
    return [item for item in obj.items if isinstance(item, cryspy.Pd)][0]


def take_object(s_l: float = None, d_l: float = None, ttheta_max: float = 60.):
    obj = cryspy.file_to_globaln(F_PBSO4)
    pd = take_pd(obj)
    # the low angles, where the axial divergence matters
    pd.range.ttheta_max = 60.
    obj.fix_variables()
    if s_l is not None:
        pd.add_items([cryspy.PdInstrAxialDivergence(
            s_l=s_l, d_l=d_l, ttheta_max=ttheta_max)])
    return obj


def calc_signal(obj):
    cryspy.rhochi_no_refinement(obj)
    pd_proc = take_pd(obj).pd_proc
    return numpy.array(pd_proc.intensity_plus_net, dtype=float) + \
        numpy.array(pd_proc.intensity_minus_net, dtype=float)


def test_zero_axial_divergence_gives_symmetric_pattern():
    signal = calc_signal(take_object())
    numpy.testing.assert_allclose(calc_signal(take_object(0., 0.03)), signal)
    signal_fcj = calc_signal(take_object(0.02, 0.03))
    assert not numpy.allclose(signal_fcj, signal)
    # the integrated intensity is kept
    assert signal_fcj.sum() == pytest.approx(signal.sum(), rel=2e-3)


def test_axial_divergence_is_refined_back():
    obj = take_object(0.02, 0.03)
    signal = calc_signal(obj)
    pd = take_pd(obj)
    pd_meas = pd.pd_meas
    ttheta = numpy.array(pd.pd_proc.ttheta, dtype=float)
    background = numpy.array(pd.pd_proc.intensity_bkg_calc, dtype=float)
    intensity = numpy.interp(numpy.array(pd_meas.ttheta, dtype=float), ttheta, signal + background)
    pd_meas.numpy_intensity = intensity
    pd_meas.numpy_intensity_sigma = numpy.sqrt(numpy.abs(intensity)) + 1.
    pd_meas.numpy_to_items()

    fit = pd.pd_instr_axial_divergence
    fit.s_l, fit.d_l = 0.01, 0.045
    fit.s_l_refinement, fit.d_l_refinement = True, True
    cryspy.rhochi_rietveld_refinement(obj)
    # S/L and H/L enter the model symmetrically: compare the pair
    found = sorted([float(fit.s_l), float(fit.d_l)])
    assert found == pytest.approx([0.02, 0.03], rel=0.03)
