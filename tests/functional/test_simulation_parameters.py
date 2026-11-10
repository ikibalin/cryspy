"""Simulation procedures take the wavelength and the 2theta range."""
import numpy
import pytest

import cryspy
from cryspy.exceptions import CrysPyError


@pytest.fixture
def crystal_object():
    cryspy_object = cryspy.str_to_globaln("global_")
    cryspy.simulation_paramagnetic_crystal(cryspy_object)
    return cryspy_object


def test_unpolarized_powder_on_given_points(crystal_object):
    cryspy.simulation_unpolarized_neutron_powder_diffraction(
        crystal_object, wavelength=2.4, ttheta_min=10., ttheta_max=70., ttheta_step=0.5)
    pd = crystal_object.pd_unpd
    ttheta = numpy.array(pd.pd_meas.ttheta)
    assert ttheta.size == 121
    assert ttheta[0] == pytest.approx(10.) and ttheta[-1] == pytest.approx(70.)
    assert pd.setup.wavelength == pytest.approx(2.4)
    # Bragg peaks above the flat background
    assert numpy.array(pd.pd_meas.intensity).max() > 2. * pd.pd_background.items[0].intensity
    # the first peak 111 of the cubic cell: 2 asin(lambda sqrt(3) / (2 a))
    a = crystal_object.crystal_ho2Ti2o7.cell.length_a
    ttheta_111 = numpy.degrees(2. * numpy.arcsin(2.4 * numpy.sqrt(3.) / (2. * a)))
    assert numpy.array(pd.pd_peak_ho2ti2o7.ttheta).min() == pytest.approx(ttheta_111, abs=1e-3)


def test_polarized_powder_defaults_kept(crystal_object):
    cryspy.simulation_polarized_neutron_powder_diffraction(crystal_object)
    ttheta = numpy.array(crystal_object.pd_pnpd.pd_meas.ttheta)
    assert ttheta[0] == pytest.approx(4.) and ttheta[-1] == pytest.approx(80.)


def test_single_crystal_reflections_in_range(crystal_object):
    cryspy.simulation_single_crystal(crystal_object, wavelength=0.84, ttheta_min=20., ttheta_max=40.)
    refln = crystal_object.diffrn_sc.diffrn_refln
    index_hkl = numpy.array([refln.index_h, refln.index_k, refln.index_l], dtype=float)
    a = crystal_object.crystal_ho2Ti2o7.cell.length_a
    ttheta = numpy.degrees(2. * numpy.arcsin(0.84 * numpy.linalg.norm(index_hkl, axis=0) / (2. * a)))
    assert ttheta.size > 0 and ttheta.min() >= 20. - 1e-6 and ttheta.max() <= 40. + 1e-6
    assert crystal_object.diffrn_sc.setup.wavelength == pytest.approx(0.84)


def test_wrong_ranges_and_satellites_without_k(crystal_object):
    with pytest.raises(CrysPyError):
        cryspy.simulation_single_crystal(crystal_object, ttheta_min=50., ttheta_max=40.)
    with pytest.raises(CrysPyError):
        cryspy.simulation_unpolarized_neutron_powder_diffraction(crystal_object, wavelength=-1.)
    with pytest.raises(CrysPyError):
        cryspy.simulation_single_crystal(crystal_object, satellites_only=True)
