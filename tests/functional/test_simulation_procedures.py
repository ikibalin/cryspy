"""Every simulation procedure runs on the paramagnetic and on the magnetic
example crystal, also when parameters of the unit cell are marked as refined
(the derivatives by the unit cell were taken then, and the flipping ratio of
a paramagnetic crystal failed with KeyError 'eq_css')."""
import numpy
import pytest

import cryspy

CRYSTALS = {
    "paramagnetic": cryspy.simulation_paramagnetic_crystal,
    "magnetic": cryspy.simulation_magnetic_crystal,
}


def take_crystal_object(crystal_kind: str, refined_cell: bool):
    cryspy_object = cryspy.str_to_globaln("global_")
    CRYSTALS[crystal_kind](cryspy_object)
    if refined_cell:
        crystal = [item for item in cryspy_object.items if isinstance(item, cryspy.Crystal)][0]
        crystal.cell.length_a_refinement = True
    return cryspy_object


def take_experiments(cryspy_object, experiment_class):
    return [item for item in cryspy_object.items if isinstance(item, experiment_class)]


def check_flip_ratio(cryspy_object):
    diffrn = take_experiments(cryspy_object, cryspy.Diffrn)[0]
    fr = numpy.array(diffrn.diffrn_refln.fr, dtype=float)
    fr_calc = numpy.array(diffrn.diffrn_refln.fr_calc, dtype=float)
    assert fr.size > 0 and numpy.all(numpy.isfinite(fr_calc))
    # the simulated data are the calculated ones, and the field polarizes the moments
    assert numpy.allclose(fr, fr_calc, atol=1e-5)
    assert numpy.abs(fr - 1.).max() > 1e-3


def check_single_crystal(cryspy_object):
    diffrn = take_experiments(cryspy_object, cryspy.Diffrn)[0]
    intensity = numpy.array(diffrn.diffrn_refln.intensity, dtype=float)
    assert intensity.size > 0 and numpy.all(numpy.isfinite(intensity))
    assert intensity.max() > 0.


def check_powder(cryspy_object):
    pd = take_experiments(cryspy_object, cryspy.Pd)[0]
    intensity = numpy.array(pd.pd_proc.intensity_plus_net, dtype=float) + \
        numpy.array(pd.pd_proc.intensity_minus_net, dtype=float)
    assert intensity.size > 0 and numpy.all(numpy.isfinite(intensity))
    assert intensity.max() > 0.


SIMULATIONS = {
    "simulation_fliping_ratio": check_flip_ratio,
    "simulation_single_crystal": check_single_crystal,
    "simulation_unpolarized_neutron_powder_diffraction": check_powder,
    "simulation_polarized_neutron_powder_diffraction": check_powder,
}


@pytest.mark.parametrize("refined_cell", [False, True])
@pytest.mark.parametrize("crystal_kind", list(CRYSTALS))
@pytest.mark.parametrize("simulation", list(SIMULATIONS))
def test_simulation(simulation, crystal_kind, refined_cell):
    cryspy_object = take_crystal_object(crystal_kind, refined_cell)
    getattr(cryspy, simulation)(cryspy_object)
    SIMULATIONS[simulation](cryspy_object)


@pytest.mark.parametrize("refined_cell", [False, True])
@pytest.mark.parametrize("crystal_kind", list(CRYSTALS))
def test_simulation_current_experiments(crystal_kind, refined_cell):
    cryspy_object = take_crystal_object(crystal_kind, refined_cell)
    cryspy.simulation_fliping_ratio(cryspy_object)
    diffrn = take_experiments(cryspy_object, cryspy.Diffrn)[0]
    fr_before = numpy.array(diffrn.diffrn_refln.fr, dtype=float)
    for item in diffrn.diffrn_refln.items:
        item.fr = 1.
    cryspy.simulation_current_experiments(cryspy_object)
    assert numpy.allclose(numpy.array(diffrn.diffrn_refln.fr, dtype=float), fr_before, atol=1e-5)
