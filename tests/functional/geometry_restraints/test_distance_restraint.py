"""A distance restraint pulls the refined distance S-O1 of PbSO4 to its target."""
import os

import pytest

import cryspy

F_PBSO4 = os.path.join(
    os.path.dirname(__file__), os.pardir, os.pardir, os.pardir, "examples",
    "rhochi_unpolarized_powder_1d_PbSO4", "rhochi_unpolarized_powder_1d_PbSO4.rcif")


def take_refined_crystal(target: str = ""):
    rcif = cryspy.load_file(F_PBSO4)
    crystal = [item for item in rcif.items if type(item).__name__ == "Crystal"][0]
    for item in crystal.atom_site.items:
        if item.label in ("S", "O1"):
            item.fract_x_refinement = True
            item.fract_z_refinement = True
    if target != "":
        crystal.add_restraint("S O1", target=target, sigma=0.002)
    cryspy.rhochi_rietveld_refinement(rcif)
    return crystal


def test_distance_restraint_pulls_to_target():
    distance_free = take_refined_crystal().calc_distance_or_angle("S O1")["value"]
    crystal = take_refined_crystal("1.40")
    distance = crystal.calc_distance_or_angle("S O1")["value"]
    assert abs(distance - 1.40) < abs(distance_free - 1.40) - 0.01
    # diff is updated at the end of the refinement
    assert crystal.restr_distance.items[0].diff == pytest.approx(distance - 1.40, abs=1e-6)
