import os
import numpy
import cryspy



def test_refine_mcif_2k():
    dir = os.path.dirname(__file__)
    f_name = os.path.join(dir, "2k.rcif")
    rhochi = cryspy.file_to_globaln(f_name)
    d_out = cryspy.rhochi_no_refinement(rhochi)
    cell = rhochi.crystal_s2k.cell
    assert cell.type_cell.startswith("o")


def test_refine_mcif_1k():
    dir = os.path.dirname(__file__)
    f_name = os.path.join(dir, "1k.rcif")
    rhochi = cryspy.file_to_globaln(f_name)
    d_out = cryspy.rhochi_no_refinement(rhochi)
    crystal = rhochi.crystal_s1k
    assert numpy.isclose(crystal.atom_site_moment["Gd_2"].crystalaxis_y, 0.) 
    d_out = cryspy.rhochi_rietveld_refinement(rhochi)
    assert numpy.isclose(crystal.atom_site_moment["Gd_2"].crystalaxis_x, -5.867102366251087) 
    assert d_out["chi_sq"] <=  57741.
