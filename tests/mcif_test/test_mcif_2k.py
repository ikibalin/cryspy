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
