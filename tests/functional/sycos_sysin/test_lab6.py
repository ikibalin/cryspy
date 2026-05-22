import os
import numpy
import cryspy


def test_lab6():
    dir = os.path.dirname(__file__)
    f_name = os.path.join(dir, "LaB6.rcif")
    rhochi = cryspy.file_to_globaln(f_name)
    d_out = cryspy.rhochi_rietveld_refinement(rhochi)
    assert d_out["chi_sq"] <= 13681849
