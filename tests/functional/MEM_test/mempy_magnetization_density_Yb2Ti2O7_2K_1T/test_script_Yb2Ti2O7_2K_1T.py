import os

import cryspy
import pytest

# MEM of these tests needs several GB: one worker runs them one after the other
# (pytest -n 2 --dist loadgroup), two at once exhaust the memory of the CI runner
pytestmark = pytest.mark.xdist_group("mem_large")

def test_mem_magnetization_density_for_chi():
    f_name = os.path.join(os.path.dirname(__file__), "mempy_magnetization_density_Yb2Ti2O7_2K_1T.rcif")
    rcif_obj = cryspy.load_file(f_name)
    # cryspy.rhochi_rietveld_refinement(rcif_obj)
    res = cryspy.mempy_magnetization_density_reconstruction(rcif_obj)
    assert res["chi_sq"] < 7.5
    # cryspy.mempy_spin_density_reconstruction(rcif_obj)
