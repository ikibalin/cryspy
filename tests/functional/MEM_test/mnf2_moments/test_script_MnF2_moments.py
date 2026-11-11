import os
import cryspy
import pytest

# MEM of these tests needs several GB: one worker runs them one after the other
# (pytest -n 2 --dist loadgroup), two at once exhaust the memory of the CI runner
pytestmark = pytest.mark.xdist_group("mem_large")


def test_mem_magnetization_density_for_momemts():
    f_name = os.path.join(os.path.dirname(__file__), "main.rcif")
    rcif_obj = cryspy.load_file(f_name)
    # cryspy.rhochi_rietveld_refinement(rcif_obj)
    res = cryspy.mempy_magnetization_density_reconstruction(rcif_obj)
    # cryspy.mempy_spin_density_reconstruction(rcif_obj)
    assert res["chi_sq"] < 2.82


def test_mem_magnetization_density_for_momemts_unpol_SC():
    f_name = os.path.join(os.path.dirname(__file__), "mnf2_unpol_SC.rcif")
    rcif_obj = cryspy.load_file(f_name)
    res = cryspy.rhochi_rietveld_refinement(rcif_obj)
    assert res["chi_sq"] < 80600
    res = cryspy.mempy_magnetization_density_reconstruction(rcif_obj)
    assert res["chi_sq"] < 80


def test_mem_nuclear_density_for_momemts_unpol_SC():
    f_name = os.path.join(os.path.dirname(__file__), "mnF2_80K.rcif")
    rcif_obj = cryspy.load_file(f_name)
    res = cryspy.mempy_nuclear_density_reconstruction(rcif_obj)
    assert res["chi_sq"] < 622


def test_unpol_SC_to_html():
    f_name = os.path.join(os.path.dirname(__file__), "mnf2_unpol_SC.rcif")
    rcif_obj = cryspy.load_file(f_name)
    try:
        rcif_obj.diffrn_exp1.report_html()
    except Exception as e:
        print(f"Error occurred: {e}")
        assert False
    assert True


def test_mem_magnetization_density_for_momemts_multipole_prior():
    f_name = os.path.join(
        os.path.dirname(__file__), "MnF2_multipole_prior.rcif"
    )
    rcif_obj = cryspy.load_file(f_name)
    rcif_obj.crystal_mnf2.get_dictionary()
    rcif_obj.crystal_mnf2.save_atom_rho_multipole_density_to_den("")
    # cryspy.rhochi_rietveld_refinement(rcif_obj)
    res = cryspy.mempy_magnetization_density_reconstruction(rcif_obj)
    # cryspy.mempy_spin_density_reconstruction(rcif_obj)
    assert res["chi_sq"] < 2.83


def test_mem_nuclear_density_with_absorption():
    """The transmission of the crystal multiplies the intensities of MEM:
    a strong absorption without a refined scale spoils the fit."""
    f_name = os.path.join(os.path.dirname(__file__), "mnF2_80K.rcif")
    rcif_obj = cryspy.load_file(f_name)
    diffrn = [item for item in rcif_obj.items if isinstance(item, cryspy.Diffrn)][0]
    diffrn.add_items([cryspy.ExptlAbsorpt(correction_type="sphere", coefficient_mu=3.),
                      cryspy.ExptlCrystal(size_rad=1.)])
    res = cryspy.mempy_nuclear_density_reconstruction(rcif_obj)
    assert res["chi_sq"] > 2 * 622
