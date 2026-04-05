import os 

import cryspy

def test_mem_magnetization_density_for_momemts():
    f_name = os.path.join(os.path.dirname(__file__), "main.rcif")
    rcif_obj = cryspy.load_file(f_name)
    # cryspy.rhochi_rietveld_refinement(rcif_obj)
    res =cryspy.mempy_magnetization_density_reconstruction(rcif_obj)
    # cryspy.mempy_spin_density_reconstruction(rcif_obj)
    assert res["chi_sq"] < 2.76


def test_mem_magnetization_density_for_momemts_unpol_SC():
    f_name = os.path.join(os.path.dirname(__file__), "mnf2_unpol_SC.rcif")
    rcif_obj = cryspy.load_file(f_name)
    res = cryspy.rhochi_rietveld_refinement(rcif_obj)
    assert res["chi_sq"] < 80600
    res = cryspy.mempy_magnetization_density_reconstruction(rcif_obj)
    assert res["chi_sq"] < 80


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
    f_name = os.path.join(os.path.dirname(__file__), "MnF2_multipole_prior.rcif")
    rcif_obj = cryspy.load_file(f_name)
    dict_crystal = rcif_obj.crystal_mnf2.get_dictionary()
    multipol_density = rcif_obj.crystal_mnf2.save_atom_rho_multipole_density_to_den("")
    # cryspy.rhochi_rietveld_refinement(rcif_obj)
    res =cryspy.mempy_magnetization_density_reconstruction(rcif_obj)
    # cryspy.mempy_spin_density_reconstruction(rcif_obj)
    assert res["chi_sq"] < 2.78
