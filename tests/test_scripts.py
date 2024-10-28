# %%
import pytest
import os
import cryspy

L_RIETVELD_REFINEMENT_FILE = [
    os.path.join("rhochi_polarized_powder_1d_Ho2Ti2O7_5K_1T", "rhochi_polarized_powder_1d_Ho2Ti2O7_5K_1T.rcif"),
    os.path.join("rhochi_polarized_powder_2d_DyAl_5K_5T", "rhochi_polarized_powder_2d_DyAl_5K_5T.rcif"),
    os.path.join("rhochi_single_susceptibility_Ho2Ti2O7", "rhochi_single_susceptibility_Ho2Ti2O7.rcif"),
    os.path.join("rhochi_unpolarized_powder_1d_mcif_URu0.96Rh0.04Si2", "rhochi_unpolarized_powder_1d_mcif_URu0.96Rh0.04Si2.rcif"),
    os.path.join("rhochi_unpolarized_powder_1d_PbSO4", "rhochi_unpolarized_powder_1d_PbSO4.rcif"),
    os.path.join("rhochi_unpolarized_powder_1d_xrays_PbSO4", "rhochi_unpolarized_powder_1d_xrays_PbSO4.rcif"),
    os.path.join("rhochi_unpolarized_tof_powder_CeCuAl", "rhochi_unpolarized_tof_powder_CeCuAl.rcif"),
    os.path.join("only_for_test", "CuDSON.rcif"),
    os.path.join("only_for_test", "emso.rcif"),
    os.path.join("only_for_test", "Mg2Sb3Yb3O14.rcif"),
    os.path.join("only_for_test", "Mg2Sb3Ho3O14.rcif"),
    os.path.join("only_for_test", "mnbi4te7_single_lambdaover2.rcif"),
]

L_RIETVELD_REFINEMENT_CHI_SQ = [
    3273.82,
    67620.70,
    1135.42,
    2171.20,
    21016.14,
    6493.18,
    10716.91,
    104444.83,
     721829, 
     101.68,
     571.93,
     410.51,
]


   
def test_atom_multiplicity():
    f_dir_examples = os.path.join(os.path.dirname(os.path.dirname(__file__)), "examples")
    f_name = os.path.join(
        f_dir_examples,
        os.path.join("only_for_test",
                     "TbCo2Ni3.rcif"))
    rcif_obj = cryspy.load_file(f_name)
    cryspy.rhochi_no_refinement(rcif_obj)
    l_mult = rcif_obj.crystal_phase1.atom_site.multiplicity
    l_res = [1, 3, 3, 2, 2]
    for mult, res in zip(l_mult, l_res):
        assert mult == res
    


def test_rietveld_refinement():
    f_dir_examples = os.path.join(os.path.dirname(os.path.dirname(__file__)), "examples")
    # %%
    for f_file, chi_sq in zip(L_RIETVELD_REFINEMENT_FILE, L_RIETVELD_REFINEMENT_CHI_SQ):
        f_name = os.path.join(f_dir_examples, f_file)

        rcif_obj = cryspy.load_file(f_name)

        res = cryspy.rhochi_rietveld_refinement(rcif_obj)
        
        assert res["chi_sq"] < chi_sq



def test_mem():
    f_dir_examples = os.path.join(os.path.dirname(os.path.dirname(__file__)), "examples")
    f_name = os.path.join(
        f_dir_examples,
        os.path.join("mempy_spin_density_YTiO3",
                     "mempy_spin_density_YTiO3.rcif"))
    try:
        rcif_obj = cryspy.load_file(f_name)
        cryspy.mempy_spin_density_reconstruction(rcif_obj)
        assert True
    except:
        assert False

    f_name = os.path.join(
        f_dir_examples,
        os.path.join("mempy_magnetization_density_Yb2Ti2O7_2K_1T",
                     "mempy_magnetization_density_Yb2Ti2O7_2K_1T.rcif"))
    try:
        rcif_obj = cryspy.load_file(f_name)
        cryspy.mempy_magnetization_density_reconstruction(rcif_obj)
        assert True
    except:
        assert False
