"""Powder integrated intensities are reused between refinement steps.

The cache key was misspelt ("iint_minu"), so iint_plus/iint_minus of a
paramagnetic phase were recalculated at every step. Now a step that changes
nothing the intensities depend on takes them from dict_in_out, and a step
that changes such a parameter recalculates them.
"""
import os

import numpy
import pytest

import cryspy
from cryspy.procedure_rhochi import rhochi_pd, rhochi_pd2d, rhochi_tof

DIR = os.path.dirname(os.path.abspath(__file__))
CASES = [
    (os.path.join(DIR, "Ho2Ti2O7_5K_1T_powder_1d_test", "main.rcif"),
     "pd_", rhochi_pd, "calc_powder_iint_1d_para",
     rhochi_pd.calc_chi_sq_for_pd_by_dictionary),
    (os.path.join(DIR, "Fe3O4_150K_6T_powder_2d_test", "Fe3O4_150K_6.0T.rcif"),
     "pd2d_", rhochi_pd2d, "calc_powder_iint_2d_para",
     rhochi_pd2d.calc_chi_sq_for_pd2d_by_dictionary),
    (os.path.join(DIR, "test_tof", "model_pV.rcif"),
     "tof_", rhochi_tof, "calc_powder_iint_1d_para",
     rhochi_tof.calc_chi_sq_for_tof_by_dictionary),
]


@pytest.mark.parametrize("f_name, prefix, module, iint_name, calc_chi_sq", CASES)
def test_iint_taken_from_cache(monkeypatch, f_name, prefix, module, iint_name, calc_chi_sq):
    obj = cryspy.file_to_globaln(f_name)
    obj.fix_variables()
    global_dict = obj.get_dictionary()
    dict_exp = [value for key, value in global_dict.items() if key.startswith(prefix)][0]
    dict_crystals = [value for key, value in global_dict.items()
                     if key.startswith("crystal_")]

    n_call = []
    iint_function = getattr(module, iint_name)

    def counted(*args, **kwargs):
        n_call.append(1)
        return iint_function(*args, **kwargs)

    monkeypatch.setattr(module, iint_name, counted)
    dict_in_out = {}
    chi_sq_1 = calc_chi_sq(dict_exp, dict_crystals, dict_in_out=dict_in_out)[0]
    n_first = len(n_call)
    assert n_first > 0
    chi_sq_2 = calc_chi_sq(dict_exp, dict_crystals, dict_in_out=dict_in_out,
                           flag_use_precalculated_data=True)[0]
    assert len(n_call) == n_first
    assert chi_sq_2 == pytest.approx(chi_sq_1, rel=1e-12)

    # a refined structure parameter changes the intensities: no cache
    dict_crystals[0]["flags_atom_b_iso"][...] = False
    dict_crystals[0]["flags_atom_b_iso"].flat[0] = True
    dict_crystals[0]["atom_b_iso"] = numpy.array(dict_crystals[0]["atom_b_iso"]) + 0.3
    chi_sq_3 = calc_chi_sq(dict_exp, dict_crystals, dict_in_out=dict_in_out,
                           flag_use_precalculated_data=True)[0]
    assert len(n_call) > n_first
    assert chi_sq_3 != pytest.approx(chi_sq_1, rel=1e-9)
