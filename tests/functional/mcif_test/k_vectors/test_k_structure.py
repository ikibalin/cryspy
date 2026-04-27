import os
import numpy
import cryspy

from cryspy.D_functions_item_loop.function_1_calc_for_magcrystal import calc_f_mag
"""
def test_k_structure():
    dir = os.path.dirname(__file__)
    f_name = os.path.join(dir, "k_structure.rcif")
    rhochi = cryspy.file_to_globaln(f_name)
    d_out = cryspy.rhochi_no_refinement(rhochi)
    chi_sq = d_out["chi_sq"]
    n_points = d_out["n_point"]
    assert chi_sq < 3036.0
    assert int(n_points) == 4901


dir = os.path.dirname(__file__)
f_name = os.path.join(dir, "k_structure.rcif")
rhochi = cryspy.file_to_globaln(f_name)
import numpy as np

def rational_approx_array(x, tol=1e-6, max_den=10000):
    x = np.asarray(x, dtype=float)
    flat = x.ravel()

    q = np.arange(1, max_den + 1)

    p = np.round(flat[:, None] * q[None, :])

    err = np.abs(flat[:, None] - p / q)

    idx = np.argmin(err, axis=1)
    best_err = err[np.arange(len(flat)), idx]

    mask = best_err < tol

    result = []
    for i, ok in enumerate(mask):
        if ok:
            qi = q[idx[i]]
            pi = int(p[i, idx[i]])
            result.append((pi, int(qi)))
        else:
            result.append(None)

    if x.ndim == 0:
        return result[0]

    return np.array(result, dtype=object).reshape(x.shape + (2,))

dict_global = rhochi.get_dictionary()
d_crystal = dict_global["crystal_testmcif"]
d_diffrn = dict_global["diffrn_testexp"]
magnetic_structure_k_vector = d_crystal["magnetic_structure_k_vector"]

k_nd = rational_approx_array(magnetic_structure_k_vector, tol=1e-6, max_den=10000)
d_out = cryspy.rhochi_no_refinement(rhochi)
pass
"""