"""Module realized some mathematical functions.

Functions
---------
    - calc_m_sigma
    - calc_scalar_product_by_vectors
    - calc_scalar_product_by_complex_vectors
    - calc_modulus_sq_by_complex_vector

"""
import numpy
import copy

from numpy.linalg import LinAlgError

def estimate_inversed_hessian_matrix(func, param_0):
    """Estimate inversed Hessian matrix."""
    n_param = len(param_0)
    np_hessian = numpy.zeros(shape=(n_param, n_param), dtype=float)
    np_first_der = numpy.zeros(shape=(n_param,), dtype=float)
    chi_sq = func(param_0)
    perc = 0.01
    
    for i_p_1, p_1 in enumerate(param_0):
        delta_p_1 = perc * numpy.abs(p_1)
        if delta_p_1 < 1e-5:
            delta_p_1 = 1e-5
        param_p = copy.deepcopy(param_0)
        param_m = copy.deepcopy(param_0)
        param_p[i_p_1] += delta_p_1
        param_m[i_p_1] -= delta_p_1
        # chi_sq_p = func(param_p)
        # chi_sq_m = func(param_m)
        # der_second = (chi_sq_p + chi_sq_m - 2.*chi_sq)/(delta_p_1**2)
        # np_hessian[i_p_1, i_p_1] = der_second

        # der_first = (chi_sq_p - chi_sq_m) / (2.*delta_p_1)
        # np_first_der[i_p_1] = der_first

        param_pp = copy.deepcopy(param_0)
        param_pm = copy.deepcopy(param_0)
        param_mp = copy.deepcopy(param_0)
        param_mm = copy.deepcopy(param_0)

        param_pp[i_p_1] += delta_p_1
        param_pm[i_p_1] += delta_p_1
        param_mp[i_p_1] -= delta_p_1
        param_mm[i_p_1] -= delta_p_1
        for i_p_2, p_2 in enumerate(param_0[:i_p_1+1]):
            delta_p_2 = perc * numpy.abs(p_2)
            if delta_p_2 < 1e-5:
                delta_p_2 = 1e-5
            param_pp[i_p_2] += delta_p_2
            param_pm[i_p_2] -= delta_p_2
            param_mp[i_p_2] += delta_p_2
            param_mm[i_p_2] -= delta_p_2
            
            chi_sq_pp = func(param_pp)
            chi_sq_pm = func(param_pm)
            chi_sq_mp = func(param_mp)
            chi_sq_mm = func(param_mm)

            der_second = (chi_sq_pp + chi_sq_mm - chi_sq_pm - chi_sq_mp) / (
                4. * delta_p_1 * delta_p_2) 
            np_hessian[i_p_1, i_p_2] = der_second
            np_hessian[i_p_2, i_p_1] = der_second

            param_pp[i_p_2] -= delta_p_2
            param_pm[i_p_2] += delta_p_2
            param_mp[i_p_2] -= delta_p_2
            param_mm[i_p_2] += delta_p_2

    try:
        np_hessian_inv = numpy.linalg.inv(np_hessian)
    except LinAlgError:
        np_hessian_inv = None
    func(param_0)
    return np_hessian_inv, np_first_der

