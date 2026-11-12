"""Significance of the parameters added to a least-squares model.

Hamilton's R-factor ratio test (W. C. Hamilton, Acta Cryst. 18 (1965)
502) asks whether a decrease of the weighted sum of squares chi_sq that
follows the refinement of more parameters is larger than the noise of the
data gives by chance. For nested models (the larger model keeps all
parameters of the smaller one) it is the F-test

    F = ((chi_sq_0 - chi_sq_1) / b) / (chi_sq_1 / (n - m)),

where chi_sq_0 and chi_sq_1 are the sums of the smaller and the larger
model, b the number of the added parameters, m the number of the
parameters of the larger model and n the number of points; F follows the
F distribution with (b, n - m) degrees of freedom when the added
parameters are zero in truth. Dividing by chi_sq_1 / (n - m) makes the
test independent of the scale of the sigmas of the data (goodness of fit
different from 1).
"""
import math

from scipy.stats import f as f_distribution


def calc_f_test_added_parameters(
        chi_sq_before: float, chi_sq_after: float, n_point: int, n_parameter: int,
        n_added: int = 1) -> tuple:
    """F statistic and p-value of the parameters added to a model.

    Parameters
    ----------
    chi_sq_before : float
        Weighted sum of squares of the model without the added parameters.
    chi_sq_after : float
        Weighted sum of squares of the model with them (refined).
    n_point : int
        Number of the measured points (reflections, flip ratios, profile points).
    n_parameter : int
        Number of the refined parameters of the model with the added ones.
    n_added : int
        Number of the added parameters.

    Returns
    -------
    tuple
        (F, p_value): the p-value is the probability to get such a decrease
        of chi_sq by chance; (0, 1) when chi_sq does not decrease, (nan,
        nan) when there are no degrees of freedom left (n_point <= n_parameter).
    """
    n_freedom = int(n_point) - int(n_parameter)
    if n_freedom <= 0 or n_added <= 0 or not (math.isfinite(chi_sq_before) and math.isfinite(chi_sq_after)):
        return math.nan, math.nan
    if chi_sq_after >= chi_sq_before:
        return 0., 1.
    if chi_sq_after <= 0.:
        return math.inf, 0.
    f_value = ((chi_sq_before - chi_sq_after) / n_added) / (chi_sq_after / n_freedom)
    return f_value, float(f_distribution.sf(f_value, n_added, n_freedom))
