"""Ikeda-Carpenter moderator pulse convoluted with a pseudo-Voigt (TOF powder).

The moderator emits neutrons of one wavelength with the time distribution
of Ikeda and Carpenter (Nucl. Instrum. Methods A 239, 536 (1985)),
with tau >= 0 the delay in microseconds:

    f(tau) = alpha/2 * [(1-R) (alpha tau)**2 exp(-alpha tau)
             + 2 R alpha**2 beta / (alpha-beta)**3
               * (exp(-beta tau) - exp(-alpha tau)
                  * (1 + (alpha-beta) tau + (alpha-beta)**2 tau**2 / 2))].

Its area is 1. The peak is f convoluted with the Thompson-Cox-Hastings
pseudo-Voigt of the sample and the instrument (FullProf Npr=13, Mantid
IkedaCarpenterPV). Every term of f is tau**n exp(-a tau) with n = 0, 1, 2,
and the convolution of tau**n exp(-a tau) is (-d/da)**n of the convolution
of exp(-a tau), which is known in closed form for a Gaussian (erfc) and a
Lorentzian (exponential integral E1). So the profile is exact, without the
finite-difference approximation of tau**2 exp(-alpha tau) used by Mantid.

Parametrisation as in Mantid IkedaCarpenterPV, with the wavelength
lambda = 2 d sin(theta_bank):

    alpha = 1 / (alpha0 + alpha1 * lambda),
    beta = 1 / beta0,
    R = exp(-81.799 / (kappa * lambda**2)).
"""
import numpy
import scipy.special

na = numpy.newaxis

# E (meV) = 81.799 / lambda**2 (Angstrom) for a neutron
NEUTRON_ENERGY_MEV_ANGSTROM_SQ = 81.799

# |Re z| above which exp(z) E1(z) is taken from its asymptotic series:
# exp(z) or E1(z) would overflow, and the series is exact to 1e-12 there.
_Z_ASYMPTOTIC = 600.


def calc_ikeda_carpenter_alpha_beta_r(
        wavelength, alpha0: float, alpha1: float, beta0: float, kappa: float):
    """Decay rates and the switching function of the Ikeda-Carpenter pulse.

    Parameters
    ----------
    wavelength : numpy.ndarray
        Neutron wavelength lambda, Angstrom.
    alpha0 : float
        Constant part of 1/alpha, microseconds.
    alpha1 : float
        Part of 1/alpha linear in lambda, microseconds per Angstrom.
    beta0 : float
        1/beta, decay time of the slow (storage) term, microseconds.
    kappa : float
        Moderator temperature parameter, meV; kappa <= 0 gives R = 0
        (no slow term).

    Returns
    -------
    alpha, beta, r : numpy.ndarray
        alpha and beta in 1/microseconds, r (dimensionless) of the shape
        of ``wavelength``.

    Reference: Mantid IkedaCarpenterPV fit function.
    """
    wavelength = numpy.asarray(wavelength, dtype=float)
    alpha = 1. / (alpha0 + alpha1 * wavelength)
    beta = numpy.full_like(wavelength, 1. / beta0)
    if kappa > 0.:
        r = numpy.exp(-NEUTRON_ENERGY_MEV_ANGSTROM_SQ / (
            kappa * numpy.square(wavelength)))
    else:
        r = numpy.zeros_like(wavelength)
    return alpha, beta, r


def _calc_gauss_exponential_terms(a, sigma, delta):
    """Convolution of exp(-a tau) (tau >= 0) with a unit-area Gaussian.

    Returns ``(F, -dF/da, d2F/da2)``: the convolutions of exp(-a tau),
    tau exp(-a tau) and tau**2 exp(-a tau). ``a`` and ``sigma`` are columns
    broadcast against ``delta`` = t - t_hkl (microseconds).
    """
    sqrt2 = numpy.sqrt(2.)
    big_a = a * numpy.square(sigma) - delta  # a sigma**2 - delta
    y = big_a / (sigma * sqrt2)
    # exp(u - y**2) = exp(-delta**2 / (2 sigma**2)): no overflow
    gauss = numpy.exp(-0.5 * numpy.square(delta / sigma))
    # exp(u) erfc(y), with u = a (a sigma**2 - 2 delta) / 2
    with numpy.errstate(over="ignore", invalid="ignore"):
        u = 0.5 * a * (a * numpy.square(sigma) - 2. * delta)
        e_u_erfc = numpy.where(
            y >= 0., gauss * scipy.special.erfcx(numpy.abs(y)),
            numpy.exp(numpy.minimum(u, 0.)) * scipy.special.erfc(y))
    c = sigma * numpy.sqrt(2. / numpy.pi)
    f_0 = 0.5 * e_u_erfc
    f_1 = -0.5 * (big_a * e_u_erfc - c * gauss)
    f_2 = 0.5 * ((numpy.square(big_a) + numpy.square(sigma)) * e_u_erfc
                 - c * big_a * gauss)
    return f_0, f_1, f_2


def _calc_exp_e1(z):
    """exp(z) E1(z) for complex z, safe for large |Re z|."""
    far = numpy.abs(z.real) > _Z_ASYMPTOTIC
    z_near = numpy.where(far, 1., z)
    with numpy.errstate(over="ignore", invalid="ignore"):
        phi = numpy.exp(z_near) * scipy.special.exp1(z_near)
    inv_z = 1. / numpy.where(far, z, 1.)
    phi_far = inv_z * (1. - inv_z * (1. - inv_z * (2. - 6. * inv_z)))
    return numpy.where(far, phi_far, phi), far, inv_z


def _calc_lorentz_exponential_terms(a, fwhm, delta):
    """Convolution of exp(-a tau) (tau >= 0) with a unit-area Lorentzian.

    Same return as ``_calc_gauss_exponential_terms``; ``fwhm`` is the full
    width at half maximum of the Lorentzian (microseconds).

    With z = -a delta + i a fwhm / 2 and phi = exp(z) E1(z), F = -Im(phi)/pi;
    dphi/dz = phi - 1/z and dz/da = z/a give the derivatives over a.
    """
    z = -a * delta + 0.5j * a * fwhm
    phi, far, inv_z = _calc_exp_e1(z)
    # z phi - 1 and z**2 phi - z + 1 lose digits for large |z|: use series
    z_phi_1 = numpy.where(
        far, -inv_z * (1. - inv_z * (2. - 6. * inv_z)), z * phi - 1.)
    z2_phi = numpy.where(
        far, 2. * inv_z * (1. - inv_z * (3. - 12. * inv_z)),
        numpy.square(z) * phi - z + 1.)
    f_0 = -phi.imag / numpy.pi
    f_1 = (z_phi_1 / a).imag / numpy.pi
    f_2 = -(z2_phi / numpy.square(a)).imag / numpy.pi
    return f_0, f_1, f_2


def _combine_ikeda_carpenter(terms_alpha, f_beta_0, alpha, beta, r):
    """Weighted sum of the convoluted exponential terms of f(tau)."""
    f_0, f_1, f_2 = terms_alpha
    alpha_beta = alpha - beta
    fast = 0.5 * (1. - r) * alpha**3 * f_2
    slow = r * alpha**3 * beta / alpha_beta**3 * (
        f_beta_0 - f_0 - alpha_beta * f_1
        - 0.5 * numpy.square(alpha_beta) * f_2)
    return fast + slow


def tof_ikeda_carpenter(  # noqa: PLR0913
        alpha, beta, r, sigma, gamma, time, time_hkl, cutoff_fwhm: float = 0.):
    """Ikeda-Carpenter pulse convoluted with a pseudo-Voigt, TOF powder.

    The Gaussian (standard deviation ``sigma``) and the Lorentzian (FWHM
    ``gamma``) are combined into one Thompson-Cox-Hastings pseudo-Voigt
    with the common FWHM H and the mixing eta, the same convention as
    ``tof_Jorgensen_VonDreele`` (FullProf/CrysFML). The pulse f(tau) is
    given in the module docstring; its tail goes to longer times.

    Parameters
    ----------
    alpha, beta, r, sigma, gamma : numpy.ndarray
        Shape (n_time,): alpha, beta in 1/microseconds (alpha > beta),
        r dimensionless in [0, 1], sigma and gamma in microseconds,
        taken at the d of each time point.
    time : numpy.ndarray
        Time of flight of the points, shape (n_time,), microseconds.
    time_hkl : numpy.ndarray
        Time of flight of the reflections, shape (n_hkl,), microseconds.
    cutoff_fwhm : float
        Half window of a peak in FWHM (with 1/alpha + 1/beta added);
        0 means no cutoff.

    Returns
    -------
    numpy.ndarray
        Profile of unit area in time, shape (n_time, n_hkl), 1/microseconds.

    References: S. Ikeda, J. M. Carpenter, Nucl. Instrum. Methods A 239,
    536 (1985); FullProf manual, Npr=13; Mantid IkedaCarpenterPV.
    """
    from .powder_diffraction_cutoff import cutoff_place, cutoff_select
    from .powder_diffraction_tof import calc_hpv_eta, inv_8ln2
    cutoff_fwhm = numpy.inf if cutoff_fwhm <= 0. else cutoff_fwhm
    h_pv, eta = calc_hpv_eta(sigma * numpy.sqrt(8. * numpy.log(2.)), gamma)
    sigma_c = h_pv * numpy.sqrt(inv_8ln2)
    delta_2d = time[:, na] - time_hkl[na, :]
    half_width = cutoff_fwhm * (h_pv + 1. / alpha + 1. / beta)
    keep, (alpha, beta, r, eta, sigma_c, h_pv), delta_2d = cutoff_select(
        delta_2d, cutoff_fwhm, half_width, (alpha, beta, r, eta, sigma_c, h_pv))
    alpha, beta, r = alpha[:, na], beta[:, na], r[:, na]
    eta, sigma_c, h_pv = eta[:, na], sigma_c[:, na], h_pv[:, na]

    profile_g = _combine_ikeda_carpenter(
        _calc_gauss_exponential_terms(alpha, sigma_c, delta_2d),
        _calc_gauss_exponential_terms(beta, sigma_c, delta_2d)[0],
        alpha, beta, r)
    profile_l = _combine_ikeda_carpenter(
        _calc_lorentz_exponential_terms(alpha, h_pv, delta_2d),
        _calc_lorentz_exponential_terms(beta, h_pv, delta_2d)[0],
        alpha, beta, r)
    profile = (1. - eta) * profile_g + eta * profile_l
    return cutoff_place(keep, profile)
