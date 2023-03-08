import numpy

from scipy.special import erfc, exp1

def calc_n(alpha, beta):
    y = 0.5 * alpha * beta / (alpha + beta)
    return y

def calc_eta(h_l, h_com):
    x = h_l/h_com
    y = 1.36603 * x - 0.47719 * numpy.square(x) + 0.11116 * numpy.power(x, 3)
    return y

def calc_h_com(h_g, h_l):
    y = numpy.power(
        numpy.power(h_g, 5)
        + 2.69269 * numpy.power(h_g, 4) * h_l
        + 2.42843 * numpy.power(h_g, 3) * numpy.square(h_l)
        + 4.47163 * numpy.square(h_g) * numpy.power(h_l, 3)
        + 0.07842 * h_g * numpy.power(h_l, 4)
        + numpy.power(h_l, 5), 0.2)
    return y

def calc_h_g(sigma):
    y = numpy.sqrt(8*numpy.log(2)) * sigma
    return y

def calc_gaussian(delta_t, sigma):
    y = numpy.exp(-0.5 * numpy.square(delta_t/sigma)) / (numpy.sqrt(2 * numpy.pi) * sigma)
    return y

def calc_lorentz(delta_t, h_l):
    y = h_l / (2 * numpy.pi * (numpy.square(delta_t) + numpy.square(0.5*h_l)))
    return y

def calc_u(delta_t, sigma, alpha):
    y = 0.5 * alpha * (alpha * numpy.square(sigma) + 2*delta_t)
    return y

def calc_v(delta_t, sigma, beta):
    y = 0.5 * beta * (beta * numpy.square(sigma) - 2*delta_t)
    return y

def calc_y(delta_t, sigma, alpha):
    y = (alpha * numpy.square(sigma) + delta_t) / (sigma * numpy.sqrt(2))
    return y

def calc_z(delta_t, sigma, beta):
    y = (beta * numpy.square(sigma) - delta_t) / (sigma * numpy.sqrt(2))
    return y

def calc_e_beta_g(delta_t, sigma, alpha, beta):
    n = calc_n(alpha, beta)
    u = calc_u(delta_t, sigma, alpha)
    v = calc_v(delta_t, sigma, beta)
    y = calc_y(delta_t, sigma, alpha)
    z = calc_z(delta_t, sigma, beta)

    exp_u = numpy.exp(u)
    exp_u[numpy.isinf(exp_u)] = 1e200
    exp_v = numpy.exp(v)
    exp_v[numpy.isinf(exp_v)] = 1e200

    res = n * (exp_u * erfc(y) + exp_v * erfc(z))
    return res

def calc_p(delta_t, h_l, alpha):
    y = alpha * delta_t + 0.5* 1j * alpha * h_l
    return y

def calc_q(delta_t, h_l, beta):
    y = - beta * delta_t + 0.5* 1j * beta * h_l
    return y

def calc_e_beta_l(delta_t, h_l, alpha, beta):
    n = calc_n(alpha, beta)
    p = calc_p(delta_t, h_l, alpha)
    q = calc_q(delta_t, h_l, beta)
    exp1_p = exp1(p)
    exp1_p[numpy.isinf(exp1_p)] = 0
    exp1_p[numpy.isnan(exp1_p)] = 0
    exp_p = numpy.exp(p)
    exp_p[numpy.isinf(exp_p)] = 1e200
    exp1_q = exp1(q)
    exp1_q[numpy.isinf(exp1_q)] = 0
    exp1_q[numpy.isnan(exp1_q)] = 0
    exp_q = numpy.exp(q)
    exp_q[numpy.isinf(exp_q)] = 1e200
    res = -2 * n / numpy.pi * (
        numpy.imag(exp_p * exp1_p) +
        numpy.imag(exp_q * exp1_q)
    )
    return res

def calc_profile(delta_t, sigma, h_l, alpha, beta_0, beta_1, r_0):
    h_g = calc_h_g(sigma)
    h_com = calc_h_com(h_g, h_l)
    eta = calc_eta(h_l, h_com)
    e_beta_0_g = calc_e_beta_g(delta_t, sigma, alpha, beta_0)
    e_beta_1_g = calc_e_beta_g(delta_t, sigma, alpha, beta_1)
    e_beta_0_l = calc_e_beta_l(delta_t, h_l, alpha, beta_0)
    e_beta_1_l = calc_e_beta_l(delta_t, h_l, alpha, beta_1)
    pv_0 = eta * e_beta_0_l + (1. - eta) * e_beta_0_g
    pv_1 = eta * e_beta_1_l + (1. - eta) * e_beta_1_g
    y = r_0 * pv_0 + (1.-r_0) * pv_1
    return y


# Z-code models

def calc_sigma_square(d, sigma_0_sq, sigma_1_sq, sigma_2_sq):
    y = sigma_0_sq + sigma_1_sq * numpy.square(d) + sigma_2_sq * numpy.power(d, 4)
    return y


def calc_h_l(d, gamma_0, gamma_1, gamma_2):
    y = gamma_0 + gamma_1 * d + gamma_2 * numpy.square(d)
    return y


def calc_r_0(d, r_01, r_02, r_03):
    """Switching function
    """
    y = r_01 * numpy.exp(-r_02 * numpy.power(d, -r_03))
    return y


def calc_alpha_prime(d, alpha_1, alpha_2):
    y = 1./(alpha_1 + alpha_2/d)
    return y

def calc_beta_0_prime(d, beta_00, beta_01):
    y = beta_00 + beta_01/d
    return y


def calc_beta_1_prime(d, beta_10):
    y = beta_10
    return y


def calc_profile_by_zcode_parameters(
        delta_t, d,
        sigma_0_sq, sigma_1_sq, sigma_2_sq,
        gamma_0, gamma_1, gamma_2,
        r_01, r_02, r_03,
        alpha_1, alpha_2,
        beta_00, beta_01,
        beta_10):
    r_0 = calc_r_0(d, r_01, r_02, r_03)
    alpha_prime = calc_alpha_prime(d, alpha_1, alpha_2)
    beta_0_prime = calc_beta_0_prime(d, beta_00, beta_01)
    beta_1_prime = calc_beta_1_prime(d, beta_10)
    sigma = numpy.sqrt(calc_sigma_square(d, sigma_0_sq, sigma_1_sq, sigma_2_sq))
    h_l = calc_h_l(d, gamma_0, gamma_1, gamma_2)

    h_g = calc_h_g(sigma)
    h_com = calc_h_com(h_g, h_l)
    alpha = 1./(alpha_prime * h_com)
    beta_0 = 1./(beta_0_prime * h_com)
    beta_1 = 1./(beta_1_prime * h_com)
    """
    alpha = alpha_prime
    beta_0 = beta_0_prime
    beta_1 = beta_1_prime
    """
    res = calc_profile(delta_t, sigma, h_l, alpha, beta_0, beta_1, r_0)
    return res