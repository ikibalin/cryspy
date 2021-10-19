import numpy
import scipy
import scipy.special


def calc_rho_normalized(radius, coeff, zeta, n, kappa=1):
    """
    Give normalized Slatter type function.

        R = ((2*zeta0)**(n0+0.5)/((2*n0)!)**0.5)*r**(n0-1)*exp(-zeta0*r)
        (see Clementi_1974)
    """
    zeta_ang = kappa*zeta/0.52918
    hh = numpy.power(2.*zeta_ang, n.astype(float)+0.5)
    coeff_norm = hh / numpy.sqrt(scipy.special.factorial(2*n))
    r_nm1 = numpy.power(
        numpy.expand_dims(radius, axis=1),
        numpy.expand_dims(n-1, axis=0))
    norm_density = (numpy.expand_dims(coeff*coeff_norm, axis=0)*r_nm1*
                    numpy.exp(-numpy.expand_dims(zeta_ang, axis=0)*
                              numpy.expand_dims(radius, axis=1))).sum(axis=1)
    return norm_density


def calc_density_spherical(radius, shell_population, shell_coeff, shell_zeta, shell_n, kappa):
    density_spherical = numpy.zeros_like(radius)
    for population, coeff, zeta, n in zip(shell_population, shell_coeff, shell_zeta, shell_n):
        rho_normalized = calc_rho_normalized(radius, coeff, zeta, n, kappa=kappa)
        density_spherical += 4.*numpy.pi*population*numpy.square(rho_normalized)
    return density_spherical