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


def calc_transs(l_max, nn, zeta, sthovl):
    """ integral( r**nn * exp(-zeta r) * j_l(Qr) dr) / Q**l
    Q = 4 pi sint / lambda
    # ATTENTION: in one of my note it was written factor r**2 but it looks that it was mistake
    """
    Q = 4*numpy.pi*sthovl
    a = [0. for h in range(17)]
    ff = [0 for h in range(l_max+1)]
    d = Q**2 + zeta**2
    a[0] = 0.0
    a[1] = 1./d
    n = nn-1
    tz = 2.*zeta
    ts = 4.*numpy.pi
    for l in range(l_max+1):
        ll = l+1
        if (ll != 1):
            a[ll] = a[ll-1]*ts*l*1./d
            a[ll-1] = numpy.zeros_like(sthovl) 
        for nx in range(ll, n+1):
            i1 = nx
            i2 = nx+1
            i3 = i2+1
            a[i3-1] = (tz*nx*a[i2-1] - (nx+l)*(nx-ll)*a[i1-1])*1./d
        ff[ll-1] = a[i3-1]
    return ff


def calc_jl_per_2sthovll(sthovl, coeff, n, zeta, kappa=1, l_max = 3):
    c_ij = numpy.expand_dims(coeff, axis=1) * numpy.expand_dims(coeff, axis=0)
    zeta_ang = kappa*zeta/0.529177
    hh = numpy.power(2.*zeta_ang, n.astype(float)+0.5)
    coeff_norm = hh / numpy.sqrt(scipy.special.factorial(2*n))
    coeff_norm_ij = numpy.expand_dims(coeff_norm, axis=1) * numpy.expand_dims(coeff_norm, axis=0)
    
    nn_ij = numpy.expand_dims(n, axis=1) + numpy.expand_dims(n, axis=0) # FIXME: why not -2?
    zeta_ij = numpy.expand_dims(zeta_ang, axis=1) + numpy.expand_dims(zeta_ang, axis=0)

    res = numpy.zeros(sthovl.shape+(l_max+1,))

    for c, cc, nn, zz in zip(c_ij.flatten(), coeff_norm_ij.flatten(), nn_ij.flatten(), zeta_ij.flatten()):
        ff = calc_transs(l_max, nn, zz, sthovl)
        ff_l = numpy.stack(ff, axis=1)
        
        val = c*cc*ff_l
        res += val
    return res


def calc_jl(sthovl, coeff, n, zeta, kappa=1, l_max = 3):
    # q_l = numpy.power(numpy.expand_dims(4*numpy.pi*sthovl, axis=1), numpy.expand_dims(numpy.arange(l_max+1), axis=0))
    q_l = numpy.power(numpy.expand_dims(2.*sthovl, axis=1), numpy.expand_dims(numpy.arange(l_max+1), axis=0))
    jl_per_2sthovll = calc_jl_per_2sthovll(sthovl, coeff, n, zeta, kappa=kappa, l_max=l_max)
    jl = jl_per_2sthovll*q_l
    return jl