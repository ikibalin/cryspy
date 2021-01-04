"""
Functions:
    - transs
    - calc_GCF

"""
import numpy
import math


def transs(lmax: int, nn: int, zeta: float, sthovl: numpy.ndarray):
    """
    Calculate integral.

    $$
    ( r**nn * exp(-zeta r) * j_l(Qr) r**2 dr) / Q**l
    $$

    for given l

    Q = 4 pi sint / lambda

    """
    Q = 4. * math.pi * sthovl
    a = [0.*Q for h in range(17)]
    ff = numpy.zeros(shape=(sthovl.size, lmax+1), dtype=float)
    d = Q**2+zeta**2
    a[0] = 0.0*Q
    a[1] = 1./d
    n = nn-1
    tz = 2.*zeta
    ts = 4.*math.pi
    for l in range(lmax+1):
        ll = l+1
        if (ll != 1):
            a[ll] = a[ll-1] * ts * l * 1./d
            a[ll-1] = 0.0*Q
        for nx in range(ll, n+1):
            i1 = nx
            i2 = nx+1
            i3 = i2+1
            a[i3-1] = (tz*nx*a[i2-1] - (nx+l)*(nx-ll)*a[i1-1]) * 1./d
        ff[:, ll-1] = a[i3-1]
    return ff

def calc_GCF(ln1: list, ldzeta1: list, lcoeff1: list, kappa1: float,
             ln2: list, ldzeta2: list, lcoeff2: list, kappa2: float,
             sthovl: float, lmax: float):
    GCF = numpy.zeros(shape=(sthovl.size, lmax+1), dtype=float)
    for n1, dzeta1, coeff1 in zip(ln1, ldzeta1, lcoeff1):
        # transformation from atomic units to inverse angstrems
        zetaA1 = dzeta1*kappa1/0.529177 
        Norm1 = math.sqrt(((2.*zetaA1)**(2*n1+1))/math.factorial(2*n1))
        for n2, dzeta2, coeff2 in zip(ln2, ldzeta2, lcoeff2):
            # transformation from atomic units to inverse angstrems
            zetaA2 = dzeta2*kappa2/0.529177
            Norm2 = math.sqrt(((2.*zetaA2)**(2*n2+1))/math.factorial(2*n2))

            nn, zeta = (n1+n2), (zetaA1+zetaA2)
            GCF += coeff1*coeff2*Norm1*Norm2*transs(lmax, nn, zeta, sthovl)
    return GCF

