"""Functions to find roots."""
import math

def rootsearch(f, a: float, b: float, dx: float):
    """Seacrh one root in range [a, b]."""
    x1 = a
    f1 = f(a)
    x2 = a + dx
    f2 = f(x2)
    if x2 > b:
        x2 = b
        f2 = f(x2)
    while f1*f2 > 0.0:
        if x1 >= b:
            return None, None
        x1 = x2
        f1 = f2
        x2 = x1 + dx
        f2 = f(x2)
        if x2 > b:
            x2 = b
            f2 = f(x2)
    return x1, x2


def bisect(f, x1: float, x2: float, switch: int = 0, epsilon: float = 1.0e-12):
    """Bisecting."""
    f1 = f(x1)
    if f1 == 0.0:
        return x1
    f2 = f(x2)
    if f2 == 0.0:
        return x2
    if f1*f2 > 0.0:
        print('Root is not bracketed')
        return None
    n = int(math.ceil(math.log(abs(x2 - x1)/epsilon)/math.log(2.0)))
    for i in range(n):
        x3 = 0.5*(x1 + x2)
        f3 = f(x3)
        if (switch == 1) and (abs(f3) > abs(f1)) and (abs(f3) > abs(f2)):
            return None
        if f3 == 0.0:
            return x3
        if f2*f3 < 0.0:
            x1 = x3
            f1 = f3
        else:
            x2 = x3
            f2 = f3
    return (x1 + x2)/2.0


def calc_roots(f, a: float, b: float, eps: float = 0):
    """
    Find all roots.

    Function f is given in range [a,b].
    Roots should not be close than eps.

    Give also estimation of mistake if mistake of zero is defined.
    """
    if eps == 0:
        eps = abs(b-a)*1./100
    l_res = []
    x1, x2 = rootsearch(f, a, b, eps)
    while x1 is not None:
        a2 = x2
        root = bisect(f, x1, x2, 1)
        if root is not None:
            # res.append( (round(root,-int(math.log(eps, 10)))))
            l_res.append(root)
        x1, x2 = rootsearch(f, a2, b, eps)
    return l_res
