"""Domains of a single crystal: structural twins (jana_twin) combined with
magnetic domains (magnetic_domain).

A domain is the whole crystal (atoms and ordered moments) transformed by a
magnetic operation (T, theta): T gives the Miller indices of the domain for
the measured ones, h_d = T h, and theta = -1 reverses the ordered moments.
Every structural twin has all the magnetic domains: the domain (s, m) has
the matrix T_m T_s, the time reversal theta_s theta_m and the volume
fraction f_s f_m.
"""
import numpy

IDENTITY = numpy.eye(3).flatten()


def take_domain_set(dict_diffrn: dict, prefix: str):
    """Matrices (9, n), time reversals (n,) and fractions (n,) of one set
    of domains with the identity first (its fraction is 1 - sum of the
    others); None when the experiment has none.

    prefix : 'twin' (jana_twin) or 'magnetic_domain' (magnetic_domain).
    """
    key_fraction = f"{prefix:}_fraction"
    if key_fraction not in dict_diffrn:
        return None
    fraction = numpy.atleast_1d(numpy.asarray(dict_diffrn[key_fraction], dtype=float))
    matrices = numpy.asarray(dict_diffrn[f"{prefix:}_matrices"], dtype=float).reshape(9, -1)
    time_reversal = numpy.asarray(dict_diffrn.get(
        f"{prefix:}_time_reversal", numpy.ones(fraction.size)), dtype=float)
    return (numpy.concatenate((IDENTITY[:, numpy.newaxis], matrices), axis=1),
            numpy.concatenate(([1., ], time_reversal)),
            numpy.concatenate(([1. - fraction.sum(), ], fraction)))


def form_domains_by_dictionary(dict_diffrn: dict):
    """All domains of a diffrn experiment: structural twins times magnetic
    domains.

    Parameters
    ----------
    dict_diffrn : dict
        Dictionary of the experiment with "twin_fraction", "twin_matrices"
        (jana_twin) and/or "magnetic_domain_fraction",
        "magnetic_domain_matrices", "magnetic_domain_time_reversal".

    Returns
    -------
    dict or None
        None without domains, otherwise
        "matrices" (9, n_domain) T on hkl, the first is the identity;
        "time_reversal" (n_domain,); "fraction" (n_domain,) = f_s f_m;
        "fraction_twin" (n_twin,), "fraction_magnetic" (n_magnetic,) with
        the identity first. The domain d = s * n_magnetic + m.
    """
    twin = take_domain_set(dict_diffrn, "twin")
    magnetic = take_domain_set(dict_diffrn, "magnetic_domain")
    if twin is None and magnetic is None:
        return None
    one = (IDENTITY[:, numpy.newaxis], numpy.ones(1), numpy.ones(1))
    twin = one if twin is None else twin
    magnetic = one if magnetic is None else magnetic
    m_twin = twin[0].transpose().reshape(-1, 3, 3)
    m_magnetic = magnetic[0].transpose().reshape(-1, 3, 3)
    # T_m T_s for every pair (s, m), s outer
    matrices = numpy.einsum("mij,sjk->smik", m_magnetic, m_twin).reshape(-1, 9).transpose()
    return {
        "matrices": matrices,
        "time_reversal": numpy.outer(twin[1], magnetic[1]).flatten(),
        "fraction": numpy.outer(twin[2], magnetic[2]).flatten(),
        "fraction_twin": twin[2], "fraction_magnetic": magnetic[2]}


def calc_dder_domain_fractions(iint_domains, fraction_twin, fraction_magnetic):
    """Derivatives of the intensities summed over the domains by the refined
    fractions (the first fraction of each set is 1 - sum of the others).

    Parameters
    ----------
    iint_domains : numpy.ndarray, shape (n_hkl, n_twin * n_magnetic)
        Intensities of every domain, domain d = s * n_magnetic + m.
    fraction_twin, fraction_magnetic : numpy.ndarray
        Fractions with the identity first.

    Returns
    -------
    tuple of numpy.ndarray
        (n_hkl, n_twin - 1) by the twin fractions and (n_hkl,
        n_magnetic - 1) by the magnetic domain fractions.
    """
    iint = iint_domains.reshape(iint_domains.shape[0], fraction_twin.size, fraction_magnetic.size)
    # sum over the magnetic domains of each twin, sum over the twins of each magnetic domain
    iint_twin = (iint * fraction_magnetic[numpy.newaxis, numpy.newaxis, :]).sum(axis=2)
    iint_magnetic = (iint * fraction_twin[numpy.newaxis, :, numpy.newaxis]).sum(axis=1)
    return (iint_twin[:, 1:] - iint_twin[:, 0:1], iint_magnetic[:, 1:] - iint_magnetic[:, 0:1])
