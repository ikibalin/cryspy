---
title: Twins
keywords: twin, twins, twinning, domain, twin law, jana, time reversal, antiferromagnetic domain, magnetic domain, 180 degree domain
categories: jana_twin, magnetic_domain
---
Twins of a single crystal are described like in Jana2006 by the `jana_twin`
loop of the `diffrn` experiment. Each row is one more domain:

- `id`;
- `matrix[1][1] ... matrix[3][3]` transform Miller indices of the first
  domain into the indices of this domain;
- `fraction` is the refinable volume fraction of this domain; the first domain
  has 1 − sum of the fractions.

- `time_reversal` (optional, +1 by default or -1) makes the twin law a
  magnetic operation: the ordered moments of this domain are reversed.

A domain is the whole crystal (atoms and ordered moments in its own axes)
transformed by the twin law. The magnetic field and the polarization of the
beam are laboratory vectors: they are taken in the axes of every domain.
Field-induced moments (`atom_site_susceptibility`) follow the field of the
domain and do not depend on `time_reversal`.

In a field the ferromagnetic component of a domain is usually along the
field: when the net ordered moment of a domain is antiparallel to its field
CrysPy warns and `time_reversal` -1 may be needed. For an antiferromagnet
the net moment is zero and the data decide: the assistant tool
`check_twin_time_reversal` calculates chi-square for +1 and -1 of every
domain of the flip ratio experiments and keeps the better one.

Magnetic domains are written separately, in the loop `magnetic_domain` of
the experiment: the crystal transformed by a magnetic operation of the
paramagnetic group that is not in the magnetic group. Each row is one more
domain (the first one, the crystal itself, is always present):

- `id`;
- `operation_xyz` is the magnetic operation as in mCIF, e.g. `x,y,z,-1`
  (the 180° domain: all ordered moments reversed) or `-y,x,z,+1`;
- `fraction` is the refinable volume fraction; the first domain has
  1 − sum of the fractions.

Structural twins (`jana_twin`) and magnetic domains are combined: every
twin has all the magnetic domains, the domain (twin s, magnetic m) has the
matrix T_m T_s on hkl (T = R^T of the operation), the time reversal
theta_s theta_m and the fraction f_s f_m. Method `define_magnetic_domains`
of the experiment writes the magnetic domains of the crystal (cosets of the
magnetic group in the grey paramagnetic group, Friedel pairs merged), for
k = 0 magnetic space groups and for propagation vectors; older files with
the magnetic domains in `jana_twin` keep loading, and this method moves
them to `magnetic_domain`.
