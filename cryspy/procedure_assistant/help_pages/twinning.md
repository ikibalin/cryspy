---
title: Twins
keywords: twin, twins, twinning, domain, twin law, jana, time reversal, antiferromagnetic domain
categories: jana_twin
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
