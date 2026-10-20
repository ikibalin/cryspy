---
title: Twins
keywords: twin, twins, twinning, domain, twin law, jana, двойники, двойникование, домены
categories: jana_twin
---
Twins of a single crystal are described like in Jana2006 by the `jana_twin`
loop of the `diffrn` experiment. Each row is one more domain:

- `id`;
- `matrix[1][1] ... matrix[3][3]` transform Miller indices of the first
  domain into the indices of this domain;
- `fraction` is the refinable volume fraction of this domain; the first domain
  has 1 − sum of the fractions.

Known limitation: with two or more twin laws the reflections are combined
incorrectly in the current version; one twin law works.
