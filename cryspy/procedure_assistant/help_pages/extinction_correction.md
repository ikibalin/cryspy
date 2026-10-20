---
title: Extinction
keywords: extinction, mosaicity, domain radius, gauss, lorentz, becker coppens, экстинкция, мозаичность
categories: extinction
---
Extinction is set by the `extinction` item of a single crystal (or powder)
experiment:

- `model`: `gauss` or `lorentz` (Becker-Coppens type), or `qani`;
- for `gauss` and `lorentz`: `mosaicity` in minutes and `radius` in micrometers;
- for `qani`: the anisotropic parameters `q_hh, q_kk, q_ll, q_hk, q_hl, q_kl`.

```
_extinction_model gauss
_extinction_mosaicity 100.0()
_extinction_radius 50.0()
```

Mosaicity and radius are often strongly correlated; refine one of them first.
