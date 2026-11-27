---
title: Conversion between FullProf files (.pcr) and rcif
keywords: fullprof, pcr, prf, dat, convert, conversion, import, export, cross-check, compare, check in fullprof, rietveld, codes
categories: pd, crystal, phase, pd_instr_resolution, pd_instr_reflex_asymmetry, pd_instr_axial_divergence
---
A FullProf refinement of a powder pattern can be checked in CrysPy and a
CrysPy model in FullProf:

```
obj = cryspy.load_fullprof("pbso4.pcr")         # data file from the .pcr
cryspy.add_fullprof_files(rcif, "pbso4.pcr")     # into an existing object
cryspy.save_fullprof_files(rcif, "fp/pbso4.pcr") # writes fp/pbso4.pcr and fp/pbso4.dat
```

The assistant tools are `add_fullprof_files` and `save_fullprof_files`.

Converted (both directions): nuclear phases (cell, space group, atoms,
occupancies, Biso or anisotropic displacements: betas of FullProf,
`Uani` of CrysPy), the scale of each phase, U V W X Y and the size
broadening (`igsize` = GauSiz), the zero shift (`offset_2theta`), the
transparency SySin, the background (points, or the polynomial of FullProf
written as 25 points), excluded regions, the empirical asymmetry
`p1 ... p4`, the axial divergence `S_L`, `D_L` with `AsyLim` as
`ttheta_max`, the wavelength, the radiation and Cthm. Parameters refined in
FullProf (code ≠ 0) are refined in CrysPy, and back. Data files: Ins = 0,
6 (D1A/D2B) and 10 (XYDATA); CrysPy writes Ins = 10.

Conventions:

- FullProf `Occ` is the occupancy times the multiplicity of the site over
  the general multiplicity; CrysPy has the plain occupancy.
- The X-ray scale of CrysPy is two times the scale of FullProf (CrysPy
  has the polarization factor with K = 0.5, FullProf 1 + Cthm cos² 2θ).
- U V W X Y of the first phase are the resolution of the instrument
  (`pd_instr_resolution`); other phases get the differences in the
  `phase` loop.
- The convention of the empirical asymmetry `p1 ... p4` differs from
  FullProf (issue #50): refine `p1 ... p4` after a conversion before
  comparing the fits.

Not converted: magnetic phases, TOF, single crystals, polarized neutrons
(the sum up + down is written to FullProf), other profiles than
pseudo-Voigt, strain (Shape1 and others), texture (March-Dollase of
FullProf and the texture of CrysPy are different functions), SyCos,
`magnetic_broadening`, `stacking_fault`, `texture_sh`. A .pcr file with
what CrysPy does not calculate is refused with a message; values that are
only approximated give a warning.
