---
title: Conversion between FullProf files (.pcr) and rcif
keywords: fullprof, pcr, prf, dat, int, single crystal, convert, conversion, import, export, cross-check, compare, check in fullprof, rietveld, codes
categories: pd, diffrn, crystal, atom_site_moment_fourier, extinction, phase, pd_instr_resolution, pd_instr_reflex_asymmetry, pd_instr_axial_divergence
---
A FullProf refinement of a powder pattern or of single crystal intensities can be checked in CrysPy and a
CrysPy model in FullProf:

```
obj = cryspy.load_fullprof("pbso4.pcr")         # data file from the .pcr
cryspy.add_fullprof_files(rcif, "pbso4.pcr")     # into an existing object
cryspy.save_fullprof_files(rcif, "fp/pbso4.pcr") # writes fp/pbso4.pcr and fp/pbso4.dat
obj = cryspy.load_fullprof("hote3.pcr")          # single crystal, Cry = 1, hote3.int
```

The assistant tools are `add_fullprof_files` and `save_fullprof_files`.

Converted (both directions): nuclear phases (cell, space group, atoms,
occupancies, Biso or anisotropic displacements: betas of FullProf,
`Uani` of CrysPy), the scale of each phase, U V W X Y and the size
broadening (`igsize` = GauSiz), the zero shift (`offset_2theta`), the
transparency SySin, the background (points, or the polynomial (Nba = 0) and the
Chebyshev polynomials (Nba = -5) of FullProf written as 25 points), excluded regions, the empirical asymmetry
`p1 ... p4`, the axial divergence `S_L`, `D_L` with `AsyLim` as
`ttheta_max`, the wavelength, the radiation and Cthm. Parameters refined in
FullProf (code ≠ 0) are refined in CrysPy, and back. Data files: Ins = 0,
6 (D1A/D2B) and 10 (XYDATA); CrysPy writes Ins = 10.

Single crystals (Cry = 1, integrated intensities in an .int file, both
directions): a nuclear phase, or a magnetic phase with a propagation
vector given by basis functions (Jbt = 1, Isy = -2, `BASR`/`BASI`). The
magnetic phase becomes a crystal in P1 with one atom for each operator and
the Fourier components of the moments (`atom_site_moment_fourier`), the
M-type form factor of FullProf (<j0>) a Lande factor 2. CrysPy writes a
magnetic crystal as such a phase in P1 with the basis functions along a, b,
c and the satellites only (Mag_Only). The scale, the cell, the wavelength
and the extinction are converted.

Conventions:

- FullProf `Occ` is the occupancy times the multiplicity of the site over
  the general multiplicity; CrysPy has the plain occupancy.
- The X-ray scale of CrysPy is two times the scale of FullProf (CrysPy
  has the polarization factor with K = 0.5, FullProf 1 + Cthm cos² 2θ).
- U V W X Y of the first phase are the resolution of the instrument
  (`pd_instr_resolution`); other phases get the differences in the
  `phase` loop.
- Single crystal: the scale of CrysPy is half the scale of FullProf
  (CrysPy adds the intensities of both spin states).
- Single crystal: the extinction model 1 of FullProf,
  y = (1 + 0.001 Ext1 F² λ³ / sin 2θ)^(-1/2), is the model `qani` with
  q_ii = Ext1/2 g*_ii and q_ij = Ext1 g*_ij (the reciprocal metric
  tensor). CrysPy applies it to each spin channel, so for magnetic
  reflections the intensities differ by up to about 1 %.
- Single crystal, magnetic: FullProf gives m = S exp(-2πi k·R) + c.c. at
  the unreduced position op(r); for half-integer k without c.c. and with a
  real S.
- The convention of the empirical asymmetry `p1 ... p4` differs from
  FullProf (issue #50): refine `p1 ... p4` after a conversion before
  comparing the fits.

Not converted: magnetic phases of powders, TOF, flipping ratios
(itypdata = 2) and polarized single crystals, the refinement of the
coefficients C_n of the basis functions (refine the Fourier components in
CrysPy), other extinction models of FullProf, twins and magnetic domains,
magnetic space groups (MCIF) for single crystals, polarized powders
(the sum up + down is written to FullProf), other profiles than
pseudo-Voigt, strain (Shape1 and others), texture (March-Dollase of
FullProf and the texture of CrysPy are different functions), SyCos,
`magnetic_broadening`, `stacking_fault`, `texture_sh`. A .pcr file with
what CrysPy does not calculate is refused with a message; values that are
only approximated give a warning.
