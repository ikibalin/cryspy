---
title: Powder experiment (pd, pd2d, tof)
keywords: powder, pd, rietveld, profile, background, resolution, u v w, asymmetry, axial divergence, finger cox jephcoat, fcj, s_l, d_l, texture, 2theta, offset, phase scale, pd2d, tof, time of flight
categories: pd_meas, pd_background, pd_instr_resolution, pd_instr_reflex_asymmetry, pd_instr_axial_divergence, phase, texture, range, exclude, setup, tof_parameters, tof_profile, tof_background, pd2d_meas
---
A constant wavelength 1d powder experiment is a `pd` data block:

- `setup`: `wavelength` (Å), `offset_2theta` (degrees, refinable), `field` (T),
  `radiation`, `K` (0 for neutrons, 0.5 for characteristic X-rays) and
  `cthm` (cos² 2θ of the monochromator) for the Lorentz-polarization factor;
- `pd_instr_resolution`: Caglioti `U`, `V`, `W` and Lorentz `X`, `Y`;
- `pd_instr_reflex_asymmetry`: `p1 ... p4`;
- `pd_instr_axial_divergence`: `s_l`, `d_l` (refinable) and `ttheta_max`
  (degrees, optional): the physical asymmetry of Finger, Cox and Jephcoat
  from the axial divergence, S/L and H/L as `S_L`, `D_L` of FullProf;
  peaks below 90 degrees get a tail to low angles, above 90 degrees to high
  angles; reflections above `ttheta_max` stay symmetric;
- `phase` loop: crystal `label`, `scale`, `igsize` (and optional `U ... Y`
  per phase);
- `pd_background` loop: points `2theta`, `intensity` (intensities refinable);
- `pd_meas` loop: `2theta` and `intensity`, `intensity_sigma` (unpolarized) or
  `intensity_plus`, `intensity_minus` with sigmas (polarized);
- optional `range`, `exclude`, `texture` (modified March function `g_1`, `g_2`,
  axis `h_ax k_ax l_ax`), `chi2`, `extinction`, `diffrn_radiation`,
  `magnetic_broadening` (extra Lorentzian width of the magnetic reflections
  of a phase, 0 for the nuclear width: magnetic domain size, see the help
  page on the extra broadening of magnetic reflections).
- optional `stacking_fault` (faults of the stacking of layers, Hendricks-
  Teller model with independent or correlated faults: only the reflections
  that see the fault shift are broad or asymmetric; see the help page on
  stacking faults and the procedure `diagnose_peak_widths`).

`pd2d` is the 2d version for polarized neutrons (`gamma_nu_...` data),
`tof` is a time-of-flight pattern with `tof_parameters` (zero, dtt1, ...),
`tof_profile` and `tof_background`.

Profile matching (Le Bail fit, free or fixed integrated intensities of
`pd_peak`) is the procedure `rhochi_profile_matching`, see its help page.
