---
title: Random stacking faults in powder: only some reflections broad or asymmetric
keywords: stacking faults, stacking disorder, layer disorder, polytypes, polytypism, OD structure, displacement faults, asymmetric peak, peak with a tail, broad hkl peaks, selective broadening, hkl dependent broadening, some peaks broad others sharp, diagnose peak widths, Hendricks Teller, Warren, rods
categories: stacking_fault
---
A crystal built of layers can have random faults in their stacking: from
time to time all the layers above a plane are shifted by an extra vector
t (polytypes, OD structures, twinning by layers). The reflections H with
H.t an integer do not see the faults and stay sharp. The others are spread
along the normal of the layers: in a powder they give broad, nearly
symmetric peaks (H with a component along the normal) or a sharp front
with a long tail to high angles (H in the plane of the layers, e.g. 110 of
layers normal to c). Example: K3Fe(CN)6 (Pbcn, layers normal to c, two per
cell): all reflections with k odd are spread (t = b/2), 00l, h0l and the
reflections with k even are sharp.

## Diagnosis

`diagnose_peak_widths(experiment, phase_label)` fits every group of
reflections with a sharp peak (the profile of the phase) and two broad
ones with free intensities (no structure is needed) and gives the part of
every group in the broad peaks. If the broad groups follow a fault rule
(k odd, h + k odd, h - k not 3n, ...) it says so (SELECTIVE_WIDTH), finds
the layers and the probability that describe the lines best
(STACKING_FAULTS) and writes the call of `add_stacking_fault`. If all
reflections are broad (ALL_BROAD) the reason is the size or strain of the
phase, not faults. `diagnose_fit` gives the SELECTIVE_WIDTH line too. Run
the profile matching (with the zero shift and the cell) first, so that
the peak positions are right.

## The loop stacking_fault

```
loop_
_stacking_fault_phase_label
_stacking_fault_normal_h
_stacking_fault_normal_k
_stacking_fault_normal_l
_stacking_fault_shift_x
_stacking_fault_shift_y
_stacking_fault_shift_z
_stacking_fault_probability
  Pbcn  0 0 2  0.0 0.5 0.0  0.178(1)
```

- `normal_h normal_k normal_l`: Miller indices of the layer planes, the
  distance between the layers is 1/|N|: "0 0 2" are two layers per c.
- `shift_x shift_y shift_z`: the fault vector t (fractions of the cell),
  in the plane of the layers.
- `probability`: probability alpha of a fault between two neighbour
  layers (refinable, 0 ... 1).

## The calculation

For a reflection H the correlation of two layers m apart is multiplied
by beta^m, beta = 1 - alpha + alpha exp(2 pi i H.t), so the Bragg point
becomes a Lorentzian along the normal with the half width
gamma = -ln|beta| / (2 pi d_L) and the shift -arg(beta) / (2 pi d_L)
(1/Å, s = 1/d units; d_L the distance between the layers; S. B. Hendricks
and E. Teller, J. Chem. Phys. 10 (1942) 147; B. E. Warren, X-ray
Diffraction, ch. 13). For t = b/2 and k odd beta = 1 - 2 alpha. Every
equivalent reflection has its own rod (the faults break the symmetry).
The powder average of the Lorentzian is calculated exactly on a fine
grid (with the 1/sqrt(s² - g²) front of the reflections in the plane of
the layers) and convolved with the instrumental profile (resolution and
asymmetry p1-p4 of the experiment). The integrated intensity of every
reflection is kept, so the same works in the Rietveld refinement
(structure factors) and in the profile matching (free intensities).

The Lorentzian form holds while gamma is smaller than the distance
between the Bragg points along the rod (alpha up to about 0.3 for
t = b/2). Only one kind of fault per phase, 1D pd only.

## Procedures

- `diagnose_peak_widths(experiment, phase_label)`: the diagnosis above.
- `add_stacking_fault(normal, shift, probability, phase_label,
  experiment)`: add the row, the probability is refined.
- `report_stacking_fault(experiment)`: the layers, the fault shift, the
  probability with its error, the correlation length of the stacking and
  the first spread reflections.

The recipe `stacking_faults` does these steps.
