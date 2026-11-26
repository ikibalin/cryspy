---
title: Stacking faults in powder: only some reflections broad or asymmetric (Hendricks-Teller model)
keywords: stacking faults, stacking disorder, layer disorder, polytypes, polytypism, OD structure, displacement faults, asymmetric peak, peak with a tail, broad hkl peaks, selective broadening, hkl dependent broadening, some peaks broad others sharp, diagnose peak widths, Hendricks Teller, DIFFaX, Warren, rods, correlated faults, twin domains by stacking
categories: stacking_fault
---
A crystal built of layers can have faults in their stacking: the next
layer is put on the previous one with an extra in-plane shift t
(polytypes, OD structures, twinning by layers). The reflections that do
not see the shift stay sharp. The others are spread along rods normal to
the layers: in a powder they give broad peaks (H with a component along
the normal) or a sharp front with a long tail to high angles (H in the
plane of the layers, e.g. 110 of layers normal to c), asymmetric where the
structure factor of the layer changes along the rod. Example: K3Fe(CN)6
(Pbcn, layers normal to c, two per cell related by the c-glide): all
reflections with k odd are spread (t = b/2), 00l, h0l and the reflections
with k even are sharp.

## Diagnosis

`diagnose_peak_widths(experiment, phase_label)` fits every group of
reflections with a sharp peak (the profile of the phase) and two broad
ones with free intensities (no structure is needed) and gives the part of
every group in the broad peaks. If the broad groups follow a fault rule
(k odd, h + k odd, h - k not 3n, ...) it says so (SELECTIVE_WIDTH), finds
the layers (also several layers per cell, e.g. "0 0 2") and the
probability that describe the lines best (STACKING_FAULTS) and writes the
call of `add_stacking_fault`. If all reflections are broad (ALL_BROAD) the
reason is the size or strain of the phase, not faults. `diagnose_fit`
gives the SELECTIVE_WIDTH line too. Run the profile matching (with the
zero shift and the cell) first, so that the peak positions are right.

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
_stacking_fault_probability_repeat
_stacking_fault_layer_cut
  Pbcn  0 0 2  0.0 0.5 0.0  0.196(13)  0.41(12)  .
```

- `normal_h normal_k normal_l`: Miller indices N of the layer planes, the
  distance between the layers is 1/|N|: "0 0 2" are two layers per c.
  The layers must be related by a symmetry operation of the crystal (a
  translation, a screw axis or a glide plane with R^T N = N and N.t = 1),
  it is found in the space group.
- `shift_x shift_y shift_z`: the fault vector t in the plane of the layers
  (N.t = 0), fractions of the cell, in the frame of the layer below. Its
  images by the symmetry of the layer are equally probable variants.
- `probability`: probability alpha of a fault after a normal contact
  (refinable, 0 ... 1).
- `probability_repeat`: probability beta of a fault after a fault
  (refinable; not given: equal to alpha, independent faults). beta > alpha:
  the faults come in runs (domains of the other stacking, twins by
  stacking); beta < alpha: they avoid each other.
- `layer_cut`: N.r of the cut between two layers (not given: the widest
  gap between the atoms). The cut decides which atoms move together.

## The calculation

The layer is cut from the structure; the ideal stacking operation
S0 = (R0, t0) is taken from the space group, a fault is S0 with the extra
shift. The operations form a Markov chain (alpha, beta) and the intensity
along every rod Q = H + tau P (N = g P, g layers per cell) is the
Hendricks-Teller sum of the correlations of the layers
(S. B. Hendricks and E. Teller, J. Chem. Phys. 10 (1942) 147; M. M. J.
Treacy, J. M. Newsam and M. W. Deem, Proc. R. Soc. Lond. A 433 (1991) 499,
DIFFaX), calculated with a transfer matrix of the states (orientation of
the layer, last operation), the structure factor of the layer and its
displacement parameters. For identical layers stacked by a translation
and independent faults it is the Lorentzian of B. E. Warren (X-ray
Diffraction, ch. 13) with beta_H = 1 - alpha + alpha exp(2 pi i H.t).

The points of the rods are given to the nearest reflection of the rod (a
spread reflection keeps its own line in the profile matching), put on a
fine 2theta grid with the Lorentz factor and convolved with the
instrumental profile (resolution and asymmetry p1-p4). Rietveld: the
intensity of the rods comes from the structure (the faults move
intensity between the reflections of a rod; the magnetic part of a
reflection keeps its Bragg intensity and is spread in the same way);
profile matching: the shapes come from the model, the intensities are
free. Neutrons and X-rays, 1D pd only, one kind of fault per phase.

A wave of the chain that no fault changes (|lambda| = 1, for example an
equivalent of a reflection that does not see the shift) stays a sharp
Bragg peak of that reflection; a reflection goes to the rods when any of
its equivalents sees the faults (the faults lower the symmetry). Polarized
neutrons: I+ and I- get the same line shapes, so I+ + I- is the
unpolarized pattern; the nuclear-magnetic interference is spread with the
lines.

## Procedures

- `diagnose_peak_widths(experiment, phase_label)`: the diagnosis above.
- `add_stacking_fault(normal, shift, probability, phase_label,
  experiment, probability_repeat, layer_cut)`: add the row (checks the
  layers and the shift), the probabilities are refined.
- `report_stacking_fault(experiment)`: the layers, the cut, the stacking
  operation, the fault shift and its variants, the probabilities with
  their errors, the correlation length of the stacking and the first
  spread reflections.

The recipe `stacking_faults` does these steps.
