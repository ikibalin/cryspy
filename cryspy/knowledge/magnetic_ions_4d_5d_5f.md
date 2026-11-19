---
title: Magnetic 4d, 5d and actinide (5f) ions
keywords: 4d, 5d, 5f, actinide, uranium, ruthenium, iridium, osmium, spin-orbit, jeff, magnetic ion, valence
---
# Magnetic 4d, 5d and actinide (5f) ions

4d and 5d ions have extended orbitals: strong covalency, small ordered
moments (often 0.3–1 μB) and, for 5d ions, a strong spin-orbit coupling
(Ir4+ jeff = 1/2). Their form factors fall off slowly with sinθ/λ. The
CrysPy table has only neutral and 1+ ions of the 4d metals (Zr, Nb, Mo, Tc,
Ru, Rh, Pd) and no 5d ions: give the coefficients of the ion in the
atom_type_scat loop (neutron_magnetic_j0_a1 ... j2_e).

Actinides: the 5f shell is between the localized 4f and the itinerant 3d
cases. Uranium is U4+ (5f2) in oxides (UO2) and close to U3+ (5f3) in many
intermetallics, where the moments are reduced (URu2Si2: 0.03 μB). The
orbital moment is large and opposite to the spin one, so the form factor
differs strongly from a spin-only one. The Lande g of the actinides below is
that of Hund's rules (Russell-Saunders coupling), a first approximation:
the intermediate coupling of 5f ions changes it; refine it or give it.
The 4d and 5d ions get the spin-only value 2.

| Ion | Usual | Configuration | Spin S | Lande g | Moment (μB) | Notes |
|---|---|---|---|---|---|---|
| Ru4+ | yes | 4d4 | 1 | 2 | 2 | Ruthenates: small moments. No form factor in the CrysPy table. |
| Ru3+ | | 4d5 | 1/2 | 2 | 1 | Low spin (α-RuCl3). No form factor in the CrysPy table. |
| Rh4+ | yes | 4d5 | 1/2 | 2 | 1 | No form factor in the CrysPy table. |
| Mo5+ | yes | 4d1 | 1/2 | 2 | 1 | No form factor in the CrysPy table. |
| Mo3+ | | 4d3 | 3/2 | 2 | 3 | No form factor in the CrysPy table. |
| Os5+ | yes | 5d3 | 3/2 | 2 | 3 | Strong spin-orbit coupling. No form factor in the CrysPy table. |
| Ir4+ | yes | 5d5 | 1/2 | 2 | 1 | jeff = 1/2 Mott insulators (Sr2IrO4, Na2IrO3); moments 0.2–0.5 μB. No form factor in the CrysPy table. |
| U4+ | yes | 5f2 | 1 | 4/5 | 3.20 | Oxides (UO2). |
| U3+ | | 5f3 | 3/2 | 8/11 | 3.27 | Close to the 5f state of many intermetallics, with reduced moments (URu2Si2: 0.03 μB). |
| U5+ | | 5f1 | 1/2 | 6/7 | 2.14 | |
| Np4+ | yes | 5f3 | 3/2 | 8/11 | 3.27 | Radioactive. |
| Np3+ | | 5f4 | 2 | 3/5 | 2.40 | Radioactive. |
| Pu3+ | yes | 5f5 | 5/2 | 2/7 | 0.71 | Radioactive; strong absorber (use 242Pu); Pu metal is not magnetic. |
| Pu4+ | | 5f4 | 2 | 3/5 | 2.40 | Radioactive. |
