---
title: Strong neutron absorbers and isotopes for substitution
keywords: absorption, absorber, neutron absorption, cross section, isotope, isotopic substitution, enriched isotope, hot neutrons, short wavelength, Gd, Sm, Eu, Cd, B, Dy
---
# Strong neutron absorbers and isotopes for substitution

The absorption cross section of most nuclei is proportional to the time the
neutron spends near the nucleus, i.e. to 1/v: σ_abs(λ) = σ_abs(1.798 Å) ·
λ / 1.798 Å. The tabulated values are for 2200 m/s neutrons (λ = 1.798 Å;
Sears, Neutron News 3 (1992) 26, the values of the CrysPy database). Gd,
Sm, Eu, Cd and B have nuclear resonances at thermal energies: for them the
1/v law is only a guide, and hot neutrons (λ ≈ 0.5 Å, D9 at ILL) reduce the
absorption much more.

Absorption matters when μ·d is not small (μ the linear absorption
coefficient, d the path in the crystal): CrysPy calculates μ from the
composition (calc_absorption_coefficient) and estimates whether a
correction changes the fit (estimate_absorption_correction). Ways to deal
with a strong absorber, best first: an isotope with a small cross section
(column Isotope), a short wavelength, a thin plate or a small crystal, and
an absorption correction of the intensities.

| Element | σ_abs natural (barn) | Isotope | σ_abs of the isotope (barn) | Notes |
|---|---|---|---|---|
| Gd | 49700 | 160Gd | 0.77 | 157Gd (259000 b) and 155Gd absorb; 158Gd (2.2 b) is used too. The strongest natural absorber. |
| Sm | 5922 | 154Sm | 8.4 | 149Sm (42080 b) absorbs; 154Sm is the usual enriched isotope (144Sm, 0.7 b, is rare). |
| Eu | 4530 | 153Eu | 312 | 151Eu (9100 b) absorbs; 153Eu still absorbs: use a short wavelength as well. |
| Cd | 2520 | 114Cd | 0.34 | 113Cd (20600 b) absorbs. |
| B | 767 | 11B | 0.0055 | 10B (3835 b) absorbs; 11B is the usual substitution in borides and borates. |
| Dy | 994 | 162Dy | 194 | 164Dy (2840 b) absorbs most; a short wavelength is the usual remedy. |
| Li | 70.5 | 7Li | 0.0454 | 6Li (940 b) absorbs. |
| Ir | 425 | 193Ir | 111 | 191Ir (954 b) absorbs. |
| Hg | 372.3 | 204Hg | 0.43 | |
| In | 193.8 | 113In | 12 | 115In (202 b) is 96 % of indium: substitution is expensive. |
| Er | 159 | 168Er | 2.74 | |
| Rh | 144.8 | | | One isotope (103Rh). |
| Hf | 104.1 | 180Hf | 13.04 | |
| Tm | 100 | | | One isotope (169Tm). |
| Au | 98.65 | | | One isotope (197Au). |
| Lu | 74 | 175Lu | 21 | 176Lu (2065 b) absorbs. |
| Ho | 64.7 | | | One isotope (165Ho). |
| Ag | 63.3 | 107Ag | 37.6 | |
| Nd | 50.5 | | | Moderate: rarely substituted. |
| Co | 37.18 | | | One isotope (59Co). |
| Yb | 34.8 | | | Moderate: rarely substituted. |
| Cl | 33.5 | 37Cl | 0.433 | 35Cl (44.1 b) absorbs; matters in chlorides with few heavy atoms. |
| Tb | 23.4 | | | One isotope (159Tb); weak, an absorption correction is rarely needed. |
| Mn | 13.3 | | | One isotope (55Mn). |
