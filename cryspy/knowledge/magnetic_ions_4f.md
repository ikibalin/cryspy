---
title: Magnetic rare-earth (4f) ions
keywords: rare earth, lanthanide, 4f, magnetic ion, valence, oxidation state, Lande factor, g_J, gJ, free-ion moment, Hund's rules, crystal field
---
# Magnetic rare-earth (4f) ions

In almost all compounds the rare earths are trivalent (R3+): the 4f shell
lies inside the 5s and 5p shells, so the ion keeps its free-ion ground
multiplet (Hund's rules) and only the crystal field splits it. The usual
ion is therefore R3+ also in intermetallics and in compounds whose formal
charge balance fails, such as the polytellurides RTe3 (TbTe3: Tb3+, the Te
layers carry the remaining charge). Exceptions: Eu and Yb are often
divalent (Eu2+ 4f7, Yb2+ 4f14 nonmagnetic), Ce can be tetravalent (Ce4+ 4f0,
nonmagnetic) or of intermediate valence, Sm is sometimes divalent.

The largest ordered moment of an ion is g_J·J (free-ion saturation moment,
column Moment); the crystal field usually reduces it, so a refined moment
clearly above g_J·J means an unphysical model (wrong scale, wrong
representation or free amplitudes of several waves). The Lande factor g_J
enters the dipole form factor f = <j0> + (2/g_J − 1)<j2>; CrysPy takes g_J of
this table (Lande g) when lande is 0. The values follow Hund's rules
(Jensen and Mackintosh, Rare Earth Magnetism, 1991, table 1.1).

| Ion | Usual | Configuration | Ground term | Lande g | Moment (μB) | Notes |
|---|---|---|---|---|---|---|
| Ce3+ | yes | 4f1 | 2F5/2 | 6/7 | 2.14 | Kondo and heavy-fermion compounds: ordered moments often well below 1 μB. No Ce3+ form factor in the CrysPy table: give atom_type_scat. |
| Ce4+ | | 4f0 | 1S0 | – | 0 | Not magnetic. |
| Ce2+ | | 4f2 | 3H4 | 4/5 | 3.20 | Rare. |
| Pr3+ | yes | 4f2 | 3H4 | 4/5 | 3.20 | Non-Kramers ion: a singlet crystal-field ground state may not order at all. |
| Pr4+ | | 4f1 | 2F5/2 | 6/7 | 2.14 | |
| Nd3+ | yes | 4f3 | 4I9/2 | 8/11 | 3.27 | Kramers ion. |
| Nd2+ | | 4f4 | 5I4 | 3/5 | 2.40 | |
| Pm3+ | yes | 4f4 | 5I4 | 3/5 | 2.40 | Radioactive. No form factor in the CrysPy table. |
| Sm3+ | yes | 4f5 | 6H5/2 | 2/7 | 0.71 | Strong neutron absorber (use 154Sm); the J-multiplets mix, so the dipole form factor is poor (its maximum is not at sinθ/λ = 0). |
| Sm2+ | | 4f6 | 7F0 | – | 0 | J = 0: Van Vleck paramagnet, no ordered moment of its own. |
| Eu2+ | yes | 4f7 | 8S7/2 | 2 | 7.0 | Strong neutron absorber (use 153Eu or hot neutrons). |
| Eu3+ | | 4f6 | 7F0 | – | 0 | J = 0: Van Vleck paramagnet, no ordered moment of its own. |
| Gd3+ | yes | 4f7 | 8S7/2 | 2 | 7.0 | Strongest natural neutron absorber (use 160Gd or hot neutrons); S-state ion, weak anisotropy. |
| Tb3+ | yes | 4f8 | 7F6 | 3/2 | 9.0 | Non-Kramers ion. |
| Tb4+ | | 4f7 | 8S7/2 | 2 | 7.0 | |
| Dy3+ | yes | 4f9 | 6H15/2 | 4/3 | 10.0 | Kramers ion; noticeable absorption (use a short wavelength). |
| Dy2+ | | 4f10 | 5I8 | 5/4 | 10.0 | |
| Ho3+ | yes | 4f10 | 5I8 | 5/4 | 10.0 | Non-Kramers ion. |
| Er3+ | yes | 4f11 | 4I15/2 | 6/5 | 9.0 | Kramers ion. |
| Tm3+ | yes | 4f12 | 3H6 | 7/6 | 7.0 | Non-Kramers ion. |
| Tm2+ | | 4f13 | 2F7/2 | 8/7 | 4.0 | |
| Yb3+ | yes | 4f13 | 2F7/2 | 8/7 | 4.0 | Kramers ion; often an effective spin 1/2 with an anisotropic g-tensor. |
| Yb2+ | | 4f14 | 1S0 | – | 0 | Not magnetic. |

La3+ (4f0) and Lu3+ (4f14) are not magnetic.
