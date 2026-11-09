---
title: Ordered magnetic structure (magnetic space group, mcif)
keywords: magnetic structure, propagation vector, k vector, satellites, magnetic domains, irreducible representation, incommensurate, вектор распространения, сателлиты, магнитные домены, ordered moment, moment, mcif, magnetic space group, shubnikov, antiferromagnet, ferromagnet, symmform, магнитная структура, магнитный момент, магнитная группа
categories: atom_site_moment, space_group_symop_magn_operation, space_group_symop_magn_centering, magnetic_structure_k, atom_site_moment_fourier, space_group_symop_magn_k, space_group_symop_magn_k_centering, parent_space_group
---
An ordered (commensurate) magnetic structure is described in the crystal data
block in the mcif way:

- instead of `space_group`, give the magnetic space group by the loops
  `space_group_symop_magn_operation` (`xyz` such as `-y,x-y+1/2,z,-1`) and
  `space_group_symop_magn_centering`;
- `atom_site_moment` gives the moment of each magnetic ion: `crystalaxis_x`,
  `crystalaxis_y`, `crystalaxis_z` in μB along the crystal axes (refinable) and
  `symmform` (for example `0,0,mz` or `mx,mx,0`) with the symmetry constraints;
- `atom_site_scat` gives the magnetic form factor.

## Magnetic structures with a propagation vector k

A structure with k (commensurate such as k = (1/2, 0, 0), or incommensurate)
is described in the crystal data block by:

- `space_group_symop_magn_k` (`operation_xyz`, `moment_uvw`, `phase` in units
  of 2π) and `space_group_symop_magn_k_centering`: the operations of the group
  of k. They replace `space_group` and describe the nuclear and the magnetic
  structure together, so every moment sits on an atom of `atom_site`;
- `magnetic_structure_k`: the propagation vectors;
- `atom_site_moment_fourier`: the Fourier components S of the magnetic atoms
  (`real_x`, ..., `imag_z` in μB along a, b, c, `phase`), FullProf convention
  m(R) = Σ_{±k} S exp(-2πi k·R), S_-k = S_k*;
- `parent_space_group`: the paramagnetic group, kept for the domains.

The easiest way from a paramagnetic crystal with `space_group`:

1. method `report_irreps_of_k` of the crystal: the group of k, the characters
   of its one-dimensional irreducible representations Γ1, Γ2, ... and the
   Fourier components each of them allows for the magnetic atoms;
2. method `add_propagation_vector` (k, magnetic atoms, number of the
   representation, start moment): converts the crystal; atoms whose orbit
   splits in the group of k become several atoms (O1, O1_2, ...);
3. method `define_magnetic_domains` of a single-crystal experiment: the
   magnetic domains as `jana_twin` with equal fractions
   (`report_magnetic_domains` of the crystal lists them);
4. procedure `simulation_single_crystal`: simulated intensities of the
   nuclear reflections and of the satellites. The recipe `magnetic_structure_k`
   does these steps.

Single-crystal reflections of the satellites are given by fractional indices
H = h ± k in `diffrn_refln` (e.g. `0.25 0.25 1`). There is no ambiguity of ±k:
for a general k only one of H - k and H + k is integer; when 2k is integer,
k and -k are the same vector and the sum has one term. Powder patterns
(unpolarized, constant wavelength) and unpolarized single-crystal intensities
are calculated; polarized data are not supported yet.
