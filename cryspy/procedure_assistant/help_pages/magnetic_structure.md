---
title: Ordered magnetic structure (magnetic space group, mcif)
keywords: magnetic structure, superspace group, mcif, MAGNDATA, msCIF, суперпространственная группа, propagation vector, k vector, satellites, magnetic domains, irreducible representation, incommensurate, вектор распространения, сателлиты, магнитные домены, ordered moment, moment, mcif, magnetic space group, shubnikov, antiferromagnet, ferromagnet, symmform, магнитная структура, магнитный момент, магнитная группа
categories: atom_site_moment, space_group_symop_magn_operation, space_group_symop_magn_centering, cell_wave_vector, atom_site_moment_fourier, atom_site_fourier_wave_vector, space_group_symop_magn_ssg_operation, space_group_symop_magn_ssg_centering, parent_space_group
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
is described as in msCIF (the files of the Bilbao MAGNDATA database and of
Jana2006 load directly):

- `cell_wave_vector` (`seq_id`, `x`, `y`, `z`): the propagation vector;
- `space_group_symop_magn_ssg_operation` and
  `space_group_symop_magn_ssg_centering` (`algebraic` such as
  `-x1,x2+1/2,-x3,-x4+1/4,+1`): the magnetic superspace group. It replaces
  `space_group` and describes the nuclear and the magnetic structure
  together, so every moment sits on an atom of `atom_site`;
- `atom_site_moment_fourier`: one row per magnetic atom and axis
  (`atom_site_label`, `axis` x, y or z along a, b, c, `param.cos`,
  `param.sin` in μB, refinable, and `param.cos_symmform`,
  `param.sin_symmform` such as `mxs1`, `-mxs1`, `0` with the symmetry
  constraints): M(x4) = cos·cos(2π x4) + sin·sin(2π x4), x4 = k·r. Without
  symmforms the constraints are calculated from the site symmetry;
- `parent_space_group`: the paramagnetic group (cell metric, domains).

A magnetic space group with `atom_site_moment` is the special case k = 0:
method `convert_to_superspace` of the crystal writes it in the categories
above (the moments become the cos components).

The easiest way from a paramagnetic crystal with `space_group`:

1. method `report_irreps_of_k` of the crystal: the group of k, the characters
   of its one-dimensional irreducible representations Γ1, Γ2, ..., the
   superspace group of each and the Fourier components (symmforms) it
   allows for the magnetic atoms;
2. method `add_propagation_vector` (k, magnetic atoms, number of the
   representation, start moment): converts the crystal; atoms whose orbit
   splits in the superspace group become several atoms (O1, O1_2, ...);
   method `report_magnetic_moments` then lists the moments of all magnetic
   atoms of the cell: the range of |M| over the cells and the kind of
   modulation (amplitude-modulated, helix, cycloid, spiral);
3. method `define_magnetic_domains` of a single-crystal experiment: the
   magnetic domains as `jana_twin` with equal fractions
   (`report_magnetic_domains` of the crystal lists them);
4. procedure `simulation_single_crystal`: simulated intensities of the
   nuclear reflections and of the satellites. The recipe `magnetic_structure_k`
   does these steps.

Single-crystal reflections of the satellites are given by fractional indices
H = h ± k in `diffrn_refln` (e.g. `0.25 0.25 1`). There is no ambiguity of ±k:
for a general k only one of H - k and H + k is integer; when 2k is integer,
both terms scatter at the same H and are added. Powder patterns
(constant wavelength, 1D) and single-crystal intensities and flip ratios
(with the chiral term) are calculated; time-of-flight and 2D powder
experiments give an error. Fourier parameters the data cannot define (the
global phase of the modulation without operations mapping k to -k, or
combinations of cos and sin when 2k is integer) are kept fixed
automatically, with a warning naming them.

k is given as fractions ('1/3') or decimals; a component within 0.001 of a
fraction with a denominator up to 12 is taken as that fraction. A cycloid
or a helix made of two modes is built with `irrep` = '1,2' (Γ1 ⊕ Γ2).
