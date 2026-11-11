---
title: Ordered magnetic structure (magnetic space group, mcif)
keywords: magnetic structure, second propagation vector, multi-k, two k vectors, второй вектор распространения, superspace group, mcif, MAGNDATA, msCIF, суперпространственная группа, propagation vector, k vector, satellites, magnetic domains, irreducible representation, incommensurate, вектор распространения, сателлиты, магнитные домены, ordered moment, moment, mcif, magnetic space group, shubnikov, representation analysis, choose representation by data, выбор представления, chi_sq of representations, antiferromagnet, ferromagnet, symmform, магнитная структура, магнитный момент, магнитная группа
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
above (the moments become the cos components). Refinement, flip ratios
and MEM give the same results for both descriptions. A commensurate
structure with k ≠ 0 is written in its magnetic supercell by
`make_magnetic_supercell` (P1 with all atoms, or with `use_symmetry True`
its magnetic space group and the atoms of the asymmetric unit).

The easiest way from a paramagnetic crystal with `space_group`:

0. procedure `find_propagation_vectors`: the propagation vectors of the
   measured satellites (fractional indices of the single-crystal
   reflections, k and -k merged) with the number of their satellites in
   every experiment; with several vectors they are added one after the
   other (multi-k structure);
1. method `report_irreps_of_k` of the crystal: the group of k, the characters
   of its irreducible representations Γ1, Γ2, ... (one-dimensional first,
   then those of dimension 2, 3, ...), the superspace group of each and the
   Fourier components (symmforms) it allows for the magnetic atoms. A
   multidimensional representation has several modes: the isotropy
   subgroups (order-parameter directions up to domains), each with its own
   superspace group, the largest first;
   With measured data (the nuclear structure fitted in the paramagnetic
   space group), procedure `rhochi_rank_magnetic_representations` (k,
   magnetic atoms, start moment in μB) adds k with every representation
   and every mode in a copy of the object, refines the allowed Fourier
   components against all experiments and lists chi_sq/n, R_F2 of the
   satellites of this k (not changed by the satellites of other vectors),
   the number of refined Fourier parameters, the magnetic domains and the
   moments, the best first. The sums of two one-dimensional representations
   (Γ1⊕Γ2) are tried too (`try_superpositions`), and single-crystal
   experiments get the magnetic domains of every model with equal fractions
   (`with_domains`): the arms of the star of k change the satellites. The
   recommended model has the fewest Fourier parameters among those with
   chi_sq/n within 5 % of the best; `keep_best` keeps it (with its domains)
   in the object. The other parameters keep their refinement flags: fix
   the nuclear structure, the scale and the extinction from the nuclear
   refinement, and give the experiments of one crystal the same scale
   (equal marks, `40.5()s`), otherwise the size of the moments is lost;
2. method `add_propagation_vector` (k, magnetic atoms, number of the
   representation, start moment, and `mode` for a multidimensional
   representation): converts the crystal; a magnetic atom without
   `atom_site_scat` (crystal from a nuclear refinement) gets it from
   `magnetic_ion` (as Tb3+) or from a charged `type_symbol`, with `lande`
   or, when it is 0, g_J of Hund's rules for 4f ions (Tb3+: 1.5) and 2 for
   the others; atoms whose orbit
   splits in the superspace group become several atoms (O1, O1_2, ...);
   method `report_magnetic_moments` then lists the moments of all magnetic
   atoms of the cell: the range of |M| over the cells and the kind of
   modulation (amplitude-modulated, helix, cycloid, spiral); with several
   propagation vectors also the total moment of every position (all waves
   added), to compare with g·J of the free ion (9 μB for Tb3+). A crystal
   with k gets a second vector by the same method (multi-k): the common
   group is lower and may free components of the first wave;
3. method `define_magnetic_domains` of a single-crystal experiment: the
   magnetic domains as `jana_twin` with equal fractions
   (`report_magnetic_domains` of the crystal lists them);
4. procedure `simulation_single_crystal`: simulated intensities of the
   nuclear reflections and of the satellites at the given `wavelength` (Å)
   in the range `ttheta_min`..`ttheta_max` (degrees); `satellites_only`
   leaves out the nuclear reflections.
   `simulation_unpolarized_neutron_powder_diffraction` and
   `simulation_polarized_neutron_powder_diffraction` give a powder profile
   with the satellites (also `ttheta_step`). The recipe
   `magnetic_structure_k` does these steps.

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

A second propagation vector (multi-k structure, 3+2 superspace group) is
added to a crystal that has one already with the same method
`add_propagation_vector`: `report_irreps_of_k` with the new k lists the
representations of the group of the new k within the present superspace
group, the chosen one adds the wave of the new k (`cell_wave_vector` row 2,
`atom_site_fourier_wave_vector`, `wave_vector_seq_id` of the Fourier rows),
and the waves of the first k are kept. The moment of an atom is the sum of
the waves: M = Σ_n (cos_n cos(2π k_n·r) + sin_n sin(2π k_n·r)). k equal to
±k of the first wave is refused; operations of the superspace group that
exchange the two vectors are not supported.

A commensurate structure (k with denominators up to 24) is exported to its
magnetic supercell with the method `make_magnetic_supercell` of the
crystal: a new data block with the supercell in P1 and the moment of
every atom (the assistant adds it to the object).
