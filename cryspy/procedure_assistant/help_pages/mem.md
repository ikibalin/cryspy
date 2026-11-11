---
title: Maximum entropy method (MEM)
keywords: mem, maximum entropy, superspace, supercell, propagation vector, spin density, magnetization density, reconstruction, den file, prior, метод максимальной энтропии, спиновая плотность, плотность намагниченности, восстановление
categories: mem_parameters, atom_electron_configuration
---
MEM reconstructs the spin or magnetization density from polarized neutron
data. The parameters are in the `mem_parameters` item placed after `global_`:

- `points_a`, `points_b`, `points_c`: grid of the unit cell (48 by default);
- `channel_col True`: spin density with moments along the field;
  `channel_ani True`: magnetization density with induced moments χ·H or
  ordered moments; `channel_nucl True`: nuclear density;
- `magnetization_plus`, `magnetization_minus`: total magnetization of the cell
  for spins along and opposite to the field;
- `prior_density`: `uniform` or `core` (spherical or multipole density from
  the crystal, see `atom_electron_configuration` and `atom_rho_multipole`);
- `use_asymmetry`, `GoF_desired` (stop criterion);
- `spin_density_to_den_file`, `magnetization_density_to_den_file`: output `.den` files.

Run `cryspy.mempy_spin_density_reconstruction(obj)` or
`cryspy.mempy_magnetization_density_reconstruction(obj)`.

The magnetic structure may be given by a magnetic space group
(`space_group_symop_magn_operation`, `atom_site_moment`) or by a magnetic
superspace group (after `convert_to_superspace`, or MAGNDATA msCIF). For
k = 0 MEM takes the magnetic space group and the ordered moments from the
superspace operations and the cos components and gives the same result.
For k ≠ 0, commensurate or incommensurate, MEM reconstructs in the cell the
density of every magnetic atom (the shape of its moment, the magnetic form
factor); the moment of the atom in the cell L is Re(A exp(-2πi k·L)) with
its Fourier amplitude A, so the satellites H = h ± k (fractional indices)
are fitted. The density is not written in `.den` files then (the directions
of the moments change from cell to cell). For a commensurate k the
structure may also be written in its magnetic supercell with its magnetic
space group (method `make_magnetic_supercell` of the crystal with
`use_symmetry True`); MEM in the supercell needs the reflections in the
indices of the supercell. Start MEM of a modulated structure from the
`core` prior (`atom_electron_configuration`): from the uniform prior the
fit of unpolarized intensities may stop in a false minimum.
