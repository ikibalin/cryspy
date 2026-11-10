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
superspace group with k = 0 (after `convert_to_superspace`): MEM takes the
magnetic space group and the ordered moments from the superspace
operations and the cos components, and gives the same result. A density
with k ≠ 0 is not periodic in the cell: for a commensurate k write the
structure in its magnetic supercell (method `make_magnetic_supercell` of
the crystal) and give the reflections in the indices of the supercell.
