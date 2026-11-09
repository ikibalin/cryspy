---
title: Maximum entropy method (MEM)
keywords: mem, maximum entropy, spin density, magnetization density, reconstruction, den file, prior, метод максимальной энтропии, спиновая плотность, плотность намагниченности, восстановление
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
