# The rcif format

CrysPy keeps the model, the measured data, the refinement settings and the
results in one text file with the extension `.rcif`. Its syntax is the
syntax of CIF (STAR); the categories and their tags are the ones of CrysPy,
close to the IUCr Core CIF dictionary where such a category exists. The
complete list of categories, generated from the code, is on the page
[rcif categories](rcif-categories.md).

## Structure of a file

```text
global_rhochi                      # whole object; the name is free

_punishment_function "(a-0.4)**2"  # items of the whole object (optional)

data_PbSO4                         # crystal
_cell_length_a 8.4797(3)
...

data_pd                            # experiment that uses the crystal
loop_
_phase_label
_phase_scale
PbSO4  0.0145(2)
...
```

A file starts with `global_` followed by an optional name. Items right after
it belong to the whole object (for example `mem_parameters` for a MEM
reconstruction, `inversed_hessian` written after a refinement). Then follow
data blocks `data_<name>`; the name is used to refer to the block (an
experiment refers to a crystal through `_phase_label`).

The kind of a data block is not written in the file: CrysPy recognises it by
the categories it contains.

| Data block | Mandatory categories | Meaning |
| --- | --- | --- |
| crystal | `cell`, `atom_site` | crystal (or magnetic) structure |
| diffrn | `setup`, `diffrn_orient_matrix`, `diffrn_refln` | single-crystal experiment |
| pd | `setup`, `pd_instr_resolution`, `phase`, `pd_background`, `pd_meas` | 1D powder experiment |
| pd2d | `setup`, `pd2d_instr_resolution`, `phase`, `diffrn_radiation`, `range`, `pd2d_background`, `pd2d_meas` | 2D powder experiment (polarized neutrons) |
| tof | `tof_parameters`, `tof_profile`, `phase`, `tof_meas` | time-of-flight powder experiment |

A block without these categories is still read and written back, as a
generic data block.

## Items and loops

A single value is written as a tag and a value on one line. The tag is
`_<category>_<attribute>`; the dot form `_<category>.<attribute>` of
DDLm is accepted too, but use one form within a category.

```text
_cell_length_a     8.4797
_cell_length_b     5.3984
_setup_wavelength  1.912
```

A table is a `loop_`: the list of tags, then one row per entry, values
separated by spaces.

```text
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_adp_type
_atom_site_B_iso_or_equiv
Pb Pb 0.1876(1) 0.25 0.1670(1) Biso 1.4()
S  S  0.0654    0.25 0.6845    Biso 0.4
```

Text with spaces is put in quotes (`'...'` or `"..."`), multi-line text
between lines that start with `;`. `#` starts a comment. Tags are not case
sensitive.

`.` and `?` mean "no value": the attribute stays undefined.

Units are those of the category description (Å and degrees for the cell,
fractions of the cell for atomic coordinates, Å² for `B_iso`, μB for magnetic
moments, ...); they are listed on the [categories page](rcif-categories.md).

## Refinement flags, standard deviations and constraints

Whether a parameter is refined is written in the value itself:

| Written | Meaning |
| --- | --- |
| `5.4321` | fixed value |
| `5.4321()` | refined value, no standard deviation yet |
| `5.4321(12)` | refined value; 12 is its standard deviation in the last digits, i.e. 0.0012 |

After a refinement CrysPy writes the new value with its standard deviation,
rounded to two significant digits of the deviation.

A letter after the parentheses is a **mark**. Refined values with the same
mark are refined as one parameter; a mark that ends with `-` changes the sign:

```text
_atom_site_moment_crystalaxis_x  2.1()a
_atom_site_moment_crystalaxis_y  2.1()a     # always equal to the x component
_atom_site_moment_crystalaxis_z -2.1()a-    # always the opposite
```

A linear constraint between different marks, for example a sum of two
occupancies equal to 1, is written in the loop `constraint` of the global
object; the equation holds exactly during the refinement:

```text
global_

loop_
_constraint_id
_constraint_expression
1 "x+y=1"

data_fe_co
...
loop_
_atom_site_label
_atom_site_occupancy
Fe1  0.7()x
Co1  0.3()y
```

An expression may use numbers, marks, `+`, `-`, `*` and `/` by numbers and
parentheses (`"2*a-b=0"`).

A constraint with `_constraint_sigma` is soft (a restraint): `f = g` adds
`((f-g)/sigma)**2` to chi-square and may be any formula of the marks
(`**`, `sqrt`, `exp`, `log`, `abs`, `sin`, `cos`, `tan`, `sind`, `cosd`),
for example `"sqrt(a**2+b**2)=2.5"` with sigma `0.05`; a hard constraint has
`.` as its sigma. In CrysPy Editor the right-click menu of the object tree
has "Add constraint...", with the check box "Soft constraint" for a soft one.

Marks are also the names of the parameters in a punishment function: a
Python expression in `_punishment_function` that is added to chi-square
(for example a soft restraint `"(a-2.0)**2"`).

Parameters fixed by symmetry (a coordinate on a mirror plane, `chi_22 =
chi_11` in a cubic group, ...) cannot be refined; CrysPy applies these
constraints itself.

## Magnetic structures with a propagation vector

A magnetic structure with a propagation vector `k` is described as in msCIF
(J. M. Perez-Mato et al., J. Phys.: Condens. Matter 24 (2012) 163201), the
format of the Bilbao MAGNDATA database and of Jana2006; MAGNDATA files with
one or several wave vectors load directly. The moment of an atom at the position `r`
(fractional, the lattice translation included) is

`M(x4) = cos * cos(2 pi x4) + sin * sin(2 pi x4)`, `x4 = k.r`,

along the unit vectors of a, b, c in μB. `cell_wave_vector` gives `k`; the
magnetic superspace group (`space_group_symop_magn_ssg_operation` and
`_centering`, operations `x1,x2,x3,x4,theta`) is the only symmetry of the
crystal: it generates the nuclear and the magnetic structure together, so
`space_group` is not given and every moment belongs to an atom of
`atom_site`. The magnetic form factor comes from `atom_site_scat`.
`atom_site_moment_Fourier` has one row per atom and axis with the refinable
`cos` and `sin` and their symmetry constraints (`symmform`: equal symbols are
equal parameters, `0` is zero, symbols are local to the atom); without
symmforms the constraints are calculated from the site symmetry.
`parent_space_group` (optional) gives the metric of the cell and the
magnetic domains.

```text
_parent_space_group.name_H-M_alt "I m m m"

loop_
_cell_wave_vector_seq_id
_cell_wave_vector_x
_cell_wave_vector_y
_cell_wave_vector_z
1 0.5 0.0 0.5

loop_
_space_group_symop_magn_ssg_operation.id
_space_group_symop_magn_ssg_operation.algebraic
1 x1,x2,x3,x4,+1
2 -x1,x2,-x3,-x1-x3+x4,-1
3 -x1,-x2,-x3,-x1-x3+x4,+1
4 x1,-x2,x3,x4,-1

loop_
_space_group_symop_magn_ssg_centering.id
_space_group_symop_magn_ssg_centering.algebraic
1 x1,x2,x3,x4,+1
2 x1+1/2,x2+1/2,x3+1/2,x4,+1
3 x1,x2,x3,x4+1/2,-1
4 x1+1/2,x2+1/2,x3+1/2,x4+1/2,-1

loop_
_atom_site_moment_Fourier.atom_site_label
_atom_site_moment_Fourier.axis
_atom_site_moment_Fourier.wave_vector_seq_id
_atom_site_moment_Fourier_param.cos
_atom_site_moment_Fourier_param.sin
Ho x 1 -0.0745() 0.1009()
Ho y 1 0.0 0.0
Ho z 1 -5.342() 7.236()
```

The full example (Ho2BaNiO5, the FullProf example `hobk`) is in
`examples/rhochi_unpolarized_powder_1d_k_Ho2BaNiO5`. Unpolarized 1D powder
diffraction with constant wavelength, single-crystal intensities and
polarized flip ratios are calculated; magnetic satellites are not
yet listed in the `pd_peak` loop, and time-of-flight and 2D powder
experiments give an error.

**Parameters not defined by the data.** When `2k` is integer (for example
`(1/2, 0, 1/2)`), the atoms see only a few values of `x4`, so some
combinations of `cos` and `sin` change no structure factor. Without
operations mapping k to -k (epsilon = -1), a shift of the modulation along
`x4` (the global phase) changes no intensity. CrysPy finds these parameters
and keeps them fixed during the refinement, with a warning naming them.

**k = 0.** A magnetic space group (`space_group_symop_magn_operation`,
`_centering`) with `atom_site_moment` is the special case `k = 0`:
`Crystal.convert_to_superspace()` writes it in the categories above (the
moments become the `cos` components) and gives the same unpolarized
intensities.

**Single crystals.** The satellites are given in `diffrn_refln` by fractional
indices `H = h ± k` (for example `0.25 0.25 1`), the nuclear reflections by
integer ones. There is no ambiguity between `+k` and `-k`: for a general `k`
only one of `H - k` and `H + k` is integer, and when `2k` is integer both
terms scatter at the same `H` and are added. (FullProf `.int` files give
integer `hkl` with the number of `±k` instead.) Magnetic domains are the
twins of `jana_twin`: the intensity is `sum_d f_d I(M_d H)`. The example
`examples/rhochi_unpolarized_single_crystal_k_Cs3Cr2Br9` uses the CrysFML
data of Cs3Cr2Br9, `k = (1/4, 1/4, 0)`.

**From a paramagnetic crystal.** The methods of `Crystal`
`report_irreps_of_k(k_1, k_2, k_3, magnetic_atoms)` (the group of `k`, its
irreducible representations, their superspace groups and the symmforms of
the magnetic atoms) and `add_propagation_vector(k_1, k_2, k_3,
magnetic_atoms, irrep, moment, mode)` convert a crystal with `space_group`:
the operations of the group of `k` get `tau = k.t + arg(chi) / 2pi` of the
representation (for a multidimensional representation, the operations of
the isotropy subgroup `mode` with `D(g) eta = lambda_g eta` for its
order-parameter direction `eta`, and `tau = k.t + arg(lambda_g) / 2pi`;
the representations are found by the method of Dixon, Numer. Math. 10
(1967) 446, the isotropy subgroups as in ISODISTORT), an operation mapping `k` to `-k` is added when the
representation allows it, atoms whose orbit splits become several atoms,
and the space group is kept as `parent_space_group`.
Called on a crystal that has a superspace group already,
`add_propagation_vector` adds a second propagation vector (a (3+2)
superspace group): the operations of the present group that keep the new
`k` get the phases of the chosen representation in a new `x5`, the atoms
whose orbit splits become several atoms, the Fourier rows of the new wave
have `wave_vector_seq_id` 2 (`atom_site_fourier_wave_vector` gives
`q1_coeff`, `q2_coeff`), and the waves of the first `k` are kept. The moment
is the sum of the waves, `M = sum_n (cos_n cos(2 pi x(3+n)) + sin_n sin(2 pi
x(3+n)))`, `x(3+n) = k_n.r`. Superspace operations that exchange the
propagation vectors (a non-diagonal epsilon block) are not supported.
`Crystal.make_magnetic_supercell()` writes a commensurate structure as a
P1 magnetic supercell with the moments of all atoms.
`Crystal.report_magnetic_domains` and `Diffrn.define_magnetic_domains(crystal)`
give the domains (cosets of the point group of the superspace group;
domains related by the inversion give the same unpolarized intensities and
are merged). `simulation_single_crystal` adds the satellites to the
simulated reflections. The assistant recipe `magnetic_structure_k` does
these steps.

## Errors

When a value cannot be understood (`_cell_length_a abc`, a non-integer
Miller index, broken quotes) CrysPy raises `cryspy.RcifFormatError`, which
names the tag and the value. A value outside the allowed list (for example an
unknown `_atom_site_adp_type`) is ignored with a `cryspy.CrysPyWarning`.

## Where to start

Take the [example](examples/index.md) closest to your experiment and change
it. In CrysPy Editor the right-click menu of the object tree adds data
blocks, loops and items.
