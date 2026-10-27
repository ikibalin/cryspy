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
parentheses (`"2*a-b=0"`). In CrysPy Editor the right-click menu of the
object tree has "Add constraint...".

Marks are also the names of the parameters in a punishment function: a
Python expression in `_punishment_function` that is added to chi-square
(for example a soft restraint `"(a-2.0)**2"`).

Parameters fixed by symmetry (a coordinate on a mirror plane, `chi_22 =
chi_11` in a cubic group, ...) cannot be refined; CrysPy applies these
constraints itself.

## Errors

When a value cannot be understood (`_cell_length_a abc`, a non-integer
Miller index, broken quotes) CrysPy raises `cryspy.RcifFormatError`, which
names the tag and the value. A value outside the allowed list (for example an
unknown `_atom_site_adp_type`) is ignored with a `cryspy.CrysPyWarning`.

## Where to start

Take the [example](examples/index.md) closest to your experiment and change
it. In CrysPy Editor the right-click menu of the object tree adds data
blocks, loops and items.
