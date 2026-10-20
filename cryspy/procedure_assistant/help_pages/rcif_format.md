---
title: Structure of an rcif file
keywords: rcif, cif, file, format, global, data block, loop, item, refine, sigma, mark, структура файла, формат, блок данных
categories: crystal, diffrn, pd, pd2d, tof
---
CrysPy reads and writes `.rcif` files: CIF syntax with CrysPy categories.

A file starts with `global_`. Items placed right after it belong to the whole
object (for example `mem_parameters` or `inversed_hessian`). Then follow data
blocks `data_<name>`. The kind of a data block is recognised by its content:

| Data block | Mandatory categories | Meaning |
| --- | --- | --- |
| crystal | `cell`, `atom_site` | crystal (or magnetic) structure |
| diffrn | `setup`, `diffrn_orient_matrix`, `diffrn_refln` | single crystal experiment |
| pd | `setup`, `pd_instr_resolution`, `phase`, `pd_background`, `pd_meas` | 1d powder experiment |
| pd2d | `setup`, `pd2d_instr_resolution`, `phase`, `diffrn_radiation`, `range`, `pd2d_background`, `pd2d_meas` | 2d powder experiment (polarized neutrons) |
| tof | `tof_parameters`, `tof_profile`, `phase`, `tof_meas` | time-of-flight powder experiment |

An experiment refers to crystals through the `phase` category: `_phase_label`
is the name of the crystal data block.

A single item is written as `_<category>_<attribute> value`, a list of items
as a `loop_` with one column per attribute.

Refinement flags are written in the value itself:

- `5.4321` the value is fixed;
- `5.4321()` the value is refined;
- `5.4321(12)` the value is refined and 12 is its standard deviation in the last digits.

To start a new object, take one of the examples distributed with CrysPy
(`examples/` folder of the repository) and change it. In CrysPy Editor the
right-click menu of the object tree adds data blocks, loops and items.
