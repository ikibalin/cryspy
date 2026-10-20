# CLAUDE.md

Guidance for Claude Code (and human contributors) working in this repository.

## Project

CrysPy is a Python library for crystallographic analysis of neutron (and X-ray)
diffraction data: polarized/unpolarized single crystal and powder refinement
(RhoChi), TOF, magnetic space groups (mcif), and MEM reconstruction of
magnetization/spin density. Input and output are CIF-like `.rcif` files parsed
with `pycifstar`. Runtime dependencies: `numpy`, `scipy`, `pycifstar`,
`matplotlib`. Supported Python: 3.9+.

## Commands

```bash
pixi run unit-tests        # tests/basic_operations/
pixi run functional-tests  # tests/functional/ (parallel, -n 2)
pixi run examples-tests    # tests/test_scripts.py, runs examples/ and checks chi_sq
pixi run test              # all of the above
```

Without pixi: `pip install -e .[dev]`, then `python -m pytest <path>`.
Run a single test with `python -m pytest tests/path/test_x.py::test_name -v`.

## Architecture

The package is layered; lower letters must not import from higher ones.

| Layer | Contents |
|---|---|
| `A_functions_base` | Pure numerical functions on numpy arrays (unit cell, symmetry, structure factor, profiles, extinction, form factors). No classes from other layers. |
| `B_parent_classes` | Generic `ItemN`, `LoopN`, `DataN`, `GlobalN` (CIF item / loop / data block / global). |
| `C_item_loop_classes` | Concrete CIF categories (`cl_1_cell.py`, `cl_1_atom_site.py`, ...). |
| `D_functions_item_loop` | Functions operating on item/loop objects. |
| `E_data_classes` | Data blocks: `Crystal`, `MagCrystal`, `Diffrn`, `Pd`, `Pd2d`, `TOF`, ... |
| `F_functions_data`, `H_functions_global` | Functions on data blocks and on the global container. |
| `procedure_rhochi`, `procedure_mempy`, `procedure_simulation` | Top-level procedures (refinement, MEM, simulation). |

Refinement goes through the dictionary API (`procedure_rhochi/rhochi_by_dictionary.py`):
objects are converted to plain dicts (`crystal_*`, `pd_*`, `diffrn_*`, ...) and
the numerical core works on those dicts. New physics goes into `A_functions_base`
first, then gets wired into the dict pipeline, then into the CIF classes.

Public API is re-exported from `cryspy/__init__.py`.

## Code conventions

- Match the surrounding code: `import numpy` (not `np`), `na = numpy.newaxis`,
  snake_case names, `calc_<quantity>_by_<inputs>` for numerical functions.
- Vectorize with numpy. Arrays carry the component axis first
  (e.g. `index_hkl` has shape `(3, n_hkl)`, a 3x3 matrix is stored as 9
  components along axis 0); extra axes are for reflections, atoms, symmetry
  elements. Avoid Python loops over reflections or points.
- Derivatives: numerical functions take `flag_<param>: bool` arguments and
  return `(value, dder)`, where `dder` is a dict of derivatives keyed by
  parameter name (`"m_real"`, `"m_imag"` for complex inputs). Keep this
  contract when adding functions so they compose with existing ones.
- Coordinate systems must be explicit. The Cartesian system is x || a*, z || c.
  State in the docstring which system (fractional, Cartesian, reciprocal) and
  which units (Å, Å⁻¹, μB, degrees vs radians) every argument uses.
- Every new public function gets a docstring with parameters, shapes and units,
  and a reference (paper, International Tables, FullProf/Jana convention) when
  implementing a published formula.

## Software engineering practices

**Scope and changes**
- Keep each change focused on one purpose; don't mix refactoring with physics
  changes in one commit or PR.
- Read the existing implementation before writing new code; extend existing
  helpers (`matrix_operations.py`, `symmetry_elements.py`, `unit_cell.py`, ...)
  instead of duplicating them.
- Prefer the simplest correct solution. Don't add abstractions, options or
  configuration that no current caller needs.
- Don't change public function signatures or `.rcif` formats without keeping
  backward compatibility; old input files must keep loading.

**Correctness and numerics**
- Validate inputs at the boundaries (CIF parsing, public API) with clear error
  messages; trust internal calls.
- Compare floats with tolerances (`numpy.allclose`, `pytest.approx`), never
  with `==`. Choose tolerances deliberately and document why.
- Watch for degenerate cases: zero-length vectors, sin θ/λ = 0, empty
  reflection lists, special positions, complex vs real dtypes.
- When a formula has an analytical derivative, check it against a finite
  difference in a test.
- Never silence an exception or warning to make a test pass; fix the cause.

**Testing**
- Every bug fix comes with a test that fails before the fix and passes after.
- Every new feature comes with fast unit tests in `tests/basic_operations/`
  on small arrays, plus a functional test in `tests/functional/` with a
  minimal `.rcif` when it touches the refinement pipeline.
- Verify against an independent reference where possible (analytical special
  case, FullProf, Jana, published data), not only against the code's own output.
- Keep tests deterministic: fixed random seeds, no network, no dependence on
  the working directory (build paths from `os.path.dirname(__file__)`).
- Run the relevant test suite before declaring a task done; report failures
  honestly with their output.

**Performance**
- Measure before optimizing. Vectorize first; cache results that are reused
  across refinement iterations (see `dict_in_out` and
  `flag_use_precalculated_data`).
- Don't trade readability for speed without a measured benefit.

**Readability**
- Names should state the physical quantity (`sthovl`, `f_m_perp`,
  `unit_cell_parameters`), not generic `x`, `tmp`, `data2`.
- Comments explain *why* (physics, convention, non-obvious trick), not *what*.
- Keep functions short and single-purpose; split long formulas into named
  intermediate quantities.
- Delete dead code instead of commenting it out; git keeps history.

**Security and hygiene**
- Never commit secrets, credentials or personal data.
- Don't commit generated or scratch files. Use `tmp/` (git-ignored) for
  experiments, outputs and throwaway scripts.
- Don't add dependencies without a strong reason; declare any new one in
  `pyproject.toml` (the single source of truth), then update `pixi.lock`.

## Git workflow

The full contributor workflow (labels, version bumps, releases) is in
`CONTRIBUTING.md`; the essentials:

- `master` is the release branch; `develop` is the integration branch.
  Branch off `develop` and open PRs into `develop`.
- Branch names are kebab-case, optionally prefixed by the issue number
  (`66-error-during-refinement-of-anisotropic-thermal-parameters`).
- Commit messages: short imperative summary line ("Fix TOF Lorentz half-angle
  scale convention"), with details in the body when needed.
- Every PR needs at least one label (`[scope] bug`, `[scope] enhancement`,
  `[scope] documentation`, `[scope] maintenance`, `[scope] significant`);
  CI checks this.
- CI (`.github/workflows/ci.yml`) runs unit, functional and examples tests
  via pixi on every push and PR; keep them green.
- Versions come from git tags via `versioningit`; don't hardcode versions.
