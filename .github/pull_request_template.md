## What and why

<!-- One purpose per PR. Link the issue: "Fixes #123". -->

## Checklist

- [ ] The PR goes into `develop` and has a `[scope] …` label
      (bug, enhancement, maintenance, documentation or significant).
- [ ] A bug fix has a test that fails before the fix; a new feature has
      unit tests in `tests/basic_operations/` (and a functional test in
      `tests/functional/` if it touches the refinement).
- [ ] Old `.rcif` files still load; public function signatures are
      unchanged or stay backward compatible.
- [ ] New public functions have docstrings with shapes, units and the
      coordinate system.
- [ ] `ruff check .` and the tests pass locally.
