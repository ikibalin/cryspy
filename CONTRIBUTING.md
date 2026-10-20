# Contributing to CrysPy

Thank you for helping to improve CrysPy. This page describes how changes
get from a branch into a release. Code conventions (numpy style, the
`(value, dder)` derivative contract, the A → H layer order) are in
[CLAUDE.md](CLAUDE.md); please read it before writing code.

## Reporting a problem

Open an [issue](https://github.com/ikibalin/cryspy/issues) with the CrysPy
version (`python -c "import cryspy; print(cryspy.__version__)"`), the Python
version, a minimal `.rcif` file or script that shows the problem, and what
you expected instead. For a mismatch with FullProf, Jana or another program,
attach the input of both programs.

## Development setup

The developer environment is managed by [pixi](https://pixi.sh):

```bash
git clone https://github.com/ikibalin/cryspy
cd cryspy
pixi run test              # unit + functional + examples tests
```

Without pixi:

```bash
pip install -e .[dev]
python -m pytest tests/
```

Dependencies are declared only in `pyproject.toml`; after changing them,
update `pixi.lock` with `pixi install`. The version is computed from git
tags by `versioningit`; never write a version number by hand.

## Branches

| Branch | Purpose |
|---|---|
| `master` | Released code. Changed only by release and hotfix PRs. |
| `develop` | Integration branch. All regular work goes here. |
| feature branches | One per change, branched off `develop`. |

Name feature branches in kebab-case, optionally with the issue number:
`66-error-during-refinement-of-anisotropic-thermal-parameters`.
Delete the branch after its PR is merged.

## Pull requests

1. Branch off the latest `develop` and open the PR into `develop`.
2. Keep one purpose per PR; do not mix refactoring with physics changes.
3. A bug fix comes with a test that fails before the fix. A new feature
   comes with unit tests in `tests/basic_operations/` and, if it touches
   the refinement pipeline, a functional test in `tests/functional/`.
4. Old `.rcif` files must keep loading.
5. Add the label that describes the change. CI refuses a PR without
   one, and the label decides the next version number:

   | Label | Use for | Version bump |
   |---|---|---|
   | `[scope] significant` | breaking or major new capability | major |
   | `[scope] enhancement` | new feature | minor |
   | `[scope] bug` | bug fix | patch |
   | `[scope] maintenance` | refactoring, CI, dependencies | patch |
   | `[scope] documentation` | documentation only | post |

   `[bot] release` and `[bot] backmerge` are reserved for the automated PRs.
6. The PR title becomes a line of the release notes, so write it as a short
   imperative sentence ("Fix TOF Lorentz half-angle scale convention").
7. CI runs the unit, functional and examples tests; keep them green.

## Releases

Releases are made by the maintainer.

1. Run the **Release PR (develop → master)** workflow from the Actions tab
   and merge the PR it opens.
2. The **Release draft update** workflow keeps a draft GitHub release with
   the next version and notes grouped into Added / Fixed / Changed.
   Review and publish it; publishing creates the tag.
3. The tag triggers **Publish to PyPI**, and the published release triggers
   **Backmerge**, which opens a `master → develop` PR. Merge it.

A hotfix branches off `master`, goes into `master` by PR, and reaches
`develop` through the backmerge.
