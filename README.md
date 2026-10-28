# CrysPy

[![PyPI](https://img.shields.io/pypi/v/cryspy.svg)](https://pypi.org/project/cryspy/)
[![Python](https://img.shields.io/pypi/pyversions/cryspy.svg)](https://pypi.org/project/cryspy/)
[![Tests](https://github.com/ikibalin/cryspy/actions/workflows/ci.yml/badge.svg?branch=develop)](https://github.com/ikibalin/cryspy/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/178411703.svg)](https://doi.org/10.5281/zenodo.4271776)

**CrysPy** is a crystallographic Python library for the analysis of
polarized and unpolarized neutron diffraction data, and of X-ray diffraction
data. It handles powders and single crystals: nuclear and magnetic structure
refinement, local magnetic susceptibility, and the reconstruction of spin and
magnetization densities by the maximum entropy method.

## Features

- Refinement (RhoChi) of single-crystal and powder data: polarized and
  unpolarized neutrons, X-rays, 1D and 2D powder patterns, time of flight.
- Magnetic structures in magnetic space groups (`mcif`), local magnetic
  anisotropy by the susceptibility approach.
- Spin and magnetization density reconstruction by the maximum entropy
  method (MEMPy).
- Input and output in `.rcif`, a CIF-like format built on the IUCr Core CIF
  dictionary, so data and results can be archived and exchanged.
- Works in scripts and Jupyter notebooks on Windows, macOS and Linux.
  [CrysPy Editor](https://github.com/ikibalin/cryspy_editor) is the
  graphical user interface for CrysPy.

## Installation

```bash
pip install cryspy          # calculations only
pip install cryspy[plot]    # + matplotlib for the plot_* methods (Jupyter)
pip install cryspy[all]     # every optional extra
```

## Example

Refine a powder pattern of PbSO4 from [`examples/`](examples/):

```python
import cryspy

rcif = cryspy.load_file("examples/rhochi_unpolarized_powder_1d_PbSO4/"
                        "rhochi_unpolarized_powder_1d_PbSO4.rcif")
result = cryspy.rhochi_rietveld_refinement(rcif)   # refines the values marked with ()
print(result["chi_sq"])
print(rcif.crystal_pbso4.cell)

with open("refined.rcif", "w") as fid:
    fid.write(rcif.to_cif())
```

More examples, one folder per experiment type, are in [`examples/`](examples/).

## Links

- Documentation and tutorials: [cryspy.fr](https://www.cryspy.fr)
- Try CrysPy in a web browser, without installing it:
  [cryspy_webapp on Binder](https://mybinder.org/v2/gh/ikibalin/cryspy_webapp/HEAD?urlpath=voila%2Frender%2Fcryspy_webapp.ipynb)
- Graphical user interface: [CrysPy Editor](https://github.com/ikibalin/cryspy_editor)
- Bug reports and feature requests: [GitHub issues](https://github.com/ikibalin/cryspy/issues)
- Contributing: [CONTRIBUTING.md](CONTRIBUTING.md)

## License

MIT License, Copyright (c) 2018-2026 Iurii Kibalin. See [LICENSE](LICENSE).
