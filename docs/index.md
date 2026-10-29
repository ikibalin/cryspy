# CrysPy

**CrysPy** is a Python library for the analysis of polarized and unpolarized
neutron diffraction data, and of X-ray diffraction data, of powders and single
crystals:

- refinement (RhoChi) of single-crystal and powder data: polarized and
  unpolarized neutrons, X-rays, 1D and 2D powder patterns, time of flight;
- magnetic structures in magnetic space groups, local magnetic anisotropy by
  the susceptibility approach;
- spin and magnetization density reconstruction by the maximum entropy
  method (MEMPy).

Data, model and results live in one [`.rcif` file](rcif-format.md), a CIF-like
format. [CrysPy Editor](https://github.com/ikibalin/cryspy_editor) is the
graphical user interface for CrysPy.

## Installation

```bash
pip install cryspy          # calculations only
pip install cryspy[plot]    # + matplotlib for the plot_* methods (Jupyter)
pip install cryspy[view3d]  # + py3Dmol: crystal.view_structure() in 3D
pip install cryspy[all]     # every optional extra
```

CrysPy needs Python 3.10 or newer, numpy, scipy and pycifstar.
Plots of data blocks, items and loops (`plots()`, `plot_*`) need matplotlib,
an optional extra: without it they raise an error that says what to install.

## First refinement

```python
import cryspy

rcif = cryspy.load_file("rhochi_unpolarized_powder_1d_PbSO4.rcif")
result = cryspy.rhochi_rietveld_refinement(rcif)   # refines the values marked with ()
print(result["chi_sq"])

with open("refined.rcif", "w") as fid:
    fid.write(rcif.to_cif())
```

The file is one of the [examples](examples/index.md). The
[guides](guides/refinement.md) explain each part of the model; the
[API reference](api/procedure_rhochi.md) lists all functions and classes.

## Messages, progress and errors

CrysPy writes nothing to the console by itself. It reports progress through
the standard `logging` module (logger `"cryspy"`), so a script decides what to
show:

```python
import logging
logging.basicConfig(level=logging.INFO, format="%(message)s")
```

A refinement can report every iteration and be stopped from a callback:

```python
def callback(info):
    print(info["iteration"], info["chi_sq_per_n"])
    return True        # return False to stop; the last iteration is kept

cryspy.rhochi_rietveld_refinement(rcif, callback=callback)
```

Errors raised by CrysPy derive from `cryspy.CrysPyError`: `RcifFormatError`
(a value or the text of an rcif file cannot be understood),
`RcifFileNotFoundError` and `RefinementError` (a procedure cannot run on the
given object). Ignored or replaced input gives a `cryspy.CrysPyWarning`.

## Citation and license

CrysPy is distributed under the MIT license. To cite it, use the
"Cite this repository" button on [GitHub](https://github.com/ikibalin/cryspy).
