"""matplotlib is an optional extra (cryspy[plot]): cryspy imports and
calculates without it, and plots say what to install."""
import subprocess
import sys

import pytest

from cryspy.A_functions_base.optional_import import (
    MissingOptionalDependency, optional_module)

# Run in a fresh interpreter where importing matplotlib and the packages
# it brings along fails (cryspy must not use them indirectly either).
BLOCK_MATPLOTLIB = """
import sys
BLOCKED = ("matplotlib", "pyparsing", "cycler", "kiwisolver", "PIL",
           "dateutil", "fontTools", "contourpy")
class Block:
    def find_spec(self, name, path=None, target=None):
        if name.split(".")[0] in BLOCKED:
            raise ImportError("blocked " + name)
sys.meta_path.insert(0, Block())
"""


def run_without_matplotlib(code):
    return subprocess.run(
        [sys.executable, "-c", BLOCK_MATPLOTLIB + code],
        capture_output=True, text=True, timeout=600)


def test_import_cryspy_without_matplotlib():
    result = run_without_matplotlib(
        "import cryspy\n"
        "cell = cryspy.Cell(length_a=5., length_b=5., length_c=5.)\n"
        "print(round(cell.volume, 3))\n")
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "125.0"


def test_plot_without_matplotlib_says_what_to_install():
    result = run_without_matplotlib(
        "from cryspy.C_item_loop_classes.cl_1_diffrn_refln import "
        "plot_scatters_with_labels\n"
        "plot_scatters_with_labels([0., 1.], [0., 1.])\n")
    assert result.returncode != 0
    assert "pip install cryspy[plot]" in result.stderr


def test_optional_module_error():
    module = optional_module("cryspy_no_such_module", extra="plot")
    with pytest.raises(MissingOptionalDependency, match=r"cryspy\[plot\]"):
        module.anything


def test_optional_module_submodule():
    matplotlib = optional_module("matplotlib")
    assert matplotlib.colors.Normalize(vmin=0, vmax=2)(1.) == 0.5
