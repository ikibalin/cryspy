"""python -m cryspy runs without the removed RhoChi and MEM classes."""
import os
import subprocess
import sys

import cryspy

DIR_ROOT = os.path.dirname(os.path.dirname(cryspy.__file__))


def run_console(*l_arg, cwd=None):
    env = dict(os.environ)
    env["PYTHONPATH"] = DIR_ROOT + os.pathsep + env.get("PYTHONPATH", "")
    return subprocess.run(
        [sys.executable, "-m", "cryspy", *l_arg], cwd=cwd, env=env,
        capture_output=True, text=True, timeout=300)


def test_main_missing_file(tmp_path):
    res = run_console(cwd=tmp_path)
    assert res.returncode == 1, res.stderr
    assert "main.rcif" in res.stdout


def test_main_no_refined_parameters(tmp_path):
    f_name = tmp_path / "main.rcif"
    f_name.write_text("""data_Fe3O4
_cell_length_a 8.5
_cell_length_b 8.5
_cell_length_c 8.5
_cell_angle_alpha 90.
_cell_angle_beta 90.
_cell_angle_gamma 90.
_space_group_name_H-M_alt "F d -3 m"
""")
    res = run_console(cwd=tmp_path)
    assert res.returncode == 0, res.stderr
    assert "no refined parameters" in res.stdout
