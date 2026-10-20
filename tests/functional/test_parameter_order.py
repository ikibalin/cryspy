"""The order of refined parameters must not depend on PYTHONHASHSEED.

rhochi_calc_chi_sq_by_dictionary collected the parameter names with
list(set(...)): the order of the gradient and of the Hessian changed with
the hash seed of each Python run, and with a flat minimum the refinement
reached a different point from one run to the next.
"""
import json
import os
import subprocess
import sys

F_NAME = os.path.join(os.path.dirname(__file__), "twin_tbte3", "main.rcif")
DIR_CRYSPY = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

SCRIPT = """
import json, sys
import cryspy
from cryspy.procedure_rhochi.rhochi_by_dictionary import rhochi_calc_chi_sq_by_dictionary
global_dict = cryspy.load_file(sys.argv[1]).get_dictionary()
# several parameters of the crystal and of the experiment, as in a refinement
for key in ("flags_atom_b_iso", "flags_atom_fract_xyz"):
    global_dict["crystal_tbte3"][key][...] = True
for key in ("flags_phase_scale", "flags_extinction_radius", "flags_twin_fraction"):
    global_dict["diffrn_t1p8nuc"][key][...] = True
names = rhochi_calc_chi_sq_by_dictionary(
    global_dict, dict_in_out={})[-1]
print(json.dumps([list(name) if isinstance(name, tuple) else name for name in names], default=str))
"""


def take_parameter_names(hash_seed: str) -> list:
    env = dict(os.environ, PYTHONHASHSEED=hash_seed)
    env["PYTHONPATH"] = DIR_CRYSPY + os.pathsep + env.get("PYTHONPATH", "")
    result = subprocess.run(
        [sys.executable, "-c", SCRIPT, F_NAME], env=env, capture_output=True,
        text=True, check=True)
    return json.loads(result.stdout.strip().splitlines()[-1])


def test_parameter_order_does_not_depend_on_hash_seed():
    l_names = [take_parameter_names(seed) for seed in ("1", "2", "11", "19")]
    assert len(l_names[0]) > 2
    for names in l_names[1:]:
        assert names == l_names[0]
