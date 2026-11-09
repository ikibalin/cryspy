# -*- coding: utf-8 -*-
"""Console mode of cryspy library.

python -m cryspy [file.rcif]

Loads the file (main.rcif by default), runs the RhoChi refinement when
the file has refined parameters and writes the result back to the file.
"""
import os
import sys

from cryspy.H_functions_global.function_1_cryspy_objects import \
    file_to_globaln
from cryspy.procedure_rhochi.rhochi import rhochi_rietveld_refinement


def main(l_arg=None) -> int:
    if l_arg is None:
        l_arg = sys.argv[1:]
    print("Console mode of cryspy library")
    f_name = l_arg[0] if len(l_arg) > 0 else "main.rcif"
    if not os.path.isfile(f_name):
        print(f"File '{f_name}' is not found")
        return 1
    obj_global = file_to_globaln(f_name)
    if len(obj_global.get_variable_names()) == 0:
        print("The file has no refined parameters.")
        return 0
    print("Refinement is running")
    res = rhochi_rietveld_refinement(obj_global)
    print("Result of refinement: \n ", res)
    with open(f_name, "w") as fid:
        fid.write(obj_global.to_cif())
    print("The file is rewritten")
    return 0


if __name__ == "__main__":
    sys.exit(main())
