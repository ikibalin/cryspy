# -*- coding: utf-8 -*-
"""Console mode of cryspy library."""
import os
import sys

from cryspy.H_functions_global.function_1_cryspy_objects import \
    file_to_globaln, RhoChi, MEM

if __name__ == "__main__":
    print("Console mode of cryspy library")
    l_arg = sys.argv
    if len(l_arg) > 1:
        f_name = l_arg[1]
        if not os.path.isfile(f_name):
            f_name = "main.rcif"
    if not os.path.isfile(f_name):
        print("File is not found")
    else:
        obj_global = file_to_globaln(f_name)
        if isinstance(obj_global, RhoChi):
            print("RhoChi object is detected")
            print("Refienement is running")
            res = obj_global.refine()
            print("Result of refinement: \n ", res)
            obj_global.save_to_file(f_name)
            print("The file is rewritten")
        elif isinstance(obj_global, MEM):
            print("MEM object is detected")
            print("The action algorithm is not written.")
        else:
            print("The found file is not defined as any predescribed class.")
