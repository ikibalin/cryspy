import os
import os.path

import sys

from cryspy.scripts.rhochi.cl_rhochi import  create_temporary

def main(l_arg= []):
    f_name_in = os.path.join(os.getcwd(), "main.rcif")
    flag = True
    if os.path.isfile(f_name_in):
        print(f"The file: \n\n'{f_name_in:}'\n\nis already exist.")
        answ = input("\nDo you want to rewrite it? (yes or no)\n").strip().lower()
        if not("yes" in answ):
            flag = False
    if flag:
        create_temporary(f_name_in=f_name_in, exp_type="4")

if __name__ == '__main__':
    l_arg = sys.argv
    main(l_arg)