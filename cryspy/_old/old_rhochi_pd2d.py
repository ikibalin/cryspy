import os
import os.path

import sys

from cryspy.scripts.rhochi.cl_rhochi import  create_temporary

def main(l_arg= []):
    f_name_in = os.path.join(os.getcwd(), "main.rcif")
    create_temporary(f_name_in=f_name_in, exp_type="3")

if __name__ == '__main__':
    l_arg = sys.argv
    main(l_arg)