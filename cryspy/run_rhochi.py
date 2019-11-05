import sys

from cryspy import rhochi_refinement

def main(l_arg= []):
    rhochi_refinement(f_name_in=None, f_name_out=None)

if __name__ == '__main__':
    l_arg = sys.argv
    main(l_arg)