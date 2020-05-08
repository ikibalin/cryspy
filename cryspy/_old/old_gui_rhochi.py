import sys

import cryspy.scripts.rhochi.rhochi_viewer as rv

def main(l_arg= []):
    rv.main(l_arg)

if __name__ == '__main__':
    l_arg = sys.argv
    main(l_arg)