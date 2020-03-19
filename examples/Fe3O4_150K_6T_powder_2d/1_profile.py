import os
import sys
import pstats
#import cProfile

def main(l_arg):
    l_func = []
    s_dir = os.path.dirname(l_arg[0])
    if len(l_arg) >= 2:
        s_py = l_arg[1]
        if os.path.isfile(s_py):
            s_name = f"python -m cProfile -o out.prof -s time {s_py:} 5"
            #cProfile.run("_obj.refine()")
            os.system(s_name)
            if len(l_arg) >= 3:
                l_func = l_arg[2:]
        else:
            l_func = l_arg[1:]
        
    p = pstats.Stats("out.prof")
    p.strip_dirs()
    if len(l_func) != 0:
        for _func in l_func:
            p.sort_stats("time").print_callees(_func)
            p.sort_stats("time").print_callers(_func)
    p.sort_stats("time").print_stats(10)
    p.sort_stats("cumtime").print_stats(30)

main(sys.argv)

