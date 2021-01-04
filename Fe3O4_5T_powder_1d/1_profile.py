import os
import sys
import pstats
#import cProfile

def main(l_arg):
    s_dir = os.path.dirname(l_arg[0])
    s_prof = "out.prof"
    _i_func = 1
    flag_out = False
    if len(l_arg) >= _i_func+1:
        if (os.path.isfile(l_arg[_i_func]) & l_arg[_i_func].endswith(".py")):
            s_py = l_arg[_i_func]
            _i_func += 1
            flag_out = True
    if len(l_arg) >= _i_func + 1:
        if l_arg[_i_func].endswith(".prof"):
            s_prof = l_arg[_i_func]
            _i_func += 1
    if len(l_arg) >= _i_func+1:
        l_func = l_arg[_i_func:]
    else:
        l_func = []

    if flag_out:
        s_name = f"python -m cProfile -o {s_prof:} -s time {s_py:} 5"
        # cProfile.run("_obj.refine()")
        os.system(s_name)

    p = pstats.Stats(s_prof)
    p.strip_dirs()
    if len(l_func) != 0:
        for _func in l_func:
            p.sort_stats("time").print_callees(_func)
            p.sort_stats("time").print_callers(_func)
    p.sort_stats("time").print_stats(10)
    p.sort_stats("cumtime").print_stats(30)

main(sys.argv)

