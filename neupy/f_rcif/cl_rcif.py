"""
transform inforamtion from cif file to class phase
"""
import sys

import os

class RCif(object):
    """
    class to store rcif data
    """

    def __init__(self):
        self.glob = {"data":[], "name":""}
        self._p_file_dir = ""
        self._p_file_name = ""

    def __repr__(self):
        ls_out = self.save_to_str()
        return "\n".join(ls_out)

    def load_from_str(self, lcontent):
        p_glob = convert_rcif_to_glob(lcontent)
        #self.data = p_glob["data"]
        self.glob = p_glob

    def load_from_file(self, fname):
        self._p_file_dir = os.path.dirname(fname)
        self._p_file_name = os.path.basename(fname)
        
        fid = open(fname, "r")
        lcontent = fid.readlines()
        fid.close()
        self.load_from_str(lcontent)
        
    def save_to_str(self):
        lstr = rcif_to_str(self.glob)
        return lstr

    def save_to_file(self, f_name):
        lstr = self.save_to_str()
        fid = open(f_name, "w")
        fid.write("\n".join(lstr)+"\n")
        fid.close()


def rcif_to_str(p_glob):
    #structure of cif is data blokcs and loop blocks, 
    #nested loops are not supported
    ls_out = []
    l_key_g = p_glob.keys()
    if "name" in l_key_g:
        ls_out.append("global_{:}\n".format(p_glob["name"].strip()))
    else:
        ls_out.append("global_\n")
    for lab in sorted(l_key_g):
        if lab.startswith("_"):
            if ((len(p_glob[lab].strip().split())>1)&
                 (not(p_glob[lab].strip().startswith("'")))):
                ls_out.append("{:} '{:}'".format(lab, p_glob[lab].strip()))
            else:
                ls_out.append("{:} {:}".format(lab, p_glob[lab].strip()))
    
    l_data = p_glob["data"]
    for ddata in l_data:    
        l_key = ddata.keys()    
        if "name" in l_key:
            ls_out.append("\ndata_{:}".format(ddata["name"].strip()))
        else:
            ls_out.append("\ndata_")
            
        for lab in sorted(l_key):
            if lab.startswith("_"):
                if ((len(ddata[lab].strip().split())>1)&
                    (not(ddata[lab].strip().startswith("'")))):
                    ls_out.append("{:} '{:}'".format(lab, ddata[lab].strip()))
                else:
                    ls_out.append("{:} {:}".format(lab, ddata[lab].strip()))

        if "loops" in l_key:
            for d_loop in ddata["loops"]:
                l_key_loop = d_loop.keys()
                if d_loop != {}:
                    if "name" in sorted(l_key_loop):
                        ls_out.append("\nloop_{:}".format(d_loop["name"].strip()))
                    else:
                        ls_out.append("\nloop_")

                    for lab in sorted(l_key_loop):
                        ls_out.append(lab)

                    n_line = len(d_loop[lab])
                    for i_line in range(n_line):
                        line = []
                        for lab in sorted(l_key_loop):
                            line.append("{:}".format(d_loop[lab][i_line].strip()))
                        ls_out.append(" "+" ".join(line))
    return ls_out





def from_dict_to_obj(dict_i, llab_d, obj, llab_o, ltype):
    lkey = dict_i.keys()
    lnumb = [ihh for ihh, hh in enumerate(llab_d) if (hh in lkey)]
    d_args = {}
    for numb in lnumb:
        if ltype[numb] == "val":
            #val = [dict_i[llab_d[numb]], False, ""]
            val = dict_i[llab_d[numb]]
            val = val
        elif ltype[numb] == "text":
            val = dict_i[llab_d[numb]]
        elif ltype[numb] == "list":
            val = dict_i[llab_d[numb]]
        elif ltype[numb] == "logic":
            val = dict_i[llab_d[numb]]
        else:
            print("mistake in type variable of 'from_dict_to_obj' function")
            val = None
        d_args.update({llab_o[numb]:val})
    obj.set_val(**d_args)
    return





def smart_spleet(line, l_name):
    """
    split string like:
    "C C 0.0033 0.0016 'International Tables Vol C Tables 4.2.6.8 and 6.1.1.4'"
    in the list like:
    ['C', 'C', '0.0033', '0.0016', 'International Tables Vol C Tables 4.2.6.8 and 6.1.1.4']
    """
    flag_in = False
    lval, val = [], []
    for hh in line.strip():
        if (hh == " ") & (not flag_in):
            if val != []:
                lval.append("".join(val))
            val = []
        elif (hh == " ") & (flag_in):
            val.append(hh)
        elif hh == "'":
            flag_in = not flag_in
        else:
            val.append(hh)
    if val != []:
        lval.append("".join(val))
    return lval




def lines_in_block(label, lcontent):
    lnumb_p = [ihh for ihh, hh in enumerate(lcontent) if hh.startswith(label)]
    if len(lnumb_p) > 1:
        lnumb_b = lnumb_p
        lnumb_e = [hh for hh in lnumb_p[1:]]
        lnumb_e.append(len(lcontent))
    elif len(lnumb_p) == 1:
        lnumb_b = lnumb_p
        lnumb_e = [len(lcontent)]
    else:
        return [], [], []
    lname = [lcontent[numb_b][len(label):].strip() for numb_b in lnumb_b]
    lcont = [lcontent[(numb_b + 1):numb_e] for numb_b, numb_e in zip(lnumb_b, lnumb_e)]
    return lcont, lname, lnumb_b


def convert_lines_to_global(lcont):
    res = {"data": []}
    ldata_lstr, ldata_name, lnumb_b = lines_in_block("data_", lcont)
    if len(lnumb_b) != 0:
        val_lstr = lcont[:lnumb_b[0]]
        res_val = convert_lines_to_vals(val_lstr)
        res.update(res_val)
    for data_lstr, data_name in zip(ldata_lstr, ldata_name):
        res_data = convert_lines_to_vals(data_lstr)
        res_data["name"] = data_name
        res["data"].append(res_data)

    return res


def convert_lines_to_vals(lcont):
    res = {}
    lflag_loops = []
    flag, flag_trigger = False, False
    for line in lcont:
        if line.startswith("loop_"):
            flag = True
            flag_trigger = False
            lflag_loops.append(flag)
        elif line.startswith("_"):
            if flag_trigger:
                flag = False
                flag_trigger = False
            lflag_loops.append(flag)
        else:
            lflag_loops.append(flag)
            if flag:
                flag_trigger = True
    del flag, flag_trigger
    lcont_loops = [hh for hh, flag in zip(lcont, lflag_loops) if flag]

    lloop_lstr, lloop_name, lloop_numb = lines_in_block("loop_", lcont_loops)

    lres_loops = []
    for loop_lstr in lloop_lstr:
        res_loops = convert_lines_to_loop(loop_lstr)
        lres_loops.append(res_loops)
    if len(lres_loops) > 0:
        res["loops"] = lres_loops

    lcont_vals = [hh.strip() for hh, flag in zip(lcont, lflag_loops) if not flag]

    del lflag_loops
    lval_numb = [iline for iline, line in enumerate(lcont_vals) if line.startswith("_")]
    lcommon = []
    if len(lval_numb) > 1:
        lcommon = [range(inumb, lval_numb[ival + 1]) for ival, inumb in enumerate(lval_numb[:-1])]
    if len(lval_numb) > 0:
        lcommon.append(range(lval_numb[-1], len(lcont_vals)))
    else:
        return res

    for lnumb in lcommon:
        numb1 = lnumb[0]
        line = lcont_vals[numb1]
        name = line.split()[0]
        value = line[len(name):].strip()
        if len(lnumb) > 1:
            value += " ".join([lcont_vals[numb] for numb in lnumb[1:]])
        res[name] = value
    return res


def convert_lines_to_loop(lcont):
    lname = []
    for iline, line in enumerate(lcont):
        if line.startswith("_"):
            lname.append(line)
        else:
            break
    res = {}
    for name in lname:
        res[name] = []
    for line in lcont[iline:]:
        lval = smart_spleet(line, lname)
        for name, val in zip(lname, lval):
            res[name].append(val)
    return res


def convert_rcif_to_glob(lcontent):
    lcont = [hh[:hh.find("#")] if hh.find("#") != -1 else hh for hh in lcontent]
    lcont = [hh.strip() for hh in lcont if hh.strip() != ""]

    lglob_lstr, lglob_name, lnumb = lines_in_block("global_", lcont)
    if len(lglob_lstr) == 0:
        glob_lstr = lcont
        glob_name = ""
    else:
        glob_lstr = lglob_lstr[0]
        glob_name = lglob_name[0]

    res_glob = convert_lines_to_global(glob_lstr)
    res_glob["name"] = glob_name
    return res_glob

    

def main(larg):
    if len(larg) > 1:
        fname = "powder.rcif"
        fname = larg[1]
    else:
        fname = input("What is the name of the 'rcif' file?\n")
        
    rcif = RCif()
    rcif.load_from_file(fname)
    print(rcif)
    

if __name__ == "__main__":
    main(sys.argv)
