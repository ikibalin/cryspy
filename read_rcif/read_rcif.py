"""
transform inforamtion from cif file to class phase
"""
import sys


import os
"""
print sys.argv
fdir = os.path.split(os.path.dirname(os.path.abspath(sys.argv[0])))[0]
print "fdir \n\n\n", fdir
fdir_1 = os.path.join(fdir, "cwidgets")

sys.path.append(fdir_1)
"""
import ccore


class crcif(object):
    """
    class to store rcif data
    """

    def __init__(self):
        self.data = []

    def load_from_str(self, lcontent):
        p_glob = convert_rcif_to_glob(lcontent)
        self.data = p_glob["data"]

    def load_from_file(self, fname):
        fid = open(fname, "r")
        lcontent = fid.readlines()
        fid.close()
        self.load_from_str(lcontent)
        
    def save_to_str(self):
        lstr = rcif_to_str(self.data)
        return lstr

    def save_to_file(self, fname):
        lstr = self.save_to_str()
        fid = open(fname, "w")
        fid.write("\n".join(lstr))
        fid.close()

    
    def take_from_ccore(self, ccore):
        def temp_func(obj, ddata, l_rcif, l_core, t_core):
            lkey = obj.__dict__.keys()
            if l_core in lkey:
                if t_core == "logic":
                    sval = obj.__dict__[l_core]
                    ddata[l_rcif] = sval 
                elif t_core == "text":
                    ddata[l_rcif] = obj.__dict__[l_core]
                elif t_core == "val":
                    sval = obj.__dict__[l_core][0]
                    ddata[l_rcif] = sval 
            return
            
        drel = rcif_core_relation()
        lab_rcif = drel["lab_rcif_ref"]
        lab_core = drel["lab_core_ref"]
        type_core = drel["type_ref"]
        d_ref = {"name": "ref"}
        for l_rcif, l_core, t_core in zip(lab_rcif, lab_core, type_core):
            obj = ccore.ref[0]
            temp_func(obj, d_ref, l_rcif, l_core, t_core)
        
        l_par_ref = ccore.take_param_ref()
        l_par_con = ccore.take_param_const()
        if ((l_par_ref != [])|(l_par_con != [])):
            d_ref["loops"] = []
        if l_par_ref != []:
            l_link_rcif = [get_link_rcif(hh, ccore) for hh in l_par_ref]
            d_par_ref = {"_refinement_param1": l_link_rcif}
            d_ref["loops"].append(d_par_ref)
        if l_par_con != []:
            l_link_rcif_1 = [get_link_rcif(hh, ccore) for hh in l_par_con]

            d_par_con = {"_constraint_param1": l_link_rcif_1}
            
            l_link_rcif_2, l_coeff = [], []
            lhelp = [ccore.set_val_by_link(hh, None)[0][2] for hh in l_par_con]
            for hh in lhelp:
                coeff = float(hh.split("[")[0].strip().replace("x1","").replace("*",""))
                link_core_2 = hh.split("[")[1].strip()[:-1]
                link_rcif_2 = get_link_rcif(link_core_2, ccore)
                l_coeff.append(coeff)
                l_link_rcif_2.append(link_rcif_2)
                
            d_par_con["_constraint_param2"] = l_link_rcif_2
            d_par_con["_constraint_coeff"] = l_coeff

            d_ref["loops"].append(d_par_con)

        self.data.append(d_ref)
        

        
        for obj in ccore.ph:
            d_ph = {"name": obj.name}
            lab_rcif = drel["lab_rcif_sp"]
            lab_core = drel["lab_core_sp"]
            type_core = drel["type_sp"]
            for l_rcif, l_core, t_core in zip(lab_rcif, lab_core, type_core):
                temp_func(obj, d_ph, l_rcif, l_core, t_core)
            lab_rcif = drel["lab_rcif_ph"]
            lab_core = drel["lab_core_ph"]
            type_core = drel["type_ph"]
            for l_rcif, l_core, t_core in zip(lab_rcif, lab_core, type_core):
                temp_func(obj, d_ph, l_rcif, l_core, t_core)
                
            lloop_d = []
            
            d_at = {}
            lab_rcif = drel["lab_rcif_at"]
            lab_core = drel["lab_core_at"]
            type_core = drel["type_at"]
            lab_rcif_used = []
            l_dhelp = []
            for obj_2 in obj.atom:
                lkey = obj_2.__dict__.keys()
                dhelp = {}
                for l_rcif, l_core, t_core in zip(lab_rcif, lab_core, type_core):
                    temp_func(obj_2, dhelp, l_rcif, l_core, t_core)
                    if ((not(l_rcif in lab_rcif_used))&(l_core in lkey)):
                        lab_rcif_used.append(l_rcif)
                l_dhelp.append(dhelp) 
            for l_rcif in lab_rcif_used:
                d_at[l_rcif] = [hh[l_rcif] for hh in l_dhelp]
            lloop_d.append(d_at)
            
            d_bscat = {}
            lab_rcif = drel["lab_rcif_bscat"]
            lab_core = drel["lab_core_bscat"]
            type_core = drel["type_bscat"]
            lab_rcif_used = []
            l_dhelp = []
            for obj_2 in obj.atom:
                dhelp = {}
                lkey = obj_2.__dict__.keys()
                for l_rcif, l_core, t_core in zip(lab_rcif, lab_core, type_core):
                    temp_func(obj_2, dhelp, l_rcif, l_core, t_core)
                    if ((not(l_rcif in lab_rcif_used))&(l_core in lkey)):
                        lab_rcif_used.append(l_rcif)
                l_dhelp.append(dhelp) 
            lab_uniq = "_atom_site_type_symbol"
            l_val = [hh[lab_uniq] for hh in l_dhelp]
            s_uniq = set(l_val)
            lind_uniq = [l_val.index(hh) for hh in s_uniq]
            for l_rcif in lab_rcif_used:
                d_bscat[l_rcif] = [l_dhelp[ind][l_rcif] for ind in  lind_uniq]
            lloop_d.append(d_bscat)

            d_beta = {}
            lab_rcif = drel["lab_rcif_beta"]
            lab_core = drel["lab_core_beta"]
            type_core = drel["type_beta"]
            lab_rcif_used = []
            l_dhelp = []
            for obj_2 in obj.atom:
                lkey = obj_2.__dict__.keys()
                dhelp = {}
                for l_rcif, l_core, t_core in zip(lab_rcif, lab_core, type_core):
                    temp_func(obj_2, dhelp, l_rcif, l_core, t_core)
                    if ((not(l_rcif in lab_rcif_used))&(l_core in lkey)):
                        lab_rcif_used.append(l_rcif)
                l_dhelp.append(dhelp) 
            for l_rcif in lab_rcif_used:
                d_beta[l_rcif] = [hh[l_rcif] for hh in l_dhelp]
            lloop_d.append(d_beta)      
            
            d_chi = {}
            lab_rcif = drel["lab_rcif_chi"]
            lab_core = drel["lab_core_chi"]
            type_core = drel["type_chi"]
            lab_rcif_used = []
            l_dhelp = []
            for obj_2 in obj.atom:
                if obj_2.modemagn:
                    lkey = obj_2.__dict__.keys()
                    dhelp = {}
                    for l_rcif, l_core, t_core in zip(lab_rcif, lab_core, type_core):
                        temp_func(obj_2, dhelp, l_rcif, l_core, t_core)
                        if ((not(l_rcif in lab_rcif_used))&(l_core in lkey)):
                            lab_rcif_used.append(l_rcif)
                    l_dhelp.append(dhelp) 
            for l_rcif in lab_rcif_used:
                d_chi[l_rcif] = [hh[l_rcif] for hh in l_dhelp]
            lloop_d.append(d_chi)              
            
            if lloop_d != []:
                d_ph["loops"] = lloop_d
                
            self.data.append(d_ph)

        for obj in ccore.exp:
            d_exp = {"name": obj.name}
            if obj.powder:
                lab_rcif = drel["lab_rcif_exp_pd"]
                lab_core = drel["lab_core_exp_pd"]
                type_core = drel["type_exp_pd"]
            else:
                lab_rcif = drel["lab_rcif_exp_sd"]
                lab_core = drel["lab_core_exp_sd"]
                type_core = drel["type_exp_sd"]
                
            for l_rcif, l_core, t_core in zip(lab_rcif, lab_core, type_core):
                temp_func(obj, d_exp, l_rcif, l_core, t_core)
                
            d_eph = {}
            if obj.powder:
                lab_rcif = drel["lab_rcif_eph_pd"]
                lab_core = drel["lab_core_eph_pd"]
                type_core = drel["type_eph_pd"]
            else:
                lab_rcif = drel["lab_rcif_eph_sd"]
                lab_core = drel["lab_core_eph_sd"]
                type_core = drel["type_eph_sd"]
                
            lab_rcif_used = []
            l_dhelp = []
            for obj_2 in obj.exp_ph:
                lkey = obj_2.__dict__.keys()
                dhelp = {}
                for l_rcif, l_core, t_core in zip(lab_rcif, lab_core, type_core):
                    temp_func(obj_2, dhelp, l_rcif, l_core, t_core)
                    if ((not(l_rcif in lab_rcif_used))&(l_core in lkey)):
                        lab_rcif_used.append(l_rcif)
                l_dhelp.append(dhelp) 
            for l_rcif in lab_rcif_used:
                d_eph[l_rcif] = [hh[l_rcif] for hh in l_dhelp]
            lloop_d.append(d_eph)
                
            if lloop_d != []:
                d_exp["loops"] = lloop_d
                
            self.data.append(d_exp)

    def trans_to_ccore(self):
        res_type = []
        lphase = []
        lexp = []
        lref = []
        for data in self.data:
            lab_phase = "_cell"
            lab_exp_pd = "_pd"
            lab_exp_sd = "_sd"
            lab_ref = "_refinement"
            lkey = data.keys()
            flag_ph = any([hh.startswith(lab_phase) for hh in lkey])
            flag_exp_pd = any([hh.startswith(lab_exp_pd) for hh in lkey])
            flag_exp_sd = any([hh.startswith(lab_exp_sd) for hh in lkey])
            flag_ref = any([hh.startswith(lab_ref) for hh in lkey])
            if flag_ph:
                phase = trans_to_cphase(data)
                lphase.append(phase)
            elif flag_exp_pd:
                exp = trans_pd_to_cexp(data)
                lexp.append(exp)
            elif flag_exp_sd:
                exp = trans_sd_to_cexp(data)
                lexp.append(exp)
            elif flag_ref:
                ref, data_ref, data_con = trans_to_cref(data)
                lref.append(ref)

        core = ccore.ccore()
        core.ph = lphase
        core.exp = lexp
        core.ref = lref
        if data_ref != {}:
            lparam = data_ref["_refinement_param1"]
            put_ref(core, lparam)
        if data_con != {}:
            lparam1 = data_con["_constraint_param1"]
            lparam2 = data_con["_constraint_param2"]
            lcoeff = data_con["_constraint_coeff"]
            put_con(core, lparam1, lparam2, lcoeff)
        return core

def get_link_rcif(link_core, ccore):
    print ""
    link_core_2 = "".join([hh for hh in link_core.split()])[1:-1]
    
    lhelp = link_core_2.split(')(')
    type_1 = lhelp[0].split(",")[0]
    name_1 = lhelp[0].split(",")[1]
    link_rcif = "${:}".format(name_1)
    if len(lhelp) == 3:
        name_2 = lhelp[1].split(",")[1]
        link_rcif += "_${:}".format(name_2)
    name_core = lhelp[-1]
    
    drel = rcif_core_relation()
    if type_1 == "exp":
        lname_exp = [hh.name for hh in ccore.exp]
        ind = lname_exp.index(name_1)
        if ccore.exp[ind].powder:
            l_lab_rcif = ["lab_rcif_exp_pd", "lab_rcif_eph_pd"]
            l_lab_core = ["lab_core_exp_pd", "lab_core_eph_pd"]
        else:
            l_lab_rcif = ["lab_rcif_exp_sd", "lab_rcif_eph_sd"]
            l_lab_core = ["lab_core_exp_sd", "lab_core_eph_sd"]
    else:
        l_lab_rcif = ["lab_rcif_ph", "lab_rcif_sp", "lab_rcif_at", "lab_rcif_bscat", 
                  "lab_rcif_beta", "lab_rcif_chi"]
        l_lab_core = ["lab_core_ph", "lab_core_sp", "lab_core_at", "lab_core_bscat", 
                  "lab_core_beta", "lab_core_chi"]
    for lab_rcif, lab_core in zip(l_lab_rcif, l_lab_core):
        if name_core in drel[lab_core]:
            ind = drel[lab_core].index(name_core)
            name_rcif = drel[lab_rcif][ind]
            link_rcif += name_rcif 
            break
    return link_rcif

def get_link(obj, drel, param):
    lname_exp = [hh.name for hh in obj.exp]
    lname_ph = [hh.name for hh in obj.ph]
    lhelp = param.split("_")
    flag_1 = (lhelp[0][0] == "$")
    if flag_1:
        name = lhelp[0][1:]
    else:
        print "Mistake in link '{:}'".format(param)
        return ""
    flag_exp = (name in lname_exp)
    flag_ph = (name in lname_ph)
    slink = ""
    ind = -1
    if (flag_exp == (not flag_ph)):
        flag_2 = (lhelp[1][0] == "$")
        if flag_2:
            name_2 = lhelp[1][1:]
        if flag_exp:
            #experiment
            slink += "(exp,{:})".format(name)
            if flag_2:
                slink += "(exp_ph,{:})".format(name_2)
                param_s = "_"+"_".join(lhelp[2:])
                if param_s in drel["lab_rcif_eph_pd"]:
                    ind = drel["lab_rcif_eph_pd"].index(param_s)
                    slink += "({:})".format(drel["lab_core_eph_pd"][ind])
                elif param_s in drel["lab_rcif_eph_sd"]:
                    ind = drel["lab_rcif_eph_sd"].index(param_s)
                    slink += "({:})".format(drel["lab_core_eph_sd"][ind])
            else:
                
                param_s = "_"+"_".join(lhelp[1:])
                """
                print "param_s: ", param_s
                print drel["lab_rcif_exp"]
                """
                if param_s in drel["lab_rcif_exp_pd"]:
                    ind = drel["lab_rcif_exp_pd"].index(param_s)
                    slink += "({:})".format(drel["lab_core_exp_pd"][ind])
                elif param_s in drel["lab_rcif_exp_sd"]:
                    ind = drel["lab_rcif_exp_sd"].index(param_s)
                    slink += "({:})".format(drel["lab_core_exp_sd"][ind])
        else:
            #phase
            slink += "(ph,{:})".format(name)
            if flag_2:
                slink += "(atom,{:})".format(name_2)
                param_s = "_"+"_".join(lhelp[2:])
                if param_s in drel["lab_rcif_at"]:
                    ind = drel["lab_rcif_at"].index(param_s)
                    slink += "({:})".format(drel["lab_core_at"][ind])
                elif param_s in drel["lab_rcif_bscat"]:
                    ind = drel["lab_rcif_bscat"].index(param_s)
                    slink += "({:})".format(drel["lab_core_bscat"][ind])
                elif param_s in drel["lab_rcif_beta"]:
                    ind = drel["lab_rcif_beta"].index(param_s)
                    slink += "({:})".format(drel["lab_core_beta"][ind])
                elif param_s in drel["lab_rcif_chi"]:
                    ind = drel["lab_rcif_chi"].index(param_s)
                    slink += "({:})".format(drel["lab_core_chi"][ind])
            else:
                param_s = "_"+"_".join(lhelp[1:])
                if param_s in drel["lab_rcif_ph"]:
                    ind = drel["lab_rcif_ph"].index(param_s)
                    slink += "({:})".format(drel["lab_core_ph"][ind])
                elif param_s in drel["lab_rcif_sp"]:
                    ind = drel["lab_rcif_sp"].index(param_s)
                    slink += "({:})".format(drel["lab_core_sp"][ind])
    else:
        print "Mistake with naming of phase or experiment in link '{:}'".format(param)
        return ""
    if ind == -1:
        print "Code word in '{:}' is not defined".format(param)
        return ""
    return slink 


def rcif_to_str(drcif):
    drel = rcif_core_relation()


    def temp_func(dd2, llab):
        lstr = []
        lkey2 = dd2.keys()    
        llab_used = []
        for lab in llab:
            if lab in lkey2:
                llab_used.append(lab)
                line = "{:}".format(lab)
                lstr.append(line)
        nval = len(dd2[llab_used[0]])
        lline = [[] for hh in range(nval)]
        for lab in llab_used:
            for val, line in zip(dd2[lab], lline):
                if isinstance(val, str):
                    if len(val.split()) > 1:
                        line.append("'{:}'".format(val))
                    else:
                        line.append("{:}".format(val))
                else:
                    line.append(str(val))
        for line in lline:
            lstr.append(" "+"  ".join(line))
        return lstr
                            

    lstr = ["#\#RCIF_1.0", "global_filerhochi", ""]
    for dd in drcif:
        line = "data_{:}".format(dd["name"])
        lstr.append(line)
        
        lab_phase = "_cell"
        lab_exp_pd = "_pd"
        lab_exp_sd = "_sd"
        lab_ref = "_refinement"
        lkey = dd.keys()
        flag_ph = any([hh.startswith(lab_phase) for hh in lkey])
        flag_exp_pd = any([hh.startswith(lab_exp_pd) for hh in lkey])
        flag_exp_sd = any([hh.startswith(lab_exp_sd) for hh in lkey])
        flag_ref = any([hh.startswith(lab_ref) for hh in lkey])
        if flag_ph:
            llab = drel["lab_rcif_ph"]
            for lab in llab:
                if ((lab in lkey)&(lab != "name")):
                    if isinstance(dd[lab], str):
                        line = "{:} '{:}'".format(lab, dd[lab])
                    else:
                        line = "{:} {:}".format(lab, dd[lab])
                    lstr.append(line)
            lstr.append("\n")
            llab = drel["lab_rcif_sp"]
            for lab in llab:
                if lab in lkey:
                    if isinstance(dd[lab], str):
                        line = "{:} '{:}'".format(lab, dd[lab])
                    else:
                        line = "{:} {:}".format(lab, dd[lab])
                    lstr.append(line)
                    
            if "loops" in dd.keys():
                lab_at = "_atom_site_fract_x"
                lab_bscat = "_atom_site_bscat"
                lab_beta = "_atom_site_aniso_beta_11"
                lab_chi = "_atom_site_susceptibility_aniso_chi_11"
                for dd2 in dd["loops"]:
                    lstr.append("")  
                    lkey2 = dd2.keys()
                    flag_at = any([hh.startswith(lab_at) for hh in lkey2])
                    flag_bscat = any([hh.startswith(lab_bscat) for hh in lkey2])
                    flag_beta = any([hh.startswith(lab_beta) for hh in lkey2])
                    flag_chi = any([hh.startswith(lab_chi) for hh in lkey2])
                    if flag_at:
                        line = "loop_"
                        lstr.append(line)
                        llab = drel["lab_rcif_at"]
                        lstr_1 = temp_func(dd2, llab)
                        lstr.extend(lstr_1)
                            
                    elif flag_bscat:
                        line = "loop_"
                        lstr.append(line)
                        llab = drel["lab_rcif_bscat"]
                        lstr_1 = temp_func(dd2, llab)
                        lstr.extend(lstr_1)
                        
                    elif flag_beta:
                        line = "loop_"
                        lstr.append(line)
                        llab = drel["lab_rcif_beta"]
                        lstr_1 = temp_func(dd2, llab)
                        lstr.extend(lstr_1)

                    elif flag_chi:
                        line = "loop_"
                        lstr.append(line)
                        llab = drel["lab_rcif_chi"]
                        lstr_1 = temp_func(dd2, llab)
                        lstr.extend(lstr_1)
                    lstr.append("")                        
                    
        elif flag_exp_pd:
            llab = drel["lab_rcif_exp_pd"]
            for lab in llab:
                if ((lab in lkey)&(lab != "name")):
                    if isinstance(dd[lab], str):
                        line = "{:} '{:}'".format(lab, dd[lab])
                    else:
                        line = "{:} {:}".format(lab, dd[lab])
                    lstr.append(line)
                    
            if "loops" in dd.keys():
                lab_eph = "_pd_phase_scale"
                for dd2 in dd["loops"]:
                    lstr.append("")                      
                    lkey2 = dd2.keys()
                    flag_eph = any([hh.startswith(lab_eph) for hh in lkey2])
                    if flag_eph:
                        line = "loop_"
                        lstr.append(line)
                        llab = drel["lab_rcif_eph_pd"]
                        lstr_1 = temp_func(dd2, llab)
                        lstr.extend(lstr_1)
                    lstr.append("")                        
                    
        elif flag_exp_sd:
            llab = drel["lab_rcif_exp_sd"]
            for lab in llab:
                if ((lab in lkey)&(lab != "name")):
                    if isinstance(dd[lab], str):
                        line = "{:} '{:}'".format(lab, dd[lab])
                    else:
                        line = "{:} {:}".format(lab, dd[lab])
                    lstr.append(line)
                    
            if "loops" in dd.keys():
                lab_eph = "_sd_phase_name"
                for dd2 in dd["loops"]:
                    lkey2 = dd2.keys()
                    flag_eph = any([hh.startswith(lab_eph) for hh in lkey2])
                    if flag_eph:
                        lstr.append("")  
                        line = "loop_"
                        lstr.append(line)
                        llab = drel["lab_rcif_eph_sd"]
                        lstr_1 = temp_func(dd2, llab)
                        lstr.extend(lstr_1)
                        lstr.append("")                        
                    
        elif flag_ref:
            llab = drel["lab_rcif_ref"]
            for lab in llab:
                if lab in lkey:
                    if isinstance(dd[lab], str):
                        line = "{:} '{:}'".format(lab, dd[lab])
                    else:
                        line = "{:} {:}".format(lab, dd[lab])
                    lstr.append(line)
            if "loops" in dd.keys():
                lab_par_ref = "_refinement_param1"
                lab_par_con = "_constraint_param1"
                for dd2 in dd["loops"]:
                    lstr.append("")                      
                    lkey2 = dd2.keys()
                    flag_par_ref = any([hh.startswith(lab_par_ref) for hh in lkey2])
                    flag_par_con = any([hh.startswith(lab_par_con) for hh in lkey2])
                    if flag_par_ref:
                        line = "loop_"
                        lstr.append(line)
                        llab = [lab_par_ref]
                        lstr_1 = temp_func(dd2, llab)
                        lstr.extend(lstr_1)
                    elif flag_par_con:
                        line = "loop_"
                        lstr.append(line)
                        llab = ["_constraint_param1", "_constraint_param2", "_constraint_coeff"]
                        
                        lstr_1 = temp_func(dd2, llab)
                        lstr.extend(lstr_1)
                    lstr.append("")                        
        lstr.append("")

    return lstr

def put_ref(obj, lparam):
    drel = rcif_core_relation()
    for param in lparam:
        slink = get_link(obj, drel, param)
        val_1, message = obj.set_val_by_link(slink, None)
        val_1[1] = True
        obj.set_val_by_link(slink, val_1)


def put_con(obj, lparam1, lparam2, lcoeff):
    drel = rcif_core_relation()
    for param1, param2, coeff in zip(lparam1, lparam2, lcoeff):
        slink_obj = get_link(obj, drel, param1)
        slink_sub = get_link(obj, drel, param2)
        s_constr = " {:}*x1 [{:}]".format(coeff, slink_sub)
        val_1, message = obj.set_val_by_link(slink_obj, None)
        val_1[1] = False
        val_1[2] = s_constr
        obj.set_val_by_link(slink_obj, val_1)    


def rcif_core_relation():
    #phase    
    llab_rcif_ph = ["name", "_cell_length_a", "_cell_length_b", "_cell_length_c",
                 "_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma"]
    llab_core_ph = ["name", "cella", "cellb", "cellc",
                  "cellalpha", "cellbeta", "cellgama"]
    ltype_ph = ["text", "val", "val", "val", "val", "val", "val"]

    llab_rcif_sp = ["_space_group_name_H-M_alt", '_space_group_it_coordinate_system_code']
    llab_core_sp = ['spgr', 'spgr_code']
    ltype_sp = ["text", "text"]

    llab_rcif_at = ["_atom_site_label", '_atom_site_type_symbol',
                "_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z",
                "_atom_site_occupancy"]
    llab_core_at = ["name", "type_n", "coordx", "coordy", "coordz", "occ"]
    ltype_at = ["text", "text", "val", "val", "val", "val"]

    llab_rcif_bscat = ["_atom_site_type_symbol", "_atom_site_bscat"]
    llab_core_bscat = ["type_n", "bscat"]
    ltype_bscat = ["text", "val"]

    llab_rcif_beta = ["_atom_site_aniso_label", '_atom_site_aniso_beta_11',
             "_atom_site_aniso_beta_22", "_atom_site_aniso_beta_33", "_atom_site_aniso_beta_12",
             "_atom_site_aniso_beta_13", "_atom_site_aniso_beta_23"]
    llab_core_beta = ["name", "beta11", "beta22", "beta33", "beta12", "beta13", "beta23"]
    ltype_beta = ["text", "val", "val", "val", "val", "val", "val"]


    llab_rcif_chi = ["_atom_site_magnetism_aniso_label", '_atom_site_magnetism_type_symbol',
             "_atom_site_magnetism_kappa",
             "_atom_site_magnetism_lfactor", "_atom_site_susceptibility_aniso_chi_11",
             "_atom_site_susceptibility_aniso_chi_22", "_atom_site_susceptibility_aniso_chi_33",
             "_atom_site_susceptibility_aniso_chi_23", "_atom_site_susceptibility_aniso_chi_13",
             "_atom_site_susceptibility_aniso_chi_12"]
    llab_core_chi = ["name", "type", "kappa", "lfactor", "chi11", "chi22", "chi33", "chi12", "chi13", "chi23"]
    ltype_chi = ["text", "text", "val", "val", "val", "val", "val", "val", "val", "val"]


    #experiment
    llab_rcif_exp_pd = ["name", "_pd_file_name_bkgr", "_pd_file_name_input", "_pd_file_name_output",
                 "_pd_chi2_up", "_pd_chi2_down", "_pd_chi2_diff", "_pd_exl2theta",
                 "_pd_beam_polarization_up", "_pd_beam_polarization_down", "_pd_resolution_u",
                 "_pd_resolution_v", "_pd_resolution_w", "_pd_resolution_x", "_pd_resolution_y",
                 "_pd_shift_const", "_pd_shift_coeff",
                 "_pd_reflex_asymetry_p1", "_pd_reflex_asymetry_p2", "_pd_reflex_asymetry_p3",
                 "_pd_reflex_asymetry_p4"]
    llab_core_exp_pd = ["name", "bkgr", "input", "output", "modechi2_up", "modechi2_down",
                "modechi2_diff", "excl2theta", "pup", "pdown", "U", "W", "V", "x", "y",
                "zeroshift", "zshift_a", "p_asym_1", "p_asym_2", "p_asym_3", "p_asym_4"]
    ltype_exp_pd = ["text", "text", "text", "text", "logic", "logic", "logic", "text", "val", "val",
             "val", "val", "val", "val", "val", "val", "val", "val", "val", "val", "val"]


    llab_rcif_eph_pd = ["_pd_phase_name", "_pd_phase_scale", "_pd_phase_igsize",
                     "_pd_phase_extinction_radius", "_pd_phase_extinction_mosaicity"]
    llab_core_eph_pd = ["name", "scale", "igsize", "ext11", "ext22"]
    ltype_eph_pd = ["text", "val", "val", "val", "val"]

    llab_rcif_exp_sd = ["name", "_sd_file_name_input", "_sd_file_name_output",
                 "_sd_beam_polarization_up", "_sd_beam_polarization_down"]
    llab_core_exp_sd = ["name", "input", "output", "pup", "pdown"]
    ltype_exp_sd = ["text", "text", "text", "val", "val"]

    llab_rcif_eph_sd = ["_sd_phase_name", "_sd_phase_extinction_radius", "_sd_phase_extinction_mosaicity"]
    llab_core_eph_sd = ["name", "ext11", "ext22"]
    ltype_eph_sd = ["text", "val", "val"]



    #refinement
    llab_rcif_ref = ["_refinement_file_name_output", "_refinement_flag_errorbars",
                 "_refinement_flag_refinement"]
    llab_core_ref = ["output", "refin", "sigmas"]
    ltype_ref = ["text", "logic", "logic"]

    drel = {}

    drel["lab_rcif_ph"] = llab_rcif_ph
    drel["lab_core_ph"] = llab_core_ph
    drel["type_ph"] = ltype_ph
    drel["lab_rcif_sp"] = llab_rcif_sp
    drel["lab_core_sp"] = llab_core_sp
    drel["type_sp"] = ltype_sp
    drel["lab_rcif_at"] = llab_rcif_at
    drel["lab_core_at"] = llab_core_at
    drel["type_at"] = ltype_at
    drel["lab_rcif_bscat"] = llab_rcif_bscat
    drel["lab_core_bscat"] = llab_core_bscat
    drel["type_bscat"] = ltype_bscat
    drel["lab_rcif_beta"] = llab_rcif_beta
    drel["lab_core_beta"] = llab_core_beta
    drel["type_beta"] = ltype_beta
    drel["lab_rcif_chi"] = llab_rcif_chi
    drel["lab_core_chi"] = llab_core_chi
    drel["type_chi"] = ltype_chi
    drel["lab_rcif_exp_pd"] = llab_rcif_exp_pd
    drel["lab_core_exp_pd"] = llab_core_exp_pd
    drel["type_exp_pd"] = ltype_exp_pd
    drel["lab_rcif_eph_pd"] = llab_rcif_eph_pd
    drel["lab_core_eph_pd"] = llab_core_eph_pd
    drel["type_eph_pd"] = ltype_eph_pd
    drel["lab_rcif_exp_sd"] = llab_rcif_exp_sd
    drel["lab_core_exp_sd"] = llab_core_exp_sd
    drel["type_exp_sd"] = ltype_exp_sd
    drel["lab_rcif_eph_sd"] = llab_rcif_eph_sd
    drel["lab_core_eph_sd"] = llab_core_eph_sd
    drel["type_eph_sd"] = ltype_eph_sd
    drel["lab_rcif_ref"] = llab_rcif_ref
    drel["lab_core_ref"] = llab_core_ref
    drel["type_ref"] = ltype_ref
    
    return drel


def trans_to_cphase(data):
    phase = ccore.ccore_ph()
    drel = rcif_core_relation()
    
    llab_rcif = drel["lab_rcif_ph"]
    llab_phase = drel["lab_core_ph"]
    ltype = drel["type_ph"]

    from_dict_to_obj(data, llab_rcif, phase, llab_phase, ltype)

    llab_rcif = drel["lab_rcif_sp"]
    llab_phase = drel["lab_core_sp"]
    ltype = drel["type_sp"]
    from_dict_to_obj(data, llab_rcif, phase, llab_phase, ltype)

    lkey = data.keys()
    if (not ("loops" in lkey)):
        return phase

    lab_xyz = "_atom_site_fract"
    lab_beta = "_atom_site_aniso_beta"
    lab_chi = "_atom_site_susceptibility"
    lab_bscat = "_atom_site_bscat"
    lab_j0 = "_atom_site_magnetism_j0"
    lab_j2 = "_atom_site_magnetism_j2"
    llabs = [lab_xyz, lab_beta, lab_chi, lab_bscat, lab_j0, lab_j2]
    lnumb = [[] for hh in llabs]

    for iloop, loop in enumerate(data["loops"]):
        lkey = loop.keys()
        for hh, numb in zip(llabs, lnumb):
            flag = any([True if hh1.startswith(hh) else False for hh1 in lkey])
            if flag:
                numb.append(iloop)

    obj_lat_xyz = ccore.ccore_ph_at()
    obj_lat_beta = ccore.ccore_ph_at()
    obj_lat_chi = ccore.ccore_ph_at()
    obj_lat_bscat = ccore.ccore_ph_at()
    if lnumb[0] != []:
        llab_rcif = drel["lab_rcif_at"]
        llab_phase = drel["lab_core_at"]
        ltype = len(llab_rcif) * ["list"]
        for hh in llab_phase:
            obj_lat_xyz.__dict__[hh] = []
        for numb in lnumb[0]:
            from_dict_to_obj(data["loops"][numb], llab_rcif, obj_lat_xyz, llab_phase, ltype)
    if lnumb[1] != []:
        llab_rcif = drel["lab_rcif_beta"]
        llab_phase = drel["lab_core_beta"]
        ltype = len(llab_rcif) * ["list"]
        for hh in llab_phase:
            obj_lat_beta.__dict__[hh] = []
        for numb in lnumb[1]:
            from_dict_to_obj(data["loops"][numb], llab_rcif, obj_lat_beta, llab_phase, ltype)
    if lnumb[2] != []:
        llab_rcif = drel["lab_rcif_chi"]
        llab_phase = drel["lab_core_chi"]
        ltype = len(llab_rcif) * ["list"]
        for hh in llab_phase:
            obj_lat_chi.__dict__[hh] = []
        for numb in lnumb[2]:
            from_dict_to_obj(data["loops"][numb], llab_rcif, obj_lat_chi, llab_phase, ltype)
    if lnumb[3] != []:
        llab_rcif = drel["lab_rcif_bscat"]
        llab_phase = drel["lab_core_bscat"]
        ltype = len(llab_rcif) * ["list"]
        for hh in llab_phase:
            obj_lat_bscat.__dict__[hh] = []
        for numb in lnumb[3]:
            from_dict_to_obj(data["loops"][numb], llab_rcif, obj_lat_bscat, llab_phase, ltype)
    if lnumb[4] != []:
        pass
    if lnumb[5] != []:
        pass

    nat = len(obj_lat_xyz.name)
    if nat == 0:
        return phase
    latom = []
    for name in obj_lat_xyz.name:
        atom = ccore.ccore_ph_at()
        atom.name = name
        latom.append(atom)
    llab_phase = drel["lab_core_at"]
    ltype = drel["type_at"]
    for lab, typ in zip(llab_phase, ltype):
        for hh, atom in zip(obj_lat_xyz.__dict__[lab], latom):
            if typ == "text":
                atom.__dict__[lab] = hh
            elif typ == "val":
                atom.__dict__[lab] = [hh, False, ""]

    if len(obj_lat_beta.name) != 0:
        for atom in latom:
            if atom.name in obj_lat_beta.name:
                ind = obj_lat_beta.name.index(atom.name)
                llab_phase = drel["lab_core_beta"]
                ltype = drel["type_beta"]
                for lab, typ in zip(llab_phase, ltype):
                    if typ == "text":
                        atom.__dict__[lab] = obj_lat_beta.__dict__[lab][ind]
                    elif typ == "val":
                        atom.__dict__[lab] = [obj_lat_beta.__dict__[lab][ind], False, ""]

    if len(obj_lat_chi.name) != 0:
        for atom in latom:
            if atom.name in obj_lat_chi.name:
                atom.modemagn = True
                ind = obj_lat_chi.name.index(atom.name)
                llab_phase = drel["lab_core_chi"]
                ltype = drel["type_chi"]
                for lab, typ in zip(llab_phase, ltype):
                    if typ == "text":
                        atom.__dict__[lab] = obj_lat_chi.__dict__[lab][ind]
                    elif typ == "val":
                        atom.__dict__[lab] = [obj_lat_chi.__dict__[lab][ind], False, ""]

    if len(obj_lat_bscat.type_n) != 0:
        for atom in latom:
            if atom.type_n in obj_lat_bscat.type_n:
                ind = obj_lat_bscat.type_n.index(atom.type_n)
                llab_phase = drel["lab_core_bscat"]
                ltype = drel["type_bscat"]
                for lab, typ in zip(llab_phase, ltype):
                    if typ == "text":
                        atom.__dict__[lab] = obj_lat_bscat.__dict__[lab][ind]
                    elif typ == "val":
                        atom.__dict__[lab] = [obj_lat_bscat.__dict__[lab][ind], False, ""]
    phase.atom = latom
    return phase


def from_dict_to_obj(dict_i, llab_d, obj, llab_o, ltype):
    lkey = dict_i.keys()
    lnumb = [ihh for ihh, hh in enumerate(llab_d) if (hh in lkey)]
    for numb in lnumb:
        if ltype[numb] == "val":
            val = [dict_i[llab_d[numb]], False, ""]
        elif ltype[numb] == "text":
            val = dict_i[llab_d[numb]]
        elif ltype[numb] == "list":
            val = dict_i[llab_d[numb]]
        elif ltype[numb] == "logic":
            val = dict_i[llab_d[numb]]
        else:
            print "mistake in type variable of 'from_dict_to_obj' function"
            val = None
        obj.__dict__[llab_o[numb]] = val
    return


def trans_pd_to_cexp(data):
    exp = ccore.ccore_exp()
    exp.powder = True
    drel = rcif_core_relation()
    
    llab_rcif = drel["lab_rcif_exp_pd"]
    llab_exp = drel["lab_core_exp_pd"]
    ltype = drel["type_exp_pd"]

    from_dict_to_obj(data, llab_rcif, exp, llab_exp, ltype)

    lkey = data.keys()
    if (not ("loops" in lkey)):
        return exp

    lab_exp = "_pd_phase_name"

    obj_ph = ccore.ccore_exp_ph()
    llab_rcif = drel["lab_rcif_eph_pd"]
    llab_exp = drel["lab_core_eph_pd"]
    ltype = len(llab_rcif) * ["list"]
    for hh in llab_exp:
        obj_ph.__dict__[hh] = []
    for data_ph in data["loops"]:
        lkey = data_ph.keys()
        if lab_exp in lkey:
            from_dict_to_obj(data_ph, llab_rcif, obj_ph, llab_exp, ltype)
            break
    if len(obj_ph.name) != 0:
        lph_exp = []
        for name in obj_ph.name:
            ph_exp = ccore.ccore_exp_ph()
            lph_exp.append(ph_exp)
        llab_exp = drel["lab_core_eph_pd"]
        ltype = drel["type_eph_pd"]
        for lab, typ in zip(llab_exp, ltype):
            for hh, ph_exp in zip(obj_ph.__dict__[lab], lph_exp):
                if typ == "text":
                    ph_exp.__dict__[lab] = hh
                elif typ == "val":
                    ph_exp.__dict__[lab] = [hh, False, ""]
        exp.exp_ph = lph_exp
    return exp

def trans_sd_to_cexp(data):
    exp = ccore.ccore_exp()
    exp.powder = False
    drel = rcif_core_relation()
    
    llab_rcif = drel["lab_rcif_exp_sd"]
    llab_exp = drel["lab_core_exp_sd"]
    ltype = drel["type_exp_sd"]

    from_dict_to_obj(data, llab_rcif, exp, llab_exp, ltype)

    lkey = data.keys()
    if (not ("loops" in lkey)):
        return exp

    lab_exp = "_sd_phase_name"

    obj_ph = ccore.ccore_exp_ph()
    llab_rcif = drel["lab_rcif_eph_sd"]
    llab_exp = drel["lab_core_eph_sd"]
    ltype = len(llab_rcif) * ["list"]
    for hh in llab_exp:
        obj_ph.__dict__[hh] = []
    for data_ph in data["loops"]:
        lkey = data_ph.keys()
        if lab_exp in lkey:
            from_dict_to_obj(data_ph, llab_rcif, obj_ph, llab_exp, ltype)
            break
    if len(obj_ph.name) != 0:
        lph_exp = []
        for name in obj_ph.name:
            ph_exp = ccore.ccore_exp_ph()
            lph_exp.append(ph_exp)
        llab_exp = drel["lab_core_eph_sd"]
        ltype = drel["type_eph_sd"]
        for lab, typ in zip(llab_exp, ltype):
            for hh, ph_exp in zip(obj_ph.__dict__[lab], lph_exp):
                if typ == "text":
                    ph_exp.__dict__[lab] = hh
                elif typ == "val":
                    ph_exp.__dict__[lab] = [hh, False, ""]
        exp.exp_ph = lph_exp
    return exp

def trans_to_cref(data):
    ref = ccore.ccore_ref()
    drel = rcif_core_relation()
    
    llab_rcif = drel["lab_rcif_ref"]
    llab_ref = drel["lab_core_ref"]
    ltype = drel["type_ref"]
    
    from_dict_to_obj(data, llab_rcif, ref, llab_ref, ltype)
    

    lkey = data.keys()
    data_ref, data_con = {}, {}
    if (not ("loops" in lkey)):
        return ref, data_ref, data_con
    lab_ref = "_refinement_param1"
    lab_con = "_constraint_param1"
    for data in data["loops"]:
        if lab_ref in data.keys():
            data_ref = {}
            data_ref[lab_ref] = data[lab_ref]
        if lab_con in data.keys():
            data_con = {}
            data_con[lab_con] = data[lab_con]
            data_con["_constraint_param2"] = data["_constraint_param2"]
            data_con["_constraint_coeff"] = data["_constraint_coeff"]
    return ref, data_ref, data_con


def smart_spleet(str):
    """
    split string like:
    "C C 0.0033 0.0016 'International Tables Vol C Tables 4.2.6.8 and 6.1.1.4'"
    in the list like:
    ['C', 'C', '0.0033', '0.0016', 'International Tables Vol C Tables 4.2.6.8 and 6.1.1.4']
    """
    flag_in = False
    lval, val = [], []
    for hh in str.strip():
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
    lval_2 = []
    for val in lval:
        
        val_2 = conv_str_to_text_float_logic(val)

        lval_2.append(val_2)
    return lval_2


def smart_numb(str):
    """
    convert string of the format '5.7(3)' into number 5.7
    """
    ind = str.find("(")
    if ind == -1:
        try:
            numb = float(str)
        except:
            numb = 0.
    else:
        numb = float(str[:ind])
    return numb


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

    lcont_vals = [hh for hh, flag in zip(lcont, lflag_loops) if not flag]

    del lflag_loops
    lval_numb = [iline for iline, line in enumerate(lcont_vals) if line.startswith("_")]
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
        value = line[len(name):]
        if len(lnumb) > 1:
            value += " ".join([lcont_vals[numb] for numb in lnumb[1:]])
        res[name] = conv_str_to_text_float_logic(value)
    return res

def conv_str_to_text_float_logic(sval):
    if (not (isinstance(sval,str))):
        return sval
    try:
        val = float(sval)
        return val
    except:
        val_1 = sval.strip()
        if len(val_1) == 0:
            val = None
        if (((val_1[0]=="\"")|(val_1[0]=="'"))&((val_1[-1]=="\"")|(val_1[-1]=="'"))):
            val_1=val_1[1:-1]
        if len(val_1) == 0:
            val = None
        elif val_1.lower() == "true":
            val = True
        elif val_1.lower() == "false":
            val = False
        else:
            val = val_1
    return val

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
        lval = smart_spleet(line)
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

def get_core_from_file(fname):
    rcif = crcif()
    rcif.load_from_file(fname)
    core = rcif.trans_to_ccore()
    return core
    
def save_core_to_file(core, fname):
    rcif_2 = crcif()
    rcif_2.take_from_ccore(core)
    rcif_2.save_to_file(fname)
    return 

def main(larg):
    if len(larg) > 1:
        fname = "powder.rcif"
        fname = larg[1]
    else:
        fname = raw_input("What is the name of the 'rcif' file?\n")
        
    rcif = crcif()
    rcif.load_from_file(fname)

    core = rcif.trans_to_ccore()
    
    rcif_2 = crcif()
    rcif_2.take_from_ccore(core)
    rcif_2.save_to_file("powder_2.rcif")
    
    

if __name__ == "__main__":
    main(sys.argv)
