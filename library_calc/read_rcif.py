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
from crystal import *
from calculated_data import *
from setup_powder_1d import *
from setup_powder_2d import *
from observed_data import *
from experiment import *
from model import *
from variable import *


class RCif(object):
    """
    class to store rcif data
    """

    def __init__(self):
        self.data = []
        self._p_file_dir = ""
        self._p_file_name = ""

    def load_from_str(self, lcontent):
        p_glob = convert_rcif_to_glob(lcontent)
        self.data = p_glob["data"]

    def load_from_file(self, fname):
        self._p_file_dir = os.path.dirname(fname)
        self._p_file_name = os.path.basename(fname)
        
        fid = open(fname, "r")
        lcontent = fid.readlines()
        fid.close()
        self.load_from_str(lcontent)
        
    def save_to_str(self):
        lstr = rcif_to_str(self.data)
        return lstr

    def save_to_file(self, f_name):
        lstr = self.save_to_str()
        fid = open(f_name, "w")
        fid.write("\n".join(lstr)+"\n")
        fid.close()

    
    def take_from_model(self, model):
        self.data = []
        def temp_func(obj, ddata, label_rcif, label_mod, type_mod):
            if type_mod == "logic":
                sval = obj.get_val(label_mod)
                ddata[label_rcif] = sval 
            elif type_mod == "text":
                ddata[label_rcif] = obj.get_val(label_mod)
            elif type_mod == "val":
                val = obj.get_val(label_mod)
                #if val is None:
                #    print(label_rcif, label_mod, obj,"\n\n")
                if isinstance(val, Variable):
                    sval = "{:}({:})".format(val.value, val.sigma)
                else:
                    sval = "{:}".format(val)
                ddata[label_rcif] = sval 
            return
            
        drel = rcif_model_relation()
        

        
        for obj in model._list_crystal:
            d_ph = {"name": obj.get_val("name")}
            llab_rcif = drel["lab_rcif_sp"]
            llab_arg = drel["lab_arg_sp"]
            ltype_arg = drel["type_sp"]
            
            space_groupe = obj.get_val("space_groupe")
            for lab_rcif, lab_arg, type_arg in zip(llab_rcif, llab_arg, 
                                                   ltype_arg):
                temp_func(space_groupe, d_ph, lab_rcif, lab_arg, type_arg)

            llab_rcif = drel["lab_rcif_cell"]
            llab_arg = drel["lab_arg_cell"]
            ltype_arg = drel["type_cell"]

            cell = obj.get_val("cell")
            for lab_rcif, lab_arg, type_arg in zip(llab_rcif, llab_arg, 
                                                   ltype_arg):
                temp_func(cell, d_ph, lab_rcif, lab_arg, type_arg)
                
            lloop_d = []
            
            d_at = {}
            llab_rcif = drel["lab_rcif_at"]
            llab_arg = drel["lab_arg_at"]
            ltype_arg = drel["type_at"]

            atom_site = obj.get_val("atom_site")
            l_atom_type = atom_site._list_atom_type

            l_lab_rcif_used = []
            l_dhelp = []
            for obj_2 in l_atom_type:
                dhelp = {}
                for lab_rcif, lab_arg, type_arg in zip(llab_rcif, llab_arg, 
                                                   ltype_arg):
                    temp_func(obj_2, dhelp, lab_rcif, lab_arg, type_arg)
                    if (not(lab_rcif in l_lab_rcif_used)): #&(l_core in lkey)
                        l_lab_rcif_used.append(lab_rcif)
                l_dhelp.append(dhelp) 
            for l_rcif in l_lab_rcif_used:
                d_at[l_rcif] = [hh[l_rcif] for hh in l_dhelp]
            lloop_d.append(d_at)
            

            d_bscat = {}
            llab_rcif = drel["lab_rcif_bscat"]
            llab_arg = drel["lab_arg_bscat"]
            ltype_arg = drel["type_bscat"]

            l_lab_rcif_used = []
            l_dhelp = []
            for obj_2 in l_atom_type:
                dhelp = {}
                for lab_rcif, lab_arg, type_arg in zip(llab_rcif, llab_arg, 
                                                   ltype_arg):
                    temp_func(obj_2, dhelp, lab_rcif, lab_arg, type_arg)
                    if (not(lab_rcif in l_lab_rcif_used)): #&(l_core in lkey)
                        l_lab_rcif_used.append(lab_rcif)
                l_dhelp.append(dhelp) 
            lab_uniq = "_atom_site_type_symbol"
            l_val = [hh[lab_uniq] for hh in l_dhelp]
            s_uniq = set(l_val)
            l_ind_uniq = [l_val.index(hh) for hh in s_uniq]
            for l_rcif in l_lab_rcif_used:
                d_bscat[l_rcif] = [l_dhelp[ind][l_rcif] for ind in l_ind_uniq]
            lloop_d.append(d_bscat)


            d_adp = {}
            llab_rcif = drel["lab_rcif_adp"]
            llab_arg = drel["lab_arg_adp"]
            ltype_arg = drel["type_adp"]

            l_lab_rcif_used = []
            l_dhelp = []
            for obj_2 in l_atom_type:
                dhelp = {}
                flag_beta = obj_2.get_val("adp_type") == "uani"
                if flag_beta:
                    for lab_rcif, lab_arg, type_arg in zip(llab_rcif, llab_arg, 
                                                           ltype_arg):
                        temp_func(obj_2, dhelp, lab_rcif, lab_arg, type_arg)
                        if (not(lab_rcif in l_lab_rcif_used)): #&(l_core in lkey)
                            l_lab_rcif_used.append(lab_rcif)
                    l_dhelp.append(dhelp) 
            for l_rcif in l_lab_rcif_used:
                d_adp[l_rcif] = [hh[l_rcif] for hh in l_dhelp]
            lloop_d.append(d_adp)
            

            d_chi = {}
            llab_rcif = drel["lab_rcif_chi"]
            llab_arg = drel["lab_arg_chi"]
            ltype_arg = drel["type_chi"]

            l_lab_rcif_used = []
            l_dhelp = []
            for obj_2 in l_atom_type:
                dhelp = {}
                flag_m = obj_2.get_val("flag_m")
                if flag_m:
                    for lab_rcif, lab_arg, type_arg in zip(llab_rcif, llab_arg, 
                                                           ltype_arg):
                        temp_func(obj_2, dhelp, lab_rcif, lab_arg, type_arg)
                        if (not(lab_rcif in l_lab_rcif_used)): #&(l_core in lkey)
                            l_lab_rcif_used.append(lab_rcif)
                    l_dhelp.append(dhelp) 
            for l_rcif in l_lab_rcif_used:
                d_chi[l_rcif] = [hh[l_rcif] for hh in l_dhelp]
            lloop_d.append(d_chi)
            
            if lloop_d != []:
                d_ph["loops"] = lloop_d
                
            self.data.append(d_ph)

        for obj in model._list_experiment:
            d_exp = {"name": obj.get_val("name")}
            if isinstance(obj, ExperimentSingle):
                llab_rcif = drel["lab_rcif_experiment_sd"]
                llab_arg = drel["lab_arg_experiment_sd"]
                ltype_arg = drel["type_experiment_sd"]
                for lab_rcif, lab_arg, type_arg in zip(llab_rcif, llab_arg, 
                                                                    ltype_arg):
                    temp_func(obj, d_exp, lab_rcif, lab_arg, type_arg)
                d_exp["_sd_file_name_output"] = os.path.basename(
                                                d_exp["_sd_file_name_output"])
                
                observed_data = obj.get_val("observed_data")
                d_exp["_sd_file_name_input"] = observed_data.get_val("file_name")
                
                setup = obj.get_val("setup")
                beam_polarization = setup.get_val("beam_polarization")
                llab_rcif = drel["lab_rcif_beam_polarization_sd"]
                llab_arg = drel["lab_arg_beam_polarization_sd"]
                ltype_arg = drel["type_beam_polarization_sd"]
                for lab_rcif, lab_arg, type_arg in zip(llab_rcif, llab_arg, 
                                                                    ltype_arg):
                    temp_func(beam_polarization, d_exp, lab_rcif, lab_arg, 
                              type_arg)

            if isinstance(obj, ExperimentPowder1D):
                llab_rcif = drel["lab_rcif_experiment_1d"]
                llab_arg = drel["lab_arg_experiment_1d"]
                ltype_arg = drel["type_experiment_1d"]
                for lab_rcif, lab_arg, type_arg in zip(llab_rcif, llab_arg, 
                                                                    ltype_arg):
                    temp_func(obj, d_exp, lab_rcif, lab_arg, type_arg)
                d_exp["_pd_file_name_output"] = os.path.basename(
                                                d_exp["_pd_file_name_output"])
                observed_data = obj.get_val("observed_data")
                d_exp["_pd_file_name_input"] = observed_data.get_val("file_name")
                
                setup = obj.get_val("setup")
                llab_rcif = drel["lab_rcif_zero_shift_1d"]
                llab_arg = drel["lab_arg_zero_shift_1d"]
                ltype_arg = drel["type_zero_shift_1d"]
                for lab_rcif, lab_arg, type_arg in zip(llab_rcif, llab_arg, 
                                                                    ltype_arg):
                    temp_func(setup, d_exp, lab_rcif, lab_arg, type_arg)

                beam_polarization = setup.get_val("beam_polarization")
                llab_rcif = drel["lab_rcif_beam_polarization_1d"]
                llab_arg = drel["lab_arg_beam_polarization_1d"]
                ltype_arg = drel["type_beam_polarization_1d"]
                for lab_rcif, lab_arg, type_arg in zip(llab_rcif, llab_arg, 
                                                                    ltype_arg):
                    temp_func(beam_polarization, d_exp, lab_rcif, lab_arg, 
                              type_arg)
                
                resolution = setup.get_val("resolution")
                llab_rcif = drel["lab_rcif_resolution_1d"]
                llab_arg = drel["lab_arg_resolution_1d"]
                ltype_arg = drel["type_resolution_1d"]
                for lab_rcif, lab_arg, type_arg in zip(llab_rcif, llab_arg, 
                                                                    ltype_arg):
                    temp_func(resolution, d_exp, lab_rcif, lab_arg, type_arg)
                
                asymmetry = setup.get_val("asymmetry")
                llab_rcif = drel["lab_rcif_asymmetry_1d"]
                llab_arg = drel["lab_arg_asymmetry_1d"]
                ltype_arg = drel["type_asymmetry_1d"]
                for lab_rcif, lab_arg, type_arg in zip(llab_rcif, llab_arg, 
                                                                    ltype_arg):
                    temp_func(asymmetry, d_exp, lab_rcif, lab_arg, type_arg)
                
            if isinstance(obj, ExperimentPowder2D):
                llab_rcif = drel["lab_rcif_experiment_2d"]
                llab_arg = drel["lab_arg_experiment_2d"]
                ltype_arg = drel["type_experiment_2d"]
                for lab_rcif, lab_arg, type_arg in zip(llab_rcif, llab_arg, 
                                                                    ltype_arg):
                    temp_func(obj, d_exp, lab_rcif, lab_arg, type_arg)
                d_exp["_2dpd_file_name_output"] = os.path.basename(
                                                d_exp["_2dpd_file_name_output"])
                observed_data = obj.get_val("observed_data")
                d_exp["_2dpd_file_name_input"] = observed_data.get_val("file_name")
                
                setup = obj.get_val("setup")
                llab_rcif = drel["lab_rcif_zero_shift_2d"]
                llab_arg = drel["lab_arg_zero_shift_2d"]
                ltype_arg = drel["type_zero_shift_2d"]
                for lab_rcif, lab_arg, type_arg in zip(llab_rcif, llab_arg, 
                                                                    ltype_arg):
                    temp_func(setup, d_exp, lab_rcif, lab_arg, type_arg)

                beam_polarization = setup.get_val("beam_polarization")
                llab_rcif = drel["lab_rcif_beam_polarization_2d"]
                llab_arg = drel["lab_arg_beam_polarization_2d"]
                ltype_arg = drel["type_beam_polarization_2d"]
                for lab_rcif, lab_arg, type_arg in zip(llab_rcif, llab_arg, 
                                                                    ltype_arg):
                    temp_func(beam_polarization, d_exp, lab_rcif, lab_arg, 
                              type_arg)
                
                resolution = setup.get_val("resolution")
                llab_rcif = drel["lab_rcif_resolution_2d"]
                llab_arg = drel["lab_arg_resolution_2d"]
                ltype_arg = drel["type_resolution_2d"]
                for lab_rcif, lab_arg, type_arg in zip(llab_rcif, llab_arg, 
                                                                    ltype_arg):
                    temp_func(resolution, d_exp, lab_rcif, lab_arg, type_arg)
                
                asymmetry = setup.get_val("asymmetry")
                llab_rcif = drel["lab_rcif_asymmetry_2d"]
                llab_arg = drel["lab_arg_asymmetry_2d"]
                ltype_arg = drel["type_asymmetry_2d"]
                for lab_rcif, lab_arg, type_arg in zip(llab_rcif, llab_arg, 
                                                                    ltype_arg):
                    temp_func(asymmetry, d_exp, lab_rcif, lab_arg, type_arg)
                
                
            d_eph = {}
            lloop_d = []
            l_calculated_data = obj._list_calculated_data
            if isinstance(obj, ExperimentSingle):
                crystal = model._list_crystal[0]
                extinction = crystal.get_val("extinction")
                llab_rcif = drel["lab_rcif_extinction"]
                llab_arg = drel["lab_arg_extinction"]
                ltype_arg = drel["type_extinction"]
                
                l_dhelp = []
                for calculated_data in l_calculated_data:
                    name = calculated_data.get_val("name")
                    dhelp = {"_sd_phase_name": name}
                    for lab_rcif, lab_arg, type_arg in zip(llab_rcif, llab_arg, 
                                                                    ltype_arg):
                        temp_func(extinction, dhelp, lab_rcif, lab_arg, type_arg)
                    l_dhelp.append(dhelp) 
                for lab_rcif in ["_sd_phase_name"]+llab_rcif:
                    d_eph[lab_rcif] = [hh[lab_rcif] for hh in l_dhelp]
                lloop_d.append(d_eph)

            if isinstance(obj, ExperimentPowder1D):
                llab_rcif = ["_pd_phase_name", "_pd_phase_scale"]
                llab_arg = ["name", "scale"]
                ltype_arg = ["text", "val"]

                l_dhelp = []
                for calculated_data in l_calculated_data:
                    dhelp = {}
                    for lab_rcif, lab_arg, type_arg in zip(llab_rcif, llab_arg, 
                                                                    ltype_arg):
                        temp_func(calculated_data, dhelp, lab_rcif, lab_arg, type_arg)


                    name = calculated_data.get_val("name")
                    ind = 0
                    for i_crystal, crystal in enumerate(model._list_crystal):
                        if crystal.get_val("name") == name:
                            ind = i_crystal
                            break

                    crystal = model._list_crystal[ind]
                    

                    temp_func(crystal, dhelp, ["_pd_phase_igsize"], ["i_g"], ["val"])
                    l_dhelp.append(dhelp) 

                for lab_rcif in llab_rcif:
                    d_eph[lab_rcif] = [hh[lab_rcif] for hh in l_dhelp]
                lloop_d.append(d_eph)

            if isinstance(obj, ExperimentPowder2D):
                llab_rcif = ["_2dpd_phase_name", "_2dpd_phase_scale"]
                llab_arg = ["name", "scale"]
                ltype_arg = ["text", "val"]

                l_dhelp = []
                for calculated_data in l_calculated_data:
                    dhelp = {}
                    for lab_rcif, lab_arg, type_arg in zip(llab_rcif, llab_arg, 
                                                                    ltype_arg):
                        temp_func(calculated_data, dhelp, lab_rcif, lab_arg, type_arg)

                    name = calculated_data.get_val("name")
                    ind = 0
                    for i_crystal, crystal in enumerate(model._list_crystal):
                        if crystal.get_val("name") == name:
                            ind = i_crystal
                            break
                    crystal = model._list_crystal[ind]


                    temp_func(crystal, dhelp, ["_2dpd_phase_igsize"], ["i_g"], ["val"])
                    l_dhelp.append(dhelp) 

                for lab_rcif in llab_rcif:
                    d_eph[lab_rcif] = [hh[lab_rcif] for hh in l_dhelp]
                lloop_d.append(d_eph)
                
            if lloop_d != []:
                d_exp["loops"] = lloop_d
                
            self.data.append(d_exp)

    def trans_to_model(self):
        l_crystal = []
        l_experiment = []
        for data in self.data:
            lab_crystal = "_cell"
            lab_exp_pd_1d = "_pd"
            lab_exp_pd_2d = "_2dpd"
            lab_exp_sd = "_sd"
            lab_ref = "_refinement"
            lkey = data.keys()
            flag_ph = any([hh.startswith(lab_crystal) for hh in lkey])
            flag_exp_pd = any([hh.startswith(lab_exp_pd_1d) for hh in lkey])
            flag_exp_2dpd = any([hh.startswith(lab_exp_pd_2d) for hh in lkey])
            flag_exp_sd = any([hh.startswith(lab_exp_sd) for hh in lkey])
            flag_ref = any([hh.startswith(lab_ref) for hh in lkey])
            if flag_ph:
                crystal, l_variable_crystal = trans_to_crystal(data)
                l_crystal.append(crystal)
            elif flag_exp_pd:
                f_dir = self._p_file_dir
                experiment, l_variable_experiment = trans_pd_to_experiment(f_dir, data, l_crystal)
                l_experiment.append(experiment)
            elif flag_exp_2dpd:
                f_dir = self._p_file_dir
                experiment, l_variable_experiment = trans_2dpd_to_experiment(f_dir, data, l_crystal)
                l_experiment.append(experiment)
            elif flag_exp_sd:
                f_dir = self._p_file_dir
                experiment, l_variable_experiment = trans_sd_to_experiment(f_dir, data, l_crystal)
                l_experiment.append(experiment)
            elif flag_ref:
                pass
                #refinement = trans_to_refinement(data)
                #l_refinement.append(refinement)
                #l_variable.extend(l_variable_experiment)


        model = Model()
        for experiment in l_experiment:
            model.add_experiment(experiment)
        for crystal in l_crystal:
            model.add_crystal(crystal) 
        
        return model

def get_link_rcif(link_core, ccore):
    """
    """
    link_core_2 = "".join([hh for hh in link_core.split()])[1:-1]
    
    lhelp = link_core_2.split(')(')
    type_1 = lhelp[0].split(",")[0]
    name_1 = lhelp[0].split(",")[1]
    link_rcif = "${:}".format(name_1)
    if len(lhelp) == 3:
        name_2 = lhelp[1].split(",")[1]
        link_rcif += "_${:}".format(name_2)
    name_core = lhelp[-1]
    
    drel = rcif_model_relation()
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
        print("Mistake in link '{:}'".format(param))
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
        print("Mistake with naming of phase or experiment in link '{:}'".format(
                param))
        return ""
    if ind == -1:
        print("Code word in '{:}' is not defined".format(param))
        return ""
    return slink 


def rcif_to_str(l_ddata):
    #structure of cif is data blokcs and loop blocks, 
    #nested loops are not supported
    ls_out = []
    for ddata in l_ddata:    
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
                        ls_out.append("\nloop_{:}".format(lab, d_loop[lab].strip()))
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

def put_ref(obj, lparam):
    drel = rcif_model_relation()
    for param in lparam:
        slink = get_link(obj, drel, param)
        val_1, message = obj.set_val_by_link(slink, None)
        val_1[1] = True
        obj.set_val_by_link(slink, val_1)


def put_con(obj, lparam1, lparam2, lcoeff):
    drel = rcif_model_relation()
    for param1, param2, coeff in zip(lparam1, lparam2, lcoeff):
        slink_obj = get_link(obj, drel, param1)
        slink_sub = get_link(obj, drel, param2)
        s_constr = " {:}*x1 [{:}]".format(coeff, slink_sub)
        val_1, message = obj.set_val_by_link(slink_obj, None)
        val_1[1] = False
        val_1[2] = s_constr
        obj.set_val_by_link(slink_obj, val_1)    


def rcif_model_relation():
    #phase    
    llab_rcif_cell = ["_cell_length_a", "_cell_length_b", "_cell_length_c",
                 "_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma"]
    llab_arg_cell = ["a", "b", "c",
                  "alpha", "beta", "gamma"]
    ltype_cell = ["val", "val", "val", "val", "val", "val"]

    llab_rcif_sp = ["_space_group_name_H-M_alt", '_space_group_it_coordinate_system_code']
    llab_arg_sp = ['spgr_given_name', 'spgr_choice']
    ltype_sp = ["text", "text"]


    llab_rcif_at = ["_atom_site_label", '_atom_site_type_symbol',
                "_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z",
                "_atom_site_occupancy", "_atom_site_b_iso_or_equiv", "_atom_site_adp_type"]
    llab_arg_at = ["name", "type_n", "x", "y", "z", "occupation", "b_iso", "adp_type"]
    ltype_at = ["text", "text", "val", "val", "val", "val", "val", "text"]

    llab_rcif_bscat = ["_atom_site_type_symbol", "_atom_site_bscat"]
    llab_arg_bscat = ["type_n", "b_scat"]
    ltype_bscat = ["text", "val"]

    llab_rcif_adp = ["_atom_site_aniso_label", '_atom_site_aniso_U_11',
             "_atom_site_aniso_U_22", "_atom_site_aniso_U_33", "_atom_site_aniso_U_12",
             "_atom_site_aniso_U_13", "_atom_site_aniso_U_23"]
    llab_arg_adp = ["name", "u_11", "u_22", "u_33", "u_12", 
                    "u_13", "u_23"]
    ltype_adp = ["text", "val", "val", "val", "val", "val", "val"]


    llab_rcif_chi = ["_atom_site_magnetism_aniso_label", '_atom_site_magnetism_type_symbol',
             "_atom_site_magnetism_kappa",
             "_atom_site_magnetism_factor_lande", "_atom_site_susceptibility_aniso_chi_11",
             "_atom_site_susceptibility_aniso_chi_22", "_atom_site_susceptibility_aniso_chi_33",
             "_atom_site_susceptibility_aniso_chi_23", "_atom_site_susceptibility_aniso_chi_13",
             "_atom_site_susceptibility_aniso_chi_12", "_atom_site_magnetism_type"]
    llab_arg_chi = ["name", "type_m", "kappa", "factor_lande", "chi_11", 
                    "chi_22", "chi_33", "chi_12", "chi_13", "chi_23", "chi_type"]
    ltype_chi = ["text", "text", "val", "val", "val", "val", "val", "val", 
                 "val", "val", "text"]




    #ExperimentPowder1D
    llab_rcif_experiment_1d = ["name", "_pd_file_name_output"]
    llab_arg_experiment_1d = ["name", "f_out"]
    ltype_experiment_1d = ["text", "text"]
    
    llab_rcif_beam_polarization_1d = ["_pd_beam_polarization_up", "_pd_beam_polarization_down"]
    llab_arg_beam_polarization_1d = ["p_u", "p_d"]
    ltype_beam_polarization_1d = ["val", "val"]

    llab_rcif_resolution_1d = ["_pd_resolution_u", "_pd_resolution_v", 
                               "_pd_resolution_w", "_pd_resolution_x", 
                               "_pd_resolution_y"]
    llab_arg_resolution_1d = ["u", "v", "w", "x", "y"]
    ltype_resolution_1d = ["val", "val", "val", "val", "val"]


    llab_rcif_asymmetry_1d = ["_pd_reflex_asymmetry_p1", 
                              "_pd_reflex_asymmetry_p2", 
                              "_pd_reflex_asymmetry_p3",
                              "_pd_reflex_asymmetry_p4"]
    llab_arg_asymmetry_1d = ["p1", "p2", "p3", "p4"]
    ltype_asymmetry_1d = ["val", "val", "val", "val"]


    llab_rcif_zero_shift_1d = ["_pd_shift_const"]
    llab_arg_zero_shift_1d = ["zero_shift"]
    ltype_zero_shift_1d = ["val"]



    #ExperimentPowder2D
    llab_rcif_experiment_2d = ["name", "_2dpd_file_name_output"]
    llab_arg_experiment_2d = ["name", "f_out"]
    ltype_experiment_2d = ["text", "text"]
    
    llab_rcif_beam_polarization_2d = ["_2dpd_beam_polarization_up", "_2dpd_beam_polarization_down"]
    llab_arg_beam_polarization_2d = ["p_u", "p_d"]
    ltype_beam_polarization_2d = ["val", "val"]


    llab_rcif_resolution_2d = ["_2dpd_resolution_u", "_2dpd_resolution_v", 
                               "_2dpd_resolution_w", "_2dpd_resolution_x", 
                               "_2dpd_resolution_y"]
    llab_arg_resolution_2d = ["u", "v", "w", "x", "y"]
    ltype_resolution_2d = ["val", "val", "val", "val", "val"]


    llab_rcif_asymmetry_2d = ["_2dpd_reflex_asymmetry_p1", 
                              "_2dpd_reflex_asymmetry_p2", 
                              "_2dpd_reflex_asymmetry_p3",
                              "_2dpd_reflex_asymmetry_p4"]
    llab_arg_asymmetry_2d = ["p1", "p2", "p3", "p4"]
    ltype_asymmetry_2d = ["val", "val", "val", "val"]


    llab_rcif_zero_shift_2d = ["_2dpd_shift_const"]
    llab_arg_zero_shift_2d = ["zero_shift"]
    ltype_zero_shift_2d = ["val"]



    #ExperimentSingle
    llab_rcif_beam_polarization_sd = ["_sd_beam_polarization_up", "_sd_beam_polarization_down"]
    llab_arg_beam_polarization_sd = ["p_u", "p_d"]
    ltype_beam_polarization_sd = ["val", "val"]

    llab_rcif_extinction = ["_sd_phase_extinction_radius", "_sd_phase_extinction_mosaicity"]
    llab_arg_extinction = ["domain_radius", "mosaicity"]
    ltype_extinction = ["val", "val"]


    llab_rcif_experiment_sd = ["name", "_sd_file_name_output"]
    llab_arg_experiment_sd = ["name", "f_out"]
    ltype_experiment_sd = ["text", "text"]


    #refinement
    llab_rcif_ref = ["_refinement_file_name_output", "_refinement_flag_errorbars",
                 "_refinement_flag_refinement"]
    llab_arg_ref = ["output", "refin", "sigmas"]
    ltype_ref = ["text", "logic", "logic"]

    drel = {}

    drel["lab_rcif_cell"] = llab_rcif_cell
    drel["lab_arg_cell"] = llab_arg_cell
    drel["type_cell"] = ltype_cell
    
    drel["lab_rcif_sp"] = llab_rcif_sp
    drel["lab_arg_sp"] = llab_arg_sp
    drel["type_sp"] = ltype_sp

    drel["lab_rcif_atom_type"] = (llab_rcif_at+llab_rcif_bscat+llab_rcif_adp+
                                  llab_rcif_chi)
    drel["lab_arg_atom_type"] = (llab_arg_at+llab_arg_bscat+llab_arg_adp+
                                 llab_arg_chi)
    drel["type_atom_type"] = (ltype_at+ltype_bscat+ltype_adp+ltype_chi)

    
    drel["lab_rcif_at"] = llab_rcif_at
    drel["lab_arg_at"] = llab_arg_at
    drel["type_at"] = ltype_at
    
    drel["lab_rcif_bscat"] = llab_rcif_bscat
    drel["lab_arg_bscat"] = llab_arg_bscat
    drel["type_bscat"] = ltype_bscat
    
    drel["lab_rcif_adp"] = llab_rcif_adp
    drel["lab_arg_adp"] = llab_arg_adp
    drel["type_adp"] = ltype_adp
    
    drel["lab_rcif_chi"] = llab_rcif_chi
    drel["lab_arg_chi"] = llab_arg_chi
    drel["type_chi"] = ltype_chi


    drel["lab_rcif_experiment_1d"] = llab_rcif_experiment_1d
    drel["lab_arg_experiment_1d"] = llab_arg_experiment_1d 
    drel["type_experiment_1d"] = ltype_experiment_1d

    drel["lab_rcif_resolution_1d"] = llab_rcif_resolution_1d
    drel["lab_arg_resolution_1d"] = llab_arg_resolution_1d 
    drel["type_resolution_1d"] = ltype_resolution_1d 

    drel["lab_rcif_asymmetry_1d"] = llab_rcif_asymmetry_1d 
    drel["lab_arg_asymmetry_1d"] = llab_arg_asymmetry_1d
    drel["type_asymmetry_1d"] = ltype_asymmetry_1d 

    drel["lab_rcif_zero_shift_1d"] = llab_rcif_zero_shift_1d
    drel["lab_arg_zero_shift_1d"] = llab_arg_zero_shift_1d 
    drel["type_zero_shift_1d"] = ltype_zero_shift_1d 



    drel["lab_rcif_experiment_2d"] = llab_rcif_experiment_2d
    drel["lab_arg_experiment_2d"] = llab_arg_experiment_2d 
    drel["type_experiment_2d"] = ltype_experiment_2d

    drel["lab_rcif_resolution_2d"] = llab_rcif_resolution_2d
    drel["lab_arg_resolution_2d"] = llab_arg_resolution_2d 
    drel["type_resolution_2d"] = ltype_resolution_2d 

    drel["lab_rcif_asymmetry_2d"] = llab_rcif_asymmetry_2d 
    drel["lab_arg_asymmetry_2d"] = llab_arg_asymmetry_2d
    drel["type_asymmetry_2d"] = ltype_asymmetry_2d 

    drel["lab_rcif_zero_shift_2d"] = llab_rcif_zero_shift_2d
    drel["lab_arg_zero_shift_2d"] = llab_arg_zero_shift_2d 
    drel["type_zero_shift_2d"] = ltype_zero_shift_2d 



    drel["lab_rcif_beam_polarization_sd"] = llab_rcif_beam_polarization_sd
    drel["lab_arg_beam_polarization_sd"] = llab_arg_beam_polarization_sd
    drel["type_beam_polarization_sd"] = ltype_beam_polarization_sd

    drel["lab_rcif_beam_polarization_1d"] = llab_rcif_beam_polarization_1d
    drel["lab_arg_beam_polarization_1d"] = llab_arg_beam_polarization_1d
    drel["type_beam_polarization_1d"] = ltype_beam_polarization_1d

    drel["lab_rcif_beam_polarization_2d"] = llab_rcif_beam_polarization_2d
    drel["lab_arg_beam_polarization_2d"] = llab_arg_beam_polarization_2d
    drel["type_beam_polarization_2d"] = ltype_beam_polarization_2d


    drel["lab_rcif_extinction"] = llab_rcif_extinction
    drel["lab_arg_extinction"] = llab_arg_extinction
    drel["type_extinction"] = ltype_extinction


    drel["lab_rcif_experiment_sd"] = llab_rcif_experiment_sd
    drel["lab_arg_experiment_sd"] = llab_arg_experiment_sd
    drel["type_experiment_sd"] = ltype_experiment_sd


    
    drel["lab_rcif_ref"] = llab_rcif_ref
    drel["lab_arg_ref"] = llab_arg_ref
    drel["type_ref"] = ltype_ref
    
    return drel




def from_dict_to_obj(dict_i, llab_d, obj, llab_o, ltype):
    lkey = dict_i.keys()
    lnumb = [ihh for ihh, hh in enumerate(llab_d) if (hh in lkey)]
    d_args = {}
    for numb in lnumb:
        if ltype[numb] == "val":
            #val = [dict_i[llab_d[numb]], False, ""]
            val = dict_i[llab_d[numb]]
            val = conv_str_to_text_float_logic(val, llab_o[numb])
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


def trans_to_crystal(data):
    """
    transform info in dictionary to class Crystal and give the list of refined parameters
    """
    l_variable = []
    
    space_groupe = SpaceGroupe()
    cell = Cell()
    atom_site = AtomSite()
    crystal = Crystal(cell=cell, atom_site=atom_site, 
                      space_groupe=space_groupe)
    
    
    crystal.set_val(name=data["name"])
    
    drel = rcif_model_relation()
    
    llab_rcif = drel["lab_rcif_sp"]
    llab_crystal = drel["lab_arg_sp"]
    ltype = drel["type_sp"]
    
    from_dict_to_obj(data, llab_rcif, space_groupe, llab_crystal, ltype)

    llab_rcif = drel["lab_rcif_cell"]
    llab_crystal = drel["lab_arg_cell"]
    ltype = drel["type_cell"]
    
    from_dict_to_obj(data, llab_rcif, cell, llab_crystal, ltype)
    
    cell.set_val(singony=space_groupe.get_val("singony"))
    

    lkey = data.keys()
    if (not ("loops" in lkey)):
        return crystal, l_variable

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
    
    l_atom = []

    if lnumb[0] != []:
        llab_rcif = drel["lab_rcif_at"]
        llab_arg = drel["lab_arg_at"]
        ltype = drel["type_at"]
        for numb in lnumb[0]:
            l_key = list(data["loops"][numb].keys())
            n_atom = len(data["loops"][numb][l_key[0]])
            for i_atom in range(n_atom):
                dd = {}
                for key in l_key:
                    dd.update({key:data["loops"][numb][key][i_atom]})
                atom_1 = AtomType()
                from_dict_to_obj(dd, llab_rcif, atom_1, llab_arg, ltype)
                l_atom.append(atom_1)
    else:
        print("Fractional coordinates of atoms are necessaire")
        return crystal, l_variable
    """    
    if lnumb[3] != []:
        llab_rcif = drel["lab_rcif_bscat"]
        llab_arg = drel["lab_arg_bscat"]
        ltype = drel["type_bscat"]
        for numb in lnumb[3]:
            l_key = list(data["loops"][numb].keys())
            n_atom = len(data["loops"][numb][l_key[0]])
            for i_atom in range(n_atom):
                dd = {}
                for key in l_key:
                    dd.update({key:data["loops"][numb][key][i_atom]})
                for atom_1 in l_atom:
                    if dd["_atom_site_type_symbol"] == atom_1.get_val("type_n"):
                        from_dict_to_obj(dd, llab_rcif, atom_1, llab_arg, ltype)
    """
    if lnumb[1] != []:
        llab_rcif = drel["lab_rcif_adp"]
        llab_arg = drel["lab_arg_adp"]
        ltype = drel["type_adp"]
        for numb in lnumb[1]:
            l_key = list(data["loops"][numb].keys())
            n_atom = len(data["loops"][numb][l_key[0]])
            for i_atom in range(n_atom):
                dd = {}
                for key in l_key:
                    dd.update({key:data["loops"][numb][key][i_atom]})
                for atom_1 in l_atom:
                    if dd["_atom_site_aniso_label"] == atom_1.get_val("name"):
                        from_dict_to_obj(dd, llab_rcif, atom_1, llab_arg, ltype)
    if lnumb[2] != []:
        llab_rcif = drel["lab_rcif_chi"]
        llab_arg = drel["lab_arg_chi"]
        ltype = drel["type_chi"]
        for numb in lnumb[2]:
            l_key = list(data["loops"][numb].keys())
            n_atom = len(data["loops"][numb][l_key[0]])
            for i_atom in range(n_atom):
                dd = {}
                for key in l_key:
                    dd.update({key:data["loops"][numb][key][i_atom]})
                for atom_1 in l_atom:
                    if dd["_atom_site_magnetism_aniso_label"] == atom_1.get_val("name"):
                        from_dict_to_obj(dd, llab_rcif, atom_1, llab_arg, ltype)
                        atom_1.set_val(flag_m=True)
            
    if lnumb[4] != []:
        pass
    if lnumb[5] != []:
        pass
    
    for atom_1 in l_atom:
        atom_site.add_atom(atom_1)
    return crystal, l_variable


def trans_pd_to_experiment(f_dir, data, l_crystal):
    """
    transform info in dictionary to class ExperimentPowder1D and give the list of refined parameters
    """
    l_variable = []
    beam_polarization = BeamPolarization()
    resolution_powder_1d = ResolutionPowder1D()
    factor_lorentz_powder_1d = FactorLorentzPowder1D()
    asymmetry_powder_1d = AsymmetryPowder1D()
    background_powder_1d = BackgroundPowder1D()
    setup_powder_1d = SetupPowder1D(resolution=resolution_powder_1d, 
            factor_lorentz=factor_lorentz_powder_1d, 
            asymmetry=asymmetry_powder_1d, beam_polarization=beam_polarization, 
            background=background_powder_1d)
        
    observed_data_powder_1d = ObservedDataPowder1D()
    
    experiment = ExperimentPowder1D(setup=setup_powder_1d, 
                                    observed_data=observed_data_powder_1d)

    
    drel = rcif_model_relation()
    
    mode_chi_sq = ""
    l_key = data.keys()
    if "_pd_chi2_up" in l_key:
        if data["_pd_chi2_up"]:
            mode_chi_sq+="up "
    if "_pd_chi2_down" in l_key:
        if data["_pd_chi2_down"]:
            mode_chi_sq+="down "
    if "_pd_chi2_diff" in l_key:
        if data["_pd_chi2_diff"]:
            mode_chi_sq+="diff "
    if "_pd_chi2_sum" in l_key:
        if data["_pd_chi2_sum"]:
            mode_chi_sq+="sum "
    if mode_chi_sq != "":
        experiment.set_val(mode_chi_sq=mode_chi_sq)


    llab_rcif = drel["lab_rcif_experiment_1d"]
    llab_arg = drel["lab_arg_experiment_1d"]
    ltype = drel["type_experiment_1d"]

    from_dict_to_obj(data, llab_rcif, experiment, llab_arg, ltype)    
    if experiment.get_val("f_out") is not None:
        f_out = experiment.get_val("f_out")
        experiment.set_val(f_out=os.path.join(f_dir, f_out))
    

    llab_rcif = drel["lab_rcif_resolution_1d"]
    llab_arg = drel["lab_arg_resolution_1d"]
    ltype = drel["type_resolution_1d"]

    from_dict_to_obj(data, llab_rcif, resolution_powder_1d, llab_arg, ltype)

    llab_rcif = drel["lab_rcif_asymmetry_1d"]
    llab_arg = drel["lab_arg_asymmetry_1d"]
    ltype = drel["type_asymmetry_1d"]

    from_dict_to_obj(data, llab_rcif, asymmetry_powder_1d, llab_arg, ltype)

    llab_rcif = drel["lab_rcif_zero_shift_1d"]
    llab_arg = drel["lab_arg_zero_shift_1d"]
    ltype = drel["type_zero_shift_1d"]

    from_dict_to_obj(data, llab_rcif, setup_powder_1d, llab_arg, ltype)


    llab_rcif = drel["lab_rcif_beam_polarization_1d"]
    llab_arg = drel["lab_arg_beam_polarization_1d"]
    ltype = drel["type_beam_polarization_1d"]

    from_dict_to_obj(data, llab_rcif, beam_polarization, llab_arg, ltype)


    f_inp = data["_pd_file_name_input"]
    observed_data_powder_1d.read_data(os.path.join(f_dir, f_inp))
    
    
    f_inp = data["_pd_file_name_bkgr"]
    background_powder_1d.read_data(os.path.join(f_dir, f_inp))
    
    wave_length = observed_data_powder_1d.get_val("wave_length")
    setup_powder_1d.set_val(wave_length=wave_length)

    field = observed_data_powder_1d.get_val("field")

    l_key = data.keys()
    if (not ("loops" in l_key)):
        return experiment, l_variable 
    
    
    
    for data_ph in data["loops"]:
        l_key = list(data_ph.keys())
        flag_scale = "_pd_phase_scale" in l_key
        if flag_scale:
            n_crystal = len(data_ph["_pd_phase_name"])
            for i_crystal in range(n_crystal):
                dd = {}
                for key in l_key:
                    dd.update({key: data_ph[key][i_crystal]})
                for crystal in l_crystal:
                    if  dd["_pd_phase_name"] == crystal.get_val("name"):
                        scale = dd["_pd_phase_scale"]
                        i_g = dd["_pd_phase_igsize"]
                        crystal.set_val(i_g=i_g)
                        name = crystal.get_val("name")
                        calculated_data = CalculatedDataPowder1D(field=field,
                                scale=scale, name=name)
                        experiment.add_calculated_data(calculated_data)
    return experiment, l_variable

def trans_2dpd_to_experiment(f_dir, data, l_crystal):
    """
    transform info in dictionary to class ExperimentPowder2D and give the list of refined parameters
    """
    l_variable = []
    beam_polarization = BeamPolarization()
    resolution_powder_2d = ResolutionPowder2D()
    factor_lorentz_powder_2d = FactorLorentzPowder2D()
    asymmetry_powder_2d = AsymmetryPowder2D()
    background_powder_2d = BackgroundPowder2D()
    setup_powder_2d = SetupPowder2D(resolution=resolution_powder_2d, 
            factor_lorentz=factor_lorentz_powder_2d, 
            asymmetry=asymmetry_powder_2d, beam_polarization=beam_polarization, 
            background=background_powder_2d)
        
    observed_data_powder_2d = ObservedDataPowder2D()
    
    experiment = ExperimentPowder2D(setup=setup_powder_2d, 
                                    observed_data=observed_data_powder_2d)

    
    drel = rcif_model_relation()
    
    mode_chi_sq = ""
    l_key = data.keys()
    if "_2dpd_chi2_up" in l_key:
        if data["_2dpd_chi2_up"]:
            mode_chi_sq+="up "
    if "_2dpd_chi2_down" in l_key:
        if data["_2dpd_chi2_down"]:
            mode_chi_sq+="down "
    if "_2dpd_chi2_diff" in l_key:
        if data["_2dpd_chi2_diff"]:
            mode_chi_sq+="diff "
    if "_2dpd_chi2_sum" in l_key:
        if data["_2dpd_chi2_sum"]:
            mode_chi_sq+="sum "
    if mode_chi_sq != "":
        experiment.set_val(mode_chi_sq=mode_chi_sq)


    llab_rcif = drel["lab_rcif_experiment_2d"]
    llab_arg = drel["lab_arg_experiment_2d"]
    ltype = drel["type_experiment_2d"]

    from_dict_to_obj(data, llab_rcif, experiment, llab_arg, ltype)    
    if experiment.get_val("f_out") is not None:
        f_out = experiment.get_val("f_out")
        experiment.set_val(f_out=os.path.join(f_dir, f_out))


    llab_rcif = drel["lab_rcif_resolution_2d"]
    llab_arg = drel["lab_arg_resolution_2d"]
    ltype = drel["type_resolution_2d"]

    from_dict_to_obj(data, llab_rcif, resolution_powder_2d, llab_arg, ltype)

    llab_rcif = drel["lab_rcif_asymmetry_2d"]
    llab_arg = drel["lab_arg_asymmetry_2d"]
    ltype = drel["type_asymmetry_2d"]

    from_dict_to_obj(data, llab_rcif, asymmetry_powder_2d, llab_arg, ltype)

    llab_rcif = drel["lab_rcif_zero_shift_2d"]
    llab_arg = drel["lab_arg_zero_shift_2d"]
    ltype = drel["type_zero_shift_2d"]

    from_dict_to_obj(data, llab_rcif, setup_powder_2d, llab_arg, ltype)


    llab_rcif = drel["lab_rcif_beam_polarization_2d"]
    llab_arg = drel["lab_arg_beam_polarization_2d"]
    ltype = drel["type_beam_polarization_2d"]

    from_dict_to_obj(data, llab_rcif, beam_polarization, llab_arg, ltype)


    f_inp=data["_2dpd_file_name_input"]
    observed_data_powder_2d.read_data(os.path.join(f_dir, f_inp))

    tth_min = data["_2dpd_tth_min"]
    tth_max = data["_2dpd_tth_max"]
    phi_min = data["_2dpd_phi_min"]
    phi_max = data["_2dpd_phi_max"]
    observed_data_powder_2d.set_val(tth_min=tth_min, tth_max=tth_max, 
                                    phi_min=phi_min, phi_max=phi_max,)
    
    f_inp=data["_2dpd_file_name_bkgr"]
    background_powder_2d.read_data(os.path.join(f_dir, f_inp))
    
    wave_length = observed_data_powder_2d.get_val("wave_length")
    setup_powder_2d.set_val(wave_length=wave_length)

    field = observed_data_powder_2d.get_val("field")

    l_key = data.keys()
    if (not ("loops" in l_key)):
        return experiment, l_variable 
    
    
    
    for data_ph in data["loops"]:
        l_key = list(data_ph.keys())
        flag_scale = "_2dpd_phase_scale" in l_key
        if flag_scale:
            n_crystal = len(data_ph["_2dpd_phase_name"])
            for i_crystal in range(n_crystal):
                dd = {}
                for key in l_key:
                    dd.update({key: data_ph[key][i_crystal]})
                for crystal in l_crystal:
                    if  dd["_2dpd_phase_name"] == crystal.get_val("name"):
                        scale = dd["_2dpd_phase_scale"]
                        i_g = dd["_2dpd_phase_igsize"]
                        crystal.set_val(i_g=i_g)
                        name = crystal.get_val("name")
                        calculated_data = CalculatedDataPowder2D(field=field,
                                scale=scale, name=name)
                        experiment.add_calculated_data(calculated_data)
    return experiment, l_variable


def trans_sd_to_experiment(f_dir, data, l_crystal):
    """
    transform info in dictionary to class ExperimentSingle and give the list of refined parameters
    """
    l_variable = []
    experiment = ExperimentSingle()
    setup = experiment.get_val("setup")
    beam_polarization = setup.get_val("beam_polarization")
    observed_data = experiment.get_val("observed_data")


    experiment.set_val(name=data["name"])

    drel = rcif_model_relation()


    llab_rcif = drel["lab_rcif_experiment_sd"] 
    llab_arg = drel["lab_arg_experiment_sd"] 
    ltype = drel["type_experiment_sd"] 
    
    from_dict_to_obj(data, llab_rcif, experiment, llab_arg, ltype)
    if experiment.get_val("f_out") is not None:
        f_out = experiment.get_val("f_out")
        experiment.set_val(f_out=os.path.join(f_dir, f_out))
    
    
    llab_rcif = drel["lab_rcif_beam_polarization_sd"]
    llab_arg = drel["lab_arg_beam_polarization_sd"]
    ltype = drel["type_beam_polarization_sd"]

    from_dict_to_obj(data, llab_rcif, beam_polarization, llab_arg, ltype)


    f_inp=data["_sd_file_name_input"]
    observed_data.read_data(os.path.join(f_dir, f_inp))
    
    
    wave_length = observed_data.get_val("wave_length")
    setup.set_val(wave_length=wave_length)

    field = observed_data.get_val("field")
    orientation = observed_data.get_val("orientation")

    lkey = data.keys()
    if (not ("loops" in lkey)):
        return experiment, l_variable 
    
    
    lab_exp = "_sd_phase_name"
    
    #llab_rcif = drel["lab_rcif_observed_data_sd"]
    #llab_arg = drel["lab_observed_data_sd"]
    #ltype = drel["type_observed_data_sd"]
    
    llab_rcif = drel["lab_rcif_extinction"]
    llab_arg = drel["lab_arg_extinction"]
    ltype = drel["type_extinction"]
    
    
    for data_ph in data["loops"]:
        l_key = list(data_ph.keys())
        flag_extinction = "_sd_phase_extinction_radius" in l_key
        #flag_scale = "_sd_phase_scale" in lkey
        if flag_extinction:
            n_crystal = len(data_ph["_sd_phase_name"])
            for i_crystal in range(n_crystal):
                dd = {}
                for key in l_key:
                    dd.update({key: data_ph[key][i_crystal]})
                for crystal in l_crystal:
                    if  dd["_sd_phase_name"] == crystal.get_val("name"):
                        extinction = crystal.get_val("extinction")
                        from_dict_to_obj(dd, llab_rcif, extinction, llab_arg, ltype)
                        name = crystal.get_val("name")
                        #should crystal be a deepcopy or not???
                        calculated_data = CalculatedDataSingle(field=field, 
                                orientation=orientation, name=name)
                        experiment.add_calculated_data(calculated_data)
                        
    return experiment, l_variable


def trans_to_refinement(data):
    """
    transform info in dictionary to intermediate class 
    """
    ref = ccore.ccore_ref()
    drel = rcif_model_relation()
    
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
    return ref


def load_crystal_to_experiment(l_experiment, l_crystal):
    pass

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
    lval_2 = []
    for val, name in zip(lval, l_name):
        val_2 = conv_str_to_text_float_logic(val, name)
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
        res[name] = conv_str_to_text_float_logic(value, name)
    return res

def conv_str_to_text_float_logic(sval, name=""):
    if (not (isinstance(sval,str))):
        return sval
    try:
        if len(sval.strip().split("("))>1:
            l_help = sval.split("(")
            val_1 = float(l_help[0])
            val = Variable(val_1, True, name)
            pass
        else:
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
