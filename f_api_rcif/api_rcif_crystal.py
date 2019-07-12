"""
transform inforamtion from cif file to class phase
"""
import sys
import os

import f_rcif.cl_rcif 

import f_crystal.cl_crystal


import f_api_rcif.api_rcif_common

import f_rhochi_model.cl_model 
import f_common.cl_variable 


def conv_data_to_crystal(l_data):
    """
    transform info from list of dictionary intolist of Crystal  class
    """
    l_flag_crystal = []
    for data in l_data:
        l_relation = data_space_groupe_relation()
        l_keys = data.keys()
        flag_crystal = any(relation[0] in l_keys for relation in l_relation)
        l_flag_crystal.append(flag_crystal)

    l_data_crystal = [data for data, flag in zip(l_data, l_flag_crystal) if flag]
    l_crystal = []
    for data in l_data_crystal:

        space_groupe = f_crystal.cl_crystal.SpaceGroupe()
        cell = f_crystal.cl_crystal.Cell()
        atom_site = f_crystal.cl_crystal.AtomSite()

        crystal = f_crystal.cl_crystal.Crystal(cell=cell, atom_site=atom_site, 
                      space_groupe=space_groupe)

        crystal.set_val(name=data["name"])
    
        l_relation = data_space_groupe_relation()
        f_api_rcif.api_rcif_common.from_dict_to_obj(data, l_relation, space_groupe)

        l_relation = data_cell_relation()
        f_api_rcif.api_rcif_common.from_dict_to_obj(data, l_relation, cell)
    
        cell.set_val(singony=space_groupe.get_val("singony"))

        l_key = data.keys()
        if (not ("loops" in l_key)):
            return crystal

        lab_xyz = "_atom_site_fract"
        lab_beta = "_atom_site_aniso_U"
        lab_chi = "_atom_site_susceptibility"
        llabs = [lab_xyz, lab_beta, lab_chi]
        lnumb = [[] for hh in llabs]

        for iloop, loop in enumerate(data["loops"]):
            l_key = loop.keys()
            for hh, numb in zip(llabs, lnumb):
                flag = any([True if hh1.startswith(hh) else False for hh1 in l_key])
                if flag:
                    numb.append(iloop)
    
        l_atom_type = []
        if lnumb[0] != []:
            l_relation = data_atom_relation()
            for numb in lnumb[0]:
                l_key = list(data["loops"][numb].keys())
                n_atom = len(data["loops"][numb][l_key[0]])
                for i_atom in range(n_atom):
                    dd = {}
                    for key in l_key:
                        dd.update({key:data["loops"][numb][key][i_atom]})
                    atom_type = f_crystal.cl_crystal.AtomType()
                    f_api_rcif.api_rcif_common.from_dict_to_obj(dd, l_relation, atom_type)
                    l_atom_type.append(atom_type)
        else:
            print("Fractional coordinates of atoms are not found")
            return crystal
    
        if lnumb[1] != []:
            l_relation = data_adp_relation()
            for numb in lnumb[1]:
                l_key = list(data["loops"][numb].keys())
                n_atom = len(data["loops"][numb][l_key[0]])
                for i_atom in range(n_atom):
                    dd = {}
                    for key in l_key:
                        dd.update({key:data["loops"][numb][key][i_atom]})
                    for atom_type in l_atom_type:
                        if dd["_atom_site_aniso_label"] == atom_type.get_val("name"):
                            f_api_rcif.api_rcif_common.from_dict_to_obj(dd, l_relation, atom_type)
    
        if lnumb[2] != []:
            l_relation = data_chi_relation()
            for numb in lnumb[2]:
                l_key = list(data["loops"][numb].keys())
                n_atom = len(data["loops"][numb][l_key[0]])
                for i_atom in range(n_atom):
                    dd = {}
                    for key in l_key:
                        dd.update({key:data["loops"][numb][key][i_atom]})
                    for atom_type in l_atom_type:
                        if dd["_atom_site_magnetism_aniso_label"] == atom_type.get_val("name"):
                            f_api_rcif.api_rcif_common.from_dict_to_obj(dd, l_relation, atom_type)
                            atom_type.set_val(flag_m=True)

        for atom_type in l_atom_type:
            atom_site.add_atom(atom_type)
        l_crystal.append(crystal)
    return l_crystal



def conv_crystal_to_data(l_crystal):
    l_data = []
    def temp_func(ddata, relation, obj):
        lab_d, lab_o, type_val = relation[:3]
        key_d = lab_d
        
        if type_val == "logic":
            sval = "True" if obj.get_val(lab_o) else "False"
            ddata[key_d] = sval 
        elif type_val == "text":
            ddata[key_d] = obj.get_val(lab_o)
        elif type_val == "val":
            val = obj.get_val(lab_o)
            if isinstance(val, f_common.cl_variable.Variable):
                sval = "{:}".format(val.print_with_sigma())
            else:
                sval = "{:}".format(val)
            ddata[key_d] = sval 
        else:
            print(50*"ERROR ")
            print("type_val: ", type_val)
        return
        
    
    for obj in l_crystal:
        dd = {"name": obj.get_val("name")}
        l_relation = data_space_groupe_relation()
        space_groupe = obj.get_val("space_groupe")

        for relation in l_relation:
            temp_func(dd, relation, space_groupe)

        l_relation = data_cell_relation()
        cell = obj.get_val("cell")

        for relation in l_relation:
            temp_func(dd, relation, cell)
            
        lloop_d = []
        
        d_at = {}
        l_relation = data_atom_relation()

        atom_site = obj.get_val("atom_site")
        l_atom_type = atom_site._list_atom_type

        l_lab_d_used = []
        l_dhelp = []
        for obj_2 in l_atom_type:
            dhelp = {}
            for relation in l_relation:
                lab_d = relation[0]
                temp_func(dhelp, relation, obj_2)
                if (not(lab_d in l_lab_d_used)): 
                    l_lab_d_used.append(lab_d)
            l_dhelp.append(dhelp) 
        for lab_d in l_lab_d_used:
            d_at[lab_d] = [hh[lab_d] for hh in l_dhelp]
        lloop_d.append(d_at)
        
        d_bscat = {}
        l_relation = data_bscat_relation()

        l_lab_d_used = []
        l_dhelp = []
        for obj_2 in l_atom_type:
            dhelp = {}
            for relation in l_relation:
                lab_d = relation[0]
                temp_func(dhelp, relation, obj_2)
                if (not(lab_d in l_lab_d_used)):
                    l_lab_d_used.append(lab_d)
            l_dhelp.append(dhelp) 
        lab_uniq = "_atom_site_type_symbol"
        l_val = [hh[lab_uniq] for hh in l_dhelp]
        s_uniq = set(l_val)
        l_ind_uniq = [l_val.index(hh) for hh in s_uniq]
        for lab_d in l_lab_d_used:
            d_bscat[lab_d] = [l_dhelp[ind][lab_d] for ind in l_ind_uniq]
        lloop_d.append(d_bscat)

        d_adp = {}
        l_relation = data_adp_relation()
        l_lab_d_used = []
        l_dhelp = []
        for obj_2 in l_atom_type:
            dhelp = {}
            flag_beta = obj_2.get_val("adp_type") == "uani"
            if flag_beta:
                for relation in l_relation:
                    lab_d = relation[0]
                    temp_func(dhelp, relation, obj_2)
                    if (not(lab_d in l_lab_d_used)): 
                        l_lab_d_used.append(lab_d)
                l_dhelp.append(dhelp) 
        for lab_d in l_lab_d_used:
            d_adp[lab_d] = [hh[lab_d] for hh in l_dhelp]
        lloop_d.append(d_adp)
        
        d_chi = {}
        l_relation = data_chi_relation()
        l_lab_d_used = []
        l_dhelp = []
        for obj_2 in l_atom_type:
            dhelp = {}
            flag_m = obj_2.get_val("flag_m")
            if flag_m:
                for relation in l_relation:
                    lab_d = relation[0]
                    temp_func(dhelp, relation, obj_2)
                    if (not(lab_d in l_lab_d_used)): 
                        l_lab_d_used.append(lab_d)
                l_dhelp.append(dhelp)
        for lab_d in l_lab_d_used:
            d_chi[lab_d] = [hh[lab_d] for hh in l_dhelp]
        lloop_d.append(d_chi)
        
        if lloop_d != []:
            dd["loops"] = lloop_d
            
        l_data.append(dd)

    return l_data




def data_space_groupe_relation():
    l_relation = [ ("_space_group_name_H-M_alt", "spgr_given_name", "text"),
        ("_space_group_it_coordinate_system_code", "spgr_choice", "text")]
    return l_relation

def data_cell_relation():
    l_relation = [("_cell_length_a", "a", "val"), ("_cell_length_b", "b", "val"),
        ("_cell_length_c", "c", "val"), ("_cell_angle_alpha", "alpha", "val"),
        ("_cell_angle_beta", "beta", "val"), ("_cell_angle_gamma", "gamma", "val")]
    return l_relation

def data_atom_relation():
    l_relation = [("_atom_site_label", "name", "text"),
        ("_atom_site_type_symbol", "type_n", "text"),
        ("_atom_site_fract_x", "x", "val"),
        ("_atom_site_fract_y", "y", "val"),
        ("_atom_site_fract_z", "z", "val"),
        ("_atom_site_occupancy", "occupation", "val"),
        ("_atom_site_b_iso_or_equiv", "b_iso", "val"),
        ("_atom_site_adp_type", "adp_type", "text")]
    return l_relation

def data_bscat_relation():
    l_relation = [("_atom_site_type_symbol", "type_n", "text"),
        ("_atom_site_bscat", "b_scat", "val")]
    return l_relation

def data_adp_relation():
    l_relation = [("_atom_site_aniso_label", "name", "text"),
        ("_atom_site_aniso_U_11", "u_11", "val"),
        ("_atom_site_aniso_U_12", "u_12", "val"),
        ("_atom_site_aniso_U_13", "u_13", "val"),
        ("_atom_site_aniso_U_22", "u_22", "val"),
        ("_atom_site_aniso_U_23", "u_23", "val"),
        ("_atom_site_aniso_U_33", "u_33", "val")]
    return l_relation


def data_chi_relation():
    l_relation = [("_atom_site_magnetism_aniso_label", "name", "text"),
        ("_atom_site_magnetism_type_symbol", "type_m", "text"),
        ("_atom_site_magnetism_kappa", "kappa", "val"),
        ("_atom_site_magnetism_factor_lande", "factor_lande", "val"),
        ("_atom_site_magnetism_type", "chi_type", "text"),
        ("_atom_site_susceptibility_aniso_chi_11", "chi_11", "val"),
        ("_atom_site_susceptibility_aniso_chi_12", "chi_12", "val"),
        ("_atom_site_susceptibility_aniso_chi_13", "chi_13", "val"),
        ("_atom_site_susceptibility_aniso_chi_22", "chi_22", "val"),
        ("_atom_site_susceptibility_aniso_chi_23", "chi_23", "val"),
        ("_atom_site_susceptibility_aniso_chi_33", "chi_33", "val")]
    return l_relation



def data_extinction_single_domain_relation():
    l_relation = [("_sdd_phase_extinction_radius", "domain_radius", "val"),
        ("_sdd_phase_extinction_mosaicity", "mosaicity", "val")]
    return l_relation

def data_domain_single_domain_relation():
    l_relation = [("_sdd_scale_domain", "scale_domain", "list")]
    return l_relation




def main(larg):
    if len(larg) > 1:
        fname = "powder.rcif"
        fname = larg[1]
    else:
        fname = input("What is the name of the 'rcif' file?\n")
        
    
    

if __name__ == "__main__":
    main(sys.argv)
