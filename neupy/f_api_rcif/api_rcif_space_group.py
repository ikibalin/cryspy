"""
transform inforamtion from cif file to class SpaceGroup
"""
import sys
import os


from neupy.f_common.cl_fitable import Fitable

from neupy.f_crystal import (
    SpaceGroup
    )


from neupy.f_api_rcif.api_rcif_common import (
    from_dict_to_obj,
    conv_str_to_text_float_logic
    )





def conv_data_to_space_group(l_data):
    """
    transform info from list of dictionary intolist of SpaceGroup class
    """
    l_flag_space_group = []
    for data in l_data:
        l_relation = data_space_group_relation()
        l_keys = data.keys()
        flag_crystal = any(relation[0] in l_keys for relation in l_relation)
        l_flag_space_group.append(flag_crystal)

    l_data_space_group = [data for data, flag in zip(l_data, l_flag_space_group) if flag]
    l_space_group = []
    for data in l_data_space_group:
        space_group = SpaceGroup()
        l_relation = data_space_group_relation()
        from_dict_to_obj(data, l_relation, space_group)
        l_space_group.append(space_group)
    return l_space_group



def conv_space_group_to_data(l_space_group):
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
            if isinstance(val, Variable):
                sval = "{:}".format(val.print_with_sigma())
            else:
                sval = "{:}".format(val)
            ddata[key_d] = sval 
        else:
            print(50*"ERROR ")
            print("type_val: ", type_val)
        return
        
    
    for obj in l_space_group:
        dd = {"name": obj.get_val("name")}
        l_relation = data_space_group_relation()
        space_group = obj.get_val("space_group")

        for relation in l_relation:
            temp_func(dd, relation, space_group)

        l_data.append(dd)
    return l_data


def data_space_group_relation():
    l_relation = [ ("_space_group_name_H-M_alt", "spgr_name", "text"),
        ("_space_group_it_coordinate_system_code", "spgr_choice", "text")]
    return l_relation

