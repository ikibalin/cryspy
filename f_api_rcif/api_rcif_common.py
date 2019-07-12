"""
transform inforamtion from cif file to class phase
"""
import sys
import os





def from_dict_to_obj(dict_i, l_relation, obj):
    """
    l_relation is list of (lab_d, lab_o, val_type)
    """
    
    l_key = list(dict_i.keys())
    l_key_lower = [hh.lower() for hh in l_key]
    l_numb = [ihh for ihh, hh in enumerate(l_relation) if (hh[0].lower() in l_key_lower)]
    d_args = {}
    for numb in l_numb:
        lab_d, lab_o, val_type = l_relation[numb]
        key_d = l_key[l_key_lower.index(lab_d.lower())]
        if val_type == "val":
            #val = [dict_i[llab_d[numb]], False, ""]
            val = dict_i[key_d]
            val = conv_str_to_text_float_logic(val, lab_o)
        elif val_type == "text":
            val = dict_i[key_d]
        elif val_type == "list":
            val = dict_i[key_d]
        elif val_type == "logic":
            val = dict_i[key_d]
            val = conv_str_to_text_float_logic(val, lab_o)
        else:
            print("mistake in type variable of 'from_dict_to_obj' function")
            val = None
        d_args.update({lab_o: val})
    obj.set_val(**d_args)
    return





def conv_str_to_text_float_logic(sval, name=""):
    if (not (isinstance(sval,str))):
        return sval
    try:
        if len(sval.strip().split("("))>1:
            l_help = sval.split("(")
            val_1 = float(l_help[0])
            val = cl_variable.Variable(val_1, True, name)
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



def main(larg):
    if len(larg) > 1:
        fname = "powder.rcif"
        fname = larg[1]
    else:
        fname = input("What is the name of the 'rcif' file?\n")
        
    
    

if __name__ == "__main__":
    main(sys.argv)
