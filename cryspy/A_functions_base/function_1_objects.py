"""Module realized some operations with python objects

"""
from types import FunctionType


def get_functions_of_objet(obj):
    """Get list of object function.
    """
    ls_out = []
    l_method = [_1 for _1, _2 in type(obj).__dict__.items()
                if ((type(_2) == FunctionType) &
                    (not(_1.startswith("_"))))]
    for method in l_method:
        func = getattr(obj, method)
        l_param = [_ for _ in
                   func.__code__.co_varnames[:func.__code__.co_argcount]
                   if _ != "self"]
        s_par = ""
        if len(l_param) > 0:
            s_par = ", ".join(l_param)
        s_val = f"{method:}({s_par:})"
        ls_out.append(s_val)
    return sorted(ls_out)


def variable_name_to_string(variable_name: tuple, flag_html: bool = True):
    """Give string presentation of variable name.
    
    Variable name: 
        ((name, val), (name, val), ...)
    """
    ls_out = []
    for name_val in variable_name[1:]:
        name, val = name_val
        if ((val is None) | (val == "")):
            ls_out.append(f"{name_val[0]:}")
        elif isinstance(name_val[1], int):
            ls_out.append(f"{name_val[0]:}[{name_val[1]:}]")
        else:
            ls_out.append(f"{name_val[0]:}_{name_val[1]:}")
    if flag_html:
        ls_out[-1] = "<b>"+ls_out[-1]+"</b>"
    return '.'+'.'.join(ls_out)


def change_variable_name(variable_name: tuple, ending: str = ""):
    """Change the ending of variable name.
    
    Output:
        modified variable name
    """
    if ending == "":
        res = variable_name
    elif len(variable_name[-1]) == 1:
        res = variable_name[:-1] + ((f"{variable_name[-1][0]:}_{ending:}", ), )
    else:
        res = variable_name[:-1] + ((f"{variable_name[-1][0]:}_{ending:}", 
                                     variable_name[-1][1]), )
    return res


def get_table_html_for_variables(obj):
    """Get html table for variables of object
    """
    ls_html=["<table border='1'>"]
    variable_names = obj.get_variable_names()
    ls_variables = []
    if len(variable_names) != 0:
        ls_variables.append("<th>Variable</th><th>Value</th><th>Sigma</th>")
        ls_variables.extend([f"<td>{variable_name_to_string(var_name):} </td> \
<td>{obj.get_variable_by_name(var_name):}</td>\
<td>{obj.get_variable_by_name(change_variable_name(var_name, 'sigma')):}</td>"
for var_name in variable_names])

    for variables in ls_variables:
        ls_html.append(f"<tr>{variables}</tr>")
    ls_html.append("</table>")

    return " ".join(ls_html)


def form_items_by_dictionary(item_class, d_items):
    """
    Form items for loop when variables are given by dictionary
    """
    items = []
    if "items" in d_items.items():
        items  = d_items["items"]
    elif len(d_items.items()) > 0:
        i=0
        for hh in zip(*d_items.items()):
            if i==0:
                keys = hh
                i += 1
            else:
                vals = hh
        
        for val in zip(*vals):
            d_item = {}
            for k,v in zip(keys, val):
                d_item[k]=v
            item = item_class(**d_item)
            items.append(item)
        
    return items