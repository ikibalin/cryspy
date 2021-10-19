from .cl_3_data import DataN

def take_items_by_class(global_obj, l_class) -> list:
    l_res = []
    for item in global_obj.items:
        if isinstance(item, l_class):
            l_res.append(item)
        elif isinstance(item, DataN):
            l_res_data = take_items_by_class(item, l_class)
            l_res.extend(l_res_data)
    return l_res