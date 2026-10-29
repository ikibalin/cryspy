







from cryspy.E_data_classes.cl_1_crystal import Crystal
from cryspy.E_data_classes.cl_2_diffrn import Diffrn



def form_dictionary_by_crystal(data_obj) -> dict:
    """
    Based on cryspy objects form dictionary which contents whole
    information about crystal structure
    """
    if isinstance(data_obj, Crystal):
        ddict = data_obj.get_dictionary()
    else:
        ddict = {}

    return ddict



def form_dictionary_by_diffrn(data_obj) -> dict:
    """
    Based on cryspy objects form dictionary from numpy arrays needed for
    calculations

    Parameters
    ----------
    global_obj : GlobalN of DataN containers
        data container with needed information.

    Returns
    -------
    ddict : dictionary
        output keys: .

    """
    if isinstance(data_obj, Diffrn):
        ddict = data_obj.get_dictionary()
    else:
        ddict = {}

    return ddict



