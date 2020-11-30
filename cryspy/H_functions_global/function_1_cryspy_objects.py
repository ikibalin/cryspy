"""Transform data to objects of cryspy library."""
import os
import os.path
from typing import List

from pycifstar import Global, to_global, Items, Loop, Data

from cryspy.A_functions_base.function_1_strings import find_prefix

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN
from cryspy.B_parent_classes.cl_3_data import DataN
from cryspy.B_parent_classes.cl_4_global import GlobalN

from cryspy.C_item_loop_classes.cl_1_atom_site import AtomSite, \
    AtomSiteL
from cryspy.C_item_loop_classes.cl_1_atom_type import AtomType, \
    AtomTypeL
from cryspy.C_item_loop_classes.cl_1_atom_site_aniso import \
    AtomSiteAniso, AtomSiteAnisoL
from cryspy.C_item_loop_classes.cl_1_atom_site_moment import \
    AtomSiteMoment, AtomSiteMomentL
from cryspy.C_item_loop_classes.cl_1_atom_site_susceptibility import \
    AtomSiteSusceptibility, AtomSiteSusceptibilityL
from cryspy.C_item_loop_classes.cl_1_atom_type_scat import \
    AtomTypeScat, AtomTypeScatL
from cryspy.C_item_loop_classes.cl_1_atom_local_axes import \
    AtomLocalAxes, AtomLocalAxesL
from cryspy.C_item_loop_classes.cl_1_atom_electron_configuration \
    import AtomElectronConfiguration, AtomElectronConfigurationL
from cryspy.C_item_loop_classes.cl_1_cell import Cell, CellL
from cryspy.C_item_loop_classes.cl_1_chi2 import Chi2, Chi2L
from cryspy.C_item_loop_classes.cl_1_diffrn_radiation import \
    DiffrnRadiation, DiffrnRadiationL
from cryspy.C_item_loop_classes.cl_1_diffrn_refln import DiffrnRefln, \
    DiffrnReflnL
from cryspy.C_item_loop_classes.cl_1_exclude import Exclude, ExcludeL
from cryspy.C_item_loop_classes.cl_1_extinction import Extinction, \
    ExtinctionL
from cryspy.C_item_loop_classes.cl_1_inversed_hessian import InversedHessian
from cryspy.C_item_loop_classes.cl_1_mem_parameters import \
    MEMParameters, MEMParametersL
from cryspy.C_item_loop_classes.cl_1_phase import Phase, PhaseL
from cryspy.C_item_loop_classes.cl_1_pd2d_instr_reflex_asymmetry import\
    Pd2dInstrReflexAsymmetry, Pd2dInstrReflexAsymmetryL
from cryspy.C_item_loop_classes.cl_1_pd2d_instr_resolution import \
    Pd2dInstrResolution, Pd2dInstrResolutionL
from cryspy.C_item_loop_classes.cl_1_pd2d_meas import Pd2dMeas
from cryspy.C_item_loop_classes.cl_1_pd2d_proc import Pd2dProc
from cryspy.C_item_loop_classes.cl_1_pd2d_peak import Pd2dPeak, \
    Pd2dPeakL
from cryspy.C_item_loop_classes.cl_1_pd2d_background import \
    Pd2dBackground
from cryspy.C_item_loop_classes.cl_1_pd_background import \
    PdBackground, PdBackgroundL
from cryspy.C_item_loop_classes.cl_1_pd_instr_reflex_asymmetry import \
    PdInstrReflexAsymmetry, PdInstrReflexAsymmetryL
from cryspy.C_item_loop_classes.cl_1_pd_instr_resolution import \
    PdInstrResolution, PdInstrResolutionL
from cryspy.C_item_loop_classes.cl_1_pd_meas import PdMeas, PdMeasL
from cryspy.C_item_loop_classes.cl_1_pd_proc import PdProc, PdProcL
from cryspy.C_item_loop_classes.cl_1_pd_peak import PdPeak, PdPeakL
from cryspy.C_item_loop_classes.cl_1_range import Range, RangeL
from cryspy.C_item_loop_classes.cl_1_refln import Refln, ReflnL
from cryspy.C_item_loop_classes.cl_1_refln_susceptibility import \
    ReflnSusceptibility, ReflnSusceptibilityL
from cryspy.C_item_loop_classes.cl_1_refine_ls import RefineLs, \
    RefineLsL
from cryspy.C_item_loop_classes.cl_1_setup import Setup, SetupL
from cryspy.C_item_loop_classes.cl_1_space_group_symop_magn_centering import \
    SpaceGroupSymopMagnCentering, SpaceGroupSymopMagnCenteringL
from cryspy.C_item_loop_classes.cl_1_texture import Texture, TextureL
from cryspy.C_item_loop_classes.cl_1_tof_background import TOFBackground, \
    TOFBackgroundL
from cryspy.C_item_loop_classes.cl_1_tof_meas import TOFMeas, TOFMeasL
from cryspy.C_item_loop_classes.cl_1_tof_parameters import TOFParameters, \
    TOFParametersL
from cryspy.C_item_loop_classes.cl_1_tof_peak import TOFPeak, TOFPeakL
from cryspy.C_item_loop_classes.cl_1_tof_proc import TOFProc, TOFProcL
from cryspy.C_item_loop_classes.cl_1_tof_profile import TOFProfile, TOFProfileL

from cryspy.C_item_loop_classes.cl_2_atom_site_scat import \
    AtomSiteScat, AtomSiteScatL
from cryspy.C_item_loop_classes.cl_2_diffrn_orient_matrix import \
    DiffrnOrientMatrix, DiffrnOrientMatrixL
from cryspy.C_item_loop_classes.cl_2_section import Section, SectionL
from cryspy.C_item_loop_classes.cl_2_space_group import SpaceGroup
from cryspy.C_item_loop_classes.cl_2_space_group_symop_magn_operation import \
    SpaceGroupSymopMagnOperation, SpaceGroupSymopMagnOperationL

from cryspy.C_item_loop_classes.cl_3_density_point import DensityPoint, \
    DensityPointL

from cryspy.E_data_classes.cl_1_crystal import Crystal
from cryspy.E_data_classes.cl_1_mag_crystal import MagCrystal
from cryspy.E_data_classes.cl_2_diffrn import Diffrn
from cryspy.E_data_classes.cl_2_pd import Pd
from cryspy.E_data_classes.cl_2_pd2d import Pd2d
from cryspy.E_data_classes.cl_2_tof import TOF

from cryspy.G_global_classes.cl_1_rhochi import RhoChi
from cryspy.G_global_classes.cl_1_mem import MEM


L_ITEM_CLASS = []
L_LOOP_CLASS = []
L_DATA_CLASS = []
L_GLOBAL_CLASS = []
L_FUNCTION_ALL = []

L_ITEM_CLASS.extend([
    AtomElectronConfiguration, AtomLocalAxes,
    AtomSite, AtomSiteAniso, AtomSiteMoment, AtomSiteScat,
    AtomSiteSusceptibility, AtomType, AtomTypeScat,
    Cell, Chi2,
    DensityPoint, DiffrnRadiation, DiffrnOrientMatrix, DiffrnRefln,
    Exclude, Extinction, MEMParameters, InversedHessian,
    Phase, PdBackground,
    PdInstrReflexAsymmetry, PdInstrResolution, PdMeas, PdProc, PdPeak,
    Pd2dBackground, Pd2dInstrReflexAsymmetry,
    Pd2dInstrResolution, Pd2dMeas, Pd2dProc, Pd2dPeak,
    Range, RefineLs, Refln, ReflnSusceptibility,
    Setup, SpaceGroupSymopMagnCentering, SpaceGroupSymopMagnOperation,
    SpaceGroup, Section, Texture, TOFBackground, TOFMeas, TOFParameters,
    TOFPeak, TOFProc, TOFProfile])

L_LOOP_CLASS.extend([
    AtomLocalAxesL, AtomElectronConfigurationL,
    AtomSiteL, AtomSiteAnisoL, AtomSiteMomentL, AtomSiteScatL,
    AtomSiteSusceptibilityL, AtomTypeL, AtomTypeScatL,
    CellL, Chi2L,
    DensityPointL, DiffrnRadiationL, DiffrnOrientMatrixL, DiffrnReflnL,
    ExtinctionL, ExcludeL, MEMParametersL,
    PhaseL, PdBackgroundL, PdInstrReflexAsymmetryL, PdInstrResolutionL,
    PdMeasL, PdProcL, PdPeakL,
    Pd2dInstrReflexAsymmetryL, Pd2dInstrResolutionL, Pd2dPeakL,
    RangeL, RefineLsL, ReflnL, ReflnSusceptibilityL,
    SectionL, SetupL, SpaceGroupSymopMagnCenteringL,
    SpaceGroupSymopMagnOperationL, TextureL, TOFBackgroundL, TOFMeasL,
    TOFParametersL, TOFPeakL, TOFProcL, TOFProfileL])

L_DATA_CLASS.extend([Crystal, MagCrystal, Diffrn, Pd, Pd2d, TOF])

L_GLOBAL_CLASS.extend([MEM, RhoChi])


def is_in_rcif_block(rcif_block, cryspy_obj):
    """Is in rcif block."""
    ls_cryspy = cryspy_obj.to_cif.split("\n")
    l_name = [_.strip().split()[0] for _ in ls_cryspy if _.startswith("_")]
    l_flag = [rcif_block.is_value(_name) for _name in l_name]
    return any(l_flag), all(l_flag)


def file_to_globaln(file_name: str, item_classes=(), loop_classes=(),
                    data_classes=(), global_classes=()) -> GlobalN:
    """Transfer file content of cif format to cryspy objects."""
    if not(os.path.isfile(file_name)):
        raise UserWarning(f"File '{file_name:}' is not found")
        return
    s_cont = str(to_global(file_name))
    res = str_to_globaln(s_cont, item_classes=item_classes,
                         loop_classes=loop_classes, data_classes=data_classes,
                         global_classes=global_classes)
    return res


def str_to_items(s_cont: str, item_classes=(), loop_classes=()) -> list:
    """Transfer string of cif format to cryspy objects."""
    data_cif = Data()
    data_cif.take_from_string(s_cont)
    l_item_class = L_ITEM_CLASS + list(item_classes)
    l_loop_class = L_LOOP_CLASS + list(loop_classes)

    l_item = items_to_itemsn(data_cif.items, l_item_class)
    for loop_cif in data_cif.loops:
        l_loop_obj = loop_to_loopn(loop_cif, l_loop_class)
        l_item.extend(l_loop_obj)
    return l_item


def str_to_globaln(s_cont: str, item_classes=(), loop_classes=(),
                   data_classes=(), global_classes=()) -> GlobalN:
    """Transfer string of cif format to cryspy objects."""
    global_cif = Global()
    global_cif.take_from_string(s_cont)
    l_global_item = []
    l_cls_global = []
    l_item_class = L_ITEM_CLASS + list(item_classes)
    l_loop_class = L_LOOP_CLASS + list(loop_classes)
    l_data_class = L_DATA_CLASS + list(data_classes)
    l_global_class = L_GLOBAL_CLASS + list(global_classes)
    if len(global_cif.items)+len(global_cif.loops) > 0:
        l_item_obj = items_to_itemsn(global_cif.items, l_item_class)
        l_global_item.extend(l_item_obj)

        for item in l_item_obj:
            if type(item) not in l_cls_global:
                l_cls_global.append(type(item))

        for loop_cif in global_cif.loops:
            l_loop_obj = loop_to_loopn(loop_cif, l_loop_class)
            l_global_item.extend(l_loop_obj)
            for loop_obj in l_loop_obj:
                if type(loop_obj) not in l_cls_global:
                    l_cls_global.append(type(loop_obj))

    for data_cif in global_cif.datas:
        flag = False
        str_data = str(data_cif)
        for cls_data in l_data_class:
            # FIXME: it's bad solution as we loose information which are
            #        not specified as knonw in cryspy library.
            data_obj = cls_data.from_cif(str_data)
            if data_obj is not None:
                l_global_item.append(data_obj)
                if not(cls_data in l_cls_global):
                    l_cls_global.append(cls_data)
                flag = True
        if not(flag):
            l_data_item = []
            l_cls_data = []

            l_item_obj = items_to_itemsn(data_cif.items, l_item_class)
            l_data_item.extend(l_item_obj)
            for item in l_item_obj:
                if type(item) not in l_cls_data:
                    l_cls_data.append(type(item))

            for loop_cif in data_cif.loops:
                l_loop_obj = loop_to_loopn(loop_cif, l_loop_class)
                l_data_item.extend(l_loop_obj)
                for loop_obj in l_loop_obj:
                    if type(loop_obj) not in l_cls_data:
                        l_cls_data.append(type(loop_obj))

            if len(l_data_item) > 0:
                data_obj = items_to_datan(data_cif.name, l_data_item,
                                          l_data_class)
                l_global_item.append(data_obj)

    global_obj = items_to_globaln(global_cif.name, l_global_item,
                                  l_global_class)
    return global_obj


def find_name_prefix(loop_name, l_item_cls):
    """Find name prefix."""
    loop_name = loop_name.strip("_")
    l_loop_prefix = [item_cls.PREFIX for item_cls in l_item_cls]
    l_loop = [loop_prefix for loop_prefix in l_loop_prefix if
              loop_name.startswith(loop_prefix)]
    if len(l_loop) > 0:
        l_loop.sort(key=lambda x: len(x))
        for prefix in l_loop:
            name = loop_name.strip("_")[len(prefix)+1:]
            item_cls = l_item_cls[l_loop_prefix.index(prefix)]
            l_attr_cif = [attr.lower() for attr in item_cls.ATTR_CIF]
            if name.lower() in l_attr_cif:
                separator = loop_name.strip("_")[len(prefix)]
                return prefix, name, separator
        # prefix = l_loop[-1]
        # name = loop_name.strip("_")[len(prefix)+1:]
        # separator = loop_name.strip("_")[len(prefix)]
        # return prefix, name, separator
    return None, None, None


# TODO: l_item_class should be added as not necessary l_loop_class fully
#       correspond to l_item_class
def loop_to_loopn(loop_cif: Loop, l_loop_class: list) -> LoopN:
    """Transform Loop object to ItemN objects."""
    # l_item_class for future.
    l_item_class = []
    l_item_prefix = [item_cls.PREFIX for item_cls in l_item_class]
    l_item_cls_in_def_loop = [loop_cls.ITEM_CLASS for loop_cls in l_loop_class]
    l_loop_prefix = [loop_cls.ITEM_CLASS.PREFIX for loop_cls in l_loop_class]
    l_loop_name = [name.lower() for name in loop_cif.names]
    flag_id = False
    flag_label = False
    flag_type = False
    l_loop_prefix_in = []
    l_loop_name_in = []
    l_loop_separator_in = []

    for loop_name in l_loop_name:
        prefix, name, separator = find_name_prefix(loop_name,
                                                   l_item_cls_in_def_loop)

        if prefix is not None:
            l_name_allowed = [h.lower() for h in l_loop_class[
                l_loop_prefix.index(prefix)].ITEM_CLASS.ATTR_CIF]
            if not(name.lower() in l_name_allowed):
                prefix = None
                name = None

        if (prefix is None) | (name is None):
            prefix, name, separator = find_name_prefix(loop_name,
                                                       l_item_cls_in_def_loop)
            if prefix is not None:
                l_name_allowed = [h.lower() for h in l_item_class[
                    l_item_prefix.index(prefix)].ATTR_CIF]
                if not(name.lower() in l_name_allowed):
                    prefix = None
                    name = None

        if (prefix is None) | (name is None):
            if loop_name.find(".") != -1:
                prefix, name = loop_name.strip("_").split(".")[:2]
                separator = "."
            else:
                prefix = loop_name.strip("_").split("_")[0]
                name = loop_name.strip("_")[len(prefix)+1:]
                separator = "_"

        if name.lower() == "id":
            flag_id = True
            prefix_id = prefix
            loop_name_id = loop_name

        if name.lower() == "label":
            flag_label = True
            prefix_label = prefix
            loop_name_label = loop_name

        if name.lower() == "type":
            flag_type = True
            prefix_type = prefix
            loop_name_type = loop_name

        l_loop_prefix_in.append(prefix)
        l_loop_name_in.append(name)
        l_loop_separator_in.append(separator)


    s_loop_prefix_in = set(l_loop_prefix_in)
    l_loopn = []
    for prefix in s_loop_prefix_in:
        item_cls = ItemN
        loop_cls = LoopN
        if prefix in l_loop_prefix:
            loop_cls = l_loop_class[l_loop_prefix.index(prefix)]
        elif prefix in l_item_prefix:
            # for future
            item_cls = l_item_class[l_loop_prefix.index(prefix)]

        l_loop_name_cif = []
        ll_val_cif = []
        separator = "_"
        for loop_name, prefix_in, name_in, separator_in in zip(
                l_loop_name, l_loop_prefix_in, l_loop_name_in,
                l_loop_separator_in):
            if prefix == prefix_in:
                l_loop_name_cif.append(loop_name)
                ll_val_cif.append(loop_cif[loop_name])
                separator = separator_in
        if flag_id:
            if prefix_id != prefix:
                loop_new_name_id = "_"+prefix+f"{separator:}id"
                l_loop_name_cif.append(loop_new_name_id)
                ll_val_cif.append(loop_cif[loop_name_id])
        elif flag_label:
            if prefix_label != prefix:
                loop_new_name_label = "_"+prefix+f"{separator:}label"
                l_loop_name_cif.append(loop_new_name_label)
                ll_val_cif.append(loop_cif[loop_name_label])
        elif flag_type:
            if prefix_type != prefix:
                loop_new_name_type = "_"+prefix+f"{separator:}type"
                l_loop_name_cif.append(loop_new_name_type)
                ll_val_cif.append(loop_cif[loop_name_type])
        ll_val_cif_t = [[ll_val_cif[i_2][i_1] for i_2 in range(len(
            ll_val_cif))] for i_1 in range(len(ll_val_cif[0]))]
        obj = Loop(names=l_loop_name_cif, values=ll_val_cif_t)
        loopn = loop_cls.from_cif(str(obj))

        if loopn is None:
            loop_cls = LoopN
            loopn = loop_cls.from_cif(str(obj))
        # if len(loopn.items) == 0:
        #     print("obj:", loop_cls)
        #     print(obj)
        if not(loopn is None):
            l_loopn.append(loopn)
        else:
            print("l_loop_name_cif: \n", l_loop_name_cif)
    return l_loopn


def items_to_itemsn(items_cif: Items, l_item_class: list) -> List[ItemN]:
    """Transform Items object to items of ItemN classes."""

    l_item = []
    l_name_cif = list(items_cif.names)
    flag_point_separator = all([name.find(".") == -1
                                for name in l_name_cif])
    if flag_point_separator:
        separator = "_"
    else:
        separator = "."

    for cls_item in l_item_class:
        prefix = cls_item.PREFIX

        separator_loc = "_"

        l_cls_mand_cif = [f"_{prefix:}{separator_loc:}{attr_cif:}"
                          for attr_cif in cls_item.ATTR_MANDATORY_CIF]

        l_cls_all_cif = [f"_{prefix:}{separator_loc:}{attr_cif:}"
                         for attr_cif in (cls_item.ATTR_MANDATORY_CIF +
                                          cls_item.ATTR_OPTIONAL_CIF)]

        flag_mand = all([items_cif.is_value(name) for name in l_cls_mand_cif])
        flag_any = any([items_cif.is_value(name) for name in l_cls_all_cif])

        if not(flag_any):
            separator_loc = "."
            l_cls_mand_cif = [f"_{prefix:}{separator_loc:}{attr_cif:}"
                              for attr_cif in cls_item.ATTR_MANDATORY_CIF]

            l_cls_all_cif = [f"_{prefix:}{separator_loc:}{attr_cif:}"
                             for attr_cif in (cls_item.ATTR_MANDATORY_CIF +
                                              cls_item.ATTR_OPTIONAL_CIF)]

            flag_mand = all([items_cif.is_value(name) for name in
                             l_cls_mand_cif])
            flag_any = any([items_cif.is_value(name) for name in
                            l_cls_all_cif])

        item_obj = None
        if flag_mand:
            items_small = items_cif.items_with_prefix(f"_{prefix:}")
            item_obj = cls_item.from_cif(str(items_small))
        elif flag_any:
            items_small = items_cif.items_with_prefix(f"_{prefix:}")
            item_obj = ItemN.from_cif(str(items_small))

        if item_obj is not None:
            l_item.append(item_obj)
            l_name_cif_remove = [name.lower() for name in l_cls_all_cif
                                 if bool(l_name_cif.count(name.lower()))]
            for name in l_name_cif_remove:
                l_name_cif.remove(name)

    if len(l_name_cif) != 0:
        if separator == ".":
            l_prefix_cif = [(name[1:]).split(separator)[0]
                            for name in l_name_cif]
        else:
            l_prefix_cif = find_prefixes(l_name_cif)

        for prefix in l_prefix_cif:
            names = [name_cif for name_cif in l_name_cif if 
                     name_cif.startswith("_"+prefix)]
            items_small = items_cif.items_with_names(names)
            item_obj = ItemN.from_cif(str(items_small))
            l_item.append(item_obj)
    return l_item


def find_prefixes(l_name_cif):
    """Find prefix."""
    l_short_prefix_cif = list(set([(name[1:]).split("_")[0]
                                   for name in l_name_cif]))
    l_prefix_cif = []
    for short_prefix_cif in l_short_prefix_cif:
        l_name_2_cif = [name for name in l_name_cif
                        if (name[1:]).split("_")[0] == short_prefix_cif]
        prefix_cif = find_prefix(*tuple(l_name_2_cif))
        l_prefix_cif.append(prefix_cif[1:])
    return l_prefix_cif


def items_to_globaln(global_name: str, items: list, l_global_class):
    """
    Items to suitable GlobalN class.

    Parameters
    ----------
    global_name : str
        name of the GlobalN.
    items : List[Union[ItemN, LoopN, DataN]]
        List of items.
    l_global_class : list
        List of GlobalN classes to choose one of them.

    Returns
    -------
    global_obj : TYPE
        Object of GlobalN class or one of l_global_classes.

    """
    global_obj = None
    classes = set([type(item) for item in items])

    for cls_global in l_global_class:
        flag_mand = all([cls_mand in classes
                         for cls_mand in cls_global.CLASSES_MANDATORY])
        flag_opt = all([class_ in cls_global.CLASSES
                        for class_ in classes])
        if (flag_mand & flag_opt):
            global_obj = cls_global()
            global_obj.global_name = global_name
            global_obj.add_items(items)

    if global_obj is None:
        global_obj = GlobalN.make_container((), tuple(classes), "global")
        global_obj.global_name = global_name
        global_obj.add_items(items)
    return global_obj


def items_to_datan(data_name: str, items: list, l_data_class):
    """
    Items to suitable GlobalN class.

    Parameters
    ----------
    data_name : str
        name of the DataN.
    items : List[Union[ItemN, LoopN]]
        List of items.
    l_data_class : list
        List of DataN classes to choose one of them.

    Returns
    -------
    data_obj : TYPE
        Object of DataN class or one of l_data_classes.

    """
    data_obj = None
    classes = set([type(item) for item in items])
    for cls_data in l_data_class:
        flag_mand = all([cls_mand in classes
                         for cls_mand in cls_data.CLASSES_MANDATORY])
        flag_opt = all([class_ in cls_data.CLASSES
                        for class_ in classes])
        if (flag_mand & flag_opt):
            data_obj = cls_data()
            data_obj.data_name = data_name
            data_obj.add_items(items)

    if data_obj is None:
        data_obj = DataN.make_container((), tuple(classes), "data")
        data_obj.data_name = data_name
        data_obj.add_items(items)
    return data_obj

# s_cont = """
# global_
#   _unknown_fuck_d 7
#   _unknown_fuck_e 3546
#   _unknown_fuck_type_g 14
#   _test_fuck_d 7
#   _test_fuck_e 3546
#   _test_fuck_type_g 14

#   data_Fe3O4

#   _test2_fuck_d 7
#   _test2_fuck_e 3546
#   _test2_fuck_type_g 14

#   _cell_angle_alpha 90.0
#   _cell_angle_beta 90.0
#   _cell_angle_gamma 90.0
#   _cell_length_a 8.56212()
#   _cell_length_b 8.56212
#   _cell_length_c 8.56212
#   _space_group_it_coordinate_system_code 2
#   _space_group_IT_number    227
#   loop_
#   _atom_site_adp_type
#   _atom_site_B_iso_or_equiv
#   _atom_site_fract_x
#   _atom_site_fract_y
#   _atom_site_fract_z
#   _atom_site_label
#   _atom_site_occupancy
#   _atom_site_type_symbol
#   Uani 0.0 0.125 0.125 0.125 Fe3A 1.0 Fe3+
#   Uani 0.0 0.5 0.5 0.5 Fe3B 1.0 Fe3+
#   Uiso 0.0 0.25521 0.25521 0.25521 O1 1.0 O2-
#   loop_
#   _atom_type_scat_length_neutron
#   _atom_type_symbol
#     0.945 Fe3+
#   0.5803 O2-
#   loop_
#   _atom_site_aniso_U_11
#   _atom_site_aniso_U_12
#   _atom_site_aniso_U_13
#   _atom_site_aniso_U_22
#   _atom_site_aniso_U_23
#   _atom_site_aniso_U_33
#   _atom_site_aniso_label
#   0.0 0.0 0.0 0.0 0.0 0.0 Fe3A
#   0.0 0.0 0.0 0.0 0.0 0.0 Fe3B
#   loop_
#   _atom_site_scat_label
#   _atom_site_scat_lande
#   Fe3A 2.0
#   Fe3B 2.0
#   loop_
#   _atom_site_susceptibility_label
#   _atom_site_susceptibility_chi_type
#   _atom_site_susceptibility_chi_11
#   _atom_site_susceptibility_chi_12
#   _atom_site_susceptibility_chi_13
#   _atom_site_susceptibility_chi_22
#   _atom_site_susceptibility_chi_23
#   _atom_site_susceptibility_chi_33
#   Fe3A Cani -3.468(74) 0.0 0.0 -3.468 0.0 -3.468
#   Fe3B Cani 3.041      0.0 0.0  3.041 0.0  3.041
#   data_mono
#   _setup_wavelength     0.840
#   _setup_field          1.000
#   _diffrn_radiation_polarization 1.0
#   _diffrn_radiation_efficiency   1.0
#   _extinction_mosaicity 100.0
#   _extinction_radius    50.0
#   _extinction_model     gauss
#   _diffrn_orient_matrix_UB_11 6.59783
#   _diffrn_orient_matrix_UB_12 -6.99807
#   _diffrn_orient_matrix_UB_13 3.3663
#   _diffrn_orient_matrix_UB_21 2.18396
#   _diffrn_orient_matrix_UB_22 -2.60871
#   _diffrn_orient_matrix_UB_23 -9.5302
#   _diffrn_orient_matrix_UB_31 7.4657
#   _diffrn_orient_matrix_UB_32 6.94702
#   _diffrn_orient_matrix_UB_33 -0.18685
#   _phase_label  Fe3O4
#   loop_
#   _diffrn_refln_index_h
#   _diffrn_refln_index_k
#   _diffrn_refln_index_l
#   _diffrn_refln_fr
#   _diffrn_refln_fr_sigma
#   0 0 8 0.64545 0.01329
#   2 0 6 1.75682 0.04540
#   0 2 6 1.67974 0.03711
# """
# obj = str_to_globaln(s_cont)
# for var_name in obj.get_variable_names():
#     print(var_name)
# obj
