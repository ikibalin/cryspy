"""Space Group class."""

from typing import NoReturn
import numpy
from fractions import Fraction

from cryspy.A_functions_base.function_1_objects import \
    form_items_by_dictionary
from cryspy.A_functions_base.function_2_space_group import \
    get_shift_by_centring_type, \
    get_symop_pcentr_multiplicity_letter_site_symmetry_coords_xyz_2, \
    get_type_hm, auto_choose_it_coordinate_system_code, \
    get_it_number_it_coordinate_system_codes_by_name_hm_extended,\
    get_it_number_by_name_hm_full, \
    get_it_number_by_name_hm_short, \
    get_it_number_by_name_hall, \
    get_it_number_by_name_schoenflies, \
    get_it_coordinate_system_codes_by_it_number, \
    get_default_it_coordinate_system_code_by_it_number, \
    get_crystal_system_by_it_number, \
    get_name_hm_extended_by_it_number_it_coordinate_system_code, \
    get_centring_type_by_name_hm_extended, \
    get_bravais_type_by_centring_type_crystal_system, \
    get_name_hm_short_by_it_number, \
    get_lattice_type_by_name_hm_short, \
    get_name_hm_full_by_it_number, \
    get_name_hall_by_it_number, \
    get_centrosymmetry_by_name_hall, \
    get_name_schoenflies_by_it_number, \
    get_laue_class_by_name_schoenflies, \
    get_point_group_hm_short_by_name_schoenflies, \
    get_generators_by_point_group_hm, \
    get_patterson_name_hm_by_lattice_type_laue_class, \
    ACCESIBLE_IT_NUMBER, ACCESIBLE_BRAVAIS_TYPE, \
    ACCESIBLE_IT_COORDINATE_SYSTEM_CODE, ACCESIBLE_LAUE_CLASS, \
    ACCESIBLE_CENTRING_TYPE, ACCESIBLE_CRYSTAL_SYSTEM, \
    ACCESIBLE_NAME_HM_SHORT, ACCESIBLE_NAME_HM_FULL, \
    ACCESIBLE_NAME_HM_EXTENDED, ACCESIBLE_NAME_SCHOENFLIES, \
    ACCESIBLE_NAME_HALL_SHORT, ACCESIBLE_REFERENCE_SETTING

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN

from cryspy.C_item_loop_classes.cl_1_space_group_symop import \
    SpaceGroupSymop, SpaceGroupSymopL

from cryspy.C_item_loop_classes.cl_1_space_group_wyckoff import \
    SpaceGroupWyckoff, SpaceGroupWyckoffL


class SpaceGroup(ItemN):
    """Space group description.

    Contains all the data items that refer to the space group as a
    whole, such as its name, Laue group etc.

    Space-group types are identified by their number as listed in
    International Tables for Crystallography Volume A, or by their
    Schoenflies symbol. Specific settings of the space groups can
    be identified by their Hall symbol, by specifying their
    symmetry operations or generators, or by giving the
    transformation that relates the specific setting to the
    reference setting based on International Tables Volume A and
    stored in this dictionary.

    The commonly used Hermann-Mauguin symbol determines the
    space-group type uniquely but several different Hermann-Mauguin
    symbols may refer to the same space-group type. A
    Hermann-Mauguin symbol contains information on the choice of
    the basis, but not on the choice of origin.

    """

    ATTR_MANDATORY_NAMES = ()
    ATTR_MANDATORY_TYPES = ()
    ATTR_MANDATORY_CIF = ()

    ATTR_OPTIONAL_NAMES = (
        "id", "name_hm_alt", "name_hm_alt_description", "it_number",
        "name_hm_ref", "it_coordinate_system_code", "name_hm_full",
        "name_hall", "name_schoenflies", "point_group_hm", "laue_class",
        "patterson_name_hm", "centring_type", "bravais_type", "crystal_system",
        "reference_setting", "transform_pp_abc", "transform_qq_xyz")
    ATTR_OPTIONAL_TYPES = (str, str, str, int, str, str, str, str, str, str,
                           str, str, str, str, str, str, str, str)
    ATTR_OPTIONAL_CIF = (
        "id", "name_H-M_alt", "name_H-M_alt_description", "IT_number",
        "name_H-M_ref", "IT_coordinate_system_code", "name_H-M_full",
        "name_Hall", "name_Schoenflies", "point_group_H-M", "Laue_class",
        "Patterson_name_H-M", "centring_type", "Bravais_type",
        "crystal_system", "reference_setting", "transform_pp_abc",
        "transform_qq_xyz")

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ("centrosymmetry", "pcentr", "reduced_space_group_symop",
                      "full_space_group_symop", "space_group_wyckoff", "shift")
    ATTR_INT_PROTECTED_NAMES = ()

    # parameters considered are refined parameters
    ATTR_REF = ()
    ATTR_SIGMA = tuple([f"{_h:}_sigma" for _h in ATTR_REF])
    ATTR_CONSTR_FLAG = tuple([f"{_h:}_constraint" for _h in ATTR_REF])
    ATTR_REF_FLAG = tuple([f"{_h:}_refinement" for _h in ATTR_REF])
    ATTR_CONSTR_MARK = tuple([f"{_h:}_mark" for _h in ATTR_REF])

    ACCESIBLE_NAME_HM_ALT = frozenset(
        set(ACCESIBLE_NAME_HM_SHORT) | set(ACCESIBLE_NAME_HM_FULL) |
        set(ACCESIBLE_NAME_HM_EXTENDED))

    # constraints on the parameters
    D_CONSTRAINTS = {
        "it_number": ACCESIBLE_IT_NUMBER,
        "bravais_type": ACCESIBLE_BRAVAIS_TYPE,
        "it_coordinate_system_code": ACCESIBLE_IT_COORDINATE_SYSTEM_CODE,
        "laue_class": ACCESIBLE_LAUE_CLASS,
        "centring_type": ACCESIBLE_CENTRING_TYPE,
        "crystal_system": ACCESIBLE_CRYSTAL_SYSTEM,
        "name_hm_alt": ACCESIBLE_NAME_HM_ALT,
        "name_hm_ref": ACCESIBLE_NAME_HM_SHORT,
        "name_schoenflies": ACCESIBLE_NAME_SCHOENFLIES,
        "name_hall": ACCESIBLE_NAME_HALL_SHORT,
        "reference_setting": ACCESIBLE_REFERENCE_SETTING,
        "name_hm_full": ACCESIBLE_NAME_HM_FULL}

    # default values for the parameters
    D_DEFAULT = {}
    for key in ATTR_SIGMA:
        D_DEFAULT[key] = 0.
    for key in (ATTR_CONSTR_FLAG + ATTR_REF_FLAG):
        D_DEFAULT[key] = False
    for key in ATTR_CONSTR_MARK:
        D_DEFAULT[key] = ""

    PREFIX = "space_group"

    def __init__(self, **kwargs) -> NoReturn:
        super(SpaceGroup, self).__init__()

        # defined for any integer and float parameters
        D_MIN = {}

        # defined for ani integer and float parameters
        D_MAX = {}

        self.__dict__["D_MIN"] = D_MIN
        self.__dict__["D_MAX"] = D_MAX
        flag_form_object = False
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)
            if key in ("it_number", "name_hm_alt", "name_hm_ref",
                       "name_schoenflies", "name_hall", "name_hm_full"):
                flag_form_object = True

        if flag_form_object:
            self.form_object()

    def form_object(self) -> NoReturn:
        """Form object."""
        self.form_object_by_it_number_it_coordinate_system_code()
        self.__dict__["shift"] = get_shift_by_centring_type(self.centring_type)

    def define_it_number_and_it_coordinate_system_code(self):
        """Define it number and it coordinate system code."""
        flag = True
        it_coordinate_system_code = None
        if self.is_attribute("it_coordinate_system_code"):
            it_coordinate_system_code = self.it_coordinate_system_code
        
        if (self.is_attribute("it_number") &
                self.is_attribute("it_coordinate_system_code")):
            it_number = self.it_number
            it_coordinate_system_code = self.it_coordinate_system_code
        elif ((self.is_attribute("name_hm_alt")) &
              (it_coordinate_system_code is not None)):
            t_res = get_type_hm(self.name_hm_alt)
            if "full" in t_res:
                it_number = get_it_number_by_name_hm_full(self.name_hm_alt)
            elif "short" in t_res:
                it_number = get_it_number_by_name_hm_short(self.name_hm_alt)
            elif "extended" in t_res:
                it_number, it_coordinate_system_codes = \
                    get_it_number_it_coordinate_system_codes_by_name_hm_extended(
                        self.name_hm_alt)
                if (not(it_coordinate_system_code in
                        it_coordinate_system_codes)):
                    it_coordinate_system_code = \
                        auto_choose_it_coordinate_system_code(
                            it_number, it_coordinate_system_codes)
        elif self.is_attribute("name_hm_alt"):
            t_res = get_type_hm(self.name_hm_alt)
            if "extended" in t_res:
                it_number, it_coordinate_system_codes = \
                    get_it_number_it_coordinate_system_codes_by_name_hm_extended(
                        self.name_hm_alt)
                if (not(it_coordinate_system_code in
                        it_coordinate_system_codes)):
                    it_coordinate_system_code = \
                        auto_choose_it_coordinate_system_code(
                            it_number, it_coordinate_system_codes)
            elif "full" in t_res:
                it_number = get_it_number_by_name_hm_full(self.name_hm_alt)
            elif "short" in t_res:
                it_number = get_it_number_by_name_hm_short(self.name_hm_alt)
        elif self.is_attribute("name_hm_ref"):
            it_number = get_it_number_by_name_hm_short(self.name_hm_ref)
        elif self.is_attribute("name_hall"):
            it_number = get_it_number_by_name_hall(self.name_hall)
        elif self.is_attribute("it_number"):
            it_number = self.it_number
            # it_coordinate_system_code = self.it_coordinate_system_code
        elif self.is_attribute("name_schoenflies"):
            it_number = get_it_number_by_name_schoenflies(
                self.name_schoenflies)
        else:
            flag = False
        if flag:
            it_coordinate_system_codes = \
                get_it_coordinate_system_codes_by_it_number(it_number)
            self.__dict__["it_number"] = it_number
            if (not(it_coordinate_system_code in it_coordinate_system_codes)):
                it_coordinate_system_code = \
                    auto_choose_it_coordinate_system_code(
                        it_number, it_coordinate_system_codes)
            self.__dict__["it_coordinate_system_code"] = \
                it_coordinate_system_code
        return flag

    def form_object_by_it_number_it_coordinate_system_code(self):
        """Form object by it number it coordinate system code.

        TODO: Solve the problem with centring_type for hexagonal systems
              (Ex.: 166 space group).
        """
        bravais_type, laue_class, patterson_name_hm, centring_type, \
            crystal_system = None, None, None, None, None
        name_hm_extended, name_hm_full, name_hm_short = None, None, None
        name_hall, name_schoenflies, point_group_hm = None, None, None
        pcentr = None
        reduced_space_group_symop, full_space_group_symop = None, None
        space_group_wyckoff = None

        lattice_type = None
        generators = ()
        self.define_it_number_and_it_coordinate_system_code()
        try:
            it_number = self.it_number
        except AttributeError:
            it_number = None
        try:
            it_c_s_c = self.it_coordinate_system_code
        except AttributeError:
            it_c_s_c = None
        crystal_system = get_crystal_system_by_it_number(it_number)
        name_hm_extended = \
            get_name_hm_extended_by_it_number_it_coordinate_system_code(
                it_number, it_c_s_c)
        if (name_hm_extended is not None):
            centring_type = get_centring_type_by_name_hm_extended(
                name_hm_extended)
        name_hm_short = get_name_hm_short_by_it_number(it_number)
        if (name_hm_short is not None):
            lattice_type = get_lattice_type_by_name_hm_short(name_hm_short)
            if centring_type is None:
                centring_type = get_centring_type_by_name_hm_extended(
                    name_hm_short)
        name_hm_full = get_name_hm_full_by_it_number(it_number)
        if ((name_hm_full is not None) & (centring_type is None)):
            centring_type = get_centring_type_by_name_hm_extended(name_hm_full)
        if ((centring_type is not None) & (crystal_system is not None)):
            bravais_type = get_bravais_type_by_centring_type_crystal_system(
                centring_type, crystal_system)
        name_hall = get_name_hall_by_it_number(it_number)
        if name_hall is not None:
            centrosymmetry = get_centrosymmetry_by_name_hall(name_hall)
        name_schoenflies = get_name_schoenflies_by_it_number(it_number)
        if name_schoenflies is not None:
            laue_class = get_laue_class_by_name_schoenflies(name_schoenflies)
            point_group_hm = get_point_group_hm_short_by_name_schoenflies(
                name_schoenflies)
            if point_group_hm is not None:
                generators = get_generators_by_point_group_hm(point_group_hm)
            if ((lattice_type is not None) & (laue_class is not None)):
                patterson_name_hm = \
                    get_patterson_name_hm_by_lattice_type_laue_class(
                        lattice_type, laue_class)
        if ((centring_type is not None) and (it_c_s_c is not None)):
            if centring_type.startswith("R") and it_c_s_c.startswith("r"):
                centring_type = "P"

        pcentr, reduced_space_group_symop, full_space_group_symop, \
            space_group_wyckoff = \
            self.form_reduced_full_space_group_symop_space_group_wyckoff(
                it_number, it_c_s_c, centrosymmetry, centring_type)

        if pcentr is not None:
            self.__dict__["pcentr"] = pcentr
        if reduced_space_group_symop is not None:
            self.__dict__["reduced_space_group_symop"] = \
                reduced_space_group_symop
        if full_space_group_symop is not None:
            self.__dict__["full_space_group_symop"] = full_space_group_symop
        if space_group_wyckoff is not None:
            self.__dict__["space_group_wyckoff"] = space_group_wyckoff

        if centrosymmetry is not None:
            self.__dict__["centrosymmetry"] = centrosymmetry
        # if name_hm_extended is not None:
        #     self.__dict__["name_hm_alt"] = name_hm_extended
        #     self.__dict__["name_hm_alt_description"] = \
        #         "The extended Hermann-Mauguin space-group symbol."
        if name_hm_full is not None:
            self.__dict__["name_hm_full"] = name_hm_full
        if name_hm_short is not None:
            self.__dict__["name_hm_ref"] = name_hm_short
            self.__dict__["name_hm_alt"] = name_hm_short
            self.__dict__["name_hm_alt_description"] = \
                "The short Hermann-Mauguin space-group symbol."
        if name_hall is not None:
            self.__dict__["name_hall"] = name_hall
        if name_schoenflies is not None:
            self.__dict__["name_schoenflies"] = name_schoenflies

        if point_group_hm is not None:
            self.__dict__["point_group_hm"] = point_group_hm
        if laue_class is not None:
            self.__dict__["laue_class"] = laue_class
        if patterson_name_hm is not None:
            self.__dict__["patterson_name_hm"] = patterson_name_hm
        if centring_type is not None:
            self.__dict__["centring_type"] = centring_type
        if bravais_type is not None:
            self.__dict__["bravais_type"] = bravais_type
        if crystal_system is not None:
            self.__dict__["crystal_system"] = crystal_system

    def form_reduced_full_space_group_symop_space_group_wyckoff(
            self, it_number: int, it_coordinate_system_code: str,
            centrosymmetry: bool, centring_type: str):
        """Form reduced full space group symop space group Wyckoff."""
        pcentr, reduced_space_group_symop, \
            full_space_group_symop = None, None, None

        symop, pcentr, _multiplicity, _letter, _site_symmetry, _l_coords_xyz_2\
            = get_symop_pcentr_multiplicity_letter_site_symmetry_coords_xyz_2(
                it_number, it_coordinate_system_code)

        item_reduced = []
        for _i_symop, _symop in enumerate(symop):
            _item = SpaceGroupSymop(id=f"{_i_symop+1:}", operation_xyz=_symop,
                                    operation_description="reduced",
                                    sg_id=it_number)
            item_reduced.append(_item)
        reduced_space_group_symop = SpaceGroupSymopL()
        reduced_space_group_symop.items = item_reduced
        
        item_full = []
        for _item in item_reduced:
            _symop = _item.operation_xyz
            _item_full = SpaceGroupSymop(operation_xyz=_symop)
            item_full.append(_item_full)
        if centrosymmetry:
            item_add = []
            for _item in item_full:
                _item_add = _item.get_symop_inversed(pcentr)
                item_add.append(_item_add)
            item_full.extend(item_add)
        item_new = []
        for _item in item_full:
            _items_new = _item.get_symops_by_centring_type(centring_type)
            item_new.extend(_items_new)
        for _i, _item in enumerate(item_new):
            _item.id = f"{_i+1:}"
            _item.sg_id = f"{it_number:}"
            _item.operation_description = "full"
        full_space_group_symop = SpaceGroupSymopL()
        full_space_group_symop.items = item_new

        _i_numb = 0
        item = []
        for _1, _2, _3, _4 in zip(_multiplicity, _letter, _site_symmetry,
                                  _l_coords_xyz_2):
            _i_numb += 1
            _item = SpaceGroupWyckoff(id=f"{_i_numb:}", multiplicity=_1,
                                      coord_xyz=_4[0], letter=_2,
                                      site_symmetry=_3)
            _item.it_coord_xyz = _4
            _item.centring_type = centring_type
            _item.form_object()
            item.append(_item)
        space_group_wyckoff = SpaceGroupWyckoffL()
        space_group_wyckoff.items = item
        return pcentr, reduced_space_group_symop, full_space_group_symop, \
            space_group_wyckoff

    @staticmethod
    def get_it_coordinate_system_codes_by_it_number(it_number: int) -> str:
        """Get it coordinate system codes by it number."""
        return get_it_coordinate_system_codes_by_it_number(it_number)

    @staticmethod
    def get_default_it_coordinate_system_code_by_it_number(it_number: int) \
            -> str:
        """Get default it coordinate system code by it number."""
        return get_default_it_coordinate_system_code_by_it_number(it_number)

    def calc_hkl_equiv(self, index_h, index_k, index_l):
        """Gequivalent reflections of hkl and its multiplicity."""
        r_11 = numpy.array(self.reduced_space_group_symop.r_11, dtype=Fraction)
        r_12 = numpy.array(self.reduced_space_group_symop.r_12, dtype=Fraction)
        r_13 = numpy.array(self.reduced_space_group_symop.r_13, dtype=Fraction)
        r_21 = numpy.array(self.reduced_space_group_symop.r_21, dtype=Fraction)
        r_22 = numpy.array(self.reduced_space_group_symop.r_22, dtype=Fraction)
        r_23 = numpy.array(self.reduced_space_group_symop.r_23, dtype=Fraction)
        r_31 = numpy.array(self.reduced_space_group_symop.r_31, dtype=Fraction)
        r_32 = numpy.array(self.reduced_space_group_symop.r_32, dtype=Fraction)
        r_33 = numpy.array(self.reduced_space_group_symop.r_33, dtype=Fraction)

        h_s = r_11*index_h + r_21*index_k + r_31*index_l
        k_s = r_12*index_h + r_22*index_k + r_32*index_l
        l_s = r_13*index_h + r_23*index_k + r_33*index_l

        hkl_s = numpy.vstack([h_s, k_s, l_s])
        hkl_s = numpy.hstack([hkl_s, -1*hkl_s])
        hkl_s = hkl_s.astype(int)
        hkl_s_un = numpy.unique(hkl_s, axis=1)
        multiplicity = int(round(hkl_s.shape[1]*1./hkl_s_un.shape[1]))
        h_s, k_s, l_s = hkl_s_un[0, :], hkl_s_un[1, :], hkl_s_un[2, :]
        return h_s, k_s, l_s, multiplicity

    def calc_xyz_mult(self, x, y, z):
        """Give unique x,y,z elements and calculate multiplicity.

        It's done for given x,y,z fract.
        """
        wyckoff = self.space_group_wyckoff.get_wyckoff_for_fract(x, y, z)
        np_r = numpy.array(wyckoff.full_r, dtype=float)
        np_b = numpy.array(wyckoff.full_b, dtype=float)
        x_s = (np_r[:, 0, 0]*x + np_r[:, 0, 1]*y + np_r[:, 0, 2]*z +
               np_b[:, 0]) % 1
        y_s = (np_r[:, 1, 0]*x + np_r[:, 1, 1]*y + np_r[:, 1, 2]*z +
               np_b[:, 1]) % 1
        z_s = (np_r[:, 2, 0]*x + np_r[:, 2, 1]*y + np_r[:, 2, 2]*z +
               np_b[:, 2]) % 1

        # l_shift = self.shift
        # l_x, l_y, l_z = [], [], []
        # for _shift in l_shift:
        #    l_x.extend(numpy.mod(x_s+_shift[0],1))
        #    l_y.extend(numpy.mod(y_s+_shift[1],1))
        #    l_z.extend(numpy.mod(z_s+_shift[2],1))
        multiplicity = wyckoff.multiplicity
        # x_out = numpy.array(l_x, dtype=float)
        # y_out = numpy.array(l_y, dtype=float)
        # z_out = numpy.array(l_z, dtype=float)

        return x_s, y_s, z_s, multiplicity

    def calc_symop_for_xyz(self, x_in, y_in, z_in):
        """Calculate symop for xyz."""
        x, y, z = x_in % 1., y_in % 1., z_in % 1.

        symop = self.full_space_group_symop
        e_11 = numpy.array(symop.r_11, dtype=float)
        e_12 = numpy.array(symop.r_12, dtype=float)
        e_13 = numpy.array(symop.r_13, dtype=float)
        e_1 = numpy.array(symop.b_1, dtype=float)

        e_21 = numpy.array(symop.r_21, dtype=float)
        e_22 = numpy.array(symop.r_22, dtype=float)
        e_23 = numpy.array(symop.r_23, dtype=float)
        e_2 = numpy.array(symop.b_2, dtype=float)

        e_31 = numpy.array(symop.r_31, dtype=float)
        e_32 = numpy.array(symop.r_32, dtype=float)
        e_33 = numpy.array(symop.r_33, dtype=float)
        e_3 = numpy.array(symop.b_3, dtype=float)

        x_s = numpy.round(numpy.mod(e_11*x + e_12*y + e_13*z + e_1, 1), decimals=5)
        y_s = numpy.round(numpy.mod(e_21*x + e_22*y + e_23*z + e_2, 1), decimals=5)
        z_s = numpy.round(numpy.mod(e_31*x + e_32*y + e_33*z + e_3, 1), decimals=5)

        xyz_s = numpy.vstack([x_s, y_s, z_s])

        xyz_s_un, unique_inverse = numpy.unique(xyz_s, return_inverse=True,
                                                axis=1)
        x_s, y_s, z_s = xyz_s_un[0, :], xyz_s_un[1, :], xyz_s_un[2, :]
        ind = (numpy.where((x-x_s)**2+(y-y_s)**2+(z-z_s)**2 < 0.00001))[0][0]

        flag = unique_inverse == ind

        item_out = [_item for _item, _flag in zip(symop.items, flag) if _flag]
        symop_out = SpaceGroupSymopL()
        symop_out.items = item_out
        return symop_out

    def calc_el_symm_for_xyz(self, x_in, y_in, z_in):
        """FIXME: should be deleted."""
        x, y, z = x_in % 1., y_in % 1., z_in % 1.

        symop = self.full_space_group_symop
        e_11 = numpy.array(symop.r_11, dtype=float)
        e_12 = numpy.array(symop.r_12, dtype=float)
        e_13 = numpy.array(symop.r_13, dtype=float)
        e_1 = numpy.array(symop.b_1, dtype=float)

        e_21 = numpy.array(symop.r_21, dtype=float)
        e_22 = numpy.array(symop.r_22, dtype=float)
        e_23 = numpy.array(symop.r_23, dtype=float)
        e_2 = numpy.array(symop.b_2, dtype=float)

        e_31 = numpy.array(symop.r_31, dtype=float)
        e_32 = numpy.array(symop.r_32, dtype=float)
        e_33 = numpy.array(symop.r_33, dtype=float)
        e_3 = numpy.array(symop.b_3, dtype=float)

        x_s = numpy.round(numpy.mod(e_11*x + e_12*y + e_13*z + e_1, 1), decimals=5)
        y_s = numpy.round(numpy.mod(e_21*x + e_22*y + e_23*z + e_2, 1), decimals=5)
        z_s = numpy.round(numpy.mod(e_31*x + e_32*y + e_33*z + e_3, 1), decimals=5)

        xyz_s = numpy.vstack([x_s, y_s, z_s])

        xyz_s_un, unique_inverse = numpy.unique(xyz_s, return_inverse=True,
                                                axis=1)
        # Output for unique_inverse in version of numpy 2.0.0 and 1.23.5 is different. 
        # This line is to unified two versions
        unique_inverse = unique_inverse.flatten()
        x_s, y_s, z_s = xyz_s_un[0, :], xyz_s_un[1, :], xyz_s_un[2, :]
        ind = (numpy.where((x-x_s)**2+(y-y_s)**2+(z-z_s)**2 < 0.00001))[0][0]
        flag = unique_inverse == ind
        o_11, o_12, o_13 = e_11[flag], e_12[flag], e_13[flag]
        o_21, o_22, o_23 = e_21[flag], e_22[flag], e_23[flag]
        o_31, o_32, o_33 = e_31[flag], e_32[flag], e_33[flag]
        o_1, o_2, o_3 = e_1[flag], e_2[flag], e_3[flag]
        return o_11, o_12, o_13, o_21, o_22, o_23, o_31, o_32, o_33, o_1, \
            o_2, o_3

    def calc_f_hkl_by_f_hkl_as(self, index_h, index_k, index_l, f_hkl_as):
        """Calculate the structure factor.

        The structure factor in asymmetric unit cell are given.
        """
        np_h = numpy.array(index_h, dtype=float)
        np_k = numpy.array(index_k, dtype=float)
        np_l = numpy.array(index_l, dtype=float)
        shift = numpy.array(self.shift, dtype=float)
        centr = self.centrosymmetry

        orig_x, orig_y, orig_z = shift[:, 0], shift[:, 1], shift[:, 2]

        # orig_x = [hh[0] for hh in shift]
        # orig_y = [hh[1] for hh in shift]
        # orig_z = [hh[2] for hh in shift]

        np_h_2d, np_orig_x_2d = numpy.meshgrid(index_h, orig_x, indexing="ij")
        np_k_2d, np_orig_y_2d = numpy.meshgrid(index_k, orig_y, indexing="ij")
        np_l_2d, np_orig_z_2d = numpy.meshgrid(index_l, orig_z, indexing="ij")

        hh = (2*numpy.pi*1j*(np_h_2d*np_orig_x_2d + np_k_2d*np_orig_y_2d +
                             np_l_2d*np_orig_z_2d)).astype(complex)
        np_orig_as = numpy.exp(hh)
        _hh = np_orig_as.sum(axis=1)
        if len(f_hkl_as.shape) == 2:
            _hh = _hh[:, numpy.newaxis]
        f_hkl_1 = f_hkl_as*_hh*1./len(shift)

        if (centr):
            p_centr = numpy.array(self.pcentr, dtype=float)
            hh = (2. * 2. * numpy.pi * 1j * (
                np_h*p_centr[0] + np_k*p_centr[1] + np_l*p_centr[2])
                ).astype(complex)
            if len(f_hkl_as.shape) == 2:
                hh = hh[:, numpy.newaxis]
            f_hkl = 0.5*(f_hkl_1+f_hkl_1.conjugate()*numpy.exp(hh))
        else:
            f_hkl = f_hkl_1
        return f_hkl

    def calc_asymmetric_cell(self, n_a, n_b, n_c):
        """Give the numbers in asymmetric cell.

         :n_a: the number of points along a axis
         :n_b: the numper of points along b axis
         :n_c: the numper of points along c axis

        na, n_b, nc should be divided on 24: 8 and 3

        :output: - l_coord is a list of coordinates in asymmetric cell
                   (frac_x = n_x/n_a and so on)
                 - l_symm contains a list of symmetry given as
                   (n_symm, centr, n_orig)
        """
        n_a_new = int(round(n_a/24))*24
        # n_b_new = int(round(n_b/24))*24
        # n_c_new = int(round(n_c/24))*24

        # l_el_symm = self.el_symm
        # f_centr = self.centr
        # p_centr = self.p_centr
        # l_orig = self.orig
        l_coord = []

        spgr_choice = self.spgr_choice
        # spgr_name = self.spgr_name
        spgr_number = self.spgr_number

        if (spgr_number == 227) & (spgr_choice == "2"):
            n_a_new = int(round(n_a/8))*8
            for n_x in range(-n_a_new//8, 3*n_a_new//8+1):
                for n_y in range(-n_a_new//8, 0+1):
                    for n_z in range(-n_a_new//4, 0+1):
                        cond_1 = (n_y < min([n_a_new//4-n_x, n_x]))
                        cond_2 = (n_z >= -n_y-n_a_new//4)
                        cond_3 = (n_z <= n_y)
                        if (cond_1 & cond_2 & cond_3):
                            coord_x = float(n_x)/float(n_a_new)
                            coord_y = float(n_y)/float(n_a_new)
                            coord_z = float(n_z)/float(n_a_new)
                            l_coord.append((coord_x, coord_y, coord_z))
        return l_coord

    def calc_rotated_matrix_for_position(self, m_chi, x, y, z):
        """Calc rotated matrix for position."""
        symop = self.full_space_group_symop
        e_11 = numpy.array(symop.r_11, dtype=float)
        e_12 = numpy.array(symop.r_12, dtype=float)
        e_13 = numpy.array(symop.r_13, dtype=float)
        e_1 = numpy.array(symop.b_1, dtype=float)

        e_21 = numpy.array(symop.r_21, dtype=float)
        e_22 = numpy.array(symop.r_22, dtype=float)
        e_23 = numpy.array(symop.r_23, dtype=float)
        e_2 = numpy.array(symop.b_2, dtype=float)

        e_31 = numpy.array(symop.r_31, dtype=float)
        e_32 = numpy.array(symop.r_32, dtype=float)
        e_33 = numpy.array(symop.r_33, dtype=float)
        e_3 = numpy.array(symop.b_3, dtype=float)

        x_s = numpy.round(numpy.mod(e_11*x + e_12*y + e_13*z + e_1, 1), decimals=5)
        y_s = numpy.round(numpy.mod(e_21*x + e_22*y + e_23*z + e_2, 1), decimals=5)
        z_s = numpy.round(numpy.mod(e_31*x + e_32*y + e_33*z + e_3, 1), decimals=5)
        # o_11, o_12, o_13, o_21, o_22, o_23, o_31, o_32, o_33, o_1, o_2, o_3 \
        #    = self.calc_el_symm_for_xyz(x, y, z)
        # np_x, np_y, np_z, mult = self.calc_xyz_mult(x, y, z)

        l_ind, l_xyz = [], []
        _ind = 0
        for _x, _y, _z in zip(x_s, y_s, z_s):
            if (_x, _y, _z) not in l_xyz:
                l_ind.append(_ind)
                l_xyz.append((_x, _y, _z))
            _ind += 1
        l_res = []
        for _ind, _xyz in zip(l_ind, l_xyz):
            _11, _12, _13 = e_11[_ind], e_12[_ind], e_13[_ind]
            _21, _22, _23 = e_21[_ind], e_22[_ind], e_23[_ind]
            _31, _32, _33 = e_31[_ind], e_32[_ind], e_33[_ind]
            # _1, _2, _3 = e_1[_ind], e_2[_ind], e_3[_ind]
            matrix_r = numpy.array([[_11, _12, _13], [_21, _22, _23],
                                    [_31, _32, _33]], dtype=float)
            matrix_rt = matrix_r.transpose()
            r_chi = numpy.matmul(matrix_r, m_chi)
            matrix_chi_rot = numpy.matmul(r_chi, matrix_rt)
            l_res.append((_xyz, matrix_chi_rot))
        return l_res

    def report(self) -> str:
        return self.report_space_group()

    def report_space_group(self) -> str:
        """Make a report about space group in string format."""
        if self.is_defined:
            self.form_object
        else:
            return ""
        ls_out = []
        ls_out.append("# Space group")

        width_left, width_right = 30, 40

        if self.is_attribute("it_number"):
            ls_out.append("|IT_number: ".rjust(width_left) + "|" +
                          f"{self.it_number:}|".ljust(width_right))
        if self.is_attribute("name_hm_alt"):
            ls_out.append(
                "|Name H-M alt: ".rjust(width_left) + "|" +
                f"\"{self.name_hm_alt:}\"|".ljust(width_right))
        if self.is_attribute("name_hm_full"):
            ls_out.append(
                "|Name H-M full: ".rjust(width_left) + "|" +
                f"\"{self.name_hm_full:}\"|".ljust(width_right))
        if self.is_attribute("name_hm_ref"):
            ls_out.append(
                "|Name H-M ref: ".rjust(width_left) + "|" +
                f"\"{self.name_hm_ref:}\"|".ljust(width_right))
        if self.is_attribute("name_hall"):
            ls_out.append("|Name Hall short: ".rjust(width_left) + "|" +
                          f"\"{self.name_hall:}\"|".ljust(width_right))
        if self.is_attribute("name_schoenflies"):
            ls_out.append(
                "|Name Schoenflies: ".rjust(width_left) + "|" +
                f"\"{self.name_schoenflies:}\"|".ljust(width_right))
        if self.is_attribute("it_coordinate_system_code"):
            ls_out.append("|IT_coordinate_system_code: ".rjust(width_left) + "|" +
                          f"\"{self.it_coordinate_system_code:}\"|".ljust(
                              width_right))
        ls_out.append("")
        if self.is_attribute("point_group_hm"):
            ls_out.append(
                "|Point group H-M: ".rjust(width_left) + "|" +
                f"\"{self.point_group_hm:}\"|".ljust(width_right))
        if self.is_attribute("laue_class"):
            ls_out.append("|Laue class: ".rjust(width_left) + "|" +
                          f"\"{self.laue_class:}\"|".ljust(width_right))
        if self.is_attribute("patterson_name_hm"):
            ls_out.append(
                "|Patterson name H-M: ".rjust(width_left) + "|" +
                f"\"{self.patterson_name_hm:}\"|".ljust(width_right))
        if self.is_attribute("centring_type"):
            ls_out.append(
                "|Centring type: ".rjust(width_left) + "|" +
                f"\"{self.centring_type:}\"|".ljust(width_right))
        if self.is_attribute("bravais_type"):
            ls_out.append(
                "|Bravais type: ".rjust(width_left) + "|" +
                f"\"{self.bravais_type:}\"|".ljust(width_right))
        if self.is_attribute("crystal_system"):
            ls_out.append(
                "|Crystal system: ".rjust(width_left) + "|" +
                f"\"{self.crystal_system:}\"|".ljust(width_right))
        ls_out.append("")
        if self.is_attribute("centrosymmetry"):
            ls_out.append(
                "|Centrosymmetry: ".rjust(width_left) + "|" +
                f"{'Yes' if self.centrosymmetry else 'No':}|".ljust(
                    width_right))

        # if generators != (): print(
        #    f"Generators: ".rjust(width_left) +
        #    ", ".join([f"\"{_}\"" for _ in generators]).ljust(width_right))

        # if symop is not None:
        #    print("Symop: ")  # pcentr
        #    print_long_list([f"\"{_:}\"" for _ in symop])

        if self.is_attribute("space_group_wyckoff"):
            s_g_w = self.space_group_wyckoff
            ls_out.append("\n## Special positions:\n")
            ls_out.append("|Position | Multiplicity| Character| Group|")
            ls_out.append("|---------|-------------|----------|------|")
            for _1, _2, _3, _4 in zip(s_g_w.multiplicity, s_g_w.letter,
                                      s_g_w.site_symmetry, s_g_w.coord_xyz):
                ls_out.append(
                    f"|{_4.rjust(12):12}| {_1:3}| {_2.rjust(3):3}| \
{_3.rjust(7):7}|")

        if self.is_attribute("reduced_space_group_symop"):
            r_s_g_s = self.reduced_space_group_symop
            ls_out.append("\n## Reduced space group symop:\n")
            line = []
            for _i,  _1 in enumerate(r_s_g_s.operation_xyz):
                line.append(f" {_1.rjust(12):12}")
                if _i % 3 == 2:
                    ls_out.append("|"+"|".join(line)+"|")
                    line = []

        if self.is_attribute("shift"):
            shift = self.shift
            ls_out.append("\n## Shift:\n")
            ls_out.append("|along a | along b| along c|")
            ls_out.append("|--------|--------|--------|")
            for _pos in shift:
                ls_out.append(f"|{float(_pos[0]):9.5f}| {float(_pos[1]):9.5f}| \
{float(_pos[2]):9.5f}|")
        return "\n".join(ls_out)

    def to_cif(self, separator: str = "_", flag_all_attributes: bool = False,
               flag_minimal: bool = True) -> str:
        """Print information about object in string in STAR format.

        Arguments
        ---------
            prefix is a prefix in front of label of attribute
            separator is a separator between prefix and attribute ("_" or ".")
            flag if it's True the value "." will be printed for undefined
            attributes flag_minimal if it's True the minimal set of object will
            be printed

        Returns
        -------
            A string in STAR/CIF format
        """
        if flag_minimal:
            ls_out = []
            prefix = self.PREFIX
            attributes = ["name_hm_alt", "it_coordinate_system_code"]
            related_attributes = ["name_H-M_alt", "IT_coordinate_system_code"]
            for _attr, cif_attr in zip(attributes, related_attributes):
                if self.is_attribute(_attr):
                    _val = getattr(self, _attr)
                    s_val = str(_val)
                    if len(s_val.split(" ")) > 1:
                        ls_out.append(
                            f"_{prefix:}{separator:}{cif_attr:} \"{s_val:}\"")
                    else:
                        ls_out.append(
                            f"_{prefix:}{separator:}{cif_attr:} {s_val:}")
                else:
                    ls_out.append(f"_{prefix:}{separator:}{cif_attr:} .")
            s_out = "\n".join(ls_out)
        else:
            s_out = super(SpaceGroup, self).to_cif(
                separator=separator, flag_all_attributes=flag_all_attributes)
        return s_out


class SpaceGroupL(LoopN):
    """SpaceGroupL class.

    Contains information about Wyckoff positions of a space group.
    Only one site can be given for each special position but the
    remainder can be generated by applying the symmetry operations
    stored in _space_group_symop.operation_xyz.
    """

    ITEM_CLASS = SpaceGroup
    ATTR_INDEX = "id"

    def __init__(self, loop_name: str = None, **kwargs) -> NoReturn:
        super(SpaceGroupL, self).__init__()
        self.__dict__["items"] = form_items_by_dictionary(self.ITEM_CLASS, kwargs)
        self.__dict__["loop_name"] = loop_name

