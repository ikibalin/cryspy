"""Settings of the space groups (ITA Table 4.3.2.1) and their symbols.

The expected operators are taken from International Tables A (2016) and
checked against gemmi. The last test compares all 530 settings with gemmi
when it is installed (a development check, gemmi is not a dependency).
"""
from fractions import Fraction

import pytest

import cryspy
from cryspy.A_functions_base.function_1_strings import transform_string_to_r_b
from cryspy.A_functions_base.function_2_space_group import (
    REFERENCE_TABLE_IT_NUMBER_IT_COORDINATE_SYSTEM_CODE_NAME_HALL,
    get_it_number_it_coordinate_system_code_by_name_hall,
    get_it_number_it_coordinate_system_code_by_name_hm)


def normalize_operators(operators):
    """Set of operators (rotation, translation modulo 1)."""
    res = set()
    for operator in operators:
        r, b = transform_string_to_r_b(operator, ("x", "y", "z"))
        res.add((tuple(Fraction(_) for _ in r.flatten()),
                 tuple(Fraction(_) % 1 for _ in b)))
    return frozenset(res)


def get_operators(space_group):
    return normalize_operators(
        [item.operation_xyz for item in space_group.full_space_group_symop.items])


OPERATORS_BMAB = (
    "x,y,z", "-x+1/2,-y+1/2,z", "-x+1/2,-y,-z+1/2", "-x+1/2,y+1/2,-z",
    "-x+1/2,y,z+1/2", "-x,-y+1/2,z+1/2", "-x,-y,-z", "-x,y+1/2,-z+1/2",
    "-x,y,z", "x+1/2,-y+1/2,z", "x+1/2,-y,-z+1/2", "x+1/2,y+1/2,-z",
    "x+1/2,y,z+1/2", "x,-y+1/2,z+1/2", "x,-y,-z", "x,y+1/2,-z+1/2")


@pytest.mark.parametrize("name", ["B m a b", "Bmab", "B m e b", "Bmeb"])
def test_cmce_setting_a_cb(name):
    space_group = cryspy.SpaceGroup(name_hm_alt=name)
    assert space_group.it_number == 64
    assert space_group.it_coordinate_system_code == "a-cb"
    assert space_group.centring_type == "B"
    assert get_operators(space_group) == normalize_operators(OPERATORS_BMAB)


def test_setting_is_kept_in_cif():
    space_group = cryspy.SpaceGroup(name_hm_alt="B m a b")
    cif = space_group.to_cif()
    assert '"B m e b"' in cif
    space_group_2 = cryspy.SpaceGroup.from_cif(cif)
    assert get_operators(space_group_2) == normalize_operators(OPERATORS_BMAB)


@pytest.mark.parametrize("name, it_number, code", [
    ("C m c a", 64, "abc"), ("Cmca", 64, "abc"), ("C m c e", 64, "abc"),
    ("A b m a", 64, "cab"), ("A b m 2", 39, "abc"), ("A b a 2", 41, "abc"),
    ("C m m a", 67, "abc"), ("C c c a", 68, "2abc"), ("C c c a:1", 68, "1abc"),
    ("B b a b:1", 68, "1a-cb"),
    ("P 21/n", 14, "b2"), ("P2_1/n", 14, "b2"), ("P 1 1 21/n", 14, "c2"),
    ("I 1 2/a 1", 15, "b3"), ("C 2/n 1 1", 15, "a2"),
    ("F d -3 m", 227, "2"), ("F d -3 m:1", 227, "1"), ("Fd-3m:2", 227, "2"),
    ("R -3 m", 166, "h"), ("R -3 m:R", 166, "r"), ("R -3 m:H", 166, "h"),
    ("P 4/n m m:1", 129, "1"), ("P23", 195, None), ("P 1", 1, None),
])
def test_name_hm(name, it_number, code):
    assert get_it_number_it_coordinate_system_code_by_name_hm(name) == \
        (it_number, code)
    space_group = cryspy.SpaceGroup(name_hm_alt=name)
    assert space_group.it_number == it_number
    if code is not None:
        assert space_group.it_coordinate_system_code == code


@pytest.mark.parametrize("name, it_number, code", [
    ("-B 2ab 2", 64, "a-cb"), ("P 3*", 146, "r"), ("-P 2yn", 14, "b2"),
    ("-C 2a 2ac", 68, "2abc"), ("P 4ab 2ab -1ab", 129, "1"),
])
def test_name_hall(name, it_number, code):
    assert get_it_number_it_coordinate_system_code_by_name_hall(name) == \
        (it_number, code)
    space_group = cryspy.SpaceGroup(name_hall=name)
    assert space_group.it_number == it_number
    assert space_group.it_coordinate_system_code == code


def test_hall_table_has_all_settings():
    table = REFERENCE_TABLE_IT_NUMBER_IT_COORDINATE_SYSTEM_CODE_NAME_HALL
    assert len(table) == 530
    assert {_[0] for _ in table} == set(range(1, 231))


@pytest.mark.parametrize("code, n_operators", [
    ("b1", 8), ("b2", 8), ("b3", 8), ("-b1", 8), ("-b2", 8), ("-b3", 8),
    ("c1", 8), ("c2", 8), ("c3", 8), ("-c1", 8), ("-c2", 8), ("-c3", 8),
    ("a1", 8), ("a2", 8), ("a3", 8), ("-a1", 8), ("-a2", 8), ("-a3", 8)])
def test_monoclinic_cell_choices(code, n_operators):
    """Cell choices 2 and 3 crashed before (no card for them)."""
    space_group = cryspy.SpaceGroup(it_number=15, it_coordinate_system_code=code)
    assert len(space_group.full_space_group_symop.items) == n_operators


def test_monoclinic_cell_choice_2_operators():
    space_group = cryspy.SpaceGroup(it_number=15, it_coordinate_system_code="b2")
    assert space_group.centring_type == "A"
    assert get_operators(space_group) == normalize_operators((
        "x,y,z", "-x+1/2,y,-z+1/2", "-x,-y,-z", "x+1/2,-y,z+1/2",
        "x,y+1/2,z+1/2", "-x+1/2,y+1/2,-z", "-x,-y+1/2,-z+1/2",
        "x+1/2,-y+1/2,z"))


def test_p31_operators():
    """P 31 had the operators of P 32."""
    space_group = cryspy.SpaceGroup(it_number=144)
    assert get_operators(space_group) == normalize_operators(
        ("x,y,z", "-y,x-y,z+1/3", "-x+y,-x,z+2/3"))


def test_origin_choice_1_with_axes_permutation():
    """The inversion centre was not transformed with the axes (P n c b:1)."""
    space_group = cryspy.SpaceGroup(
        it_number=50, it_coordinate_system_code="1cab")
    assert get_operators(space_group) == normalize_operators((
        "x,y,z", "-x,-y,z", "-x,y,-z", "x,-y,-z", "-x,-y+1/2,-z+1/2",
        "x,y+1/2,-z+1/2", "x,-y+1/2,z+1/2", "-x,y+1/2,z+1/2"))


def test_pnnn_origin_choice_1():
    """Inversion centre was at (1/4, 1/4, 1/5) instead of (1/4, 1/4, 1/4)."""
    space_group = cryspy.SpaceGroup(
        it_number=48, it_coordinate_system_code="1abc")
    assert normalize_operators(("-x+1/2,-y+1/2,-z+1/2",)) <= \
        get_operators(space_group)


def test_fd3c_origin_choice_2():
    space_group = cryspy.SpaceGroup(it_number=228, it_coordinate_system_code="2")
    assert len(space_group.full_space_group_symop.items) == 192


def test_wyckoff_with_twelfths():
    """7/12 was rounded to 4/7 (denominators were limited to 10)."""
    space_group = cryspy.SpaceGroup(it_number=178)
    wyckoff_b = [item for item in space_group.space_group_wyckoff.items
                 if item.letter == "b"][0]
    assert wyckoff_b.it_coord_xyz == [
        "x,2x,1/4", "-2x,-x,7/12", "x,-x,11/12", "-x,-2x,3/4", "2x,x,1/12",
        "-x,x,5/12"]


def test_all_settings_against_gemmi():
    """Operators of all 530 settings and Wyckoff orbits (development check)."""
    gemmi = pytest.importorskip("gemmi")

    def get_code(group):
        number, qualifier, ext = group.number, group.qualifier, group.ext
        if 3 <= number <= 15:
            return qualifier + "1" if qualifier in ("a", "b", "c") else qualifier
        if 16 <= number <= 74:
            qualifier = qualifier or "abc"
            return ext + qualifier if ext in ("1", "2") else qualifier
        return {"1": "1", "2": "2", "H": "h", "R": "r"}.get(ext)

    for group in gemmi.spacegroup_table_itb():
        code = get_code(group)
        space_group = cryspy.SpaceGroup(
            it_number=group.number, it_coordinate_system_code=code)
        operators_gemmi = normalize_operators(
            [op.triplet() for op in group.operations()])
        assert get_operators(space_group) == operators_gemmi, group.xhm()
        hall = get_it_number_it_coordinate_system_code_by_name_hall(group.hall)
        assert hall[0] == group.number, group.xhm()


CIF_OLD_TAGS_BMAB = """data_La2CuO4
_space_group_IT_number           64
_symmetry_Int_Tables_number      64
_symmetry_space_group_name_Hall  '-B 2ab 2'
_symmetry_space_group_name_H-M   'B m a b'
_cell_length_a                   5.352
_cell_length_b                   5.4
_cell_length_c                   13.157
_cell_angle_alpha                90
_cell_angle_beta                 90
_cell_angle_gamma                90
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
La1 La3+ 0. 0.0068 0.3615
Cu1 Cu2+ 0. 0. 0.
O1 O2- 0. 0.4673 0.3169
O2 O2- 0.25 0.25 0.0066
"""


@pytest.mark.parametrize("cif_text", [
    CIF_OLD_TAGS_BMAB,
    CIF_OLD_TAGS_BMAB.replace("_symmetry_space_group_name_H-M   'B m a b'\n", ""),
])
def test_old_symmetry_tags_keep_setting(cif_text):
    """COD 1008541: _space_group_IT_number with the old tag 'B m a b'.

    Only the IT number was read, so the structure was taken in the
    standard setting C m c e and O2 got 16 positions instead of 8.
    """
    crystal = cryspy.str_to_globaln(cif_text).items[0]
    space_group = crystal.space_group
    assert space_group.it_number == 64
    assert space_group.it_coordinate_system_code == "a-cb"
    assert get_operators(space_group) == normalize_operators(OPERATORS_BMAB)
    multiplicity = {atom.label: atom.multiplicity
                    for atom in crystal.atom_site.items}
    assert multiplicity == {"La1": 8, "Cu1": 4, "O1": 8, "O2": 8}
