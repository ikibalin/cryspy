"""
Functions and constants to work with space group.

List of constants:
-------------------
ACCESIBLE_BRAVAIS_TYPE
ACCESIBLE_IT_COORDINATE_SYSTEM_CODE
ACCESIBLE_LAUE_CLASS
ACCESIBLE_CENTRING_TYPE
ACCESIBLE_CRYSTAL_SYSTEM
ACCESIBLE_NAME_HM_SHORT
ACCESIBLE_NAME_SCHOENFLIES
ACCESIBLE_NAME_HALL_SHORT
ACCESIBLE_REFERENCE_SETTING

DEFAULT_REFERENCE_TABLE_IT_NUMBER_NAME_HALL_NAME_SCHOENFLIES_NAME_HM_SHORT_REFERENCE_SETTING_IT_COORDINATE_SYSTEM_CODE

D_CENTRING_TYPE_SHIFT - accessible list and shift
D_CRYSTAL_FAMILY_DESCRIPTION - accessible list and description
D_BRAVAIS_TYPE_CELL_CONSTRAINT_MODE_ABC - accessible list and description constraint_mode_abc
T_BRAVAIS_TYPE_CENTRING_TYPE_CRYSTAL_SYSTEM - relation between bravais_type, centring_type, crystal_system
List of functions:
-------------------
get_crystal_system_by_it_number(it_number:int)->str
get_default_it_coordinate_system_code_by_it_number(it_number:int)->str
get_it_number_by_name_hm_short(name:str)->int
get_it_number_by_name_schoenflies(name:str)->int
get_it_number_by_name_hall(name:str)->int
get_name_hm_short_by_it_number(it_number:int)->str
get_name_schoenflies_by_it_number(it_number:int)->str
get_name_hall_by_it_number(it_number:int)->str


"""
import os
from numpy import array, transpose, zeros
from fractions import Fraction
from cryspy.A_functions_base.function_1_strings import \
    transform_string_to_r_b, transform_r_b_to_string

from typing import Tuple

F_ITABLES = os.path.join(os.path.dirname(__file__), "itables.txt")
F_WYCKOFF = os.path.join(os.path.dirname(__file__), "wyckoff.dat")

def read_el_cards():
    """
    Read information about space group from file to list of cards ldcard.

    Info in file fitables:

    1 P1               Triclinic
    choice: 1
    centr: false
    pcentr: 0, 0, 0
    symmetry: X,Y,Z

    2 P-1              Triclinic
    ...
    """
    fid = open(F_ITABLES, "r")
    lcontent = fid.readlines()
    fid.close()

    lcontent = [hh.strip() for hh in lcontent if hh.strip() != ""]
    ldcard = []
    dcard = None
    for hh in lcontent:
        lhelp = hh.split()
        if lhelp[0].isdigit():
            if dcard != None:
                ldcard.append(dcard)
            dcard = {"it_number": int(lhelp[0]), "name": lhelp[1], "singony": lhelp[2]}
        else:
            lhelp = hh.split(":")
            if (lhelp[0].strip() in dcard.keys()):
                dcard[lhelp[0].strip()].append(lhelp[1].strip())
            else:
                dcard[lhelp[0].strip()] = [lhelp[1].strip()]
    ldcard.append(dcard)
    return ldcard


EL_CARDS = read_el_cards()


def read_wyckoff():
    with open(F_WYCKOFF, "r") as fid:
        l_cont = fid.readlines()
    l_numb_b, l_numb_e = [], []
    for _i_line, _line in enumerate(l_cont):
        l_h = _line.strip().split()
        for _i, _ in enumerate(l_h):
            if not (_.isdigit()):
                break
        if _i >= 4:
            l_numb_b.append(_i_line)
        if len(l_h) == 0:
            l_numb_e.append(_i_line)
    l_data = []
    for _numb_b, _numb_e in zip(l_numb_b, l_numb_e):
        l_param = l_cont[_numb_b].strip().split()[:5]
        hm_full = ""
        flag = False
        for _char in l_cont[_numb_b].strip():
            if _char.isalpha():
                flag = True
            if flag:
                hm_full += _char
        data = {"it_number": int(l_param[0]), "choice": int(l_param[1]), "centr_000": int(l_param[3] == 1),
                "hm_full": hm_full.strip(), "wyckoff": []}
        l_cont_2 = l_cont[(_numb_b + 1):_numb_e]
        l_wyckoff_symop = []
        l_d_card = []
        d_card = None
        for _line in l_cont_2:
            l_h = _line.strip().split()
            if l_h[0].isdigit():
                if d_card is not None:
                    l_d_card.append(d_card)
                d_card = {"multiplicity": int(l_h[0]), "letter": l_h[1], "site_symmetry": l_h[2], "symop": []}
            else:
                d_card["symop"].extend(l_h)
        l_d_card.append(d_card)
        data["wyckoff"].extend(l_d_card)
        l_data.append(data)
    return l_data


WYCKOFF = read_wyckoff()


def get_crystal_system_by_it_number(it_number: int) -> str:
    if it_number is None:
        return None
    if (it_number >= 1) & (it_number <= 2):
        res = "triclinic"
    elif (it_number >= 3) & (it_number <= 15):
        res = "monoclinic"
    elif (it_number >= 16) & (it_number <= 74):
        res = "orthorhombic"
    elif (it_number >= 75) & (it_number <= 142):
        res = "tetragonal"
    elif (it_number >= 143) & (it_number <= 167):
        res = "trigonal"
    elif (it_number >= 168) & (it_number <= 194):
        res = "hexagonal"
    elif (it_number >= 195) & (it_number <= 230):
        res = "cubic"
    else:
        res = None
    return res


ACCESIBLE_IT_NUMBER_TRICLINIC_SYSTEM = tuple(range(1, 3))
ACCESIBLE_IT_NUMBER_MONOCLINIC_SYSTEM = tuple(range(3, 16))
ACCESIBLE_IT_NUMBER_ORTHORHOMBIC_SYSTEM = tuple(range(16, 75))
ACCESIBLE_IT_NUMBER_TETRAGONAL_SYSTEM = tuple(range(7, 143))
ACCESIBLE_IT_NUMBER_TRIGONAL_SYSTEM = tuple(range(143, 168))
ACCESIBLE_IT_NUMBER_HEXAGONAL_SYSTEM = tuple(range(168, 195))
ACCESIBLE_IT_NUMBER_CUBIC_SYSTEM = tuple(range(195, 231))

ACCESIBLE_IT_NUMBER_MONOCLINIC_SYSTEM_TRIPLE_CHOICE = (5, 7, 8, 9, 12, 13, 14, 15)
ACCESIBLE_IT_NUMBER_ORTHORHOMBIC_SYSTEM_DOUBLE_CHOICE = (48, 50, 59, 68, 70)
ACCESIBLE_IT_NUMBER_TETRAGONAL_SYSTEM_DOUBLE_CHOICE = (85, 86, 88, 125, 126, 129, 130, 133, 134, 137, 138, 141, 142)
ACCESIBLE_IT_NUMBER_TRIGONAL_SYSTEM_DOUBLE_AXES = (146, 148, 155, 160, 161, 166, 167)
ACCESIBLE_IT_NUMBER_CUBIC_SYSTEM_DOUBLE_CHOICE = (201, 203, 222, 224, 227, 228)

ACCESIBLE_IT_NUMBER = (ACCESIBLE_IT_NUMBER_TRICLINIC_SYSTEM +
                       ACCESIBLE_IT_NUMBER_MONOCLINIC_SYSTEM +
                       ACCESIBLE_IT_NUMBER_ORTHORHOMBIC_SYSTEM +
                       ACCESIBLE_IT_NUMBER_TETRAGONAL_SYSTEM +
                       ACCESIBLE_IT_NUMBER_TRIGONAL_SYSTEM +
                       ACCESIBLE_IT_NUMBER_HEXAGONAL_SYSTEM +
                       ACCESIBLE_IT_NUMBER_CUBIC_SYSTEM)


def get_default_it_coordinate_system_code_by_it_number(it_number: int) -> str:
    crystal_system = get_crystal_system_by_it_number(it_number)
    if crystal_system == "triclinic":
        it_coordinate_system_code = None
    elif crystal_system == "monoclinic":
        it_coordinate_system_code = "b1"
    elif crystal_system == "orthorhombic":
        if it_number in ACCESIBLE_IT_NUMBER_ORTHORHOMBIC_SYSTEM_DOUBLE_CHOICE:
            it_coordinate_system_code = "2abc"
        else:
            it_coordinate_system_code = "abc"
    elif crystal_system == "tetragonal":
        if it_number in ACCESIBLE_IT_NUMBER_TETRAGONAL_SYSTEM_DOUBLE_CHOICE:
            it_coordinate_system_code = "2"
        else:
            it_coordinate_system_code = "1"
    elif crystal_system == "trigonal":
        if it_number in ACCESIBLE_IT_NUMBER_TRIGONAL_SYSTEM_DOUBLE_AXES:
            it_coordinate_system_code = "h"
        else:
            it_coordinate_system_code = "r"
    elif crystal_system == "hexagonal":
        it_coordinate_system_code = "h"
    elif crystal_system == "cubic":
        if it_number in ACCESIBLE_IT_NUMBER_CUBIC_SYSTEM_DOUBLE_CHOICE:
            it_coordinate_system_code = "2"
        else:
            it_coordinate_system_code = "1"
    else:
        it_coordinate_system_code = None
    return it_coordinate_system_code


def get_it_coordinate_system_codes_by_it_number(it_number: int) -> str:
    crystal_system = get_crystal_system_by_it_number(it_number)
    if crystal_system == "triclinic":
        it_coordinate_system_codes = ()
    elif crystal_system == "monoclinic":
        it_coordinate_system_codes = (
        "b1", "c1", "a1", "b2", "c2", "a2", "b3", "c3", "a3", "-b1", "-c1", "-a1", "-b2", "-c2", "-a2", "-b3", "-c3",
        "-a3")
    elif crystal_system == "orthorhombic":
        if it_number in ACCESIBLE_IT_NUMBER_ORTHORHOMBIC_SYSTEM_DOUBLE_CHOICE:
            it_coordinate_system_codes = ("1abc", "1ba-c", "1cab", "1-cba",
                                          "1bca", "1a-cb", "2abc", "2ba-c", "2cab", "2-cba", "2bca", "2a-cb")
        else:
            it_coordinate_system_codes = ("abc", "ba-c", "cab", "-cba", "bca", "a-cb")
    elif crystal_system == "tetragonal":
        if it_number in ACCESIBLE_IT_NUMBER_TETRAGONAL_SYSTEM_DOUBLE_CHOICE:
            it_coordinate_system_codes = ("2", "1")
        else:
            it_coordinate_system_codes = ("1",)
    elif crystal_system == "trigonal":
        if it_number in ACCESIBLE_IT_NUMBER_TRIGONAL_SYSTEM_DOUBLE_AXES:
            it_coordinate_system_codes = ("h", "r")
        else:
            it_coordinate_system_codes = ("h",) 
    elif crystal_system == "hexagonal":
        it_coordinate_system_codes = ("h",)
    elif crystal_system == "cubic":
        if it_number in ACCESIBLE_IT_NUMBER_CUBIC_SYSTEM_DOUBLE_CHOICE:
            it_coordinate_system_codes = ("2", "1")
        else:
            it_coordinate_system_codes = ("1",)
    else:
        it_coordinate_system_codes = ()
    return it_coordinate_system_codes


ACCESIBLE_IT_COORDINATE_SYSTEM_CODE = ("b1", "b2", "b3", "-b1", "-b2", "-b3", "c1", "c2", "c3", "-c1", "-c2", "-c3",
                                       "a1", "a2", "a3", "-a1", "-a2", "-a3", "abc", "ba-c", "cab", "-cba", "bca",
                                       "a-cb", "1abc", "1ba-c", "1cab", "1-cba",
                                       "1bca", "1a-cb", "2abc", "2ba-c", "2cab", "2-cba", "2bca", "2a-cb", "1", "2",
                                       "h", "r")

ACCESIBLE_CRYSTAL_SYSTEM = ("triclinic", "monoclinic", "orthorhombic", "tetragonal", "trigonal", "hexagonal", "cubic")


def get_it_coordinate_system_codes_by_crystal_system(crystal_system: str) -> str:
    if crystal_system.startswith("tric"):
        it_coordinate_system_codes = ()
    elif crystal_system.startswith("m"):
        it_coordinate_system_codes = ("b1", "b2", "b3", "-b1", "-b2", "-b3", "c1", "c2", "c3", "-c1", "-c2", "-c3",
                                      "a1", "a2", "a3", "-a1", "-a2", "-a3")
    elif crystal_system.startswith("o"):
        it_coordinate_system_codes = ("abc", "ba-c", "cab", "-cba", "bca", "a-cb", "1abc", "1ba-c", "1cab", "1-cba",
                                      "1bca", "1a-cb", "2abc", "2ba-c", "2cab", "2-cba", "2bca", "2a-cb")
    elif crystal_system.startswith("te"):
        it_coordinate_system_codes = ("1", "2")
    elif crystal_system.startswith("trig"):
        it_coordinate_system_codes = ("h", "r")
    elif crystal_system.startswith("h"):
        it_coordinate_system_codes = ("h",)
    elif crystal_system.startswith("c"):
        it_coordinate_system_codes = ("1", "2")
    else:
        it_coordinate_system_codes = ()
    return it_coordinate_system_codes


ACCESIBLE_LAUE_CLASS = ("-1", "2/m", "mmm", "4/m", "4/mmm", "-3", "-3m", "6/m", "6/mmm", "m-3", "m-3m")

ACCESIBLE_CENTRING_TYPE = ("P", "A", "B", "C", "F", "I", "R", "Rrev", "H")

ACCESIBLE_NAME_HM_SHORT = ("P 1", "P -1", "P 2", "P 21", "C 2", "P m", "P c", "C m", "C c", "P 2/m", "P 21/m", "C 2/m",
                           "P 2/c", "P 21/c", "C 2/c", "P 2 2 2", "P 2 2 21", "P 21 21 2", "P 21 21 21", "C 2 2 21",
                           "C 2 2 2", "F 2 2 2", "I 2 2 2",
                           "I 21 21 21", "P m m 2", "P m c 21", "P c c 2", "P m a 2", "P c a 21", "P n c 2", "P m n 21",
                           "P b a 2", "P n a 21", "P n n 2",
                           "C m m 2", "C m c 21", "C c c 2", "A m m 2", "A e m 2", "A m a 2", "A e a 2", "F m m 2",
                           "F d d 2", "I m m 2", "I b a 2", "I m a 2",
                           "P m m m", "P n n n", "P c c m", "P b a n", "P m m a", "P n n a", "P m n a", "P c c a",
                           "P b a m", "P c c n", "P b c m", "P n n m",
                           "P m m n", "P b c n", "P b c a", "P n m a", "C m c m", "C m c e", "C m m m", "C c c m",
                           "C m m e", "C c c e", "F m m m", "F d d d",
                           "I m m m", "I b a m", "I b c a", "I m m a", "P 4", "P 41", "P 42", "P 43", "I 4", "I 41",
                           "P -4", "I -4", "P 4/m", "P 42/m", "P 4/n",
                           "P 42/n", "I 4/m", "I 41/a", "P 4 2 2", "P 4 21 2", "P 41 2 2", "P 41 21 2", "P 42 2 2",
                           "P 42 21 2", "P 43 2 2", "P 43 21 2", "I 4 2 2",
                           "I 41 2 2", "P 4 m m", "P 4 b m", "P 42 c m", "P 42 n m", "P 4 c c", "P 4 n c", "P 42 m c",
                           "P 42 b c", "I 4 m m", "I 4 c m", "I 41 m d",
                           "I 41 c d", "P -4 2 m", "P -4 2 c", "P -4 21 m", "P -4 21 c", "P -4 m 2", "P -4 c 2",
                           "P -4 b 2", "P -4 n 2", "I -4 m 2", "I -4 c 2",
                           "I -4 2 m", "I -4 2 d", "P 4/m m m", "P 4/m c c", "P 4/n b m", "P 4/n n c", "P 4/m b m",
                           "P 4/m n c", "P 4/n m m", "P 4/n c c", "P 42/m m c",
                           "P 42/m c m", "P 42/n b c", "P 42/n n m", "P 42/m b c", "P 42/m n m", "P 42/n m c",
                           "P 42/n c m", "I 4/m m m", "I 4/m c m", "I 41/a m d",
                           "I 41/a c d", "P 3", "P 31", "P 32", "R 3", "P -3", "R -3", "P 3 1 2", "P 3 2 1", "P 31 1 2",
                           "P 31 2 1", "P 32 1 2", "P 32 2 1", "R 3 2",
                           "P 3 m 1", "P 3 1 m", "P 3 c 1", "P 3 1 c", "R 3 m", "R 3 c", "P -3 1 m", "P -3 1 c",
                           "P -3 m 1", "P -3 c 1", "R -3 m", "R -3 c", "P 6", "P 61",
                           "P 65", "P 62", "P 64", "P 63", "P -6", "P 6/m ", "P 63/m", "P 6 2 2", "P 61 2 2",
                           "P 65 2 2", "P 62 2 2", "P 64 2 2", "P 63 2 2", "P 6 m m",
                           "P 6 c c", "P 63 c m", "P 63 m c", "P -6 m 2", "P -6 c 2", "P -6 2 m", "P -6 2 c",
                           "P 6/m m m", "P 6/m c c", "P 63/m c m", "P 63/m m c", "P 2 3",
                           "F 2 3", "I 2 3", "P 21 3", "I 21 3", "P m -3", "P n -3", "F m -3", "F d -3", "I m -3",
                           "P a -3", "I a -3", "P 4 3 2", "P 42 3 2", "F 4 3 2",
                           "F 41 3 2", "I 4 3 2", "P 43 3 2", "P 41 3 2", "I 41 3 2", "P -4 3 m", "F -4 3 m",
                           "I -4 3 m", "P -4 3 n", "F -4 3 c", "I -4 3 d", "P m -3 m",
                           "P n -3 n", "P m -3 n", "P n -3 m", "F m -3 m", "F m -3 c", "F d -3 m", "F d -3 c",
                           "I m -3 m", "I a -3 d")

ACCESIBLE_NAME_HM_FULL = ("P 1", "P -1", "P 2", "P 21", "C 2", "P m", "P c", "C m", "C c", "P 2/m", "P 21/m", "C 2/m",
                          "P 2/c", "P 21/c", "C 2/c", "P 2 2 2", "P 2 2 21", "P 21 21 2", "P 21 21 21", "C 2 2 21",
                          "C 2 2 2", "F 2 2 2", "I 2 2 2",
                          "I 21 21 21", "P m m 2", "P m c 21", "P c c 2", "P m a 2", "P c a 21", "P n c 2", "P m n 21",
                          "P b a 2", "P n a 21", "P n n 2",
                          "C m m 2", "C m c 21", "C c c 2", "A m m 2", "A e m 2", "A m a 2", "A e a 2", "F m m 2",
                          "F d d 2", "I m m 2", "I b a 2", "I m a 2",
                          "P m m m", "P n n n", "P c c m", "P b a n", "P m m a", "P n n a", "P m n a", "P c c a",
                          "P b a m", "P c c n", "P b c m", "P n n m",
                          "P m m n", "P b c n", "P b c a", "P n m a", "C m c m", "C m c e", "C m m m", "C c c m",
                          "C m m e", "C c c e", "F m m m", "F d d d",
                          "I m m m", "I b a m", "I b c a", "I m m a", "P 4", "P 41", "P 42", "P 43", "I 4", "I 41",
                          "P -4", "I -4", "P 4/m", "P 42/m", "P 4/n",
                          "P 42/n", "I 4/m", "I 41/a", "P 4 2 2", "P 4 21 2", "P 41 2 2", "P 41 21 2", "P 42 2 2",
                          "P 42 21 2", "P 43 2 2", "P 43 21 2", "I 4 2 2",
                          "I 41 2 2", "P 4 m m", "P 4 b m", "P 42 c m", "P 42 n m", "P 4 c c", "P 4 n c", "P 42 m c",
                          "P 42 b c", "I 4 m m", "I 4 c m", "I 41 m d",
                          "I 41 c d", "P -4 2 m", "P -4 2 c", "P -4 21 m", "P -4 21 c", "P -4 m 2", "P -4 c 2",
                          "P -4 b 2", "P -4 n 2", "I -4 m 2", "I -4 c 2",
                          "I -4 2 m", "I -4 2 d", "P 4/m m m", "P 4/m c c", "P 4/n b m", "P 4/n n c", "P 4/m b m",
                          "P 4/m n c", "P 4/n m m", "P 4/n c c", "P 42/m m c",
                          "P 42/m c m", "P 42/n b c", "P 42/n n m", "P 42/m b c", "P 42/m n m", "P 42/n m c",
                          "P 42/n c m", "I 4/m m m", "I 4/m c m", "I 41/a m d",
                          "I 41/a c d", "P 3", "P 31", "P 32", "R 3", "P -3", "R -3", "P 3 1 2", "P 3 2 1", "P 31 1 2",
                          "P 31 2 1", "P 32 1 2", "P 32 2 1", "R 3 2",
                          "P 3 m 1", "P 3 1 m", "P 3 c 1", "P 3 1 c", "R 3 m", "R 3 c", "P -3 1 m", "P -3 1 c",
                          "P -3 m 1", "P -3 c 1", "R -3 m", "R -3 c", "P 6", "P 61",
                          "P 65", "P 62", "P 64", "P 63", "P -6", "P 6/m ", "P 63/m", "P 6 2 2", "P 61 2 2", "P 65 2 2",
                          "P 62 2 2", "P 64 2 2", "P 63 2 2", "P 6 m m",
                          "P 6 c c", "P 63 c m", "P 63 m c", "P -6 m 2", "P -6 c 2", "P -6 2 m", "P -6 2 c",
                          "P 6/m m m", "P 6/m c c", "P 63/m c m", "P 63/m m c", "P 2 3",
                          "F 2 3", "I 2 3", "P 21 3", "I 21 3", "P m -3", "P n -3", "F m -3", "F d -3", "I m -3",
                          "P a -3", "I a -3", "P 4 3 2", "P 42 3 2", "F 4 3 2",
                          "F 41 3 2", "I 4 3 2", "P 43 3 2", "P 41 3 2", "I 41 3 2", "P -4 3 m", "F -4 3 m", "I -4 3 m",
                          "P -4 3 n", "F -4 3 c", "I -4 3 d", "P m -3 m",
                          "P n -3 n", "P m -3 n", "P n -3 m", "F m -3 m", "F m -3 c", "F d -3 m", "F d -3 c",
                          "I m -3 m", "I a -3 d")

ACCESIBLE_NAME_SCHOENFLIES = (
"C1.1", "Ci.1", "C2.1", "C2.2", "C2.3", "Cs.1", "Cs.2", "Cs.3", "Cs.4", "C2h.1", "C2h.2", "C2h.3", "C2h.4",
"C2h.5", "C2h.6", "D2.1", "D2.2", "D2.3", "D2.4", "D2.5", "D2.6", "D2.7", "D2.8", "D2.9", "C2v.1", "C2v.2", "C2v.3",
"C2v.4", "C2v.5",
"C2v.6", "C2v.7", "C2v.8", "C2v.9", "C2v.10", "C2v.11", "C2v.12", "C2v.13", "C2v.14", "C2v.15", "C2v.16", "C2v.17",
"C2v.18", "C2v.19",
"C2v.20", "C2v.21", "C2v.22", "D2h.1", "D2h.2", "D2h.3", "D2h.4", "D2h.5", "D2h.6", "D2h.7", "D2h.8", "D2h.9", "D2h.10",
"D2h.11", "D2h.12",
"D2h.13", "D2h.14", "D2h.15", "D2h.16", "D2h.17", "D2h.18", "D2h.19", "D2h.20", "D2h.21", "D2h.22", "D2h.23", "D2h.24",
"D2h.25", "D2h.26",
"D2h.27", "D2h.28", "C4.1", "C4.2", "C4.3", "C4.4", "C4.5", "C4.6", "S4.1", "S4.2", "C4h.1", "C4h.2", "C4h.3", "C4h.4",
"C4h.5", "C4h.6",
"D4.1", "D4.2", "D4.3", "D4.4", "D4.5", "D4.6", "D4.7", "D4.8", "D4.9", "D4.10", "C4v.1", "C4v.2", "C4v.3", "C4v.4",
"C4v.5", "C4v.6", "C4v.7",
"C4v.8", "C4v.9", "C4v.10", "C4v.11", "C4v.12", "D2d.1", "D2d.2", "D2d.3", "D2d.4", "D2d.5", "D2d.6", "D2d.7", "D2d.8",
"D2d.9", "D2d.10",
"D2d.11", "D2d.12", "D4h.1", "D4h.2", "D4h.3", "D4h.4", "D4h.5", "D4h.6", "D4h.7", "D4h.8", "D4h.9", "D4h.10", "D4h.11",
"D4h.12", "D4h.13",
"D4h.14", "D4h.15", "D4h.16", "D4h.17", "D4h.18", "D4h.19", "D4h.20", "C3.1", "C3.2", "C3.3", "C3.4", "C3i.1", "C3i.2",
"D3.1", "D3.2", "D3.3",
"D3.4", "D3.5", "D3.6", "D3.7", "C3v.1", "C3v.2", "C3v.3", "C3v.4", "C3v.5", "C3v.6", "D3d.1", "D3d.2", "D3d.3",
"D3d.4", "D3d.5", "D3d.6",
"C6.1", "C6.2", "C6.3", "C6.4", "C6.5", "C6.6", "C3h.1", "C6h.1", "C6h.2", "D6.1", "D6.2", "D6.3", "D6.4", "D6.5",
"D6.6", "C6v.1", "C6v.2",
"C6v.3", "C6v.4", "D3h.1", "D3h.2", "D3h.3", "D3h.4", "D6h.1", "D6h.2", "D6h.3", "D6h.4", "T.1", "T.2", "T.3", "T.4",
"T.5", "Th.1", "Th.2",
"Th.3", "Th.4", "Th.5", "Th.6", "Th.7", "O.1", "O.2", "O.3", "O.4", "O.5", "O.6", "O.7", "O.8", "Td.1", "Td.2", "Td.3",
"Td.4", "Td.5", "Td.6",
"Oh.1", "Oh.2", "Oh.3", "Oh.4", "Oh.5", "Oh.6", "Oh.7", "Oh.8", "Oh.9", "Oh.10")

ACCESIBLE_NAME_HALL_SHORT = (
"P 1", "-P 1", "P 2y", "P 2yb", "C 2y", "P -2y", "P -2yc", "C -2y", "C -2yc", "-P 2y", "-P 2yb", "-C 2y", "-P 2yc",
"-P 2ybc",
"-C 2yc", "P 2 2", "P 2c 2", "P 2 2ab", "P 2ac 2ab", "C 2c 2", "C 2 2", "F 2 2", "I 2 2", "I 2b 2c", "P 2 -2",
"P 2c -2", "P 2 -2c", "P 2 -2a", "P 2c -2ac",
"P 2 -2bc", "P 2ac -2", "P 2 -2ab", "P 2c -2n", "P 2 -2n", "C 2 -2", "C 2c -2", "C 2 -2c", "A 2 -2", "A 2 -2b",
"A 2 -2a", "A 2 -2ab", "F 2 -2",
"F 2 -2d", "I 2 -2", "I 2 -2c", "I 2 -2a", "-P 2 2", "-P 2ab 2bc", "-P 2 2c", "-P 2ab 2b", "-P 2a 2a", "-P 2a 2bc",
"-P 2ac 2", "-P 2a 2ac", "-P 2 2ab",
"-P 2ab 2ac", "-P 2c 2b", "-P 2 2n", "-P 2ab 2a", "-P 2n 2ab", "-P 2ac 2ab", "-P 2ac 2n", "-C 2c 2", "-C 2ac 2",
"-C 2 2", "-C 2 2c",
"-C 2a 2", "-C 2a 2ac", "-F 2 2", "-F 2uv 2vw", "-I 2 2", "-I 2 2c", "-I 2b 2c", "-I 2b 2", "P 4", "P 4w", "P 4c",
"P 4cw", "I 4", "I 4bw", "P -4", "I -4",
"-P 4", "-P 4c", "-P 4a", "-P 4bc", "-I 4", "-I 4ad", "P 4 2", "P 4ab 2ab", "P 4w 2c", "P 4abw 2nw", "P 4c 2",
"P 4n 2n", "P 4cw 2c", "P 4nw 2abw", "I 4 2",
"I 4bw 2bw", "P 4 -2", "P 4 -2ab", "P 4c -2c", "P 4n -2n", "P 4 -2c", "P 4 -2n", "P 4c -2", "P 4c -2ab", "I 4 -2",
"I 4 -2c", "I 4bw -2", "I 4bw -2c",
"P -4 2", "P -4 2c", "P -4 2ab", "P -4 2n", "P -4 -2", "P -4 -2c", "P -4 -2ab", "P -4 -2n", "I -4 -2", "I -4 -2c",
"I -4 2", "I -4 2bw", "-P 4 2",
"-P 4 2c", "-P 4a 2b", "-P 4a 2bc", "-P 4 2ab", "-P 4 2n", "-P 4a 2a", "-P 4a 2ac", "-P 4c 2", "-P 4c 2c", "-P 4ac 2b",
"-P 4ac 2bc", "-P 4c 2ab",
"-P 4n 2n", "-P 4ac 2a", "-P 4ac 2ac", "-I 4 2", "-I 4 2c", "-I 4bd 2", "-I 4bd 2c", "P 3", "P 31", "P 32", "R 3",
"-P 3", "-R 3",
"P 3 2", "P 3 2\"", "P 31 2 (0 0 4)", "P 31 2\"", "P 32 2 (0 0 2)", "P 32 2\"", "R 3 2\"", "P 3 -2\"", "P 3 -2",
"P 3 -2\"c", "P 3 -2c", "R 3 -2\"", "R 3 -2\"c",
"-P 3 2", "-P 3 2c", "-P 3 2\"", "-P 3 2\"c", "-R 3 2\"", "-R 3 2\"c", "P 6", "P 61", "P 65", "P 62", "P 64", "P 6c",
"P -6", "-P 6", "-P 6c", "P 6 2",
"P 61 2 (0 0 5)", "P 65 2 (0 0 1)", "P 62 2 (0 0 4)", "P 64 2 (0 0 2)", "P 6c 2c", "P 6 -2", "P 6 -2c", "P 6c -2",
"P 6c -2c", "P -6 2", "P -6c 2", "P -6 -2",
"P -6c -2c", "-P 6 2", "-P 6 2c", "-P 6c 2", "-P 6c 2c", "P 2 2 3", "F 2 2 3", "I 2 2 3", "P 2ac 2ab 3", "I 2b 2c 3",
"-P 2 2 3", "-P 2ab 2bc 3",
"-F 2 2 3", "-F 2uv 2vw 3", "-I 2 2 3", "-P 2ac 2ab 3", "-I 2b 2c 3", "P 4 2 3", "P 4n 2 3", "F 4 2 3", "F 4d 2 3",
"I 4 2 3", "P 4acd 2ab 3", "P 4bd 2ab 3",
"I 4bd 2c 3", "P -4 2 3", "F -4 2 3", "I -4 2 3", "P -4n 2 3", "F -4a 2 3", "I -4bd 2c 3", "-P 4 2 3", "-P 4a 2bc 3",
"-P 4n 2 3", "-P 4bc 2bc 3",
"-F 4 2 3", "-F 4a 2 3", "-F 4vw 2vw 3", "-F 4ud 2vw 3", "-I 4 2 3", "-I 4bd 2c 3")

ACCESIBLE_REFERENCE_SETTING = tuple(
    [f"{str(_1).zfill(3):}: {_2:}" for _1, _2 in zip(range(1, 231), ACCESIBLE_NAME_HALL_SHORT)])

DEFAULT_REFERENCE_TABLE_IT_NUMBER_NAME_HALL_NAME_SCHOENFLIES_NAME_HM_SHORT_REFERENCE_SETTING_IT_COORDINATE_SYSTEM_CODE = tuple(
    [
        (_1, _2, _3, _4, _5, get_default_it_coordinate_system_code_by_it_number(_1)) for _1, _2, _3, _4, _5 in
        zip(range(1, 231), ACCESIBLE_NAME_HALL_SHORT, ACCESIBLE_NAME_SCHOENFLIES, ACCESIBLE_NAME_HM_SHORT,
            ACCESIBLE_REFERENCE_SETTING)
    ])


def get_it_number_by_name_hm_short(name: str) -> int:
    if name in ACCESIBLE_NAME_HM_SHORT:
        it_number = ACCESIBLE_NAME_HM_SHORT.index(name) + 1
    else:
        it_number = None
    return it_number


def get_it_number_by_name_schoenflies(name: str) -> int:
    if (name in ACCESIBLE_NAME_SCHOENFLIES):
        it_number = ACCESIBLE_NAME_SCHOENFLIES.index(name) + 1
    else:
        it_number = None
    return it_number


def get_it_number_by_name_hall(name: str) -> int:
    if (name in ACCESIBLE_NAME_HALL_SHORT):
        it_number = ACCESIBLE_NAME_HALL_SHORT.index(name) + 1
    else:
        it_number = None
    return it_number


def get_name_hm_short_by_it_number(it_number: int) -> str:
    if (it_number in ACCESIBLE_IT_NUMBER):
        name = ACCESIBLE_NAME_HM_SHORT[it_number - 1]
    else:
        name = None
    return name


def get_name_schoenflies_by_it_number(it_number: int) -> str:
    if it_number in ACCESIBLE_IT_NUMBER:
        name = ACCESIBLE_NAME_SCHOENFLIES[it_number - 1]
    else:
        name = None
    return name


def get_name_hall_by_it_number(it_number: int) -> str:
    if it_number in ACCESIBLE_IT_NUMBER:
        name = ACCESIBLE_NAME_HALL_SHORT[it_number - 1]
    else:
        name = None
    return name

#FIXME it should be checked
REFERENCE_TABLE_TRICLINIC_IT_COORDINATE_SYSTEM_CODE_NAME_HM_EXTENDED = (
    (1, None, "P 1 1 1"), (2, None, "P 1 1 1")
)

# from IT A Table 4.3.2.1
REFERENCE_TABLE_MONOCLINIC_IT_COORDINATE_SYSTEM_CODE_NAME_HM_EXTENDED = (
    (3, "b1", "P 1 2 1"), (3, "-b1", "P 1 2 1"), (3, "c1", "P 1 1 2"), (3, "-c1", "P 1 1 2"), (3, "a1", "P 2 1 1"),
    (3, "-a1", "P 2 1 1"),
    (4, "b1", "P 1 21 1"), (4, "-b1", "P 1 21 1"), (4, "c1", "P 1 1 21"), (4, "-c1", "P 1 1 21"), (4, "a1", "P 21 1 1"),
    (4, "-a1", "P 21 1 1"),
    (5, "b1", "C 1 2 1"), (5, "-b1", "A 1 2 1"), (5, "c1", "A 1 1 2"), (5, "-c1", "B 1 1 2"), (5, "a1", "B 2 1 1"),
    (5, "-a1", "C 2 1 1"),
    (5, "b2", "A 1 2 1"), (5, "-b2", "C 1 2 1"), (5, "c2", "B 1 1 2"), (5, "-c2", "A 1 1 2"), (5, "a2", "C 2 1 1"),
    (5, "-a2", "B 2 1 1"),
    (5, "b3", "I 1 2 1"), (5, "-b3", "I 1 2 1"), (5, "c3", "I 1 1 2"), (5, "-c3", "I 1 1 2"), (5, "a3", "I 2 1 1"),
    (5, "-a3", "I 2 1 1"),
    (6, "b1", "P 1 m 1"), (6, "-b1", "P 1 m 1"), (6, "c1", "P 1 1 m"), (6, "-c1", "P 1 1 m"), (6, "a1", "P m 1 1"),
    (6, "-a1", "P m 1 1"),
    (7, "b1", "P 1 c 1"), (7, "-b1", "P 1 a 1"), (7, "c1", "P 1 1 a"), (7, "-c1", "P 1 1 b"), (7, "a1", "P b 1 1"),
    (7, "-a1", "P c 1 1"),
    (7, "b2", "P 1 n 1"), (7, "-b2", "P 1 n 1"), (7, "c2", "P 1 1 n"), (7, "-c2", "P 1 1 n"), (7, "a2", "P n 1 1"),
    (7, "-a2", "P n 1 1"),
    (7, "b3", "P 1 a 1"), (7, "-b3", "P 1 c 1"), (7, "c3", "P 1 1 b"), (7, "-c3", "P 1 1 a"), (7, "a3", "P c 1 1"),
    (7, "-a3", "P b 1 1"),
    (8, "b1", "C 1 m 1"), (8, "-b1", "A 1 m 1"), (8, "c1", "A 1 1 m"), (8, "-c1", "B 1 1 m"), (8, "a1", "B m 1 1"),
    (8, "-a1", "C m 1 1"),
    (8, "b2", "A 1 m 1"), (8, "-b2", "C 1 m 1"), (8, "c2", "B 1 1 m"), (8, "-c2", "A 1 1 m"), (8, "a2", "C m 1 1"),
    (8, "-a2", "B m 1 1"),
    (8, "b3", "I 1 m 1"), (8, "-b3", "I 1 m 1"), (8, "c3", "I 1 1 m"), (8, "-c3", "I 1 1 m"), (8, "a3", "I m 1 1"),
    (8, "-a3", "I m 1 1"),
    (9, "b1", "C 1 c 1"), (9, "-b1", "A 1 a 1"), (9, "c1", "A 1 1 a"), (9, "-c1", "B 1 1 b"), (9, "a1", "B b 1 1"),
    (9, "-a1", "C c 1 1"),
    (9, "b2", "A 1 n 1"), (9, "-b2", "C 1 n 1"), (9, "c2", "B 1 1 n"), (9, "-c2", "A 1 1 n"), (9, "a2", "C n 1 1"),
    (9, "-a2", "B n 1 1"),
    (9, "b3", "I 1 a 1"), (9, "-b3", "I 1 c 1"), (9, "c3", "I 1 1 b"), (9, "-c3", "I 1 1 a"), (9, "a3", "I c 1 1"),
    (9, "-a3", "I b 1 1"),
    (10, "b1", "P 1 2/m 1"), (10, "-b1", "P 1 2/m 1"), (10, "c1", "P 1 1 2/m"), (10, "-c1", "P 1 1 2/m"),
    (10, "a1", "P 2/m 1 1"), (10, "-a1", "P 2/m 1 1"),
    (11, "b1", "P 1 21/m 1"), (11, "-b1", "P 1 21/m 1"), (11, "c1", "P 1 1 21/m"), (11, "-c1", "P 1 1 21/m"),
    (11, "a1", "P 21/m 1 1"), (11, "-a1", "P 21/m 1 1"),
    (12, "b1", "C 1 2/m 1"), (12, "-b1", "A 1 2/m 1"), (12, "c1", "A 1 1 2/m"), (12, "-c1", "B 1 1 2/m"),
    (12, "a1", "B 2/m 1 1"), (12, "-a1", "C 2/m 1 1"),
    (12, "b2", "A 1 2/m 1"), (12, "-b2", "C 1 2/m 1"), (12, "c2", "B 1 1 2/m"), (12, "-c2", "A 1 1 2/m"),
    (12, "a2", "C 2/m 1 1"), (12, "-a2", "B 2/m 1 1"),
    (12, "b3", "I 1 2/m 1"), (12, "-b3", "I 1 2/m 1"), (12, "c3", "I 1 1 2/m"), (12, "-c3", "I 1 1 2/m"),
    (12, "a3", "I 2/m 1 1"), (12, "-a3", "I 2/m 1 1"),
    (13, "b1", "P 1 2/c 1"), (13, "-b1", "P 1 2/a 1"), (13, "c1", "P 1 1 2/a"), (13, "-c1", "P 1 1 2/b"),
    (13, "a1", "P 2/b 1 1"), (13, "-a1", "P 2/c 1 1"),
    (13, "b2", "P 1 2/n 1"), (13, "-b2", "P 1 2/n 1"), (13, "c2", "P 1 1 2/n"), (13, "-c2", "P 1 1 2/n"),
    (13, "a2", "P 2/n 1 1"), (13, "-a2", "P 2/n 1 1"),
    (13, "b3", "P 1 2/a 1"), (13, "-b3", "P 1 2/c 1"), (13, "c3", "P 1 1 2/b"), (13, "-c3", "P 1 1 2/a"),
    (13, "a3", "P 2/c 1 1"), (13, "-a3", "P 2/b 1 1"),
    (14, "b1", "P 1 21/c 1"), (14, "-b1", "P 1 21/a 1"), (14, "c1", "P 1 1 21/a"), (14, "-c1", "P 1 1 21/b"),
    (14, "a1", "P 21/b 1 1"), (14, "-a1", "P 21/c 1 1"),
    (14, "b2", "P 1 21/n 1"), (14, "-b2", "P 1 21/n 1"), (14, "c2", "P 1 1 21/n"), (14, "-c2", "P 1 1 21/n"),
    (14, "a2", "P 21/n 1 1"), (14, "-a2", "P 21/n 1 1"),
    (14, "b3", "P 1 21/a 1"), (14, "-b3", "P 1 21/c 1"), (14, "c3", "P 1 1 21/b"), (14, "-c3", "P 1 1 21/a"),
    (14, "a3", "P 21/c 1 1"), (14, "-a3", "P 21/b 1 1"),
    (15, "b1", "C 1 2/c 1"), (15, "-b1", "A 1 2/a 1"), (15, "c1", "A 1 1 2/a"), (15, "-c1", "B 1 1 2/b"),
    (15, "a1", "B 2/b 1 1"), (15, "-a1", "C 2/c 1 1"),
    (15, "b2", "A 1 2/n 1"), (15, "-b2", "C 1 2/n 1"), (15, "c2", "B 1 1 2/n"), (15, "-c2", "A 1 1 2/n"),
    (15, "a2", "C 2/n 1 1"), (15, "-a2", "B 2/n 1 1"),
    (15, "b3", "I 1 2/a 1"), (15, "-b3", "I 1 2/c 1"), (15, "c3", "I 1 1 2/b"), (15, "-c3", "I 1 1 2/a"),
    (15, "a3", "I 2/c 1 1"), (15, "-a3", "I 2/b 1 1"))

# from IT A Table 4.3.2.1
REFERENCE_TABLE_ORTHORHOMBIC_IT_COORDINATE_SYSTEM_CODE_NAME_HM_EXTENDED = (
    (16, "abc", "P 2 2 2"), (16, "ba-c", "P 2 2 2"), (16, "cab", "P 2 2 2"), (16, "-cba", "P 2 2 2"),
    (16, "bca", "P 2 2 2"), (16, "a-cb", "P 2 2 2"),
    (17, "abc", "P 2 2 21"), (17, "ba-c", "P 2 2 21"), (17, "cab", "P 21 2 2"), (17, "-cba", "P 21 2 2"),
    (17, "bca", "P 2 21 2"), (17, "a-cb", "P 2 21 2"),
    (18, "abc", "P 21 21 2"), (18, "ba-c", "P 21 21 2"), (18, "cab", "P 2 21 21"), (18, "-cba", "P 2 21 21"),
    (18, "bca", "P 21 2 21"), (18, "a-cb", "P 21 2 21"),
    (19, "abc", "P 21 21 21"), (19, "ba-c", "P 21 21 21"), (19, "cab", "P 21 21 21"), (19, "-cba", "P 21 21 21"),
    (19, "bca", "P 21 21 21"), (19, "a-cb", "P 21 21 21"),
    (20, "abc", "C 2 2 21"), (20, "ba-c", "C 2 2 21"), (20, "cab", "A 21 2 2"), (20, "-cba", "A 21 2 2"),
    (20, "bca", "B 2 21 2"), (20, "a-cb", "B 2 21 2"),
    (21, "abc", "C 2 2 2"), (21, "ba-c", "C 2 2 2"), (21, "cab", "A 2 2 2"), (21, "-cba", "A 2 2 2"),
    (21, "bca", "B 2 2 2"), (21, "a-cb", "B 2 2 2"),
    (22, "abc", "F 2 2 2"), (22, "ba-c", "F 2 2 2"), (22, "cab", "F 2 2 2"), (22, "-cba", "F 2 2 2"),
    (22, "bca", "F 2 2 2"), (22, "a-cb", "F 2 2 2"),
    (23, "abc", "I 2 2 2"), (23, "ba-c", "I 2 2 2"), (23, "cab", "I 2 2 2"), (23, "-cba", "I 2 2 2"),
    (23, "bca", "I 2 2 2"), (23, "a-cb", "I 2 2 2"),
    (24, "abc", "I 21 21 21"), (24, "ba-c", "I 21 21 21"), (24, "cab", "I 21 21 21"), (24, "-cba", "I 21 21 21"),
    (24, "bca", "I 21 21 21"), (24, "a-cb", "I 21 21 21"),
    (25, "abc", "P m m 2"), (25, "ba-c", "P m m 2"), (25, "cab", "P 2 m m"), (25, "-cba", "P 2 m m"),
    (25, "bca", "P m 2 m"), (25, "a-cb", "P m 2 m"),
    (26, "abc", "P m c 21"), (26, "ba-c", "P c m 21"), (26, "cab", "P 21 m a"), (26, "-cba", "P 21 a m"),
    (26, "bca", "P b 21 m"), (26, "a-cb", "P m 21 b"),
    (27, "abc", "P c c 2"), (27, "ba-c", "P c c 2"), (27, "cab", "P 2 a a"), (27, "-cba", "P 2 a a"),
    (27, "bca", "P b 2 b"), (27, "a-cb", "P b 2 b"),
    (28, "abc", "P m a 2"), (28, "ba-c", "P b m 2"), (28, "cab", "P 2 m b"), (28, "-cba", "P 2 c m"),
    (28, "bca", "P c 2 m"), (28, "a-cb", "P m 2 a"),
    (29, "abc", "P c a 21"), (29, "ba-c", "P b c 21"), (29, "cab", "P 21 a b"), (29, "-cba", "P 21 c a"),
    (29, "bca", "P c 21 b"), (29, "a-cb", "P b 21 a"),
    (30, "abc", "P n c 2"), (30, "ba-c", "P c n 2"), (30, "cab", "P 2 n a"), (30, "-cba", "P 2 a n"),
    (30, "bca", "P b 2 n"), (30, "a-cb", "P n 2 b"),
    (31, "abc", "P m n 21"), (31, "ba-c", "P n m 21"), (31, "cab", "P 21 m n"), (31, "-cba", "P 21 n m"),
    (31, "bca", "P n 21 m"), (31, "a-cb", "P m 21 n"),
    (32, "abc", "P b a 2"), (32, "ba-c", "P b a 2"), (32, "cab", "P 2 c b"), (32, "-cba", "P 2 c b"),
    (32, "bca", "P c 2 a"), (32, "a-cb", "P c 2 a"),
    (33, "abc", "P n a 21"), (33, "ba-c", "P b n 21"), (33, "cab", "P 21 n b"), (33, "-cba", "P 21 c n"),
    (33, "bca", "P c 21 n"), (33, "a-cb", "P n 21 a"),
    (34, "abc", "P n n 2"), (34, "ba-c", "P n n 2"), (34, "cab", "P 2 n n"), (34, "-cba", "P 2 n n"),
    (34, "bca", "P n 2 n"), (34, "a-cb", "P n 2 n"),
    (35, "abc", "C m m 2"), (35, "ba-c", "C m m 2"), (35, "cab", "A 2 m m"), (35, "-cba", "A 2 m m"),
    (35, "bca", "B m 2 m"), (35, "a-cb", "B m 2 m"),
    (36, "abc", "C m c 21"), (36, "ba-c", "C c m 21"), (36, "cab", "A 21 m a"), (36, "-cba", "A 21 a m"),
    (36, "bca", "B b 21 m"), (36, "a-cb", "B m 21 b"),
    (37, "abc", "C c c 2"), (37, "ba-c", "C c c 2"), (37, "cab", "A 2 a a"), (37, "-cba", "A 2 a a"),
    (37, "bca", "B b 2 b"), (37, "a-cb", "B b 2 b"),
    (38, "abc", "A m m 2"), (38, "ba-c", "B m m 2"), (38, "cab", "B 2 m m"), (38, "-cba", "C 2 m m"),
    (38, "bca", "C m 2 m"), (38, "a-cb", "A m 2 m"),
    (39, "abc", "A e m 2"), (39, "ba-c", "B m e 2"), (39, "cab", "B 2 e m"), (39, "-cba", "C 2 m e"),
    (39, "bca", "C m 2 e"), (39, "a-cb", "A e 2 m"),
    (40, "abc", "A m a 2"), (40, "ba-c", "B b m 2"), (40, "cab", "B 2 m b"), (40, "-cba", "C 2 c m"),
    (40, "bca", "C c 2 m"), (40, "a-cb", "A m 2 a"),
    (41, "abc", "A e a 2"), (41, "ba-c", "B b e 2"), (41, "cab", "B 2 e b"), (41, "-cba", "C 2 c e"),
    (41, "bca", "C c 2 e"), (41, "a-cb", "A e 2 a"),
    (42, "abc", "F m m 2"), (42, "ba-c", "F m m 2"), (42, "cab", "F 2 m m"), (42, "-cba", "F 2 m m"),
    (42, "bca", "F m 2 m"), (42, "a-cb", "F m 2 m"),
    (43, "abc", "F d d 2"), (43, "ba-c", "F d d 2"), (43, "cab", "F 2 d d"), (43, "-cba", "F 2 d d"),
    (43, "bca", "F d 2 d"), (43, "a-cb", "F d 2 d"),
    (44, "abc", "I m m 2"), (44, "ba-c", "I m m 2"), (44, "cab", "I 2 m m"), (44, "-cba", "I 2 m m"),
    (44, "bca", "I m 2 m"), (44, "a-cb", "I m 2 m"),
    (45, "abc", "I b a 2"), (45, "ba-c", "I b a 2"), (45, "cab", "I 2 c b"), (45, "-cba", "I 2 c b"),
    (45, "bca", "I c 2 a"), (45, "a-cb", "I c 2 a"),
    (46, "abc", "I m a 2"), (46, "ba-c", "I b m 2"), (46, "cab", "I 2 m b"), (46, "-cba", "I 2 c m"),
    (46, "bca", "I c 2 m"), (46, "a-cb", "I m 2 a"),
    (47, "abc", "P m m m"), (47, "ba-c", "P m m m"), (47, "cab", "P m m m"), (47, "-cba", "P m m m"),
    (47, "bca", "P m m m"), (47, "a-cb", "P m m m"),
    (48, "1abc", "P n n n"), (48, "2abc", "P n n n"), (48, "1ba-c", "P n n n"), (48, "2ba-c", "P n n n"),
    (48, "1cab", "P n n n"), (48, "2cab", "P n n n"),
    (48, "1-cba", "P n n n"), (48, "2-cba", "P n n n"), (48, "1bca", "P n n n"), (48, "2bca", "P n n n"),
    (48, "1a-cb", "P n n n"), (48, "2a-cb", "P n n n"),
    (49, "abc", "P c c m"), (49, "ba-c", "P c c m"), (49, "cab", "P m a a"), (49, "-cba", "P m a a"),
    (49, "bca", "P b m b"), (49, "a-cb", "P b m b"),
    (50, "1abc", "P b a n"), (50, "2abc", "P b a n"), (50, "1ba-c", "P b a n"), (50, "2ba-c", "P b a n"),
    (50, "1cab", "P n c b"), (50, "2cab", "P n c b"),
    (50, "1-cba", "P n c b"), (50, "2-cba", "P n c b"), (50, "1bca", "P c n a"), (50, "2bca", "P c n a"),
    (50, "1a-cb", "P c n a"), (50, "2a-cb", "P c n a"),
    (51, "abc", "P m m a"), (51, "ba-c", "P m m b"), (51, "cab", "P b m m"), (51, "-cba", "P c m m"),
    (51, "bca", "P m c m"), (51, "a-cb", "P m a m"),
    (52, "abc", "P n n a"), (52, "ba-c", "P n n b"), (52, "cab", "P b n n"), (52, "-cba", "P c n n"),
    (52, "bca", "P n c n"), (52, "a-cb", "P n a n"),
    (53, "abc", "P m n a"), (53, "ba-c", "P n m b"), (53, "cab", "P b m n"), (53, "-cba", "P c n m"),
    (53, "bca", "P n c m"), (53, "a-cb", "P m a n"),
    (54, "abc", "P c c a"), (54, "ba-c", "P c c b"), (54, "cab", "P b a a"), (54, "-cba", "P c a a"),
    (54, "bca", "P b c b"), (54, "a-cb", "P b a b"),
    (55, "abc", "P b a m"), (55, "ba-c", "P b a m"), (55, "cab", "P m c b"), (55, "-cba", "P m c b"),
    (55, "bca", "P c m a"), (55, "a-cb", "P c m a"),
    (56, "abc", "P c c n"), (56, "ba-c", "P c c n"), (56, "cab", "P n a a"), (56, "-cba", "P n a a"),
    (56, "bca", "P b n b"), (56, "a-cb", "P b n b"),
    (57, "abc", "P b c m"), (57, "ba-c", "P c a m"), (57, "cab", "P m c a"), (57, "-cba", "P m a b"),
    (57, "bca", "P b m a"), (57, "a-cb", "P c m b"),
    (58, "abc", "P n n m"), (58, "ba-c", "P n n m"), (58, "cab", "P m n n"), (58, "-cba", "P m n n"),
    (58, "bca", "P n m n"), (58, "a-cb", "P n m n"),
    (59, "1abc", "P m m n"), (59, "2abc", "P m m n"), (59, "1ba-c", "P m m n"), (59, "2ba-c", "P m m n"),
    (59, "1cab", "P n m m"), (59, "2cab", "P n m m"),
    (59, "1-cba", "P n m m"), (59, "2-cba", "P n m m"), (59, "1bca", "P m n m"), (59, "2bca", "P m n m"),
    (59, "1a-cb", "P m n m"), (59, "2a-cb", "P m n m"),
    (60, "abc", "P b c n"), (60, "ba-c", "P c a n"), (60, "cab", "P n c a"), (60, "-cba", "P n a b"),
    (60, "bca", "P b n a"), (60, "a-cb", "P c n b"),
    (61, "abc", "P b c a"), (61, "ba-c", "P c a b"), (61, "cab", "P b c a"), (61, "-cba", "P c a b"),
    (61, "bca", "P b c a"), (61, "a-cb", "P c a b"),
    (62, "abc", "P n m a"), (62, "ba-c", "P m n b"), (62, "cab", "P b n m"), (62, "-cba", "P c m n"),
    (62, "bca", "P m c n"), (62, "a-cb", "P n a m"),
    (63, "abc", "C m c m"), (63, "ba-c", "C c m m"), (63, "cab", "A m m a"), (63, "-cba", "A m a m"),
    (63, "bca", "B b m m"), (63, "a-cb", "B m m b"),
    (64, "abc", "C m c e"), (64, "ba-c", "C c m e"), (64, "cab", "A e m a"), (64, "-cba", "A e a m"),
    (64, "bca", "B b e m"), (64, "a-cb", "B m e b"),
    (65, "abc", "C m m m"), (65, "ba-c", "C m m m"), (65, "cab", "A m m m"), (65, "-cba", "A m m m"),
    (65, "bca", "B m m m"), (65, "a-cb", "B m m m"),
    (66, "abc", "C c c m"), (66, "ba-c", "C c c m"), (66, "cab", "A m a a"), (66, "-cba", "A m a a"),
    (66, "bca", "B b m b"), (66, "a-cb", "B b m b"),
    (67, "abc", "C m m e"), (67, "ba-c", "C m m e"), (67, "cab", "A e m m"), (67, "-cba", "A e m m"),
    (67, "bca", "B m e m"), (67, "a-cb", "B m e m"),
    (68, "1abc", "C c c e"), (68, "2abc", "C c c e"), (68, "1ba-c", "C c c e"), (68, "2ba-c", "C c c e"),
    (68, "1cab", "A e a a"), (68, "2cab", "A e a a"),
    (68, "1-cba", "A e a a"), (68, "2-cba", "A e a a"), (68, "1bca", "B b e b"), (68, "2bca", "B b e b"),
    (68, "1a-cb", "B b e b"), (68, "2a-cb", "B b e b"),
    (69, "abc", "F m m m"), (69, "ba-c", "F m m m"), (69, "cab", "F m m m"), (69, "-cba", "F m m m"),
    (69, "bca", "F m m m"), (69, "a-cb", "F m m m"),
    (70, "1abc", "F d d d"), (70, "2abc", "F d d d"), (70, "1ba-c", "F d d d"), (70, "2ba-c", "F d d d"),
    (70, "1cab", "F d d d"), (70, "2cab", "F d d d"),
    (70, "1-cba", "F d d d"), (70, "2-cba", "F d d d"), (70, "1bca", "F d d d"), (70, "2bca", "F d d d"),
    (70, "1a-cb", "F d d d"), (70, "2a-cb", "F d d d"),
    (71, "abc", "I m m m"), (71, "ba-c", "I m m m"), (71, "cab", "I m m m"), (71, "-cba", "I m m m"),
    (71, "bca", "I m m m"), (71, "a-cb", "I m m m"),
    (72, "abc", "I b a m"), (72, "ba-c", "I b a m"), (72, "cab", "I m c b"), (72, "-cba", "I m c b"),
    (72, "bca", "I c m a"), (72, "a-cb", "I c m a"),
    (73, "abc", "I b c a"), (73, "ba-c", "I c a b"), (73, "cab", "I b c a"), (73, "-cba", "I c a b"),
    (73, "bca", "I b c a"), (73, "a-cb", "I c a b"),
    (74, "abc", "I m m a"), (74, "ba-c", "I m m b"), (74, "cab", "I b m m"), (74, "-cba", "I c m m"),
    (74, "bca", "I m c m"), (74, "a-cb", "I m a m")
)

# from IT A Table 4.3.2.1
REFERENCE_TABLE_TETRAGONAL_IT_COORDINATE_SYSTEM_CODE_NAME_HM_EXTENDED = (
    (79, "1", "I 4"), (80, "1", "I 41"), (87, "1", "I 4/m"), (88, "1", "I 41/a"), (88, "2", "I 41/a"),
    (89, "1", "P 4 2 2"),
    (90, "1", "P 4 21 2"), (91, "1", "P 41 2 2"), (92, "1", "P 41 21 2"), (93, "1", "P 42 2 2"), (94, "1", "P 42 21 2"),
    (95, "1", "P 43 2 2"),
    (96, "1", "P 43 21 2"), (97, "1", "I 4 2 2"), (98, "1", "I 41 2 2"), (99, "1", "P 4 m m"), (100, "1", "P 4 b m"),
    (101, "1", "P 42 c m"),
    (102, "1", "P 42 n m"), (103, "1", "P 4 c c"), (104, "1", "P 4 n c"), (105, "1", "P 42 m c"),
    (106, "1", "P 42 b c"), (107, "1", "I 4 m m"),
    (108, "1", "I 4 c e"), (109, "1", "I 41 m d"), (110, "1", "I 41 c d"), (111, "1", "P -4 2 m"),
    (112, "1", "P -4 2 c"), (113, "1", "P -4 21 m"),
    (114, "1", "P -4 21 c"), (115, "1", "P -4 m 2"), (116, "1", "P -4 c 2"), (117, "1", "P -4 b 2"),
    (118, "1", "P -4 n 2"), (119, "1", "I -4 m 2"),
    (120, "1", "I -4 c 2"), (121, "1", "I -4 2 m"), (122, "1", "I -4 2 d"), (123, "1", "P 4/m 2/m 2/m"),
    (124, "1", "P 4/m 2/c 2/c"), (125, "1", "P 4/n 2/b 2/m"),
    (125, "2", "P 4/n 2/b 2/m"), (126, "1", "P 4/n 2/n 2/c"), (126, "2", "P 4/n 2/n 2/c"), (127, "1", "P 4/m 21/b 2/m"),
    (128, "1", "P 4/m 21/n 2/c"), (129, "1", "P 4/n 21/m 2/m"),
    (129, "2", "P 4/n 21/m 2/m"), (130, "1", "P 4/n 21/c 2/c"), (130, "2", "P 4/n 21/c 2/c"),
    (131, "1", "P 42/m 2/m 2/c"), (132, "1", "P 42/m 2/c 2/m"), (133, "1", "P 42/n 2/b 2/c"),
    (133, "2", "P 42/n 2/b 2/c"), (134, "1", "P 42/n 2/n 2/m"), (134, "2", "P 42/n 2/n 2/m"),
    (135, "1", "P 42/m 21/b 2/c"), (136, "1", "P 42/m 21/n 2/m"), (137, "1", "P 42/n 21/m 2/c"),
    (137, "2", "P 42/n 21/m 2/c"), (138, "1", "P 42/n 21/c 2/m"), (138, "2", "P 42/n 21/c 2/m"),
    (139, "1", "I 4/m 21/m 2/m"), (140, "1", "I 4/m 2/c 2/m"), (141, "1", "I 41/a 2/m 2/d"),
    (141, "2", "I 41/a 2/m 2/d"), (142, "1", "I 41/a 2/c 2/d"), (142, "2", "I 41/a 2/c 2/d")
)

# from IT A Table 4.3.2.1
REFERENCE_TABLE_TRIGONAL_IT_COORDINATE_SYSTEM_CODE_NAME_HM_EXTENDED = (
    (146, "r", "R 3"), (146, "h", "R 3"), (148, "r", "R -3"), (148, "h", "R -3"), (149, "r", "P 3 1 2"),
    (150, "r", "P 3 2 1"),
    (151, "r", "P 31 1 2"), (152, "r", "P 31 2 1"), (153, "r", "P 32 1 2"), (154, "r", "P 32 2 1"), (155, "r", "R 3 2"),
    (155, "h", "R 3 2"),
    (156, "r", "P 3 m 1"), (157, "r", "P 3 1 m"), (158, "r", "P 3 c 1"), (159, "r", "P 3 1 c"), (160, "r", "R 3 m"),
    (160, "h", "R 3 m"),
    (161, "r", "R 3 c"), (161, "h", "R 3 c"), (162, "r", "P -3 1 2/m"), (163, "r", "P -3 1 2/c"),
    (164, "r", "P -3 2/m 1"), (165, "r", "P -3 2/c 1"),
    (166, "r", "R -3 2/m"), (166, "h", "R -3 2/m"), (167, "r", "R -3 2/c"), (167, "h", "R -3 2/c")
)

# from IT A Table 4.3.2.1
REFERENCE_TABLE_HEXAGONAL_IT_COORDINATE_SYSTEM_CODE_NAME_HM_EXTENDED = (
    (177, "h", "P 6 2 2"), (178, "h", "P 61 2 2"), (179, "h", "P 65 2 2"), (180, "h", "P 62 2 2"),
    (181, "h", "P 64 2 2"), (182, "h", "P 63 2 2"),
    (183, "h", "P 6 m m"), (184, "h", "P 6 c c"), (185, "h", "P 63 c m"), (186, "h", "P 63 m c"),
    (187, "h", "P -6 m 2"), (188, "h", "P -6 c 2"),
    (189, "h", "P -6 2 m"), (190, "h", "P -6 2 c"), (191, "h", "P 6/m 2/m 2/m"), (192, "h", "P 6/m 2/c 2/c"),
    (193, "h", "P 63/m 2/c 2/m"), (194, "h", "P 63/m 2/m 2/c")
)

REFERENCE_TABLE_CUBIC_IT_COORDINATE_SYSTEM_CODE_NAME_HM_EXTENDED = (
    (196, "1", "F 2 3"), (197, "1", "I 2 3"), (199, "1", "I 21 3"), (202, "1", "F 2/m -3"), (203, "1", "F 2/d -3"),
    (203, "2", "F 2/d -3"),
    (204, "1", "I 2/m -3"), (206, "1", "I 21/a -3"), (207, "1", "P 4 3 2"), (208, "1", "P 42 3 2"),
    (209, "1", "F 4 3 2"), (210, "1", "F 41 3 2"),
    (211, "1", "I 4 3 2"), (212, "1", "P 43 3 2"), (213, "1", "P 41 3 2"), (214, "1", "I 41 3 2"),
    (215, "1", "P -4 3 m"), (216, "1", "F -4 3 m"),
    (217, "1", "I -4 3 m"), (218, "1", "P -4 3 n"), (219, "1", "F -4 3 c"), (220, "1", "I -4 3 d"),
    (221, "1", "P 4/m -3 2/m"), (222, "1", "P 4/n -3 2/n"),
    (222, "2", "P 4/n -3 2/n"), (223, "1", "P 42/m -3 2/n"), (224, "1", "P 42/n -3 2/m"), (224, "2", "P 42/n -3 2/m"),
    (225, "1", "F 4/m -3 2/m"), (226, "1", "F 4/m -3 2/c"),
    (227, "1", "F 41/d -3 2/m"), (227, "2", "F 41/d -3 2/m"), (228, "1", "F 41/d -3 2/n"), (228, "2", "F 41/d -3 2/n"),
    (229, "1", "I 4/m -3 2/m"), (230, "1", "I 41/a -3 2/d")
)

REFERENCE_TABLE_IT_COORDINATE_SYSTEM_CODE_NAME_HM_EXTENDED = (
            REFERENCE_TABLE_TRICLINIC_IT_COORDINATE_SYSTEM_CODE_NAME_HM_EXTENDED+
            REFERENCE_TABLE_MONOCLINIC_IT_COORDINATE_SYSTEM_CODE_NAME_HM_EXTENDED +
            REFERENCE_TABLE_ORTHORHOMBIC_IT_COORDINATE_SYSTEM_CODE_NAME_HM_EXTENDED +
            REFERENCE_TABLE_TETRAGONAL_IT_COORDINATE_SYSTEM_CODE_NAME_HM_EXTENDED +
            REFERENCE_TABLE_TRIGONAL_IT_COORDINATE_SYSTEM_CODE_NAME_HM_EXTENDED +
            REFERENCE_TABLE_HEXAGONAL_IT_COORDINATE_SYSTEM_CODE_NAME_HM_EXTENDED +
            REFERENCE_TABLE_CUBIC_IT_COORDINATE_SYSTEM_CODE_NAME_HM_EXTENDED)

ACCESIBLE_NAME_HM_EXTENDED = frozenset([_[2] for _ in REFERENCE_TABLE_IT_COORDINATE_SYSTEM_CODE_NAME_HM_EXTENDED])


def get_it_number_it_coordinate_system_codes_by_name_hm_extended(name: str) -> int:
    flag = True
    it_number = None
    it_coordinate_system_codes = []
    for _it_number, _it_coordinate_system_code, _name in REFERENCE_TABLE_IT_COORDINATE_SYSTEM_CODE_NAME_HM_EXTENDED:
        if name == _name:
            if it_number is not None:
                flag &= it_number == _it_number
            it_number = _it_number
            it_coordinate_system_codes.append(_it_coordinate_system_code)
    if not (flag):
        print(f"For some reason for hm_name_extended \"{name:}\" it_number is not unique")
    return it_number, tuple(it_coordinate_system_codes)


def get_name_hm_extended_by_it_number_it_coordinate_system_code(it_number: int, it_coordinate_system_code) -> str:
    name_hm_extended = None
    for _it_number, _it_coordinate_system_code, _name in REFERENCE_TABLE_IT_COORDINATE_SYSTEM_CODE_NAME_HM_EXTENDED:
        if ((it_number == _it_number) & (it_coordinate_system_code == _it_coordinate_system_code)):
            name_hm_extended = _name
            break
    return name_hm_extended


# IT A Table 12.3.4.1. Standard space-group symbols
REFERENCE_TABLE_TRICLINIC_IT_NUMBER_NAME_HM_FULL = (
    (1, "P 1"), (2, "P -1")
)

REFERENCE_TABLE_MONOCLINIC_IT_COORDINATE_SYSTEM_CODE_NAME_HM_FULL = (
    (3, "b1", "P 1 2 1"), (3, "c1", "P 1 1 2"), (4, "b1", "P 1 21 1"), (4, "c1", "P 1 1 21"),
    (5, "b1", "C 1 2 1"), (5, "c1", "A 1 1 2"), (5, "-c1", "B 1 1 2"), (6, "b1", "P 1 m 1"), (6, "c1", "P 1 1 m"),
    (7, "b1", "P 1 c 1"),
    (7, "c1", "P 1 1 a"), (7, "-c1", "P 1 1 b"), (8, "b1", "C 1 m 1"), (8, "c1", "A 1 1 m"), (8, "-c1", "B 1 1 m"),
    (9, "b1", "C 1 c 1"),
    (9, "c1", "A 1 1 a"), (9, "-c1", "B 1 1 b"), (10, "b1", "P 1 2/m 1"), (10, "c1", "P 1 1 2/m"),
    (11, "b1", "P 1 21/m 1"), (11, "c1", "P 1 1 21/m"),
    (12, "b1", "C 1 2/m 1"), (12, "c1", "A 1 1 2/m"), (12, "-c1", "B 1 1 2/m"), (13, "b1", "P 1 2/c 1"),
    (13, "c1", "P 1 1 2/a"), (13, "-c1", "P 1 1 2/b"),
    (14, "b1", "P 1 21/c 1"), (14, "c1", "P 1 1 21/a"), (14, "-c1", "P 1 1 21/b"), (15, "b1", "C 1 2/c 1"),
    (15, "c1", "A 1 1 2/a"), (15, "-c1", "B 1 1 2/b"))

REFERENCE_TABLE_ORTHORHOMBIC_IT_NUMBER_NAME_HM_FULL = (
    (16, "P 2 2 2"), (17, "P 2 2 21"), (18, "P 21 21 2"), (19, "P 21 21 21"), (20, "C 2 2 21"), (21, "C 2 2 2"),
    (22, "F 2 2 2"), (23, "I 2 2 2"), (24, "I 21 21 21"), (25, "P m m 2"), (26, "P m c 21"), (27, "P c c 2"),
    (28, "P m a 2"), (29, "P c a 21"), (30, "P n c 2"), (31, "P m n 21"), (32, "P b a 2"), (33, "P n a 21"),
    (34, "P n n 2"), (35, "C m m 2"), (36, "C m c 21"), (37, "C c c 2"), (38, "A m m 2"), (39, "A e m 2"),
    (40, "A m a 2"), (41, "A e a 2"), (42, "F m m 2"), (43, "F d d 2"), (44, "I m m 2"), (45, "I b a 2"),
    (46, "I m a 2"), (47, "P 2/m 2/m 2/m"), (48, "P 2/n 2/n 2/n"), (49, "P 2/c 2/c 2/m"), (50, "P 2/b 2/a 2/n"),
    (51, "P 21/m 2/m 2/a"),
    (52, "P 2/n 21/n 2/a"), (53, "P 2/m 2/n 21/a"), (54, "P 21/c 2/c 2/a"), (55, "P 21/b 21/a 2/m"),
    (56, "P 21/c 21/c 2/n"), (57, "P 2/b 21/c 21/m"),
    (58, "P 21/n 21/n 2/m"), (59, "P 21/m 21/m 2/n"), (60, "P 21/b 2/c 21/n"), (61, "P 21/b 21/c 21/a"),
    (62, "P 21/n 21/m 21/a"), (63, "C 2/m 2/c 21/m"),
    (64, "C 2/m 2/c 21/e"), (65, "C 2/m 2/m 2/m"), (66, "C 2/c 2/c 2/m"), (67, "C 2/m 2/m 2/e"), (68, "C 2/c 2/c 2/e"),
    (69, "F 2/m 2/m 2/m"),
    (70, "F 2/d 2/d 2/d"), (71, "I 2/m 2/m 2/m"), (72, "I 2/b 2/a 2/m"), (73, "I 21/b 21/c 21/a"),
    (74, "I 21/m 21/m 21/a")
)

REFERENCE_TABLE_TETRAGONAL_IT_NUMBER_NAME_HM_FULL = (
    (75, "P 4"), (76, "P 41"), (77, "P 42"), (78, "P 43"), (79, "I 4"), (80, "I 41"),
    (81, "P -4"), (82, "I -4"), (83, "P 4/m"), (84, "P 42/m"), (85, "P 4/n"), (86, "P 42/n"),
    (87, "I 4/m"), (88, "I 41/a"), (89, "P 4 2 2"), (90, "P 4 21 2"), (91, "P 41 2 2"), (92, "P 41 21 2"),
    (93, "P 42 2 2"), (94, "P 42 21 2"), (95, "P 43 2 2"), (96, "P 43 21 2"), (97, "I 4 2 2"), (98, "I 41 2 2"),
    (99, "P 4 m m"), (100, "P 4 b m"), (101, "P 42 c m"), (102, "P 42 n m"), (103, "P 4 c c"), (104, "P 4 n c"),
    (105, "P 42 m c"), (106, "P 42 b c"), (107, "I 4 m m"), (108, "I 4 c m"), (109, "I 41 m d"), (110, "I 41 c d"),
    (111, "P -4 2 m"), (112, "P -4 2 c"), (113, "P -4 21 m"), (114, "P -4 21 c"), (115, "P -4 m 2"), (116, "P -4 c 2"),
    (117, "P -4 b 2"), (118, "P -4 n 2"), (119, "I -4 m 2"), (120, "I -4 c 2"), (121, "I -4 2 m"), (122, "I -4 2 d"),
    (123, "P 4/m 2/m 2/m"), (124, "P 4/m 2/c 2/c"), (125, "P 4/n 2/b 2/m"), (126, "P 4/n 2/n 2/c"),
    (127, "P 4/m 21/b 2/m"), (128, "P 4/m 21/n 2/c"),
    (129, "P 4/n 21/m 2/m"), (130, "P 4/n 21/c 2/c"), (131, "P 42/m 2/m 2/c"), (132, "P 42/m 2/c 2/m"),
    (133, "P 42/n 2/b 2/c"), (134, "P 42/n 2/n 2/m"),
    (135, "P 42/m 21/b 2/c"), (136, "P 42/m 21/n 2/m"), (137, "P 42/n 21/m 2/c"), (138, "P 42/n 21/c 2/m"),
    (139, "I 4/m 21/m 2/m"), (140, "I 4/m 2/c 2/m"),
    (141, "I 41/a 2/m 2/d"), (142, "I 41/a 2/c 2/d")
)

REFERENCE_TABLE_TRIGONAL_IT_NUMBER_NAME_HM_FULL = (
    (143, "P 3"), (144, "P 31"), (145, "P 32"), (146, "R 3"), (147, "P -3"), (148, "R -3"),
    (149, "P 3 1 2"), (150, "P 3 2 1"), (151, "P 31 1 2"), (152, "P 31 2 1"), (153, "P 32 1 2"), (154, "P 32 2 1"),
    (155, "R 3 2"), (156, "P 3 m 1"), (157, "P 3 1 m"), (158, "P 3 c 1"), (159, "P 3 1 c"), (160, "R 3 m"),
    (161, "R 3 c"), (162, "P -3 1 2/m"), (163, "P -3 1 2/c"), (164, "P -3 2/m 1"), (165, "P -3 2/c 1"),
    (166, "R -3 2/m"),
    (167, "R -3 2/c")
)

REFERENCE_TABLE_HEXAGONAL_IT_NUMBER_NAME_HM_FULL = (
    (168, "P 6"), (169, "P 61"), (170, "P 65"), (171, "P 62"), (172, "P 64"), (173, "P 63"),
    (174, "P -6"), (175, "P 6/m "), (176, "P 63/m"), (177, "P 6 2 2"), (178, "P 61 2 2"), (179, "P 65 2 2"),
    (180, "P 62 2 2"), (181, "P 64 2 2"), (182, "P 63 2 2"), (183, "P 6 m m"), (184, "P 6 c c"), (185, "P 63 c m"),
    (186, "P 63 m c"), (187, "P -6 m 2"), (188, "P -6 c 2"), (189, "P -6 2 m"), (190, "P -6 2 c"),
    (191, "P 6/m 2/m 2/m"), (192, "P 6/m 2/c 2/c"), (193, "P 63/m 2/c 2/m"), (194, "P 63/m 2/m 2/c")
)

REFERENCE_TABLE_CUBIC_IT_NUMBER_NAME_HM_FULL = (
    (195, "P 23"), (196, "F 23"), (197, "I 23"), (198, "P 21 3"), (199, "I 21 3"),
    (200, "P 2/m -3"), (201, "P 2/n -3"), (202, "F 2/m -3"), (203, "F 2/d -3"), (204, "I 2/m -3"), (205, "P 21/a -3"),
    (206, "I 21/a -3"),
    (207, "P 4 3 2"), (208, "P 42 3 2"), (209, "F 4 3 2"), (210, "F 41 3 2"), (211, "I 4 3 2"), (212, "P 43 3 2"),
    (213, "P 41 3 2"), (214, "I 41 3 2"),
    (215, "P -4 3 m"), (216, "F -4 3 m"), (217, "I -4 3 m"), (218, "P -4 3 n"), (219, "F -4 3 c"), (220, "I -4 3 d"),
    (221, "P 4/m -3 2/m"), (222, "P 4/n -3 2/n"), (223, "P 42/m -3 2/n"), (224, "P 42/n -3 2/m"), (225, "F 4/m -3 2/m"),
    (226, "F 4/m -3 2/c"),
    (227, "F 41/d -3 2/m"), (228, "F 41/d -3 2/c"), (229, "I 4/m -3 2/m"), (230, "I 41/a -3 2/d")
)

REFERENCE_TABLE_IT_NUMBER_NAME_HM_FULL = (REFERENCE_TABLE_TRICLINIC_IT_NUMBER_NAME_HM_FULL +
                                          REFERENCE_TABLE_ORTHORHOMBIC_IT_NUMBER_NAME_HM_FULL +
                                          REFERENCE_TABLE_TETRAGONAL_IT_NUMBER_NAME_HM_FULL +
                                          REFERENCE_TABLE_TRIGONAL_IT_NUMBER_NAME_HM_FULL +
                                          REFERENCE_TABLE_HEXAGONAL_IT_NUMBER_NAME_HM_FULL +
                                          REFERENCE_TABLE_CUBIC_IT_NUMBER_NAME_HM_FULL)

ACCESIBLE_NAME_HM_FULL = frozenset([_[1] for _ in REFERENCE_TABLE_IT_NUMBER_NAME_HM_FULL])


def get_it_number_by_name_hm_full(name: str) -> int:
    _l = [_it_number for _it_number, _name in REFERENCE_TABLE_IT_NUMBER_NAME_HM_FULL if (name == _name)]
    if len(_l) == 0:
        it_number = None
    else:
        it_number = _l[0]
    return it_number


def get_name_hm_full_by_it_number(it_number: int) -> str:
    _l = [_name for _it_number, _name in REFERENCE_TABLE_IT_NUMBER_NAME_HM_FULL if (it_number == _it_number)]
    if len(_l) == 0:
        name_hm_full = None
    else:
        name_hm_full = _l[0]
    return name_hm_full


REFERENCE_TABLE_CENTRING_TYPE_SHIFT = (
    ("P", ((Fraction(0, 1), Fraction(0, 1), Fraction(0, 1)),)),
    ("A", ((Fraction(0, 1), Fraction(0, 1), Fraction(0, 1)),
           (Fraction(0, 2), Fraction(1, 2), Fraction(1, 2)))),
    ("B", ((Fraction(0, 1), Fraction(0, 1), Fraction(0, 1)),
           (Fraction(1, 2), Fraction(0, 2), Fraction(1, 2)))),
    ("C", ((Fraction(0, 1), Fraction(0, 1), Fraction(0, 1)),
           (Fraction(1, 2), Fraction(1, 2), Fraction(0, 2)))),
    ("F", ((Fraction(0, 1), Fraction(0, 1), Fraction(0, 1)),
           (Fraction(0, 2), Fraction(1, 2), Fraction(1, 2)),
           (Fraction(1, 2), Fraction(0, 2), Fraction(1, 2)),
           (Fraction(1, 2), Fraction(1, 2), Fraction(0, 2)))),
    ("I", ((Fraction(0, 1), Fraction(0, 1), Fraction(0, 1)),
           (Fraction(1, 2), Fraction(1, 2), Fraction(1, 2)))),
    ("R", ((Fraction(0, 1), Fraction(0, 1), Fraction(0, 1)),
           (Fraction(2, 3), Fraction(1, 3), Fraction(1, 3)),
           (Fraction(1, 3), Fraction(2, 3), Fraction(2, 3)))),
    ("Rrev", ((Fraction(0, 1), Fraction(0, 1), Fraction(0, 1)),
              (Fraction(1, 3), Fraction(2, 3), Fraction(1, 3)),
              (Fraction(2, 3), Fraction(1, 3), Fraction(2, 3)))),
    ("H", ((Fraction(0, 1), Fraction(0, 1), Fraction(0, 1)),
           (Fraction(2, 3), Fraction(1, 3), Fraction(0, 3)),
           (Fraction(1, 3), Fraction(2, 3), Fraction(0, 3))))
)

ACCESIBLE_CENTRING_TYPE = frozenset([_[0] for _ in REFERENCE_TABLE_CENTRING_TYPE_SHIFT])


def get_shift_by_centring_type(centring_type: str):
    shift = ()
    for _1, _2 in REFERENCE_TABLE_CENTRING_TYPE_SHIFT:
        if _1 == centring_type:
            shift = _2
            break
    return shift


def get_centring_type_by_name_hm_extended(hm_extended: str) -> str:
    centring_type = hm_extended[0]  # it is not correct for Rrev
    if not (centring_type in ACCESIBLE_CENTRING_TYPE):
        centring_type = None
    return centring_type


ACCESIBLE_LATTICE_TYPE = ("P", "C", "I", "F", "R")


def get_lattice_type_by_name_hm_short(hm_short: str) -> str:
    lattice_type = hm_short[0]
    if not (lattice_type in ACCESIBLE_LATTICE_TYPE):
        lattice_type = None
    return lattice_type


REFERENCE_TABLE_PATTERSON_NAME_HM_LATTICE_TYPE_LAUE_CLASS = (
    ("P -1", "P", "-1"), ("P 2/m", "P", "2/m"), ("C 2/m", "C", "2/m"), ("P m m m", "P", "mmm"), ("C m m m", "C", "mmm"),
    ("I m m m", "I", "mmm"),
    ("F m m m", "F", "mmm"), ("P 4/m", "P", "4/m"), ("I 4/m", "I", "4/m"), ("P 4/m m m", "P", "4/mmm"),
    ("I 4/m m m", "I", "4/mmm"), ("P -3", "P", "-3"),
    ("R -3", "R", "-3"), ("P -3 m 1", "P", "-3m1"), ("R -3 m", "R", "-3m"), ("P -3 1 m", "P", "-31m"),
    ("P 6/m", "P", "6/m"), ("P 6/m m m", "P", "6/mmm"),
    ("P m -3", "P", "m-3"), ("I m -3", "I", "m-3"), ("F m -3", "F", "m-3"), ("P m -3 m", "P", "m-3m"),
    ("I m -3 m", "I", "m-3m"), ("F m -3 m", "F", "m-3m")
)

ACCESIBLE_PATTERSON_NAME_HM = frozenset([_[0] for _ in REFERENCE_TABLE_PATTERSON_NAME_HM_LATTICE_TYPE_LAUE_CLASS])


def get_patterson_name_hm_by_lattice_type_laue_class(lattice_type: str, laue_class: str) -> str:
    patterson_name_hm = None
    for _1, _2, _3 in REFERENCE_TABLE_PATTERSON_NAME_HM_LATTICE_TYPE_LAUE_CLASS:
        if ((_2 == lattice_type) & (_3 == laue_class)):
            patterson_name_hm = _1
            break
    return patterson_name_hm


REFERENCE_TABLE_BRAVAIS_TYPE_CENTRING_TYPE_CRYSTAL_SYSTEM = (
    ("aP", "P", "triclinic"),
    ("mP", "P", "monoclinic"),
    ("mS", "A", "monoclinic"),
    ("mS", "B", "monoclinic"),
    ("mS", "C", "monoclinic"),
    ("mI", "I", "monoclinic"),
    ("oP", "P", "orthorhombic"),
    ("oS", "A", "orthorhombic"),
    ("oS", "B", "orthorhombic"),
    ("oS", "C", "orthorhombic"),
    ("oI", "I", "orthorhombic"),
    ("oF", "F", "orthorhombic"),
    ("tP", "P", "tetragonal"),
    ("tI", "I", "tetragonal"),
    ("hP", "P", "hexagonal"),
    ("hP", "P", "trigonal"),  # FIXME: not sure
    ("hR", "R", "trigonal"),
    ("hR", "Rrev", "trigonal"),
    ("hR", "H", "trigonal"),
    ("cP", "P", "cubic"),
    ("cI", "I", "cubic"),
    ("cF", "F", "cubic")
)

ACCESIBLE_BRAVAIS_TYPE = frozenset([_[0] for _ in REFERENCE_TABLE_BRAVAIS_TYPE_CENTRING_TYPE_CRYSTAL_SYSTEM])
ACCESIBLE_CRYSTAL_SYSTEM = frozenset([_[2] for _ in REFERENCE_TABLE_BRAVAIS_TYPE_CENTRING_TYPE_CRYSTAL_SYSTEM])


def get_bravais_type_by_centring_type_crystal_system(centring_type: str, crystal_system: str) -> str:
    bravais_type = None
    for _bravais_type, _centring_type, _crystal_system in REFERENCE_TABLE_BRAVAIS_TYPE_CENTRING_TYPE_CRYSTAL_SYSTEM:
        if ((_centring_type == centring_type) & (_crystal_system == crystal_system)):
            bravais_type = _bravais_type
    return bravais_type


def get_bravais_types_by_crystal_system(crystal_system: str) -> str:
    if crystal_system.startswith("tric"):
        bravais_types = ("aP",)
    elif crystal_system.startswith("m"):
        bravais_types = ("mP", "mS")
    elif crystal_system.startswith("o"):
        bravais_types = ("oP", "oS", "oI", "oF")
    elif crystal_system.startswith("te"):
        bravais_types = ("tI", "tF")
    elif crystal_system.startswith("h"):
        bravais_types = ("hP",)
    elif crystal_system.startswith("trig"):
        bravais_types = ("hR",)
    elif crystal_system.startswith("c"):
        bravais_types = ("cP", "cI", "cF")
    else:
        bravais_types = ()
    return bravais_types


def get_centring_types_by_bravais_type(_bravais_type: str) -> str:
    if _bravais_type.endswith("P"):
        centring_types = ("P",)
    elif _bravais_type.endswith("I"):
        centring_types = ("I",)
    elif _bravais_type.endswith("F"):
        centring_types = ("F",)
    elif _bravais_type.endswith("S"):
        centring_types = ("A", "B", "C")
    elif _bravais_type.endswith("R"):
        centring_types = ("R", "Rrev", "H",)
    return centring_types


def get_crystal_system_by_bravais_type(_bravais_type: str) -> str:
    crystal_system = None
    if _bravais_type.startswith("a"):
        crystal_system = "triclinic"
    elif _bravais_type.startswith("m"):
        crystal_system = "monoclinic"
    elif _bravais_type.startswith("o"):
        crystal_system = "orthorhombic"
    elif _bravais_type.startswith("t"):
        crystal_system = "tetragonal"
    elif _bravais_type == "hP":
        crystal_system = "hexagonal"
    elif _bravais_type == "hR":
        crystal_system = "trigonal"
    elif _bravais_type.startswith("c"):
        crystal_system = "cubic"
    return crystal_system


def get_type_hm(_name: str) -> str:
    l_res = []
    if _name in ACCESIBLE_NAME_HM_SHORT:
        l_res.append("short")
    if _name in ACCESIBLE_NAME_HM_FULL:
        l_res.append("full")
    if _name in ACCESIBLE_NAME_HM_EXTENDED:
        l_res.append("extended")
    return tuple(l_res)


def get_notation(_name: str) -> str:
    res = None
    if len(get_type_hm(_name)) != 0: res = "Hermann-Mauguin"
    if _name in ACCESIBLE_NAME_HALL_SHORT: res = "Hall"
    if _name in ACCESIBLE_NAME_SCHOENFLIES: res = "Schoenflies"
    return res


# IT A: Table 8.3.5.1. Sequence of generators for the crystal classes
# The space-group generators differ from those listed here by their glide or screw
# components. The generator 1 is omitted, except for crystal class 1. The
# subscript of a symbol denotes the characteristic direction of that operation,
# where necessary. The subscripts z, y, 110, 1-10, 10-1 and 111 refer to the
# directions [001], [010], [110], [1-10], [10-1] and [111], respectively. For mirror
# reflections m, the ‘direction of m’ refers to the normal to the mirror plane. The
# subscripts may be likewise interpreted as Miller indices of that plane
# Hermann–Mauguin symbol of crystal class
# Generators Gi (sequence left to right)
REFERENCE_TABLE_POINT_GROUP_HM_SYMBOL_GENERATORS = (
    ("1", ("1",)), ("-1", ("-1",)), ("2", ("2",)), ("m", ("m",)), ("2/m", (2, -1)), ("222", ("2z", "2y")),
    ("mm2", ("2z", "my")),
    ("mmm", ("2z", "2y", "-1")), ("4", ("2z", "4")), ("-4", ("2z", "-4")), ("4/m", ("2z", "4", "-1")),
    ("422", ("2z", "4", "2y")),
    ("4mm", ("2z", "4", "my")), ("-42m", ("2z", "-4", "2y")), ("-4m2", ("2z", "-4", "my")),
    ("4/mmm", ("2z", "4", "2y", "-1")), ("3", ("3",)),
    ("-3", ("3", "-1")), ("321", ("3", "2110")), ("321:r", ("3111", "210-1")), ("312", ("3", "21-10")),
    ("3m1", ("3", "m110")),
    ("3m1:r", ("3111", "m10-1")),
    ("31m", ("3", "m1-10")),
    ("-3m1", ("3", "2110", "-1")),
    ("-3m1:r", ("3111", "210-1", -1)),
    ("-31m", ("3", "21-10", "-1")),
    ("6", ("3", "2z")),
    ("-6", ("3", "mz")),
    ("6-m", ("3", "2z", "-1")),
    ("622", ("3", "2z", "2110")),
    ("6mm", ("3", "2z", "m110")),
    ("-6m2", ("3", "mz", "m110")),
    ("-62m", ("3", "mz", "2110")),
    ("6/mmm", ("3", "2z", "2110", "-1")),
    ("23", ("2z", "2y", "3111")),
    ("m-3", ("2z", "2y", "3111", "-1")),
    ("432", ("2z", "2y", "3111", "2110")),
    ("-43m", ("2z", "2y", "3111", "m1-10")),
    ("m-3m", ("2z", "2y", "3111", "2110", "-1"))
)


def get_generators_by_point_group_hm(name: str) -> Tuple[str]:
    generators = ()
    for _1, _2 in REFERENCE_TABLE_POINT_GROUP_HM_SYMBOL_GENERATORS:
        if _1 == name:
            generators = _2
            break
    return generators


# IT A: Table 10.1.2.4. Names and symbols of the 32 crystal classes
REFERENCE_TABLE_LAUE_CLASS_SHORT_FULL_POINT_GROUP_HM_SYMBOL_SHORT_FULL_SCHOENFLIES = (
    ("-1", "-1", "1", "1", "C1"),
    ("-1", "-1", "-1", "-1", "Ci"),
    ("2/m", "2/m", "2", "2", "C2"),
    ("2/m", "2/m", "m", "m", "Cs"),
    ("2/m", "2/m", "2/m", "2/m", "C2h"),
    ("mmm", "2/m2/m2/m", "222", "222", "D2"),
    ("mmm", "2/m2/m2/m", "mm2", "mm2", "C2v"),
    ("mmm", "2/m2/m2/m", "mmm", "2/m2/m2/m", "D2h"),
    ("4/m", "4/m", "4", "4", "C4"),
    ("4/m", "4/m", "-4", "-4", "S4"),
    ("4/m", "4/m", "4/m", "4/m", "C4h"),
    ("4/mmm", "4/m2/m2/m", "422", "422", "D4"),
    ("4/mmm", "4/m2/m2/m", "4mm", "4mm", "C4v"),
    ("4/mmm", "4/m2/m2/m", "-42m", "-42m", "D2d"),
    ("4/mmm", "4/m2/m2/m", "4/mmm", "4/m2/m2/m", "D4h"),
    ("-3", "-3", "3", "3", "C3"),
    ("-3", "-3", "-3", "-3", "C3i"),
    ("-3m", "-32/m", "32", "32", "D3"),
    ("-3m", "-32/m", "3m", "3m", "C3v"),
    ("-3m", "-32/m", "-3m", "-32/m", "D3d"),
    ("6/m", "6/m", "6", "6", "C6"),
    ("6/m", "6/m", "-6", "-6", "C3h"),
    ("6/m", "6/m", "6/m", "6/m", "C6h"),
    ("6/mmm", "6/m2/m2/m", "622", "622", "D6"),
    ("6/mmm", "6/m2/m2/m", "6mm", "6mm", "D6v"),
    ("6/mmm", "6/m2/m2/m", "-62m", "-62m", "D3h"),
    ("6/mmm", "6/m2/m2/m", "6/mmm", "6/m2/m2/m", "D6h"),
    ("m-3", "2/m-3", "23", "23", "T"),
    ("m-3", "2/m-3", "m-3", "2/m-3", "Th"),
    ("m-3m", "4/m-32/m", "432", "432", "O"),
    ("m-3m", "4/m-32/m", "-43m", "-43m", "Td"),
    ("m-3m", "4/m-32/m", "m-3m", "4/m-32/m", "Oh")
)

ACCESIBLE_LAUE_CLASS_FULL = frozenset(
    [_[1] for _ in REFERENCE_TABLE_LAUE_CLASS_SHORT_FULL_POINT_GROUP_HM_SYMBOL_SHORT_FULL_SCHOENFLIES])
ACCESIBLE_POINT_GROUP_SYMBOL_SHORT = frozenset(
    [_[2] for _ in REFERENCE_TABLE_LAUE_CLASS_SHORT_FULL_POINT_GROUP_HM_SYMBOL_SHORT_FULL_SCHOENFLIES])
ACCESIBLE_POINT_GROUP_SYMBOL_FULL = frozenset(
    [_[3] for _ in REFERENCE_TABLE_LAUE_CLASS_SHORT_FULL_POINT_GROUP_HM_SYMBOL_SHORT_FULL_SCHOENFLIES])


def get_laue_class_by_name_schoenflies(name: str) -> str:
    laue_class = None
    symb = name.split(".")[0]
    for _1, _2, _3, _4, _5 in REFERENCE_TABLE_LAUE_CLASS_SHORT_FULL_POINT_GROUP_HM_SYMBOL_SHORT_FULL_SCHOENFLIES:
        if _5 == symb:
            laue_class = _1
            break
    return laue_class


def get_point_group_hm_full_by_name_schoenflies(name: str) -> str:
    point_group_hm_full = None
    symb = name.split(".")[0]
    for _1, _2, _3, _4, _5 in REFERENCE_TABLE_LAUE_CLASS_SHORT_FULL_POINT_GROUP_HM_SYMBOL_SHORT_FULL_SCHOENFLIES:
        if _5 == symb:
            point_group_hm_full = _4
            break
    return point_group_hm_full


def get_point_group_hm_short_by_name_schoenflies(name: str) -> str:
    point_group_hm_short = None
    symb = name.split(".")[0]
    for _1, _2, _3, _4, _5 in REFERENCE_TABLE_LAUE_CLASS_SHORT_FULL_POINT_GROUP_HM_SYMBOL_SHORT_FULL_SCHOENFLIES:
        if _5 == symb:
            point_group_hm_short = _3
            break
    return point_group_hm_short


def get_name_schoenfliess_by_laue_class(laue_class: str) -> Tuple[str]:
    l_symb = [_5 for _1, _2, _3, _4, _5 in
              REFERENCE_TABLE_LAUE_CLASS_SHORT_FULL_POINT_GROUP_HM_SYMBOL_SHORT_FULL_SCHOENFLIES if _1 == laue_class]
    l_res = []
    for symb in l_symb:
        for _name_schoenflies in ACCESIBLE_NAME_SCHOENFLIES:
            _symb = _name_schoenflies.split(".")[0]
            if symb == _symb:
                l_res.append(_name_schoenflies)
    return tuple(l_res)


def get_name_schoenfliess_by_point_group_hm_short(point_group: str) -> Tuple[str]:
    l_symb = [_5 for _1, _2, _3, _4, _5 in
              REFERENCE_TABLE_LAUE_CLASS_SHORT_FULL_POINT_GROUP_HM_SYMBOL_SHORT_FULL_SCHOENFLIES if _3 == point_group]
    l_res = []
    for symb in l_symb:
        for _name_schoenflies in ACCESIBLE_NAME_SCHOENFLIES:
            _symb = _name_schoenflies.split(".")[0]
            if symb == _symb:
                l_res.append(_name_schoenflies)
    return tuple(l_res)


def get_centrosymmetry_by_name_hall(name: str) -> str:
    centrosymmetry = name.startswith("-")
    if not (name in ACCESIBLE_NAME_HALL_SHORT):
        centrosymmetry = None
    return centrosymmetry


def separate_notation_it_coordinate_system_code(name: str):
    l_h = name.strip().split(":")
    notation = l_h[0].strip()
    if notation.isdigit():
        notation = int(notation)
    if len(l_h) == 1:
        it_coordinate_system_code = None
    else:
        it_coordinate_system_code = l_h[1].strip()
        if not (it_coordinate_system_code in ACCESIBLE_IT_COORDINATE_SYSTEM_CODE):
            it_coordinate_system_code = None
    return notation, it_coordinate_system_code


def get_symop_pcentr_multiplicity_letter_site_symmetry_coords_xyz_2(
        it_number: int, it_coordinate_system_code: str):
    """
    FIXME: HOW it works for 166 space group
    crystal system should be trigonal or hexagonal
    """
    crystal_system = get_crystal_system_by_it_number(it_number)
    if it_coordinate_system_code is None:
        choice = "1"
    elif "3" in it_coordinate_system_code:
        choice = "3"
    elif "2" in it_coordinate_system_code:
        choice = "2"
    elif "1" in it_coordinate_system_code:
        choice = "1"
    elif "h" in it_coordinate_system_code:
        # FIXME: IT SHOULD BE CHECKED
        # if crystal_system.startswith("trigonal"):
        #     choice = "2"
        # else:  # hexagonal
        choice = "1"
    elif "r" in it_coordinate_system_code:
        choice = "2"
    else:
        choice = "1"
    symop, p_centr = None, None
    for _el_card in EL_CARDS:
        if ((_el_card["it_number"] == it_number) & (_el_card["choice"][0] == choice)):
            symop = tuple(_el_card["symmetry"])
            # FIXME: It looks that definition of p_centr is inversed one compare with parameters given in card (at least for 227 choise 1)
            # It should be checked for all space group. It is why here factor -1.
            p_centr = -1*array([Fraction(_).limit_denominator(10) for _ in _el_card["pcentr"][0].split(",")],
                             dtype=Fraction)
            break
    _s_name, _choice = get_transform_pp_abc_choice_by_it_number_it_coordinate_system_code(it_number,
                                                                                          it_coordinate_system_code)
    Q, p = transform_string_to_r_b(_s_name, ("a", "b", "c"))
    P = transpose(Q)
    q = -1 * mult_matrix_vector(Q, p)
    p_centr_new = p_centr + q
    symop_2 = [transform_symop_operation_xyz_by_pp_abc(_symop, P, p) for _symop in symop]

    for _el_card in WYCKOFF:
        if ((_el_card["it_number"] == it_number) & (_el_card["choice"] == int(choice))):
            wyckoff = _el_card["wyckoff"]
            break
    l_multiplicity = [_h["multiplicity"] for _h in wyckoff]
    l_letter = [_h["letter"] for _h in wyckoff]
    l_site_symmetry = [_h["site_symmetry"] for _h in wyckoff]
    l_coord_xyz = [_h["symop"] for _h in wyckoff]
    l_coord_xyz_2 = [[transform_symop_operation_xyz_by_pp_abc(_coord_xyz, P, p) for _coord_xyz in coord_xyz] for
                     coord_xyz in l_coord_xyz]

    return symop_2, p_centr_new, l_multiplicity, l_letter, l_site_symmetry, l_coord_xyz_2


def transform_symop_operation_xyz_by_pp_abc(symop_operation_xyz: str, P, p) -> str:
    Q = transpose(P)  # TODO: here is proposed that Q^T = Q**-1, but I am not sure that it is true.
    q = -1 * mult_matrix_vector(Q, p)

    r_xyz, b_xyz = transform_string_to_r_b(symop_operation_xyz, ("x", "y", "z"))
    b_new = zeros(shape=(3), dtype=float)
    r_new = zeros(shape=(3, 3), dtype=float)
    QW = mult_matrixes(Q, r_xyz)
    QWP = mult_matrixes(QW, P)
    QWp = mult_matrix_vector(QW, p)
    Qw = mult_matrix_vector(Q, b_xyz)
    r_new = QWP
    b_new = QWp + Qw + q
    _s = transform_r_b_to_string(r_new, b_new, ("x", "y", "z"))
    return _s


def transform_symop_operation_xyz_by_Qq_xyz(symop_operation_xyz: str, Q, q) -> str:
    P = transpose(q)  # TODO: here is proposed that Q^T = Q**-1, but I am not sure that it is true.
    p = -1 * mult_matrix_vector(P, q)
    _s = transform_symop_operation_xyz_by_pp_abc(symop_operation_xyz, P, p)
    return _s


def mult_matrix_vector(a, v):
    cond_1 = isinstance(v[0], Fraction)
    cond_2 = isinstance(a[0, 0], Fraction)
    if (cond_1 & cond_2):
        p_0 = a[0, 0]*v[0] + a[0, 1]*v[1] + a[0, 2]*v[2]
        p_1 = a[1, 0]*v[0] + a[1, 1]*v[1] + a[1, 2]*v[2]
        p_2 = a[2, 0]*v[0] + a[2, 1]*v[1] + a[2, 2]*v[2]
        b = array([p_0, p_1, p_2], dtype=Fraction)
    else:
        p_0 = float(a[0, 0])*float(v[0]) + float(a[0, 1])*float(v[1]) + float(a[0, 2])*float(v[2])
        p_1 = float(a[1, 0])*float(v[0]) + float(a[1, 1])*float(v[1]) + float(a[1, 2])*float(v[2])
        p_2 = float(a[2, 0])*float(v[0]) + float(a[2, 1])*float(v[1]) + float(a[2, 2])*float(v[2])
        b = array([p_0, p_1, p_2], dtype=float)
    return b

 
def mult_matrixes(a, b):
    c = 0. * a
    for _i in range(3):
        for _j in range(3):
            c[_i, _j] = sum(a[_i, :] * b[:, _j])
    return c



def auto_choose_it_coordinate_system_code(it_number:int, it_coordinate_system_codes:list)->str:
    if len(it_coordinate_system_codes) == 0:
        it_coordinate_system_code = None
    elif len(it_coordinate_system_codes) > 1:
        print(f"Several values of it_coordinate_system_code have been defined:")
        print_long_list(it_coordinate_system_codes)
        default_i_c_s_c = get_default_it_coordinate_system_code_by_it_number(it_number)
        if default_i_c_s_c in it_coordinate_system_codes:
            it_coordinate_system_code = default_i_c_s_c
            print(f"The default value has been choosen:'{it_coordinate_system_code:}'.")
        else:
            l_1 = [_ for _ in it_coordinate_system_codes if not ("-" in _)]
            if len(l_1) != 0:
                _choice = l_1[0]
            else:
                _choice = it_coordinate_system_codes[0]
            it_coordinate_system_code = _choice
            print(f"The \"{it_coordinate_system_code:}\" has been choosen.")
    else:
        it_coordinate_system_code = it_coordinate_system_codes[0]
    return it_coordinate_system_code


def get_transform_pp_abc_choice_by_it_number_it_coordinate_system_code(it_number: int,
                                                                       it_coordinate_system_code: str) -> Tuple:
    # TODO: not sure about -b1, c1, -c1, a1, -a1
    if it_coordinate_system_code in ("b1", "b2", "b3", "abc", "1abc", "2abc", "1", "2", "h", "r", None):
        transform_pp_abc = "a,b,c"
    elif it_coordinate_system_code in ("-a1", "-a2", "-a3", "ba-c", "1ba-c", "2ba-c"):
        transform_pp_abc = "b,a,-c"
    elif it_coordinate_system_code in ("c1", "c2", "c3", "cab", "1cab", "2cab"):
        transform_pp_abc = "c,a,b"
    elif it_coordinate_system_code in ("-b1", "-b2", "-b3", "-cba", "1-cba", "2-cba"):
        transform_pp_abc = "-c,b,a"
    elif it_coordinate_system_code in ("a1", "a2", "a3", "bca", "1bca", "2bca"):
        transform_pp_abc = "b,c,a"
    elif it_coordinate_system_code in ("-c1", "-c2", "-c3", "a-cb", "1bca", "2a-cb"):
        transform_pp_abc = "a,-c,b"

    if it_coordinate_system_code is None:
        choice = 1
    elif "2" in it_coordinate_system_code:
        choice = 2
    elif "h" in it_coordinate_system_code:
        crystal_system = get_crystal_system_by_it_number(it_number)
        if crystal_system.startswith("trigonal"):
            choice = 2
        else:  # hexagonal
            choice = 1
    elif "3" in it_coordinate_system_code:
        choice = 3
    else:
        choice = 1
    return transform_pp_abc, choice


def print_long_list(ll):
    ls_out, s_line = [], []
    max_size = max([len(str(_)) for _ in ll])
    length_size = 80
    number_per_line = int(length_size // max_size)
    _i = Fraction(1, number_per_line)
    for _ in ll:
        s_line.append(str(_).rjust(max_size))
        if _i.denominator == 1:
            ls_out.append(", ".join(s_line))
            s_line = []
        _i += Fraction(1, number_per_line)
    ls_out.append(", ".join(s_line))
    print("\n".join(ls_out).rstrip())
    return


def devide(l_a, b, dev):
    if dev is not None:
        l_a_o = [a / dev for a in l_a]
        b_o = b / dev
    else:
        l_a_o = [a for a in l_a]
        b_o = b
    return l_a_o, b_o


def one_line(l_a, b, l_ind_exclude):
    l_ind_non_zeros = [i_a for i_a, a in enumerate(l_a) if (not (i_a in l_ind_exclude) & (a != 0))]
    l_ind_non_zeros = []
    for i_a, a in enumerate(l_a):
        flag_1 = not (i_a in l_ind_exclude)
        flag_2 = (a != 0)
        if (flag_1 & flag_2):
            l_ind_non_zeros.append(i_a)
    ind_1, dev_1 = None, None
    if len(l_ind_non_zeros) != 0:
        ind_1 = l_ind_non_zeros[0]
        dev_1 = l_a[ind_1]
    l_a_o, b_o = devide(l_a, b, dev_1)
    return l_a_o, b_o, ind_1, dev_1


def is_solution_a_b(ll_a, l_b):
    if all([b == 0 for b in l_b]):
        return True
    l_ind_exclude = []

    l_a_in_1, b_in_1 = ll_a[0], (l_b[0])%1
    l_a_in_2, b_in_2 = ll_a[1], (l_b[1])%1
    l_a_in_3, b_in_3 = ll_a[2], (l_b[2])%1
    l_a_1, b_1, ind_1, dev_1 = one_line(l_a_in_1, b_in_1, l_ind_exclude)
    if ind_1 is not None:
        val_2 = l_a_in_2[ind_1]
        l_a_in_2 = [_1 - val_2 * _2 for _1, _2 in zip(l_a_in_2, l_a_1)]
        b_in_2 = (b_in_2 - val_2 * b_1) % 1
        val_3 = l_a_in_3[ind_1]
        l_a_in_3 = [_1 - val_3 * _2 for _1, _2 in zip(l_a_in_3, l_a_1)]
        b_in_3 = (b_in_3 - val_3 * b_1) % 1
        l_ind_exclude.append(ind_1)
    elif b_in_1 != 0:
        return False

    l_a_2, b_2, ind_2, dev_2 = one_line(l_a_in_2, b_in_2, l_ind_exclude)
    if ind_2 is not None:
        val_3 = l_a_in_3[ind_2]
        l_a_in_3 = [_1 - val_3 * _2 for _1, _2 in zip(l_a_in_3, l_a_2)]
        b_in_3 = (b_in_3 - val_3 * b_2) % 1
        l_ind_exclude.append(ind_2)
    elif b_in_2 != 0:
        return False

    l_a_3, b_3, ind_3, dev_3 = one_line(l_a_in_3, b_in_3, l_ind_exclude)
    if ind_3 is not None:
        l_ind_exclude.append(ind_3)
    elif b_in_3 != 0:
        return False
    return True


def is_good_for_mask(r, b, fract_x, fract_y, fract_z):
    b_1 = array([(fract_x - b[0]) % 1, (fract_y - b[1]) % 1, (fract_z - b[2]) % 1], dtype=Fraction)
    flag_1 = is_solution_a_b(r, b_1)
    return flag_1

