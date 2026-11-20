"""CIF/STAR categories that CrysPy does not know get a class of their own
(make_item_class, make_loop_class): described by the ATTR_* tuples like
every item class, so their attributes are descriptors and the objects are
copied and pickled as the others."""
import copy
import pickle

import cryspy
from cryspy.B_parent_classes.attributes import take_item_classes
from cryspy.B_parent_classes.cif_attributes import CifValue
from cryspy.B_parent_classes.cl_1_item import ItemN, is_generic, make_item_class
from cryspy.B_parent_classes.cl_2_loop import LoopN, make_loop_class

LOOP = """loop_
_publ_author_name
_publ_author_address
'Tuilier, M H' 'Bordeaux'
'Chevalier, B' ?
"""

ITEM = """_journal_name_full 'Physica C'
_journal_volume 200
_journal_page_first 113(2)
"""


def test_unknown_loop_has_its_class():
    loop = LoopN.from_cif(LOOP)
    assert is_generic(loop) and isinstance(loop, LoopN) and type(loop) is not LoopN
    item_class = type(loop).ITEM_CLASS
    assert issubclass(item_class, ItemN) and item_class.PREFIX == "publ_author"
    assert item_class.ATTR_NAMES == ("name", "address")
    assert isinstance(item_class.__dict__["name"], CifValue)
    assert loop.items[0].name == "Tuilier, M H" and loop.items[1].address is None
    assert loop.name == ["Tuilier, M H", "Chevalier, B"]
    assert loop.get_name() == "publ_author"
    # one class for one category
    assert type(LoopN.from_cif(LOOP)) is type(loop)
    assert make_loop_class(item_class) is type(loop)


def test_unknown_item_has_its_class():
    item = ItemN.from_cif(ITEM)
    assert is_generic(item) and type(item) is not ItemN
    assert type(item).PREFIX == "journal"
    # numbers are refinable values, integers stay integers, the other
    # values text
    assert type(item).ATTR_REF == ("page_first", )
    assert item.volume == 200 and isinstance(item.volume, int)
    assert item.page_first == 113.
    assert item.page_first_sigma == 2. and item.page_first_refinement
    assert item.name_full == "Physica C"
    assert "_journal_page_first 113(2)" in item.to_cif()
    assert "_journal_volume 200\n" in item.to_cif()


def test_types_of_unknown_loop_columns():
    """A column of numbers is refinable as in the categories of CrysPy
    (a copy of atom_site_susceptibility in an experiment block), a column
    of integers is written back as integers, '007' is text."""
    loop = LoopN.from_cif("""loop_
_dummy_parameter_label
_dummy_parameter_chi
_dummy_parameter_order
_dummy_parameter_code
A 3.329(48) 1 007
B 3.32908 2 008
C . 3 ?
""")
    item_class = type(loop).ITEM_CLASS
    assert item_class.ATTR_REF == ("chi", )
    assert item_class.ATTR_TYPES == (str, float, int, str)
    assert loop.items[0].chi_refinement and loop.items[0].chi_sigma == 0.048
    assert loop.items[2].chi is None and loop.items[2].code is None
    assert loop.items[1].code == "008"
    text = loop.to_cif()
    assert "A 3.329(48) 1 007" in " ".join(text.split()).replace("  ", " ")


def test_copy_and_pickle_of_unknown_categories():
    for obj in (LoopN.from_cif(LOOP), ItemN.from_cif(ITEM)):
        for copied in (copy.deepcopy(obj), pickle.loads(pickle.dumps(obj))):
            assert type(copied) is type(obj)
            assert copied.to_cif() == obj.to_cif()
    item_class = make_item_class("my_category", ("a", "b"), attr_ref=("a", ))
    assert pickle.loads(pickle.dumps(item_class)) is item_class
    loop_class = make_loop_class(item_class)
    assert pickle.loads(pickle.dumps(loop_class)) is loop_class


def test_unknown_categories_are_not_classes_of_the_library():
    ItemN.from_cif(ITEM)
    assert all(not is_generic(item_class) for item_class in take_item_classes())
    assert not is_generic(cryspy.Cell) and not is_generic(cryspy.Cell())


def test_define_by_items():
    item_class = make_item_class("my_category", ("a", "b"))
    l_item = [item_class(), item_class()]
    l_item[0].a, l_item[1].a = "x", "y"
    loop = LoopN.define_by_items(l_item, loop_name="first")
    assert type(loop) is make_loop_class(item_class) and loop.loop_name == "first"
    assert loop.a == ["x", "y"]
