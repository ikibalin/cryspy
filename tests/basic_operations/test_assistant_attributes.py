"""Attributes of any ItemN/LoopN class: allowed values, refinable tags,
generic access of the assistant tools (also for classes defined here)."""
import cryspy
from cryspy.B_parent_classes.attributes import take_item_classes
from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN
from cryspy.B_parent_classes.cl_3_data import DataN
from cryspy.B_parent_classes.cl_4_global import GlobalN


class NewMode(ItemN):
    """New category defined outside CrysPy."""
    ATTR_MANDATORY_NAMES = ("label", )
    ATTR_MANDATORY_TYPES = (str, )
    ATTR_MANDATORY_CIF = ("label", )
    ATTR_OPTIONAL_NAMES = ("kind", "level", "width")
    ATTR_OPTIONAL_TYPES = (str, int, float)
    ATTR_OPTIONAL_CIF = ("kind", "level", "width")
    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF
    ATTR_INT_NAMES = ()
    ATTR_INT_PROTECTED_NAMES = ()
    ATTR_REF = ("width", )
    ATTR_SIGMA = ("width_sigma", )
    ATTR_CONSTR_FLAG = ("width_constraint", )
    ATTR_REF_FLAG = ("width_refinement", )
    ATTR_CONSTR_MARK = ("width_mark", )
    D_CONSTRAINTS = {"kind": ["fast", "slow", "off"], "level": {1, 2, 3}}
    D_DEFAULT = {"width_sigma": 0., "width_constraint": False,
                 "width_refinement": False, "width_mark": ""}
    PREFIX = "new_mode"

    def __init__(self, **kwargs):
        super(NewMode, self).__init__()
        self.__dict__["D_MIN"] = {}
        self.__dict__["D_MAX"] = {}
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)


class NewModeL(LoopN):
    ITEM_CLASS = NewMode
    ATTR_INDEX = "label"

    def __init__(self, loop_name: str = None, **kwargs):
        super(NewModeL, self).__init__()
        self.__dict__["items"] = form_items(kwargs)
        self.__dict__["loop_name"] = loop_name


def form_items(kwargs):
    return [NewMode(**d_item) for d_item in kwargs.get("items", [])]


def take_new_object():
    data = DataN.make_container((), (NewModeL, ), "new_block")
    data.data_name = "blk"
    loop = NewModeL()
    loop.items.extend([NewMode(label="A", kind="fast", level=1, width=0.5),
                       NewMode(label="B", kind="off", level=2, width=0.7)])
    data.add_items([loop])
    obj = GlobalN.make_container((), (DataN, ), "global")
    obj.add_items([data])
    return obj


def test_value_choices_from_class():
    assert cryspy.take_value_choices(cryspy.AtomSiteSusceptibilityL) == {
        "_atom_site_susceptibility_chi_type": ("Ciso", "Cani")}
    assert cryspy.take_value_choices(cryspy.Extinction()) == {
        "_extinction_model": ("gauss", "lorentz", "qani")}
    # a set or a long list (space group names) is not a choice
    d_choice = cryspy.take_value_choices(cryspy.SpaceGroup)
    assert "_space_group_name_h-m_alt" not in d_choice
    assert "_space_group_it_number" not in d_choice
    assert cryspy.take_value_choices(NewModeL) == {"_new_mode_kind": ("fast", "slow", "off")}
    assert cryspy.take_value_choices(LoopN()) == {}


def test_next_value_choice():
    choices = ("Ciso", "Cani")
    assert cryspy.next_value_choice(choices, "Ciso") == "Cani"
    assert cryspy.next_value_choice(choices, "cani") == "Ciso"
    assert cryspy.next_value_choice(choices, "unknown") == "Ciso"


def test_refinable_tags_and_descriptions():
    assert "_cell_length_a" in cryspy.take_refinable_cif_tags(cryspy.Cell)
    assert cryspy.take_refinable_cif_tags(NewMode()) == {"_new_mode_width"}
    d_attribute = {d["name"]: d for d in cryspy.describe_attributes(NewModeL)}
    assert d_attribute["kind"]["allowed_values"] == ("fast", "slow", "off")
    assert d_attribute["level"]["allowed_values"] == (1, 2, 3)
    assert d_attribute["label"]["mandatory"]
    assert d_attribute["width"]["refinable"]
    assert NewMode in take_item_classes()


def call(obj, tool, **arguments):
    return cryspy.call_assistant_tool(obj, tool, arguments)


def test_tools_reach_new_classes():
    obj = take_new_object()
    text = call(obj, "list_attributes", pattern="kind")["text"]
    assert "new_block_blk.new_mode[A].kind | fast | fast, slow, off" in text
    result = call(obj, "set_attribute", path="new_block_blk.new_mode[B].kind", value="SLOW")
    assert not result["flag_error"] and result["flag_object_changed"]
    assert obj.items[0].items[0].items[1].kind == "slow"
    result = call(obj, "set_attribute", path="new_mode[A].kind", value="medium")
    assert result["flag_error"] and "fast, slow, off" in result["text"]
    result = call(obj, "set_attribute", path="new_mode[A].level", value="3")
    assert obj.items[0].items[0].items[0].level == 3
    assert "blk.new_mode[A].width" in call(obj, "list_parameters")["text"]
    assert "B" in call(obj, "show_item", path="blk.new_mode")["text"]
    assert "B" in call(obj, "show_item", path="new_block_blk.new_mode")["text"]
    text = call(obj, "show_item", path="new_mode[A]")["text"]
    assert "_new_mode_kind fast" in text and "off" not in text
    text = call(None, "describe_category", name="new_mode")["text"]
    assert "_new_mode_kind (str): fast, slow, off" in text
    assert "_new_mode_width (float, refinable)" in text


def test_value_of_attribute_is_explained():
    obj = cryspy.str_to_globaln(
        "global_\n\ndata_fe\n_cell_length_a 5.0\n_cell_angle_alpha 90\n")
    text = call(obj, "list_parameters", pattern="lorentz")["text"]
    assert "extinction.model" in text and "set_attribute" in text
