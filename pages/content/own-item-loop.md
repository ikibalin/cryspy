---
layout: page
permalink: /content/own-item-loop/
---
[Back to Content][content]

**The template of python script for new item and new loop objects**

{% highlight python %}
from cryspy import ItemN, LoopN

class MyItem(ItemN):
    """MyItem class.

    Attributes
    ----------
        - param_1, param_2 (mandatory) 
        - param_3 (optional) 
    """
    # mandatory attributes
    ATTR_MANDATORY_NAMES = ("param_1", "param_2")
    # types of mandatory attributes: int, float, str, complex
    ATTR_MANDATORY_TYPES = (int, float)
    # definition of mandatory attributes in cif
    ATTR_MANDATORY_CIF = ("param_1", "param_2")

    # optional attributes
    ATTR_OPTIONAL_NAMES = ("param_3", )
    ATTR_OPTIONAL_TYPES = (str, )
    ATTR_OPTIONAL_CIF = ("param_3", )

    ATTR_NAMES = ATTR_MANDATORY_NAMES + ATTR_OPTIONAL_NAMES
    ATTR_TYPES = ATTR_MANDATORY_TYPES + ATTR_OPTIONAL_TYPES
    ATTR_CIF = ATTR_MANDATORY_CIF + ATTR_OPTIONAL_CIF

    ATTR_INT_NAMES = ()
    ATTR_INT_PROTECTED_NAMES = ()

    # parameters considered as refined parameters
    ATTR_REF = ("param_2", )
    ATTR_SIGMA = tuple([f"{_h:}_sigma" for _h in ATTR_REF])
    ATTR_CONSTR_FLAG = tuple([f"{_h:}_constraint" for _h in ATTR_REF])
    ATTR_REF_FLAG = tuple([f"{_h:}_refinement" for _h in ATTR_REF])

    # formats if cif format
    D_FORMATS = {'param_1': "{:4}", 'param_2': "{:.4f}"}

    # constraints on the parameters
    D_CONSTRAINTS = {}

    # default values for the parameters
    D_DEFAULT = {"param_1": 1}
    for key in ATTR_SIGMA:
        D_DEFAULT[key] = 0.
    for key in (ATTR_CONSTR_FLAG + ATTR_REF_FLAG):
        D_DEFAULT[key] = False

    PREFIX = "my_item" # prefix for cif file

    def __init__(self, **kwargs) -> NoReturn:
        super(MyItem, self).__init__()

        # defined for any integer and float parameters
        D_MIN = {"param_1": 0, "param_2": 0.}

        # defined for ani integer and float parameters
        D_MAX = {"param_1": 5}

        self.__dict__["D_MIN"] = D_MIN
        self.__dict__["D_MAX"] = D_MAX
        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)

	def my_method():
		pass


# loop class
class MyObjectL(LoopN):
    """Description of MyObject in loop."""

    ITEM_CLASS = MyObject
    ATTR_INDEX = None

    def __init__(self, loop_name=None) -> NoReturn:
        super(MyObjectL, self).__init__()
        self.__dict__["items"] = []
        self.__dict__["loop_name"] = loop_name
	
	def my_method():
		pass
{% endhighlight %}


[content]: /cryspy/content
