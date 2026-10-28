"""Methods of the data blocks, loops and items of a CrysPy object.

Every public method of a class derived from DataN, LoopN or ItemN (as
Pd.define_background_points or PdBackgroundL.define_points) can be called
by an AI assistant or another program, without a list of methods anywhere:
the methods are found in the classes and described from their signatures
and their docstrings in the numpydoc format (procedures.parse_docstring),
as the procedures of procedures.py:

    >>> import cryspy
    >>> obj = cryspy.load_file("main.rcif")
    >>> cryspy.call_assistant_tool(obj, "list_methods", {"pattern": "background"})
    >>> cryspy.call_assistant_tool(obj, "call_method", {
    ...     "path": "pnd", "method": "define_background_points", "step_ttheta": 2})

A parameter of a method is a basic value (float, int, bool, str; the type
is taken from the type hint or from the default value) or an item, loop or
data block of the object given by its class in the type hint (as
pd_meas: PdMeasL); the latter is found in the same data block, then in the
whole object, and is asked from the user when there are several. Methods
of the parent classes (to_cif, add_items, ...), private, deprecated,
plotting and internal ones (get_*, form_*, take_*, ...) are left out, and
so are methods without a docstring (the docstring is their description)
and methods with parameters that cannot be given (arrays).
"""
import inspect
import re
from typing import Callable, Dict, Iterator, List, Optional, Tuple

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN
from cryspy.B_parent_classes.cl_3_data import DataN
from cryspy.B_parent_classes.cl_4_global import GlobalN
from cryspy.procedure_assistant.procedures import (
    D_BASIC_TYPE, form_question, parse_docstring, take_annotation_types)

BASE_CLASSES = (ItemN, LoopN, DataN, GlobalN)
# internal methods: the dictionary pipeline, plots, flags, checks
SKIP_PREFIXES = ("_", "get_", "set_", "take_", "form_", "plot", "view", "is_", "to_", "numpy_",
                 "apply_constraints", "copy")


def take_classes_with_methods() -> List[type]:
    """All classes derived from DataN, LoopN and ItemN (also the classes
    of the loops)."""
    l_class, l_todo = [], list(BASE_CLASSES[:3])
    while l_todo:
        for sub_class in l_todo.pop(0).__subclasses__():
            if sub_class not in l_class:
                l_class.append(sub_class)
                l_todo.append(sub_class)
    return l_class


def take_method_parameters(function: Callable) -> Optional[List[tuple]]:
    """(name, kind, types, default) of the parameters of a method without
    self; kind: "basic" (types: (float, ), (int, ), ...) or "cif" (types:
    classes of items, loops or data blocks; a base class as ItemN is
    found by the name of the parameter). None when the method cannot be
    called with such values."""
    try:
        signature = inspect.signature(function)
    except (TypeError, ValueError):
        return None
    l_out = []
    for ind, param in enumerate(signature.parameters.values()):
        if ind == 0 and param.name in ("self", "cls"):
            continue
        if param.kind in (param.VAR_POSITIONAL, param.VAR_KEYWORD):
            continue
        types = take_annotation_types(param.annotation)
        default = param.default
        flag_default = default is not inspect.Parameter.empty
        if param.annotation is inspect.Parameter.empty and flag_default and default is not None:
            types = (type(default), )
        if len(types) == 1 and types[0] in D_BASIC_TYPE:
            l_out.append((param.name, "basic", types, default))
        elif all(isinstance(item, type) and issubclass(item, (ItemN, LoopN, DataN))
                 for item in types):
            l_out.append((param.name, "cif", types, default))
        elif flag_default:
            continue  # given by its default
        else:
            return None
    return l_out


def is_method_for_assistant(name: str, function) -> bool:
    if not inspect.isfunction(function) or name.startswith(SKIP_PREFIXES):
        return False
    if getattr(function, "__deprecated__", False) or not (function.__doc__ or "").strip():
        return False  # deprecated, or without a description
    return take_method_parameters(function) is not None


def take_methods(cls: type) -> Dict[str, Callable]:
    """{name: function} of the methods of the class (and of its parents
    except DataN, LoopN, ItemN, GlobalN) that an assistant may call."""
    d_method = {}
    for klass in reversed(cls.__mro__):
        if klass in BASE_CLASSES or klass is object or not issubclass(klass, BASE_CLASSES):
            continue
        for name, function in vars(klass).items():
            if is_method_for_assistant(name, function):
                d_method[name] = function
            elif name in d_method:  # redefined as something else
                del d_method[name]
    return {name: d_method[name] for name in sorted(d_method)}


def describe_method(function: Callable) -> dict:
    """Description from the signature and the docstring: {"name",
    "summary", "description", "returns", "parameters": [{"name", "kind",
    "type" (JSON schema type of "basic"), "classes" (of "cif"),
    "default" (when there is one), "description", "question"}]}."""
    d_doc = parse_docstring(function.__doc__)
    l_param = []
    for param, kind, types, default in take_method_parameters(function):
        description = d_doc["parameters"].get(param, {}).get("description", "")
        d_param = {"name": param, "kind": kind, "description": description}
        if kind == "basic":
            d_param["type"] = D_BASIC_TYPE[types[0]]
        else:
            d_param["classes"] = [item.__name__ for item in types]
        d_param["question"] = form_question(description, f"the {param.replace('_', ' '):}")
        if default is not inspect.Parameter.empty:
            d_param["default"] = default
        l_param.append(d_param)
    summary = d_doc["summary"]
    return {"name": function.__name__, "summary": summary, "description": d_doc["description"],
            "returns": d_doc["returns"], "parameters": l_param}


def format_method(d_method: dict) -> str:
    """'define_points(pd_meas: PdMeasL, step_ttheta: number=10.0): summary'."""
    ls_param = []
    for d_param in d_method["parameters"]:
        s_type = d_param.get("type") or "|".join(d_param["classes"])
        s_default = f"={d_param['default']!r:}" if "default" in d_param else ""
        ls_param.append(f"{d_param['name']:}: {s_type:}{s_default:}")
    return f"{d_method['name']:}({', '.join(ls_param):}): {d_method['summary']:}"


def iterate_method_blocks(cryspy_object) -> Iterator[Tuple[str, object]]:
    """(path, object) of the data blocks, loops and single items of the
    object (not the items of the loops)."""
    l_obj = cryspy_object.items if isinstance(cryspy_object, GlobalN) else [cryspy_object]
    for obj in l_obj:
        if isinstance(obj, DataN):
            yield obj.get_name(), obj
            for block in obj.items:
                if isinstance(block, (ItemN, LoopN)):
                    yield f"{obj.get_name():}.{block.get_name():}", block
        elif isinstance(obj, (ItemN, LoopN)):
            yield obj.get_name(), obj


def find_method_blocks(cryspy_object, path: str = "", method: str = "") -> List[Tuple[str, object]]:
    """(path, object) of the blocks given by the path that have the
    method (any method when it is ""). The path is that of a data block,
    loop or item ('pnd', 'pnd.pd_background'; the data block may be given
    by its data name), or a class or category for all such blocks ('Pd',
    'pd_background'); empty for all blocks."""
    from cryspy.procedure_assistant.attributes import take_block
    l_block = list(iterate_method_blocks(cryspy_object))
    path_sh = path.strip()
    if path_sh != "":
        block = take_block(cryspy_object, path_sh)
        word = re.sub(r"[\s\-]+", "_", path_sh.lower())
        if block is not None:
            l_block = [(block_path, obj) for block_path, obj in l_block if obj is block]
        else:
            l_block = [(block_path, obj) for block_path, obj in l_block
                       if word in (type(obj).__name__.lower(), block_path.split(".")[-1].lower(),
                                   getattr(obj, "PREFIX", "").lower())]
    if method != "":
        l_block = [(block_path, obj) for block_path, obj in l_block if method in take_methods(type(obj))]
    return l_block


def take_data_block(cryspy_object, obj) -> Optional[DataN]:
    """The data block of the object which contains obj (or is obj)."""
    if isinstance(obj, DataN):
        return obj
    l_obj = cryspy_object.items if isinstance(cryspy_object, GlobalN) else [cryspy_object]
    for data in l_obj:
        if isinstance(data, DataN) and any(block is obj for block in data.items):
            return data
    return None


def find_parameter_blocks(cryspy_object, obj, param: str, types: tuple) -> Dict[str, object]:
    """{path: block} for a parameter given by classes: in the data block of
    obj first, then in the whole object. A base class as ItemN is found by
    the name of the parameter ('cell: ItemN' -> the cell)."""
    def fits(block_path: str, block) -> bool:
        if block is obj or not isinstance(block, types):
            return False
        if any(item in BASE_CLASSES for item in types):
            return block_path.split(".")[-1].lower() == param.lower()
        return True

    data = take_data_block(cryspy_object, obj)
    l_block = list(iterate_method_blocks(cryspy_object))
    if data is not None:
        d_found = {block_path: block for block_path, block in iterate_method_blocks(data)
                   if fits(block_path, block)}
        if d_found:
            return d_found
    return {block_path: block for block_path, block in l_block if fits(block_path, block)}
