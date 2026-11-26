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
import math
import re
from typing import Callable, Dict, Iterator, List, Optional, Tuple

from cryspy.B_parent_classes.cl_1_item import ItemN
from cryspy.B_parent_classes.cl_2_loop import LoopN
from cryspy.B_parent_classes.cl_3_data import DataN
from cryspy.B_parent_classes.cl_4_global import GlobalN
from cryspy.procedure_assistant.procedures import (
    D_BASIC_TYPE, form_question, parse_docstring, take_annotation_types, take_procedures)

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
            # the classes of unknown categories of a read file have no methods
            if sub_class not in l_class and not sub_class.__dict__.get("IS_GENERIC", False):
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


# words of a request that say nothing about the method
STOP_WORDS = {
    "the", "a", "an", "of", "to", "in", "on", "for", "by", "with", "and", "or",
    "this", "that", "these", "it", "its", "my", "me", "please", "crystal",
    "data", "block", "method", "through", "into", "from", "as", "be", "is",
    "do", "make", "can", "how", "what", "all"}
# abbreviations the user may type for words of the docstrings
WORD_ALIASES = {"ssg": "superspace", "msg": "magnetic", "mcif": "magnetic",
                "irrep": "representation", "irreps": "representation",
                "k": "propagation", "sc": "single"}


def take_request_words(text: str) -> list:
    """Meaningful words (stems of 5 letters) of a request to find methods."""
    l_word = []
    for word in re.split(r"[^a-z0-9]+", text.lower()):
        word = WORD_ALIASES.get(word, word)
        if len(word) < 2 or word in STOP_WORDS:
            continue
        stem = word[:5]
        if stem not in l_word:
            l_word.append(stem)
    return l_word


def score_method(function: Callable, l_word: list) -> int:
    """Number of request words found in the name and the docstring of a
    method (stems, so 'symmetry operations' finds 'symmetry operation')."""
    return score_text(function.__name__, inspect.getdoc(function) or "", l_word)


def score_text(name: str, text: str, l_word: list) -> int:
    """Number of request words (stems) found in a name and its description
    (a method and its docstring, a tool and its description)."""
    stems = take_stems(name, text)
    return sum(1 for word in l_word if word in stems)


# verbs of method names that do not tell what the method is about: 'Add
# the texture' must not find add_propagation_vector by the word 'add'
GENERIC_NAME_WORDS = {"add", "set", "remov", "delet", "defin", "creat", "calc", "show",
                      "chang", "updat", "run", "estim", "refin", "conve", "repor"}


def count_name_words(name: str, l_word: list) -> Tuple[int, int]:
    """Numbers of request words (stems) in a method or procedure name: the
    meaningful ones and the generic verbs of GENERIC_NAME_WORDS (which only
    break ties: 'Add the propagation vector' prefers add_propagation_vector
    to find_propagation_vectors)."""
    l_name = {word[:5] for word in name.lower().split("_") if word}
    l_found = [word for word in l_word if word in l_name]
    n_generic = sum(1 for word in l_found if word in GENERIC_NAME_WORDS)
    return len(l_found) - n_generic, n_generic


# a candidate is given when its score is at least this part of the score of
# the best one (the weak ones confuse a small model)
RELATIVE_SCORE_MIN = 0.8


def take_stems(name: str, text: str) -> set:
    """Stems (5 letters) of the words of a name and its description."""
    text = name.replace("_", " ") + " " + text
    return {word[:5] for word in re.split(r"[^a-z0-9]+", text.lower()) if word}


def take_description(function: Callable) -> str:
    """What a procedure or a method does: its docstring without the
    numpydoc sections (Parameters, Returns, Examples, ...), whose words
    describe the arguments, not the operation ('form factor' of the
    argument magnetic_ion of add_magnetic_model must not make it an
    answer to 'add the magnetic form factor through multipoles')."""
    doc = inspect.getdoc(function) or ""
    return re.split(r"\n[^\n]+\n-{3,}[ \t]*(?:\n|$)", "\n" + doc, maxsplit=1)[0].strip()


def take_category_text(item_class: type) -> str:
    """Description of a category: the docstring of its item class."""
    return inspect.getdoc(item_class) or ""


def take_category_classes() -> List[type]:
    """Classes of the categories (items and rows of loops) of CrysPy."""
    from cryspy.B_parent_classes.attributes import take_item_classes
    return take_item_classes()


_FREQUENCY_CACHE = {}
# constants of the BM25 ranking of text search: saturation of the number
# of words and the part of the length normalization
BM25_K1, BM25_B = 1.2, 0.75
# the words of a name tell what a procedure, method or category is for
# better than those of its description
NAME_WEIGHT = 1.5
# part of the weight of the request words that a candidate should have
COVERAGE_MIN = 0.35


def take_word_frequencies() -> Tuple[int, dict, float]:
    """Number of the procedures, methods and categories of CrysPy, the
    number of them with every stem in the name or the description, and the
    mean number of stems of a description (computed once)."""
    if "frequencies" not in _FREQUENCY_CACHE:
        l_stems = [take_stems(name, take_description(function))
                   for name, function in take_procedures().items()]
        d_method = {}
        for cls in take_classes_with_methods():
            d_method.update(take_methods(cls))
        l_stems.extend(take_stems(name, take_description(function))
                       for name, function in d_method.items())
        l_stems.extend(take_stems(item_class.PREFIX, take_category_text(item_class))
                       for item_class in take_category_classes())
        d_frequency = {}
        for stems in l_stems:
            for stem in stems:
                d_frequency[stem] = d_frequency.get(stem, 0) + 1
        n_mean = sum(len(stems) for stems in l_stems) / max(len(l_stems), 1)
        _FREQUENCY_CACHE["frequencies"] = (len(l_stems), d_frequency, n_mean)
    return _FREQUENCY_CACHE["frequencies"]


def take_word_weights(l_word: list) -> Dict[str, float]:
    """Weight of every request word: high for a rare word of the names and
    descriptions ('multipole'), low for a common one ('magnetic', 'atom',
    'add'), as the inverse document frequency of text search (BM25). A
    word of no description weighs as the rarest ones."""
    n_candidate, d_frequency, _ = take_word_frequencies()
    l_weight = []
    for word in l_word:
        n_word = max(d_frequency.get(word, 0), 1)
        l_weight.append(math.log(1. + (n_candidate - n_word + 0.5) / (n_word + 0.5)))
    return dict(zip(l_word, l_weight))


# names of the procedures and methods that only report, and the words of a
# request that asks for a report (a question)
REPORT_NAME_PREFIXES = ("report_", "show_", "list_", "print_", "summarize_")
REPORT_REQUEST_WORDS = {"report", "show", "list", "print", "summarize", "what", "which",
                        "how", "why", "is", "are", "does", "do", "tell"}
# factor of the score of a report for an order or of an operation for a
# question
KIND_MISMATCH = 0.7


def asks_for_report(request: str) -> bool:
    """True when the request is a question or asks for a report ('What is
    the distance between N1 and C1?', 'Show the magnetic broadening'),
    False for an order ('Set magnetic domain size')."""
    l_word = re.split(r"[^a-z0-9]+", request.lower().strip())
    return "?" in request or (len(l_word) > 0 and l_word[0] in REPORT_REQUEST_WORDS)


def covers_request(l_found: list, l_word: list, d_weight: dict) -> bool:
    """True when the found request words weigh at least COVERAGE_MIN of
    the request words that some description has (a label or a value as
    'Cu' or '0.5' says nothing about the operation)."""
    d_frequency = take_word_frequencies()[1]
    weight_known = sum(d_weight[word] for word in l_word if word in d_frequency or word in l_found)
    weight_found = sum(d_weight[word] for word in l_found)
    return weight_found >= COVERAGE_MIN * weight_known


def score_candidate(name: str, text: str, l_word: list, d_weight: dict,
                    flag_report: bool = False) -> Optional[float]:
    """Score of a procedure, method or category for a request.

    The sum of: the BM25 score of the request words in its name and
    description (a long docstring that has many words by chance weighs
    less); the weights of the words in its name times NAME_WEIGHT (half for
    the generic verbs of GENERIC_NAME_WORDS, which mostly break ties: 'Add
    the propagation vector' prefers add_propagation_vector to
    find_propagation_vectors); the weight of the first request word when
    the description starts with it (a numpydoc summary starts with the
    verb of the operation: 'Restrain the distance' prefers add_restraint,
    'Restrain a distance ...', to update_restraints). The score is
    multiplied by KIND_MISMATCH when a report (report_, show_, list_)
    answers an order or an operation answers a question (flag_report of
    asks_for_report): 'Set magnetic domain size' prefers
    add_magnetic_broadening to report_magnetic_broadening. None when it does
    not fit: fewer than two request words (one for a one-word request) in
    the name and the description, no meaningful one in the name (or,
    when there is none in the name, two in the summary, the first
    sentence), or the found words weigh too little (covers_request)."""
    n_name, _ = count_name_words(name, l_word)
    if n_name == 0:
        # a short name ('atom_site_scat') may not have the words: two of
        # them in its summary then say that it is for the request
        n_summary = count_name_words("_".join(take_stems("", take_summary(text))), l_word)[0]
        n_name = 1 if n_summary >= 2 else 0
    stems = take_stems(name, text)
    l_found = [word for word in l_word if word in stems]
    n_min = 1 if len(l_word) == 1 else 2
    if n_name == 0 or len(l_found) < n_min:
        return None
    # the words of no description (labels, values: 'Cu', '0.5') say
    # nothing about the operation; of the others, the found ones should
    # weigh enough: 'add' and 'atom' do not make add_atom an answer to
    # 'Add the magnetic form factor through multipoles for the atom Cu'
    if not covers_request(l_found, l_word, d_weight):
        return None
    n_mean = take_word_frequencies()[2]
    norm = (BM25_K1 + 1.) / (1. + BM25_K1 * (1. - BM25_B + BM25_B * len(stems) / n_mean))
    score = norm * sum(d_weight[word] for word in l_found)
    l_name = {word[:5] for word in name.lower().split("_") if word}
    for word in l_word:
        if word in l_name:
            score += d_weight[word] * NAME_WEIGHT * (0.5 if word in GENERIC_NAME_WORDS else 1.)
    l_text_word = re.split(r"[^a-z0-9]+", text.lower().strip(), maxsplit=1)
    if l_text_word[0][:5] == l_word[0]:
        score += d_weight[l_word[0]]
    if name.lower().startswith(REPORT_NAME_PREFIXES) != flag_report:
        score *= KIND_MISMATCH
    return score


def find_matching_procedures(request: str, max_procedures: int = 3) -> List[dict]:
    """Procedures (action run_procedure) that fit a request, best first, by
    the rules of :func:`score_candidate`: 'Add spherical texture
    description into powder experiment' finds
    add_texture_spherical_harmonics.

    Returns
    -------
    list of dict
        {"name", "summary", "score"}.
    """
    l_word = take_request_words(request)
    if len(l_word) == 0:
        return []
    d_weight = take_word_weights(l_word)
    flag_report = asks_for_report(request)
    l_found = []
    for name, function in take_procedures().items():
        doc = (inspect.getdoc(function) or "").strip()
        score = score_candidate(name, take_description(function), l_word, d_weight, flag_report)
        if score is None:
            continue
        l_found.append({"name": name, "summary": doc.split("\n")[0] if doc else "",
                        "score": score, "function": function})
    l_found.sort(key=lambda item: (-item["score"], item["name"]))
    return l_found[:max_procedures]


def find_matching_methods(cryspy_object, request: str, max_methods: int = 3) -> List[dict]:
    """Methods of the object that fit a request of the user, best first,
    by the rules of :func:`score_candidate`: 'Convert the structure to
    superspace' finds convert_to_superspace of the crystals. A program
    with a small language model gives them with the request, so that the
    model does not answer "cannot do" for an operation that is a method
    (the model often does not look the methods up by itself).

    Returns
    -------
    list of dict
        {"name", "summary", "paths" (blocks that have it), "score"}.
    """
    l_word = take_request_words(request)
    if len(l_word) == 0:
        return []
    d_weight = take_word_weights(l_word)
    flag_report = asks_for_report(request)
    d_found = {}
    for block_path, block in iterate_method_blocks(cryspy_object):
        for name, function in take_methods(type(block)).items():
            if name in d_found:
                d_found[name]["paths"].append(block_path)
                continue
            score = score_candidate(name, take_description(function), l_word, d_weight, flag_report)
            if score is None:
                continue
            d_found[name] = {"name": name, "summary": describe_method(function)["summary"],
                             "paths": [block_path], "score": score,
                             "function": function, "block": block}
    l_method = sorted(d_found.values(), key=lambda item: (-item["score"], item["name"]))
    return l_method[:max_methods]


def find_matching_categories(cryspy_object, request: str, max_categories: int = 3) -> List[dict]:
    """Categories (items or rows of loops, action add_item) that the data
    blocks of the object may have and that fit a request, best first, by
    the rules of :func:`score_candidate` on the category name and the
    docstring of its class: 'Add description of magnetic form factor
    through multipoles for Cu' finds atom_rho_multipole. Many requests
    ('add multipoles', 'add the local axes of Fe') are answered by a new
    row of a category, not by a procedure or a method.

    Returns
    -------
    list of dict
        {"name" (category), "summary", "paths" (data blocks), "score"}.
    """
    l_word = take_request_words(request)
    if len(l_word) == 0:
        return []
    d_weight = take_word_weights(l_word)
    flag_report = asks_for_report(request)
    d_found = {}
    for block in iterate_data_blocks(cryspy_object):
        for item_class in take_block_categories(type(block)):
            name = item_class.PREFIX
            if name in d_found:
                d_found[name]["paths"].append(block.get_name())
                continue
            text = take_category_text(item_class)
            # a category is offered for add_item: 'Add multipoles for Cu'
            # has the word 'add' of its operation
            score = score_candidate(name, f"{text:} Add", l_word, d_weight, flag_report)
            if score is None:
                continue
            d_found[name] = {"name": name, "summary": take_summary(text),
                             "paths": [block.get_name()], "score": score,
                             "category_class": item_class, "block": block}
    l_category = sorted(d_found.values(), key=lambda item: (-item["score"], item["name"]))
    return l_category[:max_categories]


def iterate_data_blocks(cryspy_object) -> Iterator[DataN]:
    if isinstance(cryspy_object, GlobalN):
        for item in cryspy_object.items:
            if isinstance(item, DataN):
                yield item
    elif isinstance(cryspy_object, DataN):
        yield cryspy_object


def take_block_categories(data_class: type) -> List[type]:
    """Item classes of the items and of the rows of the loops that a data
    block of the class may have."""
    l_class = []
    for block_class in data_class.CLASSES:
        if issubclass(block_class, LoopN):
            block_class = getattr(block_class, "ITEM_CLASS", None)
        if (isinstance(block_class, type) and issubclass(block_class, ItemN)
                and isinstance(getattr(block_class, "PREFIX", None), str)
                and block_class not in l_class):
            l_class.append(block_class)
    return l_class


def take_summary(text: str) -> str:
    """First sentence of a description, on one line."""
    text = " ".join(text.split())
    match = re.match(r"(.+?\.)(\s|$)", text)
    return match.group(1) if match else text


def form_matching_methods_text(cryspy_object, request: str, max_methods: int = 3) -> str:
    """Text of :func:`find_matching_procedures`,
    :func:`find_matching_methods` and :func:`find_matching_categories` for
    the prompt of a language model ('' when nothing fits). Only the
    candidates close to the best one are given (the weak ones confuse a
    small model), the group of the best one first."""
    d_group = {
        "procedure": find_matching_procedures(request, max_methods),
        "method": find_matching_methods(cryspy_object, request, max_methods),
        "category": find_matching_categories(cryspy_object, request, max_methods)}
    l_all = [item for l_item in d_group.values() for item in l_item]
    if len(l_all) == 0:
        return ""
    score_min = RELATIVE_SCORE_MIN * max(item["score"] for item in l_all)
    d_text = {}
    for kind, l_item in d_group.items():
        l_item = [item for item in l_item if item["score"] >= score_min]
        if len(l_item) == 0:
            continue
        d_text[kind] = (l_item[0]["score"], form_candidates_text(kind, l_item))
    l_kind = sorted(d_text, key=lambda kind: -d_text[kind][0])
    text = "\n".join(d_text[kind][1] for kind in l_kind)
    return text + form_best_requirements_text(cryspy_object, d_group[l_kind[0]][0])


def form_best_requirements_text(cryspy_object, d_best: dict) -> str:
    """What the object has and misses for the best candidate when its
    docstring names the categories it needs and adds (requirements.py), so
    that the model gives only the missing values and asks about the
    others: '' otherwise."""
    from cryspy.procedure_assistant.requirements import (
        check_requirements, form_requirements_text, has_requirements)
    function = d_best.get("function")
    if function is None and d_best.get("category_class") is not None:
        return form_category_rows_text(cryspy_object, d_best)
    if function is None or not has_requirements(function):
        return ""
    d_check = check_requirements(cryspy_object, function, d_best.get("block"))
    text = ("\n" + form_requirements_text(d_check, d_best["name"], flag_short=True)
            + "\nGive only the values that are not in the object; ask the user "
              "(ask_user) when several values are possible and the request does not say which.")
    return text + form_generic_route_text(d_check)


def form_generic_route_text(d_check: dict) -> str:
    """The same rows by add_item (the category of the first keyed row the
    method adds), for a user who gives the values himself."""
    from cryspy.procedure_assistant.dependencies import take_key_attribute
    from cryspy.procedure_assistant.items import take_category_class
    for d_add in d_check["adds"]:
        if not d_add["column"]:
            continue
        try:
            item_class = take_category_class(d_add["category"])
        except ValueError:
            continue
        if take_key_attribute(item_class) == "":
            continue
        return (f"\nWithout the method: add_item {{\"category\": \"{d_add['category']}\"}} adds "
                "the rows with the rows they need and asks the values it cannot find.")
    return ""


def form_candidates_text(kind: str, l_item: List[dict]) -> str:
    """Lines of the procedures, methods or categories for the prompt."""
    if kind == "procedure":
        ls_out = ["Procedures that may do this (action run_procedure):"]
        for d_item in l_item:
            ls_out.append(f"- {d_item['name']}: {d_item['summary']} "
                          f"[run_procedure {{\"name\": \"{d_item['name']}\"}}]")
        return "\n".join(ls_out)
    if kind == "method":
        ls_out = ["Methods of the opened file that may do this (action call_method):"]
    else:
        ls_out = ["Categories that may describe this, a new row or item (action add_item):"]
    for d_item in l_item:
        l_path = d_item["paths"]
        s_path = ", ".join(l_path[:3]) + (", ..." if len(l_path) > 3 else "")
        if kind == "method":
            ls_out.append(f"- {d_item['name']} of {s_path}: {d_item['summary']} "
                          f"[call_method {{\"method\": \"{d_item['name']}\"}}]")
        else:
            ls_out.append(f"- {d_item['name']} in {s_path}: {d_item['summary']} "
                          f"[add_item {{\"category\": \"{d_item['name']}\"}}]")
    return "\n".join(ls_out)


def form_category_rows_text(cryspy_object, d_best: dict) -> str:
    """What the data block has of the best category and of the categories
    its rows need (dependencies.py): add_item adds the missing ones too."""
    from cryspy.procedure_assistant.dependencies import describe_category_rows, take_key_attribute
    item_class = d_best["category_class"]
    if take_key_attribute(item_class) == "":
        return ""
    return ("\n" + describe_category_rows(cryspy_object, d_best["block"], item_class)
            + "\nadd_item adds the rows for the given keys (one or several) with the rows they "
              "need; values not known are asked from the user.")
