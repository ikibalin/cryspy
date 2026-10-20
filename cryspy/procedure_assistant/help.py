"""Help pages of CrysPy for people and for AI assistants.

Two kinds of pages:

- guides: Markdown files in help_pages/ written by hand (how to describe
  thermal vibrations, multipoles, a powder experiment, ...);
- categories: one page per CIF category, generated from the classes of
  CrysPy (description, tags, refinable parameters, allowed values).

Search is a plain keyword search without external dependencies, so it works
the same in a script, in CrysPy Editor and with a local language model.
"""
import os
import re
from typing import Dict, List, Optional

DIR_HELP_PAGES = os.path.join(os.path.dirname(__file__), "help_pages")

# words which do not help to find a page
SET_STOP_WORDS = {
    "a", "an", "and", "are", "be", "by", "can", "do", "for", "from", "how",
    "i", "in", "is", "it", "me", "my", "of", "on", "or", "should", "the",
    "to", "use", "used", "what", "which", "with",
    "в", "и", "к", "как", "какие", "какой", "мне", "на", "надо", "не",
    "нужно", "о", "по", "с", "что", "чтобы", "это", "для", "использовать",
    "необходимо", "описать", "описания",
}


class HelpPage:
    def __init__(self, name: str, title: str, text: str,
                 keywords: tuple = (), categories: tuple = (),
                 kind: str = "guide"):
        self.name = name
        self.title = title
        self.text = text
        self.keywords = keywords
        self.categories = categories
        self.kind = kind

    def to_text(self) -> str:
        return f"# {self.title}\n\n{self.text.strip()}\n"


def split_front_matter(content: str):
    """Split '---' header with 'key: value' lines from the Markdown body."""
    d_header = {}
    if content.startswith("---"):
        part_header, _, body = content[3:].partition("\n---")
        for line in part_header.strip().splitlines():
            key, _, value = line.partition(":")
            d_header[key.strip()] = value.strip()
        content = body.lstrip("-").lstrip("\n")
    return d_header, content


def split_list(text: str) -> tuple:
    return tuple(word.strip().lower() for word in text.split(",")
                 if word.strip() != "")


def load_guide_pages() -> List[HelpPage]:
    l_page = []
    for file_name in sorted(os.listdir(DIR_HELP_PAGES)):
        if not file_name.endswith(".md"):
            continue
        with open(os.path.join(DIR_HELP_PAGES, file_name), "r",
                  encoding="utf-8") as fid:
            d_header, body = split_front_matter(fid.read())
        name = file_name[:-3]
        l_page.append(HelpPage(
            name, d_header.get("title", name), body,
            keywords=split_list(d_header.get("keywords", "")),
            categories=split_list(d_header.get("categories", "")),
            kind="guide"))
    return l_page


def form_category_text(item_class, l_data_class) -> str:
    """Reference text of one CIF category generated from its class."""
    prefix = item_class.PREFIX
    ls_out = []
    doc = (item_class.__doc__ or "").strip()
    if doc != "":
        ls_out.append(doc)
    l_parent = [cls.PREFIX for cls in l_data_class
                if any(getattr(cls_item, "ITEM_CLASS", cls_item) is item_class
                       for cls_item in cls.CLASSES)]
    if len(l_parent) > 0:
        ls_out.append(f"Used in data blocks: {', '.join(l_parent):}.")
    if len(item_class.ATTR_MANDATORY_CIF) > 0:
        ls_out.append("Mandatory tags: " + ", ".join(
            f"_{prefix:}_{name:}" for name in item_class.ATTR_MANDATORY_CIF))
    if len(item_class.ATTR_OPTIONAL_CIF) > 0:
        ls_out.append("Optional tags: " + ", ".join(
            f"_{prefix:}_{name:}" for name in item_class.ATTR_OPTIONAL_CIF))
    if len(item_class.ATTR_REF) > 0:
        ls_out.append("Refinable parameters: " + ", ".join(item_class.ATTR_REF)
                      + ". Mark them by '()' after the value.")
    d_constraints = getattr(item_class, "D_CONSTRAINTS", {})
    for name, values in d_constraints.items():
        if name.lower() in [hh.lower() for hh in item_class.ATTR_NAMES]:
            ls_out.append(f"Allowed values of {name:}: "
                          + ", ".join(str(value) for value in values))
    d_default = {key: value for key, value in
                 getattr(item_class, "D_DEFAULT", {}).items()
                 if key in item_class.ATTR_NAMES}
    if len(d_default) > 0:
        ls_out.append("Default values: " + ", ".join(
            f"{key:} = {value:}" for key, value in d_default.items()))
    return "\n\n".join(ls_out)


def load_category_pages() -> List[HelpPage]:
    from cryspy.H_functions_global.function_1_cryspy_objects import (
        L_ITEM_CLASS, L_DATA_CLASS)
    from cryspy.C_item_loop_classes.cl_1_atom_rho_multipole import \
        AtomRhoMultipole
    l_item_class = list(L_ITEM_CLASS)
    # not registered in L_ITEM_CLASS, but used by the crystal data block
    if AtomRhoMultipole not in l_item_class:
        l_item_class.append(AtomRhoMultipole)
    l_page = []
    for item_class in l_item_class:
        prefix = item_class.PREFIX.lower()
        l_page.append(HelpPage(
            prefix, f"Category {prefix:}",
            form_category_text(item_class, L_DATA_CLASS),
            keywords=tuple(prefix.split("_")) + (prefix, ),
            categories=(prefix, ), kind="category"))
    for data_class in L_DATA_CLASS:
        prefix = data_class.PREFIX.lower()
        ls_text = [(data_class.__doc__ or "").strip()]
        ls_text.append("Mandatory categories: " + ", ".join(
            take_class_prefix(cls) for cls in data_class.CLASSES_MANDATORY))
        ls_text.append("Optional categories: " + ", ".join(
            take_class_prefix(cls) for cls in data_class.CLASSES_OPTIONAL))
        l_page.append(HelpPage(
            prefix, f"Data block {prefix:}", "\n\n".join(ls_text),
            keywords=(prefix, "data", "block"), categories=(prefix, ),
            kind="category"))
    return l_page


def take_class_prefix(cls) -> str:
    if "PREFIX" in cls.__dict__:
        return cls.PREFIX
    item_class = getattr(cls, "ITEM_CLASS", None)
    if item_class is not None and hasattr(item_class, "PREFIX"):
        return item_class.PREFIX
    return cls.__name__


_L_PAGES: List[HelpPage] = []


def take_help_pages() -> List[HelpPage]:
    """Guides first, then category pages; loaded once."""
    if len(_L_PAGES) == 0:
        _L_PAGES.extend(load_guide_pages())
        _L_PAGES.extend(load_category_pages())
    return _L_PAGES


def split_words(text: str) -> List[str]:
    return [word for word in re.findall(r"\w+", text.lower())
            if word not in SET_STOP_WORDS]


def is_word_match(word: str, l_word: List[str]) -> bool:
    """Same word or the same stem (endings differ in Russian and English)."""
    for word_2 in l_word:
        if word == word_2:
            return True
        n_stem = min(len(word), len(word_2)) - 2
        if n_stem >= 4 and word[:n_stem] == word_2[:n_stem]:
            return True
    return False


def score_help_page(page: HelpPage, l_query_word: List[str],
                    query: str) -> float:
    l_keyword_word = split_words(" ".join(page.keywords))
    l_title_word = split_words(page.title)
    l_category_word = split_words(" ".join(page.categories))
    l_text_word = split_words(page.text)
    score = 0.
    for word in l_query_word:
        if is_word_match(word, l_keyword_word):
            score += 3.
        if is_word_match(word, l_title_word):
            score += 2.
        if is_word_match(word, l_category_word):
            score += 1.
        if word in l_text_word:
            score += min(l_text_word.count(word), 5) * 0.2
    # a whole keyword phrase or category name inside the query
    query_sh = query.lower()
    for keyword in page.keywords + page.categories:
        if len(keyword) > 3 and keyword in query_sh:
            score += 3.
    if page.kind == "guide":
        score *= 1.5
    return score


def search_help(query: str, n_pages: int = 3) -> List[HelpPage]:
    """Help pages that fit the query best, the best first."""
    l_query_word = split_words(query)
    if len(l_query_word) == 0:
        return []
    l_score = []
    for page in take_help_pages():
        score = score_help_page(page, l_query_word, query)
        if score > 0:
            l_score.append((score, page))
    l_score.sort(key=lambda pair: -pair[0])
    return [page for _, page in l_score[:n_pages]]


def take_help_page(name: str) -> Optional[HelpPage]:
    name_sh = name.strip().lower()
    for page in take_help_pages():
        if page.name == name_sh:
            return page
    return None


def help_topics() -> str:
    """List of help guides and categories of CrysPy."""
    l_page = take_help_pages()
    ls_out = ["Guides:"]
    ls_out.extend(f"- {page.name:}: {page.title:}" for page in l_page
                  if page.kind == "guide")
    ls_out.append("Categories (data blocks, items and loops): " + ", ".join(
        page.name for page in l_page if page.kind == "category"))
    return "\n".join(ls_out)


def get_help(query: str, max_length: int = 6000) -> str:
    """Help text for a topic name, a category or a free question.

    Examples:
        cryspy.get_help("multipoles")
        cryspy.get_help("atom_site_aniso")
        cryspy.get_help("how to describe thermal vibrations")
    """
    page = take_help_page(query)
    if page is not None:
        text = page.to_text()
    else:
        l_page = search_help(query)
        if len(l_page) == 0:
            return ("Nothing is found in CrysPy help for "
                    f"'{query:}'.\n\n{help_topics():}")
        text = "\n\n".join(page.to_text() for page in l_page)
    if len(text) > max_length:
        text = text[:max_length].rsplit("\n", 1)[0] + "\n[...]"
    return text
