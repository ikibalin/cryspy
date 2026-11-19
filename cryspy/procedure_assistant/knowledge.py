"""Search in the internal knowledge base of CrysPy (cryspy/knowledge/*.md).

The knowledge base is a folder of Markdown files inside the package:
physics that CrysPy and an AI assistant need besides the description of the
rcif format (help_pages): the magnetic ions (usual valences, Lande factors,
free-ion moments), the usual charges of the ions for the charge balance,
... Some files hold tables read by CrysPy itself
(A_functions_base/knowledge_tables.py).

The search is that of the documents of the user (documents.py): the files
are cut into fragments, ranked by BM25 and, when the program gives an
embedding function (set_document_embedding), merged with the ranking by
meaning. The index is kept in memory (the folder is small and read-only).
Every row of a table is a fragment of its own, written with the titles of
the columns and the title of the file ("Magnetic rare-earth (4f) ions:
Ion Tb3+ (Tb); Usual yes; ..."), so that a question about one ion finds its
row.
"""
import os
from typing import List

from cryspy.A_functions_base.knowledge_tables import DIR_KNOWLEDGE, split_ion, split_table_row
from cryspy.procedure_assistant import documents
from cryspy.procedure_assistant.help import split_front_matter

D_KNOWLEDGE = {"index": None}


def form_row_text(title: str, l_column: List[str], l_cell: List[str]) -> str:
    """'Title: Ion Tb3+ (Tb); Usual yes; Lande g 3/2' (empty cells left out;
    the element of an ion is added so that 'Tb' finds 'Tb3+')."""
    l_text = []
    for column, cell in zip(l_column, l_cell):
        if cell == "":
            continue
        if column in ("Ion", "Element"):
            element = split_ion(cell)[0]
            cell = cell if element == cell else f"{cell} ({element})"
        l_text.append(f"{column} {cell}")
    return f"{title}: " + "; ".join(l_text)


def split_knowledge_chunks(content: str) -> List[tuple]:
    """(line, text) fragments of a knowledge file: the text outside the
    tables as documents.split_text_chunks, every table row by itself."""
    d_header, _ = split_front_matter(content)
    l_line = content.splitlines()
    title = d_header.get("title", "")
    l_out, l_text, ind = [], [], 0
    while ind < len(l_line):
        line = l_line[ind].strip()
        line_next = l_line[ind + 1].strip() if ind + 1 < len(l_line) else ""
        if line.startswith("|") and line_next.startswith("|") and set(line_next) <= set("|-: "):
            l_column = split_table_row(line)
            ind += 2
            while ind < len(l_line) and l_line[ind].strip().startswith("|"):
                l_out.append((ind + 1, form_row_text(title, l_column, split_table_row(l_line[ind]))))
                ind += 1
            l_text.append("")
            continue
        l_text.append(l_line[ind])
        ind += 1
    for line, text in documents.split_text_chunks("\n".join(l_text)):
        l_out.append((line, text))
    return sorted(l_out)


def take_knowledge_index() -> documents.DocumentIndex:
    """Index of the knowledge base, made once (the vectors of an embedding
    model are kept in the cache folder of the user, as for documents)."""
    index = D_KNOWLEDGE["index"]
    if index is None:
        index = documents.DocumentIndex(DIR_KNOWLEDGE)
        index.d_file = {}
        for name in documents.find_document_files(DIR_KNOWLEDGE):
            with open(os.path.join(DIR_KNOWLEDGE, name), "r", encoding="utf-8") as fid:
                l_chunk = split_knowledge_chunks(fid.read())
            index.d_file[name] = {"chunks": [[0, line, text] for line, text in l_chunk],
                                  "error": ""}
        index.form_chunks()
        D_KNOWLEDGE["index"] = index
    index.load_vectors(documents.take_document_embedding_name())
    return index


def list_knowledge_files() -> List[str]:
    """Files of the knowledge base (relative to cryspy/knowledge)."""
    return sorted(take_knowledge_index().d_file)


def search_knowledge(query: str, n_results: int = 4) -> List[dict]:
    """Fragments of the knowledge base of CrysPy that fit the query best.

    Parameters
    ----------
    query : str
        The question or the key words (as "usual valence of Tb", "Lande
        factor Ho3+", "charge balance").
    n_results : int
        Number of fragments.

    Returns
    -------
    list of dict
        As documents.search_documents: {"file", "location", "text",
        "score", ...}, the best first.
    """
    d_embed = documents.D_DOCUMENTS
    return take_knowledge_index().search(query, n_results=n_results, embed=d_embed["embed"],
                                         embed_query=d_embed["embed_query"])
